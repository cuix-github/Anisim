"""No-force regression, numpy lockstep, and invalid-input rejection.

Writes separate evidence directories. Does not run a production scene.
"""
import json
import os
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "gates"), str(ROOT / "reference2d")]
from validate_3d import run_node
from validate_step import channels
from validate_level_set_node import flat_inputs, save, load
from validate_wave_zone import compare
from pfflip2d import Sim
from external_forces import apply_external

FIELDS = tuple(f"{family}_{axis}" for family in ("liquid_acceleration", "linear_drag") for axis in "uvw")
OUT = ROOT / "out/runs/force_wave_validation_20261020"


def main():
    sim = Sim(80, 32)
    sim.seed(lambda x, y: y < 12)
    sim.calibrate()
    sim.run(4, .02)
    ch = {key: value.astype(np.float32) for key, value in channels(sim).items()}
    arr = flat_inputs(sim, ch)
    dt = float(np.float32(.02))
    x, y = np.meshgrid(np.arange(81), np.arange(32) + .5, indexing="ij")
    au = (4 * np.exp(-((x - 25) / 6) ** 2) * np.cos(y / 7)).astype(np.float32)
    x, y = np.meshgrid(np.arange(80) + .5, np.arange(33), indexing="ij")
    av = (2 * np.sin(x / 9) * np.exp(-((y - 12) / 5) ** 2)).astype(np.float32)
    arrays = dict(liquid_acceleration_u=au.ravel(), liquid_acceleration_v=av.ravel(),
                  liquid_acceleration_w=np.zeros(1, np.float32),
                  linear_drag_u=np.array([.3], np.float32),
                  linear_drag_v=np.array([.1], np.float32), linear_drag_w=np.zeros(1, np.float32))
    idir = OUT / "fields_in"
    for key, val in dict(arr, **arrays).items():
        save(str(idir), key, val)
    kwargs = dict(nx=80, ny=32, nz=1, dt=dt, gravity=-9.8, preconditioner=1)
    off = OUT / "new_off"
    run_node(str(idir), str(off), extra_arrays=FIELDS, external_forces=0, **kwargs)
    plain = OUT / "new_plain"
    run_node(str(idir), str(plain), **kwargs)
    for key in ("out_u", "out_v", "out_w", "pressure", "out_positions", "out_velocities"):
        assert np.array_equal(load(str(off), key), load(str(plain), key)), key

    ref = Sim(80, 32)
    ref.pos = arr["positions"][:, :2].astype(float)
    ref.vel = arr["velocities"][:, :2].astype(float)
    ref.typ = sim.typ.copy()
    ref.rho0_face = sim.rho0_face
    um, uo = ch["u_mass"].astype(float), ch["u_mom"].astype(float)
    vm, vo = ch["v_mass"].astype(float), ch["v_mom"].astype(float)
    ref.u = np.where(um > 0, uo / np.maximum(um, 1e-30), 0)
    ref.v = np.where(vm > 0, vo / np.maximum(vm, 1e-30), 0)
    pu, pv = np.clip(ch["u_phase"].astype(float), 0, 1), np.clip(ch["v_phase"].astype(float), 0, 1)
    ref.mu, ref.mv = np.empty_like(ref.u), np.empty_like(ref.v)
    ref.phase = lambda raw: pu if raw is ref.mu else pv
    ref._enforce_walls()
    us, vs = ref.u.copy(), ref.v.copy()
    ref.v += float(np.float32(-9.8)) * dt
    ref.u = apply_external(ref.u, pu, au, float(np.float32(.3)), dt)
    ref.v = apply_external(ref.v, pv, av, float(np.float32(.1)), dt)
    ref._enforce_walls()
    ref.project(dt)
    ref.g2p(us, vs)
    ref.advect(dt)
    on = OUT / "new_on"
    stats = run_node(str(idir), str(on), extra_arrays=FIELDS, external_forces=1, **kwargs)
    good, comparison = compare("external forces before projection", str(on), ref)
    assert good, comparison
    assert np.max(np.abs(load(str(on), "out_velocities") - load(str(off), "out_velocities"))) > .001
    assert stats["final_residual"] <= 1.05e-6
    bad = dict(arrays, liquid_acceleration_u=np.array([np.nan], np.float32))
    badin = OUT / "invalid_in"
    for key, val in dict(arr, **bad).items():
        save(str(badin), key, val)
    rejected = False
    try:
        run_node(str(badin), str(OUT / "invalid_out"), extra_arrays=FIELDS, external_forces=1, **kwargs)
    except SystemExit as error:
        rejected = "invalid external" in str(error)
    assert rejected
    receipt = dict(passed=True, off_bit_identical=True, lockstep=comparison,
                   statistics=stats, nonfinite_rejected=rejected)
    (OUT / "force_gate.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
