"""The wave zone in step_3d against its spec, the 2D reference's Sim(zone=...): the
wavemaker of a numerical wave basin without a paddle, a relaxation zone at the x = 0
wall where the water's velocity is blended toward the sum of linear wave trains after
every projection (the faces) and after the FLIP update (the liquid particles). Four
checks, one-metre cells:

  off       the zone ports at 0 (and a `waves` array of nothing) are the plain step,
            every output bit for bit (one step of a still tank at nz = 1).
  lockstep  one step from one state, the node at nz = 1 against the reference's twin
            step at the same instant: a tank the zone has driven for 45 steps, one
            train of 2 m at 4 s, full from the start. Grid velocities, positions and
            velocities within the tolerances of the level-set lockstep.
  train     the dispersion check through the whole loop (the fused graph, nz = 2):
            600 x 40 x 2, 20 m deep, the zone 24 m long, one train of 2 m at 8 s
            ramped over 8 s, 50 s; the surface at 50 s from 100 to 400 m beyond the
            zone against linear theory (wavelength within 10 percent, height within 30
            percent of the train's own, the zone writing the height directly), the
            still level's drift recorded.
  oblique   the direction, in 3D under the production settings (ST, escaped
            particles, adaptive air, with the experiments' calibration ports): 160 x 20 x 96, 10 m deep, one train of 1 m at
            5 s heading 20 degrees toward +z, 15 s; the crest lines' tilt across the
            basin, from the shift between surface profiles 24 m apart, within 5 degrees
            of the train's heading, the surface finite and no higher than three heights.

Run from the repo root:  .venv\\Scripts\\python.exe gates\\validate_wave_zone.py
[off|lockstep|train|oblique]; about ten minutes on every thread. Exits nonzero on failure.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "reference2d"))
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
from pfflip2d import Sim  # noqa: E402
from validate_piston import wavelength_and_height  # noqa: E402
from validate_piston_node import surface_top  # noqa: E402
from validate_step import channels  # noqa: E402
from validate_3d import run_node  # noqa: E402
from validate_level_set_node import DT, flat_inputs, load, save  # noqa: E402
from loop3d import run_bifrost3  # noqa: E402
sys.path.insert(0, os.path.join(HERE, "..", "experiments"))
from solitary_beach import calibration  # noqa: E402

OUT = os.path.join(HERE, "..", "out", "wave_zone")
RESULTS = os.path.join(HERE, "..", "results", "3d", "wave_zone")
G = 9.8


def f32(v):
    return float(np.float32(v))


def ports(length, depth, ramp):
    return {"zone_length": f32(length), "zone_depth": f32(depth), "zone_ramp": f32(ramp)}


def waves_array(trains):
    """The node's `waves`: four floats a train (height, period, direction in degrees,
    phase)."""
    return np.asarray([v for tr in trains for v in tr], np.float32).reshape(-1)


def spec_zone(length, depth, ramp, trains):
    """The reference's zone from the node's float32 ports (its trains run along x)."""
    return dict(length=f32(length), depth=f32(depth), ramp=f32(ramp),
                waves=[(f32(H), f32(T), f32(ph)) for H, T, _, ph in trains])


def node_step(tag, sim, step_index, trains=None, zone=None, **extra):
    """One step_3d call at nz = 1 on the reference state; trains and the zone's
    scalars when given."""
    ch = channels(sim)
    arr = flat_inputs(sim, ch)
    if trains is not None:
        arr["waves"] = waves_array(trains)
    d_in, d_out = os.path.join(OUT, tag + "_in"), os.path.join(OUT, tag + "_out")
    for k, v in arr.items():
        save(d_in, k, v)
    stats = run_node(d_in, d_out, extra_arrays=("waves",) if trains is not None else (),
                     nx=sim.nx, ny=sim.ny, nz=1, dt=DT, preconditioner=1, gravity=-G,
                     step_index=step_index, **(zone or {}), **extra)
    return d_out, ch, stats


def twin(state, ch, zone, t):
    """The reference's step from the projection on, as Sim.step does it with a zone and
    nothing else on, from the node's inputs (validate_piston_node.piston_twin)."""
    ref = Sim(state.nx, state.ny, zone=zone)
    ref.t = t
    ref._zone_t = t + DT             # the inlet faces carry the trains' velocity
    arr = flat_inputs(state, ch)
    ref.pos = arr["positions"][:, :2].astype(float)
    ref.vel = arr["velocities"][:, :2].astype(float)
    ref.typ = state.typ.copy()
    ref.rho0_face = state.rho0_face
    um, uo = ch["u_mass"].astype(float), ch["u_mom"].astype(float)
    vm, vo = ch["v_mass"].astype(float), ch["v_mom"].astype(float)
    ref.u = np.where(um > 0, uo / np.maximum(um, 1e-30), 0.0)
    ref.v = np.where(vm > 0, vo / np.maximum(vm, 1e-30), 0.0)
    # the node's face phase, which the projection and the zone's faces (the inlet's
    # among them) are weighed by
    phi_u = np.clip(ch["u_phase"].astype(float), 0, 1)
    phi_v = np.clip(ch["v_phase"].astype(float), 0, 1)
    ref.mu, ref.mv = np.empty_like(ref.u), np.empty_like(ref.v)
    ref.phase = lambda raw: phi_u if raw is ref.mu else phi_v
    ref._enforce_walls()
    u_star, v_star = ref.u.copy(), ref.v.copy()
    ref.v += ref.g * DT
    ref._enforce_walls()
    ref.project(DT)
    ref._zone_faces(t + DT)
    ref._enforce_walls()
    ref.g2p(u_star, v_star)
    ref._zone_particles(t + DT)
    ref.advect(DT)
    return ref


# ------------------------------------------------------------------ off


def check_off():
    print("=== off: the zone ports at 0 are the plain step, bit for bit", flush=True)
    sim = Sim(120, 40)
    sim.seed(lambda x, y: y < 20.0)
    sim.calibrate()
    for _ in range(10):
        sim.step(DT)
    d_z, _, s_z = node_step("off_zone", sim, 10, trains=[(0.0, 0.0, 0.0, 0.0)], zone=ports(0.0, 20.0, 0.0))
    d_p, _, s_p = node_step("off_plain", sim, 10)
    worst = 0.0
    for name in ("out_u", "out_v", "out_w", "pressure", "out_positions", "out_velocities"):
        a, b = load(d_z, name), load(d_p, name)
        worst = max(worst, float(np.abs(a.astype(np.float64) - b).max()) if a.shape == b.shape else np.inf)
    same = worst == 0.0 and s_z["iterations_used"] == s_p["iterations_used"]
    print(f"    one step, {len(sim.pos):,} particles: max |diff| {worst:.1e} over the grid, pressure, positions "
          f"and velocities; iterations {s_z['iterations_used']:.0f} and {s_p['iterations_used']:.0f}  "
          + ("OK" if same else "FAIL"), flush=True)
    return same, {"max_diff": worst}


# ------------------------------------------------------------------ lockstep


def check_lockstep():
    print("=== lockstep: one step against the reference's twin, nz = 1", flush=True)
    nx, ny, h = 120, 40, 20.0
    length, ramp = 24.0, 0.0
    trains = [(2.0, 4.0, 0.0, 0.0)]
    zone = spec_zone(length, h, ramp, trains)
    sim = Sim(nx, ny, zone=zone)
    sim.seed(lambda x, y: y < h)
    sim.calibrate()
    n = 45
    dt32 = f32(DT)
    for _ in range(n):
        sim.step(dt32)          # the node's time is step_index times its float32 dt
    d_out, ch, stats = node_step("lock", sim, n, trains=trains, zone=ports(length, h, ramp))
    t = n * dt32
    ref = twin(sim, ch, zone, t)
    speed = float(np.abs(ref.vel[ref.typ == 1]).max())
    good, rec = compare(f"zone 24 m, 2 m at 4 s, t {t + dt32:.2f} s", d_out, ref)
    print(f"      the water's fastest particle {speed:.2f} m/s; iterations {stats['iterations_used']:.0f}", flush=True)
    return good, dict(rec, t=round(t + dt32, 4), max_speed=round(speed, 4), iterations=stats["iterations_used"])


def compare(name, d_out, ref):
    """The node's grid, positions and velocities against the twin's, the level-set
    lockstep's tolerances (no solid here, so its penetration counts are left out)."""
    nx, ny = ref.nx, ref.ny
    pos = load(d_out, "out_positions")[:, :2]
    vel = load(d_out, "out_velocities")[:, :2]
    du = max(float(np.abs(load(d_out, "out_u").reshape(nx + 1, ny) - ref.u).max()),
             float(np.abs(load(d_out, "out_v").reshape(nx, ny + 1) - ref.v).max()))
    dp = float(np.abs(pos - ref.pos).max())
    dv = float(np.abs(vel - ref.vel).max())
    good = dp < 2e-3 and dv < 5e-2 and du < 1e-3
    print(f"    {name:34s} grid {du:.1e}, positions {dp:.1e}, velocities {dv:.1e}  {'OK' if good else 'FAIL'}")
    return good, {"grid": du, "positions": dp, "velocities": dv}


# ------------------------------------------------------------------ loops


def check_train():
    print("=== train: the wave train through the loop, 600 x 40 x 2, 50 s", flush=True)
    nx, ny, nz, h = 600, 40, 2, 20.0
    length, ramp, H, T = 24.0, 8.0, 2.0, 8.0
    dt, substeps, every = 0.05, 1000, 111
    t0 = time.time()
    d, typ, diags = run_bifrost3("zone_train", nx, ny, nz, lambda x, y, z: y < h, substeps, dt, fused=True,
                                 dump_every=every, extra_ports=dict(ports(length, h, ramp), gravity=-G),
                                 waves=waves_array([(H, T, 0.0, 0.0)]))
    wall = time.time() - t0
    k = Sim.wavenumber(T, h, G)
    L_lin = 2 * np.pi / k
    lo, hi = int(length + 100), int(length + 400)
    steps = list(range(0, substeps, every))
    levels = []
    for kk in steps:
        pos = np.load(os.path.join(d, f"frame.{kk:04d}.npy"))
        top_s = surface_top(pos[typ == 1], nx) - h
        levels.append(round(float(np.median(top_s[lo:hi])), 3))
    eta = top_s[lo:hi]
    xw = np.arange(lo, hi) + 0.5
    detr = eta - np.polyval(np.polyfit(xw, eta, 1), xw)
    L, Hm, n = wavelength_and_height(detr)
    finite = bool(np.isfinite(pos).all())
    ok = abs(L / L_lin - 1) <= 0.10 and abs(Hm / H - 1) <= 0.30 and n >= 1 and finite
    print(f"    at {(steps[-1] + 1) * dt:.0f} s over {lo}..{hi} m: L {L:.1f} m ({L / L_lin:.3f} of {L_lin:.1f}), "
          f"H {Hm:.2f} m ({Hm / H:.3f} of the train's {H:.1f}) over {n} waves; the window's level every "
          f"{every * dt:.2f} s {levels}; {len(typ):,} particles, iterations {diags[:, 0].mean():.0f}, "
          f"peak speed {diags[:, 3].max():.1f} m/s, {wall:.0f} s  " + ("OK" if ok else "FAIL"), flush=True)
    return ok, {"grid": [nx, ny, nz], "H": H, "T": T, "dt": dt, "L": round(L, 2), "L_ratio": round(L / L_lin, 4),
                "H_measured": round(Hm, 3), "H_ratio": round(Hm / H, 4), "waves": n, "level_drift_m": levels,
                "peak_speed": round(float(diags[:, 3].max()), 2), "wall_s": round(wall, 1)}


def surface_map(pos, nx, nz):
    """The highest liquid particle of each (x, z) column; -inf where none."""
    ci = np.clip(pos[:, 0].astype(np.int64), 0, nx - 1)
    ck = np.clip(pos[:, 2].astype(np.int64), 0, nz - 1)
    top = np.full((nx, nz), -np.inf)
    np.maximum.at(top, (ci, ck), pos[:, 1])
    return top


def crest_angle(top, x0, x1, z_rows, dz):
    """The heading of the crest lines from the shift between surface profiles dz
    apart: a train heading theta toward +z has its crests x cos + z sin = const, so the
    profile at z + dz is the profile at z moved by -dz tan(theta). The shift that best
    correlates the detrended profiles, averaged over the rows, in degrees."""
    shifts = []
    for z in z_rows:
        a, b = top[x0:x1, z], top[x0:x1, z + dz]
        if not (np.isfinite(a).all() and np.isfinite(b).all()):
            continue
        a, b = a - a.mean(), b - b.mean()
        best, best_c = 0, -np.inf
        for s in range(-int(0.8 * dz), int(0.8 * dz) + 1):
            if s >= 0:
                c = float(np.dot(a[s:], b[:len(b) - s]))
            else:
                c = float(np.dot(a[:s], b[-s:]))
            c /= max(1, len(a) - abs(s))
            if c > best_c:
                best, best_c = s, c
        shifts.append(best)
    if not shifts:
        return float("nan"), shifts
    # the best s has b(x) = a(x + s): the crest sits s cells toward smaller x at z + dz,
    # and s = dz tan(theta)
    return float(np.degrees(np.arctan2(np.mean(shifts), dz))), shifts


def check_oblique():
    print("=== oblique: a train heading 20 degrees, 160 x 20 x 96, production settings, 15 s", flush=True)
    nx, ny, nz, h = 160, 20, 96, 10.0
    length, ramp, H, T, angle = 16.0, 5.0, 1.0, 5.0, 20.0
    dt, substeps = 0.05, 300
    t0 = time.time()
    d, typ, diags = run_bifrost3("zone_oblique", nx, ny, nz, lambda x, y, z: y < h, substeps, dt, fused=True,
                                 dump_every=substeps - 1,
                                 # the production settings with their calibration (rho0_face sets the
                                 # escape threshold; at the graph's default of 0 every particle escapes)
                                 extra_ports=dict(calibration(), **ports(length, h, ramp), gravity=-G),
                                 waves=waves_array([(H, T, angle, 0.0)]), dumps=("ph",))
    wall = time.time() - t0
    k = substeps - 1
    pos = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
    ph = np.load(os.path.join(d, f"ph.{k:04d}.npy"))
    liq = pos[ph >= 0.5]
    top = surface_map(liq, nx, nz)
    kw = Sim.wavenumber(T, h, G)
    c = 2 * np.pi / T / kw          # the phase speed, cells a second
    x0, x1 = int(length + 10), int(min(nx - 10, length + c * (substeps * dt - 0.5 * ramp) * np.cos(np.radians(angle)) - 6))
    rows = list(range(16, nz - 24 - 16, 8))
    deg, shifts = crest_angle(top, x0, x1, rows, 24)
    eta = top[x0:x1, 8:nz - 8] - h
    finite = bool(np.isfinite(pos).all()) and bool(np.isfinite(eta).all())
    hmax = float(eta.max()) if finite else float("inf")
    ok = finite and abs(deg - angle) <= 5.0 and hmax <= 3.0 * H
    print(f"    crests over x {x0}..{x1}: heading {deg:.1f} degrees (the train's {angle:.0f}) from the shifts {shifts} "
          f"over 24 m; the surface at most {hmax:.2f} m above still water (the train {H:.0f} m); "
          f"{len(typ):,} particles seeded, {len(pos):,} at the end; iterations {diags[:, 0].mean():.0f}, "
          f"peak speed {diags[:, 3].max():.1f} m/s, {wall:.0f} s  " + ("OK" if ok else "FAIL"), flush=True)
    return ok, {"grid": [nx, ny, nz], "H": H, "T": T, "angle": angle, "heading_measured": round(deg, 2),
                "shifts": shifts, "eta_max": round(hmax, 3), "particles_end": int(len(pos)),
                "peak_speed": round(float(diags[:, 3].max()), 2), "wall_s": round(wall, 1)}


def main(argv):
    which = argv[0] if argv else "all"
    os.makedirs(RESULTS, exist_ok=True)
    ok, rec = True, {}
    for name, fn in (("off", check_off), ("lockstep", check_lockstep), ("train", check_train),
                     ("oblique", check_oblique)):
        if which in ("all", name):
            good, r = fn()
            ok &= good
            rec[name] = r
    with open(os.path.join(RESULTS, "receipt.json" if which == "all" else f"receipt_{which}.json"), "w") as f:
        json.dump(rec, f, indent=1, default=float)
    print("PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
