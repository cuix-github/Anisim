"""Gate for the droplet radii of the 3D step node (PF-FLIP Section 7, the node's
droplet_model), against the 2D reference's Sim(droplet=...) and the laws the model is
made of. Checks:

  off       the model off is the plain path: one step on a tank with liquid particles
            planted in the air (escaped droplets), in lockstep with the reference's
            plain twin; and bit for bit the same whether droplet_model is 0 or unset
  lockstep  the model on, nz = 1: the same tank with a wind over the water, half the
            planted droplets carrying a radius already (20 um to 2 mm) and half new
            (drawn by the hash), a few below the deletion threshold; grid, positions,
            velocities, radii and the deleted set against the reference's twin
  anchors   on the reference alone: the d^2 law exactly, the hash normal's and the
            log-normal radii's statistics, Stokes' relaxation of a 20 um drop in a
            light wind against the exponential, a 1 mm drop's terminal speed against
            the drag law's own closed form
  loop      the 3D loop (sim_3d_fused) with the model on: a splash planted over a
            tank, every rad dump positive where and only where escaped liquid is,
            droplets leaving by deletion, the final state carrying the radii and a
            resume reading them; the model off leaving every radius 0 and the count
            unchanged

Usage: python validate_droplet_node.py [off|lockstep|anchors|loop|all]
Receipt: results/3d/droplet/receipt.json (receipt_<check>.json for one check).
Exit nonzero on failure.
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "reference2d"))
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
from pfflip2d import Sim  # noqa: E402
from validate_step import channels  # noqa: E402
from validate_level_set_node import DT, flat_inputs, load, save  # noqa: E402
from validate_3d import run_node  # noqa: E402
from loop3d import load_final_state, run_bifrost3, seed3d  # noqa: E402

OUT = os.path.join(HERE, "..", "out", "droplet_node")
RESULTS = os.path.join(HERE, "..", "results", "3d", "droplet")
MPC = 0.25                       # metres per cell in the 2D checks: the breaker's cells
# the model's constants as the node's float32 ports hold them, so the reference
# computes with the same numbers
D = {k: float(np.float32(v))
     for k, v in dict(Sim.DROPLET_DEFAULTS, metres_per_cell=MPC).items()}


def pf(v):
    """A float32 port value as bifcmd reads it back exactly (positional, no exponent)."""
    return np.format_float_positional(np.float32(v), unique=True, trim="-")


PORTS = {"droplet_model": 1, "droplet_r_median": pf(D["r_median"]),
         "droplet_r_sigma": pf(D["r_sigma"]), "droplet_evap_k": pf(D["evap_k"]),
         "droplet_r_min": pf(D["r_min"]), "metres_per_cell": pf(D["metres_per_cell"]),
         "air_nu": pf(D["air_nu"])}
ESCAPE = {"escape": 1, "esc_phi": pf(0.3), "drag_droplet": pf(1.0), "drag_bubble": pf(8.0),
          "buoyancy": pf(2.0)}


# ------------------------------------------------------------------ the 2D state


def tank(nx=120, ny=40, depth=20.0, wind=0.0, seed=1):
    """A still tank of the reference, escape on, calibrated, a wind in its air."""
    sim = Sim(nx, ny, escape=True, droplet=D)
    sim.seed(lambda x, y: y < depth, seed=seed)
    sim.calibrate()
    if wind:
        sim.vel[(sim.typ == 0) & (sim.pos[:, 1] > depth + 2.0), 0] = wind
    return sim


def plant(sim, n, lo, hi, seed=2):
    """n liquid particles from the bed lifted into the air between heights lo and hi:
    droplets to the escape test. Returns their indices."""
    rng = np.random.default_rng(seed)
    idx = np.flatnonzero((sim.typ == 1) & (sim.pos[:, 1] < 4.0))[:n]
    sim.pos[idx, 1] = rng.uniform(lo, hi, len(idx))
    sim.pos[idx, 0] = rng.uniform(10.0, sim.nx - 10.0, len(idx))
    return idx


def twin(state, ch, radius_in, step_index, model):
    """validate_level_set_node.twin_step on the node's inputs, escape on, with the
    droplet model or the plain constant drag. Returns the reference after the step,
    the escape flags of the surviving particles, and the deleted mask."""
    ref = Sim(state.nx, state.ny, escape=True, droplet=D if model else None)
    arr = flat_inputs(state, ch)
    ref.pos = arr["positions"][:, :2].astype(float)
    ref.vel = arr["velocities"][:, :2].astype(float)
    ref.typ = state.typ.copy()
    ref.radius = radius_in.astype(float)
    ref.nstep = step_index
    ref.rho0_face = state.rho0_face
    um, uo = ch["u_mass"].astype(float), ch["u_mom"].astype(float)
    vm, vo = ch["v_mass"].astype(float), ch["v_mom"].astype(float)
    ref.u = np.where(um > 0, uo / np.maximum(um, 1e-30), 0.0)
    ref.v = np.where(vm > 0, vo / np.maximum(vm, 1e-30), 0.0)
    ref._enforce_walls()
    u_star, v_star = ref.u.copy(), ref.v.copy()
    ref.v += ref.g * DT
    ref._enforce_walls()
    phi_u = np.clip(ch["u_phase"].astype(float), 0, 1)
    phi_v = np.clip(ch["v_phase"].astype(float), 0, 1)
    ref.mu, ref.mv = np.empty_like(ref.u), np.empty_like(ref.v)
    ref.phase = lambda raw: phi_u if raw is ref.mu else phi_v
    ref.project(DT)
    ones_air = (ref.typ == 0).astype(float)
    fu, fv = ref._face_frames()
    au, lu = ref._splat(ref.u.shape, fu, [ones_air, 1 - ones_air])
    av, lv = ref._splat(ref.v.shape, fv, [ones_air, 1 - ones_air])
    thresh = (1.0 - ref.esc_phi) * ref.rho0_face
    fa = 0.5 * (ref._sample_faces(au, np.array([0.0, 0.5]), ref.pos) +
                ref._sample_faces(av, np.array([0.5, 0.0]), ref.pos))
    fl = 0.5 * (ref._sample_faces(lu, np.array([0.0, 0.5]), ref.pos) +
                ref._sample_faces(lv, np.array([0.5, 0.0]), ref.pos))
    drop = (ref.typ == 1) & (fa > thresh)
    bub = (ref.typ == 0) & (fl > thresh)
    esc = drop | bub
    ref.g2p(u_star, v_star, skip=esc)
    gone = ref._droplet_radii(drop, DT) if model else np.zeros(len(ref.pos), bool)
    if esc.any():
        vg = ref.sample_velocity(ref.pos)
        ref.vel[drop, 1] += ref.g * DT
        if model:
            if drop.any():
                ref._droplet_drag(drop, vg, DT)
        else:
            ref.vel[drop] += ref.drag_droplet * DT * (vg[drop] - ref.vel[drop])
        ref.vel[bub, 1] += ref.buoyancy * abs(ref.g) * DT
        ref.vel[bub] += ref.drag_bubble * DT * (vg[bub] - ref.vel[bub])
    ref.advect(DT, ballistic=esc)
    if gone.any():
        ref._delete(~gone)
    return ref, esc[~gone], gone


def node_step(tag, sim, radius_in, **ports):
    """One step_3d call at nz = 1 on the reference state, radius_in as its state."""
    ch = channels(sim)
    arr = flat_inputs(sim, ch)
    arr["radius_in"] = radius_in.astype(np.float32)
    d_in, d_out = os.path.join(OUT, tag + "_in"), os.path.join(OUT, tag + "_out")
    for k, v in arr.items():
        save(d_in, k, v)
    stats = run_node(d_in, d_out, extra_arrays=("radius_in",), extra_outputs=("radius_out",),
                     nx=sim.nx, ny=sim.ny, nz=1, dt=DT, preconditioner=1,
                     rho0_face=pf(sim.rho0_face), **ESCAPE, **ports)
    return d_out, ch, stats


def compare(name, d_out, ref, esc):
    """The node's outputs against the twin's: grid, positions, velocities, escape
    flags, the particle count, and the radii (droplets within 1e-5, the rest 0)."""
    nx, ny = ref.nx, ref.ny
    pos = load(d_out, "out_positions")[:, :2]
    if len(pos) != len(ref.pos):
        print(f"    {name:28s} particle count {len(pos)} vs the reference's {len(ref.pos)}  FAIL")
        return False, {"count_node": len(pos), "count_ref": len(ref.pos)}
    vel = load(d_out, "out_velocities")[:, :2]
    rad = load(d_out, "radius_out")
    du = max(float(np.abs(load(d_out, "out_u").reshape(nx + 1, ny) - ref.u).max()),
             float(np.abs(load(d_out, "out_v").reshape(nx, ny + 1) - ref.v).max()))
    dp = float(np.abs(pos - ref.pos).max())
    dv = float(np.abs(vel - ref.vel).max())
    mism = int(((load(d_out, "out_escaped") > 0.5) != esc).sum())
    drops = ref.radius > 0
    dr = float((np.abs(rad[drops] - ref.radius[drops]) / ref.radius[drops]).max()) if drops.any() else 0.0
    zeros = bool((rad[~drops] == 0).all()) and bool((rad[drops] > 0).all())
    good = mism <= 2 and dp < 2e-3 and dv < 5e-2 and du < 1e-3 and dr < 1e-5 and zeros
    print(f"    {name:28s} grid {du:.1e}, positions {dp:.1e}, velocities {dv:.1e}, escape"
          f" mismatches {mism}, {int(drops.sum())} droplets' radii within {dr:.1e},"
          f" the rest {'0' if zeros else 'NOT 0'}, {len(pos)} particles  {'OK' if good else 'FAIL'}")
    return good, {"grid": du, "positions": dp, "velocities": dv, "escape_mismatches": mism,
                  "droplets": int(drops.sum()), "radius_rel": dr, "count": len(pos)}


# ------------------------------------------------------------------ checks


def check_off():
    print("=== off: droplet_model 0 is the plain path, nz = 1", flush=True)
    sim = tank()
    plant(sim, 60, 26.0, 36.0)
    zero = np.zeros(len(sim.pos), np.float32)
    d0, ch, _ = node_step("off_model0", sim, zero, droplet_model=0)
    ref, esc, _ = twin(sim, ch, zero, 0, model=False)
    good, rec = compare("plain twin", d0, ref, esc)
    rec["escaped_planted"] = int(esc.sum())
    good &= esc.sum() >= 50
    # unset is 0: every output bit for bit
    d1, _, _ = node_step("off_unset", sim, zero)
    same = all(np.array_equal(load(d0, n), load(d1, n)) for n in
               ("out_u", "out_v", "out_w", "pressure", "out_positions", "out_velocities",
                "out_escaped", "radius_out"))
    allzero = bool((load(d0, "radius_out") == 0).all())
    print(f"    droplet_model unset against 0: {'bit-identical' if same else 'DIFFER'};"
          f" radius_out all 0: {allzero}  {'OK' if same and allzero else 'FAIL'}")
    rec.update(unset_identical=same, radius_all_zero=allzero)
    return good and same and allzero, rec


def check_lockstep():
    print("=== lockstep: the model on against the reference's twin, nz = 1", flush=True)
    sim = tank(wind=3.0)
    idx = plant(sim, 60, 26.0, 36.0)
    radius = np.zeros(len(sim.pos))
    # 30 droplets with a radius already, 20 um to 2 mm; 25 new (0, drawn by the hash);
    # 5 just above the deletion threshold, gone after one step of evaporation
    radius[idx[:30]] = np.logspace(np.log10(2e-5), np.log10(2e-3), 30)
    radius[idx[55:]] = 1.2e-6
    n = 7
    d_out, ch, stats = node_step("lock", sim, radius, step_index=n, **PORTS)
    ref, esc, gone = twin(sim, ch, radius, n, model=True)
    good, rec = compare("droplet twin", d_out, ref, esc)
    new = int((ref.radius > 0).sum()) - 30
    print(f"    planted {len(idx)}: {int(gone.sum())} deleted (expected 5), {new} drew a new"
          f" radius, radii now {ref.radius[ref.radius > 0].min() * 1e6:.1f} to"
          f" {ref.radius.max() * 1e6:.0f} um; iterations {stats['iterations_used']:.0f}", flush=True)
    good &= int(gone.sum()) == 5 and new >= 20
    rec.update(deleted=int(gone.sum()), new_radii=new, iterations=stats["iterations_used"])
    return good, rec


def check_anchors():
    print("=== anchors: the laws the model is made of, on the reference", flush=True)
    ok, rec = True, {}
    one = np.array([True])
    # the d^2 law: r^2 falls by evap_k / 4 per second, exactly, step after step (a
    # 100 um drop at ten times dry air's rate loses half its r^2 in these 2 s)
    s = Sim(4, 4, escape=True, droplet=dict(D, evap_k=1e-8, metres_per_cell=1.0))
    s.radius = np.array([1e-4])
    r2 = []
    for _ in range(200):
        s._droplet_radii(one, 0.01)
        r2.append(s.radius[0] ** 2)
    expect = 1e-8 - 0.25 * 1e-8 * 0.01 * np.arange(1, 201)
    err = float(np.abs(np.array(r2) - expect).max() / 1e-8)
    good = err < 1e-9
    print(f"    d^2 law over 200 steps: r^2 off its line by {err:.1e} of r0^2  {'OK' if good else 'FAIL'}")
    ok &= good
    rec["d2_law_rel"] = err
    # the hash normal and the radii it draws
    z = Sim.droplet_normal(3, np.arange(400000))
    r = D["r_median"] * np.exp(D["r_sigma"] * z)
    st = {"mean": float(z.mean()), "std": float(z.std()), "median_um": float(np.median(r) * 1e6),
          "ln_std": float(np.log(r).std())}
    good = (abs(st["mean"]) < 0.01 and abs(st["std"] - 1.0) < 0.01
            and abs(st["median_um"] / (D["r_median"] * 1e6) - 1.0) < 0.02
            and abs(st["ln_std"] / D["r_sigma"] - 1.0) < 0.02)
    print(f"    hash normal, 400k draws: mean {st['mean']:+.4f}, std {st['std']:.4f}; radii median"
          f" {st['median_um']:.1f} um (set {D['r_median'] * 1e6:.0f}), spread of ln r"
          f" {st['ln_std']:.4f} (set {D['r_sigma']})  {'OK' if good else 'FAIL'}")
    ok &= good
    rec["normal"] = st
    # Stokes: a 20 um drop at rest in a wind of 1 mm/s (Re 0.003), v/U = 1 - exp(-t/tau)
    s = Sim(4, 4, escape=True, droplet=dict(D, metres_per_cell=1.0))
    r0 = 2e-5
    tau = r0 ** 2 / (4.5 * (s.rho_g / s.rho_l) * D["air_nu"])
    s.radius, s.vel = np.array([r0]), np.zeros((1, 2))
    vg, dt = np.array([[1e-3, 0.0]]), tau / 20
    v = []
    for _ in range(60):
        s._droplet_drag(one, vg, dt)
        v.append(s.vel[0, 0] / 1e-3)
    t = dt * np.arange(1, 61)
    err = float(np.abs(np.array(v) - (1.0 - np.exp(-t / tau))).max())
    good = err < 0.03
    print(f"    Stokes relaxation of a 20 um drop (tau {tau * 1e3:.2f} ms, dt tau/20): within"
          f" {err:.3f} of 1 - exp(-t/tau) over 3 tau  {'OK' if good else 'FAIL'}")
    ok &= good
    rec["stokes"] = {"tau_ms": tau * 1e3, "max_err": err}
    # the terminal speed of a 1 mm drop under gravity in still air: the law's own
    # closed form at C_d 0.44, sqrt(8 g r / (3 C_d rho_a/rho_l)), Re above 1000 there
    s.radius, s.vel = np.array([1e-3]), np.zeros((1, 2))
    for _ in range(4000):
        s.vel[0, 1] += -9.8 * 1e-3
        s._droplet_drag(one, np.zeros((1, 2)), 1e-3)
    vt = -s.vel[0, 1]
    vt_closed = float(np.sqrt(8.0 * 9.8 * 1e-3 / (3.0 * 0.44 * (s.rho_g / s.rho_l))))
    re = 2.0 * 1e-3 * vt / D["air_nu"]
    good = abs(vt / vt_closed - 1.0) < 0.01 and re > 1000.0
    print(f"    terminal speed of a 1 mm drop: {vt:.2f} m/s against the closed form"
          f" {vt_closed:.2f} m/s, Re {re:.0f}  {'OK' if good else 'FAIL'}")
    ok &= good
    rec["terminal"] = {"v": vt, "closed_form": vt_closed, "re": re}
    return ok, rec


def check_loop():
    print("=== loop: the 3D loop with the model on and off, 40 x 40 x 8", flush=True)
    nx, ny, nz, dt = 40, 40, 8, 0.02
    cal = Sim(40, 40, escape=True)
    cal.seed(lambda x, y: np.ones_like(x, dtype=bool))
    cal.calibrate()

    def seed_fn(nx_, ny_, nz_, mask):
        # the tank's water, 400 liquid particles of its bed lifted into the air
        pos, typ = seed3d(nx_, ny_, nz_, mask)
        rng = np.random.default_rng(5)
        idx = np.flatnonzero((typ == 1) & (pos[:, 1] < 4.0))[:400]
        pos[idx, 0] = rng.uniform(4.0, nx_ - 4.0, len(idx))
        pos[idx, 1] = rng.uniform(26.0, 36.0, len(idx))
        pos[idx, 2] = rng.uniform(1.0, nz_ - 1.0, len(idx))
        return pos, typ

    base = dict(ESCAPE, rho0_face=pf(cal.rho0_face))
    # cells of a decimetre, evaporation a hundred times dry air so the droplets go
    # within the run (a 100 um drop then lives 0.4 s)
    on = {**base, **PORTS, "droplet_evap_k": pf(1e-7), "metres_per_cell": pf(0.1)}
    water = lambda x, y, z: y < 12  # noqa: E731
    d_on, _, _ = run_bifrost3("droplet_on", nx, ny, nz, water, 40, dt, extra_ports=on,
                              seed_fn=seed_fn, fused=True, dump_every=1, final_state=True,
                              out_dir=os.path.join(OUT, "loop_on"))
    ok, rec = True, {}
    counts, drops, placed = [], [], True
    for k in range(40):
        rad = np.load(os.path.join(d_on, f"rad.{k:04d}.npy"))
        esc = np.load(os.path.join(d_on, f"esc.{k:04d}.npy")) > 0.5
        liq = np.load(os.path.join(d_on, f"ph.{k:04d}.npy")) > 0.5
        placed &= bool(np.array_equal(rad > 0, esc & liq))
        counts.append(len(rad))
        drops.append(int((rad > 0).sum()))
    fin = np.load(os.path.join(d_on, "final_radius.npy"))
    last = np.load(os.path.join(d_on, "rad.0039.npy"))
    final_ok = np.array_equal(fin, last)
    shrinking = all(a >= b for a, b in zip(counts, counts[1:])) and counts[-1] < counts[0]
    good = placed and shrinking and final_ok and drops[0] >= 300
    print(f"    radii positive where and only where escaped liquid, every dump: {placed}")
    print(f"    particles {counts[0]} -> {counts[-1]} ({counts[0] - counts[-1]} droplets deleted),"
          f" droplets {drops[0]} -> {drops[-1]}, never increasing: {shrinking}")
    print(f"    final_radius.npy equals the last dump: {final_ok}  {'OK' if good else 'FAIL'}", flush=True)
    ok &= good
    rec.update(count_first=counts[0], count_last=counts[-1], droplets_first=drops[0],
               droplets_last=drops[-1], placed=placed, final_equals_last=final_ok)
    # a resume reads the radii back and carries on
    init = load_final_state(d_on)
    d_re, _, _ = run_bifrost3("droplet_resume", nx, ny, nz, water, 5, dt, extra_ports=on,
                              fused=True, dump_every=1, init_arrays=init, step_offset=40,
                              out_dir=os.path.join(OUT, "loop_resume"))
    rad0 = np.load(os.path.join(d_re, "rad.0000.npy"))
    resumed = int((rad0 > 0).sum())
    good = resumed <= drops[-1] and len(rad0) <= counts[-1] and "radius" in init
    print(f"    resume from the final state: {resumed} droplets on its first step"
          f" (from {drops[-1]})  {'OK' if good else 'FAIL'}", flush=True)
    ok &= good
    rec.update(resumed_droplets=resumed)
    # the model off: no radius, no deletion
    d_off, _, _ = run_bifrost3("droplet_off", nx, ny, nz, water, 10, dt, extra_ports=base,
                               seed_fn=seed_fn, fused=True, dump_every=1,
                               out_dir=os.path.join(OUT, "loop_off"))
    zero = all(bool((np.load(os.path.join(d_off, f"rad.{k:04d}.npy")) == 0).all()) for k in range(10))
    n_off = [len(np.load(os.path.join(d_off, f"frame.{k:04d}.npy"))) for k in range(10)]
    seeded = len(np.load(os.path.join(d_on, "positions.npy")))
    good = zero and len(set(n_off)) == 1 and n_off[0] == seeded
    print(f"    model off: radii all 0 {zero}, particle count {n_off[0]} unchanged over 10"
          f" steps {len(set(n_off)) == 1}  {'OK' if good else 'FAIL'}", flush=True)
    ok &= good
    rec.update(off_all_zero=zero, off_count=n_off[0])
    return ok, rec


def main(argv):
    which = argv[0] if argv else "all"
    os.makedirs(RESULTS, exist_ok=True)
    os.makedirs(OUT, exist_ok=True)
    ok, rec = True, {}
    for name, fn in (("off", check_off), ("lockstep", check_lockstep),
                     ("anchors", check_anchors), ("loop", check_loop)):
        if which in ("all", name):
            good, r = fn()
            ok &= good
            rec[name] = r
    with open(os.path.join(RESULTS, "receipt.json" if which == "all" else f"receipt_{which}.json"),
              "w") as f:
        json.dump(rec, f, indent=1, default=float)
    print("DROPLET NODE:", "PASS" if ok else "FAIL", flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
