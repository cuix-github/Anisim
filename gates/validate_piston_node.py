"""The piston wavemaker in step_3d against its spec, the 2D reference's plate
(reference2d/validate_piston.py): a plate on the floor, the whole depth in z, sliding
as x0 + S/2 (1 - cos 2 pi t / T), its speed on the faces it closes, in the divergence
of the cut faces it covers and in the push-out. Four checks, one-metre cells:

  rest      the node with a stroke of 0 is the node with the same plate as a static
            cut-cell solid, every output bit for bit (one step of a dam break at nz = 1).
  lockstep  one step from one state, the node at nz = 1 against the reference's twin
            step at the same instant: the plate alone in a channel, mid-stroke at its
            fastest; and the plate with a static beach folded in. Grid velocities,
            positions and velocities within the tolerances of the level-set lockstep.
  train     the dispersion check through the whole loop (the fused graph, nz = 2):
            600 x 40 x 2, 20 m deep, S = 2 m, T = 8 s, 50 s; the surface at 50 s from
            100 to 400 m in front of the plate against linear theory (wavelength within
            10 percent, height within 30 percent of the transfer function), the still
            level's drift recorded; no liquid inside or behind the plate in any dump.
  fast      the hard stroke in 3D under the production settings (ST, escaped
            particles, adaptive air): 200 x 60 x 4, S = 4 m, T = 4 s, the plate's top at
            50 m, 12 s; no liquid particle crosses the plate below its top between two
            dumps 0.1 s apart, the liquid count kept. The stroke is violent (waves of
            steepness 0.3 break at the paddle) and its spray, escaped droplets dragged by
            air at over 100 m/s, reaches the ceiling 40 m up and falls behind the plate:
            counted and reported, not a leak (the fast-air question in the open questions).

Run from the repo root:  .venv\\Scripts\\python.exe gates\\validate_piston_node.py
[rest|lockstep|train|fast]; about five minutes at eight threads. Exits nonzero on failure.
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
from solitary_wave import beach_phi  # noqa: E402
from validate_piston import transfer, wavelength_and_height, wavenumber  # noqa: E402
from validate_step import channels  # noqa: E402
from validate_level_set_node import DT, extrude, f32, flat_inputs, load, lockstep, node  # noqa: E402
from loop3d import run_bifrost3, seed3d  # noqa: E402

RESULTS = os.path.join(HERE, "..", "results", "3d", "piston")
X0, WIDTH = 8.0, 2.0


def ports(S, T, top):
    return {"piston_stroke": S, "piston_period": T, "piston_x0": X0, "piston_width": WIDTH,
            "piston_top": top}


def plate_x(S, T, t):
    return X0 + 0.5 * S * (1.0 - np.cos(2.0 * np.pi * t / T))


def seed_clear_of(top):
    """seed3d with nothing inside the plate at rest: the reference's seed drops what
    lies inside a solid, and an air particle left inside is pushed out of the plate's
    front face into the water on the first step, where a column of such bubbles rises
    along the plate and drags water up it (to 30 m over 20 m of water in 3 s, the first
    run of this gate)."""
    def fn(nx, ny, nz, mask):
        pos, typ = seed3d(nx, ny, nz, mask)
        keep = ~((pos[:, 0] > X0) & (pos[:, 0] < X0 + WIDTH) & (pos[:, 1] < top))
        return pos[keep], typ[keep]
    return fn


# ------------------------------------------------------------------ rest


def check_rest():
    print("=== rest: a stroke of 0 is the same plate as a static solid, bit for bit", flush=True)
    nx, ny, top = 160, 80, 56.5
    sim = Sim(nx, ny)
    sim.seed(lambda x, y: (x > X0 + WIDTH) & (x < 40) & (y < 56))
    sim.calibrate()
    for _ in range(20):
        sim.step(DT)
    arr = flat_inputs(sim, channels(sim))
    plate = Sim.plate_phi(nx, ny, X0, WIDTH, top)
    d_p, s_p = node("rest_piston", nx, ny, 1, arr, step_index=20, **ports(0.0, 1.0, top))
    d_s, s_s = node("rest_static", nx, ny, 1, dict(arr, solid_phi=extrude(plate, 1).astype(np.float32).ravel()),
                    solid_mode=2)
    worst = 0.0
    for name in ("out_u", "out_v", "out_w", "pressure", "out_positions", "out_velocities"):
        a, b = load(d_p, name), load(d_s, name)
        worst = max(worst, float(np.abs(a.astype(np.float64) - b).max()) if a.shape == b.shape else np.inf)
    same = worst == 0.0 and s_p["iterations_used"] == s_s["iterations_used"]
    print(f"    one step, {len(arr['positions']):,} particles: max |diff| {worst:.1e} over the grid, pressure, "
          f"positions and velocities; iterations {s_p['iterations_used']:.0f} and {s_s['iterations_used']:.0f}  "
          + ("OK" if same else "FAIL"), flush=True)
    return same, {"max_diff": worst}


# ------------------------------------------------------------------ lockstep


def piston_twin(state, ch, kw, t):
    """validate_level_set_node.twin_step with the plate placed at t (no escape)."""
    ref = Sim(state.nx, state.ny, **kw)
    ref.solid_support = False          # the node's phase is the channel's (no share)
    ref.t = t
    ref._move_piston(t, float(np.float32(DT)))
    arr = flat_inputs(state, ch)
    ref.pos = arr["positions"][:, :2].astype(float)
    ref.vel = arr["velocities"][:, :2].astype(float)
    ref.typ = state.typ.copy()
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
    ref.g2p(u_star, v_star)
    ref.advect(DT)
    ref._push_out()
    return ref


def check_lockstep():
    print("=== lockstep: one step against the reference's twin, nz = 1", flush=True)
    ok, rec = True, {}
    nx, ny, h, S, T, top = 120, 40, 20.0, 2.0, 4.0, 28.5
    for name, static in (("channel", None), ("channel and beach", f32(beach_phi(nx, ny, 60.0, 1.0 / 3.0)))):
        kw = {"piston": (S, T, X0, WIDTH, top)}
        if static is not None:
            kw.update(solid_phi=static, solid_cut=True)
        sim = Sim(nx, ny, **kw)
        sim.seed(lambda x, y: (x > X0 + WIDTH) & (y < h))
        sim.calibrate()
        # t = 0.9 s, near a quarter period (the plate near its fastest), its faces off
        # the grid nodes: at exactly 1 s its right face sits on a node column, and the
        # float32 time the node computes leaves it 35e-9 short, a face open in one and
        # closed in the other
        n = 45
        for _ in range(n):
            sim.step(DT)
        # plant particles inside the plate (and inside the beach) so the push-out
        # with the plate's speed, and without it in the static solid, both act
        rng = np.random.default_rng(3)
        xl = float(plate_x(S, T, n * float(np.float32(DT))))
        near = np.flatnonzero((sim.pos[:, 0] > xl + WIDTH) & (sim.pos[:, 0] < xl + WIDTH + 1.0)
                              & (sim.pos[:, 1] < h - 2))[:40]
        sim.pos[near, 0] = xl + rng.uniform(0.2, WIDTH - 0.2, len(near))
        if static is not None:
            deep = np.flatnonzero((sim.pos[:, 0] > 70.0) & (sim.pos[:, 0] < 80.0) & (sim.typ == 1))[:40]
            bed = (sim.pos[deep, 0] - 60.0) / 3.0
            sim.pos[deep, 1] = np.maximum(bed - rng.uniform(0.2, 1.0, len(deep)), 0.6)
        ch = channels(sim)
        arr = flat_inputs(sim, ch)
        extra = {}
        if static is not None:
            arr["solid_phi"] = extrude(static, 1).astype(np.float32).ravel()
            extra["solid_mode"] = 2
        d_out, stats = node("lock_" + name.replace(" ", "_"), nx, ny, 1, arr, step_index=n,
                            **extra, **ports(S, T, top))
        t = n * float(np.float32(DT))            # the node's time: step_index times its float32 dt
        ref = piston_twin(sim, ch, kw, t)
        p0 = arr["positions"][:, :2].astype(float)
        swept = int((ref._phi_at(p0, ref._plate_phi)[0] < 0.0).sum())
        in_static = int((ref._phi_at(p0, ref.static_phi)[0] < 0.0).sum()) if static is not None else 0
        good = lockstep(f"{name} (plate at {ref.piston_x:.2f}, {ref.solid_vel[0]:.2f} m/s)", d_out, ref,
                        None, extrude(ref.solid_phi, 1))
        print(f"      planted before the step: {swept} inside the plate, {in_static} inside the beach", flush=True)
        good &= swept >= 20 and (static is None or in_static >= 20)
        ok &= good
        rec[name] = {"plate_x": round(ref.piston_x, 4), "plate_speed": round(float(ref.solid_vel[0]), 4),
                     "in_plate": swept, "in_beach": in_static, "iterations": stats["iterations_used"]}
    return ok, rec


# ------------------------------------------------------------------ loops


def plate_audit(d, typ_of, S, T, top, dt, steps):
    """Liquid inside and behind the plate in every dump."""
    worst_in, worst_behind = 0, 0
    for k in steps:
        pos = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
        typ = typ_of(k)
        liq = pos[typ == 1]
        xl = plate_x(S, T, k * dt)             # where the step that wrote dump k had it
        inside = (liq[:, 0] > xl) & (liq[:, 0] < xl + WIDTH) & (liq[:, 1] < top)
        behind = (liq[:, 0] < xl) & (liq[:, 1] < top)
        worst_in, worst_behind = max(worst_in, int(inside.sum())), max(worst_behind, int(behind.sum()))
    return worst_in, worst_behind


def surface_top(pos, nx):
    col = np.clip(pos[:, 0].astype(np.int64), 0, nx - 1)
    top = np.full(nx, -np.inf)
    np.maximum.at(top, col, pos[:, 1])
    k = np.ones(5) / 5.0
    wet = np.isfinite(top)
    filled = np.where(wet, top, 0.0)
    sm = top.copy()
    sm[wet] = (np.convolve(filled, k, "same") / np.maximum(np.convolve(wet.astype(float), k, "same"), 1e-9))[wet]
    return sm


def check_train():
    print("=== train: the wave train through the loop, 600 x 40 x 2, 50 s", flush=True)
    nx, ny, nz, h, S, T, top = 600, 40, 2, 20.0, 2.0, 8.0, 28.5
    dt, substeps, every = 0.05, 1000, 111
    front = X0 + WIDTH
    t0 = time.time()
    d, typ, diags = run_bifrost3("piston_train", nx, ny, nz, lambda x, y, z: (x > front) & (y < h),
                                 substeps, dt, fused=True, dump_every=every, seed_fn=seed_clear_of(top),
                                 extra_ports=dict(ports(S, T, top), gravity=-9.8))
    wall = time.time() - t0
    steps = list(range(0, substeps, every))
    w_in, w_behind = plate_audit(d, lambda k: typ, S, T, top, dt, steps)
    k = wavenumber(T, h)
    L_lin, H_lin = 2 * np.pi / k, S * transfer(k, h)
    lo, hi = int(front + 100), int(front + 400)
    levels = []
    for kk in steps:
        pos = np.load(os.path.join(d, f"frame.{kk:04d}.npy"))
        top_s = surface_top(pos[typ == 1], nx) - h
        levels.append(round(float(np.median(top_s[lo:hi])), 3))
    eta = top_s[lo:hi]
    xw = np.arange(lo, hi) + 0.5
    detr = eta - np.polyval(np.polyfit(xw, eta, 1), xw)
    L, H, n = wavelength_and_height(detr)
    ok = abs(L / L_lin - 1) <= 0.10 and abs(H / H_lin - 1) <= 0.30 and n >= 2 and w_in == 0 and w_behind == 0
    print(f"    at {(steps[-1] + 1) * dt:.0f} s over {lo}..{hi} m: L {L:.1f} m ({L / L_lin:.3f} of {L_lin:.1f}), "
          f"H {H:.2f} m ({H / H_lin:.3f} of {H_lin:.2f}) over {n} waves; the window's level every "
          f"{every * dt:.2f} s {levels}; liquid inside the plate {w_in}, behind it {w_behind}; "
          f"{len(typ):,} particles, iterations {diags[:, 0].mean():.0f}, {wall:.0f} s  " + ("OK" if ok else "FAIL"),
          flush=True)
    return ok, {"grid": [nx, ny, nz], "S": S, "T": T, "dt": dt, "L": round(L, 2), "L_ratio": round(L / L_lin, 4),
                "H": round(H, 3), "H_ratio": round(H / H_lin, 4), "waves": n, "level_drift_m": levels,
                "inside": w_in, "behind": w_behind, "wall_s": round(wall, 1)}


def check_fast():
    print("=== fast: a 4 m stroke at 4 s in 3D, production settings", flush=True)
    nx, ny, nz, h, S, T, top = 200, 60, 4, 20.0, 4.0, 4.0, 50.0
    dt, substeps, every = 0.02, 600, 5
    front = X0 + WIDTH
    t0 = time.time()
    d, typ, diags = run_bifrost3("piston_fast", nx, ny, nz, lambda x, y, z: (x > front) & (y < h),
                                 substeps, dt, fused=True, dump_every=every, seed_fn=seed_clear_of(top),
                                 extra_ports=dict(ports(S, T, top), gravity=-9.8, st=1, escape=1, adapt=1),
                                 dumps=("esc", "scl", "ph"))
    wall = time.time() - t0
    # liquid never merges or splits and keeps its order, so the liquid subset of two
    # dumps is the same particles: a crossing is a particle in front of the plate and
    # below its top in one dump and behind it, still below the top, in the next
    n_liq = int((typ == 1).sum())
    prev, crossed, inside, over, high, lost = None, 0, 0, 0, 0.0, 0
    for k in range(0, substeps, every):
        pos = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
        liq = pos[np.load(os.path.join(d, f"ph.{k:04d}.npy")) > 0.5]
        lost = max(lost, abs(len(liq) - n_liq))
        xl = plate_x(S, T, k * dt)
        behind = (liq[:, 0] < xl) & (liq[:, 1] < top - 1.0)
        inside = max(inside, int(((liq[:, 0] > xl) & (liq[:, 0] < xl + WIDTH) & (liq[:, 1] < top)).sum()))
        high = max(high, float(liq[:, 1].max()))
        if prev is not None:
            assert len(prev[0]) == len(liq), (k, len(prev[0]), len(liq))
            was_front = (prev[0][:, 0] > prev[1] + WIDTH) & (prev[0][:, 1] < top - 1.0)
            crossed += int((behind & was_front).sum())
            over += int((behind & ~prev[2]).sum())
        prev = (liq, xl, behind)
    ok = crossed == 0 and lost == 0
    print(f"    12 s, a dump every {every * dt:.2f} s: liquid crossing the plate below its top {crossed}; arriving "
          f"behind it over the top {over}; most inside the plate in a dump {inside} (the ST time jitter: a "
          f"particle's position belongs to its own instant); liquid count {n_liq:,}, worst change {lost}; "
          f"highest liquid {high:.1f} m (ceiling {ny}); peak speed {diags[:, 3].max():.1f} cells/s; "
          f"iterations {diags[:, 0].mean():.0f}; {wall:.0f} s  " + ("OK" if ok else "FAIL"), flush=True)
    return ok, {"grid": [nx, ny, nz], "crossed": crossed, "over_the_top": over, "inside_max": inside,
                "liquid": n_liq, "highest_liquid": round(high, 2), "peak_speed": round(float(diags[:, 3].max()), 1),
                "wall_s": round(wall, 1)}


def main(argv):
    which = argv[0] if argv else "all"
    os.makedirs(RESULTS, exist_ok=True)
    ok, rec = True, {}
    for name, fn in (("rest", check_rest), ("lockstep", check_lockstep), ("train", check_train),
                     ("fast", check_fast)):
        if which in ("all", name):
            good, r = fn()
            ok &= good
            rec[name] = r
    with open(os.path.join(RESULTS, "receipt.json" if which == "all" else f"receipt_{which}.json"), "w") as f:
        json.dump(rec, f, indent=1, default=float)
    print("PISTON NODE:", "PASS" if ok else "FAIL", flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
