"""Gate for the extra spray emission and the droplets' own integrator of the 3D step
node (PF-FLIP Section 7 after Ihmsen et al. 2012; the node's emit_model and
spray_cfl), against the 2D reference's Sim(emit=..., spray_cfl=...) and the laws they
are made of. Checks:

  off       emission off and spray_cfl 0 are the plain droplet path: one step on a
            tank with droplets planted in the air, bit for bit the same whether the
            spray ports are 0 or unset, spray_out all 0
  lockstep  emission and the integrator on, nz = 1: a tank with a wind over the water
            and a sheared surface layer (steep gradients at the interface), droplets
            planted in the air and spray particles already flying; grid, positions,
            velocities, radii, escape and spray flags and the count against the
            reference's twin, the spray it emits included
  anchors   on the reference alone: a still tank emits nothing; the count emitted
            falls as the threshold rises; the integrator's position error on a
            20 um drop relaxing in a light wind under gravity falls fourfold when
            the sub-step halves (second order; its velocity is exact); the particles
            that are not spray are bit-identical with emission on and off, with and
            without ST and adaptivity (the spray is massless and counted nowhere);
            and 800 steps of the sheared tank with emission on stay finite, the spray
            alive bounded, the others' count unchanged
  loop      the 3D loop (sim_3d_fused, 40 x 40 x 8) resumed from a sheared state with
            emission on: spr dumps 1 exactly where the mass is 0, radii positive
            there past the emitting step, the others' count the run without
            emission's at every dump, spray emitted; the same resume with emission
            off leaves the others bit-identical (the fused transfer skips a
            massless particle); a resume reads the spray flags back

Usage: python validate_spray_node.py [off|lockstep|anchors|loop|all]
Receipt: results/3d/spray/receipt.json (receipt_<check>.json for one check).
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
from validate_droplet_node import D, ESCAPE, PORTS, pf, plant  # noqa: E402
from loop3d import load_final_state, run_bifrost3, seed3d  # noqa: E402

OUT = os.path.join(HERE, "..", "out", "spray_node")
RESULTS = os.path.join(HERE, "..", "results", "3d", "spray")
# the emission's constants as the node's float32 ports hold them; a low threshold and
# a high rate so a small tank emits within one step
E = {k: float(np.float32(v)) for k, v in dict(Sim.EMIT_DEFAULTS, grad_threshold=3.0, rate=0.5).items()}
CFL = float(np.float32(0.5))
SPRAY = {"emit_model": 1, "emit_grad_threshold": pf(E["grad_threshold"]), "emit_rate": pf(E["rate"]),
         "emit_max_per_cell": int(E["max_per_cell"]), "spray_cfl": pf(CFL)}
EXTRA = ("radius_in", "spray_in")
OUTS = ("radius_out", "spray_out")


def tank(nx=120, ny=40, depth=20.0, wind=0.0, shear=0.0, seed=1, **kw):
    """A still tank of the reference, escape and the droplet model on, calibrated; a
    wind in its air above depth + 2, a sheared layer in its top three cells of water."""
    sim = Sim(nx, ny, escape=True, droplet=D, **kw)
    sim.seed(lambda x, y: y < depth, seed=seed)
    sim.calibrate()
    if wind:
        sim.vel[(sim.typ == 0) & (sim.pos[:, 1] > depth + 2.0), 0] = wind
    if shear:
        sim.vel[(sim.typ == 1) & (sim.pos[:, 1] > depth - 3.0), 0] = shear
    return sim


def twin(state, ch, radius_in, spray_in, step_index, emit, cfl):
    """validate_level_set_node.twin_step on the node's inputs with the spray: the
    reference's step from the projection on, as Sim.step does it; `emit` the
    emission's constants or None. Returns the reference after the step and the
    particles deleted."""
    ref = Sim(state.nx, state.ny, escape=True, droplet=D, emit=emit, spray_cfl=cfl)
    arr = flat_inputs(state, ch)
    ref.pos = arr["positions"][:, :2].astype(float)
    ref.vel = arr["velocities"][:, :2].astype(float)
    ref.typ = state.typ.copy()
    ref.scale = np.ones(len(ref.pos))
    ref.tau = np.zeros(len(ref.pos))
    ref.xi = np.ones(len(ref.pos))
    ref.radius = radius_in.astype(float)
    ref.spray = spray_in >= 0.5
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
    # the escape test, the spray counted in no phase
    ones_air = (ref.typ == 0).astype(float)
    ones_liq = ((ref.typ == 1) & ~ref.spray).astype(float)
    fu, fv = ref._face_frames()
    au, lu = ref._splat(ref.u.shape, fu, [ones_air, ones_liq])
    av, lv = ref._splat(ref.v.shape, fv, [ones_air, ones_liq])
    thresh = (1.0 - ref.esc_phi) * ref.rho0_face
    fa = 0.5 * (ref._sample_faces(au, np.array([0.0, 0.5]), ref.pos) +
                ref._sample_faces(av, np.array([0.5, 0.0]), ref.pos))
    fl = 0.5 * (ref._sample_faces(lu, np.array([0.0, 0.5]), ref.pos) +
                ref._sample_faces(lv, np.array([0.5, 0.0]), ref.pos))
    drop = ((ref.typ == 1) & (fa > thresh)) | ref.spray
    bub = (ref.typ == 0) & (fl > thresh)
    inside = ref.spray & (fl > thresh)
    esc = drop | bub
    ref.escaped = esc
    lf = fl / ref.rho0_face
    ref.g2p(u_star, v_star, skip=esc)
    gone = ref._droplet_radii(drop, DT) | inside
    fly = None
    if esc.any():
        vg = ref.sample_velocity(ref.pos)
        if cfl > 0.0:
            fly = drop
        else:
            ref.vel[drop, 1] += ref.g * DT
            if drop.any():
                ref._droplet_drag(drop, vg, DT)
        ref.vel[bub, 1] += ref.buoyancy * abs(ref.g) * DT
        ref.vel[bub] += ref.drag_bubble * DT * (vg[bub] - ref.vel[bub])
    ref.advect(DT, ballistic=esc, fly=fly)
    emitted = 0
    if emit:
        emitted = ref._emit(DT, lf)
        if emitted:
            gone = np.concatenate([gone, np.zeros(emitted, bool)])
    if gone.any():
        ref._delete(~gone)
    if ref.spray.any():
        order = np.concatenate([np.flatnonzero(~ref.spray), np.flatnonzero(ref.spray)])
        for name in ("pos", "vel", "typ", "scale", "tau", "xi", "radius", "escaped", "spray"):
            setattr(ref, name, getattr(ref, name)[order])
    return ref, gone, emitted


def node_step(tag, sim, radius_in, spray_in, **ports):
    """One step_3d call at nz = 1 on the reference state, radius_in and spray_in as
    its state."""
    ch = channels(sim)
    arr = flat_inputs(sim, ch)
    arr["radius_in"] = radius_in.astype(np.float32)
    arr["spray_in"] = spray_in.astype(np.float32)
    d_in, d_out = os.path.join(OUT, tag + "_in"), os.path.join(OUT, tag + "_out")
    for k, v in arr.items():
        save(d_in, k, v)
    stats = run_node(d_in, d_out, extra_arrays=EXTRA, extra_outputs=OUTS,
                     nx=sim.nx, ny=sim.ny, nz=1, dt=DT, preconditioner=1,
                     rho0_face=pf(sim.rho0_face), **ESCAPE, **ports)
    return d_out, ch, stats


def compare(name, d_out, ref):
    """The node's outputs against the twin's: grid, positions, velocities, escape and
    spray flags, the particle count, the radii (within 1e-5 where positive)."""
    nx, ny = ref.nx, ref.ny
    pos = load(d_out, "out_positions")[:, :2]
    if len(pos) != len(ref.pos):
        print(f"    {name:28s} particle count {len(pos)} vs the reference's {len(ref.pos)}  FAIL")
        return False, {"count_node": len(pos), "count_ref": len(ref.pos)}
    vel = load(d_out, "out_velocities")[:, :2]
    rad = load(d_out, "radius_out")
    spr = load(d_out, "spray_out") > 0.5
    du = max(float(np.abs(load(d_out, "out_u").reshape(nx + 1, ny) - ref.u).max()),
             float(np.abs(load(d_out, "out_v").reshape(nx, ny + 1) - ref.v).max()))
    dp = float(np.abs(pos - ref.pos).max())
    dv = float(np.abs(vel - ref.vel).max())
    mism = int(((load(d_out, "out_escaped") > 0.5) != ref.escaped).sum())
    smism = int((spr != ref.spray).sum())
    drops = ref.radius > 0
    dr = float((np.abs(rad[drops] - ref.radius[drops]) / ref.radius[drops]).max()) if drops.any() else 0.0
    zeros = bool((rad[~drops] == 0).all()) and bool((rad[drops] > 0).all())
    good = mism <= 2 and smism == 0 and dp < 2e-3 and dv < 5e-2 and du < 1e-3 and dr < 1e-5 and zeros
    print(f"    {name:28s} grid {du:.1e}, positions {dp:.1e}, velocities {dv:.1e}, escape"
          f" mismatches {mism}, spray mismatches {smism}, {int(drops.sum())} droplets' radii"
          f" within {dr:.1e}, {int(spr.sum())} spray, {len(pos)} particles  {'OK' if good else 'FAIL'}")
    return good, {"grid": du, "positions": dp, "velocities": dv, "escape_mismatches": mism,
                  "spray_mismatches": smism, "droplets": int(drops.sum()), "radius_rel": dr,
                  "spray": int(spr.sum()), "count": len(pos)}


# ------------------------------------------------------------------ checks


def check_off():
    print("=== off: the spray ports at 0 are the droplet path, nz = 1", flush=True)
    sim = tank(wind=3.0)
    plant(sim, 60, 26.0, 36.0)
    zero = np.zeros(len(sim.pos), np.float32)
    d0, _, _ = node_step("off_zero", sim, zero, zero, emit_model=0, spray_cfl=pf(0.0), **PORTS)
    d1, _, _ = node_step("off_unset", sim, zero, zero, **PORTS)
    names = ("out_u", "out_v", "out_w", "pressure", "out_positions", "out_velocities",
             "out_escaped", "radius_out", "spray_out")
    same = all(np.array_equal(load(d0, n), load(d1, n)) for n in names)
    allzero = bool((load(d0, "spray_out") == 0).all())
    drops = int((load(d0, "radius_out") > 0).sum())
    good = same and allzero and drops >= 50
    print(f"    spray ports 0 against unset: {'bit-identical' if same else 'DIFFER'}; spray_out"
          f" all 0: {allzero}; {drops} droplets as before  {'OK' if good else 'FAIL'}")
    return good, {"unset_identical": same, "spray_all_zero": allzero, "droplets": drops}


def check_lockstep():
    print("=== lockstep: emission and the integrator on against the reference's twin, nz = 1",
          flush=True)
    sim = tank(wind=30.0, shear=10.0)
    idx = plant(sim, 60, 26.0, 36.0)
    radius = np.zeros(len(sim.pos))
    spray = np.zeros(len(sim.pos))
    radius[idx[:40]] = np.logspace(np.log10(2e-5), np.log10(2e-3), 40)
    # twenty of the planted carry the spray flag already: massless, flying on
    spray[idx[40:]] = 1.0
    radius[idx[40:]] = 1e-4
    n = 7
    # a rate that makes every candidate above the threshold emit, so one step
    # exercises the emission and the per-cell cap
    rate = float(np.float32(50.0))
    d_out, ch, stats = node_step("lock", sim, radius, spray, step_index=n, **PORTS,
                                 **dict(SPRAY, emit_rate=pf(rate)))
    ref, gone, emitted = twin(sim, ch, radius, spray, n, emit=dict(E, rate=rate), cfl=CFL)
    good, rec = compare("spray twin", d_out, ref)
    print(f"    the twin emitted {emitted}, deleted {int(gone.sum())}; spray alive"
          f" {int(ref.spray.sum())}; iterations {stats['iterations_used']:.0f}", flush=True)
    good &= emitted >= 5
    rec.update(emitted=emitted, deleted=int(gone.sum()), iterations=stats["iterations_used"])
    return good, rec


def check_anchors():
    print("=== anchors: the laws on the reference", flush=True)
    ok, rec = True, {}
    # a still tank emits nothing
    s = tank(emit=dict(E), spray_cfl=CFL)
    tot = 0
    for _ in range(20):
        s.step(DT)
        tot += s.n_emitted
    good = tot == 0
    print(f"    still tank, 20 steps: {tot} emitted  {'OK' if good else 'FAIL'}")
    ok &= good
    rec["still_emitted"] = tot
    # the count falls as the threshold rises: the others are the same whatever the
    # spray, so the candidates and their hash draws are the same run to run and a
    # higher threshold can only lose emissions
    counts = []
    for thr in (2.0, 6.0, 12.0):
        s = tank(wind=30.0, shear=10.0, emit=dict(E, grad_threshold=thr), spray_cfl=CFL)
        tot = 0
        for _ in range(40):
            s.step(DT)
            tot += s.n_emitted
        counts.append(tot)
    good = counts[0] >= counts[1] >= counts[2] >= 0 and counts[0] > counts[2]
    print(f"    emitted over 40 steps at thresholds 2, 6, 12 per second: {counts}  {'OK' if good else 'FAIL'}")
    ok &= good
    rec["by_threshold"] = counts
    # the integrator's order: a 20 um drop launched into a light wind under gravity,
    # its flight at sub-steps halving against the flight at sub-steps a hundred
    # times finer (the drag's rate changes with the slip, so no closed form: the fine
    # flight stands in for it); the position error falls fourfold. Asserted in
    # Stokes' regime (launched at 0.02 cells/s, Re 0.05); the fast launch (4 cells/s,
    # Re 10, the rate varying strongly within a sub-step) reported, its ratios
    # climbing toward 4 as the sub-steps shrink
    r0, wind, T = 2e-5, 1e-3, 0.012

    def flight(cfl, v0):
        s = Sim(8, 8, escape=True, droplet=dict(D, metres_per_cell=1.0), spray_cfl=cfl)
        s.pos, s.vel = np.array([[2.0, 4.0]]), np.array([[v0, 0.0]])
        s.typ, s.radius, s.spray = np.array([1]), np.array([r0]), np.array([False])
        s.u, s.v = np.full((9, 8), wind), np.zeros((8, 9))
        s._droplet_fly(np.array([True]), np.array([T]))
        return s.pos[0].copy()

    k = 4.5 * (Sim(4, 4).rho_g / Sim(4, 4).rho_l) * D["air_nu"] / r0 ** 2   # Stokes' rate, 1/s
    out = {}
    for name, v0, fine, cfls in (("stokes", 0.02, 1e-6, (1e-4, 5e-5, 2.5e-5)),
                                 ("fast", 4.0, 1e-4, (0.01, 0.005, 0.0025))):
        truth = flight(fine, v0)
        errs = [float(np.abs(flight(c, v0) - truth).sum()) for c in cfls]
        out[name] = {"errors": errs, "ratios": [errs[i] / max(errs[i + 1], 1e-300) for i in range(2)]}
    good = all(3.5 < r < 4.5 for r in out["stokes"]["ratios"])
    print(f"    integrator on a 20 um drop (tau {1e3 / k:.2f} ms, 2 tau of flight), position errors at"
          f" sub-steps halving: Stokes {[f'{e:.1e}' for e in out['stokes']['errors']]} cells, ratios"
          f" {[f'{r:.2f}' for r in out['stokes']['ratios']]} (4 is second order); fast"
          f" {[f'{e:.1e}' for e in out['fast']['errors']]}, ratios"
          f" {[f'{r:.2f}' for r in out['fast']['ratios']]}  {'OK' if good else 'FAIL'}")
    ok &= good
    rec["integrator"] = dict(out, tau_ms=1e3 / k)
    # the others bit-identical with emission on and off
    same_all = True
    for st, adapt in ((False, False), (True, False), (False, True), (True, True)):
        a = tank(wind=30.0, shear=10.0, spray_cfl=CFL, st=st, adapt=adapt)
        b = tank(wind=30.0, shear=10.0, spray_cfl=CFL, st=st, adapt=adapt, emit=dict(E))
        tot = 0
        for _ in range(40):
            a.step(DT)
            b.step(DT)
            tot += b.n_emitted
        prim = ~b.spray
        same = (np.array_equal(a.pos, b.pos[prim]) and np.array_equal(a.vel, b.vel[prim])
                and np.array_equal(a.u, b.u))
        same_all &= same and tot > 0
        print(f"    st {int(st)} adapt {int(adapt)}, 40 steps: {tot} emitted, the others"
              f" {'bit-identical' if same else 'DIFFER'}")
    print(f"    the particles that are not spray, emission on against off  {'OK' if same_all else 'FAIL'}")
    ok &= same_all
    rec["others_identical"] = same_all
    # the long run, against the same run without emission (whose own droplets may
    # evaporate away over 20 s: the count to match is its)
    s = tank(wind=30.0, shear=10.0, emit=dict(E), spray_cfl=CFL)
    s0 = tank(wind=30.0, shear=10.0, spray_cfl=CFL)
    n0 = len(s.pos)
    alive, tot, vmax = [], 0, 0.0
    for _ in range(800):
        s.step(DT)
        s0.step(DT)
        tot += s.n_emitted
        alive.append(int(s.spray.sum()))
        vmax = max(vmax, float(np.abs(s.vel).max()))
    finite = bool(np.isfinite(s.pos).all() and np.isfinite(s.vel).all())
    others = int((~s.spray).sum())
    same = others == len(s0.pos) and np.array_equal(s.pos[~s.spray], s0.pos)
    good = finite and same and max(alive) < n0 // 2 and vmax < 1e3
    print(f"    800 steps (20 s): emitted {tot}, spray alive at most {max(alive)} and {alive[-1]} at"
          f" the end, the others {others} (seeded {n0}; {'the same as' if same else 'NOT'} without"
          f" emission, {len(s0.pos)}), max speed {vmax:.1f} cells/s, finite {finite}"
          f"  {'OK' if good else 'FAIL'}")
    ok &= good
    rec["long"] = {"emitted": tot, "alive_max": max(alive), "alive_end": alive[-1], "others": others,
                   "seeded": n0, "same_as_without": same, "vmax": vmax, "finite": finite}
    return ok, rec


def check_loop():
    print("=== loop: the 3D loop resumed from a sheared state with emission on and off, 40 x 40 x 8",
          flush=True)
    nx, ny, nz, dt = 40, 40, 8, 0.02
    cal = Sim(40, 40, escape=True)
    cal.seed(lambda x, y: np.ones_like(x, dtype=bool))
    cal.calibrate()

    def seed_fn(nx_, ny_, nz_, mask):
        pos, typ = seed3d(nx_, ny_, nz_, mask)
        rng = np.random.default_rng(5)
        idx = np.flatnonzero((typ == 1) & (pos[:, 1] < 4.0))[:200]
        pos[idx, 0] = rng.uniform(4.0, nx_ - 4.0, len(idx))
        pos[idx, 1] = rng.uniform(26.0, 36.0, len(idx))
        pos[idx, 2] = rng.uniform(1.0, nz_ - 1.0, len(idx))
        return pos, typ

    base = dict(ESCAPE, **PORTS)
    base.update(rho0_face=pf(cal.rho0_face), metres_per_cell=pf(0.1))
    water = lambda x, y, z: y < 12  # noqa: E731
    # two steps to a state, then a wind over the water and a sheared top layer
    d0, _, _ = run_bifrost3("spray_state", nx, ny, nz, water, 2, dt, extra_ports=base,
                            seed_fn=seed_fn, fused=True, dump_every=1, final_state=True,
                            out_dir=os.path.join(OUT, "loop_state"))
    init = load_final_state(d0)
    pos, vel = init["positions"], init["velocities"].copy()
    liq = init["phase_st"] >= 0.5
    vel[(~liq) & (pos[:, 1] > 14.0), 0] = 30.0
    vel[liq & (pos[:, 1] > 9.0), 0] = 10.0
    init["velocities"] = vel.astype(np.float32)
    for k, ax in (("momx", 0), ("momy", 1), ("momz", 2)):
        init[k] = (init["mass"] * vel[:, ax]).astype(np.float32)
    on = dict(base, **dict(SPRAY, emit_grad_threshold=pf(2.0)))
    off = dict(base, spray_cfl=pf(CFL))
    d_on, _, _ = run_bifrost3("spray_on", nx, ny, nz, water, 40, dt, extra_ports=on, fused=True,
                              dump_every=1, init_arrays=init, step_offset=2, final_state=True,
                              out_dir=os.path.join(OUT, "loop_on"))
    d_off, _, _ = run_bifrost3("spray_off", nx, ny, nz, water, 40, dt, extra_ports=off, fused=True,
                               dump_every=1, init_arrays=init, step_offset=2,
                               out_dir=os.path.join(OUT, "loop_off"))
    ok, rec = True, {}
    n0 = len(pos)
    # a spray particle draws its radius at its first step as a droplet, the step
    # after its emission, and the spray sits at the tail in emission order, so at a
    # dump the spray without a radius are exactly the last ones: once a 0 appears
    # among the spray's radii every radius after it is 0; the others' count follows
    # the run without emission (its own droplets may evaporate away)
    placed, others_match, emitted_any = True, True, False
    sprays, fresh = [], []
    for k in range(40):
        spr = np.load(os.path.join(d_on, f"spr.{k:04d}.npy")) > 0.5
        rad = np.load(os.path.join(d_on, f"rad.{k:04d}.npy"))
        z = rad[spr] <= 0
        placed &= bool(not z.any() or z[int(np.argmax(z)):].all())
        fresh.append(int(z.sum()))
        n_off = len(np.load(os.path.join(d_off, f"frame.{k:04d}.npy")))
        others_match &= int((~spr).sum()) == n_off
        sprays.append(int(spr.sum()))
    emitted_any = max(sprays) > 0
    fin_mass = np.load(os.path.join(d_on, "final_mass.npy"))
    fin_spr = np.load(os.path.join(d_on, "final_spray.npy")) > 0.5
    massless = bool(np.array_equal(fin_mass == 0, fin_spr))
    good = placed and others_match and emitted_any and massless
    print(f"    spray over 40 steps: {sprays[0]} -> {sprays[-1]} (at most {max(sprays)}, at most"
          f" {max(fresh)} fresh at a dump); radii positive on every spray past its first step {placed};"
          f" the others' count (seeded {n0}) the run without emission's at every dump {others_match};"
          f" mass 0 exactly on the spray at the end {massless}  {'OK' if good else 'FAIL'}", flush=True)
    ok &= good
    rec.update(spray_first=sprays[0], spray_last=sprays[-1], spray_max=max(sprays), fresh_max=max(fresh),
               placed=placed, others_match_without=others_match, massless=massless)
    # the others bit-identical with emission off
    p_on = np.load(os.path.join(d_on, "frame.0039.npy"))
    spr = np.load(os.path.join(d_on, "spr.0039.npy")) > 0.5
    p_off = np.load(os.path.join(d_off, "frame.0039.npy"))
    v_on = np.load(os.path.join(d_on, "vel.0039.npy")) if os.path.exists(os.path.join(d_on, "vel.0039.npy")) else None
    same = len(p_off) == int((~spr).sum()) and np.array_equal(p_on[~spr], p_off)
    print(f"    the others at step 40, emission on against off: {'bit-identical' if same else 'DIFFER'}"
          f" ({int((~spr).sum())} particles)  {'OK' if same else 'FAIL'}", flush=True)
    ok &= same
    rec["others_identical_3d"] = same
    # a resume reads the spray back
    init2 = load_final_state(d_on)
    d_re, _, _ = run_bifrost3("spray_resume", nx, ny, nz, water, 3, dt, extra_ports=on, fused=True,
                              dump_every=1, init_arrays=init2, step_offset=42,
                              out_dir=os.path.join(OUT, "loop_resume"))
    spr0 = np.load(os.path.join(d_re, "spr.0000.npy")) > 0.5
    good = "spray" in init2 and int(init2["spray"].sum()) == sprays[-1] and spr0.sum() >= 1
    print(f"    resume: {int(init2['spray'].sum())} spray in the state (the last dump's {sprays[-1]}),"
          f" {int(spr0.sum())} on the resumed first step  {'OK' if good else 'FAIL'}", flush=True)
    ok &= good
    rec["resumed_spray"] = int(spr0.sum())
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
    print("SPRAY NODE:", "PASS" if ok else "FAIL", flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
