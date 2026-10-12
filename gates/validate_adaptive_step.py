"""Milestone 6.5a stage-3b gate: the adaptive solve inside step_3d, in nz = 1
lockstep against the 2D reference.

step_3d with pressure_solver = 1 or 2 runs nodes/src/PFFlipAdaptive.cpp, the
C++ port of reference2d/adaptive_projection.py: 1 with junction = "naive"
(the PF-FLIP paper's junction stencil), 2 with junction = "cleanup" (the
variational solve, then a naive-stencil projection of its leftover). At
nz = 1 its tree is the reference's quadtree, so one node step must reproduce
one reference step, for each junction. Four cases, each giving the solver
a different tree:

  dam_mid        the validate_3d dam-break mid-state (160x80, frame 12)
  dam_late       the same dam break at frame 30: a splash, a different tree
  dam_late_deep  frame 30 again with adaptive_band 1 and adaptive_levels 3:
                 a three-level tree whose corners put a junction's lateral
                 neighbor one level coarser (16 such junctions on the leaf
                 level, 7 on the first multigrid level), the cleanup
                 stencil's lateral case the default trees never build
  tank           a half-full 64x64 tank after 10 frames at rest, where the
                 naive stencil's junction error is the only thing moving

The reference is fed the float32-rounded particle positions the node's ports
carry, so both sides count the same particles per cell and build the same
tree. Per state and junction: the same pressure iteration count (cleanup:
both solves together, the node's iterations_used); corrected faces and
particle velocities within 1e-5 and positions within 1e-4 of the reference
(float32 port rounding is ~1e-6 and ~1e-5 at these magnitudes; validate_3d's
uniform lockstep allows 1e-3); and the node's divergence after projection at
solver level. The uniform path (pressure_solver 0, the
default) is covered by validate_3d.py and the existing loop gates.

Run from gates:  ..\\.venv\\Scripts\\python.exe validate_adaptive_step.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "reference2d"))
from adaptive_projection import AdaptiveProjection  # noqa: E402
from pfflip2d import Sim  # noqa: E402
from validate_3d import run_node  # noqa: E402
from validate_step import channels  # noqa: E402

OUT = os.path.join(HERE, "..", "out", "adaptive_step")
DT = 0.02


def dam(frames):
    sim = Sim(160, 80)
    sim.seed(lambda x, y: (x < 40) & (y < 56))
    sim.calibrate()
    sim.run(frames, 0.06)
    return sim


def tank():
    sim = Sim(64, 64)
    sim.seed(lambda x, y: y < 32)
    sim.calibrate()
    sim.run(10, 0.05)
    return sim


def package(sim, indir):
    """Node inputs at nz = 1; returns the float32 channels and positions."""
    os.makedirs(indir, exist_ok=True)
    nx, ny = sim.nx, sim.ny
    ch = {k: v.astype(np.float32) for k, v in channels(sim).items()}
    for k, v in ch.items():
        np.save(os.path.join(indir, k + ".npy"), v.ravel())
    for k in ("w_mass", "w_mom", "w_phase"):
        np.save(os.path.join(indir, k + ".npy"), np.zeros(nx * ny * 2, np.float32))
    pos32 = sim.pos.astype(np.float32)
    np.save(os.path.join(indir, "positions.npy"),
            np.concatenate([pos32, np.full((len(pos32), 1), 0.5, np.float32)], 1))
    np.save(os.path.join(indir, "velocities.npy"),
            np.concatenate([sim.vel, np.zeros((len(sim.pos), 1))], 1).astype(np.float32))
    np.save(os.path.join(indir, "particle_phase.npy"), sim.typ.astype(np.float32))
    return ch, pos32.astype(float)


def reference_step(sim, ch, pos, junction, tree):
    """One reference step on the node's inputs, the adaptive projection with
    the given junction in place of the uniform one (validate_3d.python_step
    otherwise)."""
    ref = Sim(sim.nx, sim.ny)
    ref.pos, ref.vel, ref.typ = pos.copy(), sim.vel.copy(), sim.typ.copy()
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
    proj = AdaptiveProjection(ref, junction=junction, **tree)
    proj(DT)
    ref.g2p(u_star, v_star)
    ref.advect(DT)
    return ref, proj.log[-1]


JUNCTIONS = (("naive", 1), ("cleanup", 2))   # reference junction, pressure_solver


def check(name, sim, **ports):
    """ports: adaptive_band / adaptive_levels for the node (absent: its
    defaults, the reference's band 5 and levels 2)."""
    tree = {"band": ports.get("adaptive_band", 5), "levels": ports.get("adaptive_levels", 2)}
    print(f"=== {name} ({sim.nx}x{sim.ny})")
    indir = os.path.join(OUT, name + "_in")
    ch, pos = package(sim, indir)
    nx, ny = sim.nx, sim.ny
    ok = True
    for junction, solver in JUNCTIONS:
        outdir = os.path.join(OUT, f"{name}_{junction}_out")
        ref, log = reference_step(sim, ch, pos, junction, tree)
        stats = run_node(indir, outdir, nx=nx, ny=ny, nz=1, dt=DT, pressure_solver=solver,
                         **ports)
        bu = np.load(os.path.join(outdir, "out_u.npy")).reshape(nx + 1, ny)
        bv = np.load(os.path.join(outdir, "out_v.npy")).reshape(nx, ny + 1)
        bp = np.load(os.path.join(outdir, "out_positions.npy"))[:, :2]
        bw = np.load(os.path.join(outdir, "out_velocities.npy"))[:, :2]
        it_node, it_ref = int(stats["iterations_used"]), log["iters"] + log["iters_cleanup"]
        rows = [("face u", np.abs(bu - ref.u).max(), 1e-5), ("face v", np.abs(bv - ref.v).max(), 1e-5),
                ("positions", np.abs(bp - ref.pos).max(), 1e-4),
                ("velocities", np.abs(bw - ref.vel).max(), 1e-5)]
        split = f" ({log['iters']} + {log['iters_cleanup']})" if junction == "cleanup" else ""
        print(f"  -- {junction} (pressure_solver {solver}): {log['unknowns']} adaptive cells for "
              f"{log['uniform']} uniform ({1 - log['unknowns'] / log['uniform']:.0%} saved); "
              f"iterations node {it_node}, reference {it_ref}{split}")
        good = it_node == it_ref
        ok &= good
        if not good:
            print("    iteration counts differ  FAIL")
        for what, err, bar in rows:
            good = err < bar
            ok &= good
            print(f"    max |diff| {what}: {err:.3e}  {'OK' if good else 'FAIL'}")
        div = stats["max_divergence_after"]
        ref_div = max(log["div_fine"], log["div_coarse"])
        good = div < 1e-4
        ok &= good
        print(f"    divergence after: node {div:.2e}, reference {ref_div:.2e}  {'OK' if good else 'FAIL'}")
    print(f"{name}:", "PASS" if ok else "FAIL")
    return ok


def main():
    late = dam(30)
    results = [check("dam_mid", dam(12)), check("dam_late", late),
               check("dam_late_deep", late, adaptive_band=1, adaptive_levels=3),
               check("tank", tank())]
    print("ALL PASS" if all(results) else "FAILURES PRESENT")
    sys.exit(0 if all(results) else 1)


if __name__ == "__main__":
    main()
