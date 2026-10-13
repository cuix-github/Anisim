"""Milestone 6.5a stage-3c gate: the adaptive solve in 3D (an octree), both
junction treatments.

There is no 3D reference to lockstep against, so the 3D path is checked by
reductions, structure and physics anchors, each on one step of step_3d
through solve_step3_test.json (multigrid for the uniform solve):

  refined     adaptive_band wide enough that every cell stays fine: the tree
              is the uniform grid, so the adaptive solve must reproduce the
              uniform one. At 1e-10 the two Krylov solves must agree to 1e-7
              of the faces, which pins the 3D faces, coefficients,
              right-hand side and multigrid coarsening.
  exact       the leaf operator, dumped through PFFLIP_ADAPTIVE_DUMP on a 3D
              dam tree: every cleanup junction face's weighted difference
              must reproduce linear pressure along x, y and z exactly - the
              two-tangent sample's whole point. Naive's two-cell junction
              faces miss by the lateral offset (s/2), reported.
  dynamic     a 48x32x32 dam with a pool and a divergence-free swirl at
              1000:1: per junction, the divergence left after projection
              must sit at solver level (the local fills close every fine
              cell), and the liquid faces are reported against uniform's.
  z_invariant the swirl's x-y part alone, the same in every z slab, at nz=32
              and at nz=1: in the middle slab, away from the z walls' bands,
              the 3D solve must reproduce the nz=1 solve (which is the 2D
              reference's, validate_adaptive_step.py) within 1e-3.
  still       a half-full 32x64x32 tank at rest after the gravity kick, and
              the same tank at nz=1: uniform and cleanup must leave it still
              to 1e-4 of the kick; naive's first-order junction error is
              reported (dt g / 3 in the 2D record).

Run from gates:  ..\\.venv\\Scripts\\python.exe validate_adaptive_3d.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
from validate_3d import run_node  # noqa: E402
from loop3d import seed3d  # noqa: E402
from pfflip2d import Sim  # noqa: E402

OUT = os.path.join(HERE, "..", "out", "adaptive_3d")
RHO_L, RHO_G, G = 1000.0, 1.0, -9.8
JUNCTIONS = (("naive", 1), ("cleanup", 2))


def face_centres(nx, ny, nz):
    """Face-centre coordinates (x, y, z grids) of the u, v and w faces."""
    def grid(ni, nj, nk, off):
        return np.meshgrid(np.arange(ni) + off[0], np.arange(nj) + off[1],
                           np.arange(nk) + off[2], indexing="ij")
    return {"u": grid(nx + 1, ny, nz, (0.0, 0.5, 0.5)),
            "v": grid(nx, ny + 1, nz, (0.5, 0.0, 0.5)),
            "w": grid(nx, ny, nz + 1, (0.5, 0.5, 0.0))}


def package(d, nx, ny, nz, mask, field, seed=1):
    """Analytic face channels for a liquid mask and a velocity field, and
    particles carrying the field: seeded 8 per cell (seed3d, as the 3D loops),
    or 4 per cell at nz = 1 (the 2D seeding the node's liquid count expects)."""
    os.makedirs(d, exist_ok=True)
    for comp, (x, y, z) in face_centres(nx, ny, nz).items():
        phi = mask(x, y, z).astype(float)
        mass = RHO_G + (RHO_L - RHO_G) * phi
        vel = field(x, y, z)["uvw".index(comp)]
        np.save(os.path.join(d, f"{comp}_mass.npy"), mass.ravel().astype(np.float32))
        np.save(os.path.join(d, f"{comp}_mom.npy"), (mass * vel).ravel().astype(np.float32))
        np.save(os.path.join(d, f"{comp}_phase.npy"), phi.ravel().astype(np.float32))
    if nz == 1:
        sim = Sim(nx, ny)
        sim.seed(lambda x, y: mask(x, y, 0.5 + 0.0 * x), seed=seed)
        pos = np.column_stack([sim.pos, np.full(len(sim.pos), 0.5)])
        typ = sim.typ
    else:
        pos, typ = seed3d(nx, ny, nz, mask, seed=seed)
    vel = np.stack(field(pos[:, 0], pos[:, 1], pos[:, 2]), 1)
    np.save(os.path.join(d, "positions.npy"), pos.astype(np.float32))
    np.save(os.path.join(d, "velocities.npy"), vel.astype(np.float32))
    np.save(os.path.join(d, "particle_phase.npy"), typ.astype(np.float32))


def faces(d, nx, ny, nz):
    return (np.load(os.path.join(d, "out_u.npy")).reshape(nx + 1, ny, nz),
            np.load(os.path.join(d, "out_v.npy")).reshape(nx, ny + 1, nz),
            np.load(os.path.join(d, "out_w.npy")).reshape(nx, ny, nz + 1))


def step(indir, outdir, nx, ny, nz, solver, dt=0.02, **ports):
    st = run_node(indir, outdir, nx=nx, ny=ny, nz=nz, dt=dt, pressure_solver=solver,
                  preconditioner=1, **ports)
    return faces(outdir, nx, ny, nz), st


def against(f0, f1, masks):
    """Max and rms-relative difference of f1 from f0 over the masked faces."""
    mx = max(float(np.abs(a - b)[m].max()) for a, b, m in zip(f0, f1, masks) if m.any())
    rel = np.sqrt(sum(((a - b)[m] ** 2).sum() for a, b, m in zip(f0, f1, masks))
                  / sum((a[m] ** 2).sum() for a, m in zip(f0, masks)))
    return mx, float(rel)


A_SWIRL, B_SWIRL = 60.0, 40.0


def swirl_xy(x, y, z):
    """The curl of (0, 0, A sin(pi x/48) sin(pi y/32)): divergence-free, no
    normal flow through the 48x32 box's walls, the same in every z slab."""
    u = A_SWIRL * np.pi / 32 * np.sin(np.pi * x / 48) * np.cos(np.pi * y / 32)
    v = -A_SWIRL * np.pi / 48 * np.cos(np.pi * x / 48) * np.sin(np.pi * y / 32)
    return u, v, 0.0 * x


def swirl(x, y, z):
    """swirl_xy plus the curl of (B sin(pi y/32) sin(pi z/32), 0, 0): a
    divergence-free 3D swirl with no normal flow through any wall of the
    48x32x32 box, up to ~6 cells per second. A real step's u* is nearly
    divergence-free in the bulk too."""
    u, v, _ = swirl_xy(x, y, z)
    v = v + B_SWIRL * np.pi / 32 * np.sin(np.pi * y / 32) * np.cos(np.pi * z / 32)
    w = -B_SWIRL * np.pi / 32 * np.cos(np.pi * y / 32) * np.sin(np.pi * z / 32)
    return u, v, w


def dam_mask(x, y, z):
    return (y < 6) | ((x < 16) & (y < 24))


def liquid(nx, ny, nz):
    cen = face_centres(nx, ny, nz)
    return [dam_mask(*cen[c]) for c in "uvw"]


def refined():
    print("=== refined (48x32x32 dam, every cell fine: adaptive must reproduce uniform)")
    nx, ny, nz = 48, 32, 32
    d = os.path.join(OUT, "dam_in")
    package(d, nx, ny, nz, dam_mask, swirl)
    ok = True
    for tol in (1e-6, 1e-10):
        res = {s: step(d, os.path.join(OUT, f"refined_ps{s}_tol{tol:g}"), nx, ny, nz, s,
                       adaptive_band=48, tolerance=tol) for s in (0, 1)}
        diff = max(float(np.abs(a - b).max()) for a, b in zip(res[0][0], res[1][0]))
        scale = max(float(np.abs(f).max()) for f in res[0][0])
        line = (f"  tolerance {tol:g}: iterations uniform {res[0][1]['iterations_used']:.0f}, "
                f"adaptive {res[1][1]['iterations_used']:.0f}; max face diff {diff:.2e} "
                f"(faces up to {scale:.1f})")
        if tol == 1e-10:
            good = diff < 1e-7 * max(scale, 1.0)
            ok &= good
            line += "  OK" if good else "  FAIL"
        print(line)
    print("refined:", "PASS" if ok else "FAIL")
    return ok


def read_dump(path):
    cells, fcs = [], []
    for line in open(path):
        t = line.split()
        if t[0] == "C":
            cells.append(tuple(int(x) for x in t[1:]))
        else:
            en = int(t[5])
            fcs.append((int(t[1]), int(t[2]), int(t[3]), float(t[4]),
                        [(int(t[6 + 2 * e]), float(t[7 + 2 * e])) for e in range(en)]))
    return cells, fcs


def exact():
    print("=== exact (the 3D dam tree's leaf operator against linear pressure)")
    nx, ny, nz = 48, 32, 32
    d = os.path.join(OUT, "dam_in")
    ok = True
    for name, solver in JUNCTIONS:
        path = os.path.join(OUT, f"dump_ps{solver}.txt")
        os.environ["PFFLIP_ADAPTIVE_DUMP"] = path
        try:
            step(d, os.path.join(OUT, f"exact_ps{solver}"), nx, ny, nz, solver)
        finally:
            del os.environ["PFFLIP_ADAPTIVE_DUMP"]
        cells, fcs = read_dump(path)
        cen = np.array([[(i + 0.5) * (1 << l), (j + 0.5) * (1 << l), (k + 0.5) * (1 << l)]
                        for (l, i, j, k) in cells])
        lev = np.array([c[0] for c in cells])
        worst, n = 0.0, 0
        for axis, a, b, dist, ent in fcs:
            if not ent:
                if lev[a] == lev[b]:
                    continue
                ent = [(a, -1.0), (b, 1.0)]
            n += 1
            # sum_e w_e x_e must be the normal center distance along the axis
            err = sum(w * cen[i] for i, w in ent) - dist * np.eye(3)[axis]
            worst = max(worst, float(np.abs(err).max()))
        line = f"  {name}: {len(cells)} cells, {n} junction faces, worst miss {worst:.3g} cells"
        if name == "cleanup":
            good = worst < 1e-9
            ok &= good
            line += "  OK" if good else "  FAIL"
        else:
            line += "  (the two-cell stencil's lateral offset, s/2)"
        print(line)
    print("exact:", "PASS" if ok else "FAIL")
    return ok


def dynamic():
    print("=== dynamic (48x32x32 dam with a pool and a divergence-free swirl, 1000:1)")
    nx, ny, nz = 48, 32, 32
    d = os.path.join(OUT, "dam_in")
    f0, s0 = step(d, os.path.join(OUT, "dynamic_ps0"), nx, ny, nz, 0)
    print(f"  uniform: {s0['iterations_used']:.0f} iterations, divergence after "
          f"{s0['max_divergence_after']:.2e}")
    ok = True
    for name, solver in JUNCTIONS:
        f1, st = step(d, os.path.join(OUT, f"dynamic_ps{solver}"), nx, ny, nz, solver)
        mx, rel = against(f0, f1, liquid(nx, ny, nz))
        div = st["max_divergence_after"]
        good = div < 1e-4
        ok &= good
        print(f"  {name}: {st['iterations_used']:.0f} iterations, divergence after {div:.2e}  "
              f"{'OK' if good else 'FAIL'}; liquid faces against uniform: max {mx:.3e}, "
              f"rms relative {rel:.2e}")
    print("dynamic:", "PASS" if ok else "FAIL")
    return ok


def z_invariant():
    print("=== z_invariant (the x-y swirl in every slab: nz=32 against nz=1)")
    nx, ny = 48, 32
    runs = {}
    for nz in (1, 32):
        d = os.path.join(OUT, f"zinv{nz}_in")
        package(d, nx, ny, nz, dam_mask, swirl_xy)
        for s in (0, 1, 2):
            runs[nz, s] = step(d, os.path.join(OUT, f"zinv{nz}_ps{s}"), nx, ny, nz, s)
    ok = True
    for name, solver in JUNCTIONS:
        f3, st = runs[32, solver]
        f2 = runs[1, solver][0]
        liq2 = liquid(nx, ny, 1)
        mid = max(float(np.abs(f3[i][..., 16] - f2[i][..., 0])[liq2[i][..., 0]].max())
                  for i in range(2))
        m3 = against(runs[32, 0][0], f3, liquid(nx, ny, 32))
        m2 = against(runs[1, 0][0], f2, liq2)
        good = mid < 1e-3 and st["max_divergence_after"] < 1e-4
        ok &= good
        print(f"  {name}: middle slab against nz=1 {mid:.2e}  {'OK' if good else 'FAIL'}; "
              f"liquid faces against uniform, max / rms relative: nz=32 {m3[0]:.3f} / "
              f"{m3[1]:.2e}, nz=1 {m2[0]:.3f} / {m2[1]:.2e}")
    print("z_invariant:", "PASS" if ok else "FAIL")
    return ok


def still():
    ok = True
    for nz in (32, 1):
        nx, ny, dt = 32, 64, 0.02
        print(f"=== still ({nx}x{ny}x{nz} half-full tank at rest, one step after the gravity kick)")
        d = os.path.join(OUT, f"tank{nz}_in")
        package(d, nx, ny, nz, lambda x, y, z: y < 32, lambda x, y, z: (0.0 * x,) * 3)
        kick = abs(G) * dt
        for name, solver in (("uniform", 0),) + JUNCTIONS:
            fu, st = step(d, os.path.join(OUT, f"still{nz}_ps{solver}"), nx, ny, nz, solver, dt)
            worst = max(float(np.abs(f).max()) for f in fu)
            line = (f"  {name}: {st['iterations_used']:.0f} iterations, max face speed after "
                    f"{worst:.3e} = {worst / kick:.2e} of the kick ({kick:.3f})")
            if name == "naive":
                line += f"  (dt g / 3 = {kick / 3:.3f})"
            else:
                good = worst < 1e-4 * kick
                ok &= good
                line += "  OK" if good else "  FAIL"
            print(line)
    print("still:", "PASS" if ok else "FAIL")
    return ok


def main():
    results = [refined(), exact(), dynamic(), z_invariant(), still()]
    print("ALL PASS" if all(results) else "FAILURES PRESENT")
    sys.exit(0 if all(results) else 1)


if __name__ == "__main__":
    main()
