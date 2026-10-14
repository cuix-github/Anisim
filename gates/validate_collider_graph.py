"""Mesh colliders, stage 3: the graph from a mesh to solid_phi, through the shipped
Simulation::Common::collider node and our collider_field adapter (stock nodes only),
held to exact distances and to the node's own box path. Light: a few minutes.

  lattice   the graph's node positions are exactly origin + (i, j, k) * cell_size, in
            the node's order, for an off-grid origin and a fractional cell size
  field     a sphere and a box as meshes, in world units, against their exact signed
            distances: within one cell of the surface the field is exact to 0.02
            cells (the band the node reads for face areas and push-out), within two
            cells to 0.75 (a cut face's far end; the slack is at a box's edges and
            mid-planes); the sign is right everywhere; beyond that the stock level set
            is a swept approximation (reported, held to 0.5x..1.5x)
  settings  the shipped node's knobs move the surface as documented: thickening by
            that many voxels, offset by that many scene units, a coarser detail size
            keeps the surface within half a voxel; two colliders give the minimum of
            the two fields exactly
  step      the pillar from the stage-2 gate as a mesh, a quarter cell off the nodes
            and face centres (on them the decision is a toss-up): with
            binary faces the mesh field gives the analytic field's face areas and one
            step's grid velocities exactly; with cut cells the face areas agree within
            0.05 and the step's velocities within 0.05 cells/s. The flow runs into the
            pillar at 25 cells/s, so the step pushes particles out (the count is
            reported), and none is inside afterwards
  refuse    a field the node cannot use stops it with an error (a float64 .npy, which
            the float32 port reads as empty; a field one node short), instead of a
            step without the solid; the right field runs
  loop      still water over a submerged ball, the ball from the graph: as still as
            the plain tank (the stage-2 bar), no particle inside

Run from the repo root:  .venv\\Scripts\\python.exe gates\\validate_collider_graph.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
from colliders import box_mesh, solid_phi_from_meshes, sphere_mesh  # noqa: E402
from validate_level_set_node import (f32, flow_inputs, inside, load, node,  # noqa: E402
                                     node_fractions, node_grid3, rest_3d, sdf_box3)

TAG = "cgr_"


def grid_world(nx, ny, nz, dx, origin):
    i, j, k = node_grid3(nx, ny, nz)
    return (origin[0] + i * dx, origin[1] + j * dx, origin[2] + k * dx)


def compare(name, phi, exact_cells):
    """phi and the exact field both in cells. Within one cell of the surface, the band
    the node reads for face areas and push-out, the field must be exact to 0.02 cells;
    within two cells (a cut face's far end) to 0.5; the sign right everywhere. Beyond
    that the stock level set is a swept approximation: reported, held to 0.5x..1.5x."""
    err = np.abs(phi - exact_cells)
    d = np.abs(exact_cells)
    on = exact_cells == 0.0
    signs = (np.sign(phi) == np.sign(exact_cells)) | on
    far = d > 3.0
    ratio = (phi[far] / exact_cells[far]) if far.any() else np.array([1.0])
    e1, e2, e3 = err[d <= 1.0].max(), err[d <= 2.0].max(), err[d <= 3.0].max()
    good = (e1 <= 0.02 and e2 <= 0.75 and bool(signs.all())
            and 0.5 <= ratio.min() and ratio.max() <= 1.5)
    print(f"    {name:26s} |err| within 1 cell {e1:.4f}, within 2 {e2:.3f}, within 3 {e3:.3f};"
          f" sign wrong at {int((~signs).sum())} nodes ({int(on.sum())} on the surface);"
          f" beyond 3: {ratio.min():.2f}x..{ratio.max():.2f}x  {'OK' if good else 'FAIL'}")
    return good


def check_lattice():
    print("=== lattice: the graph's node positions")
    nx, ny, nz, dx, origin = 9, 5, 7, 0.25, (-3.0, 1.0, 2.0)
    phi, pos = solid_phi_from_meshes(TAG + "lattice", [sphere_mesh((0, 2, 3), 0.5)], nx, ny, nz,
                                     cell_size=dx, origin=origin, positions=True)
    exp = np.stack(grid_world(nx, ny, nz, dx, origin), -1)
    diff = np.abs(pos - exp).max()
    ok = diff == 0.0 and phi.shape == (nx + 1, ny + 1, nz + 1)
    print(f"    {nx}x{ny}x{nz} nodes, cell {dx}, origin {origin}: max |position error| {diff}"
          f"  {'OK' if ok else 'FAIL'}")
    return ok


def check_field():
    print("=== field: meshes against exact signed distances (in cells)")
    ok = True
    # a sphere, world units: cell 0.5, origin off the axes
    nx, ny, nz, dx, origin = 40, 32, 36, 0.5, (-2.0, 0.5, 1.0)
    c, r = (8.0, 6.0, 9.0), 3.75
    x, y, z = grid_world(nx, ny, nz, dx, origin)
    exact = (np.sqrt((x - c[0]) ** 2 + (y - c[1]) ** 2 + (z - c[2]) ** 2) - r) / dx
    phi = solid_phi_from_meshes(TAG + "sphere", [sphere_mesh(c, r)], nx, ny, nz, dx, origin)
    ok &= compare(f"sphere r {r} at cell {dx}", phi, exact)
    # a box, world units, its faces between nodes
    lo, hi = (3.1, 0.6, 4.3), (9.3, 7.1, 11.8)
    exact_b = sdf_box3(x, y, z, lo + hi) / dx
    phi_b = solid_phi_from_meshes(TAG + "box", [box_mesh(lo, hi)], nx, ny, nz, dx, origin)
    ok &= compare("box off the nodes", phi_b, exact_b)
    # the same box with its faces on nodes: the degenerate case
    lo2, hi2 = (3.0, 0.5, 4.0), (9.0, 7.0, 12.0)
    exact_b2 = sdf_box3(x, y, z, lo2 + hi2) / dx
    phi_b2 = solid_phi_from_meshes(TAG + "box_on", [box_mesh(lo2, hi2)], nx, ny, nz, dx, origin)
    ok &= compare("box with faces on nodes", phi_b2, exact_b2)
    on = exact_b2 == 0.0
    print(f"      nodes on the box's faces read {phi_b2[on].min():+.4f}..{phi_b2[on].max():+.4f}"
          " cells (exact 0)")
    # two colliders: the minimum of the two fields
    both = solid_phi_from_meshes(TAG + "both", [sphere_mesh(c, r), box_mesh(lo, hi)], nx, ny, nz,
                                 dx, origin)
    dmin = np.abs(both - np.minimum(phi, phi_b)).max()
    good = dmin <= 1e-5
    ok &= good
    print(f"    sphere and box together: max |field - min(sphere, box)| {dmin:.2e}"
          f"  {'OK' if good else 'FAIL'}")
    return ok


def crossing(line, xs):
    """Where a sampled line crosses zero (linear), first crossing from the left."""
    s = np.sign(line)
    k = int(np.argmax(s[:-1] != s[1:]))
    t = line[k] / (line[k] - line[k + 1])
    return xs[k] + t * (xs[k + 1] - xs[k])


def check_settings():
    print("=== settings: the shipped node's knobs")
    ok = True
    nx, ny, nz, dx = 48, 24, 24, 0.5
    c, r = (12.0, 6.0, 6.0), 3.0
    x = np.arange(nx + 1) * dx
    j, k = int(c[1] / dx), int(c[2] / dx)
    base = solid_phi_from_meshes(TAG + "s_base", [sphere_mesh(c, r)], nx, ny, nz, dx)
    x0 = crossing(base[:, j, k], x)
    print(f"    sphere r {r} at x {c[0]}: surface crossing at x = {x0:.4f} (exact {c[0] - r})")
    cases = (("thickening 1 (voxel 0.5)", {"thickening": 1.0}, 1.0 * dx),
             ("thickening 2", {"thickening": 2.0}, 2.0 * dx),
             ("offset 0.75 (scene units)", {"offset": 0.75}, 0.75),
             ("detail size 2 cells", {"detail_size": 2.0 * dx}, 0.0),
             ("detail size half a cell", {"detail_size": 0.5 * dx}, 0.0))
    for name, kw, shift in cases:
        phi = solid_phi_from_meshes(TAG + "s_" + name.split()[0] + name.split()[1], [sphere_mesh(c, r)],
                                    nx, ny, nz, dx, **kw)
        xc = crossing(phi[:, j, k], x)
        tol = 0.5 * kw.get("detail_size", dx) if "detail" in name else 0.05 * dx
        good = abs((x0 - xc) - shift) <= tol and inside(phi, np.array([[c[0] / dx, j, k]])) == 1
        ok &= good
        print(f"    {name:28s} surface moved out by {x0 - xc:.4f} (expected {shift:.4f},"
              f" tolerance {tol:.3f})  {'OK' if good else 'FAIL'}")
    return ok


def check_step():
    """The stage-2 pillar, shifted a quarter cell and sunk a cell into the floor, so
    no node and no face centre sits on its surface: there the decision is a genuine
    toss-up, and +-1e-7 of level-set rounding flips a face the analytic field's exact
    zero closes."""
    print("=== step: the stage-2 pillar as a mesh, a quarter cell off the nodes")
    nx, ny, nz = 40, 24, 32
    box3 = (18.25, -1.0, 12.25, 24.25, 16.25, 20.25)
    analytic = f32(sdf_box3(*node_grid3(nx, ny, nz), box3))
    mesh = f32(solid_phi_from_meshes(TAG + "pillar", [box_mesh(box3[:3], box3[3:])], nx, ny, nz))
    print(f"    field: |mesh - analytic| within 1 cell of the surface"
          f" {np.abs(mesh - analytic)[np.abs(analytic) <= 1.0].max():.4f}, everywhere"
          f" {np.abs(mesh - analytic).max():.3f}; signs agree at every node:"
          f" {bool((np.sign(mesh) == np.sign(analytic)).all())}")
    ok = True
    arr = flow_inputs(nx, ny, nz, analytic)
    os.environ["PFFLIP_SOLID_DEBUG"] = "1"   # the node reports its push-out count
    try:
        for mode_name, mode, bar in (("binary", 1, 0.0), ("cut", 2, 0.05)):
            fa = node_fractions(TAG + f"pf_a_{mode}", analytic, mode)
            fm = node_fractions(TAG + f"pf_m_{mode}", mesh, mode)
            diff = max(np.abs(fa[c] - fm[c]).max() for c in fa)
            ndiff = sum(int((fa[c] != fm[c]).sum()) for c in fa)
            # the port is float32: a float64 array reads as empty and the node refuses it
            d_a, _ = node(TAG + f"step_a_{mode}", nx, ny, nz,
                          dict(arr, solid_phi=analytic.astype(np.float32).ravel()), solid_mode=mode)
            d_m, st = node(TAG + f"step_m_{mode}", nx, ny, nz,
                           dict(arr, solid_phi=mesh.astype(np.float32).ravel()), solid_mode=mode)
            gdiff = max(np.abs(load(d_a, n) - load(d_m, n)).max()
                        for n in ("out_u", "out_v", "out_w"))
            pen = inside(mesh, load(d_m, "out_positions"))
            pushed = st.get("pushed", -1)
            good = (diff <= bar and pen == 0 and st.get("solid_on", False) and pushed > 0
                    and (gdiff == 0.0 if bar == 0.0 else gdiff <= 0.05))
            ok &= good
            print(f"    {mode_name:7s} face areas: {ndiff} faces differ, max |diff| {diff:.4f}"
                  f" (bar {bar}); one step: grid velocities max |diff| {gdiff:.2e} cells/s;"
                  f" the flow cut {pushed} particles into the pillar, all pushed out, inside"
                  f" after the step {pen}  {'OK' if good else 'FAIL'}")
    finally:
        del os.environ["PFFLIP_SOLID_DEBUG"]
    return ok


def check_refuse():
    """A field the node cannot use must stop it, not leave it running without the
    solid: a float64 .npy reads as empty through the float32 port, a field of the
    wrong size would index the wrong grid."""
    print("=== refuse: a field of the wrong dtype or size stops the node")
    nx, ny, nz = 40, 24, 32
    box3 = (18.25, -1.0, 12.25, 24.25, 16.25, 20.25)
    phi = sdf_box3(*node_grid3(nx, ny, nz), box3)       # float64
    arr = flow_inputs(nx, ny, nz, f32(phi))
    field32 = phi.astype(np.float32).ravel()            # (f32() returns float64)
    cases = (("float64 field, reads as empty", phi.ravel(), True),
             ("float32 field, one node short", field32[:-1], True),
             ("float32 field, the right size", field32, False))
    ok = True
    for name, field, want in cases:
        try:
            _, st = node(TAG + "refuse", nx, ny, nz, dict(arr, solid_phi=field), raw=True,
                         solid_mode=2)
            refused, outcome = False, f"ran, {st['iterations_used']:.0f} iterations"
        except SystemExit as e:
            refused, outcome = "refused" in str(e), str(e)
        good = refused == want
        ok &= good
        print(f"    {name}: {outcome}  {'OK' if good else 'FAIL'}")
    return ok


def check_loop():
    print("=== loop: still water over a submerged ball from the graph")
    n, steps = 32, 100
    ball = solid_phi_from_meshes(TAG + "ball", [sphere_mesh((16.0, 9.0, 16.0), 6.0)], n, n, n)
    res = {name: rest_3d(TAG + f"ball{mode}", n, 22.0, steps, f32(ball), mode)
           for name, mode in (("tank, no solid", 0), ("mesh ball, binary", 1),
                              ("mesh ball, cut cells", 2))}
    base = res["tank, no solid"][0]
    ok = True
    print(f"  {n}^3, {steps} steps: liquid rms speed, run's peak speed, most inside in a frame")
    for name, (rms, peak, pen) in res.items():
        good = pen == 0 and (rms <= 1.5 * base or name.startswith("tank"))
        ok &= good
        print(f"    {name:22s} {rms:.4f}  {peak:.3f}  {pen}  {'OK' if good else 'FAIL'}")
    print(f"    bar: within 1.5x the plain tank's {base:.4f} (stage 2 with the analytic ball:"
          " 0.0126 / 0.0126 against 0.0124)")
    return ok


CHECKS = {"lattice": check_lattice, "field": check_field, "settings": check_settings,
          "step": check_step, "refuse": check_refuse, "loop": check_loop}


def main():
    names = sys.argv[1:] or list(CHECKS)
    ok = True
    for name in names:
        ok &= bool(CHECKS[name]())
    print("COLLIDER GRAPH:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
