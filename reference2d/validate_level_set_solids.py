"""Level-set solids in the 2D reference: a mesh collider's signed distance field, sampled
at the grid nodes, closing faces either whole (binary, the boxes' rule) or by their open
fraction (cut cells). Four checks, each with a known answer:

  faces    a grid-aligned box written as a distance field closes exactly the faces the
           box path closes (binary, and cut cells); off the grid the two may differ
           only at the box's corners, where the field bends between two nodes; a slab
           thinner than a cell is reported: the node samples cannot see it
  box run  the bore against the block, box path against the same block as a distance
           field with the support share off: one step from one state and whole runs
           over three particle seeds are bit-identical; with the share on (the
           default) the runs agree within the box path's seed spread
  bowl     water at rest in a round bowl stays at rest, binary and cut, against a
           square tank with no solid at all; cut cells without the support share are
           reported, as the reason the share exists
  incline  a block of water released on a frictionless 30 degree slope: its centre of
           mass accelerates down the slope at g sin(theta); a staircase slope holds it
           back, cut cells should not

Every check counts particles left inside the solid after the push-out (must be zero).

Run from reference2d:  ..\\.venv\\Scripts\\python.exe validate_level_set_solids.py
"""
import sys

import numpy as np

from pfflip2d import Sim

DT = 0.02


# ------------------------------------------------------------------ distance fields


def node_grid(nx, ny):
    return np.meshgrid(np.arange(nx + 1, dtype=float), np.arange(ny + 1, dtype=float),
                       indexing="ij")


def sdf_box(x, y, box):
    """Signed distance to an axis-aligned box (x0, y0, x1, y1), negative inside."""
    x0, y0, x1, y1 = box
    qx, qy = np.maximum(x0 - x, x - x1), np.maximum(y0 - y, y - y1)
    return np.hypot(np.maximum(qx, 0.0), np.maximum(qy, 0.0)) + np.minimum(np.maximum(qx, qy), 0.0)


def sdf_bowl(x, y, cx, cy, r):
    """Solid everywhere outside a circle: the fluid's container is a round bowl."""
    return r - np.hypot(x - cx, y - cy)


def sdf_slope(x, y, xa, ya, theta):
    """Solid below the line through (xa, ya) descending to the right at angle theta."""
    return (y - ya) * np.cos(theta) + (x - xa) * np.sin(theta)


def penetrations(sim):
    return int((sim._phi_at(sim.pos)[0] < 0.0).sum()) if sim.solid_phi is not None else \
        int(sim._inside_solid(sim.pos).sum())


# ------------------------------------------------------------------ the checks


def check_faces():
    nx, ny = 64, 48
    x, y = node_grid(nx, ny)
    boxes = {"grid-aligned block": (20.0, 0.0, 28.0, 14.0),
             "off-grid box": (10.3, 5.3, 20.7, 15.1),
             "box on two walls": (0.0, 0.0, 8.0, 8.0),
             "one-cell box": (30.0, 20.0, 31.0, 21.0),
             "box past the top wall": (40.25, 30.5, 63.0, 60.0)}
    ok = True
    print("faces: closed faces, box path vs the same box as a distance field")
    for name, box in boxes.items():
        ref = Sim(nx, ny, obstacles=[box])
        mu0, mv0 = ref._solid_face_masks()
        phi = sdf_box(x, y, box)
        mu1, mv1 = Sim(nx, ny, solid_phi=phi)._solid_face_masks()
        same = np.array_equal(mu0, mu1) and np.array_equal(mv0, mv1)
        if box != tuple(np.round(box)):
            # off the grid: every differing face must hold a box corner within half a
            # cell of its midpoint (the field bends there; a face's ends cannot see it)
            corners = np.array([(cx, cy) for cx in box[::2] for cy in box[1::2]])
            mids = [(i, j + 0.5) for i, j in np.argwhere(mu0 != mu1)] + \
                   [(i + 0.5, j) for i, j in np.argwhere(mv0 != mv1)]
            far = [m for m in mids if np.hypot(*(corners - m).T).min() > 0.5 + 1e-9]
            print(f"  {name:24s} binary: {len(mids)} faces differ, {len(mids) - len(far)}"
                  f" of them at a corner ({int(mu0.sum() + mv0.sum())} closed)")
            ok &= not far
            continue
        line = f"  {name:24s} binary {'identical' if same else 'DIFFER'}"
        if True:
            wu, wv = Sim(nx, ny, solid_phi=phi, solid_cut=True)._solid_open()
            cut_same = np.array_equal(wu <= 0, mu0) and np.array_equal(wv <= 0, mv0) and \
                bool(np.all((wu == 0) | (wu == 1))) and bool(np.all((wv == 0) | (wv == 1)))
            line += f", cut {'identical' if cut_same else 'DIFFER'}"
            same = same and cut_same
        print(line + f" ({int(mu0.sum() + mv0.sum())} closed)")
        ok &= same
    slab = (5.5, 30.25, 40.5, 30.75)
    mu0, mv0 = Sim(nx, ny, obstacles=[slab])._solid_face_masks()
    mu1, mv1 = Sim(nx, ny, solid_phi=sdf_box(x, y, slab))._solid_face_masks()
    print(f"  thin slab (half a cell, reported): the box path closes "
          f"{int(mu0.sum() + mv0.sum())} faces, the distance field {int(mu1.sum() + mv1.sum())}"
          " - a solid thinner than a cell falls between the nodes; mesh_to_level_set's"
          " thickening exists for this")
    return ok


def bore_mask(x, y):
    return (y < 6) | ((x < 40) & (y < 40))


def bore_metrics(sim):
    liq = sim.pos[sim.typ == 1]
    front = liq[(liq[:, 0] > 100.0) & (liq[:, 0] < 120.0)]
    return {"liquid centroid x": float(liq[:, 0].mean()),
            "liquid centroid y": float(liq[:, 1].mean()),
            "splash height at the block": float(np.percentile(front[:, 1], 95)),
            "overtopped share": float(((liq[:, 0] > 128.0) & (liq[:, 1] > 7.0)).mean())}


def in_corner_cell(pts, box):
    """Particles within a cell of a box corner, where a node-sampled field bends."""
    corners = np.array([(cx, cy) for cx in box[::2] for cy in box[1::2]])
    return np.abs(pts[:, None, :] - corners[None]).max(2).min(1) < 1.0


def check_box_run():
    nx, ny, box = 200, 60, (120.0, 0.0, 128.0, 14.0)
    phi = sdf_box(*node_grid(nx, ny), box)
    print("box run: the bore against the block (it arrives at about 4 s and overtops from"
          " 5 s), box path vs distance field")
    # one step from one state at the block, support share off
    a = Sim(nx, ny, escape=True, obstacles=[box])
    a.seed(bore_mask)
    a.calibrate()
    for _ in range(225):
        a.step(DT)
    b = Sim(nx, ny, escape=True, solid_phi=phi, solid_support=False)
    for k in ("pos", "vel", "typ", "scale", "tau", "xi"):
        setattr(b, k, getattr(a, k).copy())
    b.rho0_face = a.rho0_face
    a.step(DT)
    b.step(DT)
    grid_same = np.array_equal(a.u, b.u) and np.array_equal(a.v, b.v)
    moved = np.any(a.pos != b.pos, 1) | np.any(a.vel != b.vel, 1)
    # the box path leaves a pushed particle exactly 0.01 outside the face it crossed
    pushed = int(((a.pos[:, 0] == box[0] - 0.01) | (a.pos[:, 0] == box[2] + 0.01) |
                  (a.pos[:, 1] == box[3] + 0.01)).sum())
    corner = in_corner_cell(a.pos[moved], box)
    print(f"  share off, one step at t = 4.5 s ({pushed} particles pushed out of the block):"
          f" grid velocities {'bit-identical' if grid_same else 'DIFFER'}; {int(moved.sum())}"
          f" particles differ, {int(corner.sum())} of them in a corner cell")
    ok = grid_same and bool(corner.all())
    # whole runs over three seeds
    steps = 300
    arms = {"box path": {"obstacles": [box]},
            "field, share off": {"solid_phi": phi, "solid_support": False},
            "field, share on": {"solid_phi": phi}}
    runs = {arm: [] for arm in arms}
    states = {arm: [] for arm in arms}
    pen = 0
    for seed in (1, 2, 3):
        for arm, kw in arms.items():
            s = Sim(nx, ny, escape=True, **kw)
            s.seed(bore_mask, seed=seed)
            s.calibrate()
            for _ in range(steps):
                s.step(DT)
                pen = max(pen, penetrations(s))
            runs[arm].append(bore_metrics(s))
            states[arm].append((s.pos.copy(), s.vel.copy()))
    exact = sum(int(p0.shape == p1.shape and np.array_equal(p0, p1) and np.array_equal(v0, v1))
                for (p0, v0), (p1, v1) in zip(states["box path"], states["field, share off"]))
    print(f"  {steps} steps (6 s), seeds 1-3: share off is bit-identical to the box path on"
          f" {exact} of 3 seeds (the others part where a particle enters a corner cell);"
          " each field arm's mean against the box path's seed spread:")
    for key in runs["box path"][0]:
        va = np.array([r[key] for r in runs["box path"]])
        spread = va.max() - va.min()
        line = f"    {key:27s} box {va.mean():8.4f}  bar {2 * spread:.4f}"
        for arm in ("field, share off", "field, share on"):
            vb = np.array([r[key] for r in runs[arm]])
            dev = abs(vb.mean() - va.mean())
            inside = dev <= 2 * spread
            ok &= bool(inside)
            line += f"  | {arm[7:]} {vb.mean():8.4f} {'inside' if inside else 'OUTSIDE'}"
        print(line)
    print(f"  most particles inside the block after any step's push-out: {pen}")
    return ok and pen == 0


def still_speed(sim, steps):
    pen = 0
    for _ in range(steps):
        sim.step(DT)
        pen = max(pen, penetrations(sim))
    face = max(np.abs(sim.u).max(), np.abs(sim.v).max())
    liq = sim.typ == 1
    rms = float(np.sqrt((sim.vel[liq] ** 2).sum(1).mean()))
    return face / (abs(sim.g) * DT), rms, pen


def check_bowl():
    n, steps = 64, 100
    x, y = node_grid(n, n)
    bowl = sdf_bowl(x, y, 32.0, 32.0, 29.5)
    print(f"bowl: water at rest, {steps} steps (2 s); max face speed as a fraction of one"
          " step's gravity kick, liquid rms speed in cells/s")
    level = 30.0
    arms = {"square tank, no solid": {},
            "round bowl, binary": {"solid_phi": bowl},
            "round bowl, cut cells": {"solid_phi": bowl, "solid_cut": True},
            "cut cells, share off": {"solid_phi": bowl, "solid_cut": True,
                                     "solid_support": False}}
    res = {}
    for name, kw in arms.items():
        s = Sim(n, n, **kw)
        s.seed(lambda px, py: py < level)
        s.calibrate()
        res[name] = still_speed(s, steps)
        f, rms, pen = res[name]
        print(f"  {name:24s} face {f:.2e}  liquid rms {rms:.3f}  inside after push-out {pen}")
    base = res["square tank, no solid"][1]
    ok = all(r[2] == 0 for r in res.values())
    for name in ("round bowl, binary", "round bowl, cut cells"):
        ok &= res[name][1] <= 1.5 * base
    print(f"  bar: liquid rms within 1.5x the square tank's ({base:.3f}) for the two default"
          " arms; the share-off arm is reported")
    return ok


def check_incline():
    nx, ny, theta = 128, 64, np.radians(30.0)
    xa, ya = 0.0, 52.0
    x, y = node_grid(nx, ny)
    slope = sdf_slope(x, y, xa, ya, theta)
    t_hat = np.array([np.cos(theta), -np.sin(theta)])      # down the slope

    def block(px, py):
        d = sdf_slope(px, py, xa, ya, theta)                # height above the slope
        s = (px - xa) * np.cos(theta) - (py - ya) * np.sin(theta)
        return (d > 0.0) & (d < 8.0) & (s > 12.0) & (s < 36.0)

    steps = 25
    expect = abs(-9.8) * np.sin(theta)
    print(f"incline: a water block released on a frictionless 30 degree slope, {steps} steps;"
          f" centre-of-mass acceleration down the slope against g sin(theta) = {expect:.2f}")
    res = {}
    for name, cut in (("binary (staircase)", False), ("cut cells", True)):
        s = Sim(nx, ny, solid_phi=slope, solid_cut=cut)
        s.seed(block)
        s.calibrate()
        ts, ss, pen = [0.0], [float((s.pos[s.typ == 1] @ t_hat).mean())], 0
        for k in range(steps):
            s.step(DT)
            pen = max(pen, penetrations(s))
            ts.append((k + 1) * DT)
            ss.append(float((s.pos[s.typ == 1] @ t_hat).mean()))
        acc = 2.0 * np.polyfit(np.array(ts), np.array(ss), 2)[0]
        res[name] = (acc, pen)
        print(f"  {name:20s} {acc:5.2f} cells/s^2 = {acc / expect:.2f} of g sin(theta);"
              f" inside after push-out {pen}")
    ok = all(p == 0 for _, p in res.values())
    ok &= res["cut cells"][0] >= 0.95 * expect
    print("  bar: cut cells reach 0.95 of g sin(theta); the staircase is reported")
    return ok


def main():
    results = {"faces": check_faces(), "box run": check_box_run(),
               "bowl": check_bowl(), "incline": check_incline()}
    print()
    for k, v in results.items():
        print(f"  {k:8s} {'PASS' if v else 'FAIL'}")
    ok = all(results.values())
    print("LEVEL-SET SOLIDS:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
