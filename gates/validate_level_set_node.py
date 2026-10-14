"""Level-set solids in step_3d, stage 2 of the mesh colliders: the node against the 2D
reference at nz = 1, against its own box path, and in 3D. The 2D record is
reference2d/validate_level_set_solids.py; this gate holds the node to it.

  faces  the node's open fractions, dumped through PFFLIP_SOLID_DUMP: equal to the
         reference's at nz = 1 (binary and cut cells); the same in every slab of a
         z-invariant field at nz = 4; and, with cut cells, exact for planes at any
         tilt on all three face kinds
  step   one step from one state. A grid-aligned box written as a field gives the
         box path's grid velocities and pressure bit for bit (binary and cut cells),
         at nz = 1 on the bore against the block and in 3D on a pillar; particles may
         differ only beside the box's corners and edges, where the field bends. And
         the node against the reference's twin at nz = 1: binary faces on the bore,
         cut cells on a slope
  phase  why the node needs no support share (the reference's fix for cut cells): its
         phase is the mean of the particle phases a face's kernel reaches, so a face
         the solid half covers still reads liquid, where the reference's raw-mass
         phase reads it as half air (python only)
  loops  the node's own loop (the fused graph). At nz = 1: water at rest in a round
         bowl stays at rest (binary and cut cells, against a square tank); a block
         slides down a 30 degree slope at g sin(theta) with cut cells; the bore
         against the block as a field matches the box path within its seed spread.
         In 3D: still water over a submerged ball stays as still as without it; in a
         half-full spherical bowl, where the surface meets the curved wall all round,
         cut cells raise no speed spikes against binary faces (the sliver cutoff's
         check); and a block slides on a slope tilted off every axis
Every loop counts particles inside the solid in its dumped frames (must be zero).

Run from the repo root:  .venv\\Scripts\\python.exe gates\\validate_level_set_node.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "reference2d"))
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
from pfflip2d import Sim  # noqa: E402
from validate_level_set_solids import (bore_mask, bore_metrics, in_corner_cell,  # noqa: E402
                                       node_grid, sdf_bowl, sdf_box, sdf_slope)
from validate_step import channels  # noqa: E402
from validate_3d import run_node  # noqa: E402
from loop3d import run_bifrost3, seed3d  # noqa: E402

OUT = os.path.join(HERE, "..", "out", "level_set_node")
DT = 0.02
RHO_L, RHO_G, G = 1000.0, 1.0, -9.8
MODES = (("binary", 1), ("cut", 2))


# ------------------------------------------------------------------ fields and helpers


def f32(a):
    """What the node sees: the field rounded to the float32 port."""
    return np.asarray(a, np.float32).astype(np.float64)


def extrude(phi2, nz):
    return np.repeat(np.asarray(phi2)[:, :, None], nz + 1, axis=2)


def node_grid3(nx, ny, nz):
    return np.meshgrid(np.arange(nx + 1, dtype=float), np.arange(ny + 1, dtype=float),
                       np.arange(nz + 1, dtype=float), indexing="ij")


def sdf_box3(x, y, z, box):
    """Signed distance to a box (x0, y0, z0, x1, y1, z1), negative inside."""
    lo, hi = np.array(box[:3]), np.array(box[3:])
    q = np.stack([np.maximum(lo[a] - c, c - hi[a]) for a, c in enumerate((x, y, z))])
    outside = np.sqrt((np.maximum(q, 0.0) ** 2).sum(0))
    return outside + np.minimum(q.max(0), 0.0)


def phi_at(phi, pts):
    """The node's trilinear distance at points (clamped like the node)."""
    n = np.array(phi.shape) - 1
    q = np.clip(pts, 0.0, n)
    b = np.minimum(np.floor(q).astype(np.int64), n - 1)
    f = q - b
    out = np.zeros(len(pts))
    for a in (0, 1):
        for c in (0, 1):
            for e in (0, 1):
                w = ((f[:, 0] if a else 1 - f[:, 0]) * (f[:, 1] if c else 1 - f[:, 1]) *
                     (f[:, 2] if e else 1 - f[:, 2]))
                out += w * phi[b[:, 0] + a, b[:, 1] + c, b[:, 2] + e]
    return out


def inside(phi, pos):
    return int((phi_at(phi, np.asarray(pos, np.float64)) < 0.0).sum())


def sizes(nx, ny, nz):
    return {"u": (nx + 1) * ny * nz, "v": nx * (ny + 1) * nz, "w": nx * ny * (nz + 1)}


def shapes(nx, ny, nz):
    return {"u": (nx + 1, ny, nz), "v": (nx, ny + 1, nz), "w": (nx, ny, nz + 1)}


def save(d, name, arr):
    os.makedirs(d, exist_ok=True)
    np.save(os.path.join(d, name + ".npy"), arr)


def load(d, name):
    return np.load(os.path.join(d, name + ".npy"))


def node(tag, nx, ny, nz, arrays, raw=False, **ports):
    """One step_3d call: arrays (name -> array) are written as the node's inputs."""
    d_in, d_out = os.path.join(OUT, tag + "_in"), os.path.join(OUT, tag + "_out")
    for k, v in arrays.items():
        # the field's port is float32; any other dtype reads as empty, which the
        # node refuses (raw=True sends the array as it is, to test that)
        save(d_in, k, np.asarray(v, np.float32) if k == "solid_phi" and not raw else v)
    extra = tuple(k for k in arrays if k in ("obstacles_min", "obstacles_max", "solid_phi"))
    stats = run_node(d_in, d_out, extra_arrays=extra, nx=nx, ny=ny, nz=nz, dt=DT,
                     preconditioner=1, **ports)
    return d_out, stats


def flat_inputs(sim, ch):
    """A 2D reference state and its channels, packaged for the node at nz = 1."""
    n = len(sim.pos)
    nw = sim.nx * sim.ny * 2
    arr = {k: v.astype(np.float32).ravel() for k, v in ch.items()}
    arr.update({"w_mass": np.zeros(nw, np.float32), "w_mom": np.zeros(nw, np.float32),
                "w_phase": np.zeros(nw, np.float32),
                "positions": np.column_stack([sim.pos, np.full(n, 0.5)]).astype(np.float32),
                "velocities": np.column_stack([sim.vel, np.zeros(n)]).astype(np.float32),
                "particle_phase": sim.typ.astype(np.float32)})
    return arr


# ------------------------------------------------------------------ faces


def node_fractions(tag, phi, mode):
    nx, ny, nz = (int(s) - 1 for s in phi.shape)
    sz = sizes(nx, ny, nz)
    arr = {f"{c}_{k}": np.zeros(n, np.float32) for c, n in sz.items()
           for k in ("mass", "mom", "phase")}
    arr.update({"positions": np.array([[0.5 * nx, 0.5 * ny, 0.5 * nz]], np.float32),
                "velocities": np.zeros((1, 3), np.float32),
                "particle_phase": np.zeros(1, np.float32),
                "solid_phi": phi.astype(np.float32).ravel()})
    dump = os.path.join(OUT, tag + "_open.bin")
    if os.path.exists(dump):
        os.remove(dump)
    os.environ["PFFLIP_SOLID_DUMP"] = dump
    try:
        node(tag, nx, ny, nz, arr, solid_mode=mode)
    finally:
        del os.environ["PFFLIP_SOLID_DUMP"]
    raw = np.fromfile(dump, np.float64)
    out, at = {}, 0
    for c, shp in shapes(nx, ny, nz).items():
        n = sz[c]
        out[c] = raw[at:at + n].reshape(shp)
        at += n
    assert at == len(raw), (at, len(raw))
    return out


def clip_open(a, b, c, d):
    """Exact open area of a unit square whose distance is linear: the square,
    corners (0,0) a, (1,0) b, (1,1) d, (0,1) c, clipped to the positive side."""
    pts = ((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0))
    val = (a, b, d, c)
    poly = []
    for k in range(4):
        p, q, vp, vq = pts[k], pts[(k + 1) % 4], val[k], val[(k + 1) % 4]
        if vp > 0:
            poly.append(p)
        if (vp > 0) != (vq > 0):
            t = vp / (vp - vq)
            poly.append((p[0] + t * (q[0] - p[0]), p[1] + t * (q[1] - p[1])))
    if len(poly) < 3:
        return 0.0
    xs, ys = np.array(poly).T
    return 0.5 * abs(float(np.dot(xs, np.roll(ys, -1)) - np.dot(ys, np.roll(xs, -1))))


def face_corners(phi, comp, i, j, k):
    if comp == "u":
        return phi[i, j, k], phi[i, j + 1, k], phi[i, j, k + 1], phi[i, j + 1, k + 1]
    if comp == "v":
        return phi[i, j, k], phi[i + 1, j, k], phi[i, j, k + 1], phi[i + 1, j, k + 1]
    return phi[i, j, k], phi[i + 1, j, k], phi[i, j + 1, k], phi[i + 1, j + 1, k]


def check_faces():
    print("=== faces: the node's open fractions")
    ok = True
    nx, ny = 64, 48
    x, y = node_grid(nx, ny)
    fields = {"grid-aligned block": sdf_box(x, y, (20.0, 0.0, 28.0, 14.0)),
              "off-grid box": sdf_box(x, y, (10.3, 5.3, 20.7, 15.1)),
              "round bowl": sdf_bowl(x, y, 32.0, 24.0, 21.5),
              "30 degree slope": sdf_slope(x, y, 0.0, 40.0, np.radians(30.0))}
    print("  nz = 1 against the reference's faces; nz = 4 (the field the same in every"
          " slab) against nz = 1")
    for name, phi2 in fields.items():
        line = f"    {name:20s}"
        for mode_name, mode in MODES:
            ref = Sim(nx, ny, solid_phi=f32(phi2), solid_cut=mode == 2)._solid_open()
            flat = node_fractions(f"faces_flat_{mode}", extrude(phi2, 1), mode)
            same = (np.array_equal(flat["u"][:, :, 0], ref[0]) and
                    np.array_equal(flat["v"][:, :, 0], ref[1]))
            deep = node_fractions(f"faces_deep_{mode}", extrude(phi2, 4), mode)
            dz = max(float(np.abs(deep[c] - flat[c][:, :, :1]).max()) for c in ("u", "v"))
            good = same and (dz == 0.0 if mode == 1 else dz <= 1e-12)
            ok &= good
            line += (f"  {mode_name}: {'identical' if same else 'DIFFER'} to the reference,"
                     f" slabs within {dz:.1e}")
        print(line)
    # planes: dyadic coefficients keep the nodal values exact in float32 and exactly
    # planar, so the exact clipped area is the answer to the last bit's worth
    nx, ny, nz = 12, 10, 8
    x, y, z = node_grid3(nx, ny, nz)
    planes = {"tilted (0.5, 0.75, -0.25)": 0.5 * x + 0.75 * y - 0.25 * z - 3.125,
              "tilted (-0.375, 0.5, 0.625)": -0.375 * x + 0.5 * y + 0.625 * z - 2.0625,
              "tilted (0.25, -0.5, 0.875)": 0.25 * x - 0.5 * y + 0.875 * z - 1.1875,
              "axis-aligned x = 5.3125": x - 5.3125}
    print("  cut cells on planes, every face kind against the exact clipped area;"
          " binary against the face centre's sign")
    for name, phi in planes.items():
        cut = node_fractions("faces_plane_cut", phi, 2)
        binary = node_fractions("faces_plane_bin", phi, 1)
        err, wrong, partial = 0.0, 0, 0
        for comp, shp in shapes(nx, ny, nz).items():
            for i, j, k in np.ndindex(*shp):
                a, b, c, d = face_corners(phi, comp, i, j, k)
                exact = clip_open(a, b, c, d)
                exact = 0.0 if exact < Sim.SOLID_MIN_OPEN else exact   # slivers close
                partial += 0.0 < exact < 1.0
                err = max(err, abs(cut[comp][i, j, k] - exact))
                wrong += binary[comp][i, j, k] != float(0.25 * ((a + b) + (c + d)) > 0)
        good = err <= 1e-12 and wrong == 0
        ok &= good
        print(f"    {name:28s} cut: max error {err:.1e} over {partial} partly open faces;"
              f" binary: {wrong} wrong")
    return ok


# ------------------------------------------------------------------ one step


def twin_step(state, ch, kw, rho0_face, escape):
    """The reference's step on the node's inputs: channels in place of its own
    splat, their phase in place of Eq. 7 (the node's phase is the channel's)."""
    ref = Sim(state.nx, state.ny, escape=escape, **kw)
    arr = flat_inputs(state, ch)
    ref.pos = arr["positions"][:, :2].astype(float)
    ref.vel = arr["velocities"][:, :2].astype(float)
    ref.typ = state.typ.copy()
    ref.rho0_face = rho0_face
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
    esc = None
    if escape:
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
    if esc is not None and esc.any():
        vg = ref.sample_velocity(ref.pos)
        ref.vel[drop, 1] += ref.g * DT
        ref.vel[drop] += ref.drag_droplet * DT * (vg[drop] - ref.vel[drop])
        ref.vel[bub, 1] += ref.buoyancy * abs(ref.g) * DT
        ref.vel[bub] += ref.drag_bubble * DT * (vg[bub] - ref.vel[bub])
    ref.advect(DT, ballistic=esc)
    ref._push_out()
    return ref, esc


def lockstep(name, d_out, ref, esc, phi3):
    nx, ny = ref.nx, ref.ny
    pos = load(d_out, "out_positions")[:, :2]
    vel = load(d_out, "out_velocities")[:, :2]
    du = max(float(np.abs(load(d_out, "out_u").reshape(nx + 1, ny) - ref.u).max()),
             float(np.abs(load(d_out, "out_v").reshape(nx, ny + 1) - ref.v).max()))
    dp = float(np.abs(pos - ref.pos).max())
    dv = float(np.abs(vel - ref.vel).max())
    mism = int((load(d_out, "out_escaped") > 0.5).astype(bool).__ne__(esc).sum()) \
        if esc is not None else 0
    pen_n = inside(phi3, np.column_stack([pos, np.full(len(pos), 0.5)]))
    pen_r = int((ref._phi_at(ref.pos)[0] < 0.0).sum())
    good = pen_n == 0 and pen_r == 0 and mism <= 2 and dp < 2e-3 and dv < 5e-2 and du < 1e-3
    print(f"    {name:34s} grid {du:.1e}, positions {dp:.1e}, velocities {dv:.1e},"
          f" escape mismatches {mism}; inside: node {pen_n}, reference {pen_r}"
          f"  {'OK' if good else 'FAIL'}")
    return good


def near_box_edge(pts, box):
    """Within a cell of a box edge (3D) or corner (2D box, z ignored): where the
    node-sampled field bends and the box path's least-penetration push differs."""
    lo, hi = np.array(box[:3]), np.array(box[3:])
    hit = np.zeros(len(pts), bool)
    for a, b in ((0, 1), (0, 2), (1, 2)):
        c = 3 - a - b
        along = (pts[:, c] > lo[c] - 1.0) & (pts[:, c] < hi[c] + 1.0)
        for ea in (lo[a], hi[a]):
            for eb in (lo[b], hi[b]):
                hit |= along & (np.abs(pts[:, a] - ea) < 1.0) & (np.abs(pts[:, b] - eb) < 1.0)
    return hit


def box_against_field(tag, nx, ny, nz, arr, box3, phi3, ports, corner_fn):
    """The node with the box path and with the same box as a field, binary and cut:
    grid velocities and pressure must be bit-identical; particles may differ only
    where corner_fn says the field bends."""
    lo = np.array([box3[:3]], np.float32)
    hi = np.array([box3[3:]], np.float32)
    d_box, _ = node(tag + "_box", nx, ny, nz, dict(arr, obstacles_min=lo, obstacles_max=hi),
                    **ports)
    ok = True
    pos0 = load(d_box, "out_positions")
    for mode_name, mode in MODES:
        d_f, _ = node(f"{tag}_{mode_name}", nx, ny, nz,
                      dict(arr, solid_phi=phi3.astype(np.float32).ravel()),
                      solid_mode=mode, **ports)
        grids = all(np.array_equal(load(d_box, n), load(d_f, n))
                    for n in ("out_u", "out_v", "out_w", "pressure"))
        pos1 = load(d_f, "out_positions")
        moved = (np.any(pos0 != pos1, 1) |
                 np.any(load(d_box, "out_velocities") != load(d_f, "out_velocities"), 1))
        at = corner_fn(pos0[moved])
        good = grids and bool(at.all()) and inside(phi3, pos1) == 0
        ok &= good
        print(f"    {mode_name:7s} grid velocities and pressure"
              f" {'bit-identical' if grids else 'DIFFER'}; {int(moved.sum())} particles differ,"
              f" {int(at.sum())} of them beside a corner or edge; inside after the step"
              f" {inside(phi3, pos1)}  {'OK' if good else 'FAIL'}")
    return ok


def check_step():
    print("=== step: one step from one state")
    ok = True
    # the bore against the block at nz = 1, t = 4.5 s (it arrives at about 4 s)
    nx, ny, box = 200, 60, (120.0, 0.0, 128.0, 14.0)
    phi2 = f32(sdf_box(*node_grid(nx, ny), box))
    a = Sim(nx, ny, escape=True, obstacles=[box])
    a.seed(bore_mask)
    a.calibrate()
    for _ in range(225):
        a.step(DT)
    ch = channels(a)
    arr = flat_inputs(a, ch)
    ports = {"escape": 1, "rho0_face": a.rho0_face}
    print(f"  bore at the block, nz = 1, t = 4.5 s: box path against the block as a field")
    ok &= box_against_field("bore", nx, ny, 1, arr, (box[0], box[1], -1.0, box[2], box[3], 2.0),
                            extrude(phi2, 1), ports,
                            lambda p: in_corner_cell(p[:, :2], box))
    print("  the node against the reference's twin, nz = 1")
    d_f, _ = node("bore_lock", nx, ny, 1, dict(arr, solid_phi=extrude(phi2, 1).astype(
        np.float32).ravel()), solid_mode=1, **ports)
    ref, esc = twin_step(a, ch, {"solid_phi": phi2, "solid_support": False},
                         a.rho0_face, True)
    ok &= lockstep("bore, binary faces", d_f, ref, esc, extrude(phi2, 1))
    # cut cells on a 30 degree slope, a block sliding for 10 steps
    nx, ny, theta = 128, 64, np.radians(30.0)
    slope = f32(sdf_slope(*node_grid(nx, ny), 0.0, 52.0, theta))
    s = Sim(nx, ny, solid_phi=slope, solid_cut=True, solid_support=False)

    def block(px, py):
        d = sdf_slope(px, py, 0.0, 52.0, theta)
        t = px * np.cos(theta) - (py - 52.0) * np.sin(theta)
        return (d > 0.0) & (d < 8.0) & (t > 12.0) & (t < 36.0)
    s.seed(block)
    s.calibrate()
    for _ in range(10):
        s.step(DT)
    ch = channels(s)
    d_f, _ = node("slope_lock", nx, ny, 1, dict(flat_inputs(s, ch), solid_phi=extrude(
        slope, 1).astype(np.float32).ravel()), solid_mode=2)
    ref, esc = twin_step(s, ch, {"solid_phi": slope, "solid_cut": True,
                                 "solid_support": False}, s.rho0_face, False)
    ok &= lockstep("slope, cut cells", d_f, ref, esc, extrude(slope, 1))
    # 3D: a pillar on the floor, liquid flowing at it
    nx, ny, nz = 40, 24, 32
    box3 = (18.0, 0.0, 12.0, 24.0, 16.0, 20.0)
    phi3 = sdf_box3(*node_grid3(nx, ny, nz), box3)
    print(f"  3D, {nx}x{ny}x{nz}: liquid at 25 cells/s into a pillar, box path against the"
          " pillar as a field")
    ok &= box_against_field("pillar", nx, ny, nz, flow_inputs(nx, ny, nz, phi3), box3, phi3,
                            {}, lambda p: near_box_edge(p, box3))
    return ok


def flow_inputs(nx, ny, nz, phi3, speed=25.0, level=10.0):
    """Analytic channels for liquid below level flowing at +x, and seeded particles
    outside the solid carrying the same flow."""
    arr = {}
    offs = {"u": (0.0, 0.5, 0.5), "v": (0.5, 0.0, 0.5), "w": (0.5, 0.5, 0.0)}
    for comp, shp in shapes(nx, ny, nz).items():
        x, y, z = np.meshgrid(*(np.arange(n) + o for n, o in zip(shp, offs[comp])),
                              indexing="ij")
        ph = (y < level).astype(float)
        mass = RHO_G + (RHO_L - RHO_G) * ph
        vel = speed * ph if comp == "u" else 0.0 * ph
        arr[comp + "_mass"] = mass.ravel().astype(np.float32)
        arr[comp + "_mom"] = (mass * vel).ravel().astype(np.float32)
        arr[comp + "_phase"] = ph.ravel().astype(np.float32)
    pos, typ = seed3d(nx, ny, nz, lambda x, y, z: y < level)
    keep = phi_at(phi3, pos) >= 0.0
    pos, typ = pos[keep], typ[keep]
    vel = np.zeros_like(pos)
    vel[:, 0] = speed * (typ == 1)
    arr.update({"positions": pos.astype(np.float32), "velocities": vel.astype(np.float32),
                "particle_phase": typ.astype(np.float32)})
    return arr


# ------------------------------------------------------------------ phase


def check_phase():
    print("=== phase: faces the solid partly covers, in still water in a round bowl")
    n, level = 64, 30.0
    bowl = sdf_bowl(*node_grid(n, n), 32.0, 32.0, 29.5)
    s = Sim(n, n, solid_phi=bowl, solid_cut=True)
    s.seed(lambda px, py: py < level)
    s.calibrate()
    s.p2g()
    wu, wv = s._solid_open()
    su, sv = s._solid_support_share()
    ch = channels(s)
    rows = []
    for w, mass, share, mean, off in ((wu, s.mu, su, ch["u_phase"], (0.0, 0.5)),
                                      (wv, s.mv, sv, ch["v_phase"], (0.5, 0.0))):
        jj = np.arange(w.shape[1]) + off[1]
        sub = np.broadcast_to(jj[None, :] < level - 3.0, w.shape)
        sel = (w > 0.0) & (w < 1.0) & sub & (mass > 0.0)
        rows.append((s.phase(mass)[sel], s.phase(mass / np.maximum(share, 0.1))[sel],
                     mean[sel]))
    raw, fixed, mean = (np.concatenate(c) for c in zip(*rows))
    print(f"  {len(raw)} submerged faces partly open; share reading below 0.9:")
    for name, v in (("reference, raw mass (Eq. 7)", raw),
                    ("reference, raw mass over the support share", fixed),
                    ("node, mean particle phase", mean)):
        print(f"    {name:44s} {int((v < 0.9).sum()):4d}   lowest {v.min():.3f}")
    ok = bool((mean >= 0.9).all())
    print("  bar: the node's phase reads every such face as liquid (>= 0.9)")
    return ok


# ------------------------------------------------------------------ loops


def loop(tag, nx, ny, nz, pos, typ, steps, solid=None, ports=None, every=1, seed=1):
    d, _, _ = run_bifrost3(f"lsn_{tag}", nx, ny, nz, None, steps, DT,
                           extra_ports=dict(ports or {}), seed_fn=lambda *_: (pos, typ),
                           fused=True, dump_every=every, final_state=True, solid=solid)
    frames = [load(d, f"frame.{k:04d}") for k in range(0, steps, every)]
    return d, frames


def lift(sim):
    return (np.column_stack([sim.pos, np.full(len(sim.pos), 0.5)]).astype(np.float32),
            sim.typ.astype(np.int32))


def rest_rms(d, typ):
    vel = load(d, "final_velocities")
    return float(np.sqrt((vel[typ == 1] ** 2).sum(1).mean()))


def check_loops_flat():
    print("=== loops at nz = 1 (the fused graph)")
    ok = True
    # still water in a round bowl
    n, level, steps = 64, 30.0, 100
    bowl = f32(sdf_bowl(*node_grid(n, n), 32.0, 32.0, 29.5))
    res = {}
    for name, mode in (("square tank, no solid", 0), ("round bowl, binary", 1),
                       ("round bowl, cut cells", 2)):
        s = Sim(n, n, solid_phi=bowl if mode else None)
        s.seed(lambda px, py: py < level)
        pos, typ = lift(s)
        solid = (extrude(bowl, 1), mode) if mode else None
        d, frames = loop(f"bowl{mode}", n, n, 1, pos, typ, steps, solid=solid, every=10)
        pen = max(inside(extrude(bowl, 1), f) for f in frames) if mode else 0
        res[name] = (rest_rms(d, typ), pen)
    base = res["square tank, no solid"][0]
    print(f"  bowl at rest, {steps} steps: liquid rms speed (cells/s), most inside in a frame")
    for name, (rms, pen) in res.items():
        good = pen == 0 and (rms <= 1.5 * base or name.startswith("square"))
        ok &= good
        print(f"    {name:24s} {rms:.3f}  {pen}  {'OK' if good else 'FAIL'}")
    print(f"    bar: within 1.5x the square tank's {base:.3f} (the reference: tank 0.171,"
          " binary 0.177, cut 0.192)")
    # a block on a 30 degree slope
    nx, ny, theta, steps = 128, 64, np.radians(30.0), 25
    slope = f32(sdf_slope(*node_grid(nx, ny), 0.0, 52.0, theta))
    t_hat = np.array([np.cos(theta), -np.sin(theta)])
    expect = abs(G) * np.sin(theta)

    def block(px, py):
        dd = sdf_slope(px, py, 0.0, 52.0, theta)
        t = px * np.cos(theta) - (py - 52.0) * np.sin(theta)
        return (dd > 0.0) & (dd < 8.0) & (t > 12.0) & (t < 36.0)
    print(f"  slope, {steps} steps: centre-of-mass acceleration against g sin(theta) ="
          f" {expect:.2f} (the reference: cut 0.99, staircase 0.94)")
    for name, mode in (("binary (staircase)", 1), ("cut cells", 2)):
        s = Sim(nx, ny, solid_phi=slope)
        s.seed(block)
        pos, typ = lift(s)
        d, frames = loop(f"slope{mode}", nx, ny, 1, pos, typ, steps + 1,
                         solid=(extrude(slope, 1), mode))
        ss = [float((f[typ == 1][:, :2] @ t_hat).mean()) for f in frames]
        acc = 2.0 * np.polyfit(DT * np.arange(len(ss)), np.array(ss), 2)[0]
        pen = max(inside(extrude(slope, 1), f) for f in frames)
        good = pen == 0 and (acc >= 0.95 * expect or mode == 1)
        ok &= good
        print(f"    {name:20s} {acc / expect:.2f} of g sin(theta); most inside {pen}"
              f"  {'OK' if good else 'FAIL'}")
    # the bore against the block, three seeds
    nx, ny, box, steps = 200, 60, (120.0, 0.0, 128.0, 14.0), 300
    phi2 = f32(sdf_box(*node_grid(nx, ny), box))
    print(f"  bore against the block, {steps} steps (6 s), seeds 1-3: the field arm's mean"
          " against the box path's seed spread")
    runs = {"box path": [], "field": []}
    pen = 0
    for seed in (1, 2, 3):
        for arm in runs:
            s = Sim(nx, ny, escape=True, obstacles=[box])
            s.seed(bore_mask, seed=seed)
            s.calibrate()
            pos, typ = lift(s)
            ports = {"escape": 1, "rho0_face": s.rho0_face}
            kw = ({"solid": (extrude(phi2, 1), 1)} if arm == "field" else {})
            if arm == "box path":
                d, _, _ = run_bifrost3(f"lsn_bore_box_s{seed}", nx, ny, 1, None, steps, DT,
                                       extra_ports=ports, seed_fn=lambda *_: (pos, typ),
                                       fused=True, dump_every=10, final_state=True,
                                       obstacles=(np.array([[box[0], box[1], -1.0]]),
                                                  np.array([[box[2], box[3], 2.0]])))
                frames = [load(d, f"frame.{k:04d}") for k in range(0, steps, 10)]
            else:
                d, frames = loop(f"bore_field_s{seed}", nx, ny, 1, pos, typ, steps,
                                 ports=ports, every=10, **kw)
            pen = max([pen] + [inside(extrude(phi2, 1), f) for f in frames])
            end = Sim(nx, ny)
            end.pos = load(d, "final")[:, :2].astype(float)
            end.typ = load(d, "final_phase").astype(np.int32)
            runs[arm].append(bore_metrics(end))
    for key in runs["box path"][0]:
        va = np.array([r[key] for r in runs["box path"]])
        vb = np.array([r[key] for r in runs["field"]])
        spread = va.max() - va.min()
        good = abs(vb.mean() - va.mean()) <= 2 * spread
        ok &= good
        print(f"    {key:27s} box {va.mean():8.4f}  field {vb.mean():8.4f}"
              f"  bar {2 * spread:.4f}  {'inside' if good else 'OUTSIDE'}")
    ok &= pen == 0
    print(f"    most particles inside the block in a dumped frame: {pen}")
    return ok


def rest_3d(tag, n, level, steps, phi, mode):
    pos, typ = seed3d(n, n, n, lambda px, py, pz: py < level)
    if mode:
        keep = phi_at(phi, pos) >= 0.0
        pos, typ = pos[keep], typ[keep]
    d, frames = loop(tag, n, n, n, pos.astype(np.float32), typ, steps,
                     solid=(phi, mode) if mode else None, every=10)
    peak = float(max(load(d, f"diag.{k:04d}")[-1] for k in range(steps)))
    pen = max(inside(phi, f) for f in frames) if mode else 0
    return rest_rms(d, typ), peak, pen


def check_loops_3d():
    print("=== loops in 3D")
    ok = True
    n, steps = 32, 100
    x, y, z = node_grid3(n, n, n)
    r = np.sqrt((x - 16.0) ** 2 + (y - 16.0) ** 2 + (z - 16.0) ** 2)
    # still water over a submerged ball: the same tank and free surface with and
    # without the solid, so the solid's own disturbance is all that differs
    ball = f32(np.sqrt((x - 16.0) ** 2 + (y - 9.0) ** 2 + (z - 16.0) ** 2) - 6.0)
    res = {name: rest_3d(f"ball{mode}", n, 22.0, steps, ball, mode)
           for name, mode in (("tank, no solid", 0), ("tank and ball, binary", 1),
                              ("tank and ball, cut cells", 2))}
    base = res["tank, no solid"][0]
    print(f"  still water over a submerged ball, {n}^3, {steps} steps: liquid rms speed,"
          " run's peak speed, most inside in a frame")
    for name, (rms, peak, pen) in res.items():
        good = pen == 0 and (rms <= 1.5 * base or name.startswith("tank,"))
        ok &= good
        print(f"    {name:26s} {rms:.4f}  {peak:.3f}  {pen}  {'OK' if good else 'FAIL'}")
    print(f"    bar: within 1.5x the plain tank's {base:.4f}")
    # a half-full spherical bowl: the water surface meets the curved wall all around.
    # Its rms is not comparable to a box tank's (more of its water lies near the noisy
    # surface), so cut cells answer to binary faces on the same bowl: no spikes
    bowl = f32(14.5 - r)
    res = {name: rest_3d(f"sphere{mode}", n, 16.0, steps, bowl, mode)
           for name, mode in (("spherical bowl, binary", 1), ("spherical bowl, cut cells", 2))}
    peak_bin = res["spherical bowl, binary"][1]
    print(f"  a half-full spherical bowl, {steps} steps: liquid rms speed, run's peak speed,"
          " most inside")
    for name, (rms, peak, pen) in res.items():
        good = pen == 0 and peak <= 2.0 * peak_bin
        ok &= good
        print(f"    {name:26s} {rms:.4f}  {peak:.3f}  {pen}  {'OK' if good else 'FAIL'}")
    print(f"    bar: the run's peak speed within 2x binary's {peak_bin:.3f} (cut cells without"
          " the sliver cutoff: 5.75)")
    # a slope tilted 30 degrees, its fall line 35 degrees off the x axis
    nx, ny, nz, steps = 40, 32, 36, 25
    th, az = np.radians(30.0), np.radians(35.0)
    nrm = np.array([np.sin(th) * np.cos(az), np.cos(th), np.sin(th) * np.sin(az)])
    down = np.array([np.cos(th) * np.cos(az), -np.sin(th), np.cos(th) * np.sin(az)])
    side = np.cross(nrm, down)
    p0 = np.array([6.0, 28.0, 6.0])
    x, y, z = node_grid3(nx, ny, nz)
    plane = f32((x - p0[0]) * nrm[0] + (y - p0[1]) * nrm[1] + (z - p0[2]) * nrm[2])
    expect = abs(G) * np.sin(th)

    def block(px, py, pz):
        r = np.stack([px, py, pz], -1) - p0
        h, s, c = r @ nrm, r @ down, r @ side
        return (h > 0.0) & (h < 5.0) & (s > 6.0) & (s < 18.0) & (np.abs(c) < 6.0)
    print(f"  slope off every axis, {nx}x{ny}x{nz}, {steps} steps: centre-of-mass"
          f" acceleration against g sin(theta) = {expect:.2f}")
    for name, mode in (("binary (staircase)", 1), ("cut cells", 2)):
        pos, typ = seed3d(nx, ny, nz, block)
        keep = phi_at(plane, pos) >= 0.0
        pos, typ = pos[keep], typ[keep]
        liq = pos[typ == 1]
        assert liq.min() > 1.0 and (liq.max(0) < np.array([nx, ny, nz]) - 1.0).all()
        d, frames = loop(f"slope3d{mode}", nx, ny, nz, pos.astype(np.float32), typ,
                         steps + 1, solid=(plane, mode))
        ss = [float((f[typ == 1] @ down).mean()) for f in frames]
        acc = 2.0 * np.polyfit(DT * np.arange(len(ss)), np.array(ss), 2)[0]
        pen = max(inside(plane, f) for f in frames)
        good = pen == 0 and (acc >= 0.95 * expect or mode == 1)
        ok &= good
        print(f"    {name:20s} {acc / expect:.2f} of g sin(theta); most inside {pen}"
              f"  {'OK' if good else 'FAIL'}")
    return ok


CHECKS = {"faces": check_faces, "step": check_step, "phase": check_phase,
          "loops flat": check_loops_flat, "loops 3d": check_loops_3d}


def main():
    os.makedirs(OUT, exist_ok=True)
    pick = sys.argv[1:] or list(CHECKS)
    results = {k: CHECKS[k]() for k in pick}
    print()
    for k, v in results.items():
        print(f"  {k:10s} {'PASS' if v else 'FAIL'}")
    ok = all(results.values())
    print("LEVEL-SET NODE:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
