"""Mesh colliders for the 3D loop: meshes in, solid_phi out, through the stock graph.

solid_phi_from_meshes() writes the meshes as .npy, runs collider_setup.json (the
shipped Simulation::Common::collider node, then our collider_field adapter: mesh ->
level set -> sampled at the grid nodes -> cells) and returns the field in the shape
run_bifrost3 wants for its solid argument. Static solids are computed once this way.
Also simple closed meshes to test with: a box, a UV sphere, and a bumpy "rock".
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def box_mesh(lo, hi):
    """A closed box, 8 points, 6 quads wound outward."""
    x0, y0, z0 = lo
    x1, y1, z1 = hi
    pts = np.array([[x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                    [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1]], np.float32)
    faces = [[0, 3, 2, 1], [4, 5, 6, 7], [0, 1, 5, 4], [2, 3, 7, 6], [1, 2, 6, 5], [0, 4, 7, 3]]
    fv = np.array([v for f in faces for v in f], np.int64)
    fo = np.arange(0, 4 * len(faces) + 1, 4, dtype=np.int64)
    return pts, fv, fo


def sphere_mesh(centre, radius, nu=96, nv=48, bump=None):
    """A closed UV sphere (quads, triangle fans at the poles). bump(theta, phi) may
    return a radius factor per vertex for a bumpy rock."""
    c = np.asarray(centre, np.float64)
    pts, fv, fo = [c + [0.0, radius, 0.0]], [], [0]
    for j in range(1, nv):
        th = np.pi * j / nv
        for i in range(nu):
            ph = 2.0 * np.pi * i / nu
            r = radius * (bump(th, ph) if bump else 1.0)
            pts.append(c + r * np.array([np.sin(th) * np.cos(ph), np.cos(th), np.sin(th) * np.sin(ph)]))
    pts.append(c + [0.0, -radius, 0.0])
    if bump:
        pts[0] = c + [0.0, radius * bump(0.0, 0.0), 0.0]
        pts[-1] = c + [0.0, -radius * bump(np.pi, 0.0), 0.0]
    top, bot = 0, len(pts) - 1

    def ring(j, i):
        return 1 + (j - 1) * nu + (i % nu)
    for i in range(nu):                       # top fan, wound outward
        fv += [top, ring(1, i + 1), ring(1, i)]
        fo.append(len(fv))
    for j in range(1, nv - 1):
        for i in range(nu):
            fv += [ring(j, i), ring(j, i + 1), ring(j + 1, i + 1), ring(j + 1, i)]
            fo.append(len(fv))
    for i in range(nu):
        fv += [bot, ring(nv - 1, i), ring(nv - 1, i + 1)]
        fo.append(len(fv))
    return (np.array(pts, np.float32), np.array(fv, np.int64), np.array(fo, np.int64))


def rock_mesh(centre, radius, seed=1, nu=96, nv=48, roughness=0.18):
    """A bumpy sphere: a few random low-order ripples on the radius."""
    rng = np.random.default_rng(seed)
    terms = [(rng.integers(1, 5), rng.integers(1, 5), rng.uniform(0, 2 * np.pi),
              rng.uniform(0, 2 * np.pi), rng.uniform(0.3, 1.0)) for _ in range(5)]
    norm = sum(a for *_, a in terms)

    def bump(th, ph):
        s = sum(a * np.cos(m * th + p1) * np.cos(n * ph + p2) for m, n, p1, p2, a in terms)
        return 1.0 + roughness * s / norm
    return sphere_mesh(centre, radius, nu, nv, bump)


def _save_mesh(d, suffix, mesh):
    pts, fv, fo = mesh
    paths = {}
    for name, arr in (("points", np.asarray(pts, np.float32)), ("face_vertex", np.asarray(fv, np.int64)),
                      ("face_offset", np.asarray(fo, np.int64))):
        p = os.path.join(d, name + suffix + ".npy")
        np.save(p, arr)
        paths["path_" + name + suffix] = p.replace("\\", "/")
    return paths


def solid_phi_from_meshes(tag, meshes, nx, ny, nz, cell_size=1.0, origin=(0.0, 0.0, 0.0),
                          detail_size=None, offset=0.0, thickening=0.0, min_hole_radius=0.0,
                          positions=False, out_dir=None):
    """The solid's signed distance at the (nx+1, ny+1, nz+1) grid nodes, in cells,
    negative inside, from one or two closed meshes given in world units: the grid node
    (i, j, k) sits at origin + (i, j, k) * cell_size. detail_size (world units) is the
    level set's voxel size, default one cell; offset grows the solid by that many world
    units; thickening by that many voxels (the shipped node's defaults are 0 and 1).
    Runs under out_dir, or out/collider_<tag>/. With positions, also returns the node positions
    the graph built."""
    assert 1 <= len(meshes) <= 2, "the setup graph takes one or two meshes"
    d = out_dir or os.path.join(HERE, "..", "out", "collider_" + tag)
    os.makedirs(d, exist_ok=True)
    ports = {"use_b": "true" if len(meshes) == 2 else "false",
             "detail_size": cell_size if detail_size is None else detail_size,
             "offset": offset, "thickening": thickening, "min_hole_radius": min_hole_radius,
             "nx": nx, "ny": ny, "nz": nz, "cell_size": cell_size,
             "path_solid_phi": os.path.join(d, "solid_phi.npy").replace("\\", "/")}
    if positions:
        ports["path_positions"] = os.path.join(d, "positions.npy").replace("\\", "/")
    ports.update(_save_mesh(d, "", meshes[0]))
    ports.update(_save_mesh(d, "_b", meshes[-1]))
    for f in ("solid_phi.npy", "positions.npy"):
        if os.path.exists(os.path.join(d, f)):
            os.remove(os.path.join(d, f))
    args = [os.path.join(HERE, "run_collider_setup.bat")]
    for k, v in ports.items():
        args += ["--set-port", k, str(v)]
    for k, v in zip(("origin_x", "origin_y", "origin_z"), origin):
        args += ["--set-port", k, str(v)]
    res = subprocess.run(args, capture_output=True, text=True)
    out = os.path.join(d, "solid_phi.npy")
    if res.returncode != 0 or not os.path.exists(out):
        sys.stderr.write(res.stdout[-3000:] + res.stderr[-2000:])
        raise SystemExit(f"bifcmd collider setup failed for {tag}")
    phi = np.load(out).reshape(nx + 1, ny + 1, nz + 1)
    if positions:
        return phi, np.load(os.path.join(d, "positions.npy")).reshape(nx + 1, ny + 1, nz + 1, 3)
    return phi
