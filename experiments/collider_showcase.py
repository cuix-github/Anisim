"""Mesh colliders, stage 3's showcase: the finale-class bore (256x64x128, production
config: ST at 8x dt, adaptive air, escaped droplets) breaking over three boulders that
exist only as meshes. The boulders go through the stock graph (the shipped collider
node, collider_field) once, as solid_phi; the solver sees cut cells.

  ..\\.venv\\Scripts\\python.exe collider_showcase.py smoke        # tiny plumbing check, ~2 min
  ..\\.venv\\Scripts\\python.exe collider_showcase.py cost         # 8 substeps with and without the rocks
  ..\\.venv\\Scripts\\python.exe collider_showcase.py stage 0      # ... 1, 2, 3: ~65 min each
  ..\\.venv\\Scripts\\python.exe collider_showcase.py finalize     # verdict, video
  ..\\.venv\\Scripts\\python.exe collider_showcase.py night        # all of the above, in order

The background runner caps a task at two hours, so the run goes in chunks that resume
from full-state dumps (the bore suite's scheme), one launch each. Each chunk renders
its frames to PNG and then deletes its particle dumps (about 160 MB a frame), keeping
the resume state; the video and the numbers land in out/collider_showcase/. Checks:
particles inside a boulder in any rendered frame, all and liquid, with the position
and depth of each one found (the field itself is the judge), the run's peak speed,
particle counts, honest clocks per chunk (startup included, labeled).

Ran 2026-10-14 (the README's stage-3 record carries the numbers): 20.16 s simulated in
4.13 h, 85 to 93 s a substep, no particle inside at any chunk's end, 4 of 84 frames
with one particle inside, the one with its depth recorded 0.00005 cells in; the
boulders cost nothing measurable (107.2 s against 107.18 s a substep without them).
"""
import json
import os
import shutil
import subprocess
import sys
import time

import numpy as np
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
sys.path.insert(0, os.path.join(HERE, "..", "gates"))
sys.path.insert(0, os.path.join(HERE, "..", "reference2d"))
from pfflip2d import Sim  # noqa: E402
from loop3d import run_bifrost3, seed3d, load_final_state  # noqa: E402
from colliders import rock_mesh, solid_phi_from_meshes  # noqa: E402
from run_bore_suite import scene  # noqa: E402

OUT = os.path.join(HERE, "..", "out", "collider_showcase")
PORTS = {"preconditioner": 1, "st": 1, "escape": 1, "adapt": 1}
DT = 0.12
# the finale grid, 42 substeps a chunk (5 s simulated), four chunks = 20 s
FULL = ((256, 64, 128), 42, 4, 2)      # grid, substeps per chunk, chunks, dump_every
SMOKE = ((64, 16, 32), 6, 2, 1)


def rocks(nx, ny, nz):
    """Three boulders on the bed, part buried, downstream of the reservoir."""
    h0 = 0.125 * ny
    spec = [((0.50 * nx, 0.55 * h0, 0.34 * nz), 0.165 * ny, 3),
            ((0.63 * nx, 0.45 * h0, 0.62 * nz), 0.13 * ny, 5),
            ((0.74 * nx, 0.60 * h0, 0.42 * nz), 0.19 * ny, 8)]
    return [rock_mesh(c, r, seed=s) for c, r, s in spec]


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


def field(tag, nx, ny, nz):
    """The boulders' solid_phi through the graph, cached beside the run."""
    path = os.path.join(OUT, f"solid_phi_{tag}.npy")
    if os.path.exists(path):
        return np.load(path)
    meshes = rocks(nx, ny, nz)
    t0 = time.time()
    # the setup graph takes two meshes a run: the third goes through as a second run
    # and the fields meet in a minimum, exactly what the adapter does for N colliders
    phi = np.minimum(solid_phi_from_meshes(f"{tag}_ab", meshes[:2], nx, ny, nz),
                     solid_phi_from_meshes(f"{tag}_c", meshes[2:], nx, ny, nz))
    wall = time.time() - t0
    os.makedirs(OUT, exist_ok=True)
    np.save(path, phi)
    with open(os.path.join(OUT, f"field_{tag}.json"), "w") as f:
        json.dump({"nodes": int(phi.size), "wall_s": round(wall, 1),
                   "inside_nodes": int((phi < 0).sum()),
                   "meshes": [{"points": int(len(m[0])), "faces": int(len(m[2]) - 1)}
                              for m in meshes]}, f, indent=1)
    return phi


# ------------------------------------------------------------------ rendering


class View:
    """Three-quarter view; particles paint far-to-near into a z-buffer, the boulders
    are painted as shaded quads with their own depth, and the nearer one wins."""

    def __init__(self, nx, ny, nz, meshes, scale=None):
        self.nx, self.ny, self.nz = nx, ny, nz
        self.scale = scale or max(2, int(1000 / (nx + nz)))
        yaw, pitch = np.radians(30), np.radians(18)
        self.cy, self.sy = np.cos(yaw), np.sin(yaw)
        self.cp, self.sp = np.cos(pitch), np.sin(pitch)
        self.W = int((nx + nz) * 0.80 * self.scale // 2 * 2)
        self.H = int((ny + nz) * 0.85 * self.scale // 2 * 2)
        self.rock_img, self.rock_depth = self.paint_rocks(meshes)

    def project(self, p):
        x = p[..., 0] - self.nx / 2
        y = p[..., 1] - self.ny / 2
        z = p[..., 2] - self.nz / 2
        xr = self.cy * x + self.sy * z
        zr = -self.sy * x + self.cy * z
        yr = self.cp * y - self.sp * zr
        return (xr * self.scale + self.W / 2, self.H / 2 - yr * self.scale,
                self.cp * zr + self.sp * y)       # larger depth = nearer

    def paint_rocks(self, meshes):
        img = Image.new("RGB", (self.W, self.H), (0, 0, 0))
        dep = Image.new("F", (self.W, self.H), -1e30)
        di, dd = ImageDraw.Draw(img), ImageDraw.Draw(dep)
        light = np.array([0.4, 0.8, -0.45])
        light /= np.linalg.norm(light)
        for pts, fv, fo in meshes:
            px, py, pd = self.project(pts.astype(np.float64))
            faces = [fv[fo[q]:fo[q + 1]] for q in range(len(fo) - 1)]
            depth = np.array([pd[f].mean() for f in faces])
            for q in np.argsort(depth):               # painter's order, near last
                f = faces[q]
                a, b, c = pts[f[0]], pts[f[1]], pts[f[2]]
                n = np.cross(b - a, c - a)
                nl = np.linalg.norm(n)
                if nl < 1e-12:
                    continue
                lam = max(0.0, float(n @ light) / nl)
                shade = int(70 + 120 * lam)
                col = (shade, int(shade * 0.93), int(shade * 0.82))
                poly = list(zip(px[f], py[f]))
                di.polygon(poly, fill=col)
                dd.polygon(poly, fill=float(depth[q]))
        return np.array(img), np.array(dep)

    def frame(self, pos, liquid):
        img = np.full((self.H, self.W, 3), (252, 252, 254), np.uint8)
        zbuf = np.full((self.H, self.W), -1e30)
        liq = pos[liquid]
        if len(liq) > 400000:
            liq = liq[:: len(liq) // 400000 + 1]
        x, y, depth = self.project(liq.astype(np.float64))
        order = np.argsort(depth)
        x, y, depth = x[order], y[order], depth[order]
        t = (depth - depth.min()) / max(depth.max() - depth.min(), 1e-9)
        xi = np.clip(x.astype(int), 0, self.W - 2)
        yi = np.clip(y.astype(int), 0, self.H - 2)
        col = np.stack([15 + 60 * t, 70 + 90 * t, 170 + 70 * t], 1).astype(np.uint8)
        for dx in (0, 1):
            for dy in (0, 1):
                img[yi + dy, xi + dx] = col
                zbuf[yi + dy, xi + dx] = depth
        rock = self.rock_depth > zbuf
        img[rock] = self.rock_img[rock]
        return Image.fromarray(img)


# ------------------------------------------------------------------ the run


def calibration():
    cal = Sim(64, 64, adapt=True, st=True)
    cal.seed(lambda x, y: np.ones_like(x, dtype=bool))
    cal.calibrate()
    return dict(PORTS, rho0_face=cal.rho0_face, coarse_gain=cal.coarse_gain, ws_epsilon=16)


def run_chunk(tag, grid, substeps, chunk, every, with_rocks=True, render=True):
    nx, ny, nz = grid
    mask, _, _, h1, h0, rx = scene(nx, ny, nz)
    phi = field(tag, nx, ny, nz) if with_rocks else None
    ports = calibration()

    def seed(*_):
        pos, typ = seed3d(nx, ny, nz, mask)
        if phi is not None:
            keep = phi_at(phi, pos) >= 0.0
            pos, typ = pos[keep], typ[keep]
        return pos, typ
    init = load_final_state(os.path.join(HERE, "..", "out", f"loop3_{tag}_c{chunk - 1}")) \
        if chunk else None
    t0 = time.time()
    d, typ, diag = run_bifrost3(f"{tag}_c{chunk}", nx, ny, nz, mask, substeps, DT,
                                extra_ports=ports, seed_fn=seed, init_arrays=init,
                                final_state=True, step_offset=chunk * substeps,
                                dump_every=every, solid=(phi, 2) if phi is not None else None)
    wall = time.time() - t0
    # audit and render every dumped frame, then drop the dumps (the state stays)
    frames_dir = os.path.join(OUT, f"{tag}_png")
    os.makedirs(frames_dir, exist_ok=True)
    view = View(nx, ny, nz, rocks(nx, ny, nz)) if (render and with_rocks) else None
    counts, pens, pens_liquid, peaks, offenders = [], [], [], [], []
    for k in range(0, substeps, every):
        p = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
        ph = np.load(os.path.join(d, f"ph.{k:04d}.npy"))
        counts.append(int(len(p)))
        # inside a boulder: every particle, and the liquid ones (the rendered ones;
        # air particles move far faster, past the push-out's 4.5 cells a step)
        dist = phi_at(phi, p.astype(np.float64)) if phi is not None else np.ones(len(p))
        ins = dist < 0.0
        pens.append(int(ins.sum()))
        pens_liquid.append(int((ins & (ph >= 0.5)).sum()))
        for i in np.flatnonzero(ins)[:8]:   # the evidence: where and how deep
            offenders.append({"substep": chunk * substeps + k, "position": p[i].round(3).tolist(),
                              "depth_cells": round(float(-dist[i]), 4),
                              "liquid": bool(ph[i] >= 0.5)})
        if view is not None:
            idx = (chunk * substeps + k) // every
            view.frame(p, ph >= 0.5).save(os.path.join(frames_dir, f"fr{idx:04d}.png"))
        for name in ("frame", "esc", "scl", "ph"):
            f = os.path.join(d, f"{name}.{k:04d}.npy")
            if os.path.exists(f):
                os.remove(f)
    peaks = [float(diag[k][-1]) for k in range(substeps)]
    iters = [float(diag[k][0]) for k in range(substeps)]
    entry = {"chunk": chunk, "substeps": substeps, "wall_s": round(wall, 1),
             "wall_per_substep_s": round(wall / substeps, 2),
             "sim_time_s": [round(chunk * substeps * DT, 2), round((chunk + 1) * substeps * DT, 2)],
             "particles": counts, "inside_a_boulder": pens, "inside_a_boulder_liquid": pens_liquid,
             "inside_particles": offenders, "peak_speed": max(peaks),
             "iterations_mean": round(float(np.mean(iters)), 1), "with_rocks": with_rocks}
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"{tag}_c{chunk}.json"), "w") as f:
        json.dump(entry, f, indent=1)
    print(json.dumps(entry), flush=True)
    return entry


def finalize(tag, n_chunks):
    chunks = []
    for c in range(n_chunks):
        with open(os.path.join(OUT, f"{tag}_c{c}.json")) as f:
            chunks.append(json.load(f))
    frames_dir = os.path.join(OUT, f"{tag}_png")
    video = os.path.join(OUT, f"{tag}.mp4")
    subprocess.run(["ffmpeg", "-y", "-framerate", "12", "-i", os.path.join(frames_dir, "fr%04d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "21", video],
                   check=True, capture_output=True)
    summary = {"tag": tag, "grid": chunks and FULL[0], "dt": DT, "config": PORTS,
               "simulated_s": round(sum(c["substeps"] for c in chunks) * DT, 2),
               "wall_s_total": round(sum(c["wall_s"] for c in chunks), 1),
               "wall_per_substep_s": [c["wall_per_substep_s"] for c in chunks],
               "inside_a_boulder_max": max(max(c["inside_a_boulder"]) for c in chunks),
               "inside_a_boulder_frames": sum(sum(1 for n in c["inside_a_boulder"] if n)
                                              for c in chunks),
               # the liquid split exists only for chunks audited after it was added
               "inside_a_boulder_liquid_max": max(max(c.get("inside_a_boulder_liquid", [0]))
                                                  for c in chunks),
               "peak_speed_cells_per_s": max(c["peak_speed"] for c in chunks),
               "particles_start_end": [chunks[0]["particles"][0], chunks[-1]["particles"][-1]],
               "iterations_mean": [c["iterations_mean"] for c in chunks],
               "video": os.path.basename(video), "frames": len(os.listdir(frames_dir))}
    with open(os.path.join(OUT, f"{tag}_summary.json"), "w") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1), flush=True)
    print("FINALIZED", tag, flush=True)


def cost(grid, substeps=8):
    """The honest clock: the same bore with and without the boulders, one chunk
    each, no rendering, startup included and the same in both."""
    res = {}
    for with_rocks in (False, True):
        tag = "cost_rocks" if with_rocks else "cost_plain"
        for f in os.listdir(OUT) if os.path.exists(OUT) else []:
            if f.startswith(tag):
                os.remove(os.path.join(OUT, f))
        e = run_chunk(tag, grid, substeps, 0, 1, with_rocks=with_rocks, render=False)
        res[tag] = e["wall_per_substep_s"]
        shutil.rmtree(os.path.join(HERE, "..", "out", f"loop3_{tag}_c0"), ignore_errors=True)
    res["ratio"] = round(res["cost_rocks"] / res["cost_plain"], 3)
    with open(os.path.join(OUT, "cost.json"), "w") as f:
        json.dump(dict(res, grid=grid, substeps=substeps), f, indent=1)
    print(json.dumps(res), flush=True)


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    os.makedirs(OUT, exist_ok=True)
    if cmd == "smoke":
        grid, sub, n, every = SMOKE
        for c in range(n):
            run_chunk("smoke", grid, sub, c, every)
        finalize("smoke", n)
    elif cmd == "cost":
        cost(FULL[0])
    elif cmd == "stage":
        grid, sub, n, every = FULL
        run_chunk("boulders", grid, sub, int(sys.argv[2]), every)
    elif cmd == "finalize":
        finalize("boulders", FULL[2])
    elif cmd == "night":
        grid, sub, n, every = FULL
        for c in range(n):
            run_chunk("boulders", grid, sub, c, every)
        finalize("boulders", n)
        cost(grid)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
