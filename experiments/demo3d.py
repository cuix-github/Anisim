"""Milestone 6: the first 3D run, rendered to video.

A 3D two-phase dam break (liquid column under air in a closed box) runs entirely
inside one headless Bifrost graph (sim_3d.json: three staggered stock splat chains
feeding the multigrid step node inside an iterate). Per-substep particle dumps are
rendered as a shaded point cloud from a three-quarter view and assembled into an mp4
with ffmpeg.

Run from experiments:
  ..\\.venv\\Scripts\\python.exe demo3d.py smoke   # 16x16x8 sanity run
  ..\\.venv\\Scripts\\python.exe demo3d.py         # the real scene + video
"""
import os
import subprocess
import sys

import numpy as np
from PIL import Image, ImageDraw

HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "..", "reference2d", "results", "validation")
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
# The 3D loop runner lives beside its .bat launchers; re-exported for the experiments.
from loop3d import seed3d, lattice, run_bifrost3  # noqa: E402,F401


def render_frame(pos, typ, nx, ny, nz, scale=11):
    """Shaded liquid point cloud from a three-quarter view, box edges drawn."""
    yaw, pitch = np.radians(32), np.radians(20)
    cy_, sy_ = np.cos(yaw), np.sin(yaw)
    cp_, sp_ = np.cos(pitch), np.sin(pitch)

    def project(p):
        x = p[:, 0] - nx / 2
        y = p[:, 1] - ny / 2
        z = p[:, 2] - nz / 2
        xr = cy_ * x + sy_ * z
        zr = -sy_ * x + cy_ * z
        yr = cp_ * y - sp_ * zr
        depth = cp_ * zr + sp_ * y
        return xr, yr, depth

    W = int((nx + nz) * 0.82 * scale // 2 * 2)
    H = int((ny + nz) * 0.95 * scale // 2 * 2)
    im = Image.new("RGB", (W, H), (252, 252, 254))
    d = ImageDraw.Draw(im)
    corners = np.array([[i * nx, j * ny, k * nz] for i in (0, 1) for j in (0, 1) for k in (0, 1)],
                       float)
    cx, cy2, _ = project(corners)
    px = (cx * scale + W / 2).astype(int)
    py = (H / 2 - cy2 * scale).astype(int)
    edges = [(0, 1), (0, 2), (0, 4), (1, 3), (1, 5), (2, 3), (2, 6), (3, 7),
             (4, 5), (4, 6), (5, 7), (6, 7)]
    for a, b in edges:
        d.line([(px[a], py[a]), (px[b], py[b])], fill=(205, 208, 214), width=2)

    liq = pos[typ == 1]
    if len(liq) > 220000:
        liq = liq[:: len(liq) // 220000 + 1]
    x, y, depth = project(liq)
    order = np.argsort(depth)   # far first, near painted last
    x, y, depth = x[order], y[order], depth[order]
    t = (depth - depth.min()) / max(depth.max() - depth.min(), 1e-9)
    img = np.array(im)
    xi = np.clip((x * scale + W / 2).astype(int), 0, W - 2)
    yi = np.clip((H / 2 - y * scale).astype(int), 0, H - 2)
    col = np.stack([15 + 60 * t, 70 + 90 * t, 170 + 70 * t], 1).astype(np.uint8)
    for dx in (0, 1):
        for dy in (0, 1):
            img[yi + dy, xi + dx] = col
    return Image.fromarray(img)


def main():
    smoke = len(sys.argv) > 1 and sys.argv[1] == "smoke"
    if smoke:
        d, typ, diags = run_bifrost3("smoke", 16, 16, 8,
                                     lambda x, y, z: y < 8, 4, 0.02)
        print("smoke ok:", diags[:, 0].tolist(), "max div", diags[:, 2].max())
        return

    nx, ny, nz, dt, sub = 64, 48, 32, 0.015, 480
    d, typ, diags = run_bifrost3("dam", nx, ny, nz,
                                 lambda x, y, z: (x < 20) & (y < 36), sub, dt)
    print(f"run done: CG iters mean {diags[:, 0].mean():.1f} max {diags[:, 0].max():.0f}; "
          f"worst divergence {diags[:, 2].max():.2e}; max CFL {diags[:, 3].max() * dt:.2f}")
    frames_dir = os.path.join(d, "png")
    os.makedirs(frames_dir, exist_ok=True)
    out_idx = 0
    for k in range(0, sub, 3):
        pos = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
        render_frame(pos, typ, nx, ny, nz).save(
            os.path.join(frames_dir, f"fr{out_idx:04d}.png"))
        out_idx += 1
    video = os.path.join(RESULTS, "bifrost_dam3d.mp4")
    subprocess.run(["ffmpeg", "-y", "-framerate", "24",
                    "-i", os.path.join(frames_dir, "fr%04d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", video],
                   check=True, capture_output=True)
    print("wrote", os.path.normpath(video))


if __name__ == "__main__":
    main()
