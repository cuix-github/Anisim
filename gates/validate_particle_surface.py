"""The water surface and mist (maya/particle_surface.json, particle_mist.json): the
mesh the renderer is given must sit where the particles are, and the mist must be
the particles the mesh leaves out. Light: seconds, plus the real frame when asked.

  small   a synthetic cache: a block of liquid at 8 particles a cell (16 x 6 x 8 cells),
          one lone droplet and a vertical sheet one particle thick. The surface closes
          over the block at its fill level (the top within half a cell of the
          particles' top) and reaches half a cell past its sides (the spheres'
          radius, point_size 0.5); the droplet, one voxel across, is not in the
          mesh; offset moves the surface out by its value. The mist pass takes the
          droplet, at least half the sheet, and from the block only particles along
          its convex edges (the smoothing rounds those off), under a tenth of it.
  cache   a real play cache (default out/collider_showcase/boulders512_play24): its
          last frame through the compound, the mesh's bounds against the cache's
          grid, and the wall per frame (bifcmd startup differenced out).

Run from the repo root:  .venv\\Scripts\\python.exe gates\\validate_particle_surface.py [small] [cache [<dir>]]
"""
import json
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "maya"))
from play_cache import write_frame  # noqa: E402

OUT = os.path.join(HERE, "..", "out", "surface_gate")
BAT = os.path.join(HERE, "..", "maya", "run_particle_surface_test.bat")
MIST_BAT = os.path.join(HERE, "..", "maya", "run_particle_mist_test.bat")
BLOCK = (16, 6, 8)                 # cells of liquid
DROPLET = np.array([20.0, 10.0, 4.0], np.float32)
SHEET_X = 24.0                     # a vertical sheet one particle thick, 4 x 8 cells


def surface(cache, frame, tag, **ports):
    """The compound at one frame of a cache; returns the mesh's vertices and the wall."""
    c = os.path.abspath(cache).replace("\\", "/")
    o = os.path.abspath(os.path.join(OUT, tag)).replace("\\", "/")
    os.makedirs(o, exist_ok=True)
    args = [BAT, "--set-port", "liquid_pattern", c + "/liquid.####",
            "--set-port", "frame", str(frame), "--set-port", "path_positions", o + "/pos.npy"]
    for k, v in ports.items():
        args += ["--set-port", k, str(v)]
    t0 = time.time()
    res = subprocess.run(args, capture_output=True, text=True)
    wall = time.time() - t0
    if res.returncode != 0 or "'ok_positions': (bool) true" not in res.stdout:
        sys.stderr.write(res.stdout[-2000:] + res.stderr[-2000:])
        raise SystemExit(f"bifcmd particle_surface failed for {tag}")
    return np.load(os.path.join(o, "pos.npy")), wall


def mist_pass(cache, frame, tag, **ports):
    """The mist compound at one frame: the mist positions, every particle's depth, wall."""
    c = os.path.abspath(cache).replace("\\", "/")
    o = os.path.abspath(os.path.join(OUT, tag)).replace("\\", "/")
    os.makedirs(o, exist_ok=True)
    args = [MIST_BAT, "--set-port", "liquid_pattern", c + "/liquid.####",
            "--set-port", "frame", str(frame), "--set-port", "path_positions", o + "/mist.npy",
            "--set-port", "path_depth", o + "/depth.npy"]
    for k, v in ports.items():
        args += ["--set-port", k, str(v)]
    t0 = time.time()
    res = subprocess.run(args, capture_output=True, text=True)
    wall = time.time() - t0
    if res.returncode != 0 or "'ok_depth': (bool) true" not in res.stdout:
        sys.stderr.write(res.stdout[-2000:] + res.stderr[-2000:])
        raise SystemExit(f"bifcmd particle_mist failed for {tag}")
    return np.load(os.path.join(o, "mist.npy")), np.load(os.path.join(o, "depth.npy")), wall


def small_cache():
    """Frame 2: the block at 0.5-cell spacing (centres from 0.25), the droplet and the
    sheet; returns the cache and how many particles each part has."""
    cache = os.path.join(OUT, "small_cache")
    axes = [np.arange(0.25, n, 0.5, dtype=np.float32) for n in BLOCK]
    pos = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    ys, zs = np.meshgrid(axes[1][:8], axes[2], indexing="ij")
    sheet = np.stack([np.full(ys.size, SHEET_X, np.float32), ys.ravel(), zs.ravel()], -1)
    pos = np.vstack([pos, DROPLET[None], sheet])
    vel = np.zeros_like(pos)
    write_frame(cache, 2, pos, vel, np.ones(len(pos), np.float32))
    return cache, len(pos)


def parts(pos):
    """Which part of the synthetic cache each particle (or vertex) belongs to."""
    droplet = np.abs(pos - DROPLET).max(1) < 2.0
    sheet = (np.abs(pos[:, 0] - SHEET_X) < 2.0) & ~droplet
    return ~droplet & ~sheet, droplet, sheet


def check_small():
    print("=== small: a block of liquid and a lone droplet through the compound")
    cache, n = small_cache()
    verts, wall = surface(cache, 2, "small")
    block = verts[parts(verts)[0]]
    lo, hi = block.min(0), block.max(0)
    top = float(np.median(block[block[:, 1] > BLOCK[1] - 1.0, 1]))
    near = int((np.linalg.norm(verts - DROPLET, axis=1) < 1.0).sum())
    sides = all(abs(lo[i] + 0.25) <= 0.5 and abs(hi[i] - (BLOCK[i] + 0.25)) <= 0.5
                for i in (0, 2))
    fill = abs(top - BLOCK[1]) <= 0.5
    good = len(verts) > 0 and sides and fill and near == 0
    print(f"    {n:,} particles -> {len(verts):,} vertices in {wall:.1f} s; block bounds"
          f" {np.round(lo, 2)} .. {np.round(hi, 2)} (sides within half a cell {sides});"
          f" top at {top:.2f} against the fill level {BLOCK[1]} ({fill}); droplet"
          f" vertices {near} (none: the mist's)  {'OK' if good else 'FAIL'}")
    out, _ = surface(cache, 2, "small_offset", offset=0.5)
    out = out[parts(out)[0]]
    grown = float(out.max(0)[0] - hi[0])
    good_off = 0.25 <= grown <= 0.75
    print(f"    offset 0.5: the far side moves out by {grown:.2f} cells  "
          f"{'OK' if good_off else 'FAIL'}")
    pos = np.load(os.path.join(cache, "liquid.0002.npy"))
    mist, depth, mwall = mist_pass(cache, 2, "small_mist")
    in_block, in_droplet, in_sheet = parts(pos)
    m_block, m_droplet, m_sheet = parts(mist)
    faces = np.stack([pos[:, i] < 1.25 for i in range(3)]
                     + [pos[:, i] > BLOCK[i] - 1.25 for i in range(3)], 1).sum(1)
    on_edge = faces >= 2                         # within two particles of two faces
    block_mist = (depth > 0.0) & in_block
    edges_only = bool((on_edge | ~block_mist).all())
    share = block_mist.sum() / in_block.sum()
    good_mist = (m_droplet.sum() == 1 and m_sheet.sum() >= in_sheet.sum() // 2
                 and edges_only and share < 0.1 and len(mist) == (depth > 0.0).sum())
    print(f"    mist in {mwall:.1f} s: {len(mist):,} particles; the droplet {m_droplet.sum()};"
          f" the sheet {m_sheet.sum()} of {in_sheet.sum()}; the block {m_block.sum()}"
          f" ({100 * share:.1f} percent, all on its edges {edges_only})  "
          f"{'OK' if good_mist else 'FAIL'}")
    return good and good_off and good_mist


def check_cache(cache=None):
    cache = cache or os.path.join(HERE, "..", "out", "collider_showcase", "boulders512_play24")
    print(f"=== cache: the last frame of {os.path.normpath(cache)}")
    if not os.path.exists(os.path.join(cache, "cache.json")):
        print("    no cache there  FAIL")
        return False
    with open(os.path.join(cache, "cache.json")) as f:
        info = json.load(f)
    frame, grid = info["frames"], info.get("grid")
    verts, wall = surface(cache, frame, "cache")
    _, startup = surface(small_cache()[0], 2, "startup")
    lo, hi = verts.min(0), verts.max(0)
    inside = grid is None or (np.all(lo >= -1.0) and np.all(hi <= np.array(grid) + 1.0))
    good = len(verts) > 0 and inside
    print(f"    frame {frame}: {info['liquid_particles'][-1]:,} particles -> {len(verts):,}"
          f" vertices, bounds {np.round(lo, 1)} .. {np.round(hi, 1)} within the grid {grid}"
          f" ({inside})  {'OK' if good else 'FAIL'}")
    print(f"    bifcmd wall {wall:.1f} s against {startup:.1f} s for the small cache (startup):"
          f" about {wall - startup:.1f} s to read, voxelize, smooth and contour the frame")
    return good


def main():
    args = sys.argv[1:] or ["small"]
    ok = True
    i = 0
    while i < len(args):
        if args[i] == "small":
            ok &= check_small()
        elif args[i] == "cache":
            d = args[i + 1] if i + 1 < len(args) and args[i + 1] not in ("small", "cache") else None
            if d:
                i += 1
            ok &= check_cache(d)
        else:
            raise SystemExit("unknown check " + args[i])
        i += 1
    print("PARTICLE SURFACE:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
