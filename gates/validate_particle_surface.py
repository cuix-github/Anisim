"""The water surface and mist (maya/particle_surface.json, particle_mist.json): the
mesh the renderer is given must sit where the particles are, and the mist must be
the particles the mesh leaves out. Light: seconds, plus the real frame when asked.

  small   a synthetic cache: a block of liquid at 8 particles a cell (16 x 6 x 8 cells),
          one lone droplet, a vertical sheet one particle thick and a lone ball two
          cells across. The surface closes
          over the block at its fill level (the top within half a cell of the
          particles' top) and reaches half a cell past its sides (the spheres'
          radius, point_size 0.5); the droplet, one voxel across, is not in the
          mesh; offset moves the surface out by its value. The mist pass takes the
          droplet (a droplet radius: spray, out of the surface and always mist), at
          least half the sheet, and from the block only particles along its convex
          edges (the smoothing rounds those off), under a tenth of it.
          Then the feature-preserving flow (ST-FLIP Appendix B, mcf_iterations 30 after
          no Gaussian pass): the mask itself, sampled at probes, is under 0.1 inside
          the ball and the sheet and over 0.95 in the block and far outside; under it
          the ball keeps at least seven tenths of its extent where the unmasked flow
          takes more, the block's top edge, bulk to the mask, is rounded down by at
          least a tenth of a cell (the flow acts where g is 1), the blur read at the
          ball's centre is that of a Gaussian within 15 percent of mask_sigma, and
          with the flow on the sheet (kept by the mesh) leaves the mist.
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
BALL = np.array([12.0, 14.0, 4.0], np.float32)   # a lone ball: 4 x 4 x 4 particles, 2 cells
MCF_ITERS = 30                     # ST-FLIP's k_psi
PROBE_NAMES = (("ball centre", BALL), ("ball, half a cell in", BALL + [0.5, 0, 0]),
               ("ball, inside its surface", BALL + [1.0, 0, 0]),
               ("ball, outside", BALL + [1.3, 0, 0]), ("sheet centre", (SHEET_X, 2.0, 4.0)),
               ("block centre", (8.0, 3.0, 4.0)), ("block top, inside", (8.0, 5.9, 4.0)),
               ("block top, outside", (8.0, 6.3, 4.0)), ("far outside", (8.0, 12.0, 4.0)))
PROTECTED = ("ball, half a cell in", "ball, inside its surface", "sheet centre")
MASK_SIGMA, MASK_DETAIL = 2.0, 0.5   # the mask's blur in cells and the surface voxel


def ball_blur(sigma, n=40):
    """The gate's ball (64 spheres of radius 1 voxel on a 4^3 lattice at 1-voxel
    spacing, the half-voxel ramp) blurred by a Gaussian of `sigma` voxels, read at its
    centre: the value the mask's blur should give there."""
    ax = np.arange(n) - (n - 1) / 2.0
    X, Y, Z = np.meshgrid(ax, ax, ax, indexing="ij")
    sdf = np.full(X.shape, 1e9)
    for cx in (-1.5, -0.5, 0.5, 1.5):
        for cy in (-1.5, -0.5, 0.5, 1.5):
            for cz in (-1.5, -0.5, 0.5, 1.5):
                sdf = np.minimum(sdf, np.sqrt((X - cx) ** 2 + (Y - cy) ** 2 + (Z - cz) ** 2) - 1.0)
    b = np.clip(0.5 - sdf, 0.0, 1.0)
    k = np.exp(-ax ** 2 / (2 * sigma ** 2))
    k /= k.sum()
    for axis in range(3):
        b = np.apply_along_axis(lambda v: np.convolve(v, k, mode="same"), axis, b)
    return 0.5 * float(b[n // 2, n // 2, n // 2] + b[(n - 1) // 2, (n - 1) // 2, (n - 1) // 2])


def sigma_of(blur):
    """The sigma (voxels) whose ball_blur is `blur`, by bisection."""
    lo, hi = 0.5, 12.0
    for _ in range(22):
        mid = 0.5 * (lo + hi)
        if ball_blur(mid) > blur:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
SMOOTHED = ("block centre", "block top, inside", "block top, outside", "far outside")
PROBES = np.array([p for _, p in PROBE_NAMES], np.float32)


def radii_of(cache, frame):
    """The cache's radius.#### when it has one; otherwise zeros for `frame` written
    under the gate's out (the play cache's radii are the escaped droplets')."""
    c = os.path.abspath(cache).replace("\\", "/")
    if os.path.exists(os.path.join(cache, f"radius.{frame:04d}.npy")):
        return c + "/radius.####"
    d = os.path.abspath(os.path.join(OUT, "zero_radius")).replace("\\", "/")
    os.makedirs(d, exist_ok=True)
    n = np.load(os.path.join(cache, f"liquid.{frame:04d}.npy"), mmap_mode="r").shape[0]
    np.save(os.path.join(d, f"radius.{frame:04d}.npy"), np.zeros(n, np.float32))
    return d + "/radius.####"


def surface(cache, frame, tag, probes=None, **ports):
    """The compound at one frame of a cache; returns the mesh's vertices, the wall, and
    the mask's density, its blur and the mask sampled at `probes` ((N, 3) positions in
    cells; None samples the origin alone)."""
    c = os.path.abspath(cache).replace("\\", "/")
    o = os.path.abspath(os.path.join(OUT, tag)).replace("\\", "/")
    os.makedirs(o, exist_ok=True)
    np.save(os.path.join(o, "probes.npy"),
            np.zeros((1, 3), np.float32) if probes is None else np.asarray(probes, np.float32))
    args = [BAT, "--set-port", "liquid_pattern", c + "/liquid.####",
            "--set-port", "radius_pattern", radii_of(cache, frame),
            "--set-port", "frame", str(frame), "--set-port", "path_positions", o + "/pos.npy",
            "--set-port", "path_probes", o + "/probes.npy", "--set-port", "path_density", o + "/density.npy",
            "--set-port", "path_blur", o + "/blur.npy", "--set-port", "path_mask", o + "/mask.npy"]
    for k, v in ports.items():
        args += ["--set-port", k, str(v)]
    t0 = time.time()
    res = subprocess.run(args, capture_output=True, text=True)
    wall = time.time() - t0
    if res.returncode != 0 or "'ok_mask': (bool) true" not in res.stdout:
        sys.stderr.write(res.stdout[-2000:] + res.stderr[-2000:])
        raise SystemExit(f"bifcmd particle_surface failed for {tag}")
    return (np.load(os.path.join(o, "pos.npy")), wall, np.load(os.path.join(o, "density.npy")),
            np.load(os.path.join(o, "blur.npy")), np.load(os.path.join(o, "mask.npy")))


def mist_pass(cache, frame, tag, probes=None, **ports):
    """The mist compound at one frame: the mist positions, every particle's depth, the
    wall, the mist particles' weights, and the spray density sampled at `probes`
    (an (N, 3) array of positions in cells; None samples the origin alone)."""
    c = os.path.abspath(cache).replace("\\", "/")
    o = os.path.abspath(os.path.join(OUT, tag)).replace("\\", "/")
    os.makedirs(o, exist_ok=True)
    np.save(os.path.join(o, "probes.npy"),
            np.zeros((1, 3), np.float32) if probes is None else np.asarray(probes, np.float32))
    args = [MIST_BAT, "--set-port", "liquid_pattern", c + "/liquid.####",
            "--set-port", "radius_pattern", radii_of(cache, frame),
            "--set-port", "frame", str(frame), "--set-port", "path_positions", o + "/mist.npy",
            "--set-port", "path_depth", o + "/depth.npy", "--set-port", "path_weights", o + "/weights.npy",
            "--set-port", "path_probes", o + "/probes.npy", "--set-port", "path_samples", o + "/samples.npy"]
    for k, v in ports.items():
        args += ["--set-port", k, str(v)]
    t0 = time.time()
    res = subprocess.run(args, capture_output=True, text=True)
    wall = time.time() - t0
    if res.returncode != 0 or "'ok_samples': (bool) true" not in res.stdout:
        sys.stderr.write(res.stdout[-2000:] + res.stderr[-2000:])
        raise SystemExit(f"bifcmd particle_mist failed for {tag}")
    return (np.load(os.path.join(o, "mist.npy")), np.load(os.path.join(o, "depth.npy")), wall,
            np.load(os.path.join(o, "weights.npy")), np.load(os.path.join(o, "samples.npy")))


DROPLET_R, DX, CLOUD_R, SPRAY_EPS = 1e-4, 0.25, 1.0, 1e4   # the spray checks' constants


def small_cache():
    """Frame 2: the block at 0.5-cell spacing (centres from 0.25), the droplet and the
    sheet; returns the cache and how many particles each part has."""
    cache = os.path.join(OUT, "small_cache")
    axes = [np.arange(0.25, n, 0.5, dtype=np.float32) for n in BLOCK]
    pos = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
    ys, zs = np.meshgrid(axes[1][:8], axes[2], indexing="ij")
    sheet = np.stack([np.full(ys.size, SHEET_X, np.float32), ys.ravel(), zs.ravel()], -1)
    off = np.arange(-0.75, 0.76, 0.5, dtype=np.float32)
    ball = BALL + np.stack(np.meshgrid(off, off, off, indexing="ij"), -1).reshape(-1, 3)
    n_block = len(pos)
    pos = np.vstack([pos, DROPLET[None], sheet, ball])
    vel = np.zeros_like(pos)
    # the droplet alone has a droplet radius (DROPLET_R); the rest are water
    radius = np.zeros(len(pos), np.float32)
    radius[n_block] = DROPLET_R
    write_frame(cache, 2, pos, vel, np.ones(len(pos), np.float32), radius=radius)
    return cache, len(pos)


def parts(pos):
    """Which part of the synthetic cache each particle (or vertex) belongs to."""
    droplet = np.abs(pos - DROPLET).max(1) < 2.0
    ball = np.abs(pos - BALL).max(1) < 2.5
    sheet = (np.abs(pos[:, 0] - SHEET_X) < 2.0) & ~droplet
    return ~droplet & ~sheet & ~ball, droplet, sheet, ball


def check_small():
    print("=== small: a block of liquid and a lone droplet through the compound")
    cache, n = small_cache()
    verts, wall = surface(cache, 2, "small")[:2]
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
    out, _ = surface(cache, 2, "small_offset", offset=0.5)[:2]
    out = out[parts(out)[0]]
    grown = float(out.max(0)[0] - hi[0])
    good_off = 0.25 <= grown <= 0.75
    print(f"    offset 0.5: the far side moves out by {grown:.2f} cells  "
          f"{'OK' if good_off else 'FAIL'}")
    pos = np.load(os.path.join(cache, "liquid.0002.npy"))
    # the spray density sampled on a lattice of quarter cells around the droplet, for
    # the integral below
    ax = np.arange(-CLOUD_R - 1.0, CLOUD_R + 1.0 + 1e-6, 0.25, dtype=np.float32)
    probes = DROPLET + np.stack(np.meshgrid(ax, ax, ax, indexing="ij"), -1).reshape(-1, 3)
    mist, depth, mwall, weights, samples = mist_pass(
        cache, 2, "small_mist", probes=probes, metres_per_cell=DX, cloud_radius=CLOUD_R,
        spray_detail=0.5, mist_density=1.0, default_radius=2 * DROPLET_R, spray_epsilon=SPRAY_EPS)
    in_block, in_droplet, in_sheet, _ = parts(pos)
    m_block, m_droplet, m_sheet, _ = parts(mist)
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
    # the spray drawing: a parcel of dx^3 / 8 of water as drops of radius r has a
    # cross-section of 1.5 V_p / r, in cells 1.5 dx / (8 r); the droplet's cloud must
    # integrate to it (the stored channel is the sum over spray_epsilon), and its
    # weight must be twice a radius-less mist particle's (default_radius 2 r)
    w_drop = float(weights[m_droplet][0]) if m_droplet.any() else 0.0
    w_rest = float(np.median(weights[~m_droplet])) if (~m_droplet).any() else 0.0
    integral = float(samples.sum()) * 0.25 ** 3 * SPRAY_EPS
    expect = 1.5 * DX / (8.0 * DROPLET_R)
    ratio = integral / expect
    good_spray = abs(ratio - 1.0) < 0.1 and abs(w_drop / max(w_rest, 1e-30) - 2.0) < 1e-3
    print(f"    spray: the droplet's cloud integrates to {integral:.1f} cells^2 against its"
          f" cross-section {expect:.1f} ({ratio:.3f}); its weight {w_drop:.1f}, a radius-less"
          f" mist particle's {w_rest:.1f}  {'OK' if good_spray else 'FAIL'}")
    return good and good_off and good_mist and good_spray and check_smoothing(cache, verts)


def check_smoothing(cache, control):
    """The feature-preserving flow against the raw surface and the unmasked flow."""
    print("=== smoothing: mean-curvature flow under the self-quotient mask (ST-FLIP App. B)")
    arms = {"raw": dict(smooth_iterations=0),
            "masked": dict(smooth_iterations=0, mcf_iterations=MCF_ITERS),
            "unmasked": dict(smooth_iterations=0, mcf_iterations=MCF_ITERS, mask_floor=1.0)}
    out = {"control": (control, None)}
    g = {}
    for tag, ports in arms.items():
        verts, wall, dens, blur, mask = surface(cache, 2, "smooth_" + tag, probes=PROBES, **ports)
        out[tag] = (verts, wall)
        if tag == "masked":
            print("    the mask at the probes (density psi0, its blur, g):")
            for (name, _), d, b, m in zip(PROBE_NAMES, dens, blur, mask):
                print(f"      {name:26s} {d:7.3f} {b:7.3f} {m:7.3f}")
                g[name] = float(m)
            sig = sigma_of(float(blur[0]))
    mask_ok = all(g[n] < 0.1 for n in PROTECTED) and all(g[n] > 0.95 for n in SMOOTHED)
    nominal = MASK_SIGMA / MASK_DETAIL
    sig_ok = abs(sig - nominal) <= 0.15 * nominal
    print(f"    the blur at the ball's centre is a Gaussian's of sigma {sig:.2f} voxels against"
          f" the mask's {nominal:.1f}  {'OK' if sig_ok else 'FAIL'}")
    rows = {}
    for tag, (verts, wall) in out.items():
        block, _, sheet, ball = parts(verts)
        b = verts[block]
        # the block's top: its edge (the flow rounds it) and its middle (flat already)
        edge = b[(b[:, 1] > 5.0) & (b[:, 0] > 1.4) & (b[:, 0] < 1.8) & (b[:, 2] > 2.0) & (b[:, 2] < 6.0), 1]
        mid = b[(b[:, 1] > 5.0) & (b[:, 0] > 5.0) & (b[:, 0] < 11.0) & (np.abs(b[:, 2] - 4.0) < 0.6), 1]
        edge = float(np.median(edge)) if len(edge) else float("nan")
        rough = float(mid.std()) if len(mid) else float("nan")
        ext = float((verts[ball].max(0) - verts[ball].min(0)).max()) if ball.any() else 0.0
        rows[tag] = (edge, rough, ext, int(sheet.sum()))
        print(f"    {tag:9s} top edge at {edge:.3f}, middle roughness {rough:.4f}; ball extent"
              f" {ext:.2f} cells; sheet vertices {sheet.sum():,}" + (f"; {wall:.1f} s" if wall else ""))
    kept = rows["masked"][2] >= 0.7 * rows["raw"][2]
    protects = rows["unmasked"][2] < rows["masked"][2]
    acts = rows["raw"][0] - rows["masked"][0] >= 0.1
    # the mist under the flow: the sheet, kept by the mesh now, leaves the mist; the
    # droplet stays; the ball's particles are inside
    pos = np.load(os.path.join(cache, "liquid.0002.npy"))
    mist = mist_pass(cache, 2, "smooth_mist", smooth_iterations=0, mcf_iterations=MCF_ITERS)[0]
    _, in_droplet, in_sheet, in_ball = parts(pos)
    _, m_droplet, m_sheet, m_ball = parts(mist)
    mist_ok = m_droplet.sum() == 1 and m_sheet.sum() < in_sheet.sum() // 2 and m_ball.sum() == 0
    print(f"    the mist under the flow: the droplet {m_droplet.sum()}, the sheet {m_sheet.sum()} of"
          f" {in_sheet.sum()} (the control's {in_sheet.sum() // 2} or more), the ball"
          f" {m_ball.sum()} of {in_ball.sum()}  {'OK' if mist_ok else 'FAIL'}")
    good = mask_ok and sig_ok and kept and protects and acts and mist_ok
    print(f"    the mask under 0.1 at the ball and the sheet, over 0.95 in the block and outside"
          f" ({mask_ok}); the ball keeps {rows['masked'][2] / rows['raw'][2]:.2f} of its extent under"
          f" the mask ({kept}), {rows['unmasked'][2] / rows['raw'][2]:.2f} without it (less: {protects});"
          f" the block's top edge rounded down by {rows['raw'][0] - rows['masked'][0]:.2f} cells"
          f" ({acts})  {'OK' if good else 'FAIL'}")
    return good


def check_cache(cache=None):
    cache = cache or os.path.join(HERE, "..", "out", "collider_showcase", "boulders512_play24")
    print(f"=== cache: the last frame of {os.path.normpath(cache)}")
    if not os.path.exists(os.path.join(cache, "cache.json")):
        print("    no cache there  FAIL")
        return False
    with open(os.path.join(cache, "cache.json")) as f:
        info = json.load(f)
    frame, grid = info["frames"], info.get("grid")
    verts, wall = surface(cache, frame, "cache")[:2]
    _, startup = surface(small_cache()[0], 2, "startup")[:2]
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
