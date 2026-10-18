"""The wave basin: the breaker's flume made five times wider, its piston replaced by
a wave zone at the west wall (step_3d's zone_* ports and `waves` array) that writes
two oblique wave trains of different periods into 12 m of water; they cross, run onto
the breaker's 1:10 beach and break along it. The scene, in metres: the basin
172 x 24 x 80 (along, up, across), a closed box; the zone the first 24 m, the trains'
heights ramped over 16 s; train A 3.2 m at 8 s heading 15 degrees toward the far wall
(z = 80), train B 2.4 m at 6.5 s heading 20 degrees toward the near wall (z = 0); the
beach's toe at 36 m, the shoreline at 156 m. Both trains sit in the plunging range of
the surf-similarity number on this beach (0.56 and 0.52). Production configuration as
the breaker's plain run (ST, escaped particles, adaptive air, multigrid, the fused
transfer); no droplet model and no spray emission, so the liquid particle set stays
one set from dump to dump and the video-rate blend holds, the escaped particles drawn
as droplets at render time (the paper's rule) as the breaker's were.

  ..\\.venv\\Scripts\\python.exe wave_basin.py chunk C --plan "S1:E1,S2:E2,..." [--cpm M]
      chunk C of the run: the plan lists each chunk's seconds and dump interval in
      seconds; chunk C resumes from chunk C - 1's state (the zone's clock carries on).
  ..\\.venv\\Scripts\\python.exe wave_basin.py measure [--cpm M]
      every dump so far: the highest crest beyond the zone (x, z, height), the spray
      count, a top-view height map per dump under renders/audit, summary.json.
  ..\\.venv\\Scripts\\python.exe wave_basin.py cache T0 T1 [--fps F] [--cpm M] [--radius R]
                                            [--crop SEAWARD SHOREWARD ACROSS]
      the play cache for the dumps between T0 and T1 seconds: at --fps F (24) frames
      blended between the dumps (maya/play_cache.py's Hermite blend), at --fps 0 one
      frame per dump; radius.####.npy marks the escaped particles as droplets of
      --radius R metres (0.0001); camera.json the tracking camera, one eye and target
      per frame in cells (build_surface_scene.py reads it through SURF_CAMERA_PATH),
      following the highest crest, smoothed; the beach as an OBJ beside the inputs.
      --crop (metres) keeps only the liquid within that box of the followed crest, the
      cache then named crop: the render meshes the water around the camera, not the
      whole basin. --kind K names the cache instead (stills, for a few frames apart).
      --spray-above H keeps a droplet radius only for escaped particles more than H cells
      above the bulk water (1 for maya/presets/paper_clean.ps1; 0, the default, keeps them
      all, the paper look).
  ..\\.venv\\Scripts\\python.exe wave_basin.py section KIND [--z Z] [--back B] [--up U]
                                            [--look-up L] [--look-in I] [--cpm M]
      a section of cache KIND: the liquid beyond z = Z metres (55), the cut face showing
      the wave in profile, a side camera following the crest B metres in front of the
      cut (55) and U above still water (0.75), looking at L above still water (0.75), I
      beyond the cut (0); the breaker's framing is --back 13 --up 12 --look-up 4.25
      --look-in 3.75. --spray-above H keeps a droplet radius only for escaped particles
      more than H cells above the bulk water (calm water clean; 0 keeps them all), --out K
      names the cache (section). The beach over the section's width.

Everything lands under out/runs/<tag>/ (graphs/runs.py): tag wavebasin<M>.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
sys.path.insert(0, os.path.join(HERE, "..", "reference2d"))
sys.path.insert(0, os.path.join(HERE, "..", "maya"))
from solitary_beach import calibration, write_beach_obj  # noqa: E402
from solitary_wave import beach_bed  # noqa: E402
from loop3d import DUMPS, run_bifrost3, seed3d, load_final_state  # noqa: E402
from play_cache import hermite, write_index, write_liquid  # noqa: E402
from runs import RunPaths, prepare  # noqa: E402
from PIL import Image  # noqa: E402

# the scene in metres
H, TOE, SLOPE = 12.0, 36.0, 0.1
LX, LY, LZ = 172.0, 24.0, 80.0
SHORE = TOE + H / SLOPE
ZONE, RAMP = 24.0, 8.0
# the sea (2026-10-18): two lead waves heading for the beach, nine companions around
# them and a hundred small waves in all directions, each a row of height (m), period (s),
# heading (degrees from +x toward +z) and timing (rad); the second lead's timing chosen
# so the two line up at ALIGN (x m, z m, t s), the breaking zone mid-run. A fixed seed
# replays the same sea. (The first basin, 2026-10-17, ran two trains: 3.2 m at 8 s,
# 15 degrees, and 2.4 m at 6.5 s, -20 degrees, with a 16 s ramp.)
SEA = dict(seed=7,
           leads=((3.0, 9.5, 12.0, 0.0), (1.5, 11.0, -5.0)),
           align=(100.0, 40.0, 28.0),
           companions=(9, 1.2, (8.0, 11.0), 12.0, 15.0),      # count, typical height, periods, heading, spread
           small=(100, 0.8, (2.0, 5.0), 12.0, 60.0))


def _k(T, h=H):
    om = 2 * np.pi / T
    k = om * om / 9.8
    for _ in range(50):
        th, ch = np.tanh(k * h), np.cosh(k * h)
        k -= (9.8 * k * th - om * om) / (9.8 * (th + k * h / (ch * ch)))
    return float(k)


def _scale_to(heights, typical):
    """Heights scaled so the waves they add up to are about `typical` tall (four times
    the spread of their sum)."""
    h = np.asarray(heights, float)
    return h * (typical / 4) / np.sqrt(np.sum((h / 2) ** 2 / 2))


def sea_rows(sea=SEA):
    """The sea's rows, in metres, seconds, degrees and radians: the leads, then the
    companions (random heights, periods across their range, headings within their spread),
    then the small waves (frequencies evenly spread and jittered, heights falling with
    frequency as an ocean's do, times a random factor; headings fanned, most near the
    main direction)."""
    rng = np.random.default_rng(sea["seed"])
    (H0, T0, d0, p0), (H1, T1, d1) = sea["leads"]
    xa, za, ta = sea["align"]

    def phase(T, d, ph):
        r = np.radians(d)
        return _k(T) * (xa * np.cos(r) + za * np.sin(r)) - 2 * np.pi / T * ta + ph
    p1 = (phase(T0, d0, p0) - phase(T1, d1, 0.0)) % (2 * np.pi)
    rows = [(H0, T0, d0, p0), (H1, T1, d1, float(p1))]
    n, typ, (lo, hi), mean, spread = sea["companions"]
    T = rng.uniform(lo, hi, n)
    d = mean + rng.uniform(-spread, spread, n)
    Hs = _scale_to(rng.rayleigh(1.0, n), typ)
    rows += [(float(a), float(b), float(c), float(rng.uniform(0, 2 * np.pi))) for a, b, c in zip(Hs, T, d)]
    n, typ, (lo, hi), mean, spread = sea["small"]
    f = np.linspace(1 / hi, 1 / lo, n) + rng.uniform(-0.4, 0.4, n) * (1 / lo - 1 / hi) / n
    Hs = _scale_to(f ** -2.5 * rng.rayleigh(1.0, n), typ)
    u = rng.uniform(-1, 1, n)
    d = np.clip(mean + spread * np.sign(u) * (1 - np.sqrt(1 - np.abs(u))), -60, 60)
    rows += [(float(a), float(1 / b), float(c), float(rng.uniform(0, 2 * np.pi))) for a, b, c in zip(Hs, f, d)]
    return rows
# the tracking camera's offset from the crest it follows, metres: seaward, up, toward the near wall
CAMERA_OFFSET = (-24.0, 16.0, -16.0)
# an escaped particle counts as spray in the cache only this many cells above the bulk
# surface: the escape flag also marks particles of the jittered free surface itself,
# which read as specks when every one is drawn as a droplet cloud (0 = every one)
SPRAY_ABOVE = 0.0
# crests within this of a side wall are the wall's reflection doubling the train, never followed
WALL_MARGIN = 10.0


class Run:
    def __init__(self, cpm=4.0):
        self.cpm = cpm
        self.grid = (int(LX * cpm), int(LY * cpm), int(LZ * cpm))
        self.dt = 0.1 / cpm
        self.tag = f"wavesea{cpm:g}"


def scene(run):
    m = run.cpm
    return {"m": m, "h": H * m, "toe": TOE * m, "slope": SLOPE, "shore": SHORE * m, "g": 9.8 * m,
            "zone": ZONE * m, "ramp": RAMP}


def ports(run, sc):
    return dict(calibration(), gravity=-sc["g"], zone_length=sc["zone"], zone_depth=sc["h"], zone_ramp=sc["ramp"])


def waves(run):
    """The node's `waves`: four floats a wave, the height in cells."""
    return np.asarray([(Hw * run.cpm, T, th, ph) for Hw, T, th, ph in sea_rows()], np.float32).ravel()


def beach_phi(run, sc):
    nx, ny, nz = run.grid
    i, j = np.meshgrid(np.arange(nx + 1, dtype=float), np.arange(ny + 1, dtype=float), indexing="ij")
    phi2 = (j - beach_bed(i, sc["toe"], sc["slope"])) * np.cos(np.arctan(sc["slope"]))
    return np.repeat(phi2[:, :, None], nz + 1, axis=2).astype(np.float32)


def seed(run, sc):
    """Water to h over the whole basin, air above, nothing inside the sand."""
    nx, ny, nz = run.grid
    pos, typ = seed3d(nx, ny, nz, lambda x, y, z: y < sc["h"])
    pos = pos.astype(np.float32)
    keep = pos[:, 1] > beach_bed(pos[:, 0], sc["toe"], sc["slope"])
    pos, typ = pos[keep], typ[keep]
    n = len(pos)
    zeros = np.zeros(n, np.float32)
    return {"positions": pos, "velocities": np.zeros((n, 3), np.float32),
            "particle_phase": typ.astype(np.int32), "mass": np.where(typ == 1, 1000.0, 1.0).astype(np.float32),
            "momx": zeros, "momy": zeros, "momz": zeros, "tau": zeros, "wt": np.ones(n, np.float32),
            "wtph": typ.astype(np.float32), "scale": np.ones(n, np.float32), "phase_st": typ.astype(np.float32),
            "radius": zeros, "spray": zeros}


# ------------------------------------------------------------------ measuring


def surface_map(pos, nx, nz):
    """The highest liquid particle of each (x, z) column; -inf where none."""
    ci = np.clip(pos[:, 0].astype(np.int64), 0, nx - 1)
    ck = np.clip(pos[:, 2].astype(np.int64), 0, nz - 1)
    top = np.full((nx, nz), -np.inf)
    np.maximum.at(top, (ci, ck), pos[:, 1])
    return top


def blur(a, r):
    """A box blur of radius r over both axes, -inf treated as missing."""
    wet = np.isfinite(a)
    v = np.where(wet, a, 0.0)
    k = np.ones(2 * r + 1)
    num, den = v.copy(), wet.astype(float)
    for ax in (0, 1):
        num = np.apply_along_axis(lambda s: np.convolve(s, k, "same"), ax, num)
        den = np.apply_along_axis(lambda s: np.convolve(s, k, "same"), ax, den)
    return np.where(den > 0, num / np.maximum(den, 1e-9), -np.inf)


def crests(top, sc, run):
    """The crest candidates of a surface map: local maxima of the blurred map beyond
    the zone and short of the swash, above still water, (x, z, height) in cells."""
    m = run.cpm
    sm = blur(top, int(m))
    lo, hi = int(sc["zone"] + 10 * m), int(sc["shore"] - 4 * m)
    # the side walls double an oblique train against its own reflection: a crest
    # within WALL_MARGIN metres of either is the wall's, not one to follow
    zlo, zhi = int(WALL_MARGIN * m), int(top.shape[1] - WALL_MARGIN * m)
    reg = sm[lo:hi, zlo:zhi]
    eta = reg - sc["h"]
    out = []
    if not np.isfinite(reg).any():
        return out
    best = float(np.nanmax(np.where(np.isfinite(eta), eta, -np.inf)))
    if best <= 0.0:
        return out
    i = np.arange(1, reg.shape[0] - 1)[:, None]
    k = np.arange(1, reg.shape[1] - 1)[None, :]
    c = reg[1:-1, 1:-1]
    peak = (c >= reg[:-2, 1:-1]) & (c >= reg[2:, 1:-1]) & (c >= reg[1:-1, :-2]) & (c >= reg[1:-1, 2:]) & \
        (eta[1:-1, 1:-1] >= 0.7 * best) & np.isfinite(c)
    for ii, kk in zip(*np.nonzero(peak)):
        out.append((float(lo + ii + 1 + 0.5), float(zlo + kk + 1 + 0.5), float(c[ii, kk])))
    return out


def spray_only(p, e, nx, nz, above=None):
    """The escaped particles more than `above` cells (SPRAY_ABOVE) above the bulk
    surface (the surface of the particles not flagged, blurred a cell)."""
    above = SPRAY_ABOVE if above is None else above
    top = blur(surface_map(p[~e], nx, nz), 1)
    ci = np.clip(p[:, 0].astype(np.int64), 0, nx - 1)
    ck = np.clip(p[:, 2].astype(np.int64), 0, nz - 1)
    return e & (p[:, 1] > top[ci, ck] + above)


def follow(cands_per_frame, start=None):
    """One crest per frame from the candidates: the highest at first, then the one
    nearest the previous pick (the highest where nothing is near); frames with no
    candidate repeat the last."""
    track, prev = [], start
    for cands in cands_per_frame:
        if not cands:
            track.append(prev)
            continue
        if prev is None:
            pick = max(cands, key=lambda c: c[2])
        else:
            near = [c for c in cands if (c[0] - prev[0]) ** 2 + (c[1] - prev[1]) ** 2 < 40.0 ** 2]
            pick = min(near, key=lambda c: (c[0] - prev[0]) ** 2 + (c[1] - prev[1]) ** 2) if near \
                else max(cands, key=lambda c: c[2])
        track.append(pick)
        prev = pick
    first = next((t for t in track if t is not None), None)
    return [t if t is not None else first for t in track]


def smooth_track(track, sigma):
    """A Gaussian over the frames of each coordinate (reflected at the ends)."""
    a = np.asarray(track, float)
    if sigma <= 0 or len(a) < 3:
        return a
    r = int(3 * sigma)
    k = np.exp(-0.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    k /= k.sum()
    out = np.empty_like(a)
    for c in range(a.shape[1]):
        out[:, c] = np.convolve(np.pad(a[:, c], (r, r), mode="reflect"), k, "valid")
    return out


def camera_path(track, run):
    """Eye and target per frame, in cells: the crest, and the camera CAMERA_OFFSET
    metres seaward, up and toward the near wall of it."""
    m = run.cpm
    off = np.array(CAMERA_OFFSET) * m
    frames = []
    for x, z, y in track:
        eye = (x + off[0], y + off[1], z + off[2])
        frames.append([float(v) for v in eye] + [float(x), float(y), float(z)])
    return frames


def dumps_of(run, paths):
    """Every dump of the run so far: (time, directory, index), from the chunk receipts."""
    out = []
    for c in paths.chunks():
        if not os.path.exists(paths.receipt(f"c{c}.json")):
            continue
        r = json.load(open(paths.receipt(f"c{c}.json")))
        d = paths.dumps(c)
        out += [((r["offset"] + k + 1) * run.dt, d, k) for k in range(0, r["substeps"], r["every"])]
    return out


def load_dump(d, k, with_vel=False):
    pos = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
    ph = np.load(os.path.join(d, f"ph.{k:04d}.npy"))
    e = os.path.join(d, f"esc.{k:04d}.npy")
    esc = np.load(e) if os.path.exists(e) else np.zeros(len(pos), np.float32)
    vel = np.load(os.path.join(d, f"vel.{k:04d}.npy")) if with_vel else None
    return pos, ph, esc, vel


def height_image(top, sc, run, path):
    """A top view of the surface height, the still level grey, crests light, troughs dark."""
    eta = np.where(np.isfinite(top), top - sc["h"], 0.0) / run.cpm
    img = np.clip(128 + 40 * eta, 0, 255).astype(np.uint8).T[::-1]      # z up the image, near wall at the bottom
    Image.fromarray(img).resize((img.shape[1] * 2, img.shape[0] * 2)).save(path)


def report(run, sc, paths, entries):
    nx, ny, nz = run.grid
    m = sc["m"]
    frames = paths.renders("audit")
    os.makedirs(frames, exist_ok=True)
    log, cands = [], []
    for n, (t, d, k) in enumerate(entries):
        pos, ph, esc, _ = load_dump(d, k)
        liquid = ph >= 0.5
        top = surface_map(pos[liquid], nx, nz)
        cs = crests(top, sc, run)
        cands.append(cs)
        best = max(cs, key=lambda c: c[2]) if cs else None
        rec = {"t": round(t, 3), "liquid": int(liquid.sum()), "spray_particles": int((liquid & (esc > 0.5)).sum()),
               "particles": int(len(pos))}
        if best is not None:
            rec.update({"crest_x_m": round(best[0] / m, 2), "crest_z_m": round(best[1] / m, 2),
                        "crest_height_m": round((best[2] - sc["h"]) / m, 2),
                        "depth_under_crest_m": round((sc["h"] - float(beach_bed(best[0], sc["toe"], sc["slope"]))) / m, 2),
                        "crests": len(cs)})
        log.append(rec)
        height_image(top, sc, run, os.path.join(frames, f"top{n:04d}.png"))
    track = follow(cands)
    summary = {"tag": run.tag, "grid": list(run.grid), "cells_per_metre": run.cpm, "dt": run.dt,
               "scene_m": {"h": H, "domain": [LX, LY, LZ], "zone": ZONE, "ramp_s": RAMP, "sea": SEA,
                           "toe": TOE, "slope": SLOPE, "shore": SHORE},
               "followed_crest": [None if tr is None else [round(tr[0] / m, 2), round(tr[1] / m, 2), round((tr[2] - sc["h"]) / m, 2)]
                                  for tr in track],
               "log": log}
    with open(prepare(paths.receipt("summary.json")), "w") as f:
        json.dump(summary, f, indent=1)
    for e in log[::4]:
        print("  ", json.dumps(e))
    return summary


# ------------------------------------------------------------------ the cache


def write_cache(run, sc, paths, entries, fps, radius, crop=None, kind=None, spray_above=None):
    """Frames at fps between the first and last entry (one per dump at fps 0), the
    liquid blended between the dumps around each frame, the escaped particles of the
    nearer dump marked with `radius`; then the tracking camera from the frames. With
    crop = (seaward, shoreward, across) in metres only the liquid within that box of
    the followed crest is written (the cache `crop`), the camera found first from the
    whole surface: the render then meshes the water around the camera, not the basin."""
    nx, ny, nz = run.grid
    cache = paths.cache(kind or ("crop" if crop else "full"))
    os.makedirs(cache, exist_ok=True)
    times_s, counts, cands = [], [], []
    if fps <= 0:
        plan = [(t, i, i, 0.0) for i, (t, _, _) in enumerate(entries)]
    else:
        t0, t1 = entries[0][0], entries[-1][0]
        n = int(np.floor((t1 - t0) * fps + 1e-6)) + 1
        ts = [t0 + f / fps for f in range(n)]
        times = [e[0] for e in entries]
        plan = []
        for t in ts:
            i = min(int(np.searchsorted(times, t, side="right")) - 1, len(entries) - 2)
            plan.append((t, i, i + 1, (t - times[i]) / (times[i + 1] - times[i])))
    loaded = {}

    def get(i):
        if i not in loaded:
            loaded.clear() if len(loaded) > 2 else None
            t, d, k = entries[i]
            pos, ph, esc, vel = load_dump(d, k, with_vel=True)
            liq = ph >= 0.5
            loaded[i] = (pos[liq].astype(np.float32), vel[liq].astype(np.float32), esc[liq] > 0.5)
        return loaded[i]

    def frame_of(ia, ib, u):
        pa, va, ea = get(ia)
        if ia == ib or u <= 0.0:
            return pa, va, ea
        pb, vb, eb = get(ib)
        if len(pa) != len(pb):
            raise SystemExit(f"cache: liquid counts differ between dumps {ia} and {ib} ({len(pa)} vs {len(pb)})")
        p, v = hermite(pa, va, pb, vb, entries[ib][0] - entries[ia][0], min(max(u, 0.0), 1.0))
        return p, v, (ea if u < 0.5 else eb)

    sigma = max(1.0, 0.5 * fps) if fps > 0 else 1.0
    sm = None
    if crop:
        # the first pass finds the crest to follow from the whole surface
        for f, (t, ia, ib, u) in enumerate(plan, start=1):
            p, _, _ = frame_of(ia, ib, u)
            cands.append(crests(surface_map(p, nx, nz), sc, run))
        sm = smooth_track(follow(cands), sigma)
        cands = []
    for f, (t, ia, ib, u) in enumerate(plan, start=1):
        p, v, e = frame_of(ia, ib, u)
        above = SPRAY_ABOVE if spray_above is None else spray_above
        if above > 0 and e.any():
            e = spray_only(p, e, nx, nz, above)
        if sm is None:
            cands.append(crests(surface_map(p, nx, nz), sc, run))
        else:
            cx, cz = sm[f - 1][0], sm[f - 1][1]
            m = run.cpm
            keep = ((p[:, 0] > cx - crop[0] * m) & (p[:, 0] < cx + crop[1] * m)
                    & (p[:, 2] > cz - crop[2] * m) & (p[:, 2] < cz + crop[2] * m))
            p, v, e = p[keep], v[keep], e[keep]
        r = np.where(e, np.float32(radius), np.float32(0.0))
        counts.append(write_liquid(cache, f, p, np.linalg.norm(v, axis=1), r))
        times_s.append(round(t, 4))
        if f % 24 == 0:
            print(f"    frame {f}/{len(plan)} at {t:.3f} s, {counts[-1]:,} liquid", flush=True)
    if sm is None:
        sm = smooth_track(follow(cands), sigma)
    path = camera_path(sm, run)
    with open(os.path.join(cache, "camera.json"), "w") as fh:
        json.dump({"frames": path, "crest_cells": [list(map(float, tr)) for tr in sm], "offset_m": CAMERA_OFFSET,
                   "crop_m": crop}, fh)
    write_index(cache, counts, grid=list(run.grid), times_s=times_s, source=run.tag, fps=fps,
                escaped_radius_m=radius, spray_above_cells=SPRAY_ABOVE, camera="camera.json", crop_m=crop)
    return cache, counts, times_s


def write_section(run, sc, paths, kind, zc_m, back_m=55.0, up_m=0.75, look_up_m=0.75, look_in_m=0.0,
                  spray_above=0.0, out_kind="section"):
    """A section of a cache's frames: the liquid beyond z = zc_m metres only, so the cut
    face shows the wave in profile, like a flume's glass wall; camera.json a side camera
    per frame at the followed crest's x, back_m metres in front of the cut and up_m above
    still water, looking at look_up_m above still water, look_in_m beyond the cut (the
    breaker's crop camera: 13, 12, 4.25, 3.75); the beach over the section's width.
    spray_above > 0 keeps a droplet radius only for the escaped particles more than that
    many cells above the slab's bulk water (the artist's knob for spray on calm water;
    0 keeps every escaped particle a droplet). Writes the cache `out_kind`."""
    nx, ny, nz = run.grid
    m = run.cpm
    src = paths.cache(kind)
    idx = json.load(open(os.path.join(src, "cache.json")))
    crest = json.load(open(os.path.join(src, "camera.json")))["crest_cells"]
    dst = paths.cache(out_kind)
    os.makedirs(dst, exist_ok=True)
    zc = zc_m * m
    counts, frames = [], []
    for f in range(1, idx["frames"] + 1):
        p = np.load(os.path.join(src, f"liquid.{f:04d}.npy"))
        s = np.load(os.path.join(src, f"speed.{f:04d}.npy"))
        rp = os.path.join(src, f"radius.{f:04d}.npy")
        r = np.load(rp) if os.path.exists(rp) else np.zeros(len(p), np.float32)
        keep = p[:, 2] >= zc
        p, s, r = p[keep], s[keep], r[keep]
        if spray_above > 0 and (r > 0).any():
            e = spray_only(p, r > 0, nx, nz, spray_above)
            r = np.where(e, r, np.float32(0.0))
        counts.append(write_liquid(dst, f, p, s, r))
        cx = float(crest[f - 1][0])
        frames.append([cx, float(sc["h"] + up_m * m), float(zc - back_m * m),
                       cx, float(sc["h"] + look_up_m * m), float(zc + look_in_m * m)])
    with open(os.path.join(dst, "camera.json"), "w") as fh:
        json.dump({"frames": frames, "section_z_m": zc_m, "back_m": back_m, "up_m": up_m,
                   "look_up_m": look_up_m, "look_in_m": look_in_m}, fh)
    toe, x1 = sc["toe"], float(nx)
    y1 = float((nx - toe) * sc["slope"])
    v = [(toe, 0, zc), (x1, 0, zc), (x1, y1, zc), (toe, 0, nz), (x1, 0, nz), (x1, y1, nz)]
    with open(os.path.join(dst, "beach.obj"), "w") as fh:
        fh.write("o beach\n")
        fh.writelines(f"v {a:.4f} {b:.4f} {c:.4f}\n" for a, b, c in v)
        fh.writelines("f " + " ".join(map(str, q)) + "\n" for q in ((1, 3, 2), (4, 5, 6), (1, 2, 5, 4), (2, 3, 6, 5), (3, 1, 4, 6)))
    write_index(dst, counts, grid=list(run.grid), times_s=idx.get("times_s"), source=run.tag, section_z_m=zc_m,
                from_cache=kind, camera="camera.json")
    return dst, counts


# ------------------------------------------------------------------ commands


def run_chunk(run, chunk, seconds, every_s, offset, init):
    sc = scene(run)
    paths = RunPaths(run.tag)
    nx, ny, nz = run.grid
    substeps = int(round(seconds / run.dt))
    every = max(1, int(round(every_s / run.dt)))
    print(f"{run.tag} chunk {chunk}: grid {run.grid}, dt {run.dt}, {substeps} substeps from substep {offset}, "
          f"a dump every {every}; {len(init['positions']):,} particles", flush=True)
    t0 = time.time()
    d, _, diag = run_bifrost3(f"{run.tag}_c{chunk}", nx, ny, nz, None, substeps, run.dt,
                              extra_ports=ports(run, sc), init_arrays=init, final_state=True,
                              step_offset=offset, dump_every=every, solid=(beach_phi(run, sc), 2),
                              dumps=DUMPS + ("vel",), fused=True, out_dir=paths.dumps(chunk), waves=waves(run))
    wall = time.time() - t0
    entry = {"chunk": chunk, "substeps": substeps, "every": every, "dt": run.dt, "offset": offset,
             "grid": list(run.grid), "wall_s": round(wall, 1), "wall_per_substep_s": round(wall / substeps, 2),
             "sim_time_s": [round(offset * run.dt, 3), round((offset + substeps) * run.dt, 3)],
             "iterations_mean": round(float(np.mean(diag[:, 0])), 1),
             "peak_speed_cells_per_s": round(float(diag[:, 3].max()), 2)}
    with open(prepare(paths.receipt(f"c{chunk}.json")), "w") as f:
        json.dump(entry, f, indent=1)
    print(json.dumps(entry), flush=True)
    return entry


def opt(argv, name, default, cast=float):
    return cast(argv[argv.index(name) + 1]) if name in argv else default


def main(argv):
    cmd = argv[0] if argv else ""
    run = Run(opt(argv, "--cpm", 4.0))
    sc = scene(run)
    paths = RunPaths(run.tag)
    if cmd == "chunk":
        c = int(argv[1])
        plan = [tuple(float(v) for v in p.split(":")) for p in opt(argv, "--plan", "", str).split(",")]
        offset = sum(int(round(s / run.dt)) for s, _ in plan[:c])
        init = load_final_state(paths.dumps(c - 1)) if c else seed(run, sc)
        run_chunk(run, c, plan[c][0], plan[c][1], offset, init)
    elif cmd == "measure":
        report(run, sc, paths, dumps_of(run, paths))
    elif cmd == "cache":
        t0, t1 = float(argv[1]), float(argv[2])
        entries = [e for e in dumps_of(run, paths) if t0 - 1e-9 <= e[0] <= t1 + 1e-9]
        if len(entries) < 1:
            raise SystemExit("no dumps in that window")
        crop = None
        if "--crop" in argv:
            i = argv.index("--crop")
            crop = tuple(float(v) for v in argv[i + 1:i + 4])
        cache, counts, times = write_cache(run, sc, paths, entries, opt(argv, "--fps", 24.0), opt(argv, "--radius", 1e-4),
                                           crop=crop, kind=opt(argv, "--kind", None, str),
                                           spray_above=opt(argv, "--spray-above", None))
        write_beach_obj(prepare(paths.input("beach.obj")), run, sc)
        print(f"wrote {len(counts)} frames to {cache}: t {times[0]} .. {times[-1]}; camera.json beside them")
    elif cmd == "section":
        dst, counts = write_section(run, sc, paths, argv[1], opt(argv, "--z", 55.0), opt(argv, "--back", 55.0),
                                    opt(argv, "--up", 0.75), opt(argv, "--look-up", 0.75), opt(argv, "--look-in", 0.0),
                                    opt(argv, "--spray-above", 0.0), opt(argv, "--out", "section", str))
        print(f"wrote {len(counts)} section frames to {dst}: {counts}")
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
