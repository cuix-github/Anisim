"""The particle play cache: what the Maya playback compound (particle_play.json) reads,
one numbered pair of .npy files per frame, Bifrost's frame-token convention:

  liquid.####.npy  the liquid particles' positions, (N, 3) float32, in cells
  speed.####.npy   their speeds, (N,) float32, cells per second (the viewport colour)
  cache.json       how many frames, the particle counts, where they came from

Air particles are left out: the viewport shows the water. Maya frame f reads frame f,
starting at 1. write_frame is what a run calls on each dumped frame. Run as a script it
builds a cache without simulating, two ways:

  .venv\\Scripts\\python.exe maya\\play_cache.py states <out_dir> <state_dir> [<state_dir> ...]
      one frame per full-state dump (a chunk's final*.npy)
  .venv\\Scripts\\python.exe maya\\play_cache.py blend <out_dir> <dt> <every> <fps> <dump_dir> [...]
      frames at a video rate from a run's per-substep dumps (frame/vel/ph.####.npy in
      each chunk's directory, substep dt, every Nth substep dumped), each one a cubic
      Hermite blend of the two dumps around its time (positions and velocities both
      known at each end), so a solver that stepped at 0.085 s plays at 24 frames a
      second without anything jumping. This leans on one measured fact: the solver
      keeps its liquid particles in order from one substep to the next (checked on
      the smoke scene: predicting a particle from the dump before lands within 0.05
      cells, and every dump holds the same liquid count), so particle i of one dump
      is particle i of the next. The blend refuses a pair whose liquid counts differ.

Either takes --grid NX NY NZ, recorded in cache.json for the scene builder's domain box.
"""
import json
import os
import sys

import numpy as np


def write_liquid(d, frame, p, s):
    """One frame from liquid positions (N, 3) and speeds (N,), both float32."""
    os.makedirs(d, exist_ok=True)
    np.save(os.path.join(d, f"liquid.{frame:04d}.npy"), np.ascontiguousarray(p, dtype=np.float32))
    np.save(os.path.join(d, f"speed.{frame:04d}.npy"), np.ascontiguousarray(s, dtype=np.float32))
    return int(len(p))


def write_frame(d, frame, pos, vel, phase):
    """One frame of the cache from a run's particles (all phases); returns the liquid
    count. phase >= 0.5 is liquid, as everywhere else in the solver."""
    liquid = np.asarray(phase) >= 0.5
    p = np.asarray(pos, np.float32)[liquid]
    s = np.linalg.norm(np.asarray(vel, np.float32)[liquid], axis=1)
    return write_liquid(d, frame, p, s)


def dump_times(dump_dirs, substeps, every, dt):
    """The run's dumps in order: (time, directory, index), dump k of chunk c being the
    state after substep c * substeps + k, at time (c * substeps + k + 1) * dt."""
    return [((c * substeps + k + 1) * dt, d, k)
            for c, d in enumerate(dump_dirs) for k in range(0, substeps, every)]


def load_liquid(d, k):
    """Liquid positions and velocities of dump k in directory d, float32."""
    ph = np.load(os.path.join(d, f"ph.{k:04d}.npy"))
    liquid = ph >= 0.5
    p = np.load(os.path.join(d, f"frame.{k:04d}.npy"))[liquid].astype(np.float32, copy=False)
    v = np.load(os.path.join(d, f"vel.{k:04d}.npy"))[liquid].astype(np.float32, copy=False)
    return p, v


def hermite(p0, v0, p1, v1, h, u):
    """Cubic Hermite between two states h seconds apart at fraction u: position and
    velocity (the velocity linear, enough for a colour)."""
    u = np.float32(u)
    h00, h10 = 2 * u ** 3 - 3 * u ** 2 + 1, u ** 3 - 2 * u ** 2 + u
    h01, h11 = -2 * u ** 3 + 3 * u ** 2, u ** 3 - u ** 2
    h = np.float32(h)
    p = h00 * p0 + (h10 * h) * v0 + h01 * p1 + (h11 * h) * v1
    return p, (1 - u) * v0 + u * v1


def blend(out_dir, dump_dirs, substeps, every, dt, fps, frames=None, **info):
    """Frames 1..n at `fps` from the dumps, frame 1 at the first dump's time, the last
    at or before the last dump's; `frames` caps the count. Returns the liquid counts."""
    dumps = dump_times(dump_dirs, substeps, every, dt)
    t0, t_last = dumps[0][0], dumps[-1][0]
    n = int(np.floor((t_last - t0) * fps + 1e-6)) + 1
    n = min(n, frames) if frames else n
    counts, pair, cur = [], None, None
    for f in range(1, n + 1):
        t = t0 + (f - 1) / fps
        # the dump pair around t; the last frame may sit exactly on the last dump
        i = min(int(np.searchsorted([x[0] for x in dumps], t, side="right")) - 1, len(dumps) - 2)
        if pair != i:
            a, b = dumps[i], dumps[i + 1]
            pa, va = load_liquid(a[1], a[2]) if cur is None or pair != i - 1 else cur
            pb, vb = load_liquid(b[1], b[2])
            if len(pa) != len(pb):
                raise SystemExit(f"blend: liquid counts differ between dumps {i} and {i + 1} "
                                 f"({len(pa)} vs {len(pb)}); the order assumption fails here")
            pair, cur, h = i, (pb, vb), b[0] - a[0]
        u = (t - dumps[i][0]) / h
        p, v = hermite(pa, va, pb, vb, h, min(max(u, 0.0), 1.0))
        counts.append(write_liquid(out_dir, f, p, np.linalg.norm(v, axis=1)))
    write_index(out_dir, counts, fps=fps, dt_frame_s=round(1.0 / fps, 6), time_s_first=round(t0, 6),
                blended_from=dict(dt=dt, every=every, dumps=len(dumps),
                                  dirs=[os.path.basename(os.path.normpath(d)) for d in dump_dirs]),
                **info)
    return counts


def write_index(d, counts, **info):
    """cache.json: frame count (Maya frames 1..n), liquid counts, and whatever the
    writer wants to record (times, source)."""
    with open(os.path.join(d, "cache.json"), "w") as f:
        json.dump(dict(info, frames=len(counts), liquid_particles=counts), f, indent=1)


def from_states(d, state_dirs, times=None, **info):
    """Frames 1..n from full-state dumps (final.npy, final_velocities.npy,
    final_phase.npy beside it), one per directory; info goes into the index."""
    counts = []
    for f, sd in enumerate(state_dirs, start=1):
        counts.append(write_frame(d, f, np.load(os.path.join(sd, "final.npy")),
                                  np.load(os.path.join(sd, "final_velocities.npy")),
                                  np.load(os.path.join(sd, "final_phase.npy"))))
    write_index(d, counts, source=[os.path.basename(os.path.normpath(s)) for s in state_dirs],
                **({"time_s": list(times)} if times is not None else {}), **info)
    return counts


def main(argv):
    info, a = {}, list(argv)
    if "--grid" in a:   # the grid, for the scene builder's domain box
        i = a.index("--grid")
        info["grid"], a = [int(v) for v in a[i + 1:i + 4]], a[:i] + a[i + 4:]
    if len(a) >= 3 and a[0] == "states":
        print(from_states(a[1], a[2:], **info))
    elif len(a) >= 6 and a[0] == "blend":
        dirs = a[5:]
        substeps = len([f for f in os.listdir(dirs[0]) if f.startswith("frame.")]) * int(a[3])
        print(blend(a[1], dirs, substeps, int(a[3]), float(a[2]), int(a[4]), **info))
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
