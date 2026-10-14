"""The particle play cache: what the Maya playback compound (particle_play.json) reads,
one numbered pair of .npy files per frame, Bifrost's frame-token convention:

  liquid.####.npy  the liquid particles' positions, (N, 3) float32, in cells
  speed.####.npy   their speeds, (N,) float32, cells per second (the viewport colour)
  cache.json       how many frames, the particle counts, where they came from

Air particles are left out: the viewport shows the water. Maya frame f reads frame f,
starting at 1. write_frame is what a run calls on each dumped frame; run as a script it
turns full-state dumps (a chunk's final*.npy) into frames, no simulation needed:

  .venv\\Scripts\\python.exe maya\\play_cache.py <out_dir> <state_dir> [<state_dir> ...]
"""
import json
import os
import sys

import numpy as np


def write_frame(d, frame, pos, vel, phase):
    """One frame of the cache from a run's particles (all phases); returns the liquid
    count. phase >= 0.5 is liquid, as everywhere else in the solver."""
    os.makedirs(d, exist_ok=True)
    liquid = np.asarray(phase) >= 0.5
    p = np.ascontiguousarray(np.asarray(pos, np.float32)[liquid])
    s = np.linalg.norm(np.asarray(vel, np.float32)[liquid], axis=1).astype(np.float32)
    np.save(os.path.join(d, f"liquid.{frame:04d}.npy"), p)
    np.save(os.path.join(d, f"speed.{frame:04d}.npy"), s)
    return int(liquid.sum())


def write_index(d, counts, **info):
    """cache.json: frame count (Maya frames 1..n), liquid counts, and whatever the
    writer wants to record (times, source)."""
    with open(os.path.join(d, "cache.json"), "w") as f:
        json.dump(dict(info, frames=len(counts), liquid_particles=counts), f, indent=1)


def from_states(d, state_dirs, times=None):
    """Frames 1..n from full-state dumps (final.npy, final_velocities.npy,
    final_phase.npy beside it), one per directory."""
    counts = []
    for f, sd in enumerate(state_dirs, start=1):
        counts.append(write_frame(d, f, np.load(os.path.join(sd, "final.npy")),
                                  np.load(os.path.join(sd, "final_velocities.npy")),
                                  np.load(os.path.join(sd, "final_phase.npy"))))
    write_index(d, counts, source=[os.path.basename(os.path.normpath(s)) for s in state_dirs],
                **({"time_s": list(times)} if times is not None else {}))
    return counts


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    print(from_states(sys.argv[1], sys.argv[2:]))
