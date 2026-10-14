"""The Maya particle playback (maya/particle_play.json): what the viewport is given must
be exactly what the run wrote. Light: a minute, plus a few seconds per real frame.

  small   a synthetic cache, 1,000 particles, 40 percent air: the compound's points are
          the liquid particles' positions bit for bit, in order; their colours follow
          slow + min(speed / speed_scale, 1) * (fast - slow) to float rounding; their
          sizes are point_size; changing speed_scale and point_size moves both
  cache   a real play cache (default out/collider_showcase/play_states, written by
          `collider_showcase.py play` from the run's kept states): the last frame
          through the compound, the same comparison, and the wall per frame

Run from the repo root:  .venv\\Scripts\\python.exe gates\\validate_particle_play.py [small] [cache [<dir>]]
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

OUT = os.path.join(HERE, "..", "out", "play_gate")
BAT = os.path.join(HERE, "..", "maya", "run_particle_play_test.bat")
SLOW = np.array([0.04, 0.22, 0.6], np.float32)       # the compound's defaults
FAST = np.array([0.95, 0.97, 1.0], np.float32)


def play(cache, frame, tag, **ports):
    """The compound at one frame of a cache; returns positions, colours, sizes, wall."""
    c = os.path.abspath(cache).replace("\\", "/")
    o = os.path.abspath(os.path.join(OUT, tag)).replace("\\", "/")
    os.makedirs(o, exist_ok=True)
    args = [BAT, "--set-port", "liquid_pattern", c + "/liquid.####",
            "--set-port", "speed_pattern", c + "/speed.####", "--set-port", "frame", str(frame),
            "--set-port", "path_positions", o + "/pos.npy", "--set-port", "path_colors", o + "/col.npy",
            "--set-port", "path_sizes", o + "/size.npy"]
    for k, v in ports.items():
        args += ["--set-port", k, str(v)]
    t0 = time.time()
    res = subprocess.run(args, capture_output=True, text=True)
    wall = time.time() - t0
    if res.returncode != 0 or "'ok_sizes': (bool) true" not in res.stdout:
        sys.stderr.write(res.stdout[-2000:] + res.stderr[-2000:])
        raise SystemExit(f"bifcmd particle_play failed for {tag}")
    return tuple(np.load(os.path.join(o, n + ".npy")) for n in ("pos", "col", "size")) + (wall,)


def compare(cache, frame, got, speed_scale=25.0, point_size=0.3):
    pos, col, size, _ = got
    want = np.load(os.path.join(cache, f"liquid.{frame:04d}.npy"))
    spd = np.load(os.path.join(cache, f"speed.{frame:04d}.npy"))
    f = np.minimum(spd / np.float32(speed_scale), 1.0)[:, None]
    want_col = SLOW + f * (FAST - SLOW)
    same = pos.shape == want.shape and np.array_equal(pos, want)
    cerr = float(np.abs(col - want_col).max()) if col.shape == want_col.shape else np.inf
    serr = float(np.abs(size - point_size).max()) if size.shape == (len(want),) else np.inf
    return same, cerr, serr


def small_cache(n=1000):
    """Frame 2 the one under test, frame 3 a neighbour that must not be read."""
    cache = os.path.join(OUT, "small_cache")
    rng = np.random.default_rng(0)
    pos = rng.uniform(0, 30, (n, 3)).astype(np.float32)
    vel = rng.normal(0, 40, (n, 3)).astype(np.float32)
    ph = (rng.uniform(size=n) > 0.4).astype(np.float32)
    liquid = write_frame(cache, 2, pos, vel, ph)
    write_frame(cache, 3, pos[::-1], vel, ph)
    return cache, liquid


def check_small():
    print("=== small: a synthetic cache through the compound")
    n = 1000
    cache, liquid = small_cache(n)
    ok = True
    for tag, ports in (("defaults", {}), ("ports", {"speed_scale": 20.0, "point_size": 0.5})):
        got = play(cache, 2, "small_" + tag, **ports)
        same, cerr, serr = compare(cache, 2, got, **{k: v for k, v in ports.items()})
        good = same and cerr <= 1e-6 and serr == 0.0
        ok &= good
        print(f"    {tag:8s} {liquid} liquid of {n}: positions identical {same}; colour max |err|"
              f" {cerr:.1e}; size max |err| {serr:.1e}  {'OK' if good else 'FAIL'}")
    return ok


def check_cache(cache=None):
    cache = cache or os.path.join(HERE, "..", "out", "collider_showcase", "play_states")
    print(f"=== cache: the last frame of {os.path.normpath(cache)}")
    if not os.path.exists(os.path.join(cache, "cache.json")):
        print("    no cache there (write one first, e.g. collider_showcase.py play)  FAIL")
        return False
    with open(os.path.join(cache, "cache.json")) as f:
        info = json.load(f)
    frame = info["frames"]
    got = play(cache, frame, "cache")
    same, cerr, serr = compare(cache, frame, got)
    _, _, _, small = play(small_cache()[0], 2, "startup")
    good = same and cerr <= 1e-6 and serr == 0.0
    print(f"    frame {frame}: {len(got[0]):,} points; positions identical {same}; colour max"
          f" |err| {cerr:.1e}; size max |err| {serr:.1e}  {'OK' if good else 'FAIL'}")
    print(f"    bifcmd wall {got[3]:.1f} s against {small:.1f} s for 600 points (startup):"
          f" about {got[3] - small:.1f} s to read and build the frame")
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
    print("PARTICLE PLAY:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
