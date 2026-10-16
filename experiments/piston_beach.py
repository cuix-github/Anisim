"""The first breaker from the wavemaker: a piston plate (step_3d's piston ports) makes a
wave train in 12 m of water that runs onto a 1:10 beach and plunges. The scene, in
metres: the plate 1 m thick at x = 2 m on the floor, its top at 18 m, a 4 m stroke at
a period of 8 s (linear theory: L = 76 m, H = 3.9 m, the surf-similarity number
tan(beta) / sqrt(H / L0) = 0.51, a plunging breaker); the beach's toe at 36 m, the
shoreline at 156 m, the domain 172 x 24 x 16 m. Production configuration as the
showcase (ST, escaped particles, adaptive air, multigrid, the fused transfer).

  ..\\.venv\\Scripts\\python.exe piston_beach.py preview [--cpm M] [--seconds S]
      one cell wide (nz = 1) at M cells per metre (4): the whole run in minutes, a dump
      every 0.1 s, the leading crest's log and its marks (the face steep, the lip over),
      the side-view audit frames. Finds the break before the 3D run.
  ..\\.venv\\Scripts\\python.exe piston_beach.py chunk C --plan "S1:E1,S2:E2,..." [--cpm M]
      chunk C of the 3D run: the plan lists each chunk's seconds and dump interval in
      seconds; chunk C resumes from chunk C - 1's state (the piston's clock carries on).
  ..\\.venv\\Scripts\\python.exe piston_beach.py measure [--cpm M]
      the leading crest's log over every dump of the 3D run so far, its marks, audit frames.
  ..\\.venv\\Scripts\\python.exe piston_beach.py cache T0 T1 [--cpm M]
      the play cache (liquid and speed) for the dumps between T0 and T1 seconds, for the
      surface render (maya/render_surface_preview.bat), and the beach as an OBJ beside
      the run's inputs.

Everything lands under out/runs/<tag>/ (graphs/runs.py): tag pistonbeach<M> for the 3D
run, pistonbeach<M>x1 for the preview.
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
from solitary_beach import calibration, close_up, face, side_view, surface_xz, write_beach_obj  # noqa: E402
from solitary_wave import beach_bed  # noqa: E402
from loop3d import DUMPS, run_bifrost3, seed3d, load_final_state  # noqa: E402
from play_cache import write_frame, write_index  # noqa: E402
from runs import RunPaths, prepare  # noqa: E402
from PIL import Image  # noqa: E402

# the scene in metres
H, X0, W, TOP, TOE, SLOPE = 12.0, 2.0, 1.0, 18.0, 36.0, 0.1
LX, LY, LZ, S, T = 172.0, 24.0, 16.0, 4.0, 8.0
SHORE = TOE + H / SLOPE


class Run:
    def __init__(self, cpm=4.0, flat=False):
        self.cpm = cpm
        self.grid = (int(LX * cpm), int(LY * cpm), 1 if flat else int(LZ * cpm))
        self.dt = 0.1 / cpm
        self.tag = f"pistonbeach{cpm:g}" + ("x1" if flat else "")


def scene(run):
    m = run.cpm
    return {"m": m, "h": H * m, "toe": TOE * m, "slope": SLOPE, "shore": SHORE * m, "g": 9.8 * m,
            "x0": X0 * m, "w": W * m, "top": TOP * m}


def ports(run, sc):
    return dict(calibration(), gravity=-sc["g"], piston_stroke=S * sc["m"], piston_period=T,
                piston_x0=sc["x0"], piston_width=sc["w"], piston_top=sc["top"])


def beach_phi(run, sc):
    nx, ny, nz = run.grid
    i, j = np.meshgrid(np.arange(nx + 1, dtype=float), np.arange(ny + 1, dtype=float), indexing="ij")
    phi2 = (j - beach_bed(i, sc["toe"], sc["slope"])) * np.cos(np.arctan(sc["slope"]))
    return np.repeat(phi2[:, :, None], nz + 1, axis=2).astype(np.float32)


def seed(run, sc):
    """Water to h in front of the plate, air elsewhere; nothing inside the plate or
    the sand."""
    nx, ny, nz = run.grid
    pos, typ = seed3d(nx, ny, nz, lambda x, y, z: (x > sc["x0"] + sc["w"]) & (y < sc["h"]))
    if nz == 1:
        pos[:, 2] = np.clip(pos[:, 2], 0.0, 1.0)
    in_plate = (pos[:, 0] > sc["x0"]) & (pos[:, 0] < sc["x0"] + sc["w"]) & (pos[:, 1] < sc["top"])
    keep = (pos[:, 1] > beach_bed(pos[:, 0], sc["toe"], sc["slope"])) & ~in_plate
    pos, typ = pos[keep], typ[keep]
    n = len(pos)
    zeros = np.zeros(n, np.float32)
    return {"positions": pos.astype(np.float32), "velocities": np.zeros((n, 3), np.float32),
            "particle_phase": typ.astype(np.int32), "mass": np.where(typ == 1, 1000.0, 1.0).astype(np.float32),
            "momx": zeros, "momy": zeros, "momz": zeros, "tau": zeros, "wt": np.ones(n, np.float32),
            "wtph": typ.astype(np.float32), "scale": np.ones(n, np.float32), "phase_st": typ.astype(np.float32)}


# ------------------------------------------------------------------ measuring


def leading_crest(top, sc):
    """The rightmost crest over a metre above still water, short of the swash."""
    m = sc["m"]
    eta = np.where(np.isfinite(top), top - sc["h"], -1e9)
    k = np.ones(int(2 * m) + 1) / (int(2 * m) + 1)
    sm = np.convolve(np.where(eta > -1e8, eta, 0.0), k, "same")
    i = np.arange(1, len(sm) - 1)
    cand = i[(sm[i] > 1.0 * m) & (sm[i] >= sm[i - 1]) & (sm[i] >= sm[i + 1]) &
             (i > 20 * m) & (i < sc["shore"] - 8 * m)]
    if len(cand) == 0:
        return None, None
    c = int(cand.max())
    lo, hi = max(0, c - int(3 * m)), c + int(3 * m) + 1
    j = lo + int(np.argmax(np.where(np.isfinite(top[lo:hi]), top[lo:hi], -1e9)))
    return j + 0.5, float(top[j])


def measure(pos, ph, esc, run, sc, t):
    nx, ny, nz = run.grid
    m = sc["m"]
    liquid = ph >= 0.5
    top = surface_xz(pos, liquid, nx, nz)
    x, y = leading_crest(top, sc)
    spray = int((liquid & (esc > 0.5)).sum()) if esc is not None else 0
    rec = {"t": round(t, 3), "spray_particles": spray, "liquid": int(liquid.sum())}
    if x is None:
        return rec
    steep, over = face(pos, liquid, sc, x, nz)
    rec.update({"crest_x_m": round(x / m, 2), "crest_height_m": round((y - sc["h"]) / m, 2),
                "depth_under_crest_m": round((sc["h"] - float(beach_bed(x, sc["toe"], sc["slope"]))) / m, 2),
                "face_steepness": round(steep, 2), "overhang_m": round(over / m, 2)})
    return rec


def marks(log):
    steep = next((e for e in log if e.get("face_steepness", 0) >= 2.5 and e.get("crest_height_m", 0) >= 1.0), None)
    over = next((e for e in log if e.get("overhang_m", 0) >= 1.5), None)
    return steep, over


def dumps_of(run, paths):
    """Every dump of the run so far: (time, directory, index), from the chunk receipts."""
    out, offset, c = [], 0, 0
    while os.path.exists(paths.receipt(f"c{c}.json")):
        r = json.load(open(paths.receipt(f"c{c}.json")))
        d = paths.dumps(c)
        out += [((offset + k + 1) * run.dt, d, k) for k in range(0, r["substeps"], r["every"])]
        offset += r["substeps"]
        c += 1
    return out


def audit(run, sc, paths, entries, log):
    """Side views of every dump, and the close-ups around the marks."""
    frames = paths.renders("audit")
    os.makedirs(frames, exist_ok=True)
    for n, (t, d, k) in enumerate(entries):
        p = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
        ph = np.load(os.path.join(d, f"ph.{k:04d}.npy"))
        side_view(p, ph, run, sc, os.path.join(frames, f"fr{n:04d}.png"), scale=max(1, int(4 / run.cpm)))
    steep, over = marks(log)
    centre = over or steep
    if centre is None:
        return None
    i = log.index(centre)
    shots = []
    for j in range(max(0, i - 6), min(len(log), i + 10), 2):
        t, d, k = entries[j]
        p = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
        ph = np.load(os.path.join(d, f"ph.{k:04d}.npy"))
        cx = (log[j].get("crest_x_m") or centre["crest_x_m"]) * sc["m"]
        shots.append(close_up(p, ph, run, cx, 30 * sc["m"], scale=max(1, int(round(6 / sc["m"])))))
    w, h = shots[0].size
    sheet = Image.new("RGB", (w, h * len(shots)), "white")
    for j, im in enumerate(shots):
        sheet.paste(im, (0, j * h))
    path = prepare(paths.receipt("breaker.png"))
    sheet.save(path)
    return path


def report(run, sc, paths, entries):
    log = []
    for t, d, k in entries:
        p = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
        ph = np.load(os.path.join(d, f"ph.{k:04d}.npy"))
        e = os.path.join(d, f"esc.{k:04d}.npy")
        log.append(measure(p, ph, np.load(e) if os.path.exists(e) else None, run, sc, t))
    steep, over = marks(log)
    strip = audit(run, sc, paths, entries, log)
    summary = {"tag": run.tag, "grid": list(run.grid), "cells_per_metre": run.cpm, "dt": run.dt,
               "scene_m": {"h": H, "plate": [X0, W, TOP], "stroke": S, "period": T, "toe": TOE, "slope": SLOPE,
                           "shore": SHORE, "domain": [LX, LY, LZ]},
               "steep_front": steep, "turned_over": over, "breaker_closeups": strip, "log": log}
    with open(prepare(paths.receipt("summary.json")), "w") as f:
        json.dump(summary, f, indent=1)
    for e in log[::5]:
        print("  ", json.dumps(e))
    print("steep:", json.dumps(steep))
    print("over: ", json.dumps(over))
    return summary


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
                              dumps=DUMPS + ("vel",), fused=True, out_dir=paths.dumps(chunk))
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
    cpm = opt(argv, "--cpm", 4.0)
    if cmd == "preview":
        run = Run(cpm, flat=True)
        sc = scene(run)
        seconds = opt(argv, "--seconds", 24.0)
        run_chunk(run, 0, seconds, 0.1, 0, seed(run, sc))
        report(run, sc, RunPaths(run.tag), dumps_of(run, RunPaths(run.tag)))
    elif cmd == "chunk":
        run = Run(cpm)
        sc = scene(run)
        c = int(argv[1])
        plan = [tuple(float(v) for v in p.split(":")) for p in opt(argv, "--plan", "", str).split(",")]
        paths = RunPaths(run.tag)
        offset = sum(int(round(s / run.dt)) for s, _ in plan[:c])
        init = load_final_state(paths.dumps(c - 1)) if c else seed(run, sc)
        run_chunk(run, c, plan[c][0], plan[c][1], offset, init)
    elif cmd == "measure":
        run = Run(cpm)
        report(run, scene(run), RunPaths(run.tag), dumps_of(run, RunPaths(run.tag)))
    elif cmd == "cache":
        run = Run(cpm, flat="--flat" in argv)
        paths = RunPaths(run.tag)
        t0, t1 = float(argv[1]), float(argv[2])
        cache = paths.cache("full")
        os.makedirs(cache, exist_ok=True)
        counts, times = [], []
        for t, d, k in dumps_of(run, paths):
            if t0 - 1e-9 <= t <= t1 + 1e-9:
                counts.append(write_frame(cache, len(counts) + 1, np.load(os.path.join(d, f"frame.{k:04d}.npy")),
                                          np.load(os.path.join(d, f"vel.{k:04d}.npy")),
                                          np.load(os.path.join(d, f"ph.{k:04d}.npy"))))
                times.append(round(t, 3))
        write_index(cache, counts, grid=list(run.grid), times_s=times, source=run.tag)
        write_beach_obj(prepare(paths.input("beach.obj")), run, scene(run))
        print(f"wrote {len(counts)} frames to {cache}: t {times[0] if times else '-'} .. {times[-1] if times else '-'}")
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
