"""A solitary wave onto a plane beach in 3D: the 2D control (reference2d/solitary_wave.py)
run through the Bifrost loop, one variable at a time. The same scene in metres: 20 m of
still water, a 6 m hump (H/h = 0.3) started at x = 80 m with the solitary wave's
velocity field on its liquid particles, a 1:15 beach from x = 160 m to the shoreline at
x = 460 m as a cut-cell level set, the width a slice in z with free-slip walls. Cells of
one metre at the default grid, so the wave is six cells high, which is what the ocean
plan asks of the solver and what this run measures: where the front face turns steep,
where it turns over, against the 2D reference and Grilli's breaking depth.

  ..\\.venv\\Scripts\\python.exe solitary_beach.py smoke          # a tiny grid through the graph, ~1 min
  ..\\.venv\\Scripts\\python.exe solitary_beach.py stage 0        # one chunk; 1, 2, ... resume
  ..\\.venv\\Scripts\\python.exe solitary_beach.py finalize       # the wave's log, verdict, preview video
  ..\\.venv\\Scripts\\python.exe solitary_beach.py night          # every chunk, then finalize
  ..\\.venv\\Scripts\\python.exe solitary_beach.py remeasure      # the log again from the dumps, then finalize
  (finalize and remeasure take --record: the summary and the breaker close-ups go to results/3d/solitary_beach/)

Options after any command: --grid NX NY NZ (520 56 32: a 32 m slice), --substeps N a
chunk (65), --chunks N (4; 26 s at dt 0.1), --every N dump every Nth substep (5, half a
second: the Hermite blend bridges it), --dt S (0.1; CFL 1.6 at the wave's 16 m/s, which
ST-FLIP takes), --tag NAME (beach, or the grid's), --stock (the stock-splat graph),
--cells-per-metre C (1; 2 doubles every length in cells and g, the same seconds).
Production config as the showcase: ST, escaped particles, adaptive air, multigrid, the
fused transfer. Everything lands under out/runs/<tag>/ (graphs/runs.py): the dumps, the
full play cache (every dump), the beach as an OBJ for Maya, the receipts, the side-view
audit frames and their video. A 520 x 56 x 32 run is about 7.5 M particles, 270 MB a
dump, 14 GB in all at every 5.

The measurements, per dump, from the liquid particles projected along z: the crest's
position and height (the median over z of each x column's highest liquid particle)
and the still-water depth under it; then the front face as the liquid's front in
one-metre bands above still water, its steepness (the largest rise over advance
across three bands) and its overhang (the largest excess of a band's front over a
lower band's). The first dump with the face steeper than 2.5 under a crest shoaled to 1.1 H, and
the first with an overhang over 1.5 m, the crest between the toe and the shoreline,
are the marks the 2D control sets
(reference2d/solitary_wave.py, its receipts in results/2d/solitary_wave/), to be read
against Grilli's vertical face at H_b/h_b 2.76 for this slope, about 2.6 m of water.
"""
import json
import os
import subprocess
import sys
import time

import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
sys.path.insert(0, os.path.join(HERE, "..", "reference2d"))
sys.path.insert(0, os.path.join(HERE, "..", "maya"))
from pfflip2d import Sim  # noqa: E402
from solitary_wave import profile, velocity, beach_bed, G  # noqa: E402
from loop3d import DUMPS, run_bifrost3, seed3d, load_final_state  # noqa: E402
from play_cache import write_frame, write_index  # noqa: E402
from runs import RunPaths, prepare  # noqa: E402

PORTS = {"preconditioner": 1, "st": 1, "escape": 1, "adapt": 1}
# the scene in metres: depth, wave height, the hump's start, the beach's toe and slope
H_WATER, H_WAVE, X0, TOE, SLOPE, SHORE = 20.0, 6.0, 80.0, 160.0, 1.0 / 15.0, 460.0
FULL = ((520, 56, 32), 65, 4, 5)          # grid, substeps a chunk, chunks, dump every
SMOKE = ((104, 16, 8), 8, 2, 2)           # a fifth of the scene in x and y, 5 m cells


class Run:
    def __init__(self, grid, substeps, chunks, every, dt=0.1, tag=None, fused=True, cpm=1.0):
        self.grid, self.substeps, self.chunks, self.every, self.dt = tuple(grid), substeps, chunks, every, dt
        self.fused, self.cpm = fused, cpm
        self.tag = tag or ("beach" if self.grid == FULL[0] else f"beach{grid[0]}x{grid[2]}")

    def __repr__(self):
        return (f"Run(grid={self.grid}, substeps={self.substeps}, chunks={self.chunks}, every={self.every}, "
                f"dt={self.dt}, tag={self.tag!r}, fused={self.fused}, cells_per_metre={self.cpm})")


def parse(argv):
    grid, sub, n, every = FULL
    dt, tag, fused, cpm, words, i = 0.1, None, True, 1.0, [], 0
    while i < len(argv):
        a = argv[i]
        if a == "--grid":
            grid, i = tuple(int(v) for v in argv[i + 1:i + 4]), i + 4
        elif a == "--substeps":
            sub, i = int(argv[i + 1]), i + 2
        elif a == "--chunks":
            n, i = int(argv[i + 1]), i + 2
        elif a == "--every":
            every, i = int(argv[i + 1]), i + 2
        elif a == "--dt":
            dt, i = float(argv[i + 1]), i + 2
        elif a == "--tag":
            tag, i = argv[i + 1], i + 2
        elif a == "--stock":
            fused, i = False, i + 1
        elif a == "--cells-per-metre":
            cpm, i = float(argv[i + 1]), i + 2
        elif a == "--record":
            i += 1
        elif a.startswith("--"):
            raise SystemExit(f"unknown option {a}\n" + __doc__)
        else:
            words, i = words + [a], i + 1
    return Run(grid, sub, n, every, dt, tag, fused, cpm), words


def scene(run):
    """The scene in cells: the grid's x spans the 520 m whether it has 520 cells or
    104, so a smoke grid is the same beach at coarser cells."""
    nx, ny, nz = run.grid
    m = nx / 520.0                            # cells per metre along x (the smoke's 0.2)
    return {"m": m, "h": H_WATER * m, "H": H_WAVE * m, "x0": X0 * m, "toe": TOE * m,
            "slope": SLOPE, "shore": SHORE * m, "g": G * m}


def beach_phi(nx, ny, nz, sc):
    i, j = np.meshgrid(np.arange(nx + 1, dtype=float), np.arange(ny + 1, dtype=float), indexing="ij")
    phi2 = (j - beach_bed(i, sc["toe"], sc["slope"])) * np.cos(np.arctan(sc["slope"]))
    return np.repeat(phi2[:, :, None], nz + 1, axis=2).astype(np.float32)


def seed_wave(run, sc):
    """The particles: liquid under the hump's surface and above the sand, air above;
    the wave's velocity on the liquid. Returns init_arrays for run_bifrost3."""
    nx, ny, nz = run.grid
    eta_of = lambda x: profile(x, sc["h"], sc["H"], sc["x0"])[0]   # noqa: E731
    pos, typ = seed3d(nx, ny, nz, lambda x, y, z: y < sc["h"] + eta_of(x))
    keep = pos[:, 1] > beach_bed(pos[:, 0], sc["toe"], sc["slope"])
    pos, typ = pos[keep], typ[keep]
    vel = np.zeros_like(pos, dtype=np.float32)
    liquid = typ == 1
    u, w, c = velocity(pos[liquid, 0], pos[liquid, 1], sc["h"], sc["H"], sc["x0"], sc["g"])
    vel[liquid, 0], vel[liquid, 1] = u, w
    n = len(pos)
    zeros = np.zeros(n, np.float32)
    return {"positions": pos.astype(np.float32), "velocities": vel, "particle_phase": typ.astype(np.int32),
            "mass": np.where(typ == 1, 1000.0, 1.0).astype(np.float32),
            "momx": zeros, "momy": zeros, "momz": zeros, "tau": zeros, "wt": np.ones(n, np.float32),
            "wtph": typ.astype(np.float32), "scale": np.ones(n, np.float32),
            "phase_st": typ.astype(np.float32)}, float(c)


def calibration():
    cal = Sim(64, 64, adapt=True, st=True)
    cal.seed(lambda x, y: np.ones_like(x, dtype=bool))
    cal.calibrate()
    return dict(PORTS, rho0_face=cal.rho0_face, coarse_gain=cal.coarse_gain, ws_epsilon=16)


# ------------------------------------------------------------------ measuring


def surface_xz(pos, liquid, nx, nz):
    """The highest liquid particle of each (x, z) column, then the median over z per
    x column; -inf where a column holds no liquid."""
    p = pos[liquid]
    ci = np.clip(p[:, 0].astype(np.int64), 0, nx - 1)
    ck = np.clip(p[:, 2].astype(np.int64), 0, nz - 1)
    top = np.full((nx, nz), -np.inf)
    np.maximum.at(top, (ci, ck), p[:, 1])
    wet = np.isfinite(top)
    med = np.full(nx, -np.inf)
    rows = wet.sum(1) > nz // 2
    med[rows] = np.array([np.median(top[i][wet[i]]) for i in np.flatnonzero(rows)])
    return med


def crest_of(top, lo):
    seg = top[lo:]
    i = int(np.argmax(seg))
    x = lo + i + 0.5
    if 0 < i < len(seg) - 1 and np.all(np.isfinite(seg[i - 1:i + 2])):
        a, b, c = seg[i - 1], seg[i], seg[i + 1]
        den = a - 2 * b + c
        if den < 0:
            x += 0.5 * (a - c) / den
    return float(x), float(seg[i])


def face(pos, liquid, sc, crest_x, nz):
    """The wave's front face, the 2D reference's rule over every z: the liquid's front
    in one-cell bands above still water, 10 m behind the crest to 40 m ahead, a band's
    front being the x with one cell's worth of particles beyond it across the width
    (so a stray particle does not move it and a long band is not read short); a band
    with fewer than five cells' worth (the tip) ends the count. Returns the steepness,
    the largest rise over advance across any three bands (cells per cell, 99 when
    vertical or beyond), and the overhang, the largest excess of a band's front over
    any lower band's (zero for a face that leans forward without turning over)."""
    m = sc["m"]
    p = pos[liquid]
    p = p[(p[:, 0] > crest_x - 10 * m) & (p[:, 0] < crest_x + 40 * m) & (p[:, 1] >= sc["h"])]
    if len(p) < 20:
        return 0.0, 0.0
    band = ((p[:, 1] - sc["h"]) / m).astype(np.int64)
    k = int(8 * m * nz)                 # the particles of one cell of a band, across the width
    fronts = []
    for b in range(int(band.max()) + 1):
        xs = np.sort(p[band == b, 0])
        if len(xs) < 5 * k:
            break
        fronts.append(float(xs[-k]))    # a cell's worth of particles stands behind it
    steep = 0.0
    for j in range(len(fronts) - 3):
        dx = fronts[j] - fronts[j + 3]
        steep = max(steep, 99.0 if dx <= 1e-6 else min(99.0, 3.0 * m / dx))
    over = max([fronts[j] - min(fronts[:j]) for j in range(1, len(fronts))] + [0.0])
    return steep, over


def measure(pos, ph, run, sc, t):
    nx, ny, nz = run.grid
    liquid = ph >= 0.5
    top = surface_xz(pos, liquid, nx, nz)
    x, y = crest_of(top, int(20 * sc["m"]))
    m = sc["m"]
    steepness, over = face(pos, liquid, sc, x, nz)
    return {"t": round(t, 2), "crest_x_m": round(x / m, 2), "crest_height_m": round((y - sc["h"]) / m, 2),
            "depth_under_crest_m": round((sc["h"] - float(beach_bed(x, sc["toe"], sc["slope"]))) / m, 2),
            "face_steepness": round(steepness, 2), "overhang_m": round(over / m, 2), "liquid": int(liquid.sum())}


def side_view(pos, ph, run, sc, path, scale=2):
    """The particles of the middle z slab, projected along z, liquid blue over air grey."""
    nx, ny, nz = run.grid
    img = np.full((nx * scale, ny * scale, 3), 245, np.uint8)
    slab = np.abs(pos[:, 2] - nz / 2) < 1.0
    for sel, col in (((ph < 0.5) & slab, (210, 210, 215)), ((ph >= 0.5) & slab, (30, 90, 200))):
        pix = (pos[sel, :2] * scale).astype(np.int64)
        pix[:, 0] = np.clip(pix[:, 0], 0, nx * scale - 1)
        pix[:, 1] = np.clip(pix[:, 1], 0, ny * scale - 1)
        img[pix[:, 0], pix[:, 1]] = col
    Image.fromarray(np.rot90(img)).save(path)


def close_up(pos, ph, run, cx, half, scale=6):
    """The middle z slab around column cx, `half` cells either side, the whole height."""
    nx, ny, nz = run.grid
    x0, x1 = int(max(0, cx - half)), int(min(nx, cx + half))
    img = np.full(((x1 - x0) * scale, ny * scale, 3), 245, np.uint8)
    slab = (np.abs(pos[:, 2] - nz / 2) < 1.0) & (pos[:, 0] >= x0) & (pos[:, 0] < x1)
    for sel, col in (((ph < 0.5) & slab, (210, 210, 215)), ((ph >= 0.5) & slab, (30, 90, 200))):
        pix = ((pos[sel, :2] - (x0, 0)) * scale).astype(np.int64)
        pix[:, 0] = np.clip(pix[:, 0], 0, img.shape[0] - 1)
        pix[:, 1] = np.clip(pix[:, 1], 0, img.shape[1] - 1)
        img[pix[:, 0], pix[:, 1]] = col
    return Image.fromarray(np.rot90(img))


def breaker_strip(run, sc, paths, log, centre, path, around=(-2, -1, 0, 1, 2)):
    """Close-ups of the dumps around the break (one dump apart), the crest centred."""
    frames = []
    for off in around:
        i = centre + off
        if not 0 <= i < len(log):
            continue
        substep = i * run.every
        chunk, k = substep // run.substeps, substep % run.substeps
        d = paths.dumps(chunk)
        p = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
        ph = np.load(os.path.join(d, f"ph.{k:04d}.npy"))
        frames.append(close_up(p, ph, run, log[i]["crest_x_m"] * sc["m"], 40 * sc["m"],
                               scale=max(1, int(round(6 / sc["m"])))))
    if frames:
        w, h = frames[0].size
        sheet = Image.new("RGB", (w, h * len(frames)), "white")
        for j, im in enumerate(frames):
            sheet.paste(im, (0, j * h))
        sheet.save(path)
    return path


def write_beach_obj(path, run, sc):
    """The beach as a closed wedge mesh in cells, for the Maya scene."""
    nx, ny, nz = run.grid
    x1, y1 = float(nx), float((nx - sc["toe"]) * sc["slope"])
    v = [(sc["toe"], 0, 0), (x1, 0, 0), (x1, y1, 0), (sc["toe"], 0, nz), (x1, 0, nz), (x1, y1, nz)]
    faces = [(1, 3, 2), (4, 5, 6), (1, 2, 5, 4), (2, 3, 6, 5), (3, 1, 4, 6)]
    with open(path, "w") as f:
        f.write("o beach\n")
        f.writelines(f"v {x:.4f} {y:.4f} {z:.4f}\n" for x, y, z in v)
        f.writelines("f " + " ".join(str(i) for i in fc) + "\n" for fc in faces)


# ------------------------------------------------------------------ the run


def run_chunk(run, chunk):
    nx, ny, nz = run.grid
    sc = scene(run)
    paths = RunPaths(run.tag)
    phi = beach_phi(nx, ny, nz, sc)
    ports = dict(calibration(), gravity=-sc["g"])      # cells/s^2: 9.8 at one-metre cells
    if chunk:
        init, c = load_final_state(paths.dumps(chunk - 1)), None
    else:
        init, c = seed_wave(run, sc)
    print(f"{run.tag} chunk {chunk}: {run}; {len(init['positions']):,} particles"
          + (f"; c = {c:.2f} cells/s" if c else ""), flush=True)
    t0 = time.time()
    d, _, diag = run_bifrost3(f"{run.tag}_c{chunk}", nx, ny, nz, None, run.substeps, run.dt,
                              extra_ports=ports, init_arrays=init, final_state=True,
                              step_offset=chunk * run.substeps, dump_every=run.every, solid=(phi, 2),
                              dumps=DUMPS + ("vel",), fused=run.fused, out_dir=paths.dumps(chunk))
    wall = time.time() - t0
    frames_dir = paths.renders("audit")
    os.makedirs(frames_dir, exist_ok=True)
    log, play = [], []
    for k in range(0, run.substeps, run.every):
        p = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
        ph = np.load(os.path.join(d, f"ph.{k:04d}.npy"))
        t = (chunk * run.substeps + k + 1) * run.dt
        log.append(measure(p, ph, run, sc, t))
        idx = (chunk * run.substeps + k) // run.every
        side_view(p, ph, run, sc, os.path.join(frames_dir, f"fr{idx:04d}.png"))
        vel = np.load(os.path.join(d, f"vel.{k:04d}.npy"))
        play.append(write_frame(paths.cache("full"), idx + 1, p, vel, ph))
    entry = {"chunk": chunk, "substeps": run.substeps, "every": run.every, "dt": run.dt, "grid": list(run.grid),
             "fused": run.fused, "wall_s": round(wall, 1), "wall_per_substep_s": round(wall / run.substeps, 2),
             "sim_time_s": [round(chunk * run.substeps * run.dt, 2), round((chunk + 1) * run.substeps * run.dt, 2)],
             "particles": [e["liquid"] for e in log], "play_liquid": play,
             "iterations_mean": round(float(np.mean(diag[:, 0])), 1),
             "peak_speed": round(float(max(diag[k][-1] for k in range(run.substeps))), 2), "wave": log}
    with open(prepare(paths.receipt(f"c{chunk}.json")), "w") as f:
        json.dump(entry, f, indent=1)
    print(json.dumps({k: v for k, v in entry.items() if k != "wave"}), flush=True)
    print("  wave:", json.dumps(log[-1]), flush=True)
    return entry


def finalize(run, record=False):
    """The verdict over every chunk, the audit video, the play cache's index, the beach
    mesh, and close-ups around the break; with record, the summary and the close-ups
    are copied into results/3d/solitary_beach/ as the committed evidence."""
    sc = scene(run)
    paths = RunPaths(run.tag)
    chunks = []
    for c in range(run.chunks):
        with open(paths.receipt(f"c{c}.json")) as f:
            chunks.append(json.load(f))
    log = [e for c in chunks for e in c["wave"]]
    on_beach = [e for e in log if TOE < e["crest_x_m"] < SHORE]
    steep = next((e for e in on_beach if e["face_steepness"] >= 2.5 and e["crest_height_m"] >= 1.1 * H_WAVE), None)
    broke = next((e for e in on_beach if e["overhang_m"] >= 1.5), None)
    mark = broke or steep
    zoom = None
    if mark is not None:
        centre = next(i for i, e in enumerate(log) if e["t"] == mark["t"])
        zoom = breaker_strip(run, sc, paths, log, centre, os.path.join(paths.renders("audit"), "breaker.png"))
    S0 = 1.521 * SLOPE / np.sqrt(H_WAVE / H_WATER)
    frames_dir = paths.renders("audit")
    video = prepare(paths.video("audit"))
    fps = max(1, round(1.0 / (run.every * run.dt)))
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i",
                    os.path.join(frames_dir, "fr%04d.png"), "-c:v", "libx264", "-pix_fmt", "yuv420p",
                    "-crf", "21", video], check=True)
    play = [n for c in chunks for n in c["play_liquid"]]
    write_index(paths.cache("full"), play, dt_frame_s=round(run.every * run.dt, 4),
                time_s_first=round(run.dt, 4), source=run.tag, grid=list(run.grid))
    write_beach_obj(prepare(paths.input("colliders.obj")), run, sc)
    summary = {"tag": run.tag, "grid": list(run.grid), "dt": run.dt, "every": run.every, "fused": run.fused,
               "config": PORTS, "cells_per_metre": sc["m"],
               "simulated_s": round(sum(c["substeps"] for c in chunks) * run.dt, 2),
               "wall_s_total": round(sum(c["wall_s"] for c in chunks), 1),
               "wall_per_substep_s": [c["wall_per_substep_s"] for c in chunks],
               "peak_speed_cells_per_s": max(c["peak_speed"] for c in chunks),
               "iterations_mean": [c["iterations_mean"] for c in chunks],
               "steep_front": steep, "turned_over": broke,
               "grilli": {"S0": round(float(S0), 3), "Hb_over_hb": round(float(0.841 * np.exp(6.421 * S0)), 2)},
               "video": os.path.basename(video), "video_fps": fps, "play_frames": len(play),
               "log_every_2s": log[::max(1, int(round(2.0 / (run.every * run.dt))))] + log[-1:]}
    summary["breaker_closeups"] = os.path.basename(zoom) if zoom else None
    with open(prepare(paths.receipt("summary.json")), "w") as f:
        json.dump(summary, f, indent=1)
    if record:
        res = os.path.join(HERE, "..", "results", "3d", "solitary_beach")
        os.makedirs(res, exist_ok=True)
        with open(os.path.join(res, f"{run.tag}.json"), "w") as f:
            json.dump(dict(summary, recorded_from=f"out/runs/{run.tag}", script="experiments/solitary_beach.py"),
                      f, indent=1)
        if zoom:
            Image.open(zoom).save(os.path.join(res, f"{run.tag}_breaker.png"))
        print("recorded in", os.path.normpath(res), flush=True)
    print(json.dumps({k: v for k, v in summary.items() if k != "log_every_2s"}, indent=1), flush=True)
    print("FINALIZED", run.tag, flush=True)


def main():
    cmd = sys.argv[1] if len(sys.argv) > 1 else "smoke"
    run, words = parse(sys.argv[2:])
    if cmd == "smoke":
        run = Run(*SMOKE, dt=run.dt, tag=run.tag if run.tag != "beach" else "beach_smoke", fused=run.fused)
        for c in range(run.chunks):
            run_chunk(run, c)
        finalize(run)
    elif cmd == "stage":
        run_chunk(run, int(words[0]))
    elif cmd == "finalize":
        finalize(run, record="--record" in sys.argv)
    elif cmd == "remeasure":      # the wave's log again from the dumps, when the rules change
        sc = scene(run)
        paths = RunPaths(run.tag)
        for c in range(run.chunks):
            with open(paths.receipt(f"c{c}.json")) as f:
                entry = json.load(f)
            d = paths.dumps(c)
            entry["wave"] = [measure(np.load(os.path.join(d, f"frame.{k:04d}.npy")),
                                     np.load(os.path.join(d, f"ph.{k:04d}.npy")), run, sc,
                                     (c * run.substeps + k + 1) * run.dt)
                             for k in range(0, run.substeps, run.every)]
            with open(paths.receipt(f"c{c}.json"), "w") as f:
                json.dump(entry, f, indent=1)
            print(f"chunk {c}: {len(entry['wave'])} dumps measured again", flush=True)
        finalize(run, record="--record" in sys.argv)
    elif cmd == "night":
        for c in range(run.chunks):
            run_chunk(run, c)
        finalize(run)
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main()
