"""Milestone 6.5a: the long three-way run in Bifrost. The uniform, naive and
cleanup pressure solves (step_3d's pressure_solver 0, 1 and 2) run the
stage-3a scenes at nz=1 through the fused loop graph, long enough to see
whether naive's air stirring levels off, and whether it shows.

Scenes, the 2D record's (reference2d/validate_adaptive_projection.py):
  tank  64x64, half full, at rest. Frames of 0.05 s, two substeps each: the
        reference's own step at rest (dt_frame / 2).
  dam   160x80, a 40x56 column. Frames of 0.06 s, four substeps of 0.015 s:
        CFL about 0.5 at the splash's ~30 cells per second.

Every arm starts from the same seeded state (pfflip2d.Sim.seed, seed 1) and
runs in chunks of CHUNK frames, each resumed from the previous chunk's full
state (bit-exact, gates/validate_resume.py) in its own out/loop3_<tag>
directory, so a chunk killed mid-way reruns alone. Per frame the particle
positions, phases and velocities are dumped; per substep the solver's
diagnostics. The report then makes, per scene:
  metrics  per frame and arm: rms and max speed of the air and of the
           liquid, and pressure iterations per substep. The air rms slope
           over the last 200 frames says whether it is still growing.
  video    the three arms side by side: liquid in blue, air coloured by its
           speed (still air nearly black, moving air bright), and over the
           two adaptive arms the tree's level boundaries, rebuilt from each
           arm's own particles with the reference's band_target and
           build_cells - the tree the next step's solve uses.
  chart    air and liquid rms and iterations against time, three arms.

From the repo root (run is the heavy part; meant for an idle machine):
  .venv\\Scripts\\python.exe experiments\\adaptive_long_run.py smoke
  .venv\\Scripts\\python.exe experiments\\adaptive_long_run.py run --scene tank --frames 1000
  .venv\\Scripts\\python.exe experiments\\adaptive_long_run.py report --scene tank --frames 1000
"""
import argparse
import json
import os
import subprocess
import sys
import time
from types import SimpleNamespace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
sys.path[:0] = [os.path.join(ROOT, "graphs"), os.path.join(ROOT, "reference2d")]
from loop3d import load_final_state, run_bifrost3  # noqa: E402
from pfflip2d import Sim  # noqa: E402
from adaptive_solve import build_cells  # noqa: E402
from adaptive_projection import band_target  # noqa: E402

ARMS = (("uniform", 0), ("naive", 1), ("cleanup", 2))   # name, pressure_solver
SCENES = {
    "tank": {"nx": 64, "ny": 64, "mask": lambda x, y: y < 32, "frame": 0.05, "sub": 2,
             "scale": 6},
    "dam": {"nx": 160, "ny": 80, "mask": lambda x, y: (x < 40) & (y < 56), "frame": 0.06,
            "sub": 4, "scale": 4},
}
CHUNK = 250                       # frames per bifcmd run
DUMPS = ("ph", "vel")             # beside frame.####; escape and scale are off
OUT = os.path.join(ROOT, "out", "adaptive_long")
RESULTS = os.path.join(ROOT, "results", "adaptive")


def seeded(scene):
    sc = SCENES[scene]
    sim = Sim(sc["nx"], sc["ny"])
    sim.seed(sc["mask"])
    pos3 = np.column_stack([sim.pos, np.full(len(sim.pos), 0.5)]).astype(np.float32)
    return pos3, sim.typ.astype(np.int32)


def tag(scene, arm, chunk):
    return f"adaptive_long_{scene}_{arm}_c{chunk:02d}"


def chunk_dir(scene, arm, chunk):
    return os.path.join(ROOT, "out", "loop3_" + tag(scene, arm, chunk))


def chunks(frames):
    """(chunk index, frames in it) covering `frames`."""
    return [(c, min(CHUNK, frames - c * CHUNK)) for c in range((frames + CHUNK - 1) // CHUNK)]


def complete(d, substeps, every):
    last = ((substeps - 1) // every) * every
    return (os.path.exists(os.path.join(d, "final.npy"))
            and os.path.exists(os.path.join(d, f"vel.{last:04d}.npy"))
            and os.path.exists(os.path.join(d, f"diag.{substeps - 1:04d}.npy")))


def run(scene, frames):
    sc = SCENES[scene]
    dt = sc["frame"] / sc["sub"]
    pos3, typ = seeded(scene)
    os.makedirs(OUT, exist_ok=True)
    progress = os.path.join(OUT, f"{scene}_progress.jsonl")
    print(f"{scene}: {sc['nx']}x{sc['ny']}, {len(pos3)} particles, {frames} frames of "
          f"{sc['frame']} s, dt {dt}", flush=True)
    for arm, solver in ARMS:
        state = None
        for c, nf in chunks(frames):
            d = chunk_dir(scene, arm, c)
            sub = nf * sc["sub"]
            if complete(d, sub, sc["sub"]):
                state = load_final_state(d)
                print(f"  {arm} chunk {c}: already complete", flush=True)
                continue
            t0 = time.time()
            d, _, diags = run_bifrost3(
                tag(scene, arm, c), sc["nx"], sc["ny"], 1, None, sub, dt,
                extra_ports={"pressure_solver": solver},
                seed_fn=lambda *_: (pos3, typ), fused=True, dump_every=sc["sub"],
                init_arrays=state, final_state=True, step_offset=c * CHUNK * sc["sub"],
                dumps=DUMPS)
            wall = time.time() - t0
            state = load_final_state(d)
            row = {"scene": scene, "arm": arm, "chunk": c, "frames": nf, "substeps": sub,
                   "wall_s_incl_startup": round(wall, 2),
                   "mean_iterations": float(diags[:, 0].mean()),
                   "max_divergence_after": float(diags[:, 2].max())}
            with open(progress, "a", newline="\n") as f:
                f.write(json.dumps(row) + "\n")
            print(f"  {arm} chunk {c}: {nf} frames in {wall:.1f} s (startup included), "
                  f"{row['mean_iterations']:.1f} iterations per substep", flush=True)


def frame_files(scene, frames):
    """Per global frame: (chunk dir, chunk-local dump index, time in seconds).
    A dump at local substep k holds the state after substep k + 1."""
    sc = SCENES[scene]
    dt = sc["frame"] / sc["sub"]
    out = []
    for c, nf in chunks(frames):
        for f in range(nf):
            k = f * sc["sub"]
            g = c * CHUNK * sc["sub"] + k
            out.append((c, k, (g + 1) * dt))
    return out


def level_map(nx, ny, pos, liquid):
    """The adaptive solve's leaf levels per fine cell, rebuilt as the node does."""
    view = SimpleNamespace(nx=nx, ny=ny, pos=pos[:, :2].astype(float),
                           typ=liquid.astype(np.int32))
    cells = build_cells(nx, ny, 2, band_target(view, 5, 2))
    lev = np.zeros((nx, ny), int)
    for (l, i, j) in cells:
        lev[i << l:(i + 1) << l, j << l:(j + 1) << l] = l
    return lev


BG = np.array([10, 12, 18], np.uint8)
LIQ = np.array([45, 125, 255], np.uint8)
AIR_LO = np.array([34, 34, 44], float)
AIR_HI = np.array([255, 236, 150], float)
SEAM = np.array([225, 60, 215], np.uint8)
AIR_CAP = 2.0                     # cells per second at full air brightness


def draw_panel(nx, ny, s, pos, liquid, speed, lev):
    W, H = nx * s, ny * s
    img = np.empty((H, W, 3), np.uint8)
    img[:] = BG
    px = np.clip((pos[:, 0] * s).astype(int), 0, W - 2)
    py = np.clip(H - 1 - (pos[:, 1] * s).astype(int), 1, H - 1)
    air = ~liquid
    t = np.clip(speed[air] / AIR_CAP, 0.0, 1.0)[:, None]
    col = (AIR_LO + t * (AIR_HI - AIR_LO)).astype(np.uint8)
    for dx, dy in ((0, 0), (1, 0), (0, -1), (1, -1)):
        img[py[air] + dy, px[air] + dx] = col
    for dx, dy in ((0, 0), (1, 0), (0, -1), (1, -1)):
        img[py[liquid] + dy, px[liquid] + dx] = LIQ
    if lev is not None:
        # level boundaries: vertical where lev changes along x, horizontal along y
        vi, vj = np.nonzero(lev[1:, :] != lev[:-1, :])
        for i, j in zip(vi + 1, vj):
            img[H - (j + 1) * s:H - j * s, i * s] = SEAM
        hi, hj = np.nonzero(lev[:, 1:] != lev[:, :-1])
        for i, j in zip(hi, hj + 1):
            img[H - j * s, i * s:(i + 1) * s] = SEAM
    return img


def report(scene, frames):
    from PIL import Image, ImageDraw, ImageFont
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        font = ImageFont.truetype("arial.ttf", 16)
        bold = ImageFont.truetype("arialbd.ttf", 18)
    except OSError:
        font = bold = ImageFont.load_default()
    sc = SCENES[scene]
    nx, ny, s = sc["nx"], sc["ny"], sc["scale"]
    files = frame_files(scene, frames)
    iters = {}
    for arm, _ in ARMS:
        iters[arm] = np.concatenate([
            np.stack([np.load(os.path.join(chunk_dir(scene, arm, c), f"diag.{k:04d}.npy"))
                      for k in range(nf * sc["sub"])])[:, 0]
            for c, nf in chunks(frames)])
    gap, head, top = 6, 52, 34
    W = 3 * nx * s + 2 * gap
    H = top + head + ny * s
    W += W % 2
    H += H % 2
    os.makedirs(OUT, exist_ok=True)
    video = os.path.join(OUT, f"{scene}_three_way.mp4")
    ff = subprocess.Popen(["ffmpeg", "-y", "-loglevel", "error", "-f", "rawvideo",
                           "-pix_fmt", "rgb24", "-s", f"{W}x{H}", "-r", f"{1 / sc['frame']:.4f}",
                           "-i", "-", "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20",
                           video], stdin=subprocess.PIPE)
    rows = {arm: [] for arm, _ in ARMS}
    for n, (c, k, t) in enumerate(files):
        canvas = Image.new("RGB", (W, H), tuple(int(x) for x in BG))
        dr = ImageDraw.Draw(canvas)
        dr.text((10, 8), f"{scene}  t = {t:6.2f} s  (frame {n + 1}/{len(files)})", fill=(235, 235, 235),
                font=bold)
        for a, (arm, solver) in enumerate(ARMS):
            d = chunk_dir(scene, arm, c)
            pos = np.load(os.path.join(d, f"frame.{k:04d}.npy"))
            liquid = np.load(os.path.join(d, f"ph.{k:04d}.npy")) > 0.5
            vel = np.load(os.path.join(d, f"vel.{k:04d}.npy"))
            speed = np.linalg.norm(vel[:, :2], axis=1)
            g = c * CHUNK * sc["sub"] + k          # the substep this dump ends
            it = float(iters[arm][max(0, g - sc["sub"] + 1):g + 1].mean())
            r = {"t": round(t, 4),
                 "rms_air": float(np.sqrt((speed[~liquid] ** 2).mean())),
                 "rms_liq": float(np.sqrt((speed[liquid] ** 2).mean())),
                 "max_air": float(speed[~liquid].max()), "max_liq": float(speed[liquid].max()),
                 "iterations": it}
            rows[arm].append(r)
            lev = level_map(nx, ny, pos, liquid) if solver else None
            x0 = a * (nx * s + gap)
            canvas.paste(Image.fromarray(draw_panel(nx, ny, s, pos, liquid, speed, lev)),
                         (x0, top + head))
            dr.text((x0 + 6, top + 2), f"{arm} (pressure_solver {solver})", fill=(235, 235, 235),
                    font=bold)
            dr.text((x0 + 6, top + 26), f"air rms {r['rms_air']:.3f}  liquid rms "
                    f"{r['rms_liq']:.3f}  iters {it:.0f}", fill=(190, 190, 200), font=font)
        ff.stdin.write(np.asarray(canvas).tobytes())
        if (n + 1) % 100 == 0:
            print(f"  rendered {n + 1}/{len(files)}", flush=True)
    ff.stdin.close()
    if ff.wait() != 0:
        raise SystemExit("ffmpeg failed")

    fig, ax = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    for arm, _ in ARMS:
        tt = [r["t"] for r in rows[arm]]
        ax[0].plot(tt, [r["rms_air"] for r in rows[arm]], label=arm)
        ax[1].plot(tt, [r["rms_liq"] for r in rows[arm]], label=arm)
        ax[2].plot(tt, [r["iterations"] for r in rows[arm]], label=arm)
    for a2, lab in zip(ax, ("air rms speed (cells/s)", "liquid rms speed (cells/s)",
                            "pressure iterations per substep")):
        a2.set_ylabel(lab)
        a2.grid(alpha=0.3)
    ax[0].legend()
    ax[2].set_xlabel("simulated time (s)")
    fig.suptitle(f"{scene}: uniform vs naive vs cleanup, nz=1 in Bifrost")
    fig.tight_layout()
    chart = os.path.join(OUT, f"{scene}_chart.png")
    fig.savefig(chart, dpi=110)

    def window(arm, key, f):          # mean over the 10 frames ending at frame f
        return float(np.mean([r[key] for r in rows[arm][max(0, f - 10):f]]))
    marks = [f for f in (100, 200, 500, 1000, len(files)) if f <= len(files)]
    summary = {"scene": scene, "frames": len(files), "frame_s": sc["frame"],
               "substeps_per_frame": sc["sub"], "arms": {}}
    for arm, solver in ARMS:
        last = rows[arm][-200:]
        slope = (np.polyfit([r["t"] for r in last], [r["rms_air"] for r in last], 1)[0]
                 * 50 * sc["frame"]) if len(last) > 10 else None
        summary["arms"][arm] = {
            "pressure_solver": solver,
            "mean_iterations_per_substep": round(float(iters[arm].mean()), 2),
            "air_rms_at": {str(f): round(window(arm, "rms_air", f), 4) for f in marks},
            "liquid_rms_at": {str(f): round(window(arm, "rms_liq", f), 4) for f in marks},
            "air_rms_slope_per_50_frames_last_200": None if slope is None else round(float(slope), 4),
            "max_liquid_speed": round(max(r["max_liq"] for r in rows[arm]), 3),
            "max_air_speed": round(max(r["max_air"] for r in rows[arm]), 3)}
    with open(os.path.join(OUT, f"{scene}_rows.json"), "w", newline="\n") as f:
        json.dump(rows, f)
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, f"long_run_{scene}.json"), "w", newline="\n") as f:
        json.dump(summary, f, indent=1)
    print(json.dumps(summary, indent=1))
    print(f"video {video}\nchart {chart}")


def smoke():
    """Light: 8 substeps per arm on the tank. The ports reach step_3d through
    the loop graph (the three arms' iteration counts differ), velocities are
    dumped, and one tree rebuild is timed for each scene."""
    sc = SCENES["tank"]
    pos3, typ = seeded("tank")
    its = {}
    for arm, solver in ARMS:
        d, _, diags = run_bifrost3(f"adaptive_long_smoke_{arm}", sc["nx"], sc["ny"], 1, None, 8,
                                   sc["frame"] / sc["sub"], extra_ports={"pressure_solver": solver},
                                   seed_fn=lambda *_: (pos3, typ), fused=True,
                                   dump_every=sc["sub"], dumps=DUMPS)
        vel = np.load(os.path.join(d, "vel.0006.npy"))
        assert vel.shape == (len(pos3), 3), vel.shape
        its[arm] = diags[:, 0].tolist()
        print(f"  {arm}: iterations per substep {its[arm]}", flush=True)
    ok = its["uniform"] != its["naive"] and its["naive"] != its["cleanup"]
    for scene in SCENES:
        p, t = seeded(scene)
        t0 = time.time()
        level_map(SCENES[scene]["nx"], SCENES[scene]["ny"], p, t == 1)
        print(f"  tree rebuild for a {scene} frame: {time.time() - t0:.2f} s", flush=True)
    print("SMOKE OK" if ok else "SMOKE FAIL: arms did not differ")
    sys.exit(0 if ok else 1)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=("run", "report", "smoke"))
    ap.add_argument("--scene", choices=tuple(SCENES))
    ap.add_argument("--frames", type=int, default=1000)
    a = ap.parse_args()
    if a.what == "smoke":
        smoke()
    elif a.scene is None:
        ap.error("--scene is required")
    elif a.what == "run":
        run(a.scene, a.frames)
    else:
        report(a.scene, a.frames)


if __name__ == "__main__":
    main()
