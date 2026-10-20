"""Small nz=1 wave basin, executed by our real fused Bifrost simulation graph.

python experiments/wave_2d.py --tag wave2d_spray_20261020
Default: 344x48x1, 0.5 m cells, 30 seconds, exact samples at 24 fps.
All wave headings are made unidirectional; phases agree at x=100 m with the
3D sea's z=40 m slice. This is a diagnostic analogue, not a slice of the 3D solution.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "graphs"))
import wave_basin as basin
from loop3d import DUMPS, load_final_state, run_bifrost3
from runs import RunPaths
from wave_features import atomic_json, exact_frame, spray_ports
import numpy as np


def flat_waves(run, z_m=40.0, height_scale=1.0):
    rows = np.asarray(basin.sea_rows(), dtype=float)
    for row in rows:
        height, period, heading, phase = row
        theta = np.radians(heading)
        row[3] = phase + basin._k(period) * (100.0 * (np.cos(theta) - 1) + z_m * np.sin(theta))
        row[0] *= run.cpm * height_scale
        row[2] = 0.0
    return rows.astype(np.float32).ravel()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--seconds", type=float, default=30.0)
    ap.add_argument("--cpm", type=float, default=2.0)
    ap.add_argument("--chunk-seconds", type=float, default=5.0)
    ap.add_argument("--wave-height-scale", type=float, default=1.0)
    ap.add_argument("--emit-grad-threshold", type=float, default=10.0)
    ap.add_argument("--emit-rate", type=float, default=0.05)
    ap.add_argument("--until", type=float, help="Pause at a chunk boundary for inspection")
    ap.add_argument("--resume", action="store_true")
    args = ap.parse_args()
    if not args.tag.replace("_", "").replace("-", "").isalnum():
        ap.error("tag must be a simple directory name")
    if not all(np.isfinite(v) and v > 0 for v in
               (args.seconds, args.cpm, args.chunk_seconds, args.wave_height_scale,
                args.emit_grad_threshold, args.emit_rate)):
        ap.error("Durations, resolution, wave scale and emitter settings must be positive")
    run = basin.Run(args.cpm)
    run.grid = (round(basin.LX * args.cpm), round(basin.LY * args.cpm), 1)
    run.dt, run.tag = 1 / 48, args.tag
    paths = RunPaths(run.tag)
    root = Path(paths.root)
    sc = basin.scene(run)
    ports = dict(basin.ports(run, sc), **spray_ports(run.cpm))
    ports.update(emit_grad_threshold=args.emit_grad_threshold, emit_rate=args.emit_rate)
    wave_array = flat_waves(run, height_scale=args.wave_height_scale)
    plan = dict(tag=run.tag, grid=run.grid, dt=run.dt, seconds=args.seconds,
                ports=ports, cells_per_metre=run.cpm, waves=wave_array.tolist(),
                wave_height_scale=args.wave_height_scale,
                interpretation="Unidirectional nz=1 diagnostic, not the 3D solution",
                source_graph="graphs/sim_3d_fused.json",
                adaptive_pressure="Uniform: adaptive pressure does not support the solid beach")
    total = round(args.seconds / run.dt)
    chunk_steps = round(args.chunk_seconds / run.dt)
    stop = round((args.until or args.seconds) / run.dt)
    if total % 2 or chunk_steps % 2 or stop > total or stop <= 0 or (
            stop != total and stop % chunk_steps):
        ap.error("Use whole cache-frame durations and pause only at chunk boundaries")
    counts, times, records, offset_start, prior_wall = [], [], [], 0, 0.0
    if args.resume:
        saved = json.loads((root / "inputs/plan.json").read_text())
        if saved != json.loads(json.dumps(plan)):
            raise ValueError("Resume configuration differs from the original plan")
        progress = json.loads((root / "receipts/progress.json").read_text())
        counts, times, records = (progress[key] for key in ("counts", "times", "chunks"))
        offset_start = records[-1]["offset"] + records[-1]["steps"]
        prior_wall = progress["wall_s"]
        if stop <= offset_start:
            raise ValueError("No new time requested")
        init = load_final_state(paths.dumps(records[-1]["chunk"]))
    else:
        root.mkdir(parents=True, exist_ok=False)
        atomic_json(root / "inputs/plan.json", plan)
        init = basin.seed(run, sc)
        init["positions"][:, 2] = 0.5
    start = time.monotonic()
    for chunk, offset in enumerate(range(offset_start, stop, chunk_steps), start=len(records)):
        nstep = min(chunk_steps, stop - offset)
        d = paths.dumps(chunk)
        t0 = time.monotonic()
        run_bifrost3(run.tag, *run.grid, None, nstep, run.dt,
                     extra_ports=ports, init_arrays=init, release_init=True,
                     final_state=True, step_offset=offset, dump_every=2,
                     solid=(basin.beach_phi(run, sc), 2), dumps=DUMPS + ("vel",),
                     fused=True, out_dir=d, waves=wave_array,
                     log_prefix=str(root / "logs" / f"c{chunk}"))
        elapsed = time.monotonic() - t0
        diags = np.stack([np.load(Path(d) / f"diag.{k:04d}.npy") for k in range(nstep)])
        if not np.isfinite(diags).all() or (diags[:, 0] < 0).any() or (diags[:, 1] > 1.05e-6).any():
            raise ValueError("Invalid solver diagnostics")
        for k in range(0, nstep, 2):
            for prefix in ("frame", "vel", "rad"):
                if not np.isfinite(np.load(Path(d) / f"{prefix}.{k:04d}.npy")).all():
                    raise ValueError(f"Nonfinite {prefix} at {offset+k}")
            counts.append(exact_frame(d, k, root / "cache" / "full", len(counts) + 1,
                                      cells_per_metre=run.cpm, layers=True))
            times.append((offset + k + 1) * run.dt)
        rec = dict(chunk=chunk, offset=offset, steps=nstep, wall_s=elapsed,
                   end_s=(offset + nstep) * run.dt,
                   max_spray=max(c["spray"] for c in counts),
                   max_droplets=max(c["droplets"] for c in counts),
                   iterations_max=float(diags[:, 0].max()),
                   residual_max=float(diags[:, 1].max()))
        records.append(rec)
        atomic_json(root / "receipts" / f"c{chunk}.json", rec)
        atomic_json(root / "status.json", dict(state="running", **rec))
        atomic_json(root / "receipts/progress.json",
                    dict(counts=counts, times=times, chunks=records,
                         wall_s=prior_wall + time.monotonic() - start))
        print(json.dumps(rec), flush=True)
        if offset + nstep < stop:
            init = load_final_state(d)
    for layer in ("full", "water", "droplets", "spray", "air"):
        d = root / "cache" / "full"
        if layer != "full":
            d /= layer
        atomic_json(d / "cache.json", dict(frames=len(counts), fps=24, times_s=times,
                    counts=[c[layer] for c in counts], units="metres",
                    grid=[basin.LX, basin.LY, 1 / run.cpm], source=run.tag,
                    interpolation="none; exact solver snapshots every 2 substeps"))
    atomic_json(root / "receipts" / "summary.json",
                dict(wall_s=prior_wall + time.monotonic() - start, chunks=records,
                     max_counts={key: max(c[key] for c in counts) for key in counts[0]},
                     frames=len(counts), cache_times_s=[times[0], times[-1]]))
    state = "complete" if stop == total else "inspection_pause"
    atomic_json(root / "status.json", dict(state=state, seconds=stop * run.dt,
                frames=len(counts), wall_s=prior_wall + time.monotonic() - start))
    print(state.upper(), root, flush=True)


if __name__ == "__main__":
    main()
