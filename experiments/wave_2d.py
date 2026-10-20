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


def flat_waves(run, z_m=40.0):
    rows = np.asarray(basin.sea_rows(), dtype=float)
    for row in rows:
        height, period, heading, phase = row
        theta = np.radians(heading)
        row[3] = phase + basin._k(period) * (100.0 * (np.cos(theta) - 1) + z_m * np.sin(theta))
        row[0] *= run.cpm
        row[2] = 0.0
    return rows.astype(np.float32).ravel()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--seconds", type=float, default=30.0)
    ap.add_argument("--cpm", type=float, default=2.0)
    ap.add_argument("--chunk-seconds", type=float, default=5.0)
    args = ap.parse_args()
    if not args.tag.replace("_", "").replace("-", "").isalnum():
        ap.error("tag must be a simple directory name")
    run = basin.Run(args.cpm)
    run.grid = (round(basin.LX * args.cpm), round(basin.LY * args.cpm), 1)
    run.dt, run.tag = 1 / 48, args.tag
    paths = RunPaths(run.tag)
    root = Path(paths.root)
    root.mkdir(parents=True, exist_ok=False)
    sc = basin.scene(run)
    ports = dict(basin.ports(run, sc), **spray_ports(run.cpm))
    plan = dict(tag=run.tag, grid=run.grid, dt=run.dt, seconds=args.seconds,
                ports=ports, cells_per_metre=run.cpm, waves=flat_waves(run).tolist(),
                interpretation="Unidirectional nz=1 diagnostic, not the 3D solution",
                source_graph="graphs/sim_3d_fused.json",
                adaptive_pressure="Uniform: adaptive pressure does not support the solid beach")
    atomic_json(root / "inputs" / "plan.json", plan)
    init = basin.seed(run, sc)
    init["positions"][:, 2] = 0.5
    total = round(args.seconds / run.dt)
    chunk_steps = round(args.chunk_seconds / run.dt)
    counts, times, records = [], [], []
    start = time.monotonic()
    for chunk, offset in enumerate(range(0, total, chunk_steps)):
        nstep = min(chunk_steps, total - offset)
        d = paths.dumps(chunk)
        t0 = time.monotonic()
        run_bifrost3(run.tag, *run.grid, None, nstep, run.dt,
                     extra_ports=ports, init_arrays=init, release_init=True,
                     final_state=True, step_offset=offset, dump_every=2,
                     solid=(basin.beach_phi(run, sc), 2), dumps=DUMPS + ("vel",),
                     fused=True, out_dir=d, waves=flat_waves(run),
                     log_prefix=str(root / "logs" / f"c{chunk}"))
        elapsed = time.monotonic() - t0
        diags = np.stack([np.load(Path(d) / f"diag.{k:04d}.npy") for k in range(nstep)])
        if not np.isfinite(diags).all() or (diags[:, 0] < 0).any():
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
                   iterations_max=float(diags[:, 0].max()))
        records.append(rec)
        atomic_json(root / "receipts" / f"c{chunk}.json", rec)
        atomic_json(root / "status.json", dict(state="running", **rec))
        print(json.dumps(rec), flush=True)
        if offset + nstep < total:
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
                dict(wall_s=time.monotonic() - start, chunks=records,
                     max_counts={key: max(c[key] for c in counts) for key in counts[0]},
                     frames=len(counts), cache_times_s=[times[0], times[-1]]))
    atomic_json(root / "status.json", dict(state="complete", seconds=args.seconds,
                                          frames=len(counts), wall_s=time.monotonic() - start))
    print("COMPLETE", root, flush=True)


if __name__ == "__main__":
    main()
