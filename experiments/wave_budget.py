"""Read-only capacity planning for a replacement force-driven wave scene.

Reads historical receipts and file sizes. Does not simulate, render, delete,
schedule, or change solver settings. The scaling scenarios are estimates,
not measured speedups or a forecast for the unimplemented force generator.

python experiments/wave_budget.py --output out/wave_budget_20261020.json
"""
import argparse
from collections import Counter
from datetime import datetime
import json
import math
from pathlib import Path
import os
import shutil


ROOT = Path(__file__).resolve().parents[1]
GIB = 1024 ** 3


def inventory():
    groups = Counter()
    top = Counter()
    files = 0
    failures = 0
    for base, _, names in os.walk(ROOT / "out"):
        parts = Path(base).relative_to(ROOT / "out").parts
        category = parts[0] if parts else "(root files)"
        for name in names:
            try:
                size = (Path(base) / name).stat().st_size
            except OSError:
                failures += 1
                continue
            top[category] += size
            files += 1
            if len(parts) >= 2 and parts[0] == "runs":
                key = (parts[1], parts[2] if len(parts) > 2 else "(root files)")
                groups[key] += size
    disk = shutil.disk_usage(ROOT)
    return {
        "meaning": "Logical file sizes; not allocated blocks or a duplicate audit",
        "free_gib": round(disk.free / GIB, 3),
        "out_gib": round(sum(top.values()) / GIB, 3),
        "files": files,
        "unreadable_files": failures,
        "out_top": {k: round(v / GIB, 3) for k, v in top.most_common(12)},
        "run_groups": [
            {"run": r, "category": c, "gib": round(n / GIB, 3)}
            for (r, c), n in groups.most_common(30)
        ],
        "wavesea4_cache_gib": round(groups["wavesea4", "cache"] / GIB, 3),
        "wavesea4_dumps_gib": round(groups["wavesea4", "dumps"] / GIB, 3),
    }


def historical_cost():
    rows = [
        json.loads((ROOT / "out/runs/wavesea4/receipts" / f"c{i}.json").read_text())
        for i in range(5)
    ]
    grid = rows[0]["grid"]
    dt = rows[0]["dt"]
    offset = 0
    for row in rows:
        if row["grid"] != grid or row["dt"] != dt or row["offset"] != offset:
            raise ValueError("Historical receipts do not describe one contiguous configuration")
        offset += row["substeps"]
    total = sum(r["wall_s"] for r in rows)
    return {
        "source": "out/runs/wavesea4/receipts/c0.json through c4.json",
        "thread_evidence": "out/night/sea_night2.log: OMP_NUM_THREADS=32",
        "threads": 32,
        "grid": grid,
        "cells": math.prod(grid),
        "dt_s": dt,
        "steps": offset,
        "simulated_s": offset * dt,
        "wall_hours": total / 3600,
        "mean_step_s": total / offset,
        "chunk_mean_step_s": [r["wall_per_substep_s"] for r in rows],
        "new_force_spray_configuration_measured": False,
    }


def candidates(base, duration=30, budget_hours=15):
    # A capacity envelope, not an approved physical scene. Ceil ensures coverage.
    domain = (180, 32, 72)
    dt = 1 / 48
    steps = round(duration / dt)
    result = []
    for dx in (0.3, 1 / 3, 0.35):
        grid = [math.ceil(length / dx - 1e-10) for length in domain]
        ratio = math.prod(grid) / base["cells"]
        scaled = base["wall_hours"] * ratio * steps / base["steps"]
        result.append({
            "nominal_domain_m": domain,
            "actual_domain_m": [n * dx for n in grid],
            "dx_m": dx,
            "grid": grid,
            "cells": math.prod(grid),
            "dt_s": dt,
            "steps": steps,
            "cell_step_scaling_hours": round(scaled, 3),
            "scenarios_including_one_hour_cache_checks": {
                str(factor): round(scaled * factor + 1, 3)
                for factor in (1, 1.25, 1.5)
            },
            "step_time_limit_with_one_hour_reserved_s": (budget_hours - 1) * 3600 / steps,
        })
    return result


def storage_envelope(particles=32_000_000, seconds=30, fps=12):
    # Existing dump schema: xyz, velocity xyz, ph/esc/scl/rad/spr = 44 bytes.
    frames = math.ceil(seconds * fps)
    raw = particles * frames * 44
    # Full state: xyz, velocity xyz and 11 scalar channels = 68 bytes.
    # Conservatively keep both chunk inputs and outputs at 2-second intervals.
    checkpoints = particles * 68 * 2 * math.ceil(seconds / 2)
    # Lightweight display copy only; raw retains every particle. One million
    # bulk samples plus allowance for one million escaped/emitted particles.
    display = 2_000_000 * frames * 28
    overhead = 60 * GIB
    total = raw + checkpoints + display + overhead
    return {
        "assumed_total_particles": particles,
        "assumed_exact_dump_fps": fps,
        "assumed_frames": frames,
        "raw_gib": round(raw / GIB, 3),
        "checkpoint_inputs_and_outputs_gib": round(checkpoints / GIB, 3),
        "display_copy_gib": round(display / GIB, 3),
        "other_allowance_gib": 60,
        "estimated_gib": round(total / GIB, 3),
        "planning_cap_gib": 650,
        "free_disk_reserve_gib": 250,
        "limitations": (
            "Particle count, spray growth and I/O speed unmeasured. "
            "The 650 GiB cap is a planning target, not an implemented guard. "
            "Playback thinning/export has not been implemented. "
            "No compression or deletion savings are assumed."
        ),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    base = historical_cost()
    disk = inventory()
    storage = storage_envelope()
    report = {
        "measured_at": datetime.now().astimezone().isoformat(),
        "status": "planning_only_no_new_simulation_or_schedule",
        "historical": base,
        "candidates": candidates(base),
        "model_limitations": (
            "Linear cell-step scaling screens candidates only. New forcing, "
            "particle occupancy, pressure iterations, timestep stability and "
            "dump cadence differ. Multipliers are sensitivity scenarios, "
            "not measured overhead or confidence intervals."
        ),
        "storage": storage,
        "inventory": disk,
        "remaining_gib_after_planning_cap_and_reserve":
            round(disk["free_gib"] - storage["planning_cap_gib"] -
                  storage["free_disk_reserve_gib"], 3),
        "delete_required_by_this_envelope": disk["free_gib"] < 900,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "historical_hours": round(base["wall_hours"], 3),
        "candidate_hours_with_25_percent_sensitivity_and_cache": [
            {"dx_m": c["dx_m"], "hours": c["scenarios_including_one_hour_cache_checks"]["1.25"]}
            for c in report["candidates"]
        ],
        "free_gib": disk["free_gib"],
        "estimated_new_gib": storage["estimated_gib"],
        "delete_required_by_this_envelope": report["delete_required_by_this_envelope"],
    }, indent=2))


if __name__ == "__main__":
    main()
