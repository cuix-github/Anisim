"""Inspect final states at chunk boundaries; does not run the solver."""
import argparse
import json
from pathlib import Path
import numpy as np


def main(root):
    plan = json.loads((root / "inputs/plan.json").read_text())
    progress = json.loads((root / "receipts/progress.json").read_text())
    checkpoints = []
    for chunk in progress["chunks"]:
        d = root / "dumps" / f"c{chunk['chunk']}"
        p = np.load(d / "final.npy") * plan["dx"]
        v = np.load(d / "final_velocities.npy") * plan["dx"]
        ph = np.load(d / "final_phase.npy")
        spray = np.load(d / "final_spray.npy") > .5
        assert all(np.isfinite(a).all() for a in (p, v, ph))
        water = (ph >= .5) & ~spray
        speed = np.linalg.norm(v[water], axis=1)
        rear = water & (p[:, 0] < 10)
        interior = water & (p[:, 0] > 68)
        checkpoints.append(dict(
            time_s=chunk["end_s"], non_emitted_liquid_particles=int(water.sum()),
            liquid_speed_p99_m_s=float(np.percentile(speed, 99)),
            liquid_speed_max_m_s=float(speed.max()),
            liquid_max_displacement_cells_per_step=float(speed.max() * plan["dt"] / plan["dx"]),
            rear_highest_liquid_m=float(p[rear, 1].max()) if rear.any() else None,
            downstream_highest_liquid_m=float(p[interior, 1].max()) if interior.any() else None,
            liquid_within_1m_of_ceiling=int(np.count_nonzero(water & (p[:, 1] > plan["domain_m"][1] - 1)))))
    receipt = dict(
        inspected="Full checkpoint states at 5 s intervals, not every substep",
        limitation="Highest particle is not crest height. Liquid includes escaped droplets. This is not a convergence test.",
        checkpoints=checkpoints,
        sampled_peak_escaped_beyond_source=max(f["escaped_beyond_source"] for f in progress["frames"]),
        sampled_peak_extra_spray_beyond_source=max(f["spray_beyond_source"] for f in progress["frames"]),
        elapsed_including_bake_export_s=progress["elapsed_s"],
        graph_calls_including_startup_io_s=sum(c["wall_s"] for c in progress["chunks"]))
    (root / "receipts/health.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    main(parser.parse_args().run.resolve())
