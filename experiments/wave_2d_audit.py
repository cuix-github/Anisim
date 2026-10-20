"""Inspect actual cached particles and plot diagnostic snapshots; no rendering."""
import argparse
import json
from pathlib import Path

import numpy as np
from wave_features import atomic_json


def audit(root):
    cache = root / "cache/full"
    plan = json.loads((root / "inputs/plan.json").read_text())
    length, height = plan.get("domain_m", [172, 24])[:2]
    indices = {layer: json.loads((cache / layer / "cache.json").read_text())
               for layer in ("water", "droplets", "spray")}
    records = []
    for frame, time_s in enumerate(indices["water"]["times_s"], start=1):
        record = dict(frame=frame, time_s=time_s)
        for layer in indices:
            p = np.load(cache / layer / f"liquid.{frame:04d}.npy")
            if not np.isfinite(p).all():
                raise ValueError(f"Nonfinite {layer} at frame {frame}")
            record[layer] = len(p)
            if layer != "water":
                # Ten metres beyond the wave zone excludes direct forcing.
                record[layer + "_beyond_zone"] = int(np.count_nonzero(
                    (p[:, 0] > 34) & (p[:, 0] < length - 10) & (p[:, 1] < height - 2)))
                record[layer + "_near_top_or_far_wall"] = int(np.count_nonzero(
                    (p[:, 0] >= length - 2) | (p[:, 1] >= height - 2)))
            else:
                region = p[(p[:, 0] > 34) & (p[:, 0] < 150)]
                record["highest_non_escaped_point_above_still_m"] = float(region[:, 1].max() - 12)
        records.append(record)
    eligible = [r for r in records if r["spray_beyond_zone"] >= 50 and r["droplets_beyond_zone"] >= 200]
    best = max(eligible or records,
               key=lambda r: min(r["spray_beyond_zone"], r["droplets_beyond_zone"]))
    summary = dict(passed=bool(eligible),
                   best_simultaneous_frame=best,
                   frames_with_both_populations=len(eligible),
                   max_near_boundary={key: max(r[key + "_near_top_or_far_wall"] for r in records)
                                      for key in ("droplets", "spray")},
                   peaks={key: max(records, key=lambda r: r[key]) for key in
                          ("spray_beyond_zone", "droplets_beyond_zone",
                           "highest_non_escaped_point_above_still_m")},
                   plan=dict(wave_height_scale=plan["wave_height_scale"],
                             emit_rate=plan["ports"]["emit_rate"],
                             emit_grad_threshold=plan["ports"]["emit_grad_threshold"]))
    atomic_json(root / "receipts/particle_audit.json", summary)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.style.use("dark_background")
    picks = sorted(set((round(24 * 24), best["frame"],
                       max(records, key=lambda r: r["spray_beyond_zone"])["frame"],
                       max(records, key=lambda r: r["droplets_beyond_zone"])["frame"])))
    fig, axes = plt.subplots(len(picks), 1, figsize=(16, 4 * len(picks)), squeeze=False)
    for ax, frame in zip(axes[:, 0], picks):
        for layer, color, size in (("water", "#298bd0", .5), ("droplets", "#a5ffff", 5),
                                    ("spray", "#ff9d28", 10)):
            p = np.load(cache / layer / f"liquid.{frame:04d}.npy")
            ax.scatter(p[:, 0], p[:, 1], s=size, c=color, label=layer, linewidths=0)
        ax.plot([36, length], [0, (length - 36) * .1], c="#d6be92", lw=2)
        ax.axhline(12, c="grey", ls="--", lw=.5)
        ax.set(xlim=(34, min(190, length)), ylim=(0, min(32, height)), ylabel="y (metres)")
        rec = records[frame - 1]
        ax.set_title(f"{rec['time_s']:.2f} s | beyond wave zone: "
                     f"{rec['droplets_beyond_zone']} escaped / {rec['spray_beyond_zone']} emitted")
        ax.legend(loc="upper left", markerscale=2)
        ax.set_aspect("equal")
    fig.suptitle("Actual solver particles; marker sizes enlarged for diagnosis, not droplet diameters")
    axes[-1, 0].set_xlabel("x (metres)")
    fig.tight_layout(rect=(0, 0, 1, .95))
    fig.savefig(root / "particle_diagnostic.png", dpi=120, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps(summary, indent=2))
    if not summary["passed"]:
        raise ValueError("The diagnostic does not yet show both populations beyond the wave zone")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run", type=Path)
    audit(parser.parse_args().run.resolve())
