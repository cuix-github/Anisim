"""Particle-only inspection figures and simple single-wave calibration receipts."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "experiments"))
import wave_basin


def audit(root):
    plan = json.loads((root / "inputs/plan.json").read_text())
    progress = json.loads((root / "receipts/progress.json").read_text())
    records = progress["frames"]
    ts = np.array([r["time_s"] for r in records])
    esc = np.array([r["escaped_beyond_source"] for r in records])
    spr = np.array([r["spray_beyond_source"] for r in records])
    picks = sorted(set([max(1, int(len(ts) * x)) for x in (.4, .65, .85)] +
                       [int(np.argmax(esc + spr)) + 1]))
    fig, axs = plt.subplots(3, 2, figsize=(14, 10), facecolor="#181c22")
    for ax in axs.flat:
        ax.set_facecolor("#11161e")
        ax.tick_params(colors="white")
        for spine in ax.spines.values():
            spine.set_color("#777")
    for ax, frame in zip(axs.flat, picks):
        for layer, color, size in (("water", "#2279c5", .8), ("droplets", "#65edff", 5), ("spray", "#ffb847", 8)):
            p = np.load(root / f"cache/full/{layer}/liquid.{frame:04d}.npy")
            if plan["grid"][2] > 1 and layer == "water":
                mid = plan["domain_m"][2] / 2
                p = p[np.abs(p[:, 2] - mid) < 6]
            ax.scatter(p[:, 0], p[:, 1], s=size, c=color, linewidths=0, label=layer)
        toe, slope = plan["beach_toe_m"], plan["beach_slope"]
        x = np.array([0, toe, plan["domain_m"][0]])
        ax.plot(x, np.maximum(x - toe, 0) * slope, c="#a58a6d")
        ax.axvspan(14, 62, color="#aa4444", alpha=.1)
        ax.set_xlim(0, plan["domain_m"][0])
        ax.set_ylim(0, 28)
        ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{ts[frame-1]:.2f} s | escaped {esc[frame-1]} | spray {spr[frame-1]}", color="white")
        ax.set_xlabel("x (m)", color="white")
    axs[2, 0].plot(ts, esc, color="#65edff", label="escaped liquid beyond source")
    axs[2, 0].plot(ts, spr, color="#ffb847", label="extra spray beyond source")
    axs[2, 0].legend(fontsize=8)
    axs[2, 0].set_xlabel("Simulation time (s)", color="white")
    axs[2, 1].plot(ts, [r["liquid_max_y_m"] for r in records], c="#7fc4ff")
    axs[2, 1].axhline(12, color="#777")
    axs[2, 1].set_title("Highest liquid particle beyond source (not wave height)", color="white")
    axs[2, 1].set_xlabel("Simulation time (s)", color="white")
    fig.suptitle("Particle diagnostic | blue: bulk, cyan: escaped water, orange: emitted spray", color="white")
    fig.tight_layout()
    fig.savefig(root / "particle_diagnostic.png", dpi=140, facecolor=fig.get_facecolor())
    plt.close(fig)
    if len(plan["rows"]) == 1:
        probes = {x: [] for x in (75, 95)}
        for chunk in progress["chunks"]:
            d = root / "dumps" / f"c{chunk['chunk']}"
            for f in sorted(d.glob("frame.*.npy")):
                p = np.load(f)
                ph = np.load(d / f.name.replace("frame.", "ph."))
                for x in probes:
                    q = p[(ph > .5) & (np.abs(p[:, 0] * plan["dx"] - x) < 2.5)]
                    probes[x].append(float(np.percentile(q[:, 1], 97.5)) * plan["dx"])
        fit = {}
        om = 2 * np.pi / plan["rows"][0][1]
        for x, values in probes.items():
            use = ts >= 20
            matrix = np.column_stack((np.cos(om * ts[use]), np.sin(om * ts[use]), np.ones(use.sum())))
            coeff = np.linalg.lstsq(matrix, np.array(values)[use], rcond=None)[0]
            fit[x] = dict(height_m=float(2 * np.hypot(*coeff[:2])), phase=float(np.arctan2(coeff[1], coeff[0])),
                          mean_particle_surface_m=float(coeff[2]))
        delta = (fit[95]["phase"] - fit[75]["phase"]) % (2 * np.pi)
        wavelength = 2 * np.pi * 20 / delta
        theory = 2 * np.pi / wave_basin._k(plan["rows"][0][1], 12)
        result = dict(probes=fit, measured_wavelength_m=wavelength, linear_theory_wavelength_m=theory,
                      wavelength_relative_error=abs(wavelength / theory - 1), fit_interval_s=[20, float(ts[-1])],
                      interpretation="97.5th height percentile in each 5 m probe; harmonic fit; not a surface mesh")
        (root / "receipts/single_wave.json").write_text(json.dumps(result, indent=2) + "\n")
        print(json.dumps(result, indent=2))
    print(root / "particle_diagnostic.png")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("run", type=Path)
    audit(ap.parse_args().run.resolve())
