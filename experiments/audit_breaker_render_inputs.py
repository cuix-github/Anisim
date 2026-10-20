"""Inspect the frozen breaker cache and saved scenes; never start Maya or a solver.

The particle plot is an orthographic data plot, not a reconstructed water surface.
Screen-size estimates use a 35 mm lens and 36 mm horizontal film aperture.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from verify_breaker_render_study import scene

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    cache = ROOT / "out/runs/wavesea4/cache/breaking_section"
    study = ROOT / "out/runs/breaker_render_study_20261019"
    files = [cache / n for n in ["cache.json", "liquid.0001.npy", "radius.0001.npy",
                                "speed.0001.npy", "beach.obj"]]
    hashes = {p.name: digest(p) for p in files}
    metadata = json.loads((cache / "cache.json").read_text())
    assert metadata["times_s"] == [24.9]
    summary = json.loads((ROOT / "out/runs/wavesea4/receipts/summary.json").read_text())
    dx = 1.0 / summary["cells_per_metre"]
    cells = np.load(cache / "liquid.0001.npy", mmap_mode="r")
    radius = np.load(cache / "radius.0001.npy", mmap_mode="r")
    assert cells.shape == (metadata["liquid_particles"][0], 3)
    assert radius.shape == (len(cells),)
    assert np.isfinite(cells).all() and np.isfinite(radius).all()
    # Exact reader() condition in generate_particle_surface.py, including equality.
    surface = radius < 1e-9
    points = cells * dx
    report = {
        "source_revision": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "worktree_status": subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True).strip(),
        "cache_hashes": hashes,
        "time_s": metadata["times_s"][0], "cell_size_m": dx,
        "particles": len(points), "surface_input_particles": int(surface.sum()),
        "excluded_from_surface_by_radius": int((~surface).sum()),
        "bounds_m": [points.min(0).tolist(), points.max(0).tolist()],
        "height_quantiles_m": {}, "saved_renders": {},
        "limits": [
            "No production mesh or level set is evaluated or measured here.",
            "Surface input membership is not a test of survival in the final mesh.",
            "The radius channel is a render classification; it is not each parcel's geometric radius.",
            "Screen size estimates assume horizontal film fit and ignore occlusion.",
            "Cell size alone does not establish equal resolved physics or screen-space detail."
        ],
    }
    for label, mask in [("surface_input", surface), ("excluded", ~surface)]:
        report["height_quantiles_m"][label] = np.quantile(points[mask, 1], [0, .1, .5, .9, 1]).tolist()
    for arm in ["A", "B", "C", "D", "E", "candidate_D"]:
        receipt = json.loads((study / arm / "receipt.json").read_text())
        assert receipt["input_sha256"] == hashes, f"Cache differs from {arm}'s recorded input"
        assert digest(study / arm / "breaker.ma") == receipt["scene_sha256"], f"Scene {arm} changed"
        nodes, graphs, connections = scene(study / arm / "breaker.ma")
        log_path = study / arm / "breaker.1.log"
        log = log_path.read_text()
        knobs = receipt["overrides"]
        eye, target = np.array([float(x) for x in knobs["SURF_CAMERA"].split()]).reshape(2, 3) * dx
        distance = float(np.linalg.norm(target - eye))
        focal_px = int(knobs["SURF_WIDTH"]) * 35 / 36
        geo = {
            "eye_m": eye.tolist(), "target_m": target.tolist(), "target_distance_m": distance,
            "surface_voxel_m": .5 * dx, "particle_sphere_radius_m": .5 * dx,
            "voxel_px_at_target_estimate": .5 * dx * focal_px / distance,
            "sphere_diameter_px_at_target_estimate": dx * focal_px / distance,
            "triangles": int(re.search(r"(\d+) tris -- input0_0_0_graph_outputs_0_mesh", log).group(1)),
            "ray_depths": {}, "scene_sha256": digest(study / arm / "breaker.ma"),
            "log_sha256": digest(log_path),
            "water_roughness": float(re.search(r'setAttr "\.specular_roughness" ([^;]+);', nodes["water"]).group(1)),
            "mcf_iterations": int(graphs["particleSurfaceShape"]["particle_surface.mcf_iterations"].rstrip("L")),
            "mist_hidden": bool(re.search(r'setAttr[^;\n]*"\.v"\s+(?:no|0|false);', nodes["particleMistShape"])),
            "water_assigned_directly": 'connectAttr "water.out" "waterSG.ss";' in connections,
        }
        for name in ["diffuse", "specular", "transmission", "volume indirect"]:
            match = re.search(re.escape(name) + r"\s+samples\s+\d+ / depth\s+(\d+)", log)
            geo["ray_depths"][name] = int(match.group(1))
        assert geo["mist_hidden"] == (arm in ["C", "D", "E", "candidate_D"]), arm
        assert geo["water_assigned_directly"] == (arm != "A"), arm
        report["saved_renders"][arm] = geo

    # Full-depth projection and a thin central slice distinguish overlapping points
    # from an actual local cross-section. Plot every selected particle, no decimation.
    fig, axes = plt.subplots(2, 1, figsize=(13, 8), constrained_layout=True)
    selections = [np.ones(len(points), dtype=bool), (points[:, 2] >= 67.75) & (points[:, 2] < 68.25)]
    titles = ["All cache particles projected from the side (18 m depth overlaps)",
              "Local side slice: z = 67.75 to 68.25 m (0.5 m thick)"]
    report["plot_selections"] = []
    for ax, selected, title in zip(axes, selections, titles):
        for member, color, label, size in [(surface, "#167f99", "Sent to surface reconstruction", .4),
                                          (~surface, "#e06b24", "Excluded from surface; hidden when mist is off", 5)]:
            part = points[selected & member]
            ax.scatter(part[:, 0], part[:, 1], s=size, c=color, linewidths=0, label=label)
        ax.set(title=title, xlabel="x (metres)", ylabel="y (metres)", xlim=(120, 166), ylim=(8, 22))
        ax.set_aspect("equal", adjustable="box")
        ax.grid(alpha=.15)
        ax.legend(loc="lower left", fontsize=8, markerscale=2)
        report["plot_selections"].append({"title": title, "particles": int(selected.sum()),
                                         "excluded": int((selected & ~surface).sum())})
    fig.suptitle("24.9 s frozen cache | particle positions only, no meshing or rendering\n"
                 "Plot marker sizes are for legibility, not physical droplet radii", fontsize=12)
    fig.savefig(args.out / "particle_side_views.png", dpi=150)
    plt.close(fig)
    assert hashes == {p.name: digest(p) for p in files}, "Source cache changed during inspection"
    report["wall_s"] = time.perf_counter() - started
    (args.out / "audit.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: report[k] for k in ["particles", "surface_input_particles", "excluded_from_surface_by_radius", "plot_selections", "wall_s"]}, indent=2))


if __name__ == "__main__":
    main()
