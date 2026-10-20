"""Verify three existing shading-study outputs without starting Maya."""
import json
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

from verify_breaker_render_study import compare, digest, scene

ROOT = Path(__file__).resolve().parents[1]


def main():
    root = ROOT / "out/runs/breaker_shading_study_20261019"
    previous = ROOT / "out/runs/breaker_render_study_20261019/D"
    cache = ROOT / "out/runs/wavesea4/cache/breaking_section"
    hashes = {p.name: digest(p) for p in sorted(cache.iterdir()) if p.is_file()}
    images, receipts = {}, {}
    for arm in "ABC":
        directory = root / arm
        receipt = json.loads((directory / "receipt.json").read_text())
        assert receipt["status"] == "passed" and receipt["exit_code"] == 0, arm
        assert receipt["input_sha256"] == hashes, arm
        assert receipt["scene_sha256"] == digest(directory / "breaker.ma"), arm
        assert receipt["image_sha256"] == digest(directory / "breaker.0001.png"), arm
        with Image.open(directory / "breaker.0001.png") as im:
            assert im.size == (960, 540), arm
            images[arm] = np.array(im.convert("RGB"))
        log = (directory / "breaker.1.log").read_text()
        assert "rendering image at 960 x 540, 4 AA samples" in log, arm
        triangles = int(re.search(r"(\d+) tris -- input0_0_0_graph_outputs_0_mesh", log).group(1))
        assert triangles == 292516, (arm, triangles)
        receipt["water_triangles"] = triangles
        receipts[arm] = receipt
    assert len({r["git_sha"] for r in receipts.values()}) == 1
    with Image.open(previous / "breaker.0001.png") as im:
        assert np.array_equal(np.asarray(im.convert("RGB")), images["A"]), "Control pixels differ"
    control = compare(previous / "breaker.ma", root / "A/breaker.ma")
    matte = compare(root / "A/breaker.ma", root / "B/breaker.ma")
    roughness = compare(root / "A/breaker.ma", root / "C/breaker.ma")
    for diff in [control, matte, roughness]:
        assert not any(diff["changed_graph_values"].values()), diff
    assert not control["changed_nodes"] and not control["added_connections"] and not control["removed_connections"]
    assert matte["changed_nodes"] == ["surfaceDiagnosticMatte"]
    assert matte["removed_connections"] == ['connectAttr "water.out" "waterSG.ss";']
    assert matte["added_connections"] == [
        'connectAttr "surfaceDiagnosticMatte.msg" ":defaultShaderList1.s" -na;',
        'connectAttr "surfaceDiagnosticMatte.out" "waterSG.ss";']
    assert roughness["changed_nodes"] == ["water"]
    assert not roughness["added_connections"] and not roughness["removed_connections"]
    original_water = scene(root / "A/breaker.ma")[0]["water"]
    rough_water = scene(root / "C/breaker.ma")[0]["water"]
    attr = r'setAttr "\.specular_roughness" ([^;]+);'
    assert abs(float(re.search(attr, original_water).group(1)) - .03) < 1e-6
    assert abs(float(re.search(attr, rough_water).group(1)) - .12) < 1e-6
    assert re.sub(attr, "ROUGHNESS", original_water) == re.sub(attr, "ROUGHNESS", rough_water)
    roi = images["B"][170:, 250:]
    clipped = float((roi.min(axis=-1) == 255).mean())
    result = {
        "status": "passed", "control_pixels_identical_to_prior_D": True,
        "scene_differences": {"prior_D_to_A": control, "A_to_B": matte, "A_to_C": roughness},
        "matte_foreground_roi": {"x_from": 250, "y_from": 170, "all_rgb_255_fraction": clipped},
        "total_launch_evaluation_render_wall_s": sum(r["wall_s"] for r in receipts.values()),
        "renders": receipts,
        "limits": ["The production mesh vertices were not exported or hashed.",
                   "Matte highlights clip under the unchanged lighting; they cannot settle fine normal detail.",
                   "No new simulation, surface-resolution arm, spray emission or extra image was run."],
    }
    (root / "verified_shading_study.json").write_text(json.dumps(result, indent=2) + "\n")
    # A comparison figure only: source pixels are ungraded and unretouched.
    fig, axes = plt.subplots(1, 3, figsize=(19.2, 4.6), constrained_layout=True)
    titles = ["A | Water control | roughness 0.03", "B | Matte geometry diagnostic", "C | Water | roughness 0.12"]
    for ax, arm, title in zip(axes, "ABC", titles):
        ax.imshow(images[arm], interpolation="nearest")
        ax.set_title(title, fontsize=11)
        ax.axis("off")
    fig.suptitle("24.9 s | same cache, camera and reconstruction settings | foam / mist hidden", fontsize=13)
    fig.text(.5, .01, "A reproduces the previous D pixel-for-pixel. B has clipped highlights. C softens some reflections but retains the glass-like appearance.",
             ha="center", fontsize=10)
    fig.savefig(root / "comparison.png", dpi=150)
    plt.close(fig)
    print(f"PASS: three receipts, unchanged inputs, intended scene differences, 292,516 triangles each; total {result['total_launch_evaluation_render_wall_s']:.3f} s")
    print(f"Matte foreground ROI: {clipped:.2%} of pixels clip in all RGB channels")


if __name__ == "__main__":
    main()
