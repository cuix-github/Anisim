"""One explicitly selected still per invocation; never runs the solver.

The 24.9 s wavesea4 study: A control, B foam off, C also mist hidden,
D C with less surface smoothing, E C with a four-times brighter environment.
Outputs, scene files, hashes, full launch logs and settings stay under out/.
Run each arm once and inspect it before continuing. The candidate accepts
the letter of the chosen preview, with optional explicit preset and camera
overrides for the presentation shot (not another single-variable comparison).
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

from PIL import Image

ROOT = Path(__file__).resolve().parents[3]
ARMS = {
    "A": {},
    "B": {"SURF_FOAM_STRENGTH": "0"},
    "C": {"SURF_FOAM_STRENGTH": "0", "SURF_MIST_VISIBLE": "0"},
    "D": {"SURF_FOAM_STRENGTH": "0", "SURF_MIST_VISIBLE": "0",
          "SURF_MCF_ITERS": "10"},
    "E": {"SURF_FOAM_STRENGTH": "0", "SURF_MIST_VISIBLE": "0",
          "SURF_SKY_COLOR": "0.12 0.14 0.18"},
}
CAMERA = "596 76 188 570 58 256"  # cells: eye (149,19,47)m, aim (142.5,14.5,64)m


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("arm", choices=ARMS)
    parser.add_argument("--candidate", action="store_true")
    parser.add_argument("--preset", default="paper_clean")
    parser.add_argument("--camera", help="Candidate only: six values in cells")
    parser.add_argument("--out", type=Path, default=ROOT / "out/runs/breaker_render_study_20261019")
    args = parser.parse_args()
    if not args.candidate and (args.preset != "paper_clean" or args.camera):
        raise SystemExit("Preview arms use their fixed preset and camera")
    if args.camera:
        if len(args.camera.split()) != 6:
            raise SystemExit("Camera requires six numeric values")
        tuple(float(v) for v in args.camera.split())
    if git("status", "--porcelain"):
        raise SystemExit("Dirty checkout: commit the intended instrument before measuring")
    cache = ROOT / "out/runs/wavesea4/cache/breaking_section"
    inputs = {p.name: digest(p) for p in sorted(cache.iterdir()) if p.is_file()}
    index = json.loads((cache / "cache.json").read_text())
    assert index["times_s"] == [24.9] and index["liquid_particles"] == [811717]
    output_root = args.out.resolve()
    label = ("candidate_" if args.candidate else "") + args.arm
    output = output_root / label
    if output.exists():
        raise SystemExit(f"Refusing to reuse an arm directory: {output}")
    if args.arm != "A" or args.candidate:
        control = json.loads((output_root / "A/receipt.json").read_text())
        assert control["status"] == "passed", "Control must pass first"
        assert control["input_sha256"] == inputs, "Cache changed since control"
    if args.candidate:
        chosen = json.loads((output_root / args.arm / "receipt.json").read_text())
        assert chosen["status"] == "passed", "Chosen preview must pass first"
    output.mkdir(parents=True)
    width, height = (1920, 1080) if args.candidate else (960, 540)
    arm_overrides = ARMS[args.arm] if args.preset == "paper_clean" else {}
    overrides = {"SURF_WIDTH": str(width), "SURF_HEIGHT": str(height),
                 "SURF_CAMERA": args.camera or CAMERA, "SURF_METRES_PER_CELL": "0.25",
                 "SURF_THREADS": "0", "SURF_GPU": "0", **arm_overrides}
    child_env = {k: v for k, v in os.environ.items() if not k.startswith("SURF_")}
    child_env.update(overrides)
    receipt = {"arm": label, "git_sha": git("rev-parse", "HEAD"),
               "cache": str(cache.relative_to(ROOT)), "time_s": 24.9,
               "input_sha256": inputs, "preset": args.preset,
               "overrides": overrides, "status": "running",
               "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    receipt_path = output / "receipt.json"
    receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
    command = ["powershell.exe", "-NoProfile", "-ExecutionPolicy", "Bypass", "-File",
               str(ROOT / "maya/render_surface_preview.ps1"), "-Cache", str(cache),
               "-Out", str(output), "-Name", "breaker", "-Frames", "1",
               "-Obj", str(cache / "beach.obj"), "-Preset", args.preset,
               "-MaxMinutes", "12"]
    started = time.perf_counter()
    print(f"START {label}: {width}x{height}, cache 24.9 s, {output}", flush=True)
    with (output / "launch.log").open("w", encoding="utf-8") as log:
        result = subprocess.run(command, cwd=ROOT, env=child_env, stdout=log,
                                stderr=subprocess.STDOUT)
    receipt["wall_s"] = round(time.perf_counter() - started, 3)
    receipt["exit_code"] = result.returncode
    try:
        assert result.returncode == 0, "Launcher failed; inspect launch.log"
        notes = (output / "breaker.notes.txt").read_text()
        assert "ERROR" not in notes and "NO IMAGE" not in notes, notes
        frame = output / "breaker.0001.png"
        with Image.open(frame) as im:
            assert im.size == (width, height), f"Unexpected size: {im.size}"
            im.verify()
        assert inputs == {p.name: digest(p) for p in sorted(cache.iterdir()) if p.is_file()}, "Cache mutated"
        receipt["image_sha256"] = digest(frame)
        receipt["scene_sha256"] = digest(output / "breaker.ma")
        receipt["notes"] = notes
        receipt["status"] = "passed"
    except Exception as exc:
        receipt["status"] = "failed"
        receipt["error"] = str(exc)
        raise
    finally:
        receipt_path.write_text(json.dumps(receipt, indent=2) + "\n")
        print(json.dumps(receipt, indent=2), flush=True)


if __name__ == "__main__":
    main()
