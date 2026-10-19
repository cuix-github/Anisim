"""Three approved stills from the frozen 24.9 s cache, one per invocation.

A reproduces the earlier D control (10 MCF, foam/mist hidden).
B changes only the assigned surface material to neutral matte.
C changes only water roughness from 0.03 to 0.12.
All use the original 960 x 540 preview camera and paper_clean lighting.
The existing runner enforces clean revisions, separate arm directories,
control-first ordering, unchanged cache hashes and complete render receipts.
"""
import sys

import breaker_render_study as runner


def main():
    baseline = {**runner.ARMS["D"], "SURF_WATER_ROUGHNESS": "0.03",
                "SURF_DIAGNOSTIC_MATTE": "0"}
    runner.ARMS = {
        "A": baseline,
        "B": {**baseline, "SURF_DIAGNOSTIC_MATTE": "1"},
        "C": {**baseline, "SURF_WATER_ROUGHNESS": "0.12"},
    }
    runner.__doc__ = __doc__
    if "--candidate" in sys.argv:
        raise SystemExit("This study only authorizes the three preview arms")
    if not any(a == "--out" or a.startswith("--out=") for a in sys.argv):
        sys.argv += ["--out", str(runner.ROOT / "out/runs/breaker_shading_study_20261019")]
    runner.main()


if __name__ == "__main__":
    main()
