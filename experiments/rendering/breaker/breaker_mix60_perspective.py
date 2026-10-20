"""One authorized 1080p perspective view of the frozen 24.9 s 60/40 foam scene."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import breaker_foam_render as render
import breaker_deep_water as base

SOURCE = base.ROOT / "out/runs/breaker_foam_mix60_20261019/C"
OUT = base.ROOT / "out/runs/breaker_mix60_perspective_20261020"
RESOLUTION = (1920, 1080)


def read(path):
    return json.loads(path.read_text())


def prepare_worker():
    import maya.standalone
    maya.standalone.initialize(name="python")
    import maya.cmds as c
    try:
        c.file(str(SOURCE / "candidate.mb"), open=True, force=True)
        before = read(SOURCE / "state.json")
        assert not c.ls(type="bifrostGraphShape")
        assert base.geometry(c) == before["geometry"]
        for node in render.NODE_KEYS:
            actual = {attr: c.getAttr(node + "." + attr) for attr in before[node]}
            assert json.loads(json.dumps(actual)) == before[node], node
        assert base.sha(SOURCE / "foam.ass") == before["foam_ass_sha256"]
        target = np.array([143.0, 12.5, 48.0])
        old_target = np.asarray(before["view"]["target_m"])
        direction = np.asarray(before["persp"]["translate"][0]) / 4 - old_target
        direction /= np.linalg.norm(direction)
        distance, width = 100.0, 56.0
        eye = target + direction * distance
        c.setAttr("persp.translate", *(eye * 4))
        c.setAttr("perspShape.orthographic", False)
        c.setAttr("perspShape.horizontalFilmAperture", 36 / 25.4)
        c.setAttr("perspShape.verticalFilmAperture", 20.25 / 25.4)
        c.setAttr("perspShape.focalLength", 36 * distance / width)
        c.setAttr("perspShape.filmFit", 1)  # horizontal fit
        c.setAttr("perspShape.centerOfInterest", distance * 4)
        c.setAttr("perspShape.depthOfField", False)
        for attr, value in (("width", RESOLUTION[0]), ("height", RESOLUTION[1]),
                            ("pixelAspect", 1), ("deviceAspectRatio", 16 / 9)):
            c.setAttr("defaultResolution." + attr, value)
        c.setAttr("foamLayerShape.dso", (OUT / "foam.ass").as_posix(), type="string")
        c.setAttr("opticalEXR.prefix", (OUT / "optics").as_posix(), type="string")
        c.setAttr("defaultArnoldRenderOptions.log_filename", (OUT / "arnold.log").as_posix(), type="string")
        state = dict(before)
        for node in ("persp", "perspShape"):
            attrs = list(before[node])
            if node == "perspShape":
                attrs += ["focalLength", "horizontalFilmAperture", "verticalFilmAperture", "filmFit", "depthOfField"]
            state[node] = {attr: c.getAttr(node + "." + attr) for attr in attrs}
        state["view"] = dict(time_s=24.9, target_m=target.tolist(), eye_m=eye.tolist(),
                             distance_m=distance, target_plane_width_m=width,
                             height_above_still_water_m=float(eye[1] - 12),
                             resolution=RESOLUTION, perspective=True,
                             extent="Existing 50x68 m reconstruction; not the entire 172x80 m basin")
        assert base.geometry(c) == before["geometry"]
        base.save(OUT / "state.json", state)
        c.file(rename=str(OUT / "candidate.mb"))
        c.file(save=True, type="mayaBinary", force=True)
        c.arnoldExportAss(f=str(OUT / "scene_check.ass"), mask=24)
        base.save(OUT / "prepare.json", dict(status="passed",
                 scene_sha256=base.sha(OUT / "candidate.mb"),
                 geometry_identical=True, shading_identical=True, foam_identical=True,
                 view=state["view"]))
        base.exit_worker(0)
    except BaseException:
        base.save(OUT / "prepare.json", dict(status="failed", error=traceback.format_exc()))
        base.exit_worker(1)


def prepare():
    assert not base.git("status", "--porcelain")
    OUT.mkdir(exist_ok=False)
    shutil.copy2(SOURCE / "foam.ass", OUT / "foam.ass")
    env = dict(os.environ, MIX60_PERSPECTIVE_WORKER="prepare")
    with (OUT / "prepare.log").open("w") as log:
        subprocess.run([render.MAYAPY, str(Path(__file__).resolve())], cwd=base.ROOT, env=env,
                       stdout=log, stderr=log, check=True, timeout=180)
    assert read(OUT / "prepare.json")["status"] == "passed"
    assert read(OUT / "state.json")["perspShape"]["orthographic"] is False


def launch():
    assert not base.git("status", "--porcelain")
    prepared = read(OUT / "prepare.json")
    assert prepared["status"] == "passed"
    assert prepared["scene_sha256"] == base.sha(OUT / "candidate.mb")
    assert not (OUT / "receipt.json").exists(), "Only one new render authorized"
    # Use the known Maya launch quoting. Allow the separate lightweight 2D viewer.
    launcher = (SOURCE / "launch.ps1").read_text()
    launcher = launcher.replace("FOAM_MIX60_SCRIPT", "MIX60_PERSPECTIVE_SCRIPT")
    launcher = "\n".join(line for line in launcher.splitlines()
                         if "Another Maya process is open" not in line)
    launcher = launcher.replace("900000", "2700000").replace("900 seconds", "2700 seconds")
    (OUT / "launch.ps1").write_text(launcher)
    env = {k: v for k, v in os.environ.items() if not k.startswith("SURF_")}
    env.update(MIX60_PERSPECTIVE_WORKER="render", MIX60_PERSPECTIVE_SCRIPT=str(Path(__file__).resolve()))
    receipt = dict(status="running", git_sha=base.git("rev-parse", "HEAD"),
                   authorized_render_count=1, started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    base.save(OUT / "receipt.json", receipt)
    t0 = time.monotonic()
    try:
        with (OUT / "launch.log").open("w") as log:
            subprocess.run(["powershell.exe", "-NoProfile", "-File", str(OUT / "launch.ps1")],
                           env=env, cwd=base.ROOT, stdout=log, stderr=log, check=True, timeout=2740)
        worker = read(OUT / "worker_result.json")
        assert worker["status"] == "passed"
        receipt.update(status="passed", worker=worker,
                       png_sha256=base.sha(OUT / "beauty.png"),
                       exr_sha256=base.sha(OUT / worker["exr"]))
    except BaseException:
        receipt.update(status="failed", error=traceback.format_exc())
        raise
    finally:
        receipt["process_s"] = time.monotonic() - t0
        base.save(OUT / "receipt.json", receipt)
    render.check(OUT, resolution=RESOLUTION)


if __name__ == "__main__":
    mode = os.environ.get("MIX60_PERSPECTIVE_WORKER")
    if mode:
        prepare_worker() if mode == "prepare" else render.render(OUT, resolution=RESOLUTION)
    else:
        parser = argparse.ArgumentParser()
        parser.add_argument("action", choices=("prepare", "render"))
        {"prepare": prepare, "render": launch}[parser.parse_args().action]()
