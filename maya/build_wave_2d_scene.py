"""Build a particle playback scene with the actual configured solver graph beside it.

mayapy maya/build_wave_2d_scene.py <run directory>
The solver compound has no connected outputs: timeline playback evaluates only
cached particles. Its inputs are complete for manual, explicitly requested evaluation.
"""
import json
import os
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CONFIGS = [ROOT / "nodes/build/PFFlipSpike-1.0.0/PFFlipSpikePackConfig.json",
           ROOT / "maya/pfflip_sim_config.json", ROOT / "maya/pfflip_graphs_config.json"]
os.environ["BIFROST_LIB_CONFIG_FILES"] = ";".join(map(str, CONFIGS))
for key in ("MAYA_DISABLE_CER", "MAYA_DISABLE_CIP", "MAYA_DISABLE_CLIC_IPM"):
    os.environ[key] = "1"


def main():
    run = Path(sys.argv[1]).resolve()
    plan = json.loads((run / "inputs/plan.json").read_text())
    length, height = plan.get("domain_m", [172, 24])[:2]
    nz = plan["grid"][2]
    basin_width = plan.get("domain_m", [172, 24, 1])[2] if nz > 1 else 0
    toe, slope = plan.get("beach_toe_m", 36), plan.get("beach_slope", .1)
    cache = run / "cache/full"
    index = json.loads((cache / "cache.json").read_text())
    import numpy as np
    # Use identity DAG transforms: Bifrost viewport glyph sizes/transforms must
    # use the same centimetre coordinates as Maya's internal representation.
    display_cache = run / "cache/maya_cm"
    for layer in ("water", "droplets", "spray", "air"):
        dest = display_cache / layer
        dest.mkdir(parents=True, exist_ok=True)
        for frame in range(1, index["frames"] + 1):
            src = cache / layer / f"liquid.{frame:04d}.npy"
            np.save(dest / src.name, np.load(src) * np.float32(100))
    inspection = run / "receipts/particle_audit.json"
    start_frame = (json.loads(inspection.read_text())["best_simultaneous_frame"]["frame"]
                   if inspection.exists() else min(600, index["frames"]))
    events = np.concatenate([np.load(cache / layer / f"liquid.{start_frame:04d}.npy")
                             for layer in ("droplets", "spray")])
    events = events[events[:, 0] > 34]
    if len(events):
        lo, hi = np.percentile(events[:, 0], [2, 98])
        detail_x, detail_width = float((lo + hi) / 2), float(np.clip(hi - lo + 20, 60, 150))
    else:
        detail_x, detail_width = 129, 90
    import maya.standalone
    maya.standalone.initialize(name="python")
    import maya.cmds as c
    c.loadPlugin("bifrostGraph")
    c.file(new=True, force=True)
    fps = index.get("fps", 24)
    c.currentUnit(linear="m", time="film" if fps == 24 else f"{fps}fps")
    for layer, color, size in (("water", (0.06, 0.4, 0.85), 0.12),
                               ("droplets", (0.35, 0.95, 1.0), 0.18),
                               ("spray", (1.0, 0.65, 0.1), 0.25),
                               ("air", (0.4, 0.4, 0.4), 0.04)):
        shape = c.createNode("bifrostGraphShape", name=layer + "ParticlesShape")
        parent = c.listRelatives(shape, parent=True)[0]
        parent = c.rename(parent, layer + "_cached_particles")
        shape = c.listRelatives(parent, shapes=True)[0]
        node = c.vnnCompound(shape, "/", addNode="PFFlipGraphs,User::PFFlip,particle_play")[0]
        c.vnnNode(shape, "/input", createOutputPort=("frame", "float"))
        c.vnnNode(shape, "/output", createInputPort=("points", "Object"))
        c.vnnConnect(shape, "/input.frame", "/" + node + ".frame")
        c.vnnConnect(shape, "/" + node + ".points", "/output.points")
        c.setAttr(shape + ".overrideEnabled", True)
        c.setAttr(shape + ".overrideRGBColors", True)
        c.setAttr(shape + ".overrideColorRGB", *color)
        for port, value in (("liquid_pattern", str(display_cache / layer / "liquid.####").replace("\\", "/")),
                            ("speed_pattern", str(cache / layer / "speed.####").replace("\\", "/")),
                            ("point_size", str(size * 100)), ("speed_scale", "20"),
                            ("slow_color", "{" + ",".join(map(str, color)) + "}"),
                            ("fast_color", "{" + ",".join(map(str, color)) + "}")):
            c.vnnNode(shape, "/" + node, setPortDefaultValues=(port, value))
        for port in ("slow_color", "fast_color"):
            actual = c.vnnNode(shape, "/" + node, queryPortDefaultValues=port)
            if not np.allclose(np.asarray(actual, float), color):
                raise RuntimeError(f"Playback color failed to roundtrip: {layer}, {port}")
        c.connectAttr("time1.outTime", shape + ".frame")
        if layer == "air":
            c.setAttr(shape + ".visibility", False)
    shape = c.createNode("bifrostGraphShape", name="waveSolverGraphShape")
    parent = c.rename(c.listRelatives(shape, parent=True)[0], "wave_solver_OPEN_IN_BIFROST")
    shape = c.listRelatives(parent, shapes=True)[0]
    node = c.vnnCompound(shape, "/", addNode="PFFlipSim,User::PFFlip,sim_3d_fused")[0]
    compound = json.loads((ROOT / "graphs/sim_3d_fused.json").read_text())["compounds"][0]
    manual = run / "manual_solver_output"
    manual.mkdir(exist_ok=True)
    values = dict(plan["ports"], nx=plan["grid"][0], ny=plan["grid"][1], nz=nz,
                  dt=plan["dt"], substeps=round(plan["seconds"] / plan["dt"]),
                  step_offset=0, dump_every=2, solid_mode=2)
    for port in compound["ports"]:
        name = port["portName"]
        if name.startswith("path_final_"):
            short = name[len("path_final_"):]
            values[name] = str(manual / ("final.npy" if short == "positions" else f"final_{short}.npy"))
        elif name.startswith("path_"):
            values[name] = str(run / "dumps/c0" / (name[5:] + ".npy"))
        elif name.endswith("_pattern"):
            prefix = "frame" if name == "pos_pattern" else name.removesuffix("_pattern")
            values[name] = str(manual / (prefix + ".####"))
    values.update(plan.get("force_paths", {}))
    for key, value in values.items():
        c.vnnNode(shape, "/" + node, setPortDefaultValues=(key, str(value).replace("\\", "/")))
    c.addAttr(shape, longName="purpose", dataType="string")
    c.setAttr(shape + ".purpose", "Actual custom solver graph. No outputs connected; playback reads cache only.",
              type="string")
    # A metre-scale beach outline provides scale without rendering a water surface.
    c.curve(name="beach_profile", d=1,
            p=[(0, 0, 0), (toe, 0, 0), (length, (length - toe) * slope, 0)])
    c.curve(name="still_water_12m", d=1, p=[(0, 12, 0), (min(length, toe + 12 / slope), 12, 0)])
    c.curve(name="domain_outline", d=1,
            p=[(0, 0, 0), (length, 0, 0), (length, height, 0), (0, height, 0), (0, 0, 0)])
    cameras = {"persp": "persp"}
    for name, center, width in (("overview", (length / 2, height / 2, 180), length * 1.08),
                                 ("breaker_detail", (detail_x, 12, 100), detail_width)):
        cam, camshape = c.camera(name=name)
        cameras[name] = cam
        c.setAttr(cam + ".translate", *center)
        c.setAttr(camshape + ".orthographic", True)
        c.setAttr(camshape + ".orthographicWidth", width)
        c.setAttr(camshape + ".farClipPlane", 1000)
        c.setAttr(camshape + ".nearClipPlane", 0.01)
    c.setAttr("persp.translate", detail_x, 12, 100)
    c.setAttr("persp.rotate", 0, 0, 0)
    c.setAttr("perspShape.orthographic", True)
    c.setAttr("perspShape.orthographicWidth", detail_width)
    c.setAttr("perspShape.farClipPlane", 1000)
    c.setAttr("perspShape.nearClipPlane", .01)
    if nz > 1:
        c.curve(name="far_beach_profile", d=1,
                p=[(0, 0, basin_width), (toe, 0, basin_width),
                   (length, (length - toe) * slope, basin_width)])
        for name, eye, target in (
                ("overview", (length * .5, 90, basin_width + 130), (length * .5, 10, basin_width / 2)),
                ("breaker_detail", (145, 50, basin_width + 55), (105, 12, basin_width / 2)),
                ("persp", (145, 50, basin_width + 55), (105, 12, basin_width / 2))):
            name = cameras[name]
            c.setAttr(name + ".translate", *eye)
            c.setAttr(c.listRelatives(name, shapes=True)[0] + ".orthographic", False)
            locator = c.spaceLocator()[0]
            c.setAttr(locator + ".translate", *target)
            constraint = c.aimConstraint(locator, name, aimVector=(0, 0, -1), upVector=(0, 1, 0),
                                        worldUpType="vector", worldUpVector=(0, 1, 0))[0]
            c.delete(constraint, locator)
    for panel in c.getPanel(type="modelPanel") or []:
        c.modelEditor(panel, edit=True, camera=cameras["breaker_detail"], grid=False,
                      displayAppearance="wireframe")
    c.playbackOptions(min=1, max=index["frames"], ast=1, aet=index["frames"],
                      playbackSpeed=1, maxPlaybackSpeed=1)
    # Force representative frames to evaluate in Maya and verify visible geometry.
    audit = []
    for frame in sorted(set((1, index["frames"] // 4, index["frames"] // 2,
                             3 * index["frames"] // 4, start_frame, index["frames"]))):
        t0 = time.perf_counter()
        c.currentTime(frame)
        layers = {}
        for layer in ("water", "droplets", "spray"):
            points = np.load(cache / layer / f"liquid.{frame:04d}.npy")
            if len(points):
                bounds = c.exactWorldBoundingBox(layer + "_cached_particles")
                expected = np.r_[points.min(axis=0), points.max(axis=0)]
                if not np.allclose(bounds, expected, atol=1e-4):
                    raise RuntimeError(f"Playback unit/bounds mismatch: {layer} frame {frame}")
                layers[layer] = dict(particles=len(points), bounds_m=bounds)
        audit.append(dict(frame=frame, evaluation_s=time.perf_counter() - t0, layers=layers))
    c.currentTime(start_frame)
    c.select(clear=True)
    viewport_setup = (
        "import maya.cmds as c\n"
        "def setup_wave_view():\n"
        "    for p in c.getPanel(type='modelPanel') or []:\n"
        "        c.modelEditor(p, edit=True, camera='persp', grid=False, "
        "displayAppearance='wireframe', displayLights='default')\n"
        "c.evalDeferred(setup_wave_view)\n")
    c.scriptNode(name="wave_viewport_setup", scriptType=1,
                 sourceType="python", beforeScript=viewport_setup)
    scene = run / (sys.argv[2] if len(sys.argv) > 2 else
                   ("scenes/wave_3d.ma" if nz > 1 else "scenes/wave_2d.ma"))
    scene.parent.mkdir(exist_ok=True)
    c.file(rename=str(scene))
    c.file(save=True, type="mayaAscii")
    (run / "receipts/maya_scene.json").write_text(json.dumps(
        dict(status="passed", scene=str(scene), graph=node, evaluations=audit,
             start_frame=start_frame, detail_center_x_m=detail_x,
             detail_width_m=detail_width, cameras=cameras), indent=2))
    launch = run / "open_in_maya.bat"
    launch.write_text('@echo off\nsetlocal\nset "BIFROST_LIB_CONFIG_FILES=' +
                     ";".join(map(str, CONFIGS)) + '"\nstart "" "C:\\Program Files\\Autodesk\\Maya2027\\bin\\maya.exe" "' +
                     str(scene) + '"\n', encoding="utf-8")
    print("SAVED", scene, flush=True)
    import ctypes
    ctypes.windll.kernel32.TerminateProcess(ctypes.windll.kernel32.GetCurrentProcess(), 0)


if __name__ == "__main__":
    main()
