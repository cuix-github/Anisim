"""Export and measure the existing breaker surface; no solver or renderer calls.

Run with the repo Python. A bounded mayapy worker loads the saved approved control,
checks its triangle count, then removes MCF and mesher smoothing in sequence.
All stages use the full source crop, so the measured interior patch has no new
reconstruction boundary. Static meshes and normals stay in ignored out/.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
CACHE = ROOT / "out/runs/wavesea4/cache/breaking_section"
SOURCE = ROOT / "out/runs/breaker_shading_study_20261019/A/breaker.ma"
OUT = ROOT / "out/runs/breaker_surface_diagnosis_20261019_v2"
STAGES = [("control", 10, .1), ("no_mcf", 0, .1), ("raw", 0, 0.)]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def worker(out):
    os.environ["BIFROST_LIB_CONFIG_FILES"] = str(ROOT / "maya/pfflip_graphs_config.json")
    for flag in ("MAYA_DISABLE_CER", "MAYA_DISABLE_CIP", "MAYA_DISABLE_CLIC_IPM"):
        os.environ[flag] = "1"
    import numpy as np
    import maya.standalone
    t0 = time.perf_counter()
    maya.standalone.initialize(name="python")
    import maya.cmds as c
    import maya.api.OpenMaya as om
    c.file(str(SOURCE), open=True, force=True)
    c.currentTime(1)
    setup_s = time.perf_counter() - t0
    graph = "particleSurfaceShape"
    original_graph = json.loads(c.getAttr(graph + ".sc"))
    converter = c.createNode("bifrostGeoToMaya", name="diagnosisConverter")
    # This converter defaults to recomputing Maya normals. Preserve Bifrost's instead.
    c.setAttr(converter + ".transferNormals", True)
    c.connectAttr(graph + ".mesh", converter + ".bifrostGeo")
    live = c.createNode("mesh", name="diagnosisLiveShape")
    c.connectAttr(converter + ".mayaMesh[0]", live + ".inMesh")
    positions = np.load(CACHE / "liquid.0001.npy")
    radii = np.load(CACHE / "radius.0001.npy")
    metres = positions * .25
    # Interior curl patch, away from the crop's z=60/78 m cut faces.
    roi = ((metres[:, 0] >= 134) & (metres[:, 0] <= 152) &
           (metres[:, 1] >= 13) & (metres[:, 1] <= 23) &
           (metres[:, 2] >= 63) & (metres[:, 2] <= 75))
    indices = np.flatnonzero(roi)
    np.savez_compressed(out / "probes.npz", indices=indices, positions=positions[roi],
                        radii=radii[roi])
    result = {"setup_open_s": setup_s, "converter_transfer_normals": c.getAttr(converter + ".transferNormals"),
              "roi_metres": [[134, 152], [13, 23], [63, 75]], "probes": len(indices), "stages": {}}

    def mesh_fn(node):
        selection = om.MSelectionList(); selection.add(node)
        dag = selection.getDagPath(0)
        if dag.apiType() == om.MFn.kTransform:
            dag.extendToShape()
        return om.MFnMesh(dag)

    def arrays(fn):
        return dict(points=np.asarray(fn.getPoints(om.MSpace.kObject))[:, :3].astype(np.float32),
                    triangles=np.asarray(fn.getTriangles()[1], dtype=np.int32).reshape(-1, 3),
                    face_counts=np.asarray(fn.getVertices()[0], dtype=np.int32),
                    face_vertices=np.asarray(fn.getVertices()[1], dtype=np.int32),
                    normals=np.asarray(fn.getNormals(om.MSpace.kObject), dtype=np.float32),
                    normal_ids=np.asarray(fn.getNormalIds()[1], dtype=np.int32))

    for stage, mcf, smoothing in STAGES:
        ts = time.perf_counter()
        print("EVALUATING", stage, flush=True)
        if stage != "control":
            c.vnnNode(graph, "/particle_surface", setPortDefaultValues=("mcf_iterations", str(mcf)))
            c.vnnNode(graph, "/particle_surface", setPortDefaultValues=("smoothing", str(smoothing)))
        stage_graph = json.loads(c.getAttr(graph + ".sc"))
        (out / (stage + "_graph.json")).write_text(json.dumps(stage_graph, indent=2))
        fn = mesh_fn(live)
        data = arrays(fn)  # Pull the actual Bifrost output through the public Maya converter.
        evaluation_s = time.perf_counter() - ts
        assert len(data["points"]) and np.isfinite(data["points"]).all()
        if stage == "control":
            assert len(data["triangles"]) == 292516, "Control disagrees with prior Arnold triangle count"
            assert stage_graph == original_graph
        duplicate = c.duplicate(c.listRelatives(live, parent=True)[0], name="baked_" + stage)[0]
        c.delete(duplicate, constructionHistory=True)
        copy_fn = mesh_fn(duplicate)
        copy_fn.lockFaceVertexNormals(np.repeat(np.arange(len(data["face_counts"])),
                                                data["face_counts"]).tolist(), data["face_vertices"].tolist())
        copy_data = arrays(copy_fn)
        for key in data:
            assert np.array_equal(data[key], copy_data[key]), (stage, key, "bake differs")
        c.sets(duplicate, edit=True, forceElement="waterSG")
        c.setAttr(duplicate + ".visibility", int(stage == "control"))
        np.savez_compressed(out / (stage + "_mesh.npz"), **data)
        # Exact closest points on the converted mesh. The signed projection below
        # is a local normal-side test, NOT a global winding-number inside test.
        tq = time.perf_counter()
        intersector = om.MMeshIntersector(); intersector.create(fn.object())
        near = np.empty((len(indices), 3), np.float32)
        normals = np.empty_like(near)
        for j, p in enumerate(positions[roi]):
            hit = intersector.getClosestPoint(om.MPoint(float(p[0]), float(p[1]), float(p[2])))
            near[j] = (hit.point.x, hit.point.y, hit.point.z)
            normals[j] = (hit.normal.x, hit.normal.y, hit.normal.z)
        np.savez_compressed(out / (stage + "_closest.npz"), points=near, normals=normals)
        result["stages"][stage] = dict(mcf_iterations=mcf, mesher_smoothing=smoothing,
            vertices=len(data["points"]), triangles=len(data["triangles"]),
            evaluation_s=evaluation_s, closest_query_s=time.perf_counter() - tq,
            wall_s=time.perf_counter() - ts)
        print("STAGE", stage, json.dumps(result["stages"][stage]), flush=True)
        (out / "worker_progress.json").write_text(json.dumps(result, indent=2))

    # Keep camera, lights, materials and three baked meshes; remove graph evaluation.
    c.delete("particleSurface", "particleMist", c.listRelatives(live, parent=True)[0], converter)
    assert not c.ls(type="bifrostGraphShape")
    c.file(rename=str(out / "baked_breaker.mb"))
    c.file(save=True, type="mayaBinary", force=True)
    tr = time.perf_counter()
    c.file(str(out / "baked_breaker.mb"), open=True, force=True)
    assert not c.ls(type="bifrostGraphShape")
    for stage, _, _ in STAGES:
        saved = np.load(out / (stage + "_mesh.npz"))
        reloaded = arrays(mesh_fn("baked_" + stage))
        for key, value in reloaded.items():
            if key in ("normals", "normal_ids"):
                continue  # Maya may compact/reorder the normal table on disk.
            assert np.array_equal(saved[key], value), (stage, key, "reload differs")
        before = saved["normals"][saved["normal_ids"]]
        after = reloaded["normals"][reloaded["normal_ids"]]
        error = float(np.max(np.abs(before - after)))
        result["stages"][stage]["reload_corner_normal_max_error"] = error
        assert np.allclose(before, after, rtol=0, atol=2e-6), (stage, error, "reload normals differ")
    result["baked_reload_verify_s"] = time.perf_counter() - tr
    result["baked_bytes"] = (out / "baked_breaker.mb").stat().st_size
    result["worker_wall_s"] = time.perf_counter() - t0
    result["status"] = "passed"
    (out / "worker_result.json").write_text(json.dumps(result, indent=2))
    print("BAKED_RELOAD_PASS", result["baked_reload_verify_s"], flush=True)
    maya.standalone.uninitialize()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker)
        return
    assert not git("status", "--porcelain"), "Commit the instrument before running"
    assert not OUT.exists(), "Never overwrite a previous experiment"
    prior = json.loads((SOURCE.parent / "receipt.json").read_text())
    hashes = {p.name: digest(p) for p in sorted(CACHE.iterdir()) if p.is_file()}
    assert hashes == prior["input_sha256"] and digest(SOURCE) == prior["scene_sha256"]
    for name in ("maya/particle_surface.json", "maya/particle_mist.json"):
        assert git("rev-parse", "HEAD:" + name) == git("rev-parse", prior["git_sha"] + ":" + name)
    OUT.mkdir(parents=True)
    receipt = {"git_sha": git("rev-parse", "HEAD"), "input_sha256": hashes,
               "scene_sha256": digest(SOURCE), "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
               "status": "running", "timeout_s": 600}
    (OUT / "receipt.json").write_text(json.dumps(receipt, indent=2))
    ts = time.perf_counter()
    with (OUT / "worker.log").open("w") as log:
        try:
            run = subprocess.run([r"C:\Program Files\Autodesk\Maya2027\bin\mayapy.exe", str(Path(__file__)),
                                  "--worker", str(OUT)], cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, timeout=600)
            receipt["exit_code"] = run.returncode
            assert run.returncode == 0, "Worker failed; inspect worker.log"
            result = json.loads((OUT / "worker_result.json").read_text())
            assert result["status"] == "passed"
            assert hashes == {p.name: digest(p) for p in sorted(CACHE.iterdir()) if p.is_file()}
            assert digest(SOURCE) == receipt["scene_sha256"]
            receipt["status"] = "passed"
        finally:
            receipt["wall_s"] = time.perf_counter() - ts
            if receipt["status"] == "running":
                receipt["status"] = "failed"
            (OUT / "receipt.json").write_text(json.dumps(receipt, indent=2))
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
