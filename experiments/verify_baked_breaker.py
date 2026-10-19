"""Independently reopen static meshes and verify geometry/normals; never reconstruct."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
DEFAULT = ROOT / "out/runs/breaker_surface_diagnosis_20261019_v3"


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def worker(out):
    import numpy as np
    import maya.standalone
    maya.standalone.initialize(name="python")
    import maya.cmds as c
    import maya.api.OpenMaya as om
    started = time.perf_counter()
    c.file(str(out / "baked_breaker.mb"), open=True, force=True)
    assert not c.ls(type="bifrostGraphShape"), "Baked scene still has live graphs"
    checks = {}
    meshes = {}
    for stage in ("control", "no_mcf", "raw"):
        saved = np.load(out / (stage + "_mesh.npz"))
        selection = om.MSelectionList(); selection.add("baked_" + stage)
        dag = selection.getDagPath(0); dag.extendToShape()
        fn = om.MFnMesh(dag)
        meshes[stage] = fn
        points = np.asarray(fn.getPoints())[:, :3].astype(np.float32)
        tri = np.asarray(fn.getTriangles()[1], dtype=np.int32).reshape(-1, 3)
        normals = np.asarray(fn.getNormals(), dtype=np.float32)
        ids = np.asarray(fn.getNormalIds()[1], dtype=np.int32)
        assert np.array_equal(points, saved["points"]), stage
        assert np.array_equal(tri, saved["triangles"]), stage
        a, b = normals[ids], saved["normals"][saved["normal_ids"]]
        error = float(np.max(np.abs(a - b)))
        assert np.allclose(a, b, rtol=0, atol=2e-6), (stage, error)
        assert bool(c.getAttr("baked_" + stage + ".visibility")) == (stage == "control")
        checks[stage] = {"positions_exact": True, "triangles_exact": True,
                         "face_corner_normal_max_error": error, "triangles": len(tri)}
    # Reuse static geometry for correspondence; this does not evaluate any graph.
    query_start = time.perf_counter()
    for source, target in (("raw", "control"), ("control", "raw")):
        positions = np.load(out / (source + "_mesh.npz"))["points"]
        intersector = om.MMeshIntersector(); intersector.create(meshes[target].object())
        closest = np.empty_like(positions)
        normals = np.empty_like(positions)
        faces = np.empty(len(positions), np.int32)
        for i, p in enumerate(positions):
            hit = intersector.getClosestPoint(om.MPoint(float(p[0]), float(p[1]), float(p[2])))
            closest[i] = (hit.point.x, hit.point.y, hit.point.z)
            normals[i] = (hit.normal.x, hit.normal.y, hit.normal.z)
            faces[i] = hit.face
        np.savez_compressed(out / (source + "_to_" + target + ".npz"),
                            closest=closest, normals=normals, faces=faces)
    optics = {attr: c.getAttr("water." + attr) for attr in
              ("metalness", "specularIOR", "specularRoughness", "transmission",
               "transmissionColor", "transmissionDepth", "transmissionScatter", "transmissionExtraRoughness")}
    result = {"status": "passed", "mesh_checks": checks, "optics": optics,
              "static_correspondence_s": time.perf_counter() - query_start,
              "open_and_verify_s": time.perf_counter() - started,
              "baked_scene_sha256": digest(out / "baked_breaker.mb")}
    (out / "static_verification_worker.json").write_text(json.dumps(result, indent=2))
    print("STATIC_VERIFICATION_PASS", json.dumps(checks), flush=True)
    import ctypes
    sys.stdout.flush(); sys.stderr.flush()
    ctypes.windll.kernel32.TerminateProcess(ctypes.windll.kernel32.GetCurrentProcess(), 0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--out", type=Path, default=DEFAULT)
    args = parser.parse_args()
    if args.worker:
        worker(args.worker)
        return
    out = args.out.resolve()
    assert not subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT).strip()
    receipt = json.loads((out / "receipt.json").read_text())
    completed = json.loads((out / "worker_result.json").read_text())
    assert completed["status"] == "passed", "No completed geometry export to verify"
    assert not (out / "static_verification.json").exists()
    env = os.environ.copy()
    env["BIFROST_LIB_CONFIG_FILES"] = str(ROOT / "maya/pfflip_graphs_config.json")
    for key in ("MAYA_DISABLE_CER", "MAYA_DISABLE_CIP", "MAYA_DISABLE_CLIC_IPM"):
        env[key] = "1"
    start = time.perf_counter()
    with (out / "static_verification.log").open("w") as log:
        run = subprocess.run([r"C:\Program Files\Autodesk\Maya2027\bin\mayapy.exe", str(Path(__file__)),
                              "--worker", str(out)], env=env, stdout=log, stderr=subprocess.STDOUT, timeout=120)
    assert run.returncode == 0, "Static verifier failed; inspect static_verification.log"
    result = json.loads((out / "static_verification_worker.json").read_text())
    assert result["status"] == "passed"
    cache = ROOT / "out/runs/wavesea4/cache/breaking_section"
    assert {p.name: digest(p) for p in cache.iterdir() if p.is_file()} == receipt["input_sha256"]
    result.update(exit_code=run.returncode, process_wall_s=time.perf_counter() - start,
                  git_sha=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                  export_launcher_status=receipt["status"], export_launcher_exit_code=receipt["exit_code"])
    (out / "static_verification.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
