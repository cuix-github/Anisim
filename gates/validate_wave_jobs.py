"""Cache birth/death, restart corruption and frozen-runtime gates; optional Bifrost resume."""
import argparse
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "experiments"), str(ROOT / "graphs"), str(ROOT / "gates")]
import numpy as np
from wave_features import exact_frame, atomic_json
from wave_resume import checkpoint, sha, validate_plan
from loop3d import FINAL_STATE, load_final_state, run_bifrost3


def small_state(d):
    d.mkdir(parents=True, exist_ok=True)
    np.save(d / "final.npy", np.array([[1, 2, .5], [2, 3, .5]], np.float32))
    for name in FINAL_STATE:
        a = np.zeros((2, 3) if name == "velocities" else (2,), np.float32)
        if name == "mass":
            a[:] = 1000
        np.save(d / f"final_{name}.npy", a)


class Files(unittest.TestCase):
    def test_partial_and_nonfinite_checkpoint_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            small_state(d)
            self.assertEqual(checkpoint(d)["particles"], 2)
            (d / "final_radius.npy").unlink()
            with self.assertRaises(FileNotFoundError):
                checkpoint(d)
            np.save(d / "final_radius.npy", np.array([0, np.nan], np.float32))
            with self.assertRaisesRegex(ValueError, "Nonfinite"):
                checkpoint(d)

    def test_birth_death_and_fresh_spray_radius(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            src, dest = d / "dumps", d / "cache"
            src.mkdir()
            for k, n in enumerate((3, 5)):
                a = dict(frame=np.arange(n * 3, dtype=np.float32).reshape(n, 3),
                         vel=np.ones((n, 3), np.float32), ph=np.ones(n, np.float32),
                         rad=np.zeros(n, np.float32), spr=np.zeros(n, np.float32),
                         esc=np.zeros(n, np.float32))
                a["spr"][-1] = a["esc"][-1] = 1
                for name, array in a.items():
                    np.save(src / f"{name}.{k:04d}.npy", array)
                rec = exact_frame(src, k, dest, k + 1, cells_per_metre=2, layers=True)
                self.assertEqual(rec["full"], n)
                self.assertEqual(rec["spray"], 1)
                self.assertEqual(rec["fresh_spray_without_radius"], 1)
                np.testing.assert_array_equal(np.load(dest / f"liquid.{k+1:04d}.npy"), a["frame"] / 2)
                self.assertEqual(np.load(dest / f"radius.{k+1:04d}.npy")[-1], 0)
                self.assertEqual(np.load(dest / f"spray.{k+1:04d}.npy")[-1], 1)

    def test_frozen_runtime_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "new_run"
            runtime = root / "inputs/runtime"
            runtime.mkdir(parents=True)
            code = runtime / "worker.py"
            code.write_text("original")
            gate = root / "gates.json"
            atomic_json(gate, {"passed": True})
            atomic_json(root / "inputs/runtime_manifest.json", {"worker.py": sha(code)})
            plan = dict(root=str(root), tag="new_run", source=str(Path(tmp) / "original"),
                        runtime=str(runtime), gates=str(gate), from_s=18, to_s=30, grid=[688, 96, 320])
            self.assertEqual(validate_plan(plan), root.resolve())
            code.write_text("changed")
            with self.assertRaisesRegex(ValueError, "Frozen runtime"):
                validate_plan(plan)


def physics(root, source):
    import validate_spray_node as spray
    spray.OUT = str(root / "spray")
    Path(spray.OUT).mkdir(parents=True, exist_ok=False)
    records = {}
    for name, fn in (("off", spray.check_off), ("lockstep", spray.check_lockstep),
                     ("loop", spray.check_loop)):
        good, record = fn()
        records[name] = record
        if not good:
            raise AssertionError(f"Spray {name} gate failed")
    # Resume with all requested settings and the actual beach/wave zone, using
    # the small nz=1 scene. This is stronger than just checking fields exist.
    import wave_basin as basin
    plan = json.loads((source / "inputs/plan.json").read_text())
    last = json.loads((source / "receipts/summary.json").read_text())["chunks"][-1]
    offset = last["offset"] + last["steps"]
    run = basin.Run(plan["cells_per_metre"])
    run.grid, run.dt = tuple(plan["grid"]), plan["dt"]
    init = load_final_state(str(source / f"dumps/c{last['chunk']}"))
    common = dict(extra_ports=plan["ports"], fused=True, final_state=True,
                  dump_every=1, solid=(basin.beach_phi(run, basin.scene(run)), 2),
                  waves=np.asarray(plan["waves"], np.float32))
    paths = {}
    for tag, steps, state in (("full", 4, init), ("first", 2, init)):
        paths[tag], _, _ = run_bifrost3(tag, *run.grid, None, steps, run.dt,
            init_arrays=state, step_offset=offset, out_dir=str(root / tag), **common)
    paths["resumed"], _, _ = run_bifrost3("resumed", *run.grid, None, 2, run.dt,
        init_arrays=load_final_state(paths["first"]), step_offset=offset + 2,
        out_dir=str(root / "resumed"), **common)
    for name in ("final",) + tuple("final_" + x for x in FINAL_STATE):
        np.testing.assert_array_equal(np.load(Path(paths["full"]) / (name + ".npy")),
                                      np.load(Path(paths["resumed"]) / (name + ".npy")), err_msg=name)
    records["beach_zone_spray_resume_bit_exact"] = True
    return records


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--physics", action="store_true")
    ap.add_argument("--resume-source", type=Path,
                    help="Completed wave_2d run used for the beach/zone restart gate")
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()
    if args.physics and not args.resume_source:
        ap.error("--physics requires --resume-source")
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=False)
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Files))
    record = dict(passed=result.wasSuccessful(), file_tests=result.testsRun)
    if record["passed"] and args.physics:
        record["physics"] = physics(args.out, args.resume_source.resolve())
    atomic_json(args.out / "receipt.json", record)
    sys.exit(0 if record["passed"] else 1)
