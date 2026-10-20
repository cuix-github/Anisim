"""Prepare and supervise a separate 18--30 s spray continuation.

prepare writes an immutable runtime snapshot and a plan; it never simulates.
run executes sequential short chunks, publishing only validated checkpoints.
Failures stop in needs_attention for the scheduled Codex check to diagnose.
No old files are removed, no rendering is chained, no solver setting is relaxed.
"""
import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from datetime import datetime

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "graphs"))
from wave_features import atomic_json, spray_ports
from loop3d import FINAL_STATE, load_final_state, run_bifrost3, DUMPS
import wave_basin as basin
import numpy as np

GIB = 1024 ** 3


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while data := f.read(8 * 1024 * 1024):
            h.update(data)
    return h.hexdigest()


def checkpoint(path, scan=True):
    """Validate complete restart state, reading bounded slices, never legacy fills."""
    path = Path(path)
    fields = ["final"] + ["final_" + x for x in FINAL_STATE]
    arrays = {key: np.load(path / (key + ".npy"), mmap_mode="r") for key in fields}
    n = len(arrays["final"])
    if not n:
        raise ValueError("Empty checkpoint")
    for key, a in arrays.items():
        shape = (n, 3) if key in ("final", "final_velocities") else (n,)
        if a.dtype != np.float32 or a.shape != shape:
            raise ValueError(f"Invalid checkpoint channel {key}: {a.shape}, {a.dtype}")
        if scan:
            for i in range(0, n, 250_000):
                if not np.isfinite(a[i:i + 250_000]).all():
                    raise ValueError(f"Nonfinite checkpoint channel {key}")
    if scan:
        for i in range(0, n, 250_000):
            spr = arrays["final_spray"][i:i + 250_000]
            mass = arrays["final_mass"][i:i + 250_000]
            rad = arrays["final_radius"][i:i + 250_000]
            if (rad < 0).any() or not np.array_equal(mass == 0, spr > 0.5):
                raise ValueError("Invalid spray mass/radius in checkpoint")
    return dict(particles=n, files={key + ".npy": {
        "bytes": (path / (key + ".npy")).stat().st_size,
        "sha256": sha(path / (key + ".npy")) if scan else None} for key in fields})


def memory():
    class M(ctypes.Structure):
        _fields_ = [("length", ctypes.c_ulong), ("load", ctypes.c_ulong)] + [
            (x, ctypes.c_ulonglong) for x in
            ("total_ram", "free_ram", "total_commit", "free_commit", "total_va", "free_va", "extended")]
    m = M()
    m.length = ctypes.sizeof(m)
    if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m)):
        raise ctypes.WinError()
    return dict(free_ram_gib=m.free_ram / GIB, free_commit_gib=m.free_commit / GIB)


def prepare(args):
    root = ROOT / "out" / "runs" / args.tag
    root.mkdir(parents=True, exist_ok=False)
    source = (ROOT / "out/runs/wavesea4/dumps/c2").resolve()
    source_check = checkpoint(source)
    snapshot = root / "inputs" / "runtime"
    for folder in ("experiments", "graphs", "reference2d", "maya"):
        for f in (ROOT / folder).rglob("*"):
            if f.is_file() and f.suffix in (".py", ".json", ".bat", ".ps1") and "__pycache__" not in f.parts:
                dest = snapshot / f.relative_to(ROOT)
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(f, dest)
    pack = Path("nodes/build/PFFlipSpike-1.0.0")
    shutil.copytree(ROOT / pack, snapshot / pack)
    manifests = {str(f.relative_to(snapshot)).replace("\\", "/"): sha(f)
                 for f in snapshot.rglob("*") if f.is_file()}
    atomic_json(root / "inputs/runtime_manifest.json", manifests)
    run = basin.Run(4)
    ports = dict(basin.ports(run, basin.scene(run)), **spray_ports(4))
    # 24 half-second chunks. Keep all inputs, dumps and checkpoints: no retention deletion.
    segments = [dict(index=i, offset=720 + i * 20, steps=20) for i in range(24)]
    plan = dict(version=1, tag=args.tag, root=str(root), source=str(source),
                source_check=source_check, runtime=str(snapshot), python=sys.executable,
                code_commit=subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                grid=list(run.grid), dt=0.025, cpm=4, from_s=18, to_s=30,
                dump_every=5, ports=ports, waves=basin.waves(run).tolist(), segments=segments,
                start_at="2026-10-20T07:00:00-04:00", threads=32,
                disk_reserve_gib=250, stall_seconds=1800, max_attempts_per_segment=3,
                gates=str((ROOT / "out/runs/wave_preflight_20261020/receipt.json").resolve()),
                storage_estimate_gib=500,
                meaning="Spray enabled from 18 s only; earlier history remains the original simulation")
    atomic_json(root / "inputs/plan.json", plan)
    atomic_json(root / "status.json", dict(state="prepared", start_at=plan["start_at"],
                                          completed_until_s=18, plan=str(root / "inputs/plan.json")))
    print(json.dumps(dict(plan=str(root / "inputs/plan.json"), estimated_gib=500,
                          source_particles=source_check["particles"]), indent=2))


def validate_plan(plan):
    root = Path(plan["root"]).resolve()
    if root.name != plan["tag"] or Path(plan["source"]).resolve().is_relative_to(root):
        raise ValueError("Invalid destination/source separation")
    if plan["from_s"] != 18 or plan["to_s"] != 30 or plan["grid"] != [688, 96, 320]:
        raise ValueError("Plan differs from the approved run")
    runtime = Path(plan["runtime"])
    manifests = json.loads((root / "inputs/runtime_manifest.json").read_text())
    for rel, digest in manifests.items():
        if sha(runtime / rel) != digest:
            raise ValueError(f"Frozen runtime changed: {rel}")
    gate = json.loads(Path(plan["gates"]).read_text())
    if gate.get("passed") is not True:
        raise ValueError("Preflight gates have not passed")
    return root


def worker(plan, segment, attempt):
    root = validate_plan(plan)
    spec = plan["segments"][segment]
    receipts = root / "receipts"
    source = Path(plan["source"]) if segment == 0 else Path(
        json.loads((receipts / f"c{segment - 1}.json").read_text())["checkpoint"])
    expected = plan["source_check"] if segment == 0 else json.loads(
        (receipts / f"c{segment - 1}.json").read_text())["validation"]
    # Validate the preceding completed state before trusting it.
    actual = checkpoint(source)
    if actual != expected:
        raise ValueError("Checkpoint no longer matches its completion receipt")
    d = root / "dumps" / f"c{segment}" / f"attempt_{attempt:02d}"
    d.mkdir(parents=True, exist_ok=False)
    run = basin.Run(plan["cpm"])
    init = load_final_state(str(source))
    t0 = time.monotonic()
    _, _, diag = run_bifrost3(plan["tag"], *plan["grid"], None, spec["steps"], plan["dt"],
        extra_ports=plan["ports"], init_arrays=init, release_init=True,
        final_state=True, step_offset=spec["offset"], dump_every=plan["dump_every"],
        solid=(basin.beach_phi(run, basin.scene(run)), 2), dumps=DUMPS + ("vel",),
        fused=True, out_dir=str(d), waves=np.asarray(plan["waves"], np.float32),
        log_prefix=str(root / "logs" / f"c{segment}_a{attempt}_bifcmd"))
    if not np.isfinite(diag).all() or (diag[:, 0] < 0).any() or (diag[:, 1] > 1.05e-6).any():
        raise ValueError("Invalid solver diagnostics")
    validation = checkpoint(d)
    record = dict(segment=segment, attempt=attempt, checkpoint=str(d), validation=validation,
                  offset=spec["offset"], steps=spec["steps"], dt=plan["dt"],
                  dump_every=plan["dump_every"], end_s=(spec["offset"] + spec["steps"]) * plan["dt"],
                  wall_s=time.monotonic() - t0, iterations_max=float(diag[:, 0].max()),
                  residual_max=float(diag[:, 1].max()), runtime_commit=plan["code_commit"])
    atomic_json(receipts / f"c{segment}.json", record)
    print(json.dumps({k: v for k, v in record.items() if k != "validation"}), flush=True)


def stop_tree(process):
    if process.poll() is None:
        subprocess.run(["taskkill", "/PID", str(process.pid), "/T", "/F"], capture_output=True, check=False)
        process.wait(timeout=30)


def owned_processes(runtime):
    """Detect an orphan from this exact frozen runtime before another launch."""
    script = ("Get-CimInstance Win32_Process | Where-Object { "
              "$_.Name -in @('python.exe','bifcmd.exe','cmd.exe') } | "
              "Select-Object ProcessId,Name,CommandLine | ConvertTo-Json -Compress")
    result = subprocess.run(["powershell.exe", "-NoProfile", "-Command", script],
                            capture_output=True, text=True, check=True)
    rows = json.loads(result.stdout) if result.stdout.strip() else []
    if isinstance(rows, dict):
        rows = [rows]
    marker = str(Path(runtime).resolve()).replace("\\", "/").lower()
    return [row for row in rows if row["ProcessId"] != os.getpid() and marker in
            (row.get("CommandLine") or "").replace("\\", "/").lower()]


def supervise(plan):
    import msvcrt
    root = validate_plan(plan)
    if datetime.now().astimezone() < datetime.fromisoformat(plan["start_at"]):
        raise ValueError("The approved start time has not arrived")
    lockfile = root / "supervisor.lock"
    lockfile.touch(exist_ok=True)
    lock = open(lockfile, "r+b")
    lock.seek(0)
    if lockfile.stat().st_size == 0:
        lock.write(b"0")
        lock.flush()
    lock.seek(0)
    msvcrt.locking(lock.fileno(), msvcrt.LK_NBLCK, 1)
    status = dict(state="running", pid=os.getpid(), completed_until_s=18, tag=plan["tag"])
    child = None
    try:
        orphans = owned_processes(plan["runtime"])
        if orphans:
            raise ValueError(f"Existing processes for this runtime; inspect before restart: {orphans}")
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000001)
        for spec in plan["segments"]:
            i = spec["index"]
            receipt = root / "receipts" / f"c{i}.json"
            if receipt.exists():
                status["completed_until_s"] = json.loads(receipt.read_text())["end_s"]
                continue
            free = shutil.disk_usage(root).free / GIB
            if free < plan["disk_reserve_gib"] + 25:
                raise ValueError(f"Insufficient disk for another chunk: {free:.1f} GiB")
            if memory()["free_commit_gib"] < 12:
                raise ValueError("Insufficient available commit before chunk")
            cdir = root / "dumps" / f"c{i}"
            previous = sorted(cdir.glob("attempt_*")) if cdir.exists() else []
            attempt = max([int(p.name.split("_")[-1]) for p in previous], default=-1) + 1
            if attempt >= plan["max_attempts_per_segment"]:
                raise ValueError(f"Attempt budget exhausted at segment {i}; investigation required")
            logs = root / "logs"
            logs.mkdir(exist_ok=True)
            env = dict(os.environ, OMP_NUM_THREADS=str(plan["threads"]), PYTHONUNBUFFERED="1")
            cmd = [plan["python"], str(Path(plan["runtime"]) / "experiments/wave_resume.py"),
                   "worker", str(root / "inputs/plan.json"), "--segment", str(i), "--attempt", str(attempt)]
            with open(logs / f"c{i}_a{attempt}.log", "w") as log:
                child = subprocess.Popen(cmd, cwd=plan["runtime"], env=env, stdout=log, stderr=log,
                                          creationflags=subprocess.CREATE_NO_WINDOW |
                                          subprocess.BELOW_NORMAL_PRIORITY_CLASS)
                status.update(segment=i, attempt=attempt, worker_pid=child.pid, error=None)
                d = cdir / f"attempt_{attempt:02d}"
                progress_time, last_step, start = time.monotonic(), -1, time.monotonic()
                while child.poll() is None:
                    diagnostics = sorted(d.glob("diag.*.npy")) if d.exists() else []
                    step = len(diagnostics)
                    if step != last_step:
                        progress_time, last_step = time.monotonic(), step
                    status.update(checked_at=datetime.now().astimezone().isoformat(),
                                  steps_in_chunk=step, memory=memory(),
                                  free_disk_gib=shutil.disk_usage(root).free / GIB)
                    atomic_json(root / "status.json", status)
                    with open(logs / "monitor.jsonl", "a") as f:
                        f.write(json.dumps(status) + "\n")
                    if status["memory"]["free_commit_gib"] < 2 or status["free_disk_gib"] < plan["disk_reserve_gib"]:
                        stop_tree(child)
                        raise ValueError("Resource reserve reached; stopped only this worker tree")
                    if time.monotonic() - progress_time > plan["stall_seconds"]:
                        stop_tree(child)
                        raise ValueError("No new completed step for 30 minutes; investigate before retry")
                    time.sleep(15)
                if child.returncode != 0 or not receipt.exists():
                    raise ValueError(f"Worker c{i} attempt {attempt} failed ({child.returncode}); see {log.name}")
            status["completed_until_s"] = json.loads(receipt.read_text())["end_s"]
            status["last_chunk_wall_s"] = time.monotonic() - start
            atomic_json(root / "status.json", status)
        status.update(state="complete", completed_until_s=30, finished_at=datetime.now().astimezone().isoformat())
    except Exception as e:
        if child is not None:
            stop_tree(child)
        status.update(state="needs_attention", error=str(e), checked_at=datetime.now().astimezone().isoformat())
        raise
    finally:
        atomic_json(root / "status.json", status)
        ctypes.windll.kernel32.SetThreadExecutionState(0x80000000)
        lock.close()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    sub = ap.add_subparsers(dest="command", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--tag", required=True)
    for command in ("run", "worker"):
        p = sub.add_parser(command)
        p.add_argument("plan")
        if command == "worker":
            p.add_argument("--segment", type=int, required=True)
            p.add_argument("--attempt", type=int, required=True)
    args = ap.parse_args()
    if args.command == "prepare":
        if not args.tag.replace("_", "").replace("-", "").isalnum():
            ap.error("tag must be a simple directory name")
        prepare(args)
    else:
        plan = json.loads(Path(args.plan).read_text())
        if args.command == "worker":
            worker(plan, args.segment, args.attempt)
        else:
            supervise(plan)


if __name__ == "__main__":
    main()
