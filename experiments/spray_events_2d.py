"""Particle-only plunge: independent NumPy / actual Bifrost C++ loop.

No Arnold, production input, or cache deletion. Each --out must be a new job.
Float32 state boundaries match Bifrost's public array ports; NumPy performs its
own transfers, pressure solve, classification, emission and flight in float64.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/"reference2d"), str(ROOT/"graphs"), str(ROOT/"gates")]
from pfflip2d import Sim
from loop3d import run_bifrost3
from validate_spray_events import f32, ports
from validate_droplet_node import ESCAPE

NX, NY = 64, 80
D = f32(dict(Sim.DROPLET_DEFAULTS, integrator=2, dt_max=.005, metres_per_cell=1))
E = f32(dict(Sim.EMIT_DEFAULTS, model=2, grad_threshold=2, rate=.6,
             cloud_size=16, max_per_cell=128, max_active=40000, seed=137))


def mask(x, y):
    return (y < 18) | ((abs(x-32) < 7) & (y > 39) & (y < 64))


class TagPhaseDiagnostic(Sim):
    """Only for isolating the existing graph's weighted-tag pressure coefficient.

    This is not Eq.7 and is not a replacement default for the reference.
    No C++ projection or particle state is used to advance this reference.
    """
    def _face_phase(self):
        assert not self.st and not self.adapt and not self.has_solids
        valid = (~self.spray).astype(float)
        fu,fv = self._face_frames()
        du,nu = self._splat(self.u.shape,fu,[valid,valid*self.typ])
        dv,nv = self._splat(self.v.shape,fv,[valid,valid*self.typ])
        eps = float(np.float32(1e-6))
        return np.clip(nu/(du+eps),0,1),np.clip(nv/(dv+eps),0,1)


def make_sim(model=2, phase_law="mass"):
    cls = Sim if phase_law=="mass" else TagPhaseDiagnostic
    s = cls(NX, NY, alpha_l=float(np.float32(.97)), alpha_g=float(np.float32(.9)),
            g=float(np.float32(-9.8)), escape=True, esc_phi=float(np.float32(.3)),
            droplet=D, emit=dict(E, model=model), spray_cfl=.5)
    s.seed(mask, seed=1)
    s.calibrate()
    s.rho0_face = float(np.float32(s.rho0_face))
    s.pos = s.pos.astype(np.float32).astype(float)
    return s


def save_py(path, s):
    np.savez_compressed(path, pos=s.pos.astype(np.float32), vel=s.vel.astype(np.float32),
                        typ=s.typ, spr=s.spray, esc=s.escaped, rad=s.radius.astype(np.float32))


def run_python(out, steps, dt, every, model, phase_law="mass"):
    out.mkdir(parents=True)
    s = make_sim(model, phase_law)
    records = []
    start = time.perf_counter()
    for j in range(steps):
        s.step(dt)
        # The C++ graph exports float32 particle state between its iterations.
        for name in ("pos", "vel", "radius", "scale", "tau"):
            setattr(s, name, getattr(s, name).astype(np.float32).astype(float))
        assert np.isfinite(s.pos).all() and np.isfinite(s.vel).all()
        speed = max(float(np.abs(s.u).max()), float(np.abs(s.v).max()))
        rec = dict(s.emission_stats if model==2 else dict(accepted=s.n_emitted, alive=int(s.spray.sum())))
        rec.update(step=j, time=(j+1)*dt, cfl=speed*dt, iterations=s.iters,
                   primary=int((~s.spray).sum()), escaped=int((s.escaped & (s.typ==1) & ~s.spray).sum()))
        records.append(rec)
        if j % every == 0:
            save_py(out/f"frame.{j:04d}.npz", s)
        if (j+1) % 50 == 0 or j==steps-1:
            print(f"Python {j+1}/{steps}, t={(j+1)*dt:.3f}, spray={int(s.spray.sum())}, "
                  f"CFL={speed*dt:.3f}, elapsed={time.perf_counter()-start:.1f}s", flush=True)
    save_py(out/"final.npz", s)
    (out/"ledger.json").write_text(json.dumps(records), encoding="utf-8")
    return dict(seconds=time.perf_counter()-start, cfl_max=max(r["cfl"] for r in records),
                emitted=sum(r["accepted"] for r in records), alive=int(s.spray.sum()))


def run_cpp(out, steps, dt, every, model, phase_law="mass"):
    if out.exists():
        raise FileExistsError("Preserve existing arm: " + str(out))
    s = make_sim(model)
    p3 = np.c_[s.pos, np.full(len(s.pos), .5)].astype(np.float32)
    start = time.perf_counter()
    os.environ["PFFLIP_SPRAY_AUDIT"] = "1"
    _, _, diags = run_bifrost3("spray_events_2d", NX, NY, 1, lambda x,y,z: mask(x,y),
        steps, dt, seed_fn=lambda *args: (p3, s.typ.copy()), fused=True,
        extra_ports=dict(ESCAPE, **ports(dict(E, model=model), D), rho0_face=s.rho0_face),
        dump_every=every, final_state=True, dumps=("esc", "ph", "rad", "spr", "vel"),
        out_dir=str(out), log_prefix=str(out/"bifcmd"))
    assert np.isfinite(diags).all()
    records = [json.loads(line[len("SPRAY_EVENT "):]) for line in
               (out/"bifcmd.stderr.log").read_text().splitlines() if line.startswith("SPRAY_EVENT ")]
    (out/"ledger.json").write_text(json.dumps(records), encoding="utf-8")
    return dict(seconds=time.perf_counter()-start, cfl_max=float(diags[:,3].max()*dt),
                emitted=sum(r["accepted"] for r in records), alive=records[-1]["alive"] if records else None)


def load_frame(out, arm, step):
    if arm=="python":
        with np.load(out/arm/f"frame.{step:04d}.npz") as a:
            return {k:a[k] for k in a.files}
    job = json.loads((out/"job.json").read_text())
    d = Path(job["comparison_cpp_source"])/arm if job.get("comparison_cpp_source") else out/arm
    return dict(pos=np.load(d/f"frame.{step:04d}.npy")[:,:2],
                vel=np.load(d/f"vel.{step:04d}.npy")[:,:2],
                typ=np.load(d/f"ph.{step:04d}.npy")>.5,
                spr=np.load(d/f"spr.{step:04d}.npy")>.5,
                esc=np.load(d/f"esc.{step:04d}.npy")>.5,
                rad=np.load(d/f"rad.{step:04d}.npy"))


def metrics(out, steps, dt, every):
    rows = []
    for j in range(0, steps, every):
        pair = [load_frame(out, arm, j) for arm in ("python", "bifrost")]
        waters = [a["pos"][(a["typ"]==1)&~a["spr"]] for a in pair]
        hist = [np.histogram2d(p[:,0],p[:,1],bins=(16,20),range=((0,NX),(0,NY)))[0]
                for p in waters]
        tv = float(.5*np.abs(hist[0]/hist[0].sum()-hist[1]/hist[1].sum()).sum())
        com = float(np.linalg.norm(waters[0].mean(0)-waters[1].mean(0)))
        row = dict(time=(j+1)*dt, water_histogram_tv=tv, water_com_error_cells=com,
                   water_counts=[len(p) for p in waters], spray_counts=[int(a["spr"].sum()) for a in pair],
                   water_height_p99=[float(np.quantile(p[:,1],.99)) for p in waters])
        rows.append(row)
    result = dict(rows=rows, max_water_tv=max(r["water_histogram_tv"] for r in rows),
                  max_com_error_cells=max(r["water_com_error_cells"] for r in rows),
                  max_p99_height_error=max(abs(r["water_height_p99"][0]-r["water_height_p99"][1]) for r in rows))
    result["bulk_gate_pass"] = (result["max_water_tv"] < .15 and
                                result["max_com_error_cells"] < 1 and result["max_p99_height_error"] < 4)
    (out/"comparison.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k!="rows"}), flush=True)
    return result


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--seconds", type=float, default=20)
    p.add_argument("--dt", type=float, default=.01)
    p.add_argument("--every", type=int, default=10)
    p.add_argument("--model", type=int, choices=(1,2), default=2)
    p.add_argument("--phase-law", choices=("mass","tag"), default="mass")
    p.add_argument("--bifrost-from", type=Path)
    p.add_argument("--phase", choices=("all","python","bifrost","compare"), default="all")
    a = p.parse_args(); out = a.out.resolve()
    dt, steps = float(np.float32(a.dt)), round(a.seconds/a.dt)
    dirty = subprocess.check_output(["git","status","--porcelain"],cwd=ROOT,text=True).strip()
    if a.phase != "compare" and dirty:
        raise RuntimeError("Commit intended experiment sources first: " + dirty)
    out.mkdir(parents=True, exist_ok=True)
    job_path = out/"job.json"
    if not job_path.exists():
        sha = subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
        dirty = subprocess.check_output(["git","status","--porcelain"],cwd=ROOT,text=True).strip()
        dll = ROOT/"nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll"
        if not dll.exists():
            dll = ROOT/"nodes/build/src/PFFlipSpikeOps.dll"
        job = dict(kind="2D spray validation", owner="codex", created="2026-10-23", commit=sha,
                   worktree=str(ROOT), dirty=dirty, reads="generated 2D seed only", nx=NX, ny=NY,
                   dt=dt, steps=steps, dump_every=a.every, droplet=D, emission=dict(E,model=a.model),
                   phase_law=a.phase_law, comparison_cpp_source=str(a.bifrost_from.resolve()) if a.bifrost_from else None,
                   binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(), command=sys.argv,
                   retention="small reproducible validation; existing 40s data untouched")
        job_path.write_text(json.dumps(job,indent=2), encoding="utf-8")
    receipt_path = out/"runtime.json"
    receipt = json.loads(receipt_path.read_text()) if receipt_path.exists() else {}
    for name, fn in (("bifrost",run_cpp),("python",run_python)):
        if a.phase in ("all",name):
            receipt[name] = fn(out/name,steps,dt,a.every,a.model,a.phase_law)
            receipt[name]["commit"] = subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip()
            receipt[name]["reference_sha256"] = hashlib.sha256((ROOT/"reference2d/pfflip2d.py").read_bytes()).hexdigest()
            receipt_path.write_text(json.dumps(receipt,indent=2))
    if a.phase in ("all","compare"):
        metrics(out,steps,dt,a.every)


if __name__=="__main__":
    main()
