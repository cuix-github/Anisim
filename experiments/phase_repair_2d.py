"""Explicit matched-phase, independent-solver 20s experiment and ablation ladder.

Phase 2 is Eq.7 plus pure-cell faces, with the documented cell-centre sampling
choice. Legacy defaults are never silently selected. No production inputs.
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
sys.path[:0] = [str(ROOT/'experiments'), str(ROOT/'reference2d'), str(ROOT/'graphs')]
from spray_events_2d import NX, NY, D, E, ESCAPE, mask, ports, save_py
from phase_contract import PhaseContractSim
from loop3d import run_bifrost3


def make(layer, model=2, tol=1e-6, dt=.0025, escape_rule=1):
    escape = escape_rule if layer != 'base' else 0
    # Preserve the original FLIP weights at the original dt; when subdividing,
    # preserve the Eq.13 viscosity proxy, rather than silently doubling damping.
    ratio = dt/.0025
    s = PhaseContractSim(NX, NY, phase_model=model, phase_alpha=1.,
        pressure_tolerance=float(np.float32(tol)), g=float(np.float32(-9.8)),
        alpha_l=float(np.float32(1-.03*ratio)), alpha_g=float(np.float32(1-.1*ratio)),
        escape=escape, esc_phi=float(np.float32(.3)),
        droplet=D if escape else None, emit=E if layer=='full' else None,
        spray_cfl=.5)
    s.seed(mask, seed=1); s.calibrate()
    s.rho0_face=float(np.float32(s.rho0_face))
    s.pos=s.pos.astype(np.float32).astype(float)
    s.escaped=np.zeros(len(s.pos), bool)
    return s


def py(out, s, steps, dt, every):
    out.mkdir(); ledger=[]; start=time.perf_counter()
    for j in range(steps):
        s.step(dt)
        for n in ('pos','vel','radius','scale','tau'):
            setattr(s,n,getattr(s,n).astype(np.float32).astype(float))
        assert np.isfinite(s.pos).all() and np.isfinite(s.vel).all()
        record=dict(s.emission_stats) if s.emit else {}
        record.update(step=j,time=(j+1)*dt,iterations=int(s.iters),
            cfl=float(max(abs(s.u).max(),abs(s.v).max())*dt),
            primary=int((~s.spray).sum()),alive=int(s.spray.sum()),
            escaped_water=int((s.escaped & (s.typ==1) & ~s.spray).sum()))
        ledger.append(record)
        if j%every==0: save_py(out/f'frame.{j:04d}.npz',s)
        if (j+1)%400==0:
            print(f'Python {j+1}/{steps}, spray={s.spray.sum()}, elapsed={time.perf_counter()-start:.1f}s',flush=True)
    save_py(out/'final.npz',s)
    (out/'ledger.json').write_text(json.dumps(ledger))
    return dict(seconds=time.perf_counter()-start,steps=len(ledger),max_cfl=max(r['cfl'] for r in ledger))


def cpp(out, s, steps, dt, every):
    out=out.resolve()  # bifcmd NPY IO requires absolute paths, also for direct API use
    if out.exists(): raise FileExistsError(out)
    config=dict(ESCAPE,**ports(E,D),phase_model=s.phase_model,phase_alpha=s.phase_alpha,
        rho0_face=s.rho0_face,tolerance=s.pressure_tolerance,transport_model=0,
        alpha_liquid=s.alpha[1],alpha_air=s.alpha[0],st=0,adapt=0,pressure_solver=0)
    config.update(escape=int(s.escape),droplet_model=2 if s.escape else 0,
                  emit_model=2 if s.emit else 0)
    p3=np.c_[s.pos,np.full(len(s.pos),.5)].astype(np.float32)
    os.environ['PFFLIP_SPRAY_AUDIT']='1'
    start=time.perf_counter()
    _,_,diags=run_bifrost3('phase_repair',NX,NY,1,lambda x,y,z:mask(x,y),steps,dt,
        seed_fn=lambda *a:(p3,s.typ.copy()),fused=True,extra_ports=config,
        final_state=True,dumps=('ph','vel','rad','spr','esc'),dump_every=every,
        out_dir=str(out),log_prefix=str(out/'bifcmd'))
    assert np.isfinite(diags).all()
    events=[json.loads(l[len('SPRAY_EVENT '):]) for l in (out/'bifcmd.stderr.log').read_text().splitlines()
            if l.startswith('SPRAY_EVENT ')]
    (out/'ledger.json').write_text(json.dumps(events))
    return dict(seconds=time.perf_counter()-start,steps=len(diags),max_cfl=float(diags[:,3].max()*dt),
                max_pressure_residual=float(diags[:,1].max()),max_iterations=float(diags[:,0].max()))


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path)
    p.add_argument('--layer',choices=('base','escape','full'),required=True)
    p.add_argument('--model',type=int,choices=(1,2),default=2)
    p.add_argument('--seconds',type=float,default=20);p.add_argument('--dt',type=float,default=.0025)
    p.add_argument('--tolerance',type=float,default=1e-6);p.add_argument('--every',type=int,default=40)
    p.add_argument('--escape-rule',type=int,choices=(1,2),default=1)
    a=p.parse_args()
    dirty=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    if dirty: raise RuntimeError('Commit experiment first: '+dirty)
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    dt=float(np.float32(a.dt));steps=round(a.seconds/a.dt)
    s=make(a.layer,a.model,a.tolerance,a.dt,a.escape_rule)
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    job=dict(kind='Independent uniform phase repair experiment',commit=subprocess.check_output(
        ['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),command=sys.argv,
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),layer=a.layer,
        phase_model=a.model,phase_alpha=1,cell_phase_sampling='separate Eq.6 cell-centre mass splat',
        dt=dt,steps=steps,dump_every=a.every,nx=NX,ny=NY,rho0_face=s.rho0_face,
        tolerance=s.pressure_tolerance,alpha_liquid=s.alpha[1],alpha_air=s.alpha[0],
        droplet=D if s.escape else None,emission=E if s.emit else None,escape=s.escape,
        esc_phi=s.esc_phi,transport_model=0,independent_momentum_port=True,
        thresholds=dict(tv=.15,com_cells=1,height_p99_cells=4),
        omp_num_threads=os.environ.get('OMP_NUM_THREADS','runtime default'),
        scope='2D uniform; no ST/adapt/solids; escape=1 count law, escape=2 phase threshold only')
    (out/'job.json').write_text(json.dumps(job,indent=2))
    np.savez_compressed(out/'initial.npz',pos=s.pos,vel=s.vel,typ=s.typ,spr=s.spray,
                        esc=s.escaped,rad=s.radius)
    receipt={}
    for name,fn in [('bifrost',cpp),('python',py)]:
        receipt[name]=fn(out/name,make(a.layer,a.model,a.tolerance,a.dt,a.escape_rule),steps,dt,a.every)
        (out/'runtime.json').write_text(json.dumps(receipt,indent=2))
    print(json.dumps(receipt,indent=2),flush=True)


if __name__=='__main__': main()
