"""Small serialized independent-fluid diagnostics; no render or production IO."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference2d'),str(ROOT/'graphs')]
from phase_contract import PhaseContractSim
from loop3d import run_bifrost3


def make(case,model,tol):
    nx,ny=(48,40) if case=='rest' else (64,80)
    s=PhaseContractSim(nx,ny,phase_model=model,g=float(np.float32(-9.8)),
        alpha_l=float(np.float32(.97)),alpha_g=float(np.float32(.9)),pressure_tolerance=tol)
    mask=(lambda x,y:y<12) if case=='rest' else (lambda x,y:(y<18)|((abs(x-32)<7)&(y>39)&(y<64)))
    s.seed(mask,seed=1);s.calibrate();s.rho0_face=float(np.float32(s.rho0_face))
    s.pos=s.pos.astype(np.float32).astype(float)
    return s


def measures(pos,vel,typ):
    v=vel.astype(float);p=pos.astype(float);m=np.where(typ==1,1000.,1.)
    result=dict(energy=float(np.sum(m*(.5*np.sum(v*v,axis=1)+float(np.float32(9.8))*p[:,1]))),
                count=len(p),max_component=float(abs(v).max()))
    for name,sel in [('water',typ==1),('air',typ==0)]:
        result[name+'_rms']=float(np.sqrt(np.mean(np.sum(v[sel]**2,axis=1))))
    result['water_com']=p[typ==1].mean(0).tolist()
    result['water_p99']=float(np.quantile(p[typ==1,1],.99))
    return result


def run_python(out,s,steps,dt,every):
    out.mkdir();records=[];start=time.perf_counter()
    for j in range(steps):
        s.step(dt)
        for n in ('pos','vel','radius','scale','tau'):
            setattr(s,n,getattr(s,n).astype(np.float32).astype(float))
        assert np.isfinite(s.pos).all() and np.isfinite(s.vel).all()
        if j%every==0:
            np.savez_compressed(out/f'frame.{j:04d}.npz',pos=s.pos.astype(np.float32),vel=s.vel.astype(np.float32),typ=s.typ)
            records.append(dict(step=j,time=(j+1)*dt,**measures(s.pos,s.vel,s.typ)))
        if (j+1)%400==0:
            print(f'Python t={(j+1)*dt:.3f}s elapsed={time.perf_counter()-start:.1f}s',flush=True)
    np.savez_compressed(out/'final.npz',pos=s.pos,vel=s.vel,typ=s.typ)
    (out/'metrics.json').write_text(json.dumps(records,indent=2))
    return dict(wall_seconds=time.perf_counter()-start,final=measures(s.pos,s.vel,s.typ))


def run_cpp(out,s,steps,dt,every):
    start=time.perf_counter();p3=np.c_[s.pos,np.full(len(s.pos),.5)].astype(np.float32)
    run_bifrost3('phase_diagnostic',s.nx,s.ny,1,lambda x,y,z:y<1,steps,dt,
        seed_fn=lambda *a:(p3,s.typ.copy()),fused=True,final_state=True,
        extra_ports=dict(phase_model=s.phase_model,phase_alpha=s.phase_alpha,rho0_face=s.rho0_face,
            escape=0,st=0,adapt=0,emit_model=0,tolerance=s.pressure_tolerance),
        dumps=('ph','vel'),dump_every=every,out_dir=str(out),log_prefix=str(out/'bifcmd'))
    records=[]
    for j in range(0,steps,every):
        records.append(dict(step=j,time=(j+1)*dt,**measures(np.load(out/f'frame.{j:04d}.npy')[:,:2],
            np.load(out/f'vel.{j:04d}.npy')[:,:2],np.load(out/f'ph.{j:04d}.npy')>.5)))
    (out/'metrics.json').write_text(json.dumps(records,indent=2))
    return dict(wall_seconds=time.perf_counter()-start,final=measures(np.load(out/'final.npy')[:,:2],
        np.load(out/'final_velocities.npy')[:,:2],np.load(out/'final_phase.npy')>.5))


def compare(out,s,steps,every,initial,case):
    stats={arm:json.loads((out/arm/'metrics.json').read_text()) for arm in ('python','bifrost')}
    rows=[]
    for j in range(0,steps,every):
        with np.load(out/'python'/f'frame.{j:04d}.npz') as a:
            py=a['pos'][a['typ']==1]
        d=out/'bifrost';cpp=np.load(d/f'frame.{j:04d}.npy')[:,:2][np.load(d/f'ph.{j:04d}.npy')>.5]
        hs=[np.histogram2d(p[:,0],p[:,1],bins=(s.nx//4,s.ny//4),range=((0,s.nx),(0,s.ny)))[0] for p in (py,cpp)]
        rows.append(dict(step=j,tv=float(.5*abs(hs[0]/hs[0].sum()-hs[1]/hs[1].sum()).sum()),
            com=float(np.linalg.norm(py.mean(0)-cpp.mean(0))),
            height_p99=float(abs(np.quantile(py[:,1],.99)-np.quantile(cpp[:,1],.99)))))
    # Sparse dumps precede the actual final step. Include the saved final state
    # in both physics and comparison checks instead of calling the last dump final.
    with np.load(out/'python'/'final.npz') as a:
        py=a['pos'][a['typ']==1]
    d=out/'bifrost';cpp=np.load(d/'final.npy')[:,:2][np.load(d/'final_phase.npy')>.5]
    hs=[np.histogram2d(p[:,0],p[:,1],bins=(s.nx//4,s.ny//4),range=((0,s.nx),(0,s.ny)))[0] for p in (py,cpp)]
    rows.append(dict(step=steps-1,tv=float(.5*abs(hs[0]/hs[0].sum()-hs[1]/hs[1].sum()).sum()),
        com=float(np.linalg.norm(py.mean(0)-cpp.mean(0))),
        height_p99=float(abs(np.quantile(py[:,1],.99)-np.quantile(cpp[:,1],.99)))))
    runtime=json.loads((out/'runtime.json').read_text())
    worst={k:max(r[k] for r in rows) for k in ('tv','com','height_p99')}
    result=dict(max_difference=worst,similarity_pass=worst['tv']<.05 and worst['com']<.2 and worst['height_p99']<1,
                states=rows,physics={})
    for arm,rec in stats.items():
        final=runtime[arm]['final']
        energy=max([r['energy'] for r in rec]+[final['energy']])/initial['energy']
        drift=abs(final['water_com'][1]-initial['water_com'][1])
        result['physics'][arm]=dict(max_energy_ratio=energy,final_water_rms=final['water_rms'],
            final_air_rms=final['air_rms'],water_com_y_drift=drift,
            rest_pass=(final['water_rms']<.1 and drift<.1 and energy<1.005) if case=='rest' else None)
    (out/'comparison.json').write_text(json.dumps(result,indent=2))
    return {k:v for k,v in result.items() if k!='states'}


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--case',choices=('rest','plunge'),required=True);p.add_argument('--model',type=int,choices=(1,2),required=True)
    p.add_argument('--seconds',type=float,default=6);p.add_argument('--dt',type=float,default=.0025)
    p.add_argument('--tolerance',type=float,default=1e-6);p.add_argument('--every',type=int,default=40)
    a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    dirty=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip();assert not dirty,dirty
    dt,tol=float(np.float32(a.dt)),float(np.float32(a.tolerance));steps=round(a.seconds/a.dt)
    s=make(a.case,a.model,tol);initial=measures(s.pos,s.vel,s.typ)
    job=dict(kind='Independent 2D phase comparison',case=a.case,phase_model=a.model,dt=dt,tolerance=tol,
        steps=steps,dump_every=a.every,nx=s.nx,ny=s.ny,initial=initial,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256((ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll').read_bytes()).hexdigest(),
        physics='escape/ST/adapt/emission off; closed domain',command=sys.argv,
        independent_momentum_port=True,includes_actual_final_state=True,
        thresholds=dict(tv=.05,com_cells=.2,height_p99_cells=1,rest_water_rms=.1,rest_com_drift=.1,rest_energy_ratio=1.005))
    (out/'job.json').write_text(json.dumps(job,indent=2))
    runtime={}
    for arm,fn in [('bifrost',run_cpp),('python',run_python)]:
        runtime[arm]=fn(out/arm,make(a.case,a.model,tol),steps,dt,a.every)
        (out/'runtime.json').write_text(json.dumps(runtime,indent=2))
    result=compare(out,s,steps,a.every,initial,a.case);print(json.dumps(result,indent=2),flush=True)
    if not result['similarity_pass'] or (a.case=='rest' and not all(r['rest_pass'] for r in result['physics'].values())):
        raise SystemExit('Validation not accepted; failed arms retained')


if __name__=='__main__':main()
