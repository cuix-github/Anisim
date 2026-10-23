"""Serialized transport ablations; no C++/production modification or rendering."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference2d'),str(ROOT/'experiments')]
from transport_contract import TransportContractSim,extend_faces,rk3_local
from phase_contract_2d import make,measures


def analytic():
    values=np.full((12,11),99.);valid=np.zeros(values.shape,bool);valid[3:7,3:7]=True
    values[valid]=2.75
    ext,known=extend_faces(values,valid,3)
    assert np.array_equal(ext[valid],values[valid]) and np.all(ext[known]==2.75)
    assert np.all(ext[~known]==99.)
    pos=np.array([[2.,3.],[-1.,.5]])
    vel=np.array([2.,-.3])
    moved=rk3_local(pos,2.,lambda p,ids:np.broadcast_to(vel,p.shape),cfl=.25)
    assert np.allclose(moved,pos+2*vel,rtol=0,atol=1e-13)
    # Linear field rotates a point: trilinear grid sampling can represent it exactly.
    errors=[]
    for dt in (.4,.2,.1):
        p=np.array([[1.,0.]])
        for _ in range(round(2/dt)):
            p=rk3_local(p,dt,lambda q,ids:np.c_[-q[:,1],q[:,0]])
        errors.append(float(np.linalg.norm(p[0]-[np.cos(2),np.sin(2)])))
    assert errors[0]/errors[1]>6 and errors[1]/errors[2]>6,errors
    return dict(rotation_errors=errors,constant_extension=True,constant_transport=True)


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--mode',type=int,required=True);p.add_argument('--seconds',type=float,default=6)
    p.add_argument('--case',choices=('plunge','rest','wide'),default='plunge')
    p.add_argument('--control',type=Path);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    tol=float(np.float32(1e-10));dt=float(np.float32(.0025));steps=round(a.seconds/.0025)
    s=make('rest' if a.case=='rest' else 'plunge',2,tol)
    s.__class__=TransportContractSim;s.transport_model=a.mode
    if a.case=='wide':
        s=TransportContractSim(112,80,phase_model=2,transport_model=a.mode,pressure_tolerance=tol,
            g=float(np.float32(-9.8)),alpha_l=float(np.float32(.97)),alpha_g=float(np.float32(.9)))
        s.seed(lambda x,y:(y<18)|((abs(x-56)<7)&(y>39)&(y<64)),seed=1)
        s.calibrate();s.rho0_face=float(np.float32(s.rho0_face));s.pos=s.pos.astype(np.float32).astype(float)
    initial=measures(s.pos,s.vel,s.typ);rows=[];restarts=[];start=time.perf_counter()
    job=dict(kind='Transport ablation',mode=a.mode,case=a.case,steps=steps,dt=dt,initial=initial,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),analytic=analytic())
    (out/'job.json').write_text(json.dumps(job,indent=2))
    for j in range(steps):
        s.step(dt)
        if getattr(s,'pressure_restarts',0):
            restarts.append(dict(step=j,count=s.pressure_restarts,true_residual=s.pressure_residual))
        for n in ('pos','vel','radius','scale','tau'):setattr(s,n,getattr(s,n).astype(np.float32).astype(float))
        assert np.isfinite(s.pos).all() and np.isfinite(s.vel).all()
        if j%40==0:
            np.savez_compressed(out/f'frame.{j:04d}.npz',pos=s.pos.astype(np.float32),vel=s.vel.astype(np.float32),typ=s.typ)
            rows.append(dict(step=j,time=(j+1)*dt,**measures(s.pos,s.vel,s.typ)))
            if a.control:
                old=np.load(a.control/f'frame.{j:04d}.npz')
                for n in ('pos','vel','typ'):assert np.array_equal(old[n],getattr(s,{'pos':'pos','vel':'vel','typ':'typ'}[n])),(j,n)
        if (j+1)%400==0:print(f'mode={a.mode} {a.case} t={(j+1)*dt:.2f} elapsed={time.perf_counter()-start:.1f}s',flush=True)
    np.savez_compressed(out/'final.npz',pos=s.pos,vel=s.vel,typ=s.typ)
    final=measures(s.pos,s.vel,s.typ)
    if a.control:
        old=np.load(a.control/'final.npz')
        for n in ('pos','vel','typ'):assert np.array_equal(old[n],getattr(s,n)),n
    (out/'metrics.json').write_text(json.dumps(rows,indent=2))
    result=dict(wall_seconds=time.perf_counter()-start,final=final,
        max_sampled_energy_ratio=max([initial['energy'],final['energy']]+[r['energy'] for r in rows])/initial['energy'],
        bit_identical_control=bool(a.control),pressure_restarts=restarts)
    (out/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)


if __name__=='__main__':main()
