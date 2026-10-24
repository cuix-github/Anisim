"""Finite 2D base-layer gate; no spray, no production, no rendering.

Control: reproduce the previous volume-only graph. One-variable arms then
compare transport-only to volume+transport. A failure stops this gate and is
retained; downstream long-impact acceptance is not launched by this script.
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

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference2d'),str(ROOT/'gates'),str(ROOT/'graphs')]
from validate_volume_contract import make,measure
from validate_transport_contract import arrays
from phase_area import area,particle_field,triangle_fraction
from loop3d import run_bifrost3,load_final_state

DT=float(np.float32(.0025))
LIMITS=dict(local_position=8e-6,local_velocity=5e-5,rest_rms=.1,
            rest_com=.1,energy_ratio=1.005,area_drift=.02,
            analytic_area_error=1e-10,seeded_geometry_relative_error=.05)


def write(p,x):p.write_text(json.dumps(x,indent=2),encoding='utf-8')


def geometry(s):return area(particle_field(s))


def calibration():
    assert float(triangle_fraction(np.array(-1.),np.array(-1.),np.array(1.),0))==.25
    assert float(triangle_fraction(np.array(-1.),np.array(1.),np.array(1.),0))==.75
    rows=[];nx,ny=48,40
    for h in (7.3,12.,20.6):
        field=np.broadcast_to(np.arange(ny)+.5-h,(nx,ny))
        a=area(field,0)['area'];expected=nx*(ny-h)
        assert abs(a-expected)<LIMITS['analytic_area_error']
        rows.append(dict(kind='analytic_plane',height=h,measured=a,expected=expected))
    s=make(0,'rest')
    for name,mask,expected in [
        ('plane12',lambda x,y:y<12,48*12),
        ('plane20',lambda x,y:y<20,48*20),
        ('disk8',lambda x,y:(x-24)**2+(y-20)**2<8**2,np.pi*8**2),
        ('disk12',lambda x,y:(x-24)**2+(y-20)**2<12**2,np.pi*12**2)]:
        s.seed(mask,jitter=0,seed=1)
        s.pos=np.maximum(.51,np.minimum([nx-.51,ny-.51],s.pos)).astype(np.float32).astype(float)
        a=geometry(s);err=abs(a['area']/expected-1)
        rows.append(dict(kind=name,**a,expected=expected,relative_error=err))
    return dict(passed=all(r.get('relative_error',0)<=LIMITS['seeded_geometry_relative_error'] for r in rows),rows=rows)


def cpp(s,out,steps,state=None,offset=0):
    assert not out.exists()
    p=dict(gravity=s.g,phase_model=2,phase_alpha=1.,rho0_face=s.rho0_face,
        tolerance=s.pressure_tolerance,alpha_liquid=s.alpha[1],alpha_air=s.alpha[0],
        volume_tau=s.volume_tau,transport_model=s.transport_model,
        escape=0,emit_model=0,droplet_model=0,st=0,adapt=0)
    _,_,diag=run_bifrost3('base_combination',s.nx,s.ny,1,lambda x,y,z:y<12,
        steps,DT,init_arrays=arrays(s) if state is None else state,step_offset=offset,
        fused=True,final_state=True,dump_every=40,dumps=('ph','vel'),extra_ports=p,
        out_dir=str(out),log_prefix=str(out/'bifcmd'))
    assert np.isfinite(diag).all() and np.max(diag[:,1])<=s.pressure_tolerance
    return float(diag[:,1].max())


def pair(out,tau,mode,case,steps):
    s=make(tau,case);s.transport_model=mode
    initial=dict(measure(s),**geometry(s));started=time.perf_counter()
    residual=cpp(s,out/'bifrost',steps)
    histories={}
    for backend in ('python','bifrost'):
        q=make(tau,case);q.transport_model=mode;h=[dict(time=0,**initial)]
        for j in range(steps):
            if backend=='python':
                q.step(DT);q.pos=q.pos.astype(np.float32).astype(float);q.vel=q.vel.astype(np.float32).astype(float)
            elif j%40==0 or j==steps-1:
                end=j==steps-1;d=out/'bifrost'
                q.pos=np.load(d/('final.npy' if end else f'frame.{j:04d}.npy'))[:,:2].astype(float)
                q.vel=np.load(d/('final_velocities.npy' if end else f'vel.{j:04d}.npy'))[:,:2].astype(float)
            if j%40==0 or j==steps-1:h.append(dict(time=(j+1)*DT,**measure(q),**geometry(q)))
        histories[backend]=h
        if backend=='python':pp=q.pos.copy();pv=q.vel.copy()
    row=dict(tau=tau,transport=mode,case=case,steps=steps,initial=initial,**histories,
        position_max=float(abs(q.pos-pp).max()),velocity_max=float(abs(q.vel-pv).max()),
        pressure_residual_max=residual,seconds=time.perf_counter()-started)
    write(out/'receipt.json',row);return row


def identical(a,b):
    names=['final.npy']+[p.name for p in a.glob('final_*.npy')]
    checks={n:bool(np.array_equal(np.load(a/n),np.load(b/n))) for n in names}
    assert len(checks)==13 and all(checks.values()),checks
    return checks


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);os.environ['OMP_NUM_THREADS']='1'
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    write(out/'job.json',dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),thresholds=LIMITS,dt=DT,
        scope=__doc__,geometry='cell-centre phi=.5 piecewise-linear area; mean of both diagonals; boundary constant extension',
        ports='phase=2; escape/emission/droplet/ST/adapt off; tau 0/.5; transport 0/3'))
    r=dict(passed=False,cases={})
    try:
        r['calibration']=calibration();write(out/'receipt.json',r)
        # Existing known state first, after a guard-only C++ change.
        s=make(.5,'density');cpp(s,out/'control',1)
        r['control']=identical(ROOT/'out/runs/volume_candidate_20261024/local/bifrost',out/'control')
        for case in ('density','rest'):
            row=pair(out/('local_'+case),.5,3,case,3);r['cases']['local_'+case]=row
            assert row['position_max']<LIMITS['local_position'] and row['velocity_max']<LIMITS['local_velocity']
        # Exact checkpoint continuity of the combined state, including momentum.
        s=make(.5,'density');s.transport_model=3
        cpp(s,out/'resume/full',80);cpp(s,out/'resume/first',40)
        cpp(s,out/'resume/second',40,load_final_state(str(out/'resume/first')),40)
        r['resume']=identical(out/'resume/full',out/'resume/second');write(out/'receipt.json',r)
        for name,tau in [('transport_only',0.),('combined',.5)]:
            row=pair(out/name,tau,3,'rest',2400);r['cases'][name]=row
            checks={}
            for backend in ('python','bifrost'):
                h=row[backend];f=h[-1];first=h[0]
                checks[backend]=dict(rms=f['rms']<=LIMITS['rest_rms'],
                    com=np.linalg.norm(np.array(f['com'])-first['com'])<=LIMITS['rest_com'],
                    energy=max(q['energy'] for q in h)/first['energy']<=LIMITS['energy_ratio'],
                    area=max(abs(q['area']/first['area']-1) for q in h)<=LIMITS['area_drift'])
                checks[backend]={k:bool(v) for k,v in checks[backend].items()}
            row['checks']=checks;write(out/'receipt.json',r)
            assert all(all(v.values()) for v in checks.values()),(name,checks)
        # A failed geometric calibration does not invalidate independent code
        # and static-motion checks, but still blocks the overall volume gate.
        r['mechanics_passed']=True
        assert r['calibration']['passed'], 'geometry calibration failed; mechanical subgates do not certify volume'
        r['passed']=True
    except BaseException as e:
        r['failure']=repr(e);raise
    finally:
        write(out/'receipt.json',r);print(json.dumps(dict(passed=r['passed'],failure=r.get('failure')),indent=2),flush=True)


if __name__=='__main__':main()
