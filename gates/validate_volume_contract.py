"""Local graph lockstep, density relaxation and fixed still-water acceptance.

Runs stop at the first failed physical gate. No escape, spray, 3D or rendering.
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
from volume_contract import VolumeContractSim,target_divergence
from audit_volume_contract import make as witness_make
from loop3d import run_bifrost3


def make(tau,case,model=0):
    nx,ny=(48,40) if case=='rest' else (24,24)
    s=VolumeContractSim(nx,ny,volume_tau=tau,volume_model=model,g=float(np.float32(-9.8)) if case=='rest' else 0,
        pressure_tolerance=float(np.float32(1e-6)),alpha_l=float(np.float32(.97)),alpha_g=float(np.float32(.9)))
    s.seed(lambda x,y:y<12,jitter=.35 if case=='rest' else 0,seed=1);s.calibrate()
    s.pos=np.maximum(.51,np.minimum([nx-.51,ny-.51],s.pos)).astype(np.float32).astype(float)
    if case=='density':s.pos=witness_make(True).pos.copy()
    s.rho0_face=float(np.float32(s.rho0_face));return s


def measure(s):
    liquid=s.typ==1
    m=s.masses()
    e=np.sum(m*(.5*np.sum(s.vel**2,axis=1)-s.g*s.pos[:,1]))
    log=np.log(s.number_density()/s.density_reference()); log-=log.mean()
    return dict(energy=float(e),rms=float(np.sqrt(np.mean(np.sum(s.vel[liquid]**2,axis=1)))),
        com=s.pos[liquid].mean(axis=0).tolist(),density_log_rms=float(np.sqrt(np.mean(log**2))))


def run(out,tau,case,steps,model=0):
    s=make(tau,case,model);initial=measure(s);dt=float(np.float32(.0025))
    p3=np.c_[s.pos,np.full(len(s.pos),.5)].astype(np.float32)
    d=out/'bifrost';t=time.perf_counter()
    _,_,diag=run_bifrost3('volume_gate',s.nx,s.ny,1,lambda x,y,z:y<12,steps,dt,
        seed_fn=lambda *a:(p3,s.typ.copy()),fused=True,final_state=True,dump_every=40,
        dumps=('ph','vel'),out_dir=str(d),log_prefix=str(d/'bifcmd'),
        extra_ports=dict(gravity=s.g,phase_model=2,phase_alpha=1.,rho0_face=s.rho0_face,
            escape=0,emit_model=0,droplet_model=0,transport_model=0,volume_tau=tau,volume_model=model,st=0,adapt=0))
    assert np.isfinite(diag).all() and np.max(diag[:,1])<=float(np.float32(1e-6))
    history=[initial]
    for j in range(steps):
        s.step(dt)
        s.pos=s.pos.astype(np.float32).astype(float);s.vel=s.vel.astype(np.float32).astype(float)
        if j%40==0:history.append(measure(s))
    history.append(measure(s))
    b=make(tau,case,model);chistory=[initial]
    for j in range(0,steps,40):
        b.pos=np.load(d/f'frame.{j:04d}.npy')[:,:2].astype(float)
        b.vel=np.load(d/f'vel.{j:04d}.npy')[:,:2].astype(float)
        chistory.append(measure(b))
    b.pos=np.load(d/'final.npy')[:,:2].astype(float);b.vel=np.load(d/'final_velocities.npy')[:,:2].astype(float)
    chistory.append(measure(b))
    row=dict(steps=steps,tau=tau,volume_model=model,initial=initial,python=history,bifrost=chistory,
        final_position_max=float(np.max(abs(b.pos-s.pos))),final_velocity_max=float(np.max(abs(b.vel-s.vel))),
        max_pressure_residual=float(diag[:,1].max()),seconds=time.perf_counter()-t)
    (out/'receipt.json').write_text(json.dumps(row,indent=2));return row


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);os.environ['OMP_NUM_THREADS']='1'
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    rec=dict(passed=False,cases={})
    (out/'job.json').write_text(json.dumps(dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),command=sys.argv,
        scope=__doc__,thresholds=dict(single_step_position=4e-6,single_step_velocity=2e-5,
            density_rms_ratio_after_one_second=.8,rest_rms=.1,rest_com=.1,rest_energy_ratio=1.005)),indent=2))
    try:
        # Independent formula truth table and invalid support. No clamp of zero density.
        d=target_divergence(np.array([.8,1,1.2]),np.ones(3),.5)
        assert d[0]<0<d[-1] and abs(d.sum())<1e-15
        for bad in (np.array([0.,1]),np.array([np.nan,1])):
            try:target_divergence(bad,np.ones(2),.5)
            except ValueError:pass
            else:raise AssertionError('invalid support accepted')
        # Post-build control against the previous binary's saved zero-gravity state.
        control=run(out/'control',0.,'density',4)
        old=ROOT/'out/runs/volume_contract_20261024/density_error'
        names=['final.npy']+[p.name for p in old.glob('final_*.npy')]
        assert all(np.array_equal(np.load(old/n),np.load(out/'control/bifrost'/n)) for n in names)
        rec['legacy_final_arrays_identical']=len(names)
        local=run(out/'local',.5,'density',1)
        assert local['final_position_max']<4e-6 and local['final_velocity_max']<2e-5
        rec['cases']['local']=local
        density=run(out/'density',.5,'density',400);rec['cases']['density']=density
        for key in ('python','bifrost'):
            ratio=density[key][-1]['density_log_rms']/density['initial']['density_log_rms']
            assert ratio<=.8,(key,'density relaxation',ratio)
        # Both arms use the same already-admissible positions, and the one change is tau.
        for name,tau in [('rest_off',0.),('rest_on',.5)]:
            row=run(out/name,tau,'rest',2400);rec['cases'][name]=row
            for key in ('python','bifrost'):
                final=row[key][-1];ratio=max(r['energy'] for r in row[key])/row['initial']['energy']
                drift=np.linalg.norm(np.array(final['com'])-row['initial']['com'])
                assert final['rms']<=.1 and drift<=.1 and ratio<=1.005,(name,key,final,drift,ratio)
        rec['passed']=True
    except BaseException as e:
        rec['failure']=str(e);raise
    finally:
        (out/'receipt.json').write_text(json.dumps(rec,indent=2))
        print(json.dumps(dict(passed=rec['passed'],failure=rec.get('failure')),indent=2),flush=True)


if __name__=='__main__':main()
