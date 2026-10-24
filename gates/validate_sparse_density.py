"""One finite-density candidate, with failures retained and no long-run admission."""
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
sys.path[:0]=[str(ROOT/p) for p in ('gates','reference2d','graphs')]
from volume_contract import finite_target_divergence, target_divergence
from validate_interface_closure import make,cpp,py_step,compare,rest,DT,LIMITS
from validate_volume_contract import run as density_run
from validate_base_combination import identical
from diagnose_interface_support import state
from validate_3d import run_node


def write(p,v):
    p.write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')


def math_checks():
    n=np.array([0,.5,1,1.5,2.]); ref=np.ones(5)
    d=finite_target_divergence(n,ref,.5)
    assert np.array_equal(d,[-2,-1,0,1,2])
    for bad,reference,tau in [([-1,1],[1,1],.5),([0,0],[1,1],.5),
        ([np.nan,1],[1,1],.5),([1,1],[0,1],.5),([1,1],[1,1],0),
        ([1e308,1],[1e-308,1],.5)]:
        try:finite_target_divergence(bad,reference,tau)
        except ValueError:pass
        else:raise AssertionError('invalid target accepted')
    eps=1e-4
    a=finite_target_divergence([1-eps,1+eps],[1,1],.5)
    b=target_divergence([1-eps,1+eps],[1,1],.5)
    assert np.max(abs(a-b))<1e-10
    return dict(target=d.tolist(),near_reference_difference=float(np.max(abs(a-b))))


def controls(out):
    rows={}
    for name,tau,steps in [('off',0.,4),('log',.5,1)]:
        density_run(out/name,tau,'density',steps)
        rows[name]=identical(ROOT/f'out/runs/sparse_density_20261024/before_r2/{name}/bifrost',out/name/'bifrost')
    return rows


def projection(out):
    """Prescribed grid isolates the new target; direct distances are the oracle."""
    nx,ny=16,12
    s=make(nx,ny,mode=0);s.typ[:]=0
    regular=s.pos.copy();pos=regular.copy()
    pos[np.linalg.norm(pos-[.5,9.5],axis=1)<1,0]+=3
    centres=np.stack(np.meshgrid(np.arange(nx)+.5,np.arange(ny)+.5,indexing='ij'),-1).reshape(-1,2)
    def count(points):
        return np.array([np.sum(np.maximum(1-np.sum((points-c)**2,axis=1),0.)**3) for c in centres]).reshape(nx,ny)
    n,ref=count(pos),count(regular)
    raw=(n/ref-1)/.5; target=raw-raw.mean()
    assert n[0,9]==0 and abs(target.sum())<1e-11
    inp=out/'input';inp.mkdir(parents=True)
    for axis,shape in [('u',(nx+1,ny)),('v',(nx,ny+1)),('w',(nx,ny,2))]:
        for suffix,val in [('mass',1.),('mom',0.),('phase',0.)]:
            np.save(inp/f'{axis}_{suffix}.npy',np.full(shape,val,np.float32).ravel())
    for name,val in dict(positions=np.c_[pos,np.full(len(pos),.5)].astype(np.float32),
        velocities=np.zeros((len(pos),3),np.float32),particle_phase=np.zeros(len(pos),np.float32)).items():
        np.save(inp/f'{name}.npy',val)
    ports=dict(nx=nx,ny=ny,nz=1,dt=DT,gravity=0,phase_model=2,escape=0,
        droplet_model=0,emit_model=0,volume_tau=.5,rho0_face=s.rho0_face,tolerance=1e-6)
    # Same input, old model: the old refusal must not silently change.
    try:run_node(str(inp),str(out/'legacy_refusal'),volume_model=0,**ports)
    except SystemExit as e:assert 'refuses empty/nonfinite particle support' in str(e),str(e)
    else:raise AssertionError('legacy law accepted empty support')
    stats=run_node(str(inp),str(out/'finite'),volume_model=1,**ports)
    u=np.load(out/'finite/out_u.npy').reshape(nx+1,ny).astype(float)
    v=np.load(out/'finite/out_v.npy').reshape(nx,ny+1).astype(float)
    div=np.diff(u,axis=0)+np.diff(v,axis=1)
    error=float(abs(div-target).max())
    assert error<5e-5 and stats['final_residual']<=1e-6
    assert np.all(u[[0,-1]]==0) and np.all(v[:,[0,-1]]==0)
    refusals={}
    for label,override in [('invalid_model',dict(volume_model=2)),('negative_tau',dict(volume_model=1,volume_tau=-1.))]:
        p=dict(ports);p.update(override)
        try:run_node(str(inp),str(out/label),**p)
        except SystemExit as e:
            assert 'volume contract' in str(e),str(e);refusals[label]=str(e)
        else:raise AssertionError('invalid configuration accepted')
    return dict(empty_centres=int((n==0).sum()),max_target=float(abs(target).max()),
        divergence_target_max_error=error,pressure_residual=stats['final_residual'],
        legacy_refusal_retained=True,invalid_configuration_refusals=refusals)


def common(out):
    rows=[]
    for step in (1280,3200,7200):
        s=state(step);s.volume_tau=.5;s.volume_model=1
        holes=int((s.number_density()==0).sum());water=int(((s.typ==1)&~s.spray).sum())
        a,res=cpp(s,out/f'step_{step}',1);py_step(s)
        row=dict(step=step,empty_centres=holes,**compare(s,a),residual=res)
        assert int(((a['particle_phase']==1)&(a['spray']<.5)).sum())==water
        row['water_inventory_unchanged']=True
        rows.append(row);write(out/'progress.json',rows)
    return rows


def sparse_and_restart(out):
    s=make(tau=.5);s.volume_model=1
    selected=np.linalg.norm(s.pos-[.5,9.5],axis=1)<1
    s.pos[selected,0]+=3
    assert np.any(s.number_density()==0)
    a,res=cpp(s,out/'sparse',1);py_step(s)
    row=dict(**compare(s,a),residual=res)
    s=make(tau=.5);s.volume_model=1
    a,_=cpp(s,out/'full',80,every=40)
    b,_=cpp(s,out/'first',40,every=40)
    c,_=cpp(s,out/'second',40,b,40,every=40)
    checks={k:bool(np.array_equal(a[k],c[k])) for k in a}
    assert all(checks.values()),checks
    return dict(sparse=row,restart=checks)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--stage',choices=['local','density','rest'],default='local')
    ap.add_argument('--requires',type=Path)
    args=ap.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    if args.stage!='local':
        assert args.requires and json.loads(args.requires.read_text())['passed'],'prerequisite must pass'
    out=args.out.resolve();out.mkdir(parents=True,exist_ok=False)
    os.environ['OMP_NUM_THREADS']='1'
    binary=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    r=dict(passed=False,stage=args.stage,source_sha=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),scope=__doc__,limits=LIMITS)
    write(out/'job.json',r);start=time.perf_counter()
    try:
        if args.stage=='local':
            for label,fn in [('controls',controls),('math',lambda p:math_checks()),
                ('projection',projection),('common',common),('sparse_restart',sparse_and_restart)]:
                (out/label).mkdir()
                r[label]=fn(out/label);write(out/'receipt.json',r)
                print('PASS '+label,flush=True)
        elif args.stage=='density':
            row=density_run(out/'finite',.5,'density',400,model=1);r['density']=row
            write(out/'receipt.json',r)
            for backend in ('python','bifrost'):
                ratio=row[backend][-1]['density_log_rms']/row['initial']['density_log_rms']
                assert ratio<=.8,(backend,'density relaxation',ratio)
        else:r['rest']=rest(out,model=1,taus=(.5,))
        r['passed']=True
    except BaseException as e:r['failure']=repr(e);raise
    finally:
        r['seconds']=time.perf_counter()-start;write(out/'receipt.json',r)
        print(json.dumps(dict(passed=r['passed'],stage=args.stage,seconds=r['seconds'],failure=r.get('failure'))),flush=True)


if __name__=='__main__':main()
