"""Trace a real graph, then isolate stages over a sequence (not acceptance).

Each arm retains its own particle state. Sharing pressure velocity deliberately
removes one independent stage to locate the origin of trajectory differences.
The traced graph must match a pre-instrumentation retained control bit for bit.
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
sys.path[:0]=[str(ROOT/'experiments'),str(ROOT/'gates')]
from phase_contract_2d import make,run_cpp
from loop3d import FINAL_STATE


def delta(a,b):
    return dict(max=float(abs(a-b).max()),different=int(np.count_nonzero(a!=b)))


def read_record(f,nx,ny):
    dims=np.fromfile(f,np.int32,3)
    assert np.array_equal(dims,[nx,ny,1]),dims
    fields=[]
    for shape in [(nx+1,ny),(nx,ny+1)]*3:
        fields.append(np.fromfile(f,np.float64,np.prod(shape)).reshape(shape))
    return fields


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--control',type=Path,required=True);p.add_argument('--seconds',type=float,default=6)
    p.add_argument('--reuse-graph',action='store_true')
    a=p.parse_args();out=a.out.resolve();control=a.control.resolve()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out.mkdir(parents=True,exist_ok=a.reuse_graph)
    dt=float(np.float32(.0025));tol=float(np.float32(1e-10));steps=round(a.seconds/.0025)
    trace=out/'stages.bin';d=out/'bifrost'
    if not a.reuse_graph:
        (out/'job.json').write_text(json.dumps(dict(kind='Causal sequence stage replacement, NOT independent validation',
            commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            binary_sha256=hashlib.sha256((ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll').read_bytes()).hexdigest(),
            steps=steps,dt=dt,tolerance=tol,control=str(control)),indent=2))
        os.environ['PFFLIP_FLAT_STAGE_TRACE']=str(trace)
        print('Running actual graph with stage trace',flush=True)
        run_cpp(d,make('plunge',2,tol),steps,dt,1)
        del os.environ['PFFLIP_FLAT_STAGE_TRACE']
    assert trace.stat().st_size==steps*(12+3*((65*80)+(64*81))*8)
    names=['final.npy']+['final_'+n+'.npy' for n in FINAL_STATE]
    for n in names:assert np.array_equal(np.load(d/n),np.load(control/n)),n
    print('Traced graph final state matches pre-instrumentation control',flush=True)
    result=dict(control_final_bit_identical=names,arms={})
    # Stop at the first mismatch, retaining both states and exact stage fields.
    for label in ('shared_pressure','shared_p2g_and_pressure'):
        s=make('plunge',2,tol);original=s.p2g;first=None;max_p2g=0.;max_phi=0.;started=time.perf_counter()
        with open(trace,'rb') as f:
            for j in range(steps):
                fields=read_record(f,64,80);stages={}
                def p2g():
                    original()
                    stages['pre_u']=delta(s.u,fields[0]);stages['pre_v']=delta(s.v,fields[1])
                    stages['phi_u']=delta(s._face_phase()[0],fields[4]);stages['phi_v']=delta(s._face_phase()[1],fields[5])
                    if label=='shared_p2g_and_pressure':s.u,s.v=[v.copy() for v in fields[:2]]
                def project(dt):s.u,s.v=[v.copy() for v in fields[2:4]]
                s.p2g=p2g;s.project=project;s.step(dt)
                s.pos=s.pos.astype(np.float32).astype(float);s.vel=s.vel.astype(np.float32).astype(float)
                pos=np.load(d/f'frame.{j:04d}.npy')[:,:2];vel=np.load(d/f'vel.{j:04d}.npy')[:,:2]
                stages['pos']=delta(s.pos,pos);stages['vel']=delta(s.vel,vel)
                max_p2g=max(max_p2g,stages['pre_u']['max'],stages['pre_v']['max'])
                max_phi=max(max_phi,stages['phi_u']['max'],stages['phi_v']['max'])
                if stages['pos']['different'] or stages['vel']['different']:
                    first=dict(step=j,time=(j+1)*dt,stages=stages)
                    np.savez_compressed(out/(label+'_first.npz'),py_pos=s.pos,py_vel=s.vel,cpp_pos=pos,cpp_vel=vel,
                        **dict(zip(('pre_u','pre_v','post_u','post_v','phi_u','phi_v'),fields)))
                    break
                if (j+1)%400==0:print(f'{label}: {j+1} steps bit exact ({time.perf_counter()-started:.1f}s)',flush=True)
        mom={}
        if first is None:
            for axis in (0,1):mom['xy'[axis]]=delta(s._port_momentum[:,axis],np.load(d/('final_mom'+'xy'[axis]+'.npy')))
        result['arms'][label]=dict(first_mismatch=first,steps_checked=j+1,max_p2g_difference=max_p2g,
            max_phi_difference=max_phi,final_momentum=mom,wall_seconds=time.perf_counter()-started)
        (out/'receipt.json').write_text(json.dumps(result,indent=2))
        print(json.dumps(result['arms'][label],indent=2),flush=True)


if __name__=='__main__':main()
