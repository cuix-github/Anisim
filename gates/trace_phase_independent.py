"""Locate the first differences with independent pressure solves each step."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import time
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'experiments'),str(ROOT/'gates')]
from phase_contract_2d import make
from trace_phase_sequence import read_record,delta


def main():
    p=argparse.ArgumentParser();p.add_argument('--trace',type=Path,required=True)
    p.add_argument('--control',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    dt=float(np.float32(.0025));tol=float(np.float32(1e-10));s=make('plunge',2,tol)
    other=make('plunge',2,tol);cpp_prev=s.pos.copy();original_p2g=s.p2g;original_project=s.project
    (out/'job.json').write_text(json.dumps(dict(kind='Independent pressure, every-step causal trace',
        trace=str(a.trace.resolve()),control=str(a.control.resolve()),dt=dt,tolerance=tol,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()),indent=2))
    first={};rows=[];started=time.perf_counter();d=a.trace/'bifrost'
    with open(a.trace/'stages.bin','rb') as f:
        for j in range(2400):
            fields=read_record(f,64,80);stages={};event=dict(step=j,time_start=j*dt)
            def p2g():
                original_p2g()
                stages['pre_velocity_max']=max(float(abs(s.u-fields[0]).max()),float(abs(s.v-fields[1]).max()))
                ph=s._face_phase();stages['phi_max']=max(float(abs(x-y).max()) for x,y in zip(ph,fields[4:]))
                # Reconstruct the cells from C++ input positions using the
                # independent splat; actual exported face fields remain primary evidence.
                other.pos=cpp_prev;other.p2g()
                diff=(s.cell_phase>=.5)!=(other.cell_phase>=.5)
                stages['cell_tags_different']=int(diff.sum())
                if diff.any() and 'cell_tag' not in first:
                    cells=np.argwhere(diff)
                    first['cell_tag']=dict(**event,cells=[dict(index=ij.tolist(),python=float(s.cell_phase[tuple(ij)]),
                        cpp=float(other.cell_phase[tuple(ij)])) for ij in cells])
                    np.savez_compressed(out/'first_cell_tag.npz',pos=s.pos,vel=s.vel,typ=s.typ,momentum=s._port_momentum,
                        cpp_pos=cpp_prev,py_pre_u=s.u,py_pre_v=s.v,py_phi_u=ph[0],py_phi_v=ph[1],
                        cpp_pre_u=fields[0],cpp_pre_v=fields[1],cpp_phi_u=fields[4],cpp_phi_v=fields[5])
            def project(dt):
                original_project(dt)
                stages['projected_velocity_max']=max(float(abs(s.u-fields[2]).max()),float(abs(s.v-fields[3]).max()))
            s.p2g=p2g;s.project=project;s.step(dt)
            s.pos=s.pos.astype(np.float32).astype(float);s.vel=s.vel.astype(np.float32).astype(float)
            cpp_pos=np.load(d/f'frame.{j:04d}.npy')[:,:2];cpp_vel=np.load(d/f'vel.{j:04d}.npy')[:,:2]
            stages['pos']=delta(s.pos,cpp_pos);stages['vel']=delta(s.vel,cpp_vel)
            for name in ('pos','vel'):
                if stages[name]['different'] and name not in first:first[name]=dict(**event,stages=stages.copy())
            if j==0:first['initial_stage_comparison']=stages.copy()
            rows.append(dict(**event,**stages));cpp_prev=cpp_pos.astype(float)
            # Verify this remains the previous independent trajectory at every
            # retained sample, so instrumentation cannot manufacture the divergence.
            old=a.control/'python'/f'frame.{j:04d}.npz'
            if old.exists():
                z=np.load(old);assert np.array_equal(s.pos,z['pos']) and np.array_equal(s.vel,z['vel']),j
            if (j+1)%400==0:print(f'Independent: {j+1}/2400 ({time.perf_counter()-started:.1f}s)',flush=True)
    final=np.load(a.control/'python'/'final.npz')
    assert np.array_equal(s.pos,final['pos']) and np.array_equal(s.vel,final['vel'])
    result=dict(first=first,control_samples_and_final_bit_identical=True,rows=rows,wall_seconds=time.perf_counter()-started)
    (out/'receipt.json').write_text(json.dumps(result,indent=2))
    # Causal counterfactual: same pre-pressure velocity and RHS; replace only
    # phase coefficients at the first class disagreement. Not a proposed fix.
    z=np.load(out/'first_cell_tag.npz');projected=[]
    for arm in ('py','cpp'):
        t=make('plunge',2,tol);t.u=z['py_pre_u'].copy();t.v=z['py_pre_v'].copy()
        t._contract_faces=(z[arm+'_phi_u'],z[arm+'_phi_v']);t.v+=t.g*dt;t._enforce_walls();t.project(dt)
        projected.append((t.u.copy(),t.v.copy()))
    result['first_tag_phase_only_projected_velocity_max']=max(float(abs(x-y).max()) for x,y in zip(*projected))
    (out/'receipt.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))


if __name__=='__main__':main()
