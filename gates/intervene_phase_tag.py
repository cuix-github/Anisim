"""Change only cell labels at a recorded first divergence; diagnostic only."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'experiments')]
from phase_contract_2d import make
from phase_contract import compress,pure_faces


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    (out/'job.json').write_text(json.dumps(dict(kind='Single cell-label intervention; diagnostic only',
        source=str(a.source.resolve()),commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()),indent=2))
    z=np.load(a.source/'first_cell_tag.npz');dt=float(np.float32(.0025));tol=float(np.float32(1e-10))
    s=make('plunge',2,tol);s.pos=z['pos'];s.vel=z['vel'];s._port_momentum=z['momentum'];s.p2g()
    assert np.array_equal(s.u,z['py_pre_u']) and np.array_equal(s.v,z['py_pre_v'])
    assert all(np.array_equal(x,z['py_phi_'+n]) for x,n in zip(s._face_phase(),'uv'))
    other=make('plunge',2,tol);other.pos=z['cpp_pos'];other.p2g()
    raw=[compress(v,s.rho0_face,s.rho_l,s.rho_g,s.phase_alpha) for v in (s.mu,s.mv)]
    labels_only=tuple(x.astype(np.float32).astype(float) for x in pure_faces(raw,other.cell_phase))
    phases={'own':s._face_phase(),'only_cpp_cell_labels':labels_only,'all_cpp_face_phase':(z['cpp_phi_u'],z['cpp_phi_v'])}
    projected={};results={}
    for name,phase in phases.items():
        t=make('plunge',2,tol);t.u=z['py_pre_u'].copy();t.v=z['py_pre_v'].copy()
        t._contract_faces=phase;t.v+=t.g*dt;t._enforce_walls();t.project(dt)
        projected[name]=(t.u.copy(),t.v.copy())
        results[name]=dict(phase_max_difference=max(float(abs(x-y).max()) for x,y in zip(phase,phases['own'])),
            projected_velocity_max_difference=max(float(abs(x-y).max()) for x,y in zip(projected[name],projected['own'])))
    result=dict(cells_different=int(np.count_nonzero((s.cell_phase>=.5)!=(other.cell_phase>=.5))),arms=results,
        labels_only_vs_all_cpp_projected_max=max(float(abs(x-y).max()) for x,y in zip(projected['only_cpp_cell_labels'],projected['all_cpp_face_phase'])))
    (out/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
