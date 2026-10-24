"""Generic classification contracts, actual C++ node, no plunge-specific tuning."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'gates'),str(ROOT/'reference2d')]
from validate_3d import run_node
from phase_contract import PhaseContractSim
from pfflip2d import Sim


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);rows=[]
    cutoff=float(np.float32(.3));nx,ny=16,12
    for kind in ('air','at_lower','above_lower','at_upper','above_upper','water','linear'):
        s=PhaseContractSim(nx,ny,escape=2,phase_model=2,esc_phi=cutoff,g=0)
        s.pos=np.array([[4.,4.],[4.,4.],[12.,4.],[12.,4.]])
        s.typ=np.array([1,0,1,0]);s.vel=np.zeros_like(s.pos)
        levels=dict(air=0,at_lower=cutoff,above_lower=float(np.nextafter(np.float32(.3),np.float32(1))),
                    at_upper=float(np.float32(1-cutoff)),above_upper=float(np.nextafter(np.float32(1-cutoff),np.float32(1))),water=1)
        u=np.zeros((nx+1,ny));v=np.zeros((nx,ny+1));w=np.zeros((nx,ny,2))
        if kind=='linear':
            u[:]=np.arange(nx+1)[:,None]/nx;v[:]=(np.arange(nx)[:,None]+.5)/nx;w[:]=(np.arange(nx)[:,None,None]+.5)/nx
        else:
            u[:]=levels[kind];v[:]=levels[kind];w[:]=levels[kind]
        s._contract_faces=(u.astype(np.float32).astype(float),v.astype(np.float32).astype(float))
        ph=s.escape_phase_samples();expected=np.where(s.typ==1,ph<cutoff,ph>1-cutoff)
        d=out/kind;ind=d/'input';ind.mkdir(parents=True)
        for name,grid in zip('uvw',(u,v,w)):
            for suffix,arr in [('mass',np.ones_like(grid)),('mom',np.zeros_like(grid)),('phase',grid)]:
                np.save(ind/f'{name}_{suffix}.npy',arr.astype(np.float32).ravel())
        np.save(ind/'positions.npy',np.c_[s.pos,np.full(4,.5)].astype(np.float32))
        np.save(ind/'velocities.npy',np.zeros((4,3),np.float32));np.save(ind/'particle_phase.npy',s.typ.astype(np.float32))
        stats=run_node(str(ind),str(d/'cpp'),nx=nx,ny=ny,nz=1,dt=.0025,gravity=0,
            phase_model=2,escape=2,esc_phi=cutoff,rho0_face=3.1015625,buoyancy=0)
        actual=np.load(d/'cpp/out_escaped.npy')>.5
        keep=~((s.typ==0)&expected)  # section 3.4 deletes escaped air, not water
        assert np.array_equal(actual,expected[keep]),(kind,actual,expected[keep])
        assert np.array_equal(np.load(d/'cpp/out_positions.npy')[:,:2],s.pos[keep])
        rows.append(dict(kind=kind,phase=ph.tolist(),expected_before_deletion=expected.tolist(),
            expected_surviving=expected[keep].tolist(),deleted_air=int((~keep).sum()),actual=actual.tolist()))
    for kw in (dict(phase_model=1),dict(esc_phi=0),dict(esc_phi=.5)):
        try:PhaseContractSim(8,8,escape=2,**kw)
        except ValueError:pass
        else:raise AssertionError(kw)
    try:Sim(8,8,escape=2)
    except ValueError:pass
    else:raise AssertionError('missing phase contract accepted')
    for i,kw in enumerate((dict(phase_model=0),dict(esc_phi=0),dict(esc_phi=.5))):
        args=dict(nx=nx,ny=ny,nz=1,dt=.0025,gravity=0,phase_model=2,escape=2,esc_phi=cutoff)
        args.update(kw)
        try:run_node(str(ind),str(out/f'refuse_{i}'),**args)
        except SystemExit as e:assert 'phase-based escape' in str(e)
        else:raise AssertionError(kw)
    (out/'receipt.json').write_text(json.dumps(dict(passed=True,cases=rows,python_rejections=4,cpp_rejections=3,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()),indent=2))
    print('PASS: seven actual-node classification fields; invalid contracts refused')


if __name__=='__main__':main()
