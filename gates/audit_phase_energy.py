"""Common-state full-graph one-step energy audit, not trajectory continuation."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'experiments'),str(ROOT/'gates')]
from phase_contract_2d import make
from analyze_phase_contract import load
from loop3d import run_bifrost3
from audit_pressure_replay import read_dump,flux,operator


def energy(pos,vel,typ):
    m=np.where(typ==1,1000.,1.)
    return dict(kinetic=float(.5*np.sum(m[:,None]*vel.astype(float)**2)),
        potential=float(np.sum(m*float(np.float32(9.8))*pos[:,1])))


def total(e):return e['kinetic']+e['potential']


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);source=a.source.resolve()
    job=json.loads((source/'job.json').read_text());dt=job['dt'];records=[]
    (out/'job.json').write_text(json.dumps(dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source=str(source),dt=dt,scope='Common-state one-step audit. Momentum reconstructed once from saved float32 velocity, shared by both implementations; not exact original trajectory replay.'),indent=2))
    for j in (1120,1200,1280,1360,1440):
        state=load(source,'bifrost',j)
        for model in (1,2):
            s=make('plunge',model,float(np.float32(1e-10)))
            s.pos=state['pos'].astype(float);s.vel=state['vel'].astype(float);s.typ=state['typ'].copy()
            n=len(s.pos);m=np.where(s.typ==1,1000.,1.).astype(np.float32)
            mom=(m[:,None]*s.vel).astype(np.float32)
            s._port_momentum=mom.astype(float)
            z=np.zeros(n,np.float32);one=np.ones(n,np.float32);phase=s.typ.astype(np.float32)
            arrays=dict(positions=np.c_[s.pos,np.full(n,.5)],particle_phase=phase,phase_st=phase,
                velocities=np.c_[s.vel,z],mass=m,momx=mom[:,0],momy=mom[:,1],momz=z,
                tau=z,wt=one,wtph=phase,scale=one,radius=z,spray=z)
            d=out/f'step_{j}_phase_{model}';d.mkdir()
            os.environ['PFFLIP_PRESSURE_DUMP']=str(d/'pressure_double.bin')
            initial=energy(s.pos,s.vel,s.typ);trace={}
            original_g2p=s.g2p
            def g2p(*args,**kw):
                original_g2p(*args,**kw)
                trace['after_g2p']=energy(s.pos,s.vel,s.typ)
            s.g2p=g2p
            s.step(dt);trace['after_advection']=energy(s.pos,s.vel,s.typ)
            run_bifrost3('energy',s.nx,s.ny,1,lambda x,y,z:y<18,1,dt,fused=True,
                init_arrays=arrays,extra_ports=dict(phase_model=model,phase_alpha=1.,rho0_face=s.rho0_face,
                    escape=0,st=0,adapt=0,emit_model=0,tolerance=float(np.float32(1e-10))),
                final_state=True,dumps=('ph','vel'),out_dir=str(d/'bifrost'),log_prefix=str(d/'bifcmd'))
            cpp_pos=np.load(d/'bifrost/final.npy')[:,:2];cpp_vel=np.load(d/'bifrost/final_velocities.npy')[:,:2]
            cpp=energy(cpp_pos,cpp_vel,state['typ'])
            coeff,rhs,pressure,before=read_dump(d/'pressure_double.bin')
            after=[v-dt*f for v,f in zip(before,flux(pressure,coeff))]
            grid=lambda fields:sum(float(np.sum(v[b>0]**2/b[b>0])) for b,v in zip(coeff,fields))*.5
            residual=float(np.linalg.norm(rhs-operator(pressure,coeff))/np.linalg.norm(rhs))
            dv=float(abs(cpp_vel-s.vel).max());de=abs(total(cpp)-total(trace['after_advection']))/total(initial)
            record=dict(step=j,time=(j+1)*dt,phase_model=model,initial=initial,python=trace,cpp=cpp,
                grid_projection_kinetic_before=grid(before),grid_projection_kinetic_after=grid(after),
                true_residual=residual,velocity_difference=dv,total_energy_relative_difference=de,
                cpp_total_energy_step_delta=total(cpp)-total(initial))
            records.append(record)
            (out/'receipt.json').write_text(json.dumps(records,indent=2))
            assert residual<1.1e-10 and grid(after)<=grid(before)*(1+1e-10),(j,model,record)
            assert dv<5e-5 and de<1e-6,(j,model,record)
            print(f'step={j} phase={model} energy_delta={record["cpp_total_energy_step_delta"]:.3f} grid_delta={grid(after)-grid(before):.3f} dv={dv:.3g}',flush=True)
    os.environ.pop('PFFLIP_PRESSURE_DUMP',None)


if __name__=='__main__':main()
