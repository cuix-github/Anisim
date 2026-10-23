"""Causal stage replacement and particle energy accounting on a shared state.

Shared-field arms isolate stages; they are NOT independent-solver acceptance.
All actual graph arrays are compared to a previously retained control first.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'experiments'),str(ROOT/'gates')]
from phase_contract_2d import make
from loop3d import run_bifrost3,FINAL_STATE
from audit_pressure_replay import read_dump,flux

def delta(a,b):
    return dict(max=float(abs(a-b).max()),different=int(np.count_nonzero(a!=b)))

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);source=a.source.resolve()
    old=source/'bifrost'
    arrays={n:np.load(old/(n+'.npy')) for n in ('positions','velocities','particle_phase',
        'phase_st','mass','momx','momy','momz','tau','wt','wtph','scale','radius','spray')}
    dt=float(np.float32(.0025));tol=float(np.float32(1e-10))
    def make_state():
        s=make('plunge',2,tol);s.pos=arrays['positions'][:,:2].astype(float)
        s.vel=arrays['velocities'][:,:2].astype(float);s.typ=arrays['particle_phase'].astype(int)
        s._port_momentum=np.stack([arrays['momx'],arrays['momy']],1).astype(float)
        return s
    s=make_state();d=out/'control';os.environ['PFFLIP_PRESSURE_DUMP']=str(out/'pressure_double.bin')
    (out/'job.json').write_text(json.dumps(dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        source=str(source),kind='Single-state causal stage audit; shared fields are diagnostics only'),indent=2))
    run_bifrost3('causality',64,80,1,lambda x,y,z:y<18,1,dt,fused=True,init_arrays=arrays,
        extra_ports=dict(phase_model=2,phase_alpha=1.,rho0_face=s.rho0_face,escape=0,st=0,adapt=0,
                         emit_model=0,tolerance=tol),final_state=True,dumps=('ph','vel'),
        out_dir=str(d),log_prefix=str(out/'bifcmd'))
    names=['final.npy']+['final_'+n+'.npy' for n in FINAL_STATE]
    for n in names:assert np.array_equal(np.load(old/n),np.load(d/n)),n
    assert (source/'pressure_double.bin').read_bytes()==(out/'pressure_double.bin').read_bytes()
    coeff,rhs,pressure,before3=read_dump(out/'pressure_double.bin')
    after3=[v-dt*f for v,f in zip(before3,flux(pressure,coeff))]
    before=[v[:,:,0] for v in before3[:2]];after=[v[:,:,0] for v in after3[:2]]
    pre=[v.copy() for v in before];pre[1][:,1:-1]-=s.g*dt
    targets=dict(pos=np.load(d/'final.npy')[:,:2],vel=np.load(d/'final_velocities.npy')[:,:2],
                 mom=np.stack([np.load(d/('final_mom'+n+'.npy')) for n in 'xy'],1))
    results=dict(control_bit_identical=names,pressure_dump_identical=True,arms={})
    for label in ('independent','shared_pressure','shared_p2g_and_pressure'):
        s=make_state();original_p2g=s.p2g;stages={}
        def p2g():
            original_p2g()
            stages['p2g_u']=delta(s.u,pre[0]);stages['p2g_v']=delta(s.v,pre[1])
            if label=='shared_p2g_and_pressure':s.u,s.v=[v.copy() for v in pre]
        s.p2g=p2g
        if label!='independent':
            def project(dt):s.u,s.v=[v.copy() for v in after]
            s.project=project
        s.step(dt)
        for n,v in [('pos',s.pos),('vel',s.vel),('mom',s._port_momentum)]:
            stages[n]=delta(v.astype(np.float32),targets[n])
        results['arms'][label]=stages
    # Exact algebraic breakdown of G2P, using identical C++ double face fields.
    s=make_state();s.p2g();m=s.masses();alpha=np.where(s.typ==1,s.alpha[1],s.alpha[0])[:,None]
    sample=lambda fields:np.stack([s._sample_faces(v,np.array(off),s.pos) for v,off in
                                   zip(fields,((0,.5),(.5,0)))],1)
    interp0=sample(pre);interp_gravity=sample(before);interp1=sample(after)
    vblend=alpha*s.vel+(1-alpha)*interp0
    vgravity=vblend+interp_gravity-interp0
    vpressure=vgravity+interp1-interp_gravity
    ph=s._face_phase();phi=.5*sum(s._sample_faces(v,np.array(off),s.pos) for v,off in zip(ph,((0,.5),(.5,0))))
    k=lambda v:.5*m*np.sum(v*v,axis=1)
    work=dict(blend=k(vblend)-k(s.vel),gravity=k(vgravity)-k(vblend),pressure=k(vpressure)-k(vgravity))
    groups={}
    for label,mask in [('water',s.typ==1),('air',s.typ==0),('water_phi_below_03',(s.typ==1)&(phi<.3)),
        ('water_phi_03_05',(s.typ==1)&(phi>=.3)&(phi<.5)),('water_phi_05_09',(s.typ==1)&(phi>=.5)&(phi<.9)),
        ('water_phi_ge_09',(s.typ==1)&(phi>=.9))]:
        groups[label]=dict(count=int(mask.sum()),**{n:float(v[mask].sum()) for n,v in work.items()},
            positive_pressure_work=float(np.maximum(work['pressure'][mask],0).sum()))
    ids=np.argsort(work['pressure'])[-40:][::-1]
    results['g2p_energy_by_group']=groups
    results['top_pressure_gain_particles']=[dict(index=int(i),type=int(s.typ[i]),phi=float(phi[i]),
        pos=s.pos[i].tolist(),old_velocity=s.vel[i].tolist(),new_velocity=vpressure[i].tolist(),
        pressure_work=float(work['pressure'][i])) for i in ids]
    results['g2p_kinetic_change']=float((k(vpressure)-k(s.vel)).sum())
    np.savez_compressed(out/'energy_fields.npz',pos=s.pos,typ=s.typ,phi=phi,pressure_work=work['pressure'],
        old_velocity=s.vel,g2p_velocity=vpressure,pre_u=pre[0],pre_v=pre[1],post_u=after[0],post_v=after[1],
        phi_u=ph[0],phi_v=ph[1])
    (out/'receipt.json').write_text(json.dumps(results,indent=2))
    print(json.dumps({k:v for k,v in results.items() if k!='top_pressure_gain_particles'},indent=2))

if __name__=='__main__':main()
