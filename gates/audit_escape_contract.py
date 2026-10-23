"""Common-state escape/phase/flight audit; never silently changes the solver.

Rebuild common momentum once from a saved float32 velocity, shared by both
implementations. This is a single-step test, not exact trajectory continuation.
The paper threshold comparison uses documented averaged MAC-face interpolation.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'experiments'),str(ROOT/'gates'),str(ROOT/'graphs')]
from phase_repair_2d import make, ports, D, E, ESCAPE
from spray_events_2d import load_frame, NX, NY
from loop3d import run_bifrost3


def common(state,step,dt):
    s=make('escape',2,1e-10)
    s.pos=state['pos'].astype(float);s.vel=state['vel'].astype(float)
    s.typ=state['typ'].astype(np.int32);s.spray=state['spr'].copy()
    s.escaped=state['esc'].copy();s.radius=state['rad'].astype(float)
    n=len(s.pos);s.scale=np.ones(n);s.tau=np.zeros(n);s.xi=np.ones(n)
    s.nstep=step;s.t=step*dt
    m=s.masses().astype(np.float32);mom=(m[:,None]*s.vel).astype(np.float32)
    s._port_momentum=mom.astype(float)
    z=np.zeros(n,np.float32);one=np.ones(n,np.float32);ph=s.typ.astype(np.float32)
    arr=dict(positions=np.c_[s.pos,np.full(n,.5)],velocities=np.c_[s.vel,z],
        particle_phase=ph,phase_st=ph,mass=m,momx=mom[:,0],momy=mom[:,1],momz=z,
        tau=z,wt=one,wtph=ph,scale=one,radius=s.radius,spray=s.spray)
    return s,arr


def distribution(x):
    return dict(count=len(x),quantiles=np.quantile(x,[0,.1,.5,.9,1]).tolist()) if len(x) else dict(count=0)


def classification(s):
    s.p2g();fu,fv=s._face_frames();v=(~s.spray).astype(float)
    au,lu=s._splat(s.u.shape,fu,[v*(s.typ==0),v*(s.typ==1)])
    av,lv=s._splat(s.v.shape,fv,[v*(s.typ==0),v*(s.typ==1)])
    def sample(u,w):
        return .5*(s._sample_faces(u,np.array([0,.5]),s.pos)+s._sample_faces(w,np.array([.5,0]),s.pos))
    air,liquid=sample(au,av)/s.rho0_face,sample(lu,lv)/s.rho0_face
    phi=sample(*s._face_phase())
    current=(s.typ==1)&~s.spray&(air>1-s.esc_phi)
    # Paper default delta_phi=.2: liquid candidate phi<.3, bubble phi>.7.
    paper=(s.typ==1)&~s.spray&(phi<.3)
    bubble=(s.typ==0)&(liquid>1-s.esc_phi)
    paper_bubble=(s.typ==0)&(phi>.7)
    row=dict(current_water_candidates=int(current.sum()),phase_water_candidates=int(paper.sum()),
        both=int((current&paper).sum()),count_only=int((current&~paper).sum()),
        phase_only=int((paper&~current).sum()),current_candidate_phi=distribution(phi[current]),
        count_threshold_margin=distribution(abs(air[(s.typ==1)&~s.spray]-(1-s.esc_phi))),
        current_bubbles=int(bubble.sum()),phase_bubbles=int(paper_bubble.sum()),
        bubble_disagreements=int((bubble!=paper_bubble).sum()))
    return row,dict(phi=phi,air=air,liquid=liquid,current=current,paper=paper,bubble=bubble)


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    job=json.loads((a.source/'job.json').read_text());dt=job['dt'];rows=[]
    (out/'job.json').write_text(json.dumps(dict(source=str(a.source.resolve()),
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope=__doc__,phase_sampling='average of bilinear interpolants on the two MAC face fields',
        pressure_tolerance=float(np.float32(1e-10)),emission='off in common-state audit; existing spray retained'),indent=2))
    for j in (1280,3200,7200):
        state=load_frame(a.source,'bifrost',j);s,arr=common(state,j+1,dt)
        report,fields=classification(s);n=len(s.pos)
        d=out/f'step_{j}';d.mkdir()
        np.savez_compressed(d/'classification.npz',pos=s.pos,typ=s.typ,spr=s.spray,**fields)
        capture={};g2p=s.g2p;radii=s._droplet_radii;flight=s._droplet_fly
        def watch_g2p(u,v,skip=None):
            capture['skip']=skip.copy();g2p(u,v,skip)
        def watch_radii(drop,h):
            capture['drop']=drop.copy();capture['radius_before']=s.radius.copy()
            velocity=s.vel.copy();gone=radii(drop,h)
            capture['radius_after']=s.radius.copy();capture['gone']=gone.copy()
            capture['velocity_changed_by_radius_assignment']=bool(np.any(s.vel!=velocity))
            return gone
        def watch_flight(drop,times):
            capture['flight_before']=s.vel[drop].copy();flight(drop,times)
            capture['flight_after']=s.vel[drop].copy()
        s.g2p=watch_g2p;s._droplet_radii=watch_radii;s._droplet_fly=watch_flight
        s.step(dt)
        assert np.array_equal(capture['skip'],fields['current']|fields['bubble']|state['spr'])
        assert not capture['velocity_changed_by_radius_assignment']
        run_bifrost3('escape_contract',NX,NY,1,lambda x,y,z:y<18,1,dt,
            init_arrays=arr,step_offset=j+1,fused=True,final_state=True,dumps=('ph','vel','rad','spr','esc'),
            extra_ports=dict(ESCAPE,**ports(dict(E,model=0),D),phase_model=2,phase_alpha=1.,
                rho0_face=s.rho0_face,tolerance=float(np.float32(1e-10))),
            out_dir=str(d/'bifrost'),log_prefix=str(d/'bifcmd'))
        actual=dict(pos=np.load(d/'bifrost/final.npy')[:,:2],vel=np.load(d/'bifrost/final_velocities.npy')[:,:2],
            typ=np.load(d/'bifrost/final_phase.npy')>.5,rad=np.load(d/'bifrost/final_radius.npy'),
            spr=np.load(d/'bifrost/final_spray.npy')>.5,esc=np.load(d/'bifrost/esc.0000.npy')>.5)
        assert len(s.pos)==len(actual['pos'])
        err={k:float(abs(actual[k].astype(float)-getattr(s,nm)).max()) for k,nm in
             [('pos','pos'),('vel','vel'),('rad','radius'),('typ','typ'),('spr','spray'),('esc','escaped')]}
        assert err['pos']<2e-5 and err['vel']<1e-4 and err['rad']<1e-9 and not (err['typ'] or err['spr'] or err['esc']),err
        mass=np.load(d/'bifrost/final_mass.npy');escaped=actual['esc']&actual['typ']&~actual['spr']
        new=capture['drop']&(capture['radius_before']<=0)&~state['spr']
        report.update(step=j,time=(j+1)*dt,actual_graph_max_errors=err,
            input_primary=int((~state['spr']).sum()),output_primary=int((~actual['spr']).sum()),
            retained_escaped_primary_mass=distribution(mass[escaped]),
            output_spray_mass=distribution(mass[actual['spr']]),
            newly_escaped_radii_metres=distribution(capture['radius_after'][new]),
            radius_assignment_velocity_continuous=True,
            all_droplet_flight_velocity_change=distribution(np.linalg.norm(capture['flight_after']-capture['flight_before'],axis=1)),
            evaporated_primary=int((capture['gone']&~state['spr']).sum()))
        rows.append(report);(out/'receipt.json').write_text(json.dumps(rows,indent=2))
        print(json.dumps(report),flush=True)


if __name__=='__main__':main()
