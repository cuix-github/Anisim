"""Five-second still-water physics anchor with phase escape + emission enabled."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'experiments'),str(ROOT/'reference2d'),str(ROOT/'graphs')]
from phase_contract import PhaseContractSim
from phase_repair_2d import D,E,ESCAPE,ports,py
from loop3d import run_bifrost3


def make():
    s=PhaseContractSim(48,40,phase_model=2,escape=2,esc_phi=float(np.float32(.3)),
        droplet=D,emit=E,spray_cfl=.5,g=float(np.float32(-9.8)),
        alpha_l=float(np.float32(.97)),alpha_g=float(np.float32(.9)))
    s.seed(lambda x,y:y<12,seed=1);s.calibrate();s.rho0_face=float(np.float32(s.rho0_face))
    s.pos=s.pos.astype(np.float32).astype(float);s.escaped=np.zeros(len(s.pos),bool)
    return s


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);s=make();dt=float(np.float32(.0025))
    initial=s.pos[s.typ==1].mean(0);p3=np.c_[s.pos,np.full(len(s.pos),.5)].astype(np.float32)
    (out/'job.json').write_text(json.dumps(dict(kind=__doc__,steps=2000,dt=dt,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        thresholds=dict(water_rms=.1,water_com_drift=.1,max_spray_births=0)),indent=2))
    os.environ['PFFLIP_SPRAY_AUDIT']='1'
    run_bifrost3('phase_escape_rest',48,40,1,lambda x,y,z:y<12,2000,dt,
        seed_fn=lambda *a:(p3,s.typ.copy()),fused=True,final_state=True,dump_every=40,
        dumps=('ph','vel','rad','spr','esc'),extra_ports=dict(ESCAPE,**ports(E,D),escape=2,
            phase_model=2,phase_alpha=1.,rho0_face=s.rho0_face),out_dir=str(out/'bifrost'),log_prefix=str(out/'bifrost/bifcmd'))
    py(out/'python',s,2000,dt,40)
    cpp_pos=np.load(out/'bifrost/final.npy')[:,:2];cpp_vel=np.load(out/'bifrost/final_velocities.npy')[:,:2]
    cpp_typ=np.load(out/'bifrost/final_phase.npy')>.5;cpp_spray=np.load(out/'bifrost/final_spray.npy')>.5
    ledger=json.loads((out/'python/ledger.json').read_text())
    result={}
    for arm,pos,vel,typ,spray in [('python',s.pos,s.vel,s.typ,s.spray),('bifrost',cpp_pos,cpp_vel,cpp_typ,cpp_spray)]:
        water=(typ==1)&~spray
        rms=float(np.sqrt(np.mean(np.sum(vel[water].astype(float)**2,axis=1))))
        drift=float(np.linalg.norm(pos[water].astype(float).mean(0)-initial))
        result[arm]=dict(water_rms=rms,com_drift=drift,spray_alive=int(spray.sum()),
                        pass_physics=rms<.1 and drift<.1 and not spray.any())
    result['python_births']=sum(r['accepted'] for r in ledger)
    cpp_events=[json.loads(line[len('SPRAY_EVENT '):]) for line in
        (out/'bifrost/bifcmd.stderr.log').read_text().splitlines() if line.startswith('SPRAY_EVENT ')]
    assert len(cpp_events)==2000
    result['cpp_births']=sum(r['accepted'] for r in cpp_events)
    result['final_position_max']=float(abs(s.pos-cpp_pos).max())
    result['final_velocity_max']=float(abs(s.vel-cpp_vel).max())
    (out/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
    assert all(result[k]['pass_physics'] for k in ('python','bifrost')) and result['python_births']==result['cpp_births']==0,result


if __name__=='__main__':main()
