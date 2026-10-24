"""Actual graph lifecycle gate at nz=1: passive coupling and exact resume.

ST/adapt on/off are tested independently of the 20s visual comparison.
All runs share the same sheared 2D seed and opt-in event parameters.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parent/"reference2d"),str(HERE.parent/"graphs")]
from pfflip2d import Sim
from loop3d import run_bifrost3, load_final_state
from validate_spray_events import D, E, ports
from validate_droplet_node import ESCAPE


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args(); out=args.out.resolve(); out.mkdir(parents=True,exist_ok=True)
    assert not subprocess.check_output(["git","status","--porcelain"],cwd=HERE.parent,text=True).strip()
    nx,ny,dt=40,32,float(np.float32(.01))
    rec=dict(passed=False, cases={})
    os.environ["PFFLIP_SPRAY_AUDIT"]="1"
    def run(name, steps, init, kwargs, offset=0):
        dest=out/name
        assert not dest.exists(), "Preserve previous runs"
        d,_,diag=run_bifrost3(name,nx,ny,1,lambda x,y,z:y<14,steps,dt,
            extra_ports=kwargs,init_arrays=init,step_offset=offset,fused=True,
            final_state=True,dumps=("spr","rad","ph","vel"),out_dir=str(dest),
            log_prefix=str(dest/"bifcmd"))
        assert np.isfinite(diag).all() and (diag[:,0]<4000).all()
        return load_final_state(d)
    try:
        for st,adapt in ((0,0),(1,0),(0,1),(1,1)):
            sim=Sim(nx,ny,escape=True,droplet=D,st=bool(st),adapt=bool(adapt),spray_cfl=.5)
            sim.seed(lambda x,y:y<14); sim.calibrate()
            pos=np.c_[sim.pos,np.full(len(sim.pos),.5)].astype(np.float32)
            typ=sim.typ.astype(np.float32); vel=np.zeros_like(pos)
            vel[(typ==0)&(pos[:,1]>16),0]=30
            vel[(typ==1)&(pos[:,1]>11),0]=10
            z=np.zeros(len(pos),np.float32); one=np.ones(len(pos),np.float32)
            mass=np.where(typ>.5,1000,1).astype(np.float32)
            initial=dict(positions=pos,particle_phase=typ,velocities=vel,mass=mass,
                momx=mass*vel[:,0],momy=z,momz=z,tau=z,wt=one,wtph=typ,
                scale=one,phase_st=typ,radius=z,spray=z)
            p=dict(ESCAPE,**ports(dict(E,rate=1,grad_threshold=2)),rho0_face=float(np.float32(sim.rho0_face)),
                   st=st,adapt=adapt)
            key=f"st{st}_adapt{adapt}"
            on=run(key+"_on",24,initial,p)
            off=run(key+"_off",24,initial,dict(p,emit_model=0))
            keep=on["spray"]<.5
            fields=("positions","velocities","mass","phase_st","scale","tau")
            same={f:bool(np.array_equal(on[f][keep],off[f])) for f in fields}
            assert all(same.values()), (key,same)
            assert int((~keep).sum())>0
            assert (on["radius"][~keep]>0).all()
            assert np.array_equal(on["mass"]==0,~keep)
            py=[]
            for enable in (False,True):
                s=Sim(nx,ny,escape=True,droplet=D,st=bool(st),adapt=bool(adapt),
                      emit=dict(E,rate=1,grad_threshold=2) if enable else None,spray_cfl=.5)
                s.seed(lambda x,y:y<14);s.calibrate();s.pos=pos[:,:2].astype(float);s.vel=vel[:,:2].astype(float)
                for _ in range(24):s.step(dt)
                py.append(s)
            assert py[1].spray.any()
            assert np.array_equal(py[0].pos,py[1].pos[~py[1].spray])
            assert np.array_equal(py[0].vel,py[1].vel[~py[1].spray])
            rec["cases"][key]=dict(primary_identical=same,spray=int((~keep).sum()),python_primary_identical=True)
            if st==0 and adapt==0:
                first=run("resume_first",12,initial,p)
                second=run("resume_second",12,first,p,12)
                resume={f:bool(np.array_equal(on[f],second[f])) for f in on}
                assert all(resume.values()),resume
                rec["resume_bit_identical"]=resume
            print(key,"PASS",rec["cases"][key],flush=True)
        rec["passed"]=True
    finally:
        (out/"receipt.json").write_text(json.dumps(rec,indent=2))


if __name__=="__main__":main()
