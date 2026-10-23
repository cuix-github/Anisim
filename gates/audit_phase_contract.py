"""Reproduce the pressure-phase contract mismatch and recheck common states."""
import argparse
import json
import os
from pathlib import Path
import sys
import numpy as np

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE),str(HERE.parent/"experiments")]
from spray_events_2d import make_sim, D, E, load_frame
from validate_spray_events import ports
import validate_spray_node as gate


def main():
    p=argparse.ArgumentParser();p.add_argument("--source",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path);a=p.parse_args()
    src,out=a.source.resolve(),a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    job=json.loads((src/"job.json").read_text());gate.DT=job["dt"];gate.OUT=str(out)
    os.environ["PFFLIP_SPRAY_AUDIT"]="1"
    rec=dict(common_states={},phase_difference={})
    for step in (0,720,1280):
        state=load_frame(src,"bifrost",step)
        s=make_sim();s.pos=state["pos"].astype(float);s.vel=state["vel"].astype(float)
        s.typ=state["typ"].astype(int);s.radius=state["rad"].astype(float)
        s.spray=state["spr"];s.escaped=state["esc"]
        s.scale=np.ones(len(s.pos));s.tau=np.zeros(len(s.pos));s.xi=np.ones(len(s.pos))
        d,ch,stats=gate.node_step(f"state_{step}",s,s.radius,s.spray,step_index=step+1,**ports(E,D))
        ref,_,_=gate.twin(s,ch,s.radius,s.spray,step+1,E,.5,droplet=D)
        good,r=gate.compare(f"common state {step}",d,ref)
        assert good and r["escape_mismatches"]==0
        rec["common_states"][str(step)]=dict(passed=True,**r)
        # Same particles: compare the two definitions before pressure solving.
        s.p2g();mass_phase=s._face_phase()
        err=np.concatenate([(mass_phase[i]-ch[k]).ravel() for i,k in enumerate(("u_phase","v_phase"))])
        rec["phase_difference"][str(step)]=dict(max=float(abs(err).max()),
            rms=float(np.sqrt(np.mean(err*err))),faces_above_point1=int((abs(err)>.1).sum()))
    (out/"receipt.json").write_text(json.dumps(rec,indent=2))
    print(json.dumps(rec["phase_difference"],indent=2))


if __name__=="__main__":main()
