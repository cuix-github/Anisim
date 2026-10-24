"""True-residual certification/restart, preserving the original failed arms."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'gates'),str(ROOT/'reference2d'),str(ROOT/'experiments')]
from validate_pressure_failure import call
from transport_contract import TransportContractSim


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=('before','after'),required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--source',type=Path,required=True)
    p.add_argument('--before',type=Path);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    job=dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),stage=a.stage)
    (out/'job.json').write_text(json.dumps(job,indent=2))
    result={};os.environ['PFFLIP_RESIDUAL_RESTART_TRACE']='1'
    for name in ('flow2d_jacobi','flow3d_mg','legacy2d','adaptive2d'):
        src=a.source/name;ports=json.loads((src/'converged/ports.json').read_text())
        result[name]={}
        for tol in (1e-12,1e-13,1e-14):
            tag=f'{tol:.0e}';d=out/name/tag
            stats,_=call(src/'input',d,dict(ports,tolerance=tol),legacy=name=='legacy2d')
            stats['restart_messages']=(d/'stderr.log').read_text().count('residual-restart:')
            result[name][tag]=stats
            if a.stage=='after' and stats['iterations_used']>=0:
                assert stats['final_residual']<=float(np.float32(tol)),(name,tag,stats)
        (out/'receipt.json').write_text(json.dumps(result,indent=2))
    # Reproduce the observed failure from the same seed, without a replay or
    # injected numerical error. Save the exact grid at the first refused step.
    s=TransportContractSim(112,80,phase_model=2,transport_model=3,
        pressure_tolerance=float(np.float32(1e-10)),g=float(np.float32(-9.8)),
        alpha_l=float(np.float32(.97)),alpha_g=float(np.float32(.9)))
    s.seed(lambda x,y:(y<18)|((abs(x-56)<7)&(y>39)&(y<64)),seed=1)
    s.calibrate();s.rho0_face=float(np.float32(s.rho0_face));s.pos=s.pos.astype(np.float32).astype(float)
    project=s.project;restarts=[];j=0
    def audited(dt):
        try:project(dt)
        except RuntimeError:
            np.savez_compressed(out/'refused_grid.npz',u=s.u,v=s.v,phi_u=s._face_phase()[0],phi_v=s._face_phase()[1])
            raise
        if getattr(s,'pressure_restarts',0):restarts.append(dict(step=j,count=s.pressure_restarts,residual=s.pressure_residual))
    s.project=audited
    failure=None
    for j in range(400):
        try:s.step(float(np.float32(.0025)))
        except RuntimeError as e:failure=dict(step=j,message=str(e),residual=s.pressure_residual);break
        for n in ('pos','vel','radius','scale','tau'):setattr(s,n,getattr(s,n).astype(np.float32).astype(float))
    result['observed_python']=dict(steps_completed=j if failure else 400,failure=failure,restarts=restarts)
    if a.stage=='before':assert failure is not None
    else:
        assert failure is None and restarts,result['observed_python']
        prior=json.loads((a.before/'receipt.json').read_text())
        recovered=[]
        for name in ('flow2d_jacobi','flow3d_mg','legacy2d','adaptive2d'):
            for tag,stats in result[name].items():
                if prior[name][tag]['iterations_used']<0 and stats['iterations_used']>=0:
                    recovered.append(dict(case=name,tolerance=tag,restarts=stats['restart_messages']))
        assert recovered,'No actual C++ failure recovered; do not claim C++ coverage'
        result['cpp_recovered']=recovered
    (out/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
