"""One bounded diagnostic of the failed area measurement, no time integration.

Same analytic circle and physical domain at h=1,.5,.25. The original h=1
calibration stays failed; refinement must not retroactively pass that gate.
"""
import argparse,json,sys,subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference2d'),str(ROOT/'gates')]
from phase_contract import PhaseContractSim
from phase_area import area,particle_field


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    a.out.mkdir(parents=True,exist_ok=False);rows=[];rec=dict(passed=False,rows=rows)
    (a.out/'job.json').write_text(json.dumps(dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        scope=__doc__,criteria='strictly decreasing relative errors; finest <2%; no correction of coarse data'),indent=2))
    try:
        for scale in (1,2,4):
            s=PhaseContractSim(48*scale,40*scale)
            s.seed(lambda x,y:(x/scale-24)**2+(y/scale-20)**2<8**2,jitter=0);s.calibrate()
            s.rho0_face=float(np.float32(s.rho0_face))
            s.pos=np.maximum(.51,np.minimum([s.nx-.51,s.ny-.51],s.pos)).astype(np.float32).astype(float)
            g=area(particle_field(s));measured=g['area']/scale**2;expected=np.pi*8**2
            rows.append(dict(cell_size=1/scale,measured=measured,expected=expected,
                relative_error=abs(measured/expected-1),diagonal_difference=g['diagonal_difference']/scale**2))
        errors=[q['relative_error'] for q in rows]
        assert errors[2]<errors[1]<errors[0] and errors[2]<.02,errors
        rec['passed']=True
    except BaseException as e:rec['failure']=repr(e);raise
    finally:
        (a.out/'receipt.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec,indent=2))


if __name__=='__main__':main()
