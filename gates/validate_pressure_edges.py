"""Zero-budget, invalid-config and actual multi-step refusal checks."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'gates'), str(ROOT/'graphs')]
from validate_pressure_failure import call
from loop3d import run_bifrost3


def main():
    p = argparse.ArgumentParser(); p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True); a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve(); out.mkdir(parents=True,exist_ok=False); source=a.source.resolve()
    (out/'job.json').write_text(json.dumps(dict(commit=subprocess.check_output(
        ['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()),indent=2))
    result={}
    src=source/'zero_rhs'
    kw=json.loads((src/'capped/ports.json').read_text())
    for solver in (1,2):
        st,_=call(src/'input',out/f'zero_adaptive_{solver}',dict(kw,pressure_solver=solver))
        assert st['iterations_used']==0 and st['final_residual']==0 and st['particle_count']>0,st
        result[f'zero_adaptive_{solver}']=st
    src=source/'flow2d_jacobi'; kw=json.loads((src/'converged/ports.json').read_text())
    for name,delta in [('zero_tol',dict(tolerance=0)),('negative_cap',dict(max_iterations=-1)),
                       ('zero_cap_nonzero_rhs',dict(max_iterations=0))]:
        st,_=call(src/'input',out/name,dict(kw,**delta))
        assert st['iterations_used']==-2 and st['particle_count']==0,st
        result[name]=dict(refused=True,particle_count=0)
    loop=out/'loop'
    try:
        run_bifrost3('pressure_failure',16,16,1,lambda x,y,z:(x<8)&(y<12),3,.01,
            fused=True,dump_every=3,final_state=True,out_dir=str(loop),
            extra_ports=dict(max_iterations=1),log_prefix=str(loop/'bifcmd'))
    except SystemExit as e:
        assert 'pressure' in str(e) and 'refused step 0' in str(e),str(e)
        result['loop_refusal']=str(e)
    else: raise AssertionError('failed chunk returned successfully')
    result['negative_step_diagnostics']=[int(np.load(q).ravel()[0]) for q in sorted(loop.glob('diag.*.npy'))]
    (out/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
