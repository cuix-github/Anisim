"""Actual-graph refusal checks for the deliberately narrow volume candidate."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'gates'),str(ROOT/'graphs')]
from validate_volume_contract import make
from loop3d import run_bifrost3


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);rows={}
    s=make(.5,'density');pos=np.c_[s.pos,np.full(len(s.pos),.5)].astype(np.float32)
    for name,overrides,nz in [('negative',dict(volume_tau=-1),1),('phase',dict(phase_model=0),1),
        ('escape',dict(escape=2),1),('3d',{},2),('population',{},1),('empty_support',{},1)]:
        q=pos.copy();typ=s.typ.copy()
        if name=='population':q=q[:-1];typ=typ[:-1]
        if name=='empty_support':q[:,:2]=12
        kwargs=dict(volume_tau=.5,gravity=0.,phase_model=2,rho0_face=s.rho0_face)
        kwargs.update(overrides);d=out/name
        try:
            run_bifrost3(name,24,24,nz,lambda x,y,z:y<12,1,.0025,
                seed_fn=lambda *args:(q,typ),fused=True,extra_ports=kwargs,
                out_dir=str(d),log_prefix=str(d/'bifcmd'))
        except SystemExit as exc:
            msg=str(exc)
            assert 'refused step' in msg and 'volume contract' in msg,msg
            rows[name]=msg
        else:raise AssertionError('unsupported volume configuration accepted: '+name)
    (out/'receipt.json').write_text(json.dumps(dict(passed=True,cases=rows),indent=2))
    print('PASS: six actual graph configurations refused')


if __name__=='__main__':main()
