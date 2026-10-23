"""Short actual-graph regression against retained pre-phase-change spray output."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'experiments'))
from spray_events_2d import run_cpp

def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);source=a.source.resolve()
    job=json.loads((source/'job.json').read_text())
    assert job['nx']==64 and job['ny']==80 and job['emission']['model']==2
    (out/'job.json').write_text(json.dumps(dict(kind='Phase mode 0 full-graph compatibility',source=str(source),
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()),indent=2))
    run_cpp(out/'bifrost',41,job['dt'],40,2)
    checked=[]
    for j in (0,40):
        for channel in ('frame','vel','ph','rad','spr','esc'):
            n=f'{channel}.{j:04d}.npy'
            assert np.array_equal(np.load(out/'bifrost'/n),np.load(source/'bifrost'/n)),n
            checked.append(n)
    (out/'receipt.json').write_text(json.dumps(dict(passed=True,bit_identical=checked),indent=2))
    print('PASS: 12 legacy full-graph arrays bit-identical')

if __name__=='__main__':main()
