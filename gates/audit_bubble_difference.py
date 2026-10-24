"""Check old/new actual C++ steps differ only by classified-air removal.

Uses retained common-input impact checkpoints. No simulation is run.
"""
import argparse,json,subprocess
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    a.out.mkdir(parents=True,exist_ok=False);rows=[];rec=dict(passed=False,cases=rows)
    before=ROOT/'out/runs/phase_repair_20261023/escape_phase_audit'
    after=ROOT/'out/runs/acceptance_closure_20261024/impact_states_after'
    try:
        for step in (1280,3200,7200):
            old=before/f'step_{step}/bifrost';new=after/f'step_{step}/bifrost'
            inputs=['positions','velocities','mass','momx','momy','momz','particle_phase','phase_st','tau','wt','wtph','scale','radius','spray']
            same_inputs={n:bool(np.array_equal(np.load(old/f'{n}.npy'),np.load(new/f'{n}.npy'))) for n in inputs}
            assert all(same_inputs.values()),same_inputs
            removed=(np.load(old/'final_phase.npy')<.5)&(np.load(old/'esc.0000.npy')>.5)
            assert removed.any()
            names=['final.npy']+[p.name for p in old.glob('final_*.npy')]
            same={n:bool(np.array_equal(np.load(old/n)[~removed],np.load(new/n))) for n in names}
            row=dict(step=step,deleted_air=int(removed.sum()),identical_inputs=same_inputs,
                unchanged_survivor_arrays=same);rows.append(row)
            assert len(same)==13 and all(same.values()),row
        rec['passed']=True
    except BaseException as e:rec['failure']=repr(e);raise
    finally:
        (a.out/'receipt.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec,indent=2))


if __name__=='__main__':main()
