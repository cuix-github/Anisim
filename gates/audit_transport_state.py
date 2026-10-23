"""Read-only common-state transport intervention and zero-flow boundary witness."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference2d'),str(ROOT/'experiments')]
from phase_contract_2d import make,measures
from transport_contract import TransportContractSim


def main():
    p=argparse.ArgumentParser();p.add_argument('--state',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    state=np.load(a.state);moved=[];result=dict(source=str(a.state.resolve()),transport={})
    (out/'job.json').write_text(json.dumps(dict(commit=subprocess.check_output(
        ['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()),indent=2))
    for mode in range(4):
        s=make('plunge',2,float(np.float32(1e-10)));s.__class__=TransportContractSim;s.transport_model=mode
        s.pos=state['pos'].copy();s.vel=state['g2p_velocity'].copy();s.typ=state['typ'].copy()
        s.u=state['post_u'].copy();s.v=state['post_v'].copy()
        s._contract_faces=(state['phi_u'].copy(),state['phi_v'].copy())
        s.advect(float(np.float32(.0025)));moved.append(s.pos.copy())
        delta=np.linalg.norm(moved[-1]-moved[0],axis=1)
        result['transport'][str(mode)]=dict(water_max_cells=float(delta[s.typ==1].max()),
            water_mean_cells=float(delta[s.typ==1].mean()),air_max_cells=float(delta[s.typ==0].max()))
    assert result['transport']['2']['air_max_cells']==0
    # Legacy seeding includes particles inside its later .51-cell clamp band.
    # Record the mismatch; do not declare it the cause of impact instability.
    s=make('plunge',2,1e-10);initial=s.pos.copy();energy=measures(s.pos,s.vel,s.typ)['energy']
    s.advect(.0025);delta=s.pos-initial;changed=np.any(delta!=0,axis=1)
    result['zero_flow_boundary_witness']=dict(moved=int(changed.sum()),water_moved=int((changed&(s.typ==1)).sum()),
        displacement_max_cells=float(abs(delta).max()),
        energy_fraction=measures(s.pos,s.vel,s.typ)['energy']/energy-1,
        interpretation='Seed-versus-clamp support mismatch; not proof of the cause of late spray.')
    (out/'receipt.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))


if __name__=='__main__':main()
