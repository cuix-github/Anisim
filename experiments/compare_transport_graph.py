"""Independent full trajectories using the same unchanged similarity thresholds."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'experiments'),str(ROOT/'reference2d'),str(ROOT/'gates')]
from phase_contract_2d import make,measures
from transport_contract import TransportContractSim
from validate_transport_contract import cpp


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--control',type=Path)
    a=p.parse_args();assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    source=a.source.resolve();job=json.loads((source/'job.json').read_text())
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    s=make('rest' if job['case']=='rest' else 'plunge',2,float(np.float32(1e-10)))
    s.__class__=TransportContractSim;s.transport_model=job['mode']
    if job['case']=='wide':
        s=TransportContractSim(112,80,phase_model=2,transport_model=job['mode'],pressure_tolerance=float(np.float32(1e-10)),
            g=float(np.float32(-9.8)),alpha_l=float(np.float32(.97)),alpha_g=float(np.float32(.9)))
        s.seed(lambda x,y:(y<18)|((abs(x-56)<7)&(y>39)&(y<64)),seed=1)
        s.calibrate();s.rho0_face=float(np.float32(s.rho0_face));s.pos=s.pos.astype(np.float32).astype(float)
    initial=measures(s.pos,s.vel,s.typ)
    (out/'job.json').write_text(json.dumps(dict(source=str(source),commit=subprocess.check_output(
        ['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),**{k:job[k] for k in ('mode','case','steps','dt')},
        binary_sha256=hashlib.sha256((ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll').read_bytes()).hexdigest(),
        thresholds=dict(tv=.05,com_cells=.2,height_p99_cells=1)),indent=2))
    d=out/'bifrost';cpp(s,d,job['steps'],job['dt'],40)
    stats=[];rows=[]
    def compare(j,py,cp,cv,ct):
        mask=ct==1;pl=py['pos'][py['typ']==1];cl=cp[mask]
        hist=lambda p:np.histogram2d(p[:,0],p[:,1],bins=(s.nx//4,s.ny//4),range=((0,s.nx),(0,s.ny)))[0]
        a,b=hist(pl),hist(cl)
        rows.append(dict(step=j,tv=float(.5*abs(a/a.sum()-b/b.sum()).sum()),
            com_cells=float(np.linalg.norm(pl.mean(0)-cl.mean(0))),
            height_p99_cells=float(abs(np.quantile(pl[:,1],.99)-np.quantile(cl[:,1],.99)))))
        stats.append(dict(step=j,time=(j+1)*job['dt'],**measures(cp,cv,ct)))
    for j in range(0,job['steps'],40):
        with np.load(source/f'frame.{j:04d}.npz') as py:
            compare(j,py,np.load(d/f'frame.{j:04d}.npy')[:,:2],np.load(d/f'vel.{j:04d}.npy')[:,:2],
                (np.load(d/f'ph.{j:04d}.npy')>.5).astype(int))
    with np.load(source/'final.npz') as py:
        compare(job['steps']-1,py,np.load(d/'final.npy')[:,:2],np.load(d/'final_velocities.npy')[:,:2],
            (np.load(d/'final_phase.npy')>.5).astype(int))
    worst={k:max(row[k] for row in rows) for k in ('tv','com_cells','height_p99_cells')}
    passed=worst['tv']<=.05 and worst['com_cells']<=.2 and worst['height_p99_cells']<=1
    checked=[]
    if a.control:
        for path in sorted(a.control.glob('final*.npy')):
            assert np.array_equal(np.load(path),np.load(d/path.name)),path.name
            checked.append(path.name)
    result=dict(similarity_pass=passed,worst=worst,bit_identical_controls=checked,
        max_sampled_energy_ratio=max([initial['energy']]+[r['energy'] for r in stats])/initial['energy'],
        final=stats[-1],comparison=rows)
    (out/'metrics.json').write_text(json.dumps(stats,indent=2));(out/'receipt.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='comparison'},indent=2))
    if not passed:raise SystemExit('Independent similarity failed; outputs retained')


if __name__=='__main__':main()
