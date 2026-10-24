"""One frozen diagnostic after the interface closure's positive-support refusal.

Distinguish pre-existing sampling holes from detached-particle exclusion. Then
disable only density feedback in the three common-input execution checks. This
does not qualify the failed integrated candidate or relax any numerical gate.
"""
import hashlib,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np
from validate_interface_closure import ROOT,make,cpp,py_step,compare,write

def state(step):
    d=ROOT/f'out/runs/acceptance_closure_20261024/impact_states_after/step_{step}/bifrost'
    s=make(64,80,transport=3,g=-9.8)
    for attr,file in [('pos','positions'),('vel','velocities'),('typ','particle_phase'),
                      ('radius','radius'),('spray','spray'),('scale','scale')]:
        a=np.load(d/f'{file}.npy')
        if attr in ('pos','vel'):a=a[:,:2]
        if attr=='spray':a=a>.5
        if attr=='typ':a=a.astype(int)
        setattr(s,attr,a.astype(float) if attr in ('pos','vel','radius','scale') else a)
    s.tau=np.zeros(len(s.pos));s.xi=np.zeros((len(s.pos),2));s.escaped=np.zeros(len(s.pos),bool)
    return s

def main():
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=Path(sys.argv[1]).resolve();out.mkdir(parents=True,exist_ok=False)
    os.environ['OMP_NUM_THREADS']='1'
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    write(out/'job.json',dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),command=sys.argv,scope=__doc__))
    r=dict(passed=False,support=[],common=[]);start=time.perf_counter()
    try:
        for step in (1280,3200,7200):
            s=state(step);q=s.pos-.5;b=np.floor(q+.5).astype(np.int64)
            all_primary=s._splat((64,80),(b,q-b),[(~s.spray).astype(float)])[0]
            active=s.number_density();s.p2g();zero=active==0;pre=all_primary==0
            r['support'].append(dict(step=step,primary=int((~s.spray).sum()),
                detached=int(((s.typ==1)&(s.radius>0)&~s.spray).sum()),
                zero_before_exclusion=int(pre.sum()),zero_after_exclusion=int(zero.sum()),
                newly_empty=int((zero&~pre).sum()),empty_water_cells=int((zero&(s.cell_phase>=.5)).sum()),
                zero_coordinates=np.argwhere(zero).tolist(),
                active_positive_min=float(active[active>0].min())))
        write(out/'receipt.json',r)
        for step in (1280,3200,7200):
            s=state(step);water=int(((s.typ==1)&~s.spray).sum())
            a,res=cpp(s,out/f'common_tau0_{step}',1);py_step(s)
            row=dict(step=step,**compare(s,a),residual=res)
            assert int(((s.typ==1)&~s.spray).sum())==water
            row['water']=water;r['common'].append(row);write(out/'receipt.json',r)
        r['passed']=True
    except BaseException as e:r['failure']=repr(e);raise
    finally:
        r['seconds']=time.perf_counter()-start;write(out/'receipt.json',r)
        print(json.dumps(r),flush=True)

if __name__=='__main__':main()
