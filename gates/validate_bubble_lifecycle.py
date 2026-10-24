"""Generic isolated-air inclusions: actual graph, Python, legacy and restart.

Before/after modes explicitly distinguish recorded old behavior from the paper
expectation. No wave geometry, damping or emission threshold is tuned here.
"""
import argparse,json,hashlib,os,subprocess,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'reference2d'),str(ROOT/'gates'),str(ROOT/'graphs')]
from phase_contract import PhaseContractSim
from validate_transport_contract import arrays
from loop3d import run_bifrost3,load_final_state


def make(rule,drop):
    s=PhaseContractSim(16,12,g=0,escape=rule,esc_phi=float(np.float32(.3)),
        droplet=dict(r_sigma=0,evap_k=0) if drop else None)
    s.seed(lambda x,y:np.ones(len(x),bool),jitter=0);s.calibrate()
    s.pos=np.maximum(.51,np.minimum([15.49,11.49],s.pos)).astype(np.float32).astype(float)
    s.rho0_face=float(np.float32(s.rho0_face))
    ids=[int(np.argmin(np.sum((s.pos-q)**2,axis=1))) for q in ([5.25,4.25],[10.75,7.75])]
    s.typ[ids]=0;s.p2g()
    assert (s.escape_phase_samples()[ids]>.7).all()
    return s


def run(s,d,steps=4,state=None,offset=0):
    _,_,diag=run_bifrost3('bubble_lifecycle',16,12,1,lambda x,y,z:y<12,steps,float(np.float32(.0025)),
        init_arrays=arrays(s) if state is None else state,fused=True,final_state=True,
        dumps=('ph','vel','esc'),dump_every=1,step_offset=offset,
        extra_ports=dict(gravity=0,phase_model=2,rho0_face=s.rho0_face,escape=s.escape,esc_phi=s.esc_phi,
            droplet_model=1 if s.droplet else 0,droplet_r_sigma=0,droplet_evap_k=0),
        out_dir=str(d),log_prefix=str(d/'bifcmd'))
    assert np.isfinite(diag).all() and (diag[:,1]<=float(np.float32(1e-6))).all()
    return load_final_state(str(d))


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--before',type=Path);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);os.environ['OMP_NUM_THREADS']='1'
    rec=dict(passed=False,mode='paper_after' if a.before else 'baseline_old_behavior',cases={})
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    (out/'job.json').write_text(json.dumps(dict(scope=__doc__,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),mode=rec['mode'],
        expected_new='two escaped air inclusions removed only for escape=2; all 766 water particles retained'),indent=2))
    try:
        for name,rule,drop in [('off',0,False),('legacy',1,False),('phase_plain',2,False),('phase_droplet',2,True)]:
            s=make(rule,drop);initial=arrays(s);actual=run(s,out/name)
            for _ in range(4):
                s.step(float(np.float32(.0025)))
                s.pos=s.pos.astype(np.float32).astype(float);s.vel=s.vel.astype(np.float32).astype(float)
            expected=766 if a.before and rule==2 else 768
            assert len(actual['positions'])==len(s.pos)==expected,(name,len(actual['positions']),len(s.pos),expected)
            dp=float(abs(actual['positions'][:,:2]-s.pos).max());dv=float(abs(actual['velocities'][:,:2]-s.vel).max())
            assert dp<8e-6 and dv<5e-5,(name,dp,dv)
            assert np.array_equal(actual['particle_phase'],s.typ)
            assert sum(actual['particle_phase']>.5)==766
            assert all(len(v)==expected for v in actual.values())
            row=dict(remaining=expected,water=766,air=int(sum(actual['particle_phase']<.5)),position_max=dp,velocity_max=dv)
            if a.before and rule!=2:
                old=load_final_state(str(a.before/name));checks={k:bool(np.array_equal(old[k],actual[k])) for k in old}
                assert all(checks.values()),checks;row['legacy_identical']=checks
            if a.before and rule==2:
                assert np.array_equal(actual['positions'],initial['positions'][initial['particle_phase']>.5])
                assert np.all(actual['mass']==1000.) and np.all(actual['radius']==0)
            rec['cases'][name]=row
        if a.before:
            s=make(2,True);first=run(s,out/'resume_first',2)
            second=run(s,out/'resume_second',2,first,2)
            full=load_final_state(str(out/'phase_droplet'))
            checks={k:bool(np.array_equal(full[k],second[k])) for k in full}
            assert all(checks.values()),checks;rec['resume_identical']=checks
        rec['passed']=True
    except BaseException as e:rec['failure']=repr(e);raise
    finally:
        (out/'receipt.json').write_text(json.dumps(rec,indent=2));print(json.dumps(rec,indent=2),flush=True)


if __name__=='__main__':main()
