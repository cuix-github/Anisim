"""Independent Eq.7 / pure-face gates through the actual fused Bifrost P2G."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'gates'),str(ROOT/'reference2d'),str(ROOT/'graphs')]
from phase_contract import PhaseContractSim, compress, pure_faces
from validate_p2g import CHANS, make_arrays
from loop3d import run_bifrost3


def frozen():
    dirty=subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    assert not dirty,dirty
    return subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()


def run_p2g(out,s,model):
    out.mkdir(parents=True,exist_ok=False)
    arrays=make_arrays(np.c_[s.pos,np.full(len(s.pos),.5)],np.c_[s.vel,np.zeros(len(s.pos))],s.masses(),s.typ)
    args=[str(ROOT/'graphs/run_p2g_test.bat')]
    kw=dict(nx=s.nx,ny=s.ny,nz=1,phase_model=model,phase_alpha=s.phase_alpha,
            rho0_face=s.rho0_face,rho_liquid=s.rho_l,rho_air=s.rho_g)
    for k,v in kw.items(): args+=['--set-port',k,str(v)]
    for k,v in arrays.items():
        f=out/(k+'_in.npy');np.save(f,v);args+=['--set-port','path_'+k,str(f)]
    for k in CHANS: args+=['--set-port','path_'+k,str(out/(k+'.npy'))]
    r=subprocess.run(args,capture_output=True,text=True)
    (out/'stdout.log').write_text(r.stdout);(out/'stderr.log').write_text(r.stderr)
    assert r.returncode==0,(r.stdout[-1000:],r.stderr[-1000:])
    return {k:np.load(out/(k+'.npy')) for k in CHANS}


def sim(nx=24,ny=32,model=2,case='mixed'):
    s=PhaseContractSim(nx,ny,phase_model=model,g=float(np.float32(-9.8)),
                      alpha_l=float(np.float32(.97)),alpha_g=float(np.float32(.9)))
    s.seed(lambda x,y: np.zeros_like(x,dtype=bool) if case=='air' else
           np.ones_like(x,dtype=bool) if case=='water' else y<ny*.45,seed=137)
    s.calibrate();s.rho0_face=float(np.float32(s.rho0_face));s.pos=s.pos.astype(np.float32).astype(float)
    return s


def algebra():
    rho0=2.6796875;floor=np.log(1000)*rho0
    raw=floor+1000*rho0*np.array([-.1,0,.0625,.25,1,4])
    assert np.allclose(compress(raw,rho0),[0,0,.25,.5,1,1],atol=1e-15)
    assert np.allclose(compress(raw,rho0,alpha=4),[0,0,.125,.25,.5,1],atol=1e-15)
    # A [air,air,water,water] strip must retain the mixed face, not blur it.
    cell=np.array([[0],[0],[1],[1.]])
    u=np.full((5,1),.37);v=np.full((4,2),.63)
    a,b=pure_faces((u,v),cell)
    assert np.array_equal(a[:,0],[.37,0,.37,1,.37]) and np.array_equal(b,v)
    for kw in [dict(phase_alpha=0),dict(phase_model=3),dict(st=True),dict(adapt=True)]:
        try: PhaseContractSim(8,8,**kw)
        except ValueError: pass
        else: raise AssertionError(kw)
    return dict(eq7=True,pure_face_truth_table=True,invalid_modes_refused=True)


def main():
    a=argparse.ArgumentParser();a.add_argument('--out',type=Path,required=True);arg=a.parse_args()
    out=arg.out.resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'job.json').write_text(json.dumps(dict(commit=frozen(),kind='phase contract static and one-step gates'),indent=2))
    receipt=dict(algebra=algebra(),fields={},steps={})
    for model in (1,2):
        for case in ('air','water','mixed'):
            s=sim(model=model,case=case)
            s.vel=np.random.default_rng(37).normal(0,2,s.vel.shape).astype(np.float32).astype(float)
            s.p2g();actual=run_p2g(out/f'field_{model}_{case}',s,model)
            err=max(float(abs(actual[n+'_phase'].reshape(p.shape)-p).max()) for n,p in zip(('u','v'),s._face_phase()))
            assert err<1e-7,(model,case,err)
            if case=='air': assert all((p==0).all() for p in s._face_phase())
            if case=='water' and model==2:
                assert all((p[2:-2,2:-2]==1).all() for p in s._face_phase())
            receipt['fields'][f'{model}_{case}']=dict(max_phase_error=err)
        s=sim(model=model);p3=np.c_[s.pos,np.full(len(s.pos),.5)].astype(np.float32)
        dt=float(np.float32(.0025));d=out/f'loop_{model}'
        run_bifrost3('phase_one',s.nx,s.ny,1,lambda x,y,z:y<14.4,1,dt,
            seed_fn=lambda *a:(p3,s.typ.copy()),fused=True,extra_ports=dict(phase_model=model,
            phase_alpha=1.,rho0_face=s.rho0_face,escape=0),final_state=True,out_dir=str(d),
            log_prefix=str(d/'bifcmd'))
        s.step(dt)
        pos=np.load(d/'final.npy')[:,:2];vel=np.load(d/'final_velocities.npy')[:,:2]
        dp=float(abs(pos-s.pos).max());dv=float(abs(vel-s.vel).max())
        assert dp<4e-6 and dv<2e-5,(dp,dv)
        receipt['steps'][str(model)]=dict(position_max=dp,velocity_max=dv)
    # A graph may not silently apply the uniform calibration to adaptive/ST states.
    for name,kw in [('st',dict(st=1)),('adapt',dict(adapt=1)),('alpha',dict(phase_alpha=0))]:
        s=sim();p3=np.c_[s.pos,np.full(len(s.pos),.5)].astype(np.float32);d=out/('refuse_'+name)
        try:
            run_bifrost3('refuse_phase',s.nx,s.ny,1,lambda x,y,z:y<14.4,1,.0025,
                seed_fn=lambda *a:(p3,s.typ.copy()),fused=True,
                extra_ports=dict(dict(phase_model=2,rho0_face=s.rho0_face),**kw),
                out_dir=str(d),log_prefix=str(d/'bifcmd'))
        except SystemExit:
            assert 'phase contract' in (d/'bifcmd.stderr.log').read_text()
        else: raise AssertionError('unsupported configuration accepted: '+name)
    receipt['graph_rejections']=3
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
