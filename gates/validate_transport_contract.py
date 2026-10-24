"""Paper-motivated transport: analytic and actual graph/Python local checks."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'experiments'),str(ROOT/'reference2d'),str(ROOT/'graphs')]
from phase_contract_2d import make
from transport_contract import TransportContractSim
from transport_contract_2d import analytic
from loop3d import run_bifrost3


def arrays(s):
    n=len(s.pos);z=np.zeros(n,np.float32);one=np.ones(n,np.float32);phase=s.typ.astype(np.float32)
    mass=s.masses().astype(np.float32)
    mom=(mass[:,None]*s.vel).astype(np.float32) if s._port_momentum is None else s._port_momentum.astype(np.float32)
    return dict(positions=np.c_[s.pos,np.full(n,.5)].astype(np.float32),velocities=np.c_[s.vel,z].astype(np.float32),
        mass=mass,momx=mom[:,0],momy=mom[:,1],momz=z,particle_phase=phase,phase_st=phase,
        tau=z,wt=one,wtph=phase,scale=one,radius=z,spray=z)


def cpp(s,out,steps=3,dt=float(np.float32(.0025)),every=1,extra=None):
    kw=dict(phase_model=2,transport_model=s.transport_model,phase_alpha=1.,rho0_face=s.rho0_face,
        tolerance=s.pressure_tolerance,alpha_liquid=s.alpha[1],alpha_air=s.alpha[0],escape=0,st=0,adapt=0,emit_model=0)
    kw.update(extra or {})
    return run_bifrost3('transport',s.nx,s.ny,1,lambda x,y,z:y<1,steps,dt,fused=True,
        init_arrays=arrays(s),extra_ports=kw,final_state=True,dumps=('ph','vel'),dump_every=every,
        out_dir=str(out),log_prefix=str(out/'bifcmd'))


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'job.json').write_text(json.dumps(dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()),indent=2))
    bat=out/'compile.bat'
    bat.write_text('@echo off\ncall "C:/Program Files/Microsoft Visual Studio/2022/Community/VC/Auxiliary/Build/vcvars64.bat" -vcvars_ver=14.38 >nul\n'
        +f'cl /nologo /EHsc /O2 /std:c++17 "{ROOT / "gates/transport_math.cpp"}" /Fe:"{out / "math.exe"}" /Fo:"{out / "math.obj"}"\nexit /b %errorlevel%\n')
    r=subprocess.run([str(bat)],capture_output=True,text=True)
    (out/'compile.log').write_text(r.stdout+r.stderr);assert r.returncode==0,r.stdout+r.stderr
    math=json.loads(subprocess.check_output([str(out/'math.exe')],text=True))
    ref=analytic();assert np.allclose(math['rotation_errors'],ref['rotation_errors'],rtol=1e-12)
    receipt=dict(cpp_math=math,python_math=ref,local={})
    for case in ('initial','impact','rest'):
        for mode in (1,2,3):
            s=make('rest' if case=='rest' else 'plunge',2,float(np.float32(1e-10)))
            s.__class__=TransportContractSim;s.transport_model=mode
            steps=3
            if case=='impact':
                z=np.load(ROOT/'out/runs/phase_causality_20261023/control_320/energy_fields.npz')
                s.pos=z['pos'].copy();s.vel=z['old_velocity'].copy();s.typ=z['typ'].copy()
                s._port_momentum=(s.masses()[:,None]*s.vel).astype(np.float32).astype(float);steps=1
            d=out/f'{case}_{mode}';cpp(s,d,steps)
            for _ in range(steps):
                s.step(float(np.float32(.0025)))
                s.pos=s.pos.astype(np.float32).astype(float);s.vel=s.vel.astype(np.float32).astype(float)
            dp=float(abs(np.load(d/'final.npy')[:,:2]-s.pos).max())
            dv=float(abs(np.load(d/'final_velocities.npy')[:,:2]-s.vel).max())
            assert dp<8e-6 and dv<5e-5,(case,mode,dp,dv)
            receipt['local'][f'{case}_{mode}']=dict(position_max=dp,velocity_max=dv)
            (out/'receipt.json').write_text(json.dumps(receipt,indent=2))
    for name,kw in [('bad_mode',dict(transport_model=4)),('phase',dict(phase_model=1)),('escape',dict(escape=1))]:
        s=make('rest',2,float(np.float32(1e-10)));s.__class__=TransportContractSim;s.transport_model=3
        d=out/('refuse_'+name)
        try:cpp(s,d,1,extra=kw)
        except SystemExit as e:assert 'unsupported experimental transport' in str(e),str(e)
        else:raise AssertionError(name+' not refused')
    receipt['refusals']=3
    (out/'receipt.json').write_text(json.dumps(receipt,indent=2));print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
