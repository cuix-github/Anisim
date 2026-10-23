"""Before/after failure-contract gate on unrelated 2D/3D pressure systems.

Capture --stage before with the unmodified binary, then --stage after in a
fresh directory. Both successful outputs and deliberately capped solves stay.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'gates'),str(ROOT/'graphs')]
from validate_adaptive_3d import package
from audit_pressure_replay import read_dump,operator

CASES=[('rest2d',24,24,1,0,1,False,False),
       ('flow2d_jacobi',24,20,1,0,0,False,False),
       ('flow3d_mg',16,16,12,0,1,False,False),
       ('flow3d_jacobi',16,16,12,0,0,False,False),
       ('cut3d',16,16,12,0,1,True,False),
       ('adaptive2d',32,32,1,1,1,False,False),
       ('adaptive3d_cleanup',24,24,16,2,1,False,False),
       ('zero_rhs',12,12,1,0,1,False,True),
       ('legacy2d',24,20,1,0,0,False,False)]


def call(indir,out,ports,legacy=False,solid=False):
    out.mkdir(parents=True,exist_ok=False)
    cmd=[str(ROOT/'graphs'/('run_solve_step.bat' if legacy else 'run_solve_step3.bat'))]
    for n,v in ports.items():
        if legacy and n in ('nz','preconditioner','pressure_solver','adaptive_band','adaptive_levels'):continue
        cmd+=['--set-port',n,str(v)]
    names=[g+n for g in ('uv' if legacy else 'uvw') for n in ('_mass','_mom','_phase')]
    names+=['positions','velocities','particle_phase']+(['solid_phi'] if solid else [])
    for n in names:cmd+=['--set-port','path_'+n,str(indir/(n+'.npy')).replace('\\','/')]
    outputs=['out_u','out_v','pressure','out_positions','out_velocities']+([] if legacy else ['out_w','out_escaped'])
    for n in outputs:cmd+=['--set-port','path_'+n,str(out/(n+'.npy')).replace('\\','/')]
    env=os.environ.copy()
    if not legacy and ports['pressure_solver']==0:env['PFFLIP_PRESSURE_DUMP']=str(out/'pressure_double.bin')
    proc=subprocess.run(cmd,capture_output=True,text=True,env=env)
    (out/'stdout.log').write_text(proc.stdout);(out/'stderr.log').write_text(proc.stderr)
    assert proc.returncode==0,(out,proc.stderr[-1000:])
    stats={}
    for line in proc.stdout.splitlines():
        if "Output port '" in line and 'ok_' not in line:
            stats[line.split("'")[1]]=float(line.rsplit(')',1)[-1].strip())
    assert 'iterations_used' in stats
    if (out/'pressure_double.bin').exists():
        coeff,b,p,_=read_dump(out/'pressure_double.bin')
        stats['independent_true_residual']=float(np.linalg.norm(b-operator(p,coeff))/max(np.linalg.norm(b),1e-300))
    stats['particle_count']=len(np.load(out/'out_positions.npy'))
    (out/'stats.json').write_text(json.dumps(stats,indent=2));(out/'ports.json').write_text(json.dumps(ports,indent=2))
    return stats,outputs


def main():
    p=argparse.ArgumentParser();p.add_argument('--stage',choices=('before','after'),required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--before',type=Path)
    a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    (out/'job.json').write_text(json.dumps(dict(kind='Generic pressure failure contract',stage=a.stage,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256((ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll').read_bytes()).hexdigest()),indent=2))
    result={}
    for name,nx,ny,nz,solver,pre,solid,zero in CASES:
        inp=out/name/'input';inp.mkdir(parents=True)
        mask=lambda x,y,z:(y<ny*.28)|((x<nx*.4)&(y<ny*.7))
        def field(x,y,z):
            if zero or name=='rest2d':return x*0,y*0,z*0
            return (2*np.sin(2*np.pi*x/nx)*np.cos(np.pi*y/ny),
                    3*np.sin(np.pi*y/ny)*np.cos(np.pi*x/nx),
                    0*z if nz==1 else np.sin(2*np.pi*z/nz)*np.cos(np.pi*x/nx))
        package(str(inp),nx,ny,nz,mask,field,seed=37)
        if solid:
            x,y,z=np.meshgrid(np.arange(nx+1),np.arange(ny+1),np.arange(nz+1),indexing='ij')
            np.save(inp/'solid_phi.npy',(y-3-.13*x-.07*z).astype(np.float32).ravel())
        kw=dict(nx=nx,ny=ny,nz=nz,dt=.01,gravity=0 if zero else -9.8,
            rho_liquid=1000,rho_air=1,alpha_liquid=.97,alpha_air=.9,
            pressure_solver=solver,preconditioner=pre,adaptive_band=2,adaptive_levels=2,tolerance=1e-6)
        if solid:kw['solid_mode']=2
        case={}
        for arm,cap in [('converged',4000),('capped',0 if zero else 1)]:
            case[arm],outputs=call(inp,out/name/arm,dict(kw,max_iterations=cap),name=='legacy2d',solid)
            st=case[arm]
            if a.stage=='after':
                if arm=='converged' or zero:
                    assert st['iterations_used']>=0 and st['final_residual']<=1.01e-6,(name,arm,st)
                    for n in outputs:
                        assert np.array_equal(np.load(out/name/arm/(n+'.npy')),np.load(a.before/name/arm/(n+'.npy'))),(name,arm,n)
                    case[arm]['bit_identical_output_arrays']=len(outputs)
                else:
                    assert st['iterations_used']==-2 and st['particle_count']==0,(name,st)
                    assert st['final_residual']>1e-6 or not np.isfinite(st['final_residual'])
                if 'independent_true_residual' in st:
                    assert np.isclose(st['final_residual'],st['independent_true_residual'],rtol=1e-4,atol=1e-11),(name,st)
        if name=='flow2d_jacobi':
            cap=int(case['converged']['iterations_used'])
            case['last_allowed_iteration'],outputs=call(inp,out/name/'last_allowed_iteration',dict(kw,max_iterations=cap))
            if a.stage=='after':
                st=case['last_allowed_iteration'];assert st['iterations_used']>=0 and st['final_residual']<=1.01e-6
                for n in outputs:assert np.array_equal(np.load(out/name/'last_allowed_iteration'/(n+'.npy')),np.load(out/name/'converged'/(n+'.npy')))
        result[name]=case;(out/'receipt.json').write_text(json.dumps(result,indent=2))
        print(name,json.dumps(case),flush=True)


if __name__=='__main__':main()
