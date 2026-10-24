"""Primary represented area must survive sub-cutoff radii when evaporation is off.

Prescribed pure-air carrier with one water primary and one air sentinel. Existing
and newly drawn sub-cutoff radii are both tested. Legacy deletion stays unchanged.
"""
import argparse,hashlib,os,subprocess,time
from pathlib import Path
import numpy as np
from validate_interface_closure import ROOT,DT,make,write
from validate_3d import run_node

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--before',action='store_true');a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);os.environ['OMP_NUM_THREADS']='1'
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    write(out/'job.json',dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),before=a.before,scope=__doc__))
    rec=dict(passed=False,cases=[]);start=time.perf_counter()
    try:
        for mode in (0,1):
            for existing in (False,True):
                d=out/f'm{mode}_existing{int(existing)}';ind=d/'input';ind.mkdir(parents=True)
                for axis,shape in [('u',(17,12)),('v',(16,13)),('w',(16,12,2))]:
                    for suffix,value in [('mass',1.),('mom',0.),('phase',.2)]:
                        np.save(ind/f'{axis}_{suffix}.npy',np.full(shape,value,np.float32).ravel())
                arrays=dict(positions=np.array([[8,6,.5],[12,6,.5]],np.float32),
                    velocities=np.zeros((2,3),np.float32),particle_phase=np.array([1,0],np.float32),
                    radius_in=np.array([5e-7 if existing else 0.,0.],np.float32))
                for k,v in arrays.items():np.save(ind/f'{k}.npy',v)
                stats=run_node(str(ind),str(d/'cpp'),nx=16,ny=12,nz=1,dt=DT,gravity=0,
                    phase_model=2,escape=2,esc_phi=.3,rho0_face=3.1015625,droplet_model=1,
                    droplet_transfer_model=mode,droplet_r_median=5e-7,droplet_r_sigma=0.,
                    droplet_r_min=1e-6,droplet_evap_k=0.,extra_arrays=('radius_in',),
                    extra_outputs=('radius_out','out_mass'))
                m=np.load(d/'cpp/out_mass.npy');r=np.load(d/'cpp/radius_out.npy')
                s=make(mode=mode);s.pos=arrays['positions'][:,:2].astype(float)
                s.typ=np.array([1,0]);s.spray=np.zeros(2,bool);s.radius=arrays['radius_in'].astype(float)
                s.droplet.update(r_median=float(np.float32(5e-7)),r_min=float(np.float32(1e-6)),r_sigma=0.)
                gone=s._droplet_radii(np.array([True,False]),DT)
                expected=1 if a.before or mode==0 else 2
                assert len(m)==len(r)==expected and int((~gone).sum())==expected
                if expected==2:
                    assert m[0]==0 and r[0]==np.float32(1e-6) and s.radius[0]==r[0]
                assert stats['final_residual']==0
                rec['cases'].append(dict(mode=mode,existing=existing,remaining=len(m),
                    water_preserved=bool(expected==2),radius=r.tolist()))
        rec['passed']=True
    except BaseException as e:rec['failure']=repr(e);raise
    finally:
        rec['seconds']=time.perf_counter()-start;write(out/'receipt.json',rec);print(rec,flush=True)

if __name__=='__main__':main()
