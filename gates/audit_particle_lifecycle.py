"""Prescribed-field lifecycle witness against PF-FLIP section 3.4.

The retained-bubble implementation failed this gate; that before result is
preserved. Its negative result was not a pressure failure. The expectation of
deleting escaped air particles is unchanged after the implementation repair.
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
sys.path[:0]=[str(ROOT/'gates')]
from validate_3d import run_node


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);os.environ['OMP_NUM_THREADS']='1'
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    r=dict(passed=False,cases={})
    (out/'job.json').write_text(json.dumps(dict(scope=__doc__,
        commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest()),indent=2))
    try:
        for name,phase,typ in [('water_control',1.,1),('escaped_drop',0.,1),('escaped_bubble',1.,0)]:
            d=out/name;ind=d/'input';ind.mkdir(parents=True)
            for axis,shape in [('u',(17,12)),('v',(16,13)),('w',(16,12,2))]:
                for suffix,value in [('mass',1.),('mom',0.),('phase',phase)]:
                    np.save(ind/f'{axis}_{suffix}.npy',np.full(shape,value,np.float32).ravel())
            # Keep one water sentinel: shipped write_Numpy rejects zero-length
            # output arrays. Missing files are not evidence of successful deletion.
            np.save(ind/'positions.npy',np.array([[8.,6.,.5],[12.,6.,.5]],np.float32))
            np.save(ind/'velocities.npy',np.zeros((2,3),np.float32))
            np.save(ind/'particle_phase.npy',np.array([typ,1],np.float32))
            stats=run_node(str(ind),str(d/'cpp'),nx=16,ny=12,nz=1,dt=.0025,gravity=0,
                phase_model=2,escape=2,esc_phi=.3,rho0_face=3.1015625,buoyancy=0,
                droplet_model=1,droplet_r_median=.0001,droplet_r_sigma=0.,droplet_evap_k=0.,
                extra_outputs=('radius_out','out_phase_state'))
            pos=np.load(d/'cpp/out_positions.npy');esc=np.load(d/'cpp/out_escaped.npy')
            rad=np.load(d/'cpp/radius_out.npy')
            labels=np.load(d/'cpp/out_phase_state.npy')
            r['cases'][name]=dict(input_phase=phase,particle_is_water=bool(typ),
                remaining=len(pos),escaped=esc.tolist(),radii=rad.tolist(),
                positions=pos.tolist(),labels=labels.tolist(),pressure_residual=stats['final_residual'])
        assert r['cases']['water_control']['remaining']==2
        assert r['cases']['escaped_drop']['remaining']==2 and all(r>0 for r in r['cases']['escaped_drop']['radii'])
        b=r['cases']['escaped_bubble']
        assert b['remaining']==1 and b['labels']==[1.], 'PF section 3.4 requires deleting the escaped bubble and preserving the water sentinel'
        assert b['positions']==[[12.,6.,.5]] and b['escaped']==[0.] and b['radii']==[0.]
        r['passed']=True
    except BaseException as e:
        r['failure']=repr(e);raise
    finally:
        (out/'receipt.json').write_text(json.dumps(r,indent=2))
        print(json.dumps(r,indent=2),flush=True)


if __name__=='__main__':main()
