"""Read-only transport ablation report; draws cached particles, no simulation."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np


def read(p):return json.loads(p.read_text())


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    root=a.root.resolve();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    result=dict(scope='Uniform 2D; escape/emission/ST/adapt off. Not full-spray or production acceptance.',runs={})
    for d in sorted(root.iterdir()):
        if not d.is_dir() or not (d/'receipt.json').exists() or not (d/'job.json').exists():continue
        job,rec=read(d/'job.json'),read(d/'receipt.json')
        if 'mode' not in job:continue
        result['runs'][d.name]=dict(job=job,receipt={k:v for k,v in rec.items() if k!='comparison'})
        if job['case']=='rest':
            initial=job.get('initial') or read(Path(job['source'])/'job.json')['initial']
            drift=abs(rec['final']['water_com'][1]-initial['water_com'][1])
            result['runs'][d.name]['still_water_gate']=dict(com_y_drift=drift,
                passed=rec['final']['water_rms']<.1 and drift<.1 and rec['max_sampled_energy_ratio']<1.005,
                limits=dict(water_rms=.1,com_y_drift=.1,energy_ratio=1.005))
    (out/'summary.json').write_text(json.dumps(result,indent=2))
    # The narrow plunge completed in all four ablation arms; display every arm.
    labels=[('control','Original'),('rk3_only','RK3 only'),('extrap_only','Extension only'),('extrap_rk3','Extension + RK3')]
    fig,axes=plt.subplots(1,2,figsize=(12,4.2),layout='constrained')
    for name,label in labels:
        d=root/name;job=read(d/'job.json');rec=read(d/'receipt.json')
        rows=read(d/'metrics.json')+[dict(time=6,**rec['final'])]
        for ax,key,den in zip(axes,('energy','water_rms'),(job['initial']['energy'],1)):
            ax.plot([r['time'] for r in rows],[r[key]/den for r in rows],label=label)
    axes[0].set_ylabel('Particle mechanical energy / initial');axes[0].axhline(1,c='grey',ls=':')
    axes[1].set_ylabel('Water RMS speed (cells/s)')
    for ax in axes:ax.set_xlabel('Simulation time (s)');ax.grid(alpha=.2)
    axes[0].legend(fontsize=9);fig.suptitle('Independent NumPy | same 64 x 80 plunge | four controlled transport arms')
    fig.savefig(out/'ablation_curves.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,4,figsize=(13.5,5.5),layout='constrained')
    for ax,(name,label,is_cpp) in zip(axes,[('control','Original / NumPy',False),('cpp_control','Original / Bifrost',True),
            ('extrap_rk3','Extension + RK3 / NumPy',False),('cpp_extrap_rk3','Extension + RK3 / Bifrost',True)]):
        d=root/name
        if is_cpp:
            pos=np.load(d/'bifrost/final.npy')[:,:2];typ=np.load(d/'bifrost/final_phase.npy')>.5
        else:
            z=np.load(d/'final.npz');pos=z['pos'];typ=z['typ']==1
        for m,c,size in [(~typ,'#cdd4dd',.15),(typ,'#174bb8',.7)]:
            ax.scatter(pos[m,0],pos[m,1],s=size,c=c,linewidths=0)
        ax.set(xlim=(0,64),ylim=(0,80),aspect='equal',title=label);ax.set_xlabel('x (cells)')
    axes[0].set_ylabel('y (cells)')
    fig.suptitle('Actual 6 s states | blue: water, grey: air | escape / spray OFF\nTransport changes are modest; independent height-tail gate still fails')
    fig.savefig(out/'before_after.png',dpi=180);plt.close(fig)
    print(json.dumps({n:dict(mode=r['job']['mode'],case=r['job']['case'],
        energy=r['receipt']['max_sampled_energy_ratio'],similarity=r['receipt'].get('similarity_pass'),
        rest=r.get('still_water_gate')) for n,r in result['runs'].items()},indent=2))


if __name__=='__main__':main()
