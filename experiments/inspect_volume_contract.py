"""Read-only before/after evidence for the bounded volume candidate."""
import argparse
import json
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'gates')]
from validate_volume_contract import make


def times(row):
    dt=float(np.float32(.0025))
    return np.r_[0,(np.arange(0,row['steps'],40)+1)*dt,row['steps']*dt]


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',type=Path,required=True);a=p.parse_args()
    root=a.run.resolve();out=root/'visuals';out.mkdir(exist_ok=True)
    density=json.loads((root/'density/receipt.json').read_text())
    initial=make(.5,'density');final=make(.5,'density')
    final.pos=np.load(root/'density/bifrost/final.npy')[:,:2].astype(float)
    fig,ax=plt.subplots(2,3,figsize=(15,8),layout='constrained')
    fig.suptitle('Volume-correction candidate | actual Bifrost graph | local tests only',fontsize=16)
    for j,s,title in [(0,initial,'Before: controlled density error'),(1,final,'After: 1 simulated second')]:
        r=s.number_density()/s.density_reference()-1
        im=ax[0,j].imshow(r.T,origin='lower',extent=(0,24,0,24),vmin=-.8,vmax=.8,cmap='coolwarm')
        ax[0,j].set(title=title,xlabel='x (cells)',ylabel='y (cells)')
        fig.colorbar(im,ax=ax[0,j],label='Number density / reference - 1')
    for name,style in [('python','--'),('bifrost','-')]:
        ax[0,2].plot(times(density),[r['density_log_rms'] for r in density[name]],style,label=name)
    ax[0,2].set(title='Controlled density relaxation',xlabel='Simulated seconds',ylabel='Centered log-density RMS')
    ax[0,2].legend()
    for name,color in [('rest_off','#a65a29'),('rest_on','#147f78')]:
        row=json.loads((root/name/'receipt.json').read_text());t=times(row)
        ax[1,0].plot(t,[r['rms'] for r in row['bifrost']],label=name,color=color)
        ax[1,1].plot(t,[r['density_log_rms'] for r in row['bifrost']],label=name,color=color)
        ax[1,2].plot(t,[r['energy']/row['initial']['energy'] for r in row['bifrost']],label=name,color=color)
    for k,title,ylabel in [(0,'Still-water motion','Water RMS (cells/s)'),
        (1,'Still-water sampling density','Centered log-density RMS'),(2,'Still-water mechanical energy','Energy / initial')]:
        ax[1,k].set(title=title,xlabel='Simulated seconds',ylabel=ylabel);ax[1,k].legend()
    ax[1,0].axhline(.1,color='grey',ls=':',label='final RMS limit')
    for axy in ax.ravel():axy.grid(alpha=.2)
    fig.savefig(out/'comparison.png',dpi=135);plt.close(fig)
    print(out/'comparison.png')


if __name__=='__main__':main()
