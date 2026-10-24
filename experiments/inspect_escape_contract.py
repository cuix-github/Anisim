"""Read-only figures from the classification audit; no synthetic fluid motion."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    a.out.mkdir(parents=True,exist_ok=False)
    with np.load(a.source/'step_7200/classification.npz') as q:
        pos=q['pos'];water=(q['typ']==1)&~q['spr'];old=q['current'];new=q['paper'];phi=q['phi']
        fig,axes=plt.subplots(1,2,figsize=(12,6),layout='constrained')
        axes[0].scatter(*pos[water].T,s=.5,color='#247fbe',label='Primary water')
        for mask,color,name in [(old&~new,'#f18a32','Count rule only'),(new&~old,'#2cbd74','Phase rule only')]:
            axes[0].scatter(*pos[mask].T,s=8,c=color,label=f'{name}: {mask.sum()}')
        axes[0].set(xlim=(0,64),ylim=(0,80),aspect='equal',xlabel='x (cells)',ylabel='y (cells)',title='Same cached state at 18.0025 s')
        axes[0].legend(fontsize=8)
        axes[1].hist(phi[old],bins=np.linspace(0,1,21),color='#f18a32',edgecolor='white')
        axes[1].axvline(.3,c='#2cbd74',ls='--',label='Paper liquid-candidate cutoff')
        axes[1].set(xlabel='Sampled corrected phase (0 air, 1 water)',ylabel='Water particles selected by count rule',
                    title='The count criterion is not the phase criterion')
        axes[1].legend(fontsize=8)
        fig.suptitle('Classifier audit: field sampling is averaged bilinear MAC-face interpolation',fontsize=12)
        fig.savefig(a.out/'classification_18s.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
