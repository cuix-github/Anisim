"""Focused unlit geometry plot for the user-circled floating diamond."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from breaker_corner_audit import SOURCE,OUT
from breaker_foam_layer import project
import breaker_deep_water as base


def main():
    r=json.loads((OUT/'analysis.json').read_text())
    particles=np.load(base.CACHE/'liquid.0001.npy',mmap_mode='r')*.25
    fig,axs=plt.subplots(1,2,figsize=(11,5));findings={}
    for ax,arm in zip(axs,['raw','control']):
        f=min(r[arm]['fragments'],key=lambda f:np.linalg.norm(np.array(f['pixel_center'])-[463,194]))
        lo,hi=np.array(f['bounds_m']);z=np.load(SOURCE/(arm+'_mesh.npz'))
        p=z['points']*.25;t=z['triangles'];inside=((p>=lo-1e-6)&(p<=hi+1e-6)).all(1)
        faces=t[inside[t].all(1)];assert len(faces)==8
        pix=project(p);center=np.array([463,194]);xy=(pix-center)*15/960*100
        ax.add_collection(PolyCollection(xy[faces],facecolors='#c3dbe3',edgecolors='#344957',linewidths=1.5,alpha=.7))
        ax.scatter(xy[inside,0],xy[inside,1],s=30,c='#29404a',label='Mesh vertices')
        nearby=((particles>lo-.05)&(particles<hi+.05)).all(1)
        pp=(project(particles[nearby])-center)*15/960*100
        ax.scatter(pp[:,0],pp[:,1],s=50,marker='+',c='#d42c24',label='Nearby cached particle centres')
        ax.set(xlim=(-15,15),ylim=(12,-12),aspect='equal',xlabel='horizontal / cm',ylabel='vertical / cm',title=arm+' / 6 vertices, 8 triangles')
        ax.legend(fontsize=8,loc='upper right');findings[arm]=f
    fig.suptitle('Floating diamond near original image pixel (463, 194) / geometry only')
    fig.tight_layout();fig.savefig(OUT/'isolated_fragment.png',dpi=150)
    findings['interpretation']='Angular silhouette exists in raw and smoothed meshes. Nearby particles are location evidence, not proof of their exact meshing contribution. This does not diagnose all attached crest corners.'
    base.save(OUT/'isolated_fragment.json',findings)
    print(json.dumps(findings,indent=2))


if __name__=='__main__':main()
