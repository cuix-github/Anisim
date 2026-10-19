"""Inspect disconnected crest fragments in existing meshes, without rebuilding."""
import json
from pathlib import Path
import time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import breaker_deep_water as base
from breaker_foam_layer import project

SOURCE=base.ROOT/'out/runs/breaker_surface_diagnosis_20261019_v3'
OUT=base.ROOT/'out/runs/breaker_corner_audit_20261019'


def components(p,t):
    parent=np.arange(len(p))
    def root(x):
        while parent[x]!=x:
            parent[x]=parent[parent[x]];x=parent[x]
        return x
    for tri in t:
        a=root(tri[0])
        for b in tri[1:]:
            b=root(b)
            if a!=b:parent[b]=a
    ids=np.array([root(i) for i in range(len(p))]);return ids,ids[t[:,0]]


def main():
    assert not base.git('status','--porcelain');assert not OUT.exists();OUT.mkdir()
    start=time.perf_counter();report={};fig,axs=plt.subplots(3,1,figsize=(15,12))
    particles=np.load(base.CACHE/'liquid.0001.npy',mmap_mode='r')*.25
    for ax,name in zip(axs,['raw','no_mcf','control']):
        z=np.load(SOURCE/(name+'_mesh.npz'));p=z['points']*.25;t=z['triangles'];pix=project(p)
        ids,faceids=components(p,t);unique,counts=np.unique(faceids,return_counts=True)
        chosen=[]
        for comp,n in zip(unique,counts):
            if n>200:continue
            verts=np.flatnonzero(ids==comp);pc=p[verts];px=project(pc);center=px.mean(0)
            if not (280<center[0]<720 and 100<center[1]<310):continue
            lo=pc.min(0);hi=pc.max(0);near=((particles>lo-.05)&(particles<hi+.05)).all(1)
            chosen.append(dict(component=int(comp),vertices=len(verts),triangles=int(n),pixel_center=center.tolist(),
                               bounds_m=[lo.tolist(),hi.tolist()],extent_m=(hi-lo).tolist(),
                               nearby_particles_in_padded_AABB=int(near.sum())))
        chosen.sort(key=lambda r:r['pixel_center'][1])
        visible=((pix[t].mean(1)[:,0]>280)&(pix[t].mean(1)[:,0]<720)&(pix[t].mean(1)[:,1]<330))
        theta=np.radians(12);depth=-p[:,0]*np.sin(theta)+p[:,2]*np.cos(theta)
        triangles=np.flatnonzero(visible);triangles=triangles[np.argsort(depth[t[triangles]].mean(1))[::-1]]
        ax.add_collection(PolyCollection(pix[t[triangles]],facecolors='#789fa8',edgecolors='#24333b',linewidths=.25))
        for i,r in enumerate(chosen):
            x,y=r['pixel_center'];ax.annotate(f"{i}: {r['triangles']} tris",(x,y),(x+8,y-8),fontsize=8,color='#b31c24')
        ax.set(xlim=(280,720),ylim=(330,130),aspect='equal',title=name+' / actual triangular silhouette, no lighting')
        report[name]=dict(mesh_sha256=base.sha(SOURCE/(name+'_mesh.npz')),fragments=chosen,
                         normals=z['normals'].shape[0],vertices=len(p),triangles=len(t))
    fig.tight_layout();fig.savefig(OUT/'wireframes.png',dpi=150);plt.close(fig)
    report.update(git_sha=base.git('rev-parse','HEAD'),wall_s=time.perf_counter()-start,
                  note='AABB particle counts are proximity evidence, not exact particle-to-fragment provenance.')
    base.save(OUT/'analysis.json',report)
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
