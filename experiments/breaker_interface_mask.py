"""Offline curvature mask on the frozen 24.9 s mesh; no Maya or render launch.

Curvature is a heuristic for aerated interfaces, not a measured bubble fraction.
The sphere/plane gates check the estimator before it is used on the breaker.
"""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
MESH = ROOT/'out/runs/breaker_surface_diagnosis_20261019_v3/control_mesh.npz'
DEFAULT = ROOT/'out/runs/breaker_interface_scatter_20261019'


def smoothstep(lo, hi, x):
    u=np.clip((x-lo)/(hi-lo),0,1)
    return u*u*(3-2*u)


def curvature(p, t):
    """Signed mean curvature H=(k1+k2)/2, positive on an outward sphere."""
    q=p[t]
    cross=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])
    twice_area=np.linalg.norm(cross,axis=1)
    assert np.all(twice_area>1e-15)
    area=np.zeros(len(p)); normal=np.zeros_like(p); lap=np.zeros_like(p)
    for corner in range(3):
        np.add.at(area,t[:,corner],twice_area/6)
        np.add.at(normal,t[:,corner],cross)
        i,j,k=corner,(corner+1)%3,(corner+2)%3
        cot=np.einsum('ij,ij->i',q[:,j]-q[:,i],q[:,k]-q[:,i])/twice_area
        d=(q[:,k]-q[:,j])*cot[:,None]
        np.add.at(lap,t[:,j],d); np.add.at(lap,t[:,k],-d)
    normal/=np.maximum(np.linalg.norm(normal,axis=1)[:,None],1e-15)
    h=-np.einsum('ij,ij->i',lap,normal)/(4*area)
    edges=np.concatenate([t[:,[0,1]],t[:,[1,2]],t[:,[2,0]]])
    edges.sort(axis=1)
    edges,counts=np.unique(edges,axis=0,return_counts=True)
    return h,area,normal,edges,counts,twice_area/2


def filter_field(h,p,edges,iterations=4):
    # Four normalized graph-diffusion passes. Edge weights suppress jumps
    # longer than 0.125 m; no vertex or surface normal is modified.
    a,b=edges.T
    weights=np.exp(-np.sum((p[a]-p[b])**2,axis=1)/(2*.125**2))
    den=np.bincount(a,weights,minlength=len(p))+np.bincount(b,weights,minlength=len(p))
    for _ in range(iterations):
        avg=(np.bincount(a,weights*h[b],minlength=len(p))+
             np.bincount(b,weights*h[a],minlength=len(p)))/np.maximum(den,1e-15)
        h=.5*h+.5*avg
    return h


def analytic_gates():
    p=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)
    t=np.array([[0,2,4],[2,1,4],[1,3,4],[3,0,4],
                [2,0,5],[1,2,5],[3,1,5],[0,3,5]])
    for _ in range(4):
        pts=p.tolist(); lookup={}; faces=[]
        def middle(a,b):
            key=tuple(sorted((a,b)))
            if key not in lookup:
                v=(p[a]+p[b])*.5; v/=np.linalg.norm(v)
                lookup[key]=len(pts); pts.append(v.tolist())
            return lookup[key]
        for a,b,c in t:
            ab,bc,ca=middle(a,b),middle(b,c),middle(c,a)
            faces.extend([[a,ab,ca],[ab,b,bc],[ca,bc,c],[ab,bc,ca]])
        p=np.array(pts);t=np.array(faces)
    h,area,_,_,counts,_=curvature(p,t)
    mean=float(np.average(h,weights=area)); median=float(np.median(h))
    assert np.all(counts==2) and abs(mean-1)<.02 and abs(median-1)<.02,(mean,median)
    x,z=np.meshgrid(np.arange(8,dtype=float),np.arange(8,dtype=float),indexing='ij')
    p=np.column_stack([x.ravel(),np.zeros(x.size),z.ravel()]); tris=[]
    for i in range(7):
        for j in range(7):
            a=i*8+j; tris.extend([[a,a+1,a+8],[a+1,a+9,a+8]])
    ph,*_=curvature(p,np.array(tris))
    assert np.max(np.abs(ph))<1e-12
    return dict(unit_sphere_area_mean_H=mean,unit_sphere_median_H=median,
                plane_max_abs_H=float(np.max(np.abs(ph))),status='passed')


def diagrams(p,t,mask,out,center=(145.5,13.5,70),width_m=15):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    color=plt.get_cmap('inferno')(mask[t].mean(axis=1))
    theta=np.radians(12); u=-p[:,0]*np.cos(theta)-p[:,2]*np.sin(theta)
    depth=-p[:,0]*np.sin(theta)+p[:,2]*np.cos(theta)
    center_u=-center[0]*np.cos(theta)-center[2]*np.sin(theta)
    half_height=width_m*540/960/2
    views=[(f'Side / exact {width_m:g} m render framing',np.column_stack([u,p[:,1]]),depth,
            (center_u-width_m/2,center_u+width_m/2),(center[1]-half_height,center[1]+half_height)),
           ('Top / full crop; artificial edges must remain dark',p[:,[0,2]],-p[:,1],
            (p[:,0].min(),p[:,0].max()),(p[:,2].min(),p[:,2].max()))]
    fig,axs=plt.subplots(2,1,figsize=(12,10),facecolor='#171b21')
    for ax,(title,xy,d,xlim,ylim) in zip(axs,views):
        # Far triangles first; these are unlit data plots, not Arnold previews.
        order=np.argsort(d[t].mean(axis=1))[::-1]
        ax.add_collection(PolyCollection(xy[t[order]],facecolors=color[order],edgecolors='none'))
        ax.set(xlim=xlim,ylim=ylim,aspect='equal',title=title)
        ax.set_facecolor('black');ax.tick_params(colors='white');ax.title.set_color('white')
    fig.suptitle('Curvature mask: dark = 0, bright = 1 (offline diagnostic)',color='white')
    fig.tight_layout();fig.savefig(out/'mask_views.png',dpi=140);plt.close(fig)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,default=DEFAULT)
    args=parser.parse_args();out=args.out.resolve()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    assert not out.exists();out.mkdir(parents=True)
    start=time.perf_counter();gates=analytic_gates()
    with np.load(MESH) as z:
        p=z['points'].astype(float)*.25;t=z['triangles']
    h,area,normals,edges,counts,ta=curvature(p,t)
    assert np.all(counts==2),'Frozen mesh has open or nonmanifold edges'
    filtered=filter_field(h,p,edges)
    # Inverse radius 2..4 /m corresponds to 0.5..0.25 m radius on a sphere.
    curved=smoothstep(2,4,np.abs(filtered))
    lo=p.min(axis=0);hi=p.max(axis=0)
    cut_distance=np.minimum.reduce([p[:,0]-lo[0],hi[0]-p[:,0],p[:,2]-lo[2],hi[2]-p[:,2]])
    cut=smoothstep(.6,1.0,cut_distance)
    # This experiment targets the curl above the basin body. This spatial gate
    # is explicit, not a claim of a physically derived foam generation model.
    height=smoothstep(12.25,13.0,p[:,1])
    mask=curved*cut*height
    tm=mask[t].mean(axis=1)
    calm=(np.abs(filtered[t]).max(axis=1)<1)
    cut_tri=(cut_distance[t].min(axis=1)<.5)
    record=dict(status='numeric_checks_passed_visual_review_pending',gates=gates,
                mesh_sha256=hashlib.file_digest(MESH.open('rb'),'sha256').hexdigest(),
                git_sha=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                total_area_m2=float(ta.sum()),selected_area_m2=float(ta[tm>.05].sum()),
                strong_area_m2=float(ta[tm>.5].sum()),
                weighted_area_fraction=float(np.sum(ta*tm)/ta.sum()),
                selected_area_fraction=float(ta[tm>.05].sum()/ta.sum()),
                calm_max_weight=float(tm[calm].max()),cut_max_weight=float(tm[cut_tri].max()),
                selected_triangles=int((tm>.01).sum()),
                filter=dict(iterations=4,edge_gaussian_sigma_m=.125,blend=.5),
                curvature_ramp_per_m=[2,4],cut_exclusion_m=[.6,1.0],height_gate_m=[12.25,13.0],
                H_percentiles=np.percentile(filtered,[0,1,25,50,75,99,100]).tolist())
    assert record['calm_max_weight']==0 and record['cut_max_weight']==0,record
    assert 0<record['selected_area_fraction']<.2,record
    np.savez_compressed(out/'mask.npz',points_m=p,triangles=t,normals=normals,
                        curvature=filtered,raw_curvature=h,mask=mask,cut_distance_m=cut_distance)
    np.save(out/'support.0000.npy',(p[mask>.001]*4).astype(np.float32))
    diagrams(p,t,mask,out)
    record['wall_s']=time.perf_counter()-start
    (out/'mask_receipt.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
