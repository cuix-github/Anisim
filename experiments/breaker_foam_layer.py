"""Static, render-only surface foam: build and inspect native Arnold volume points.

No solver, no ray tracing, no new physical bubbles. The two arms share positions,
radii and total integrated extinction; only its spatial distribution changes.
"""
import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import breaker_deep_water as base
from breaker_interface_mask import smoothstep

MASK=base.ROOT/'out/runs/breaker_interface_scatter_20261019'
OUT=base.ROOT/'out/runs/breaker_foam_layer_20261019'
RADIUS=.0225  # metres; 4.5 cm total thickness, not physical bubble diameter
OFFSET=.0075 # centre relative to the sampled water surface, metres
AREA=.014**2 # integration area per quadrature point, square metres
TAU=.8       # target normal optical depth at full mask before overlaps/curvature
STEP=.005    # 5 mm in metres; converted to scene units only on export
MAYA_METERS_PER_UNIT=.01  # Match the frozen Maya ASS; simulation metre is 4 scene units.


def ass_units(path):
    # Read only the header, never the potentially huge base85 geometry lines.
    with Path(path).open() as f:
        for _ in range(12):
            line=f.readline()
            if line.startswith('### meters_per_unit:'):
                return float(line.split(':',1)[1])
    raise AssertionError('Missing ASS unit metadata: '+str(path))


def noise(p,length,seed):
    rng=np.random.default_rng(seed)
    v=rng.normal(size=(12,3));v/=np.linalg.norm(v,axis=1)[:,None]
    v*=rng.uniform(.7,1.3,(12,1))*2*np.pi/length
    return np.sin(p@v.T+rng.uniform(0,2*np.pi,12)).sum(1)/np.sqrt(6)


def project(p):
    theta=np.radians(12);u=-p[:,0]*np.cos(theta)-p[:,2]*np.sin(theta)
    cu=-145.5*np.cos(theta)-70*np.sin(theta)
    return np.column_stack([(u-cu+7.5)*64,(17.71875-p[:,1])*64])


def arnold():
    dll=os.add_dll_directory('C:/Program Files/Autodesk/Arnold/Maya2027/bin')
    sys.path.insert(0,'C:/Program Files/Autodesk/Arnold/Maya2027/scripts')
    import arnold as a
    return a,dll


def export_points(path,points,density):
    a,dll=arnold();a.AiBegin();u=a.AiUniverse()
    try:
        # Coordinates already use the frozen Maya scene's convention. Omitting
        # this metadata scales them another 100 times when loaded as a stand-in.
        a.AiNodeSetFlt(a.AiUniverseGetOptions(u),'meters_per_unit',MAYA_METERS_PER_UNIT)
        # Retain the previous arm's density grouping: only unit metadata changes.
        # The prior failure did not establish a fault in per-point user data.
        group=np.floor(density/max(density.max(),1e-20)*63).astype(int)
        group[density==0]=-1
        exported=np.zeros(len(points));groups=[]
        for bucket in np.unique(group):
            if bucket<0:continue
            index=np.flatnonzero(group==bucket);value=float(density[index].mean())
            name=f'surfaceFoam_{bucket:02d}';shader=f'foamVolume_{bucket:02d}'
            node=a.AiNode(u,'points',name)
            pos=np.ascontiguousarray(points[index]*4,dtype=np.float32)
            radii=np.full(len(index),RADIUS*4,np.float32)
            for attr,array,kind in [('points',pos,a.AI_TYPE_VECTOR),('radius',radii,a.AI_TYPE_FLOAT)]:
                a.AiNodeSetArray(node,attr,a.AiArrayConvert(len(array),1,kind,array.ctypes.data_as(ctypes.c_void_p)))
            a.AiNodeSetStr(node,'mode','sphere');a.AiNodeSetFlt(node,'step_size',STEP*4)
            fog=a.AiNode(u,'standard_volume',shader)
            a.AiNodeSetStr(fog,'density_channel','');a.AiNodeSetFlt(fog,'density',value*.25)
            a.AiNodeSetFlt(fog,'scatter',.98);a.AiNodeSetRGB(fog,'scatter_color',1,1,1)
            a.AiNodeSetFlt(fog,'scatter_anisotropy',.2)
            shaders=a.AiArrayAllocate(1,1,a.AI_TYPE_NODE);assert a.AiArraySetPtr(shaders,0,fog)
            a.AiNodeSetArray(node,'shader',shaders)
            exported[index]=value;groups.append((name,shader,index,value))
        assert a.AiSceneWrite(u,str(path))
        assert ass_units(path)==MAYA_METERS_PER_UNIT
        a.AiUniverseDestroy(u);u=a.AiUniverse()
        assert a.AiSceneLoad(u,str(path))
        actual_density=np.zeros(len(points))
        for name,shader,index,value in groups:
            node=a.AiNodeLookUpByName(u,name);assert node
            for attr,expected,components in [('points',np.ascontiguousarray(points[index]*4,np.float32),3),
                                              ('radius',np.full(len(index),RADIUS*4,np.float32),1)]:
                arr=a.AiNodeGetArray(node,attr);assert a.AiArrayGetNumElements(arr)==len(index)
                ptr=a.AiArrayMap(arr)
                actual=np.ctypeslib.as_array(ctypes.cast(ptr,ctypes.POINTER(ctypes.c_float)),shape=(len(index)*components,)).copy().reshape(expected.shape)
                a.AiArrayUnmap(arr);assert np.array_equal(actual,expected),attr
            fog=a.AiNodeLookUpByName(u,shader);assert not a.AiNodeIsLinked(fog,'density')
            actual_density[index]=a.AiNodeGetFlt(fog,'density')*4
            assert abs(a.AiNodeGetFlt(node,'step_size')-STEP*4)<1e-7
        error=float(np.abs(actual_density-density).sum()/density.sum())
        budget_error=float(abs(actual_density.sum()-density.sum())/density.sum())
        assert error<.08 and budget_error<1e-6,(error,budget_error)
        return dict(status='passed',meters_per_unit=ass_units(path),points=int((density>0).sum()),logical_points=len(points),groups=len(groups),
                    density_representation='64 constant-density bins; zero-density points omitted',
                    density_quantization_relative_l1=error,budget_relative_error=budget_error,
                    integrated_extinction_m2=float(actual_density.sum()*4*np.pi/3*RADIUS**3),point_radius_m=RADIUS,step_m=STEP,
                    density_range_per_m=[float(density.min()),float(density.max())],ass_sha256=base.sha(path))
    finally:
        a.AiUniverseDestroy(u);a.AiEnd();dll.close()


def diagnostic(out,p,t,mask,surface,centers,weights,factors,face_ids,area):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    fig,axes=plt.subplots(2,2,figsize=(16,10),facecolor='#171b21')
    tm=mask[t].mean(1);pixel=project(p)
    for i,arm in enumerate(['A','B']):
        f=np.ones(len(weights)) if arm=='A' else factors
        perface=np.bincount(face_ids,weights=weights*f,minlength=len(t))/np.maximum(area,1e-20)
        for j in range(2):
            ax=axes[i,j];ax.set_facecolor('#080b0f')
            xy=p[:,[0,2]] if j else pixel
            ax.add_collection(PolyCollection(xy[t],facecolors=plt.get_cmap('inferno')(np.clip(perface/2,0,1)),edgecolors='none'))
            if j:ax.set(xlim=(p[:,0].min(),p[:,0].max()),ylim=(p[:,2].min(),p[:,2].max()),xlabel='x / m',ylabel='z / m')
            else:ax.set(xlim=(230,720),ylim=(360,140),xlabel='pixel x',ylabel='pixel y')
            ax.set_aspect('equal');ax.set_title(arm+' / '+('top support' if j else 'side support'),color='white');ax.tick_params(colors='white')
    fig.suptitle('Offline support check / triangle averages only, NOT the fine density field',color='white')
    fig.tight_layout();fig.savefig(out/'support.png',dpi=120);plt.close(fig)
    # Point-level detail to avoid hiding the sub-triangle structure in averages.
    pix=project(surface);sel=(pix[:,0]>330)&(pix[:,0]<630)&(pix[:,1]>160)&(pix[:,1]<310)
    fig,axs=plt.subplots(1,2,figsize=(14,5),facecolor='#171b21')
    for ax,f,label in zip(axs,[np.ones(len(weights)),factors],['A: uniform density budget','B: clustered density, same integrated budget']):
        ax.scatter(pix[sel,0],pix[sel,1],s=.2,c=f[sel],cmap='inferno',vmin=0,vmax=3,rasterized=True)
        ax.set(xlim=(330,630),ylim=(310,160),aspect='equal',title=label);ax.set_facecolor('black')
    fig.tight_layout();fig.savefig(out/'density_detail.png',dpi=160);plt.close(fig)


def build(out):
    assert not base.git('status','--porcelain'),'Commit instrument first'
    assert not out.exists();out.mkdir(parents=True);started=time.perf_counter()
    assert json.loads((MASK/'mask_review.json').read_text())['status']=='passed'
    z=np.load(MASK/'mask.npz');p=z['points_m'];t=z['triangles'];mask=z['mask'];normals=z['normals']
    q=p[t];area=np.linalg.norm(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]),axis=1)/2
    eligible=np.flatnonzero(mask[t].max(1)>.001)
    counts=np.maximum(1,np.ceil(area[eligible]/AREA).astype(int))
    ids=np.repeat(eligible,counts);local=np.arange(len(ids))-np.repeat(np.cumsum(counts)-counts,counts)
    # Equal-area stratification, with irrational rotation within each triangle.
    u=(local+.5)/np.repeat(counts,counts);v=np.mod(local*.61803398875+ids*.754877666,1)
    su=np.sqrt(u);bary=np.column_stack([1-su,su*(1-v),su*v])
    surface=np.einsum('ni,nij->nj',bary,q[ids])
    normal=np.einsum('ni,nij->nj',bary,normals[t[ids]]);normal/=np.linalg.norm(normal,axis=1)[:,None]
    m=np.einsum('ni,ni->n',bary,mask[t[ids]])
    # Interpolating a vertex mask can leak across a triangle straddling a crop
    # exclusion. Apply the analytic exclusion again at each actual sample.
    low=p.min(0);high=p.max(0)
    cut=np.minimum.reduce([surface[:,0]-low[0],high[0]-surface[:,0],surface[:,2]-low[2],high[2]-surface[:,2]])
    m*=smoothstep(.6,1.0,cut)
    weights=area[ids]/np.repeat(counts,counts)*m
    centers=surface+normal*OFFSET
    factors=smoothstep(-.25,.5,noise(surface,.24,774))*smoothstep(-.35,.45,noise(surface,.065,991))
    factors/=np.average(factors,weights=weights)
    volume=4*np.pi/3*RADIUS**3;dA=TAU*weights/volume;dB=dA*factors
    budget=[float(np.sum(d)*volume) for d in [dA,dB]]
    assert abs(budget[0]-budget[1])/budget[0]<1e-12
    edge=np.minimum.reduce([centers[:,0]-low[0],high[0]-centers[:,0],centers[:,2]-low[2],high[2]-centers[:,2]])-RADIUS
    active=weights>1e-8
    assert edge[active].min()>.5
    assert np.max(np.linalg.norm(centers-surface,axis=1))+RADIUS<.0300001
    assert np.isfinite(dB).all() and dB.min()>=0
    np.savez_compressed(out/'layer.npz',points_m=centers,surface_m=surface,normal=normal,weights_m2=weights,
                        density_A_per_m=dA,density_B_per_m=dB,factors=factors,source_triangles=ids)
    exports={}
    for arm,d in [('A',dA),('B',dB)]:
        folder=out/arm;folder.mkdir();exports[arm]=export_points(folder/'foam.ass',centers,d)
    diagnostic(out,p,t,mask,surface,centers,weights,factors,ids,area)
    record=dict(status='numeric_passed_visual_review_pending',git_sha=base.git('rev-parse','HEAD'),
                mask_sha256=base.sha(MASK/'mask.npz'),layer_sha256=base.sha(out/'layer.npz'),
                points=len(centers),sphere_diameter_m=2*RADIUS,center_offset_m=OFFSET,
                distance_to_source_surface_upper_bound_m=OFFSET+RADIUS,
                nominal_inward_extent_m=RADIUS-OFFSET,nominal_outward_extent_m=RADIUS+OFFSET,
                crop_edge_min_m=float(edge[active].min()),normal_optical_depth_parameter=TAU,
                integrated_extinction_m2=budget,clustering_length_m=[.24,.065],
                zero_weight_fraction_B=float((factors==0).mean()),factor_percentiles=np.percentile(factors,[0,25,50,75,95,100]).tolist(),
                exports=exports,simulation=False,physical_bubble_model=False,wall_s=time.perf_counter()-started)
    base.save(out/'layer_check.json',record);print(json.dumps(record,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--out',type=Path,default=OUT)
    build(parser.parse_args().out.resolve())
