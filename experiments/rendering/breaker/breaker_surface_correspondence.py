"""Offline stock-Bifrost contour audit; no solver, Maya, or Arnold launch.

Control reproduces the frozen mesh before half-scale contouring is permitted.
The particle spheres, volume resolution and curvature flow stay unchanged.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT/'maya'))
import generate_particle_surface as gen
import breaker_interface_mask as interface

BIF = Path('C:/Program Files/Autodesk/Bifrost/Maya2027/3.1.0.8/bifrost')
CACHE = ROOT/'out/runs/wavesea4/cache/breaking_section'
FROZEN = ROOT/'out/runs/breaker_surface_diagnosis_20261019_v3/control_mesh.npz'
OUT = ROOT/'out/runs/breaker_surface_correspondence_20261019'
KNOBS = dict(point_size=.5, detail_size=.5, offset=0, smooth_iterations=0,
             smoothing=.1, mcf_iterations=10, mask_sigma=2, mask_theta=2,
             mask_zeta=5, mask_floor=0, mask_dilate=1)
CENTER = np.array([144.565, 14.69, 75.66])


def save(path, data):
    path.write_text(json.dumps(data, indent=2)+'\n')


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def clean_sha():
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    return subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()


def graph(scale, resample=False):
    surf = gen.surface('particle_surface')
    # Expose the exact production SDF, without changing its construction.
    final = next(c['source'] for c in surf['connections'] if c['target']=='v2m.volumes.vol')
    surf['ports'] += [gen.P('final_volume','output','Object')]
    surf['connections'] += [{'source':final,'target':'.final_volume'}]
    if resample:
        surf['compoundNodes'].append(dict(nodeName='refine',nodeType='Geometry::Volume::resample_volume'))
        surf['values'].append(dict(valueName='refine.detail_size_scale',valueType='float',value='0.5f'))
        surf['connections'].append(dict(source=final,target='refine.source_volume'))
        for connection in surf['connections']:
            if connection['target']=='v2m.volumes.vol':
                connection['source']='refine.resampled_volume'
    if scale != 1:
        surf['values'] += [dict(valueName='v2m.detail_size_scale',valueType='float',value=f'{scale}f'),
                           dict(valueName='v2m.enable_subdivision',valueType='bool',value='true')]
    ports = [gen.P(k,'input',t,d) for k,t,d in gen.KNOBS]
    ports += [gen.P('liquid_pattern','input','string',''),gen.P('radius_pattern','input','string',''),
              gen.P('frame','input','float','1f'),gen.P('path_probes','input','string','')]
    nodes = [dict(nodeName='surf',nodeType='particle_surface'),
             dict(nodeName='topo',nodeType='Geometry::Mesh::get_mesh_structure'),
             dict(nodeName='probe_type',valueType=gen.F3),
             dict(nodeName='read',nodeType='File::NumPy::read_NumPy'),
             dict(nodeName='far',valueType='float'),
             dict(nodeName='sample',nodeType='Geometry::Query::sample_volume')]
    values = [dict(valueName='probe_type.value',valueType=gen.F3,value=[gen.f3(0,0,0)]),
              dict(valueName='read.frame',valueType='long',value='0L'),
              dict(valueName='far.value',valueType='float',value='100f'),
              dict(valueName='sample.property',valueType='string',value='voxel_signed_distance'),
              dict(valueName='sample.sampler',valueType='Geometry::Query::SamplerType',value='kLinear')]
    conns = [dict(source='.'+p,target='surf.'+p) for p in [*KNOBS,'liquid_pattern','radius_pattern','frame']]
    conns += [dict(source=a,target=b) for a,b in [
        ('surf.mesh','topo.mesh'),('probe_type.output','read.type'),('.path_probes','read.file_path'),
        ('read.data','sample.positions'),('surf.final_volume','sample.volume'),('far.output','sample.default')]]
    for name,source in [('positions','topo.point_position'),('indices','topo.face_vertex'),
                        ('offsets','topo.face_offset'),('sdf','sample.sampled_data')]:
        if name in ('indices','offsets'):
            nodes.append(dict(nodeName='cast_'+name,nodeType='Core::Type_Conversion::to_long'))
            conns.append(dict(source=source,target='cast_'+name+'.from'))
            source='cast_'+name+'.long'
        gen.writer(nodes,values,conns,ports,name,source,None)
    return dict(gen.HEADER,compounds=[dict(name='User::PFFlip::surface_correspondence',ports=ports,
                compounds=[surf],compoundNodes=nodes,values=values,connections=conns)])


def run(arm):
    commit = clean_sha()
    if arm=='fine':
        assert json.loads((OUT/'control/receipt.json').read_text())['frozen_control_match']
    if arm=='local_refine':
        assert json.loads((OUT/'local_control/receipt.json').read_text())['local_control_match']
    out = OUT/arm
    if arm=='smoke':
        attempt=1
        while out.exists():
            attempt+=1;out=OUT/f'smoke_{attempt}'
    assert not out.exists(), out
    out.mkdir(parents=True)
    scale = .5 if arm=='fine' else 1
    cache = CACHE
    knobs = KNOBS.copy()
    source_hashes={p.name:sha(p) for p in CACHE.glob('*.npy')}
    if arm.startswith('local_'):
        # Diagnostic crop only: 4 m padding around a 20 cm fragment. Require
        # exact local mesh correspondence before interpreting its refined arm.
        cache=out/'cache';cache.mkdir()
        p=np.load(CACHE/'liquid.0001.npy');r=np.load(CACHE/'radius.0001.npy')
        keep=(np.abs(p*.25-CENTER)<4).all(1)
        np.save(cache/'liquid.0001.npy',p[keep]);np.save(cache/'radius.0001.npy',r[keep])
    if arm=='smoke':
        cache = out/'cache';cache.mkdir()
        x,y,z = np.meshgrid(np.arange(5)*.5,np.arange(5)*.5,np.arange(5)*.5,indexing='ij')
        pts = np.column_stack([x.ravel(),y.ravel(),z.ravel()]).astype(np.float32)
        np.save(cache/'liquid.0001.npy',pts);np.save(cache/'radius.0001.npy',np.zeros(len(pts),np.float32))
        knobs['mcf_iterations']=0
        probes = pts[:5]
    else:
        # Fixed 2 mm slice through the known floating fragment, in simulation metres.
        ax = np.linspace(-.2,.2,201)
        x,y = np.meshgrid(ax+ CENTER[0],ax+CENTER[1])
        probes = np.column_stack([x.ravel(),y.ravel(),np.full(x.size,CENTER[2])])*4
    np.save(out/'probes.npy',probes.astype(np.float32))
    save(out/'graph.json',graph(scale,resample=arm=='local_refine'))
    hashes = {p.name:sha(p) for p in cache.glob('*.npy')}
    args = [str(BIF/'bin/bifcmd.exe'),str(out/'graph.json')]
    inputs = dict(knobs, liquid_pattern=(cache/'liquid.####').as_posix(),
                  radius_pattern=(cache/'radius.####').as_posix(),frame=1,
                  path_probes=(out/'probes.npy').as_posix())
    for name in ['positions','indices','offsets','sdf']:
        inputs['path_'+name]=(out/(name+'.npy')).as_posix()
    for k,v in inputs.items():args += ['--set-port',k,str(v)]
    env = os.environ.copy();env['PATH']=str(BIF/'bin')+';'+str(BIF/'thirdparty/bin')+';'+env['PATH']
    record = dict(git_sha=commit,arm=arm,inputs=inputs,mesh_detail_scale=scale,
                  graph_sha256=sha(out/'graph.json'),cache_hashes_before=hashes,
                  source_cache_hashes_before=source_hashes,status='running')
    save(out/'receipt.json',record)
    start=time.perf_counter()
    with (out/'bifcmd.log').open('w') as log:
        result=subprocess.run(args,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=600)
    record.update(wall_s=time.perf_counter()-start,exit_code=result.returncode)
    save(out/'receipt.json',record)
    log=(out/'bifcmd.log').read_text()
    assert result.returncode==0,log[-4000:]
    for name in ['positions','indices','offsets','sdf']:
        assert f"'ok_{name}': (bool) true" in log,log[-4000:]
    p=np.load(out/'positions.npy');idx=np.load(out/'indices.npy');off=np.load(out/'offsets.npy')
    assert len(idx)%3==0 and np.all(np.diff(off)==3),(idx.shape,off)
    t=idx.reshape(-1,3).astype(np.int64)
    assert t.max()<len(p) and np.isfinite(p).all()
    np.savez_compressed(out/'mesh.npz',points=p,triangles=t)
    record.update(vertices=len(p),triangles=len(t),bounds_cells=[p.min(0).tolist(),p.max(0).tolist()])
    if arm=='control':
        z=np.load(FROZEN)
        record['frozen_position_max_error_cells']=float(np.max(np.abs(p-z['points']))) if p.shape==z['points'].shape else None
        # Face order/winding may differ across the Bifrost-to-Maya conversion.
        def faces(a):
            a=np.sort(a,axis=1);return a[np.lexsort(a.T[::-1])]
        record['frozen_control_match']=(record['frozen_position_max_error_cells'] is not None and
            record['frozen_position_max_error_cells']<1e-5 and np.array_equal(faces(t),faces(z['triangles'])))
        save(out/'receipt.json',record)
        assert record['frozen_control_match'],record
    if arm=='local_control':
        z=np.load(FROZEN)
        def fragment(points, triangles):
            keep=(np.abs(points*.25-CENTER)<.25).all(1)
            used=np.unique(triangles[keep[triangles].all(1)])
            points=points[used];return points[np.lexsort(points.T[::-1])]
        a=fragment(p,t);b=fragment(z['points'],z['triangles'])
        record['local_control_match']=np.array_equal(a,b)
        record['local_fragment_vertices']=len(a)
        save(out/'receipt.json',record)
        assert record['local_control_match'],record
    record['cache_hashes_after']={p.name:sha(p) for p in cache.glob('*.npy')}
    assert record['cache_hashes_after']==hashes
    record['source_cache_hashes_after']={p.name:sha(p) for p in CACHE.glob('*.npy')}
    assert record['source_cache_hashes_after']==source_hashes
    record['status']='passed';save(out/'receipt.json',record)
    print(json.dumps(record,indent=2))


def analyze():
    clean_sha()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    report=dict(gates=interface.analytic_gates(),arms={})
    fig,axs=plt.subplots(2,2,figsize=(12,10))
    for column,arm in enumerate(['control','fine']):
        out=OUT/arm;z=np.load(out/'mesh.npz');p=z['points'].astype(float)*.25;t=z['triangles']
        h,area,n,edges,counts,ta=interface.curvature(p,t)
        assert np.all(counts==2), 'Open/nonmanifold mesh'
        filtered=interface.filter_field(h,p,edges)
        cutdist=np.minimum.reduce([p[:,0]-p[:,0].min(),p[:,0].max()-p[:,0],p[:,2]-p[:,2].min(),p[:,2].max()-p[:,2]])
        exclusion=interface.smoothstep(.6,1,cutdist)*interface.smoothstep(12.25,13,p[:,1])
        mask=interface.smoothstep(2,4,np.abs(filtered))*exclusion
        np.savez_compressed(out/'mask.npz',points_m=p,triangles=t,normals=n,raw_curvature=h,
                            curvature=filtered,mask=mask,cut_distance_m=cutdist)
        tm=mask[t].mean(1)
        near=(np.abs(p-CENTER)<[.25,.25,.25]).all(1);faces=t[near[t].all(1)]
        xy=(p-CENTER)[:,:2]*100
        ax=axs[0,column]
        ax.add_collection(PolyCollection(xy[faces],facecolors='#b8dce9',edgecolors='#35566b',linewidths=.7,alpha=.65))
        sdf=np.load(out/'sdf.npy').reshape(201,201)*.25
        axis=np.linspace(-20,20,201)
        ax.contour(axis,axis,sdf,levels=[0],colors=['#cf3a24'],linewidths=2)
        ax.set(xlim=(-20,20),ylim=(-20,20),aspect='equal',title=f'{arm}: local triangles {len(faces)}',xlabel='x / cm',ylabel='y / cm')
        # Same physical window in both arms. Color is raw curvature, avoiding a
        # claim that four graph-filter passes have equal physical width after refinement.
        crest=(p[t].mean(1)[:,0]>140)&(p[t].mean(1)[:,0]<149)&(p[t].mean(1)[:,1]>12.5)&(p[t].mean(1)[:,2]>68)&(p[t].mean(1)[:,2]<77)
        chosen=t[crest];order=np.argsort(p[chosen,2].mean(1))[::-1]
        color=plt.get_cmap('inferno')(interface.smoothstep(2,4,np.abs(h))[chosen].mean(1))
        ax=axs[1,column];ax.add_collection(PolyCollection(p[chosen[order]][:,:,[0,1]],facecolors=color[order],edgecolors='none'))
        ax.set(xlim=(140,149),ylim=(12.5,16.5),aspect='equal',title='Raw curvature response: 2–4 /m',xlabel='x / m',ylabel='y / m')
        report['arms'][arm]=dict(vertices=len(p),triangles=len(t),local_triangles=len(faces),
            surface_area_m2=float(ta.sum()),selected_area_m2=float(ta[tm>.05].sum()),
            strong_area_m2=float(ta[tm>.5].sum()),mask_area_weighted=float((tm*ta).sum()/ta.sum()),
            cut_max=float(tm[(cutdist[t].min(1)<.5)].max()),
            curvature_filter_note='Legacy 4 graph passes; physical smoothing width differs with mesh resolution; not a matched optical A/B.')
        interface.diagrams(p,t,mask,out)
    a=np.load(OUT/'control/sdf.npy');b=np.load(OUT/'fine/sdf.npy')
    report['final_field_sample_max_difference_cells']=float(np.max(np.abs(a-b)))
    assert np.array_equal(a,b), 'Contouring experiment changed its input field'
    fig.suptitle('Offline geometry / curvature audit — no Arnold render\nRed: same sampled SDF zero contour in a fixed z slice; mesh is projected, not a slice')
    fig.tight_layout();fig.savefig(OUT/'comparison.png',dpi=140);plt.close(fig)
    save(OUT/'analysis.json',report);print(json.dumps(report,indent=2))


def analyze_local():
    commit=clean_sha()
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.collections import PolyCollection
    from breaker_corner_audit import components
    report=dict(git_sha=commit,arms={})
    fig,axs=plt.subplots(2,2,figsize=(12,10))
    frozen=np.load(OUT/'control/sdf.npy')
    for col,arm in enumerate(['local_control','local_refine']):
        out=OUT/arm;z=np.load(out/'mesh.npz');p=z['points'].astype(float)*.25;t=z['triangles']
        ids,face_ids=components(p,t)
        nearest=int(np.argmin(np.linalg.norm(p-CENTER,axis=1)))
        connected=ids==ids[nearest];faces=t[face_ids==ids[nearest]]
        h,area,n,edges,counts,ta=interface.curvature(p,t)
        filtered=interface.filter_field(h,p,edges)
        sdf=np.load(out/'sdf.npy')
        assert np.max(np.abs(sdf-frozen))<1e-5
        xy=(p-CENTER)[:,:2]*100
        local=(np.abs(p[t].mean(1)-CENTER)<[.5,.5,.5]).all(1)
        triangles=t[local];order=np.argsort(p[triangles,2].mean(1))[::-1]
        ax=axs[0,col]
        ax.add_collection(PolyCollection(xy[triangles[order]],facecolors='#b8dce9',edgecolors='#35566b',linewidths=.7,alpha=.65))
        ax.set(xlim=(-35,35),ylim=(-25,25),aspect='equal',title=arm+' / projected local mesh',xlabel='x / cm',ylabel='y / cm')
        color=plt.get_cmap('inferno')(interface.smoothstep(2,4,np.abs(filtered))[triangles].mean(1))
        ax=axs[1,col];ax.add_collection(PolyCollection(xy[triangles[order]],facecolors=color[order],edgecolors='#333',linewidths=.2))
        ax.set(xlim=(-35,35),ylim=(-25,25),aspect='equal',title='Legacy curvature mask: dark 0, bright 1',xlabel='x / cm',ylabel='y / cm')
        report['arms'][arm]=dict(nearest_vertex_distance_m=float(np.linalg.norm(p[nearest]-CENTER)),
            component_vertices=int(connected.sum()),component_triangles=len(faces),
            component_bounds_m=[p[connected].min(0).tolist(),p[connected].max(0).tolist()],
            component_filtered_H_range_per_m=[float(filtered[connected].min()),float(filtered[connected].max())],
            component_mean_mask=float(np.mean(interface.smoothstep(2,4,np.abs(filtered[connected])))),
            source_sdf_max_difference_cells=float(np.max(np.abs(sdf-frozen))),
            manifold=bool(np.all(counts==2)))
    report['interpretation']='The nearest connected component is a spatial correspondence check, not exact particle provenance. Resampling changes the discrete volume as well as the contour: not a paper-faithful replacement until field, topology and smoothing correspondence are established.'
    fig.suptitle('Local diagnostic only — unchanged cached particles / no Arnold render')
    fig.tight_layout();fig.savefig(OUT/'local_comparison.png',dpi=140);plt.close(fig)
    save(OUT/'local_analysis.json',report);print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['smoke','control','fine','local_control','local_refine','analyze','analyze_local'])
    args=parser.parse_args()
    if args.action=='analyze':analyze()
    elif args.action=='analyze_local':analyze_local()
    else:run(args.action)
