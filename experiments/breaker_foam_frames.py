"""Prepare later-frame counterparts of the validated foam looks, without simulation.

Stock reconstruction settings, normal export, mask, and A/B foam algorithms are
kept. Source overview caches already contain the established droplet labels.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import breaker_surface_correspondence as surface
import breaker_interface_mask as masklib
import breaker_foam_layer as foam
import breaker_deep_water as base

OUT=base.ROOT/'out/runs/breaker_foam_views_20261019/frames'
SOURCE=base.ROOT/'out/runs/wavesea4/cache/overview'
FRAMES={'control':None,'later26':4,'later28':6}


def mesh(name):
    commit=surface.clean_sha();out=OUT/name;assert not out.exists();out.mkdir(parents=True)
    frame=FRAMES[name];cache=out/'cache';cache.mkdir()
    origin=base.CACHE if frame is None else SOURCE
    number=1 if frame is None else frame
    hashes={kind:base.sha(origin/f'{kind}.{number:04}.npy') for kind in ['liquid','radius','speed']}
    pos=np.load(origin/f'liquid.{number:04}.npy')
    keep=np.ones(len(pos),bool) if frame is None else ((pos[:,0]>=480)&(pos[:,2]>=240)&(pos[:,2]<312))
    for kind in hashes:
        np.save(cache/f'{kind}.0001.npy',np.load(origin/f'{kind}.{number:04}.npy')[keep])
    cache_hashes={kind:base.sha(cache/f'{kind}.0001.npy') for kind in hashes}
    g=surface.graph(1);top=g['compounds'][0]
    top['compoundNodes'] += [dict(nodeName='read_normal',nodeType='Geometry::Properties::get_geo_property')]
    top['values'] += [dict(valueName='read_normal.property',valueType='string',value='point_normal')]
    top['connections'] += [dict(source='surf.mesh',target='read_normal.geometry'),
                           dict(source='probe_type.output',target='read_normal.type')]
    surface.gen.writer(top['compoundNodes'],top['values'],top['connections'],top['ports'],
                       'normals','read_normal.data',surface.gen.F3)
    base.save(out/'graph.json',g);np.save(out/'probes.npy',np.array([[582,54,280]],np.float32))
    args=[str(surface.BIF/'bin/bifcmd.exe'),str(out/'graph.json')]
    inputs=dict(surface.KNOBS,liquid_pattern=(cache/'liquid.####').as_posix(),
                radius_pattern=(cache/'radius.####').as_posix(),frame=1,path_probes=(out/'probes.npy').as_posix())
    for key in ['positions','indices','offsets','sdf','normals']:inputs['path_'+key]=(out/(key+'.npy')).as_posix()
    for key,value in inputs.items():args += ['--set-port',key,str(value)]
    env=os.environ.copy();env['PATH']=str(surface.BIF/'bin')+';'+str(surface.BIF/'thirdparty/bin')+';'+env['PATH']
    receipt=dict(status='running',git_sha=commit,name=name,source=str(origin),source_frame=number,
                 source_hashes=hashes,cache_hashes=cache_hashes,inputs=inputs,particles=int(keep.sum()),
                 time_s=24.9 if frame is None else json.loads((SOURCE/'cache.json').read_text())['times_s'][frame-1])
    base.save(out/'mesh_receipt.json',receipt);started=time.perf_counter()
    with (out/'bifcmd.log').open('w') as log:
        p=subprocess.run(args,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=900)
    receipt.update(exit_code=p.returncode,process_s=time.perf_counter()-started)
    base.save(out/'mesh_receipt.json',receipt)
    assert p.returncode==0,(out/'bifcmd.log').read_text()[-3000:]
    log=(out/'bifcmd.log').read_text()
    for key in ['positions','indices','offsets','sdf','normals']:assert f"'ok_{key}': (bool) true" in log
    points=np.load(out/'positions.npy');t=np.load(out/'indices.npy').reshape(-1,3).astype(np.int64)
    normals=np.load(out/'normals.npy')
    assert points.shape==normals.shape and np.isfinite(normals).all()
    assert np.max(np.abs(np.linalg.norm(normals,axis=1)-1))<1e-4
    np.savez_compressed(out/'mesh.npz',points=points,triangles=t,normals=normals)
    if name=='control':
        old=np.load(surface.FROZEN);assert np.array_equal(points,old['points'])
        old_normals=np.empty_like(points)
        old_normals[old['face_vertices']]=old['normals'][old['normal_ids']]
        error=float(np.max(np.abs(normals-old_normals)))
        receipt['normal_transfer_max_error']=error
        base.save(out/'mesh_receipt.json',receipt);assert error<1e-5,error
    else:
        assert json.loads((OUT/'control/mesh_receipt.json').read_text())['status']=='passed'
    assert hashes=={kind:base.sha(origin/f'{kind}.{number:04}.npy') for kind in hashes}
    receipt.update(status='passed',vertices=len(points),triangles=len(t),mesh_sha256=base.sha(out/'mesh.npz'))
    base.save(out/'mesh_receipt.json',receipt);print(json.dumps(receipt,indent=2))


def mask(name):
    surface.clean_sha();out=OUT/name;folder=out/'mask';assert not folder.exists();folder.mkdir()
    assert json.loads((out/'mesh_receipt.json').read_text())['status']=='passed'
    z=np.load(out/'mesh.npz');p=z['points'].astype(float)*.25;t=z['triangles']
    h,area,n,edges,counts,ta=masklib.curvature(p,t);assert np.all(counts==2)
    filtered=masklib.filter_field(h,p,edges)
    lo=p.min(0);hi=p.max(0)
    cut=np.minimum.reduce([p[:,0]-lo[0],hi[0]-p[:,0],p[:,2]-lo[2],hi[2]-p[:,2]])
    m=masklib.smoothstep(2,4,np.abs(filtered))*masklib.smoothstep(.6,1,cut)*masklib.smoothstep(12.25,13,p[:,1])
    tm=m[t].mean(1);assert tm[(cut[t].min(1)<.5)].max()==0
    np.savez_compressed(folder/'mask.npz',points_m=p,triangles=t,normals=n,curvature=filtered,
                        raw_curvature=h,mask=m,cut_distance_m=cut)
    masklib.diagrams(p,t,m,folder)
    receipt=dict(status='numeric_passed_visual_review_pending',time_s=json.loads((out/'mesh_receipt.json').read_text())['time_s'],
                 selected_area_m2=float(ta[tm>.05].sum()),strong_area_m2=float(ta[tm>.5].sum()),
                 cut_max=0,height_gate_m=[12.25,13],mesh_sha256=base.sha(out/'mesh.npz'))
    base.save(folder/'mask_receipt.json',receipt);print(json.dumps(receipt,indent=2))


def layer(name):
    folder=OUT/name;foam.MASK=folder/'mask';foam.build(folder/'foam')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['mesh','mask','layer']);parser.add_argument('name',choices=list(FRAMES))
    args=parser.parse_args();globals()[args.action](args.name)
