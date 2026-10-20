"""Bound only the density-mask field, with a paired geometry preservation gate.

The raw level set, particle selection, MCF and mesher remain stock and unchanged.
Eight metres of empty padding exceeds the three-metre Gaussian half-width.
"""
import argparse
import os
from pathlib import Path
import subprocess
import time
import numpy as np
import breaker_six_views as study
import breaker_surface_correspondence as surface
import breaker_deep_water as base


def graph(cache):
    p=np.load(cache/'liquid.0001.npy',mmap_mode='r')
    r=np.load(cache/'radius.0001.npy',mmap_mode='r')
    p=p[r<1e-9]
    step=surface.KNOBS['detail_size']
    lo=np.floor(p.min(0)/step)*step-32
    hi=np.ceil(p.max(0)/step)*step+32
    g=surface.graph(1)
    surf=g['compounds'][0]['compounds'][0]
    targets=['k_vol.min_bounds','k_vol.max_bounds']
    surf['connections']=[c for c in surf['connections'] if c['target'] not in targets]
    for key,v in zip(targets,[lo,hi]):
        surf['values'].append(dict(valueName=key,valueType='Math::float3',value=surface.gen.f3(*map(float,v))))
    return g,dict(min_cells=lo.tolist(),max_cells=hi.tolist(),padding_m=8,detail_cells=step)


def reconstruct(cache,out):
    commit=surface.clean_sha()
    assert not out.exists()
    out.mkdir(parents=True)
    g,bounds=graph(cache)
    base.save(out/'graph.json',g)
    np.save(out/'probes.npy',np.array([[582,54,280]],np.float32))
    inputs=dict(surface.KNOBS,liquid_pattern=(cache/'liquid.####').as_posix(),
                radius_pattern=(cache/'radius.####').as_posix(),frame=1,path_probes=(out/'probes.npy').as_posix())
    for k in ['positions','indices','offsets','sdf']: inputs['path_'+k]=(out/(k+'.npy')).as_posix()
    cmd=[str(surface.BIF/'bin/bifcmd.exe'),str(out/'graph.json')]
    for k,v in inputs.items():cmd+=['--set-port',k,str(v)]
    env=os.environ.copy();env['PATH']=str(surface.BIF/'bin')+';'+str(surface.BIF/'thirdparty/bin')+';'+env['PATH']
    receipt=dict(status='running',git_sha=commit,mask_bounds=bounds,inputs=inputs,
                 cache_hashes={k:base.sha(cache/(k+'.0001.npy')) for k in ['liquid','radius']})
    base.save(out/'receipt.json',receipt)
    start=time.perf_counter();minimum=study.free_gb()
    assert minimum>30
    with (out/'bifcmd.log').open('w') as log:
        process=subprocess.Popen(cmd,env=env,stdout=log,stderr=subprocess.STDOUT)
        while process.poll() is None:
            free=study.free_gb();minimum=min(minimum,free);elapsed=time.perf_counter()-start
            base.save(out/'progress.json',dict(pid=process.pid,elapsed_s=elapsed,free_gib=free,minimum_free_gib=minimum))
            if free<5 or elapsed>900:
                process.kill();process.wait()
                receipt.update(status='failed',reason='RAM/time guard',elapsed_s=elapsed,free_gib=free)
                base.save(out/'receipt.json',receipt);raise RuntimeError(receipt)
            time.sleep(2)
    receipt.update(exit_code=process.returncode,process_s=time.perf_counter()-start,minimum_free_gib=minimum)
    base.save(out/'receipt.json',receipt)
    assert process.returncode==0
    log=(out/'bifcmd.log').read_text()
    for k in ['positions','indices','offsets','sdf']:assert f"'ok_{k}': (bool) true" in log
    p=np.load(out/'positions.npy');t=np.load(out/'indices.npy').reshape(-1,3)
    assert np.isfinite(p).all() and np.all(np.diff(np.load(out/'offsets.npy'))==3)
    receipt.update(status='passed',vertices=len(p),triangles=len(t))
    base.save(out/'receipt.json',receipt)
    return p,t,receipt


def control():
    old=study.OUT/'frames/pre/tiles/1'
    out=study.OUT/'bounded_control'
    p,t,receipt=reconstruct(old/'cache',out)
    original=np.load(old/'positions.npy');faces=np.load(old/'indices.npy').reshape(-1,3)
    assert p.shape==original.shape and t.shape==faces.shape
    error=float(np.max(np.abs(p-original)))
    exact_faces=np.array_equal(t,faces)
    receipt.update(position_max_error_cells=error,exact_faces=exact_faces,
                   geometry_preserved=error<1e-5 and exact_faces)
    base.save(out/'receipt.json',receipt)
    assert receipt['geometry_preserved'],receipt
    print(receipt)


def frame(name):
    assert study.read(study.OUT/'bounded_control/receipt.json')['geometry_preserved']
    spec=study.FRAMES[name];root=study.OUT/'frames'/name
    cache=root/'bounded_cache';assert not cache.exists();cache.mkdir(parents=True)
    p=np.load(spec['source']/f"liquid.{spec['frame']:04}.npy",mmap_mode='r')
    keep=(p[:,0]>=spec['xmin']*4)&(p[:,2]>=40)&(p[:,2]<312)
    hashes={}
    for k in ['liquid','radius']:
        src=spec['source']/f"{k}.{spec['frame']:04}.npy";hashes[k]=base.sha(src)
        np.save(cache/(k+'.0001.npy'),np.load(src,mmap_mode='r')[keep])
    count=int(keep.sum());del p,keep
    p,t,receipt=reconstruct(cache,root/'bounded_reconstruction')
    np.save(root/'positions.npy',p);np.save(root/'indices.npy',t.ravel())
    env=os.environ.copy();env.update(SIX_WORKER='normals',SIX_FRAME=name)
    with (root/'normal.log').open('w') as log:
        proc=subprocess.run([study.render.MAYAPY,str(Path(study.__file__).resolve())],env=env,stdout=log,stderr=subprocess.STDOUT,timeout=180)
    assert proc.returncode==0
    finish(name)


def finish(name):
    surface.clean_sha()
    spec=study.FRAMES[name];root=study.OUT/'frames'/name
    receipt=study.read(root/'bounded_reconstruction/receipt.json')
    p=np.load(root/'positions.npy');t=np.load(root/'indices.npy').reshape(-1,3)
    n=np.load(root/'normals.npy');assert n.shape==p.shape and np.isfinite(n).all()
    bad=np.flatnonzero(np.abs(np.linalg.norm(n,axis=1)-1)>1e-4)
    if len(bad):
        # Remove only a wholly disconnected, numerically collapsed component.
        faces=np.isin(t,bad).any(1)
        assert np.isin(t[faces],bad).all()
        diameter=float(np.linalg.norm(np.ptp(p[bad],axis=0))*.25)
        assert diameter<1e-5 and len(bad)<=8
        base.save(root/'collapsed_component.json',dict(vertices=bad.tolist(),faces=int(faces.sum()),diameter_m=diameter,
                  reason='Zero-area disconnected component has undefined normals; original arrays retained.'))
        keep=np.ones(len(p),bool);keep[bad]=False
        remap=np.cumsum(keep)-1
        p,n,t=p[keep],n[keep],remap[t[~faces]]
    assert np.max(np.abs(np.linalg.norm(n,axis=1)-1))<1e-4
    np.savez_compressed(root/'mesh.npz',points=p,triangles=t,normals=n)
    hashes={k:base.sha(spec['source']/f"{k}.{spec['frame']:04}.npy") for k in ['liquid','radius']}
    count=len(np.load(root/'bounded_cache/liquid.0001.npy',mmap_mode='r'))
    for k,h in hashes.items():assert base.sha(spec['source']/f"{k}.{spec['frame']:04}.npy")==h
    receipt.update(vertices=len(p),triangles=len(t),collapsed_vertices_removed=len(bad))
    base.save(root/'mesh_receipt.json',dict(receipt,time_s=spec['time_s'],source=str(spec['source']),
              source_frame=spec['frame'],source_hashes=hashes,crop_x_min_m=spec['xmin'],crop_z_m=[10,78],
              particles=count,mesh_sha256=base.sha(root/'mesh.npz'),simulation=False,
              bounded_mask_only=True,interpolation_geometry=False))
    print(study.read(root/'mesh_receipt.json'))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['control','pre','break','shore','finish_shore'])
    a=parser.parse_args().action
    if a=='control':control()
    elif a=='finish_shore':finish('shore')
    else:frame(a)
