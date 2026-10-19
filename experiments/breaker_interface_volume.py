"""Bake a curvature-gated, 10 cm inner interface band to a reusable VDB.

Stock Bifrost builds/exports the sparse grid. Public Maya mesh queries locate
each voxel against the frozen surface. This is render preparation, not a solve.
"""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import traceback

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import breaker_deep_water as base
from breaker_interface_mask import DEFAULT, smoothstep

BIF=Path('C:/Program Files/Autodesk/Bifrost/Maya2027/3.1.0.8/bifrost')
MAYAPY='C:/Program Files/Autodesk/Maya2027/bin/mayapy.exe'
VOXEL=.125  # solver-cell world units = 0.03125 m


def graph(out,export=False):
    nodes=[]; conns=[]; values=[];ports=[]
    def node(name,kind):nodes.append(dict(nodeName=name,nodeType=kind))
    def value(name,typ,val):
        nodes.append(dict(nodeName=name,valueType=typ))
        values.append(dict(valueName=name+'.value',valueType=typ,value=val))
        return name+'.output'
    def setv(name,typ,val):values.append(dict(valueName=name,valueType=typ,value=val))
    def connect(a,b):conns.append(dict(source=a,target=b))
    def writer(name,source,path):
        node(name,'File::NumPy::write_NumPy');connect(source,name+'.data')
        setv(name+'.file_path','string',path.as_posix());setv(name+'.frame','long','0L')
        setv(name+'.enable','bool','true');setv(name+'.overwrite','bool','false')
        ports.append(dict(portName='ok_'+name,portDirection='output',portType='bool'))
        connect(name+'.success','.ok_'+name)
    f3=value('f3','array<Math::float3>',[{'x':'0f','y':'0f','z':'0f'}])
    flt=value('flt','array<float>',['0f']);zero=value('zero','float','0f')
    node('read','File::NumPy::read_NumPy');connect(f3,'read.type')
    setv('read.file_path','string',(out/'support.####').as_posix());setv('read.frame','int','0L')
    node('pts','Geometry::Points::construct_points');connect('read.data','pts.point_position')
    node('size','Geometry::Properties::set_geo_property');connect('pts.points','size.geometry')
    setv('size.property','string','point_size');setv('size.data','float','1.5f');setv('size.default','float','1.5f')
    node('vol','Geometry::Converters::points_to_volume');connect('size.out_geometry','vol.points')
    setv('vol.detail_size','float',str(VOXEL)+'f')
    setv('vol.store_level_set','bool','false');setv('vol.store_fog_density','bool','true')
    setv('vol.fog_density','float','1f');setv('vol.rebuild_level_set','bool','false')
    node('vp','Geometry::Properties::get_voxel_position');connect('vol.volume','vp.volume')
    writer('positions','vp.voxel_position',out/('export_voxels.npy' if export else 'voxels.npy'))
    if not export:
        node('get','Geometry::Properties::get_geo_property');connect('vol.volume','get.geometry')
        connect(flt,'get.type');setv('get.property','string','voxel_fog_density')
        writer('support_density','get.data',out/'support_density.npy')
    else:
        node('density','File::NumPy::read_NumPy');connect(flt,'density.type')
        setv('density.file_path','string',(out/'extinction.####').as_posix());setv('density.frame','int','0L')
        node('set','Geometry::Properties::set_geo_property_data');connect('vol.volume','set.geometry')
        setv('set.property','string','voxel_fog_density')
        connect('density.data','set.data')
        ports.append(dict(portName='ok_set',portDirection='output',portType='bool'))
        connect('set.success','.ok_set')
        node('verify','Geometry::Properties::get_geo_property');connect('set.out_geometry','verify.geometry')
        connect(flt,'verify.type');setv('verify.property','string','voxel_fog_density')
        writer('density_check','verify.data',out/'export_density.npy')
        nodes.append(dict(nodeName='volumes',nodeType='Core::Array::build_array',multiInPortNames=['v']))
        connect('set.out_geometry','volumes.first.v')
        node('write','File::Geometry::write_OpenVDB_volume')
        connect('volumes.array','write.volumes')
        setv('write.filename','string','interface.vdb');setv('write.directory','string',out.as_posix())
        setv('write.enable','bool','true');setv('write.overwrite','bool','false')
        setv('write.properties','string','voxel_fog_density')
        ports.append(dict(portName='ok_vdb',portDirection='output',portType='bool'))
        connect('write.success','.ok_vdb')
    return dict(header=dict(metadata=[dict(metaName='adskFileFormatVersion',metaValue='100L')]),
                namespaces=[],types=[],compounds=[dict(name='User::PFFlip::interface_grid',ports=ports,
                    compoundNodes=nodes,connections=conns,values=values)])


def sample(out):
    import maya.standalone
    maya.standalone.initialize(name='python')
    import maya.cmds as c
    import maya.api.OpenMaya as om
    start=time.perf_counter()
    try:
        c.file(str(base.SOURCE),open=True,force=True)
        selection=om.MSelectionList();selection.add('baked_control')
        dag=selection.getDagPath(0);dag.extendToShape();mesh=om.MFnMesh(dag)
        inter=om.MMeshIntersector();inter.create(mesh.object())
        data=np.load(out/'mask.npz');p=data['points_m']*4;t=data['triangles'];mask=data['mask']
        assert np.allclose(np.asarray(mesh.getPoints())[:,:3],p,rtol=0,atol=1e-7)
        assert np.array_equal(np.asarray(mesh.getVertices()[1]).reshape(-1,3),t)
        vox=np.load(out/'voxels.npy'); support=np.load(out/'support_density.npy')
        assert len(vox)==len(support)
        ids=np.flatnonzero(support>1e-6)
        assert len(ids)<3000000,('Unexpectedly large query set',len(ids))
        closest=np.zeros((len(ids),3)); normal=np.zeros_like(closest);faces=np.empty(len(ids),np.int32)
        for i,idx in enumerate(ids):
            hit=inter.getClosestPoint(om.MPoint(*map(float,vox[idx])))
            closest[i]=tuple(hit.point)[:3];normal[i]=tuple(hit.normal)[:3];faces[i]=hit.face
            if i%100000==0:print('closest',i,'/',len(ids),flush=True)
        delta=vox[ids]-closest
        signed=np.linalg.norm(delta,axis=1)*np.where(np.einsum('ij,ij->i',delta,normal)<0,-1,1)
        depth=-signed*.25
        tri=p[t[faces]];v0=tri[:,1]-tri[:,0];v1=tri[:,2]-tri[:,0];v2=closest-tri[:,0]
        d00=np.einsum('ij,ij->i',v0,v0);d01=np.einsum('ij,ij->i',v0,v1)
        d11=np.einsum('ij,ij->i',v1,v1);d20=np.einsum('ij,ij->i',v2,v0);d21=np.einsum('ij,ij->i',v2,v1)
        den=d00*d11-d01*d01
        b=(d11*d20-d01*d21)/den;cc=(d00*d21-d01*d20)/den
        weights=np.clip(np.column_stack([1-b-cc,b,cc]),0,1);weights/=weights.sum(axis=1)[:,None]
        w=np.sum(mask[t[faces]]*weights,axis=1)
        band=smoothstep(0,.03125,depth)*(1-smoothstep(.06875,.1,depth))
        density=np.zeros(len(vox),np.float32)
        density[ids]=(5*w*band).astype(np.float32)  # 20 /m peak extinction, 4 world units/m
        assert np.isfinite(density).all() and density.max()<=5
        assert np.all(density[ids[(depth<=0)|(depth>=.1)]]==0)
        active=density>1e-6
        assert active.any()
        # Direct voxel-center check in addition to the surface-mask check.
        xyz=vox[active]*.25;lo=data['points_m'].min(0);hi=data['points_m'].max(0)
        cut=np.minimum.reduce([xyz[:,0]-lo[0],hi[0]-xyz[:,0],xyz[:,2]-lo[2],hi[2]-xyz[:,2]])
        assert cut.min()>.45,cut.min()
        np.save(out/'extinction.0000.npy',density)
        np.savez_compressed(out/'volume_checks.npz',active_positions_cells=vox[active],density_per_cell=density[active])
        record=dict(status='passed',voxel_count=len(vox),queried_voxels=len(ids),active_voxels=int(active.sum()),
                    voxel_size_m=.03125,band_depth_m=[0,.1],peak_extinction_per_m=float(density.max()*4),
                    mean_free_path_at_peak_m=1/(float(density.max())*4),
                    min_active_distance_to_cut_m=float(cut.min()),
                    integrated_extinction_m2=float(density.sum(dtype=float)*4*.03125**3),
                    grid_bounds_cells=[vox.min(0).tolist(),vox.max(0).tolist()],
                    active_bounds_cells=[vox[active].min(0).tolist(),vox[active].max(0).tolist()],
                    wall_s=time.perf_counter()-start,geometry=base.geometry(c))
        base.save(out/'sample_receipt.json',record);base.exit_worker(0)
    except BaseException:
        base.save(out/'sample_receipt.json',dict(status='failed',error=traceback.format_exc()))
        base.exit_worker(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['grid','sample','export'])
    parser.add_argument('--out',type=Path,default=DEFAULT);args=parser.parse_args();out=args.out.resolve()
    assert not base.git('status','--porcelain');assert base.sha(base.SOURCE)==base.SOURCE_SHA
    assert json.loads((out/'mask_review.json').read_text())['status']=='passed'
    mode=args.mode;receipt=dict(status='running',mode=mode,git_sha=base.git('rev-parse','HEAD'))
    env=os.environ.copy();env.update(MAYA_DISABLE_CER='1',MAYA_DISABLE_CIP='1',MAYA_DISABLE_CLIC_IPM='1')
    if mode=='sample':
        assert not (out/'extinction.0000.npy').exists()
        env.update(INTERFACE_SAMPLE_WORKER='1',INTERFACE_OUT=str(out))
        command=[MAYAPY,str(Path(__file__).resolve())]
    else:
        assert not (out/('voxels.npy' if mode=='grid' else 'interface.vdb')).exists()
        g=out/(mode+'.json');base.save(g,graph(out,mode=='export'))
        env['PATH']=str(BIF/'bin')+';'+str(BIF/'thirdparty/bin')+';'+env['PATH']
        command=[str(BIF/'bin/bifcmd.exe'),str(g)]
    start=time.perf_counter()
    try:
        with (out/(mode+'.log')).open('w') as log:
            proc=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=env,timeout=600)
        receipt['exit_code']=proc.returncode
        assert proc.returncode==0,'Read '+mode+'.log'
        if mode=='export':
            assert np.array_equal(np.load(out/'voxels.npy'),np.load(out/'export_voxels.npy')),'Grid ordering changed'
            assert (out/'interface.vdb').stat().st_size>1024
            receipt.update(vdb_sha256=base.sha(out/'interface.vdb'),vdb_bytes=(out/'interface.vdb').stat().st_size,
                           grid_order_exact=True)
        receipt['status']='passed'
    except BaseException:
        receipt.update(status='failed',error=traceback.format_exc());raise
    finally:
        receipt['wall_s']=time.perf_counter()-start;base.save(out/(mode+'_process.json'),receipt)
        print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    if os.environ.get('INTERFACE_SAMPLE_WORKER')=='1':sample(Path(os.environ['INTERFACE_OUT']))
    else:main()
