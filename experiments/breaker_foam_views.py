"""Paired A/B viewing study requested by Xue; existing looks, cached particles only.

Preparation is offline. Rendering runs only through licensed GUI Maya, one image
at a time, followed by the established EXR/AOV/color checks.
"""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback

import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import breaker_deep_water as base
import breaker_foam_render as render
from breaker_foam_layer import check_imported, ass_units

ROOT=base.ROOT
OUT=ROOT/'out/runs/breaker_foam_views_20261019'
ORIGINAL=ROOT/'out/runs/breaker_foam_layer_20261019_v6'
FRAMES=OUT/'frames'
VIEWS={
    'wide':dict(time_s=24.9,width_m=24,rotation=[0,168,0],target_m=[145.5,13.5,70]),
    'oblique':dict(time_s=24.9,width_m=18,rotation=[-10,145,0],target_m=[145.5,13.5,70]),
    'later26':dict(time_s=26.025,width_m=18,rotation=[0,168,0],target_m=[150,13,70]),
    'later28':dict(time_s=28.025,width_m=18,rotation=[0,168,0],target_m=[156,12,70]),
}


def prepare_worker(view):
    import maya.standalone
    maya.standalone.initialize(name='python')
    import maya.cmds as c
    started=time.perf_counter();out=OUT/view;spec=VIEWS[view]
    try:
        for arm in ['A','B']:
            folder=out/arm
            source=ORIGINAL/arm
            c.file(str(source/'candidate.mb'),open=True,force=True)
            before=json.loads((source/'state.json').read_text())
            assert base.geometry(c)==before['geometry']
            assert not c.ls(type='bifrostGraphShape')
            if spec['time_s']!=24.9:
                import maya.api.OpenMaya as om
                z=np.load(FRAMES/view/'mesh.npz');points=z['points'];triangles=z['triangles'];normals=z['normals']
                c.delete('baked_control')
                transform=c.createNode('transform',name='baked_control')
                selection=om.MSelectionList();selection.add(transform);parent=selection.getDependNode(0)
                fn=om.MFnMesh();fn.create(om.MPointArray(points.tolist()),[3]*len(triangles),
                                         triangles.ravel().tolist(),parent=parent)
                fn.setVertexNormals(om.MVectorArray(normals.tolist()),list(range(len(points))))
                c.sets(transform,e=True,forceElement='waterSG')
                normal_counts,normal_ids=fn.getNormalIds()
                corners=np.asarray(fn.getNormals())[np.asarray(normal_ids)].reshape(-1,3,3)
                error=float(np.max(np.abs(corners-normals[triangles])))
                assert error<1e-5,error
                assert np.array_equal(np.asarray(fn.getPoints())[:,:3].astype(np.float32),points)
                base.save(folder/'normal_transfer.json',dict(status='passed',max_corner_normal_error=error,
                          mesh_sha256=base.sha(FRAMES/view/'mesh.npz')))
            pitch,yaw,_=np.radians(spec['rotation'])
            direction=np.array([np.sin(yaw)*np.cos(pitch),-np.sin(pitch),np.cos(yaw)*np.cos(pitch)])
            eye=np.asarray(spec['target_m'])*4+300*direction
            c.setAttr('persp.translate',*eye,type='double3')
            c.setAttr('persp.rotate',*spec['rotation'],type='double3')
            c.setAttr('perspShape.orthographic',True)
            c.setAttr('perspShape.orthographicWidth',spec['width_m']*4)
            c.setAttr('foamLayerShape.dso',(folder/'foam.ass').as_posix(),type='string')
            c.setAttr('opticalEXR.prefix',(folder/'optics').as_posix(),type='string')
            c.setAttr('defaultArnoldRenderOptions.log_filename',(folder/'arnold.log').as_posix(),type='string')
            state={node:{a:c.getAttr(node+'.'+a) for a in before[node]} for node in render.NODE_KEYS}
            state=json.loads(json.dumps(state))
            for node in render.NODE_KEYS:
                if node not in ['persp','perspShape']:
                    assert state[node]==before[node],node
            state.update(geometry=base.geometry(c),foam_ass_sha256=base.sha(folder/'foam.ass'),foam_arm=arm,
                         foam_layer_sha256=base.sha(out/'layer.npz'),mist=False,view=spec)
            base.save(folder/'state.json',state)
            c.file(rename=str(folder/'candidate.mb'));c.file(save=True,type='mayaBinary',force=True)
            c.arnoldExportAss(f=str(folder/'scene_check.ass'),mask=24)
            assert ass_units(folder/'scene_check.ass')==ass_units(folder/'foam.ass')==.01
            base.save(folder/'export_check.json',dict(status='passed',scene_sha256=base.sha(folder/'candidate.mb'),
                      ass_sha256=base.sha(folder/'scene_check.ass')))
        base.save(out/'prepare_worker.json',dict(status='passed',wall_s=time.perf_counter()-started))
        base.exit_worker(0)
    except BaseException:
        base.save(out/'prepare_worker.json',dict(status='failed',error=traceback.format_exc()));base.exit_worker(1)


def prepare(view):
    assert not base.git('status','--porcelain')
    out=OUT/view;assert not out.exists();out.mkdir(parents=True)
    source=ORIGINAL if VIEWS[view]['time_s']==24.9 else FRAMES/view/'foam'
    if VIEWS[view]['time_s']==24.9:
        assert json.loads((source/'analysis.json').read_text())['status']=='passed'
        source_receipt=source/'analysis.json'
    else:
        assert json.loads((source/'review.json').read_text())['status']=='passed'
        assert json.loads((FRAMES/view/'mesh_receipt.json').read_text())['status']=='passed'
        source_receipt=source/'layer_check.json'
    shutil.copy2(source/'layer.npz',out/'layer.npz')
    for arm in ['A','B']:
        folder=out/arm;folder.mkdir();shutil.copy2(source/arm/'foam.ass',folder/'foam.ass')
        assert base.sha(folder/'foam.ass')==base.sha(source/arm/'foam.ass')
    env=os.environ.copy();env.update(FOAM_VIEW_WORKER='prepare',FOAM_VIEW_NAME=view)
    with (out/'prepare.log').open('w') as log:
        result=subprocess.run([render.MAYAPY,str(Path(__file__).resolve())],cwd=ROOT,env=env,
                              stdout=log,stderr=subprocess.STDOUT,timeout=180)
    assert result.returncode==0,(out/'prepare.log').read_text()[-3000:]
    for arm in ['A','B']:check_imported(out/arm)
    base.save(out/'prepare.json',dict(status='passed',git_sha=base.git('rev-parse','HEAD'),view=VIEWS[view],
              source=str(source),source_receipt_sha256=base.sha(source_receipt),simulation=False))


def launch(view,arm):
    assert not base.git('status','--porcelain')
    out=OUT/view;folder=out/arm
    assert json.loads((out/'prepare.json').read_text())['status']=='passed'
    alignment=json.loads((folder/'alignment_check.json').read_text());assert alignment['status']=='passed'
    assert alignment['foam_ass_sha256']==base.sha(folder/'foam.ass')
    assert json.loads((folder/'export_check.json').read_text())['scene_sha256']==base.sha(folder/'candidate.mb')
    if arm=='B':assert json.loads((out/'A/image_check.json').read_text())['status']=='passed'
    path=folder/'receipt.json';assert not path.exists()
    # Use the already verified launch quoting from the preceding valid render.
    launcher=(ORIGINAL/arm/'launch.ps1').read_text().replace('FOAM_RENDER_SCRIPT','FOAM_VIEW_SCRIPT')
    (folder/'launch.ps1').write_text(launcher)
    env=os.environ.copy();env.update(FOAM_VIEW_WORKER='render',FOAM_VIEW_NAME=view,FOAM_VIEW_ARM=arm,
           FOAM_VIEW_SCRIPT=str(Path(__file__).resolve()),MAYA_DISABLE_CER='1',MAYA_DISABLE_CIP='1',MAYA_DISABLE_CLIC_IPM='1')
    receipt=dict(status='running',git_sha=base.git('rev-parse','HEAD'),view=VIEWS[view],arm=arm,
                 started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    base.save(path,receipt);started=time.perf_counter()
    try:
        with (folder/'launch.log').open('w') as log:
            result=subprocess.run(['powershell.exe','-NoProfile','-ExecutionPolicy','Bypass','-File',str(folder/'launch.ps1')],
                                  cwd=ROOT,env=env,stdout=log,stderr=subprocess.STDOUT,timeout=930)
        receipt['exit_code']=result.returncode;assert result.returncode==0
        worker=json.loads((folder/'worker_result.json').read_text());assert worker['status']=='passed'
        receipt.update(status='passed',worker=worker,png_sha256=base.sha(folder/'beauty.png'),
                       exr_sha256=base.sha(folder/worker['exr']))
    except BaseException:
        receipt.update(status='failed',error=traceback.format_exc());raise
    finally:
        receipt['process_s']=time.perf_counter()-started;base.save(path,receipt)
    render.check(folder)


def report():
    from PIL import Image,ImageDraw,ImageFont
    font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',22)
    small=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',16)
    results={};contact=Image.new('RGB',(960,4*310),(20,23,28));draw=ImageDraw.Draw(contact)
    labels={'wide':'24.9 s | wider side | 24 m view',
            'oblique':'24.9 s | elevated oblique | 18 m view',
            'later26':'26.025 s | side | 18 m view','later28':'28.025 s | side | 18 m view'}
    for row,view in enumerate(VIEWS):
        out=OUT/view;states=[];checks={}
        pair=Image.new('RGB',(1920,590),(20,23,28));d=ImageDraw.Draw(pair)
        for col,arm in enumerate(['A','B']):
            folder=out/arm;checks[arm]=json.loads((folder/'image_check.json').read_text())
            assert checks[arm]['status']=='passed'
            states.append(json.loads((folder/'state.json').read_text()))
            image=Image.open(folder/'beauty_view.png').convert('RGB')
            d.text((col*960+10,8),labels[view]+' | '+('A uniform thin surface' if arm=='A' else 'B clustered foam'),font=font,fill='white')
            pair.paste(image,(col*960,50));contact.paste(image.resize((480,270),Image.Resampling.LANCZOS),(col*480,row*310+40))
        differing=[k for k in states[0] if states[0][k]!=states[1][k]]
        assert set(differing)=={'foam_ass_sha256','foam_arm'},differing
        pair.save(out/'comparison.png')
        draw.text((8,row*310+4),labels[view]+'     A: uniform (left) / B: clustered (right)',font=small,fill='white')
        results[view]=dict(settings=VIEWS[view],images=checks,changed_state_keys=differing,
                           geometry=states[0]['geometry'])
    contact.save(OUT/'comparison.png')
    base.save(OUT/'analysis.json',dict(status='passed',render_count=8,simulation=False,mist=False,
              views=results,geometry_audit_candidates_used=False,
              foam_note='Per-frame static curvature mask and existing A/B appearance; no persistent foam simulation.',
              startup_failure='One owned GUI launch stopped before scene loading because runpy lacked the script import path; preserved under wide/A/failed_launch_01.'))


def show_plan():
    base.save(OUT/'plan.json',dict(views=VIEWS,render_count=8,resolution=[960,540],
              user_request='Paired A/B views at different angles, distances and cached times after 24.9 s.',
              estimated_total_minutes=[40,60],notify_before_exceeding_minutes=60,
              simulation=False,mist=False,geometry_audit_candidates_used=False))


if __name__=='__main__':
    if os.environ.get('FOAM_VIEW_WORKER'):
        view=os.environ['FOAM_VIEW_NAME']
        if os.environ['FOAM_VIEW_WORKER']=='prepare':prepare_worker(view)
        else:render.render(OUT/view/os.environ['FOAM_VIEW_ARM'])
    else:
        parser=argparse.ArgumentParser(description=__doc__)
        parser.add_argument('action',choices=['prepare','A','B','report','plan']);parser.add_argument('view',choices=list(VIEWS),nargs='?')
        args=parser.parse_args()
        if args.action=='report':report()
        elif args.action=='plan':show_plan()
        else:
            assert args.view
            prepare(args.view) if args.action=='prepare' else launch(args.view,args.action)
