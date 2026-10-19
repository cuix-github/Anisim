"""Two authorized GUI-Maya stills of static uniform/clustered surface foam.

Preparation and EXR analysis never trace rays. Mist and solver remain disabled.
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
from analyze_breaker_deep_water import pfm,convert,OIIO,MAYAPY,COMPONENTS,LAYERS
from breaker_foam_layer import ass_units,RADIUS,OFFSET

OUT=base.ROOT/'out/runs/breaker_foam_layer_20261019_v2'
CONTROL=base.ROOT/'out/runs/breaker_interface_sss_20261019_v3/A'
NODE_KEYS=['water','sky','sun','sunShape','persp','perspShape','defaultArnoldDriver',
           'defaultArnoldRenderOptions','defaultColorMgtGlobals']


def prepare(out):
    import maya.standalone
    maya.standalone.initialize(name='python')
    import maya.cmds as c
    start=time.perf_counter()
    try:
        c.file(str(base.SOURCE),open=True,force=True)
        assert not c.ls(type='bifrostGraphShape')
        baseline=json.loads((CONTROL/'state.json').read_text())
        for node in ['persp','perspShape','defaultArnoldDriver']:
            for attr,value in baseline[node].items():
                if isinstance(value,list):c.setAttr(node+'.'+attr,*value[0],type='double3')
                else:c.setAttr(node+'.'+attr,value)
        for attr,val in [('width',960),('height',540),('deviceAspectRatio',960/540)]:c.setAttr('defaultResolution.'+attr,val)
        assert c.listConnections('waterSG.surfaceShader')==['water']
        assert base.geometry(c)==baseline['geometry']
        standin=c.createNode('aiStandIn',name='foamLayerShape')
        if c.attributeQuery('overrideShaders',node=standin,exists=True):c.setAttr(standin+'.overrideShaders',False)
        c.setAttr('defaultRenderGlobals.imageFilePrefix','beauty',type='string')
        for arm in ['A','B']:
            folder=out/arm
            c.setAttr(standin+'.dso',(folder/'foam.ass').as_posix(),type='string')
            c.setAttr('opticalEXR.prefix',(folder/'optics').as_posix(),type='string')
            c.setAttr('defaultArnoldRenderOptions.log_filename',(folder/'arnold.log').as_posix(),type='string')
            state={node:{attr:c.getAttr(node+'.'+attr) for attr in baseline[node]} for node in NODE_KEYS}
            state=json.loads(json.dumps(state))  # Maya tuples vs persisted JSON lists
            for node in NODE_KEYS:assert state[node]==baseline[node],node
            state.update(geometry=base.geometry(c),mist=False,foam_ass_sha256=base.sha(folder/'foam.ass'),
                         foam_layer_sha256=base.sha(out/'layer.npz'),foam_arm=arm,control_receipt_sha256=base.sha(CONTROL/'receipt.json'))
            base.save(folder/'state.json',state)
            c.file(rename=str(folder/'candidate.mb'));c.file(save=True,type='mayaBinary',force=True)
            c.arnoldExportAss(f=str(folder/'scene_check.ass'),mask=24)
            exported=(folder/'scene_check.ass').read_text()
            assert 'foam.ass' in exported and 'foamLayer' in exported
            child_unit=ass_units(folder/'foam.ass');parent_unit=ass_units(folder/'scene_check.ass')
            assert child_unit==parent_unit,(child_unit,parent_unit)
            transform=c.listRelatives(standin,parent=True,fullPath=True)[0]
            matrix=np.asarray(c.xform(transform,q=True,worldSpace=True,matrix=True)).reshape(4,4)
            assert np.allclose(matrix,np.eye(4),atol=1e-8),matrix.tolist()
            z=np.load(out/'layer.npz');active=z['density_'+arm+'_per_m']>0
            # These point/radius values were checked against native ASS readback.
            points=(z['points_m'][active]*4).astype(np.float32);radius=np.float32(RADIUS*4)
            expected=np.r_[(points-radius).min(0),(points+radius).max(0)]
            actual=np.asarray(c.exactWorldBoundingBox(standin),float)
            error=float(np.abs(actual-expected).max())
            assert error<.01,('Imported stand-in bounds differ',actual.tolist(),expected.tolist())
            imported_m=points.astype(float)*.25
            offsets=np.linalg.norm(imported_m-z['surface_m'][active],axis=1)
            assert np.max(np.abs(offsets-OFFSET))<1e-5
            assert offsets.max()+RADIUS<.03001
            alignment=dict(status='passed',child_meters_per_unit=child_unit,parent_meters_per_unit=parent_unit,
                           import_scale=child_unit/parent_unit,maya_world_bounds=actual.tolist(),expected_bounds=expected.tolist(),
                           bbox_max_error_scene_units=error,center_to_source_surface_max_m=float(offsets.max()),
                           sphere_to_source_surface_upper_bound_m=float(offsets.max()+RADIUS),matrix=matrix.tolist(),
                           layer_sha256=base.sha(out/'layer.npz'),foam_ass_sha256=base.sha(folder/'foam.ass'))
            base.save(folder/'alignment_check.json',alignment)
            base.save(folder/'export_check.json',dict(status='passed',scene_sha256=base.sha(folder/'candidate.mb'),
                                                      ass_sha256=base.sha(folder/'scene_check.ass')))
        base.save(out/'prepare.json',dict(status='passed',wall_s=time.perf_counter()-start));base.exit_worker(0)
    except BaseException:
        base.save(out/'prepare.json',dict(status='failed',error=traceback.format_exc()));base.exit_worker(1)


def render(folder):
    import maya.cmds as c
    start=time.perf_counter()
    try:
        c.scriptEditorInfo(writeHistory=True,historyFilename=str(folder/'maya.log'))
        c.file(str(folder/'candidate.mb'),open=True,force=True)
        c.workspace(fileRule=['images',folder.as_posix()])
        assert not c.ls(type='bifrostGraphShape')
        assert c.getAttr('baked_control.visibility')
        assert not c.getAttr('baked_raw.visibility') and not c.getAttr('baked_no_mcf.visibility')
        state=json.loads((folder/'state.json').read_text());assert base.geometry(c)==state['geometry']
        assert c.listConnections('waterSG.surfaceShader')==['water']
        assert base.sha(Path(c.getAttr('foamLayerShape.dso')))==state['foam_ass_sha256']
        renderstart=time.perf_counter();c.arnoldRender(seq='1',w=960,h=540,cam='persp',srv=False)
        render_s=time.perf_counter()-renderstart
        pngs=list((folder/'tmp').rglob('*.png'));assert len(pngs)==1,pngs
        shutil.move(str(pngs[0]),folder/'beauty.png')
        exrs=list(folder.glob('*.exr'));assert len(exrs)==1,exrs
        base.save(folder/'worker_result.json',dict(status='passed',render_call_s=render_s,worker_s=time.perf_counter()-start,exr=exrs[0].name))
        base.exit_worker(0)
    except BaseException:
        base.save(folder/'worker_result.json',dict(status='failed',error=traceback.format_exc()));base.exit_worker(1)


def launch(out,mode):
    assert not base.git('status','--porcelain'),'Commit instrument first'
    assert base.sha(base.SOURCE)==base.SOURCE_SHA
    review=json.loads((out/'review.json').read_text());assert review['status']=='passed'
    layer=json.loads((out/'layer_check.json').read_text());assert layer['layer_sha256']==base.sha(out/'layer.npz')
    hashes={p.name:base.sha(p) for p in base.CACHE.iterdir() if p.is_file()}
    assert hashes==json.loads((CONTROL/'receipt.json').read_text())['cache_sha256']
    assert json.loads((CONTROL/'control_check.json').read_text())['status']=='passed'
    control_receipt=json.loads((CONTROL/'receipt.json').read_text())
    assert base.sha(CONTROL/'beauty.png')==control_receipt['png_sha256']
    assert base.sha(CONTROL/control_receipt['worker']['exr'])==control_receipt['exr_sha256']
    env={k:v for k,v in os.environ.items() if not k.startswith('SURF_')}
    env.update(FOAM_RENDER_MODE=mode,FOAM_RENDER_OUT=str(out),FOAM_RENDER_SCRIPT=str(Path(__file__).resolve()),
               MAYA_DISABLE_CER='1',MAYA_DISABLE_CIP='1',MAYA_DISABLE_CLIC_IPM='1')
    folder=out if mode=='prepare' else out/mode
    path=folder/('prepare_process.json' if mode=='prepare' else 'receipt.json');assert not path.exists()
    if mode=='prepare':
        command=[MAYAPY,str(Path(__file__).resolve())];timeout=180
    else:
        assert json.loads((out/'prepare.json').read_text())['status']=='passed'
        alignment=json.loads((folder/'alignment_check.json').read_text());assert alignment['status']=='passed'
        assert alignment['foam_ass_sha256']==base.sha(folder/'foam.ass')
        assert json.loads((folder/'export_check.json').read_text())['scene_sha256']==base.sha(folder/'candidate.mb')
        assert layer['exports'][mode]['ass_sha256']==base.sha(folder/'foam.ass')
        if mode=='B':assert json.loads((out/'A/image_check.json').read_text())['status']=='passed'
        launcher=r"""
$ErrorActionPreference = 'Stop'
if (Get-Process maya,mayapy -ErrorAction SilentlyContinue) { throw 'Another Maya process is open' }
if ((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB -lt 20) { throw 'Less than 20 GB RAM free' }
$pythonScript = $env:FOAM_RENDER_SCRIPT -replace '\\', '/'
$command = '"python(\"import runpy; runpy.run_path(r''' + $pythonScript + ''', run_name=''__main__'')\")"'
$p = Start-Process -FilePath 'C:/Program Files/Autodesk/Maya2027/bin/maya.exe' -ArgumentList @('-command',$command) -PassThru -WindowStyle Hidden
Write-Output "MAYA_PID $($p.Id)"
if (-not $p.WaitForExit(900000)) { Stop-Process -Id $p.Id -Force; throw 'Owned render timed out after 900 seconds' }
exit $p.ExitCode
"""
        (folder/'launch.ps1').write_text(launcher,encoding='utf-8')
        command=['powershell.exe','-NoProfile','-ExecutionPolicy','Bypass','-File',str(folder/'launch.ps1')];timeout=930
    receipt=dict(status='running',mode=mode,git_sha=base.git('rev-parse','HEAD'),cache_sha256=hashes,
                 source_sha256=base.SOURCE_SHA,control_receipt_sha256=base.sha(CONTROL/'receipt.json'),
                 started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    base.save(path,receipt);start=time.perf_counter()
    try:
        with (folder/(mode+'_launch.log')).open('w') as log:
            proc=subprocess.run(command,env=env,cwd=base.ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
        receipt['exit_code']=proc.returncode;assert proc.returncode==0,'Read worker log'
        assert hashes=={p.name:base.sha(p) for p in base.CACHE.iterdir() if p.is_file()}
        if mode!='prepare':
            worker=json.loads((folder/'worker_result.json').read_text());assert worker['status']=='passed'
            receipt.update(worker=worker,png_sha256=base.sha(folder/'beauty.png'),exr_sha256=base.sha(folder/worker['exr']),
                           foam_ass_sha256=base.sha(folder/'foam.ass'))
        receipt['status']='passed'
    except BaseException:
        receipt.update(status='failed',error=traceback.format_exc());raise
    finally:
        receipt['process_s']=time.perf_counter()-start;base.save(path,receipt);print(json.dumps(receipt,indent=2))


def check(folder):
    from PIL import Image
    receipt=json.loads((folder/'receipt.json').read_text());assert receipt['status']=='passed'
    assert base.sha(folder/'beauty.png')==receipt['png_sha256']
    exr=folder/receipt['worker']['exr'];assert base.sha(exr)==receipt['exr_sha256']
    info=subprocess.check_output([OIIO,'--info','-v',str(exr)],text=True)
    (folder/'exr_info.txt').write_text(info)
    assert '960 x  540, 54 channel, float openexr' in info and 'arnold/license_state: "valid"' in info
    log=(folder/'arnold.1.log').read_text();assert 'ERROR' not in log and '| render done' in log
    data={}
    for layer in LAYERS:
        ch='R,G,B' if layer=='beauty' else ','.join(layer+'.'+s for s in ('XYZ' if layer in ['N','P'] else 'RGB'))
        target=folder/(layer+'.pfm');subprocess.run([OIIO,str(exr),'--ch',ch,'-o',str(target)],check=True)
        data[layer]=pfm(target)
    subprocess.run([MAYAPY,str(Path(__file__).resolve()),'convert','--out',str(folder)],check=True)
    total=data['beauty_raw'].sum(dtype=float)
    closure=float(np.abs(sum(data[k].astype(float) for k in COMPONENTS)-data['beauty_raw']).sum()/total)
    assert closure<1e-5,closure
    png=np.asarray(Image.open(folder/'beauty.png').convert('RGB'),float)
    view=np.asarray(Image.open(folder/'beauty_view.png').convert('RGB'),float);delta=np.abs(png-view)
    assert delta.mean()<=1 and (delta.max(2)>8).mean()<.01
    shares={k:float(data[k].sum(dtype=float)/total) for k in COMPONENTS}
    result=dict(status='passed',aov_closure_relative_l1=closure,component_raw_rgb_share=shares,
                display_mean_abs_levels=float(delta.mean()),render_call_s=receipt['worker']['render_call_s'],process_s=receipt['process_s'])
    if shares['volume']<=1e-5:result.update(status='failed',error='Foam volume did not contribute')
    if shares['sss']!=0:result.update(status='failed',error='Unexpected SSS contribution')
    base.save(folder/'image_check.json',result);print(json.dumps(result,indent=2))
    assert result['status']=='passed',result


def report(out):
    from PIL import Image,ImageDraw,ImageFont
    checks={a:json.loads((out/a/'image_check.json').read_text()) for a in ['A','B']};assert all(r['status']=='passed' for r in checks.values())
    states=[json.loads((out/a/'state.json').read_text()) for a in ['A','B']]
    diffs=[k for k in states[0] if states[0][k]!=states[1][k]];assert set(diffs)=={'foam_arm','foam_ass_sha256'},diffs
    imageA=pfm(out/'A/beauty_raw.pfm');imageB=pfm(out/'B/beauty_raw.pfm');assert not np.array_equal(imageA,imageB)
    result=dict(status='passed',arms=checks,changed_state_keys=diffs,simulation_rerun=False,mesh_rebuilt=False,mist=False,
                layer=json.loads((out/'layer_check.json').read_text()),source_control=str(CONTROL),
                raw_beauty_mean_absolute_difference=float(np.abs(imageA-imageB).mean()))
    font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',20)
    sheet=Image.new('RGB',(960,1148),(20,23,28));draw=ImageDraw.Draw(sheet)
    for i,(arm,label) in enumerate([('A','A | Uniform thin surface foam / corrected baseline'),('B','B | Clustered foam / same total extinction budget')]):
        draw.text((12,i*574+6),label,font=font,fill='white');sheet.paste(Image.open(out/arm/'beauty_view.png').convert('RGB'),(0,i*574+34))
    sheet.save(out/'comparison.png')
    # A smaller display is a scale/readability aid, NOT a new camera rendering.
    sheet.resize((480,574),Image.Resampling.LANCZOS).save(out/'comparison_small.png')
    base.save(out/'analysis.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    if 'FOAM_RENDER_MODE' in os.environ:
        out=Path(os.environ['FOAM_RENDER_OUT']);mode=os.environ['FOAM_RENDER_MODE']
        if mode=='prepare':prepare(out)
        else:render(out/mode)
    else:
        p=argparse.ArgumentParser();p.add_argument('mode',choices=['prepare','A','B','checkA','checkB','report','convert']);p.add_argument('--out',type=Path,default=OUT)
        args=p.parse_args();out=args.out.resolve()
        if args.mode.startswith('check'):check(out/args.mode[-1])
        elif args.mode=='report':report(out)
        elif args.mode=='convert':convert(out)
        else:launch(out,args.mode)
