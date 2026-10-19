"""Two licensed GUI renders comparing baseline with gated subsurface scattering.

The separate VDB band failed its roundtrip gate. This is explicitly a different
render-only approximation, not a validated implementation of that volume band.
"""
import argparse
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback
sys.path.insert(0,str(Path(__file__).resolve().parent))
import breaker_deep_water as base

OUT=base.ROOT/'out/runs/breaker_interface_sss_20261019'
MASK=base.ROOT/'out/runs/breaker_interface_scatter_20261019'
MAYAPY='C:/Program Files/Autodesk/Maya2027/bin/mayapy.exe'


def prepare(out):
    import maya.standalone
    maya.standalone.initialize(name='python')
    import maya.cmds as c
    started=time.perf_counter()
    try:
        c.file(str(base.SOURCE),open=True,force=True)
        assert not c.ls(type='bifrostGraphShape')
        before=base.geometry(c)
        assert c.listConnections('waterSG.surfaceShader')==['water']
        prior=json.loads((base.PRIOR/'B/state.json').read_text())
        compare={'water':['base','metalness','specular','specularIOR','specularRoughness',
                          'transmission','transmissionDepth','transmissionColorR','transmissionColorG',
                          'transmissionColorB','transmissionScatterR','transmissionScatterG','transmissionScatterB'],
                 'sky':['intensity','colorR','colorG','colorB','camera'],
                 'sun':['visibility','rotateX','rotateY','rotateZ'],
                 'sunShape':['intensity','colorR','colorG','colorB','aiAngle']}
        for node,attrs in compare.items():
            for attr in attrs:
                assert c.getAttr(node+'.'+attr)==prior[node][attr],(node,attr)
        theta=math.radians(12);center=(582,54,280);forward=(-math.sin(theta),0,math.cos(theta))
        c.setAttr('persp.translate',*[center[i]-300*forward[i] for i in range(3)],type='double3')
        c.setAttr('persp.rotate',0,168,0,type='double3')
        c.setAttr('perspShape.orthographic',True);c.setAttr('perspShape.orthographicWidth',60)
        for attr,value in [('width',960),('height',540),('deviceAspectRatio',960/540)]:
            c.setAttr('defaultResolution.'+attr,value)
        modes=c.attributeQuery('colorManagement',node='defaultArnoldDriver',listEnum=True)[0].split(':')
        c.setAttr('defaultArnoldDriver.ai_translator','png',type='string')
        c.setAttr('defaultArnoldDriver.colorManagement',modes.index('Use View Transform'))
        c.setAttr('defaultRenderGlobals.imageFilePrefix','beauty',type='string')
        import numpy as np
        import maya.api.OpenMaya as om
        sel=om.MSelectionList();sel.add('baked_control');dag=sel.getDagPath(0);dag.extendToShape()
        mesh=om.MFnMesh(dag);weights=np.load(MASK/'mask.npz')['mask']
        c.polyColorSet('baked_control',create=True,colorSet='interfaceMask',representation='RGB')
        c.polyColorSet('baked_control',currentColorSet=True,colorSet='interfaceMask')
        colors=om.MColorArray([om.MColor((float(w),float(w),float(w),1)) for w in weights])
        mesh.setVertexColors(colors,list(range(mesh.numVertices)))
        c.setAttr(dag.fullPathName()+'.aiExportColors',True)
        read=c.shadingNode('aiUserDataColor',asUtility=True,name='interfaceMaskReader')
        c.setAttr(read+'.colorAttrName','interfaceMask',type='string')
        fog=c.shadingNode('aiStandardSurface',asShader=True,name='interfaceScatter')
        for attr,value in dict(base=1,subsurface=1,specular=1,specularIOR=1.33,
                               specularRoughness=.03,transmission=0,metalness=0,
                               subsurfaceScale=.2,subsurfaceAnisotropy=.2).items():
            c.setAttr(fog+'.'+attr,value)
        c.setAttr(fog+'.baseColor',1,1,1,type='double3')
        c.setAttr(fog+'.subsurfaceColor',.98,.98,.98,type='double3')
        c.setAttr(fog+'.subsurfaceRadius',1,1,1,type='double3')
        types=c.attributeQuery('subsurfaceType',node=fog,listEnum=True)[0].split(':')
        randomwalk=next(i for i,n in enumerate(types) if 'randomwalk' in n.lower().replace(' ',''))
        c.setAttr(fog+'.subsurfaceType',randomwalk)
        mix=c.shadingNode('aiMixShader',asShader=True,name='interfaceMix')
        c.connectAttr('water.outColor',mix+'.shader1');c.connectAttr(fog+'.outColor',mix+'.shader2')
        gain=c.shadingNode('multDoubleLinear',asUtility=True,name='interfaceEnable')
        c.connectAttr(read+'.outColorR',gain+'.input1');c.connectAttr(gain+'.output',mix+'.mix')
        c.connectAttr(mix+'.outColor','waterSG.surfaceShader',force=True)
        # Verify persisted vertex/face-vertex colors before any licensed render.
        actual=np.asarray(mesh.getVertexColors(colorSet='interfaceMask'))[:,0]
        assert np.max(np.abs(actual-weights))<1e-6
        nodes={'water':['base','metalness','specular','specularIOR','specularRoughness','transmission',
                       'transmissionColor','transmissionDepth','transmissionScatter','thinWalled'],
               'sky':['color','intensity','camera'],'sun':['visibility','rotate'],'sunShape':['color','intensity'],
               'persp':['translate','rotate'],'perspShape':['orthographic','orthographicWidth'],
               'interfaceEnable':['input2'],
               'interfaceScatter':['base','subsurface','subsurfaceType','subsurfaceColor','subsurfaceRadius',
                                   'subsurfaceScale','subsurfaceAnisotropy','specular','specularIOR','specularRoughness'],
               'defaultArnoldDriver':['colorManagement'],
               'defaultArnoldRenderOptions':['AASamples','GISpecularDepth','GITransmissionDepth','GITotalDepth','GIVolumeDepth','GIVolumeSamples'],
               'defaultColorMgtGlobals':['wsn','vtn']}
        for arm in ['A','B']:
            folder=out/arm;folder.mkdir()
            c.setAttr(gain+'.input2',int(arm=='B'))
            c.setAttr('opticalEXR.prefix',(folder/'optics').as_posix(),type='string')
            c.setAttr('defaultArnoldRenderOptions.log_filename',(folder/'arnold.log').as_posix(),type='string')
            state={node:{attr:c.getAttr(node+'.'+attr) for attr in attrs} for node,attrs in nodes.items()}
            state.update(geometry=base.geometry(c),mask_sha256=base.sha(MASK/'mask.npz'),mist=False,
                         approximation='curvature-gated randomwalk SSS, radius scale 5 cm; no separate VDB band')
            assert state['geometry']==before
            base.save(folder/'state.json',state)
            c.file(rename=str(folder/'candidate.mb'));c.file(save=True,type='mayaBinary',force=True)
        base.save(out/'prepare.json',dict(status='passed',wall_s=time.perf_counter()-started))
        base.exit_worker(0)
    except BaseException:
        base.save(out/'prepare.json',dict(status='failed',error=traceback.format_exc()));base.exit_worker(1)


def render(folder):
    import maya.cmds as c
    started=time.perf_counter();result=dict(status='running')
    try:
        c.scriptEditorInfo(writeHistory=True,historyFilename=str(folder/'maya.log'))
        c.file(str(folder/'candidate.mb'),open=True,force=True)
        c.workspace(fileRule=['images',folder.as_posix()])
        assert not c.ls(type='bifrostGraphShape')
        assert c.getAttr('baked_control.visibility')
        assert not c.getAttr('baked_raw.visibility') and not c.getAttr('baked_no_mcf.visibility')
        assert c.getAttr('interfaceEnable.input2')==int(folder.name=='B')
        state=json.loads((folder/'state.json').read_text());assert base.geometry(c)==state['geometry']
        start=time.perf_counter();c.arnoldRender(seq='1',w=960,h=540,cam='persp',srv=False)
        result['render_call_s']=time.perf_counter()-start
        pngs=list((folder/'tmp').rglob('*.png'));assert len(pngs)==1,pngs
        shutil.move(str(pngs[0]),folder/'beauty.png')
        exrs=list(folder.glob('*.exr'));assert len(exrs)==1,exrs
        result.update(status='passed',exr=exrs[0].name,worker_s=time.perf_counter()-started)
        base.save(folder/'worker_result.json',result);base.exit_worker(0)
    except BaseException:
        result.update(status='failed',error=traceback.format_exc())
        base.save(folder/'worker_result.json',result);base.exit_worker(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('mode',choices=['prepare','A','B'])
    parser.add_argument('--out',type=Path,default=OUT);args=parser.parse_args();out=args.out.resolve()
    assert not base.git('status','--porcelain'),'Commit instrument first'
    assert base.sha(base.SOURCE)==base.SOURCE_SHA
    assert json.loads((MASK/'mask_review.json').read_text())['status']=='passed'
    hashes={p.name:base.sha(p) for p in base.CACHE.iterdir() if p.is_file()}
    assert hashes==json.loads((base.PRIOR/'A/receipt.json').read_text())['cache_sha256']
    mask_hash=base.sha(MASK/'mask.npz')
    env={k:v for k,v in os.environ.items() if not k.startswith('SURF_')}
    env.update(INTERFACE_RENDER_MODE=args.mode,INTERFACE_RENDER_OUT=str(out),
               INTERFACE_RENDER_SCRIPT=str(Path(__file__).resolve()),MAYA_DISABLE_CER='1',MAYA_DISABLE_CIP='1',MAYA_DISABLE_CLIC_IPM='1')
    folder=out if args.mode=='prepare' else out/args.mode
    path=folder/('prepare_process.json' if args.mode=='prepare' else 'receipt.json');assert not path.exists()
    if args.mode=='prepare':
        assert not out.exists();out.mkdir(parents=True)
        command=[MAYAPY,str(Path(__file__).resolve())];timeout=180
    else:
        assert json.loads((out/'prepare.json').read_text())['status']=='passed'
        if args.mode=='B':
            a=json.loads((out/'A/receipt.json').read_text());assert a['status']=='passed'
            assert a['git_sha']==base.git('rev-parse','HEAD')
            assert json.loads((out/'A/control_check.json').read_text())['status']=='passed'
        launcher=r"""
$ErrorActionPreference = 'Stop'
if (Get-Process maya,mayapy -ErrorAction SilentlyContinue) { throw 'Another Maya process is open' }
if ((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB -lt 20) { throw 'Less than 20 GB RAM free' }
$pythonScript = $env:INTERFACE_RENDER_SCRIPT -replace '\\', '/'
$command = '"python(\"import runpy; runpy.run_path(r''' + $pythonScript + ''', run_name=''__main__'')\")"'
$p = Start-Process -FilePath 'C:/Program Files/Autodesk/Maya2027/bin/maya.exe' -ArgumentList @('-command',$command) -PassThru -WindowStyle Hidden
Write-Output "MAYA_PID $($p.Id)"
if (-not $p.WaitForExit(600000)) { Stop-Process -Id $p.Id -Force; throw 'Owned render timed out after 600 seconds' }
exit $p.ExitCode
"""
        (folder/'launch.ps1').write_text(launcher,encoding='utf-8')
        command=['powershell.exe','-NoProfile','-ExecutionPolicy','Bypass','-File',str(folder/'launch.ps1')];timeout=630
    receipt=dict(status='running',mode=args.mode,git_sha=base.git('rev-parse','HEAD'),cache_sha256=hashes,
                 source_sha256=base.SOURCE_SHA,mask_sha256=mask_hash,started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
    base.save(path,receipt);start=time.perf_counter()
    try:
        with (folder/(args.mode+'_launch.log')).open('w') as log:
            proc=subprocess.run(command,env=env,cwd=base.ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
        receipt['exit_code']=proc.returncode;assert proc.returncode==0,'Read worker log'
        assert hashes=={p.name:base.sha(p) for p in base.CACHE.iterdir() if p.is_file()}
        assert base.sha(base.SOURCE)==base.SOURCE_SHA and base.sha(MASK/'mask.npz')==mask_hash
        if args.mode!='prepare':
            worker=json.loads((folder/'worker_result.json').read_text());assert worker['status']=='passed'
            receipt.update(worker=worker,png_sha256=base.sha(folder/'beauty.png'),exr_sha256=base.sha(folder/worker['exr']))
        receipt['status']='passed'
    except BaseException:
        receipt.update(status='failed',error=traceback.format_exc());raise
    finally:
        receipt['process_s']=time.perf_counter()-start;base.save(path,receipt);print(json.dumps(receipt,indent=2))


if __name__=='__main__':
    if 'INTERFACE_RENDER_MODE' in os.environ:
        out=Path(os.environ['INTERFACE_RENDER_OUT']);mode=os.environ['INTERFACE_RENDER_MODE']
        if mode=='prepare':prepare(out)
        else:render(out/mode)
    else:main()
