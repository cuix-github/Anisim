"""Two approved cache-only previews: neutral water optics and finer existing mist.

Reuses the frozen mesh and existing particle_mist compound; never runs the solver.
The side view changes only the camera and frames a 15 m portion of the curl.
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

sys.path.insert(0, str(Path(__file__).resolve().parent))
import breaker_deep_water as base

ROOT = base.ROOT
OUT = ROOT / 'out/runs/breaker_neutral_mist_20261019'
MAYAPY = 'C:/Program Files/Autodesk/Maya2027/bin/mayapy.exe'
MIST = dict(surface_radius=.5, surface_detail=.5, smooth_iterations=0,
            mist_threshold=.75, spray_detail=.25, cloud_radius=.5,
            metres_per_cell=.25, default_radius=.0001, mist_density=.00025,
            spray_epsilon=10000, mcf_iterations=10, mask_sigma=2,
            mask_theta=2, mask_zeta=5, mask_floor=0, mask_dilate=1)
WATER = dict(base=0, metalness=0, specular=1, specularIOR=1.333,
             specularRoughness=.12, transmission=1, transmissionDepth=40,
             transmissionColor=(.3,.55,.6), transmissionScatter=(.15,.15,.15),
             transmissionScatterAnisotropy=.75, transmissionExtraRoughness=0,
             thinWalled=False)
SKY_INTENSITY = .44 / .171667  # Match baseline dome mean RGB; not total illumination.


def prepare(out):
    # These scene helpers set BIFROST_LIB_CONFIG_FILES at import time. Import
    # before Maya prepends its shipped libraries during standalone.initialize.
    sys.path.insert(0,str(ROOT/'maya'))
    from build_surface_scene import graph
    import maya.standalone
    maya.standalone.initialize(name='python')
    import maya.cmds as c
    started=time.perf_counter()
    try:
        for plugin in ('bifrostGraph','mtoa'):
            if not c.pluginInfo(plugin,q=True,loaded=True):
                c.loadPlugin(plugin)
        c.file(str(base.SOURCE),open=True,force=True)
        before=base.geometry(c)
        assert not c.ls(type='bifrostGraphShape')
        for attr,value in WATER.items():
            if isinstance(value,tuple): c.setAttr('water.'+attr,*value,type='double3')
            else: c.setAttr('water.'+attr,value)
        physical=c.shadingNode('aiPhysicalSky',asTexture=True,name='candidateSky')
        for attr,value in base.SKY.items(): c.setAttr(physical+'.'+attr,value)
        c.connectAttr(physical+'.outColor','sky.color',force=True)
        c.setAttr('sky.intensity',SKY_INTENSITY)
        c.setAttr('sky.aiExposure',0)
        c.setAttr('sky.resolution',2048)
        c.setAttr('sun.visibility',False)
        ports=[('liquid_pattern',(base.CACHE/'liquid.####').as_posix()),
               ('radius_pattern',(base.CACHE/'radius.####').as_posix()),*MIST.items()]
        shape=graph(c,'candidateMistShape','particle_mist','mist',ports)
        fog=c.shadingNode('aiStandardVolume',asShader=True,name='candidateMistFog')
        c.setAttr(fog+'.densityChannel','voxel_fog_density',type='string')
        c.setAttr(fog+'.density',10000)
        c.setAttr(fog+'.scatter',1)
        c.setAttr(fog+'.scatterColor',1,1,1,type='double3')
        c.setAttr(fog+'.scatterAnisotropy',.4)
        sg=c.sets(renderable=True,noSurfaceShader=True,empty=True,name='candidateMistSG')
        c.connectAttr(fog+'.outColor',sg+'.volumeShader')
        c.sets(shape,e=True,forceElement=sg)
        assert c.listConnections('waterSG.surfaceShader')==['water']
        c.currentTime(1)
        modes=c.attributeQuery('colorManagement',node='defaultArnoldDriver',listEnum=True)[0].split(':')
        c.setAttr('defaultArnoldDriver.colorManagement',modes.index('Use View Transform'))
        nodes={'water':list(WATER),'candidateSky':list(base.SKY),
               'sky':['intensity','aiExposure','resolution','camera'],'sun':['visibility'],
               'persp':['translate','rotate'],'perspShape':['focalLength','orthographic','orthographicWidth'],
               'candidateMistFog':['density','densityChannel','scatter','scatterColor','scatterAnisotropy'],
               'defaultArnoldDriver':['colorManagement'],
               'defaultArnoldRenderOptions':['AASamples','GISpecularDepth','GITransmissionDepth','GITotalDepth','GIVolumeDepth','GIVolumeSamples'],
               'defaultColorMgtGlobals':['wsn','vtn']}
        for arm in ('original','side'):
            folder=out/arm
            if arm=='side':
                # Horizontal sightline normal to the dominant +12 degree train.
                # Aim at the curling lip, rather than the entire cropped slab.
                import math
                theta=math.radians(12)
                center=(145.5*4,13.5*4,70*4)
                forward=(-math.sin(theta),0,math.cos(theta))
                c.setAttr('persp.translate',*[center[i]-300*forward[i] for i in range(3)],type='double3')
                c.setAttr('persp.rotate',0,168,0,type='double3')
                c.setAttr('perspShape.orthographic',True)
                c.setAttr('perspShape.orthographicWidth',60)
            c.workspace(fileRule=['images',folder.as_posix()])
            c.setAttr('defaultArnoldRenderOptions.log_filename',(folder/'arnold.log').as_posix(),type='string')
            c.setAttr('opticalEXR.prefix',(folder/'optics').as_posix(),type='string')
            state={node:{attr:c.getAttr(node+'.'+attr) for attr in attrs} for node,attrs in nodes.items()}
            state.update(geometry=base.geometry(c),mist_ports=MIST,
                         mist_compound='User::PFFlip::particle_mist',mist_visible=c.getAttr(shape+'.visibility'))
            assert state['geometry']==before
            base.save(folder/'state.json',state)
            c.file(rename=str(folder/'candidate.mb'))
            c.file(save=True,type='mayaBinary',force=True)
        base.save(out/'prepare.json',dict(status='passed',wall_s=time.perf_counter()-started,
                                        geometry=before,mist_ports=MIST))
        base.exit_worker(0)
    except BaseException:
        base.save(out/'prepare.json',dict(status='failed',error=traceback.format_exc()))
        base.exit_worker(1)


def render(out):
    import maya.cmds as c
    started=time.perf_counter()
    result=dict(status='running')
    try:
        c.scriptEditorInfo(writeHistory=True,historyFilename=str(out/'maya.log'))
        c.file(str(out/'candidate.mb'),open=True,force=True)
        assert c.ls(type='bifrostGraphShape')==['candidateMistShape']
        assert c.getAttr('baked_control.visibility') and c.getAttr('candidateMistShape.visibility')
        assert not c.getAttr('baked_raw.visibility') and not c.getAttr('baked_no_mcf.visibility')
        state=json.loads((out/'state.json').read_text())
        assert base.geometry(c)==state['geometry']
        result['open_s']=time.perf_counter()-started
        tick=time.perf_counter()
        c.arnoldRender(seq='1',w=960,h=540,cam='persp',srv=False)
        result['render_call_s']=time.perf_counter()-tick
        pngs=list((out/'tmp').rglob('*.png')); assert len(pngs)==1,pngs
        shutil.move(str(pngs[0]),out/'beauty.png')
        exrs=list(out.rglob('*.exr')); assert len(exrs)==1,exrs
        result.update(status='passed',exr=exrs[0].name,worker_s=time.perf_counter()-started)
        base.save(out/'worker_result.json',result)
        base.exit_worker(0)
    except BaseException:
        result.update(status='failed',error=traceback.format_exc())
        base.save(out/'worker_result.json',result)
        base.exit_worker(1)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('mode',choices=['prepare','original','side'])
    parser.add_argument('--out',type=Path,default=OUT)
    args=parser.parse_args(); out=args.out.resolve()
    assert not base.git('status','--porcelain'),'Commit before measuring'
    assert base.sha(base.SOURCE)==base.SOURCE_SHA
    baseline=json.loads((base.PRIOR/'A/receipt.json').read_text())
    assert baseline['control_check']['status']=='passed'
    hashes={p.name:base.sha(p) for p in base.CACHE.iterdir() if p.is_file()}
    assert hashes==baseline['cache_sha256']
    script=Path(__file__).resolve()
    env={k:v for k,v in os.environ.items() if not k.startswith('SURF_')}
    env.update(BIFROST_LIB_CONFIG_FILES=str(ROOT/'maya/pfflip_graphs_config.json'),
               MAYA_DISABLE_CER='1',MAYA_DISABLE_CIP='1',MAYA_DISABLE_CLIC_IPM='1',
               NEUTRAL_MODE=args.mode,NEUTRAL_OUT=str(out),NEUTRAL_SCRIPT=str(script))
    if args.mode=='prepare':
        assert not out.exists(),out
        out.mkdir(parents=True)
        for arm in ('original','side'): (out/arm).mkdir()
        command=[MAYAPY,str(script)]
        folder=out; timeout=600
    else:
        assert json.loads((out/'prepare.json').read_text())['status']=='passed'
        folder=out/args.mode; assert not (folder/'receipt.json').exists()
        if args.mode=='side':
            prior=json.loads((out/'original/receipt.json').read_text())
            assert prior['status']=='passed' and prior['git_sha']==base.git('rev-parse','HEAD')
        launcher=r"""
$ErrorActionPreference = 'Stop'
if (Get-Process maya,mayapy -ErrorAction SilentlyContinue) { throw 'Another Maya process is open' }
if ((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB -lt 20) { throw 'Less than 20 GB RAM free' }
$pythonScript = $env:NEUTRAL_SCRIPT -replace '\\', '/'
$command = '"python(\"import runpy; runpy.run_path(r''' + $pythonScript + ''', run_name=''__main__'')\")"'
$p = Start-Process -FilePath 'C:/Program Files/Autodesk/Maya2027/bin/maya.exe' -ArgumentList @('-command',$command) -PassThru -WindowStyle Hidden
Write-Output "MAYA_PID $($p.Id)"
if (-not $p.WaitForExit(480000)) { Stop-Process -Id $p.Id -Force; throw 'Owned render timed out after 480 seconds' }
exit $p.ExitCode
"""
        (folder/'launch.ps1').write_text(launcher,encoding='utf-8')
        command=['powershell.exe','-NoProfile','-ExecutionPolicy','Bypass','-File',str(folder/'launch.ps1')]
        timeout=510
    receipt=dict(mode=args.mode,status='running',git_sha=base.git('rev-parse','HEAD'),
                 started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                 cache_sha256=hashes,source_sha256=base.SOURCE_SHA)
    base.save(folder/'receipt.json',receipt)
    start=time.perf_counter()
    try:
        with (folder/'launch.log').open('w',encoding='utf-8') as log:
            proc=subprocess.run(command,env=env,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,timeout=timeout)
        receipt['exit_code']=proc.returncode
        assert proc.returncode==0,'Read launch.log / worker_result.json'
        assert hashes=={p.name:base.sha(p) for p in base.CACHE.iterdir() if p.is_file()}
        if args.mode!='prepare':
            worker=json.loads((folder/'worker_result.json').read_text()); assert worker['status']=='passed'
            receipt.update(worker=worker,png_sha256=base.sha(folder/'beauty.png'),
                           exr_sha256=base.sha(folder/worker['exr']),scene_sha256=base.sha(folder/'candidate.mb'))
        receipt['status']='passed'
    except BaseException:
        receipt.update(status='failed',error=traceback.format_exc()); raise
    finally:
        receipt['process_s']=time.perf_counter()-start
        base.save(folder/'receipt.json',receipt)
        print(json.dumps(receipt,indent=2),flush=True)


if __name__=='__main__':
    if 'NEUTRAL_MODE' in os.environ:
        out=Path(os.environ['NEUTRAL_OUT']); mode=os.environ['NEUTRAL_MODE']
        if mode=='prepare': prepare(out)
        else: render(out/mode)
    else: main()
