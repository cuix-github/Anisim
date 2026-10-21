"""Fixed-camera ocean look study. Existing geometry only; no simulation/foam generation.

The caller supplies a frozen scene and a separate new output directory. Preparation
validates camera/mesh/units, copies external volumes, and saves the effective setup.
Rendering is a separate explicit operation. Never overwrite a completed arm.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import traceback

import numpy as np

HERE = Path(__file__).resolve()
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parents[1] / 'breaker'))
import breaker_deep_water as base
from analyze_breaker_deep_water import OIIO, MAYAPY, COMPONENTS, LAYERS, pfm, convert
from breaker_foam_layer import ass_units

WATER = dict(base=0, metalness=0, specular=1, specularIOR=1.333,
             specularRoughness=.08, transmission=1, transmissionDepth=40,
             transmissionColor=(.40, .54, .55),
             transmissionScatter=(.018, .025, .024),
             transmissionScatterAnisotropy=.5, transmissionExtraRoughness=0,
             thinWalled=False)
CAMERA = dict(persp=['translate', 'rotate', 'scale'],
              perspShape=['orthographic', 'focalLength', 'horizontalFilmAperture',
                          'verticalFilmAperture', 'filmFit', 'lensSqueezeRatio',
                          'horizontalFilmOffset', 'verticalFilmOffset', 'depthOfField'])


def read(p):
    return json.loads(Path(p).read_text())


def capture(c, nodes):
    return json.loads(json.dumps({n: {a: c.getAttr(n+'.'+a) for a in attrs}
                                 for n, attrs in nodes.items()}))


def geometry_equal(a, b):
    key = 'points_normals_topology_sha256'
    return ({k: v for k, v in a.items() if k != key} ==
            {k: v for k, v in b.items() if k != key} and
            [a[key][i] for i in (0, 2, 3)] == [b[key][i] for i in (0, 2, 3)])


def environment(out):
    """Linear ACEScg radiance map; broad overcast structure, no water texture."""
    w, h = 2048, 1024
    theta = (np.arange(h)[:, None] + .5) / h * np.pi
    phi = (np.arange(w)[None, :] + .5) / w * 2*np.pi
    y = np.cos(theta) + np.zeros_like(phi)
    x, z = np.sin(theta)*np.sin(phi), np.sin(theta)*np.cos(phi)
    up = np.clip(y, 0, 1)
    axis = np.array([-.45, .65, -.61237244]); axis /= np.linalg.norm(axis)
    lobe = np.exp((axis[0]*x + axis[1]*y + axis[2]*z - 1) / .23)
    structure = .10*np.sin(3*phi+1.2)*np.sin(theta)**2 + .06*np.cos(5*phi-2*theta)
    top = .70 + .60*np.sqrt(up) + 2.2*lobe + structure
    lower = .06 + .20*np.exp(y/.12)
    blend = np.clip((y+.04)/.10, 0, 1); blend=blend*blend*(3-2*blend)
    radiance = (blend*top+(1-blend)*lower)[..., None]*np.array([.97, .99, 1.0])
    assert np.isfinite(radiance).all() and radiance.min() > 0
    p = out/'overcast_linear.pfm'
    with p.open('wb') as f:
        f.write(f'PF\n{w} {h}\n-1.0\n'.encode())
        f.write(np.flipud(radiance).astype('<f4').tobytes())
    target=out/'overcast.exr'
    subprocess.run([OIIO, str(p), '--attrib', 'oiio:ColorSpace', 'ACEScg',
                    '-d', 'half', '-o', str(target)], check=True)
    return dict(type='analytic neutral overcast', resolution=[w,h],
                color_space='linear ACEScg, read without conversion',
                radiance_min=float(radiance.min()), radiance_max=float(radiance.max()),
                radiance_mean_rgb=radiance.mean((0,1)).tolist(), sha256=base.sha(target))


def scale_mist(source, target, multiplier):
    text=source.read_text()
    pattern=r'(?m)^(\s*density\s+)([-+0-9.eE]+)(\s*)$'
    result, count=re.subn(pattern, lambda m:m[1]+format(float(m[2])*multiplier,'.12g')+m[3], text)
    assert count>0
    assert re.sub(pattern,'DENSITY',text)==re.sub(pattern,'DENSITY',result)
    target.write_text(result)
    return dict(multiplier=multiplier, material_nodes=count,
                positions_radii_and_other_parameters_unchanged=True, sha256=base.sha(target))


def prepare_worker(source, out, config):
    import maya.standalone
    maya.standalone.initialize(name='python')
    import maya.cmds as c
    try:
        c.file(str(source/'candidate.mb'), open=True, force=True)
        assert not c.ls(type='bifrostGraphShape'), 'No live fluid graph allowed'
        before=base.geometry(c); camera=capture(c,CAMERA)
        original=read(source/'state.json')
        assert geometry_equal(before,original['geometry'])
        for key in ['water','sky','sun','sunShape','defaultArnoldDriver',
                    'defaultArnoldRenderOptions','defaultColorMgtGlobals']:
            assert capture(c,{key:list(original[key])})[key]==original[key],key
        # Every look knob is explicit; the camera and surface remain untouched.
        water=dict(WATER); water.update(config.get('water',{}))
        for attr,v in water.items():
            if isinstance(v,(list,tuple)): c.setAttr('water.'+attr,*v,type='double3')
            else: c.setAttr('water.'+attr,v)
        file=c.shadingNode('file',asTexture=True,name='oceanOvercastRadiance')
        c.setAttr(file+'.fileTextureName',(out/'overcast.exr').as_posix(),type='string')
        c.setAttr(file+'.colorSpace','Raw',type='string')
        c.setAttr(file+'.ignoreColorSpaceFileRules',True)
        c.connectAttr(file+'.outColor','sky.color',force=True)
        c.setAttr('sky.intensity',config.get('sky_intensity',1.0))
        c.setAttr('sky.camera',0)
        c.setAttr('sunShape.intensity',config.get('sun_intensity',1.0))
        c.setAttr('sunShape.color',1,.985,.96,type='double3')
        if c.attributeQuery('aiAngle',node='sunShape',exists=True): c.setAttr('sunShape.aiAngle',3)
        c.setAttr('foamLayerShape.dso',(out/'foam.ass').as_posix(),type='string')
        c.setAttr('mistLayerShape.dso',(out/'mist.ass').as_posix(),type='string')
        for attr,v in [('width',1920),('height',1080),('deviceAspectRatio',16/9)]:
            c.setAttr('defaultResolution.'+attr,v)
        opts={'AASamples':4,'GISpecularDepth':3,'GITransmissionDepth':8,
              'GITotalDepth':10,'GIVolumeDepth':3,'GIVolumeSamples':2}
        for attr,v in opts.items(): c.setAttr('defaultArnoldRenderOptions.'+attr,v)
        for attr,v in [('threads_autodetect',False),('threads',16)]:
            if c.attributeQuery(attr,node='defaultArnoldRenderOptions',exists=True):
                c.setAttr('defaultArnoldRenderOptions.'+attr,v)
        c.setAttr('opticalEXR.prefix',(out/'optics').as_posix(),type='string')
        c.setAttr('defaultArnoldRenderOptions.log_filename',(out/'arnold.log').as_posix(),type='string')
        c.setAttr('defaultRenderGlobals.imageFilePrefix','beauty',type='string')
        state=capture(c,{**CAMERA,'water':list(water),'sky':['intensity','camera'],
                         'sun':['rotate','visibility'],'sunShape':['intensity','color'],
                         'defaultArnoldRenderOptions':list(opts),
                         'defaultColorMgtGlobals':['wsn','vtn'],
                         'defaultArnoldDriver':['colorManagement']})
        state.update(geometry=base.geometry(c),view=original['view'],
                     camera_before=camera,environment=read(out/'environment.json'),
                     foam_ass_sha256=base.sha(out/'foam.ass'),mist=read(out/'mist_adjustment.json'),
                     source_scene_sha256=base.sha(source/'candidate.mb'),
                     simulation=False,new_whitewater=False,new_particles=False,
                     new_surface_or_displacement=False)
        assert capture(c,CAMERA)==camera,'Camera changed'
        assert geometry_equal(before,state['geometry']),'Geometry changed'
        assert state['defaultColorMgtGlobals']==original['defaultColorMgtGlobals']
        base.save(out/'state.json',state)
        c.file(rename=str(out/'candidate.mb'));c.file(save=True,type='mayaBinary',force=True)
        c.arnoldExportAss(f=str(out/'scene_check.ass'),mask=24)
        assert ass_units(out/'scene_check.ass')==ass_units(out/'foam.ass')==ass_units(out/'mist.ass')==.01
        base.save(out/'prepare.json',dict(status='passed',camera_unchanged=True,geometry_unchanged=True,
                    foam_file_unchanged=base.sha(out/'foam.ass')==base.sha(source/'foam.ass'),
                    scene_sha256=base.sha(out/'candidate.mb')))
        base.exit_worker(0)
    except BaseException:
        base.save(out/'prepare.json',dict(status='failed',error=traceback.format_exc()));base.exit_worker(1)


def render_worker(out):
    import maya.cmds as c
    start=time.perf_counter()
    try:
        c.scriptEditorInfo(writeHistory=True,historyFilename=str(out/'maya.log'))
        c.file(str(out/'candidate.mb'),open=True,force=True)
        s=read(out/'state.json')
        assert capture(c,CAMERA)==s['camera_before']
        assert geometry_equal(base.geometry(c),s['geometry'])
        assert not c.ls(type='bifrostGraphShape')
        assert c.listConnections('waterSG.surfaceShader')==['water']
        c.workspace(fileRule=['images',out.as_posix()])
        t=time.perf_counter();c.arnoldRender(seq='1',w=1920,h=1080,cam='persp',srv=False)
        render_s=time.perf_counter()-t
        pngs=list((out/'tmp').rglob('*.png')); assert len(pngs)==1
        shutil.move(str(pngs[0]),out/'beauty.png')
        exrs=list(out.glob('optics*.exr'));assert len(exrs)==1
        base.save(out/'worker_result.json',dict(status='passed',render_call_s=render_s,
                  worker_s=time.perf_counter()-start,exr=exrs[0].name))
        base.exit_worker(0)
    except BaseException:
        base.save(out/'worker_result.json',dict(status='failed',error=traceback.format_exc()));base.exit_worker(1)


def main():
    ap=argparse.ArgumentParser();ap.add_argument('stage',choices=['prepare','render','check'])
    ap.add_argument('--source',type=Path,required=True);ap.add_argument('--out',type=Path,required=True)
    ap.add_argument('--config',type=Path);args=ap.parse_args()
    source,out=args.source.resolve(),args.out.resolve()
    if args.stage=='prepare':
        assert not out.exists(),'Use a new output directory'
        assert not base.git('status','--porcelain'),'Commit instrument first'
        out.mkdir(parents=True); config=read(args.config) if args.config else {}
        base.save(out/'requested_config.json',config)
        base.save(out/'environment.json',environment(out))
        shutil.copy2(source/'foam.ass',out/'foam.ass')
        base.save(out/'mist_adjustment.json',scale_mist(source/'mist.ass',out/'mist.ass',config.get('mist_density_multiplier',.15)))
        base.save(out/'provenance.json',dict(git_sha=base.git('rev-parse','HEAD'),instrument_sha256=base.sha(HERE),
                       source=str(source),source_scene_sha256=base.sha(source/'candidate.mb'),
                       operation='integrated look candidate, not single-variable causal experiment'))
        env=dict(os.environ,OCEAN_LOOK_WORKER='prepare',OCEAN_LOOK_SOURCE=str(source),OCEAN_LOOK_OUT=str(out))
        with (out/'prepare.log').open('w') as log:
            subprocess.run([MAYAPY,str(HERE)],env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=300)
        assert read(out/'prepare.json')['status']=='passed'
        print(json.dumps(read(out/'prepare.json')))
    elif args.stage=='render':
        assert not (out/'receipt.json').exists(),'No repeated render in an existing arm'
        assert read(out/'prepare.json')['scene_sha256']==base.sha(out/'candidate.mb')
        bootstrap=out/'bootstrap.mel'
        python_code='import runpy; runpy.run_path('+repr(HERE.as_posix())+", run_name='__main__')"
        bootstrap.write_text('python('+json.dumps(python_code)+');\n')
        launcher='''$ErrorActionPreference = 'Stop'
if (Get-Process maya,mayapy,bifcmd -ErrorAction SilentlyContinue) { throw 'Another Maya or simulation process is active' }
$p = Start-Process -FilePath 'C:/Program Files/Autodesk/Maya2027/bin/maya.exe' -ArgumentList @('-script',('\"' + $env:OCEAN_LOOK_BOOTSTRAP + '\"')) -PassThru -WindowStyle Hidden
Write-Output "MAYA_PID $($p.Id)"
if (-not $p.WaitForExit(1800000)) { Stop-Process -Id $p.Id -Force; throw 'Owned render exceeded 30 minutes' }
exit $p.ExitCode
'''
        (out/'launch.ps1').write_text(launcher)
        env=dict(os.environ,OCEAN_LOOK_WORKER='render',OCEAN_LOOK_OUT=str(out),OCEAN_LOOK_BOOTSTRAP=str(bootstrap),
                 MAYA_DISABLE_CER='1',MAYA_DISABLE_CIP='1',MAYA_DISABLE_CLIC_IPM='1')
        receipt=dict(status='running',git_sha=base.git('rev-parse','HEAD'),
                     started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()))
        base.save(out/'receipt.json',receipt);start=time.perf_counter()
        try:
            with (out/'launch.log').open('w') as log:
                subprocess.run(['powershell.exe','-NoProfile','-ExecutionPolicy','Bypass','-File',str(out/'launch.ps1')],
                               env=env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=1830)
            worker=read(out/'worker_result.json');assert worker['status']=='passed',worker
            receipt.update(status='passed',worker=worker,png_sha256=base.sha(out/'beauty.png'),
                           exr_sha256=base.sha(out/worker['exr']))
        except BaseException:
            receipt.update(status='failed',error=traceback.format_exc());raise
        finally:
            receipt['process_s']=time.perf_counter()-start;base.save(out/'receipt.json',receipt)
    else:
        receipt=read(out/'receipt.json');assert receipt['status']=='passed'
        exr=out/receipt['worker']['exr'];assert base.sha(exr)==receipt['exr_sha256']
        info=subprocess.check_output([OIIO,'--info','-v',str(exr)],text=True)
        (out/'exr_info.txt').write_text(info)
        assert '54 channel' in info and 'arnold/license_state: "valid"' in info
        data={}
        for layer in LAYERS:
            ch='R,G,B' if layer=='beauty' else ','.join(layer+'.'+s for s in ('XYZ' if layer in ['N','P'] else 'RGB'))
            target=out/(layer+'.pfm')
            subprocess.run([OIIO,str(exr),'--ch',ch,'-o',str(target)],check=True)
            data[layer]=pfm(target)
        env=dict(os.environ,OCEAN_LOOK_WORKER='convert',OCEAN_LOOK_OUT=str(out))
        subprocess.run([MAYAPY,str(HERE)],env=env,check=True)
        total=data['beauty_raw'].sum(dtype=float)
        closure=float(np.abs(sum(data[k].astype(float) for k in COMPONENTS)-data['beauty_raw']).sum()/total)
        assert closure<1e-5
        result=dict(status='passed',aov_closure_relative_l1=closure,
                    component_raw_rgb_share={k:float(data[k].sum(dtype=float)/total) for k in COMPONENTS},
                    render_call_s=receipt['worker']['render_call_s'])
        base.save(out/'image_check.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    worker=os.environ.get('OCEAN_LOOK_WORKER')
    if worker=='prepare':
        folder=Path(os.environ['OCEAN_LOOK_OUT']);prepare_worker(Path(os.environ['OCEAN_LOOK_SOURCE']),folder,read(folder/'requested_config.json'))
    elif worker=='render': render_worker(Path(os.environ['OCEAN_LOOK_OUT']))
    elif worker=='convert': convert(Path(os.environ['OCEAN_LOOK_OUT']))
    else: main()
