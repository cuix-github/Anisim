"""Two authorized stills on a validated frozen mesh; no graph or solver runs.

This is an integrated look candidate, not a one-variable causal experiment.
The original-camera image changes the material/lighting package. The side image
changes only its camera. Existing PR #39 artifacts supply the validated control.
"""
import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
PRIOR = ROOT / 'out/runs/breaker_optical_aovs_20261019_v3'
SOURCE = PRIOR / 'B/optical.mb'
SOURCE_SHA = 'd1ddd4682c46f04e76aae5586d2cf74215617e0d4012318393e26404e81e6100'
CACHE = ROOT / 'out/runs/wavesea4/cache/breaking_section'
DEFAULT = ROOT / 'out/runs/breaker_deep_water_20261019'
WATER = dict(base=0, metalness=0, specular=1, specularIOR=1.333,
             specularRoughness=.12, transmission=1, transmissionDepth=10,
             transmissionColor=(.744, .9400832057, 1),
             transmissionScatter=(.0675, .1340401471, .5),
             transmissionScatterAnisotropy=.75, transmissionExtraRoughness=0,
             thinWalled=False)
SKY = dict(elevation=30, azimuth=160, turbidity=6, sunSize=5,
           intensity=1, enableSun=True, enableSky=True)


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def exit_worker(code):
    import ctypes
    sys.stdout.flush()
    sys.stderr.flush()
    terminate = ctypes.windll.kernel32.TerminateProcess
    terminate.argtypes = [ctypes.c_void_p, ctypes.c_uint]
    terminate.restype = ctypes.c_int
    terminate(ctypes.c_void_p(-1), code)
    raise ctypes.WinError()


def geometry(c):
    import maya.api.OpenMaya as om
    import numpy as np
    selection = om.MSelectionList()
    selection.add('baked_control')
    dag = selection.getDagPath(0)
    dag.extendToShape()
    mesh = om.MFnMesh(dag)
    counts, indices = mesh.getVertices()
    arrays = [np.asarray(mesh.getPoints(om.MSpace.kWorld), dtype='<f8'),
              np.asarray(mesh.getNormals(), dtype='<f4'),
              np.asarray(counts, dtype='<i4'), np.asarray(indices, dtype='<i4')]
    return dict(vertices=mesh.numVertices, triangles=len(mesh.getTriangles()[1]) // 3,
                bounds=c.exactWorldBoundingBox('baked_control'),
                points_normals_topology_sha256=[hashlib.sha256(a.tobytes()).hexdigest() for a in arrays])


def worker(out, arm):
    import maya.cmds as c
    started = time.perf_counter()
    result = dict(arm=arm, status='running')
    try:
        c.scriptEditorInfo(writeHistory=True, historyFilename=str(out / 'maya.log'))
        c.file(str(SOURCE), open=True, force=True)
        assert not c.ls(type='bifrostGraphShape'), 'Frozen scene must not evaluate a graph'
        assert c.getAttr('baked_control.visibility')
        assert not c.getAttr('baked_raw.visibility') and not c.getAttr('baked_no_mcf.visibility')
        before = geometry(c)
        previous = json.loads((PRIOR / 'B/state.json').read_text())
        assert before['bounds'] == previous['geometry']['bounds']
        assert before['triangles'] == 292516
        assert c.listConnections('waterSG.surfaceShader') == ['water'], 'Foam must be bypassed'
        assert c.getAttr('defaultArnoldRenderOptions.GISpecularDepth') == 3
        c.currentTime(1)
        for attr, value in WATER.items():
            if isinstance(value, tuple):
                c.setAttr('water.' + attr, *value, type='double3')
            else:
                c.setAttr('water.' + attr, value)
        physical = c.shadingNode('aiPhysicalSky', asTexture=True, name='candidateSky')
        for attr, value in SKY.items():
            c.setAttr(physical + '.' + attr, value)
        c.connectAttr(physical + '.outColor', 'sky.color', force=True)
        c.setAttr('sky.intensity', 1)
        c.setAttr('sky.aiExposure', 0)
        c.setAttr('sky.resolution', 2048)
        c.setAttr('sun.visibility', False)
        # World units are solver cells: 0.25 m each. The dominant wave train
        # travels at +12 degrees in XZ. This horizontal sightline is perpendicular
        # to that direction. Other trains have different headings.
        if arm == 'side':
            x0, y0, z0, x1, y1, z1 = before['bounds']
            theta = math.radians(12)
            direction = (-math.sin(theta), 0, math.cos(theta))
            center = ((x0+x1)/2, (y0+y1)/2, (z0+z1)/2)
            eye = [center[i] - 300 * direction[i] for i in range(3)]
            c.setAttr('persp.translate', *eye, type='double3')
            c.setAttr('persp.rotate', 0, 168, 0, type='double3')
            c.setAttr('perspShape.orthographic', True)
            width = ((x1-x0)*math.cos(theta) + (z1-z0)*math.sin(theta)) * 1.08
            c.setAttr('perspShape.orthographicWidth', width)
        for attr, value in [('width', 960), ('height', 540), ('deviceAspectRatio', 960/540)]:
            c.setAttr('defaultResolution.' + attr, value)
        c.setAttr('defaultArnoldDriver.ai_translator', 'png', type='string')
        modes = c.attributeQuery('colorManagement', node='defaultArnoldDriver', listEnum=True)[0].split(':')
        c.setAttr('defaultArnoldDriver.colorManagement', modes.index('Use View Transform'))
        c.setAttr('opticalEXR.prefix', (out / 'optics').as_posix(), type='string')
        c.setAttr('defaultArnoldRenderOptions.log_filename', (out / 'arnold.log').as_posix(), type='string')
        c.workspace(fileRule=['images', out.as_posix()])
        c.setAttr('defaultRenderGlobals.imageFilePrefix', 'beauty', type='string')
        assert geometry(c) == before, 'Look candidate changed geometry or normals'
        nodes = {'water': list(WATER), 'candidateSky': list(SKY),
                 'sky': ['intensity','aiExposure','resolution','camera'],
                 'sun': ['visibility'],
                 'persp': ['translate','rotate'],
                 'perspShape': ['focalLength','orthographic','orthographicWidth'],
                 'defaultArnoldDriver': ['colorManagement'],
                 'defaultArnoldRenderOptions': ['AASamples','GISpecularDepth','GITransmissionDepth','GITotalDepth','GIVolumeDepth','GIVolumeSamples'],
                 'defaultColorMgtGlobals': ['wsn','vtn']}
        state = {node: {attr: c.getAttr(node+'.'+attr) for attr in attrs} for node,attrs in nodes.items()}
        state['geometry'] = before
        state['sky_connection'] = c.listConnections('sky.color', plugs=True)
        save(out / 'state.json', state)
        c.file(rename=str(out / 'candidate.mb'))
        c.file(save=True, type='mayaBinary', force=True)
        result['setup_s'] = time.perf_counter() - started
        render_start = time.perf_counter()
        c.arnoldRender(seq='1', w=960, h=540, cam='persp', srv=False)
        result['render_call_s'] = time.perf_counter() - render_start
        pngs = list((out / 'tmp').rglob('*.png'))
        assert len(pngs) == 1, pngs
        shutil.move(str(pngs[0]), out / 'beauty.png')
        exrs = list(out.rglob('*.exr'))
        assert len(exrs) == 1, exrs
        result.update(status='passed', exr=exrs[0].name, worker_s=time.perf_counter()-started)
        save(out / 'worker_result.json', result)
        exit_worker(0)
    except BaseException:
        result.update(status='failed', error=traceback.format_exc())
        save(out / 'worker_result.json', result)
        exit_worker(1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('arm', choices=['original', 'side'])
    parser.add_argument('--out', type=Path, default=DEFAULT)
    args = parser.parse_args()
    assert not git('status', '--porcelain'), 'Commit instrument before measuring'
    assert sha(SOURCE) == SOURCE_SHA
    control = json.loads((PRIOR / 'A/receipt.json').read_text())
    assert control['status'] == 'passed' and control['control_check']['status'] == 'passed'
    cache_hashes = {p.name: sha(p) for p in CACHE.iterdir() if p.is_file()}
    assert cache_hashes == control['cache_sha256']
    out = args.out.resolve() / args.arm
    assert not out.exists(), out
    if args.arm == 'side':
        first = json.loads((out.parent / 'original/receipt.json').read_text())
        assert first['status'] == 'passed' and first['git_sha'] == git('rev-parse','HEAD')
    out.mkdir(parents=True)
    receipt = dict(arm=args.arm, status='running', git_sha=git('rev-parse','HEAD'),
                   started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                   source_sha256=SOURCE_SHA, cache_sha256=cache_hashes,
                   validated_control=str(PRIOR), render_timeout_s=360)
    save(out / 'receipt.json', receipt)
    launcher = r"""
$ErrorActionPreference = 'Stop'
if (Get-Process maya,mayapy -ErrorAction SilentlyContinue) { throw 'Another Maya process is open' }
if ((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB -lt 20) { throw 'Less than 20 GB RAM free' }
$pythonScript = $env:BREAKER_LOOK_SCRIPT -replace '\\', '/'
$command = '"python(\"import runpy; runpy.run_path(r''' + $pythonScript + ''', run_name=''__main__'')\")"'
$p = Start-Process -FilePath 'C:/Program Files/Autodesk/Maya2027/bin/maya.exe' -ArgumentList @('-command', $command) -PassThru -WindowStyle Hidden
Write-Output "MAYA_PID $($p.Id)"
if (-not $p.WaitForExit(360000)) { Stop-Process -Id $p.Id -Force; throw 'Owned Maya timed out after 360 seconds' }
exit $p.ExitCode
"""
    (out / 'launch.ps1').write_text(launcher, encoding='utf-8')
    env = {k:v for k,v in os.environ.items() if not k.startswith('SURF_')}
    env.update(BREAKER_LOOK_WORKER='1', BREAKER_LOOK_OUT=str(out), BREAKER_LOOK_ARM=args.arm,
               BREAKER_LOOK_SCRIPT=str(Path(__file__).resolve()),
               MAYA_DISABLE_CER='1', MAYA_DISABLE_CIP='1', MAYA_DISABLE_CLIC_IPM='1')
    start = time.perf_counter()
    try:
        with (out / 'launch.log').open('w', encoding='utf-8') as log:
            proc = subprocess.run(['powershell.exe','-NoProfile','-ExecutionPolicy','Bypass',
                                   '-File',str(out/'launch.ps1')],env=env,cwd=ROOT,
                                  stdout=log,stderr=subprocess.STDOUT,timeout=390)
        receipt['exit_code'] = proc.returncode
        assert proc.returncode == 0, 'Maya failed: see worker_result.json and launch.log'
        result = json.loads((out/'worker_result.json').read_text())
        assert result['status'] == 'passed'
        assert sha(SOURCE) == SOURCE_SHA
        assert cache_hashes == {p.name:sha(p) for p in CACHE.iterdir() if p.is_file()}
        receipt.update(status='passed', worker=result, png_sha256=sha(out/'beauty.png'),
                       exr_sha256=sha(out/result['exr']), scene_sha256=sha(out/'candidate.mb'))
    except BaseException:
        receipt.update(status='failed', error=traceback.format_exc())
        raise
    finally:
        receipt['process_s'] = time.perf_counter()-start
        save(out/'receipt.json',receipt)
        print(json.dumps(receipt,indent=2),flush=True)


if __name__ == '__main__':
    if os.environ.get('BREAKER_LOOK_WORKER') == '1':
        worker(Path(os.environ['BREAKER_LOOK_OUT']),os.environ['BREAKER_LOOK_ARM'])
    else:
        main()
