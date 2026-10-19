"""Bounded optical diagnosis on already baked geometry; never runs a graph/solver.

Each invocation starts one hidden Maya render session (the subscription's Render
Sequence path). A is the conversion control; B changes only specular depth 1->3.
PNG preserves the prior display output. A separate full-float linear EXR keeps
unprocessed component AOVs and a raw beauty LPE, alongside denoised RGBA.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / 'out/runs/breaker_surface_diagnosis_20261019_v3'
PRIOR = ROOT / 'out/runs/breaker_shading_study_20261019/A'
CACHE = ROOT / 'out/runs/wavesea4/cache/breaking_section'
DEFAULT = ROOT / 'out/runs/breaker_optical_aovs_20261019'
AOVS = ['beauty_raw', 'specular', 'specular_direct', 'specular_indirect',
        'transmission', 'transmission_direct', 'transmission_indirect',
        'diffuse', 'sss', 'coat', 'sheen', 'emission', 'volume', 'background', 'N', 'P', 'Z']
# Defined before measuring: conversion may change stochastic samples/roundoff,
# but must remain within one 8-bit level on average, including the curl itself.
CONTROL_LIMITS = {'mean_abs_rgb_levels': 1.0, 'fraction_pixels_over_8_levels': .01}
CURL = [200, 175, 800, 530]


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def git(*args):
    return subprocess.check_output(['git', *args], cwd=ROOT, text=True).strip()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')


def exit_worker(code):
    # The installed Maya/MtoA can crash during plugin teardown. Exit only after
    # validated outputs and a flushed success/failure record, as the mesh exporter.
    import ctypes
    sys.stdout.flush()
    sys.stderr.flush()
    terminate = ctypes.windll.kernel32.TerminateProcess
    terminate.argtypes = [ctypes.c_void_p, ctypes.c_uint]
    terminate.restype = ctypes.c_int
    terminate(ctypes.c_void_p(-1), code)
    raise ctypes.WinError()


def snapshot(c):
    import maya.api.OpenMaya as om
    nodes = ['water', 'sky', 'sun', 'sunShape', 'persp', 'perspShape',
             'defaultArnoldRenderOptions', 'defaultArnoldDriver', 'defaultArnoldDenoiser',
             'defaultColorMgtGlobals', 'opticalEXR']
    nodes += c.ls(type='mesh', long=True)
    result = {}
    for node in nodes:
        values = {}
        for attr in c.listAttr(node, scalar=True, read=True) or []:
            if attr in ('message', 'caching', 'nodeState'):
                continue
            try:
                value = c.getAttr(node + '.' + attr)
                json.dumps(value)
                values[attr] = value
            except (RuntimeError, TypeError, ValueError):
                pass
        result[node] = values
    sel = om.MSelectionList()
    sel.add('baked_control')
    dag = sel.getDagPath(0)
    dag.extendToShape()
    fn = om.MFnMesh(dag)
    result['geometry'] = {'vertices': fn.numVertices, 'triangles': len(fn.getTriangles()[1]) // 3,
                          'bounds': c.exactWorldBoundingBox('baked_control'),
                          'shading_engines': c.listConnections(dag.fullPathName(), type='shadingEngine')}
    result['imager_connections'] = c.listConnections('defaultArnoldRenderOptions.imagers', plugs=True)
    return result


def worker(out, arm, inspect_only=False):
    import maya.cmds as c
    import mtoa.aovs as aovs
    began = time.perf_counter()
    result = {'arm': arm, 'status': 'running'}
    try:
        c.scriptEditorInfo(writeHistory=True, historyFilename=str(out / 'maya.log'))
        c.file(str(SOURCE / 'baked_breaker.mb'), open=True, force=True)
        assert not c.ls(type='bifrostGraphShape'), 'Live graph in frozen scene'
        assert c.getAttr('baked_control.visibility')
        assert not c.getAttr('baked_raw.visibility') and not c.getAttr('baked_no_mcf.visibility')
        c.currentTime(1)
        c.setAttr('defaultRenderGlobals.currentRenderer', 'arnold', type='string')
        for attr, value in {'animation': 1, 'extensionPadding': 4, 'outFormatControl': 0,
                            'putFrameBeforeExt': 1, 'startFrame': 1, 'endFrame': 1}.items():
            c.setAttr('defaultRenderGlobals.' + attr, value)
        c.setAttr('defaultRenderGlobals.imageFilePrefix', 'beauty', type='string')
        c.setAttr('defaultArnoldDriver.ai_translator', 'png', type='string')
        c.setAttr('defaultResolution.width', 960)
        c.setAttr('defaultResolution.height', 540)
        c.setAttr('defaultResolution.deviceAspectRatio', 960 / 540)
        for attr, value in {'AASamples': 4, 'GIVolumeDepth': 3, 'GIVolumeSamples': 2,
                            'threads_autodetect': 1, 'renderDevice': 0,
                            'log_verbosity': 2, 'log_to_file': 1, 'log_to_console': 0}.items():
            c.setAttr('defaultArnoldRenderOptions.' + attr, value)
        assert c.getAttr('defaultArnoldRenderOptions.GISpecularDepth') == 1
        assert c.getAttr('defaultArnoldRenderOptions.GITransmissionDepth') == 8
        if arm == 'B':
            c.setAttr('defaultArnoldRenderOptions.GISpecularDepth', 3)
        c.workspace(fileRule=['images', out.as_posix()])
        c.setAttr('defaultArnoldRenderOptions.log_filename', (out / 'arnold.log').as_posix(), type='string')
        # Preserve RGBA's original OIDN processing; other named layers remain raw.
        c.setAttr('defaultArnoldDenoiser.layerSelection', 'RGBA', type='string')
        c.setAttr('defaultArnoldDenoiser.outputSuffix', '', type='string')
        driver = c.createNode('aiAOVDriver', name='opticalEXR')
        c.setAttr(driver + '.ai_translator', 'exr', type='string')
        c.setAttr(driver + '.mergeAOVs', True)
        c.setAttr(driver + '.halfPrecision', False)
        c.setAttr(driver + '.prefix', (out / 'optics').as_posix(), type='string')
        box = c.createNode('aiAOVFilter', name='opticalBox')
        c.setAttr(box + '.ai_translator', 'box', type='string')
        rgba = aovs.AOVInterface().addAOV('RGBA').node
        c.connectAttr(driver + '.message', rgba + '.outputs[1].driver')
        c.connectAttr('defaultArnoldFilter.message', rgba + '.outputs[1].filter')
        for name in AOVS:
            node = aovs.AOVInterface().addAOV(name, 'rgba' if name == 'beauty_raw' else None).node
            c.connectAttr(driver + '.message', node + '.outputs[0].driver', force=True)
            if name not in ('N', 'P', 'Z'):
                c.connectAttr(box + '.message', node + '.outputs[0].filter', force=True)
            if name == 'beauty_raw':
                c.setAttr(node + '.lightPathExpression', 'C.*', type='string')
        state = snapshot(c)
        assert state['geometry']['triangles'] == 292516
        write_json(out / 'state.json', state)
        result['open_and_setup_s'] = time.perf_counter() - began
        if not inspect_only:
            c.file(rename=str(out / 'optical.mb'))
            c.file(save=True, type='mayaBinary', force=True)
            start = time.perf_counter()
            c.arnoldRender(seq='1', w=960, h=540, cam='persp', srv=False)
            result['render_call_s'] = time.perf_counter() - start
            images = list((out / 'tmp').glob('*.png'))
            assert len(images) == 1, [str(p) for p in images]
            shutil.move(str(images[0]), out / 'beauty.png')
            exrs = list(out.rglob('*.exr'))
            assert len(exrs) == 1, [str(p) for p in exrs]
            result['exr'] = str(exrs[0].relative_to(out))
        result['worker_wall_s'] = time.perf_counter() - began
        result['status'] = 'passed'
        write_json(out / 'worker_result.json', result)
        exit_worker(0)
    except BaseException:
        result['status'] = 'failed'
        result['error'] = traceback.format_exc()
        write_json(out / 'worker_result.json', result)
        exit_worker(1)


def control_check(frame):
    import numpy as np
    from PIL import Image
    old = np.asarray(Image.open(PRIOR / 'breaker.0001.png').convert('RGB'), dtype=np.float64)
    new = np.asarray(Image.open(frame).convert('RGB'), dtype=np.float64)
    assert old.shape == new.shape == (540, 960, 3)
    difference = np.abs(new - old)
    x0, y0, x1, y1 = CURL
    result = {'pixels_exact': bool(np.array_equal(old, new)), 'limits': CONTROL_LIMITS,
              'curl_xyxy': CURL, 'regions': {}}
    for name, diff in [('full', difference), ('curl', difference[y0:y1, x0:x1])]:
        metrics = {'mean_abs_rgb_levels': float(diff.mean()),
                   'fraction_pixels_over_8_levels': float((diff.max(axis=2) > 8).mean()),
                   'max_abs_rgb_levels': float(diff.max())}
        result['regions'][name] = metrics
    result['status'] = 'passed' if all(
        values[key] <= threshold for values in result['regions'].values()
        for key, threshold in CONTROL_LIMITS.items()) else 'failed'
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('arm', choices=['inspect', 'A', 'B'])
    parser.add_argument('--out', type=Path, default=DEFAULT)
    args = parser.parse_args()
    assert not git('status', '--porcelain'), 'Commit instrument before measuring'
    out = args.out.resolve() / args.arm
    assert not out.exists(), f'Refusing to overwrite {out}'
    baked = SOURCE / 'baked_breaker.mb'
    expected = json.loads((SOURCE / 'static_verification.json').read_text())
    assert expected['status'] == 'passed' and digest(baked) == expected['baked_scene_sha256']
    prior = json.loads((PRIOR / 'receipt.json').read_text())
    hashes = {p.name: digest(p) for p in CACHE.iterdir() if p.is_file()}
    assert hashes == prior['input_sha256'] and digest(PRIOR / 'breaker.0001.png') == prior['image_sha256']
    if args.arm == 'B':
        control = json.loads((out.parent / 'A/receipt.json').read_text())
        assert control['status'] == 'passed' and control['control_check']['status'] == 'passed'
        assert control['git_sha'] == git('rev-parse', 'HEAD'), 'Revalidate control after code changes'
    out.mkdir(parents=True)
    receipt = {'arm': args.arm, 'git_sha': git('rev-parse', 'HEAD'), 'status': 'running',
               'started_utc': time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()),
               'cache_sha256': hashes, 'baked_sha256': digest(baked),
               'prior_png_sha256': prior['image_sha256'], 'control_limits': CONTROL_LIMITS,
               'timeout_s': 360, 'requested_specular_depth': 3 if args.arm == 'B' else 1}
    write_json(out / 'receipt.json', receipt)
    env = {k: v for k, v in os.environ.items() if not k.startswith('SURF_')}
    for key in ('MAYA_DISABLE_CER', 'MAYA_DISABLE_CIP', 'MAYA_DISABLE_CLIC_IPM'):
        env[key] = '1'
    env.update(BREAKER_OPTICAL_WORKER='1', BREAKER_OPTICAL_OUT=str(out), BREAKER_OPTICAL_ARM=args.arm)
    launch = r"""
$ErrorActionPreference = 'Stop'
if (Get-Process maya,mayapy -ErrorAction SilentlyContinue) { throw 'Another Maya process is open' }
if ((Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB -lt 20) { throw 'Less than 20 GB RAM free' }
$pythonScript = $env:BREAKER_OPTICAL_SCRIPT -replace '\\', '/'
$command = '"python(\"import runpy; runpy.run_path(r''' + $pythonScript + ''', run_name=''__main__'')\")"'
$p = Start-Process -FilePath 'C:/Program Files/Autodesk/Maya2027/bin/maya.exe' -ArgumentList @('-command', $command) -PassThru -WindowStyle Hidden
Write-Output "MAYA_PID $($p.Id)"
if (-not $p.WaitForExit(360000)) { Stop-Process -Id $p.Id -Force; throw 'Owned Maya process timed out after 360 seconds' }
exit $p.ExitCode
"""
    # A script file avoids nested command-shell quoting; only our fixed path enters MEL.
    (out / 'launch.ps1').write_text(launch, encoding='utf-8')
    env['BREAKER_OPTICAL_SCRIPT'] = str(Path(__file__).resolve())
    start = time.perf_counter()
    try:
        with (out / 'launch.log').open('w', encoding='utf-8') as log:
            run = subprocess.run(['powershell.exe', '-NoProfile', '-ExecutionPolicy', 'Bypass',
                                  '-File', str(out / 'launch.ps1')], env=env, cwd=ROOT,
                                 stdout=log, stderr=subprocess.STDOUT, timeout=390)
        receipt['exit_code'] = run.returncode
        assert run.returncode == 0, 'Maya failed; see launch.log/worker_result.json'
        result = json.loads((out / 'worker_result.json').read_text())
        assert result['status'] == 'passed', result
        receipt['worker'] = result
        assert digest(baked) == receipt['baked_sha256']
        assert hashes == {p.name: digest(p) for p in CACHE.iterdir() if p.is_file()}
        if args.arm != 'inspect':
            receipt['png_sha256'] = digest(out / 'beauty.png')
            receipt['exr_sha256'] = digest(out / result['exr'])
            if args.arm == 'A':
                receipt['control_check'] = control_check(out / 'beauty.png')
                assert receipt['control_check']['status'] == 'passed', 'Conversion control disagrees; stop before B'
        receipt['status'] = 'passed'
    except BaseException:
        receipt['status'] = 'failed'
        receipt['error'] = traceback.format_exc()
        raise
    finally:
        receipt['process_wall_s'] = time.perf_counter() - start
        write_json(out / 'receipt.json', receipt)
        print(json.dumps(receipt, indent=2), flush=True)


if __name__ == '__main__':
    if os.environ.get('BREAKER_OPTICAL_WORKER') == '1':
        worker(Path(os.environ['BREAKER_OPTICAL_OUT']), os.environ['BREAKER_OPTICAL_ARM'],
               os.environ['BREAKER_OPTICAL_ARM'] == 'inspect')
    else:
        main()
