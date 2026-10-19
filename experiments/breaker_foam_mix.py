"""One authorized 24.9 s oblique preview: 25% A + 75% B volume density.

Reuse the validated frozen scene, controls and quadrature points. No simulation,
reconstruction, shader change or new foam support. Only C is ray traced.
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
sys.path.insert(0, str(Path(__file__).resolve().parent))
import breaker_deep_water as base
import breaker_foam_layer as layer
import breaker_foam_render as render

OUT = base.ROOT / 'out/runs/breaker_foam_mix_20261019'
SOURCE = base.ROOT / 'out/runs/breaker_foam_layer_20261019_v6'
VIEW = base.ROOT / 'out/runs/breaker_foam_views_20261019/oblique'
ALPHA = .75


def read(path):
    return json.loads(path.read_text())


def controls():
    assert read(SOURCE / 'analysis.json')['status'] == 'passed'
    prior = read(SOURCE / 'layer_check.json')
    assert base.sha(SOURCE / 'layer.npz') == prior['layer_sha256']
    hashes = {}
    for arm in ['A', 'B']:
        folder = VIEW / arm
        receipt = read(folder / 'receipt.json')
        assert receipt['status'] == read(folder / 'image_check.json')['status'] == 'passed'
        for name, expected in [('beauty.png', receipt['png_sha256']),
                               (receipt['worker']['exr'], receipt['exr_sha256']),
                               ('candidate.mb', read(folder / 'export_check.json')['scene_sha256']),
                               ('foam.ass', prior['exports'][arm]['ass_sha256'])]:
            assert base.sha(folder / name) == expected, (arm, name)
            hashes[f'{arm}/{name}'] = expected
    a, b = [read(VIEW / arm / 'state.json') for arm in ['A', 'B']]
    assert {k for k in a if a[k] != b[k]} == {'foam_arm', 'foam_ass_sha256'}
    return hashes


def build():
    assert not base.git('status', '--porcelain'), 'Commit instrument first'
    assert not OUT.exists()
    hashes = controls()
    OUT.mkdir(parents=True)
    (OUT / 'C').mkdir()
    with np.load(SOURCE / 'layer.npz') as z:
        arrays = {k: z[k] for k in z.files}
    a, b = arrays['density_A_per_m'], arrays['density_B_per_m']
    # Endpoint reproduction is the offline control; no duplicate render is spent.
    assert np.array_equal(a * 1 + b * 0, a)
    assert np.array_equal(a * 0 + b * 1, b)
    c = (1 - ALPHA) * a + ALPHA * b
    assert np.isfinite(c).all() and (c >= 0).all()
    assert np.array_equal(c > 0, a > 0)
    assert np.allclose(c[b == 0], (1 - ALPHA) * a[b == 0])
    assert np.all((c >= np.minimum(a, b)) & (c <= np.maximum(a, b)))
    budgets = [float(d.sum() * 4 * np.pi / 3 * layer.RADIUS**3) for d in [a, b, c]]
    assert max(abs(v / budgets[0] - 1) for v in budgets) < 1e-12
    arrays['density_C_per_m'] = c
    np.savez_compressed(OUT / 'layer.npz', **arrays)
    exported = layer.export_points(OUT / 'C/foam.ass', arrays['points_m'], c)
    # Inspect all three fields on identical source samples and the same color scale.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    s = arrays['surface_m']
    selected = (s[:, 0] > 142) & (s[:, 0] < 150) & (s[:, 1] > 12) & (s[:, 1] < 17)
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)
    for ax, d, label in zip(axes, [a, c, b], ['A uniform', 'C 25% A + 75% B', 'B clustered']):
        factor = np.divide(d, a, out=np.zeros_like(d), where=a > 0)
        scatter = ax.scatter(s[selected, 0], s[selected, 1], s=.6, c=factor[selected],
                             vmin=0, vmax=4.5, cmap='inferno', rasterized=True)
        ax.set(xlim=(142, 150), ylim=(12, 17), aspect='equal', title=label,
               xlabel='x / simulation m', ylabel='y / simulation m')
    fig.colorbar(scatter, ax=axes, label='Local density / A density; point plot, not a render')
    fig.savefig(OUT / 'density_check.png', dpi=140)
    plt.close(fig)
    base.save(OUT / 'layer_check.json', dict(status='passed', alpha_B=ALPHA,
              equation='density_C = 0.25 * density_A + 0.75 * density_B',
              endpoint_control='A and B reproduced exactly; existing scene/render hashes verified',
              source_hashes=hashes, source_layer_sha256=base.sha(SOURCE / 'layer.npz'),
              layer_sha256=base.sha(OUT / 'layer.npz'), integrated_extinction_m2=budgets,
              density_export=exported, support='Exactly A support; unchanged positions and radius',
              git_sha=base.git('rev-parse', 'HEAD'), simulation=False, render_count_authorized=1))
    print(json.dumps(read(OUT / 'layer_check.json'), indent=2))


def prepare_worker():
    import maya.standalone
    maya.standalone.initialize(name='python')
    import maya.cmds as c
    folder = OUT / 'C'
    try:
        c.file(str(VIEW / 'B/candidate.mb'), open=True, force=True)
        before = read(VIEW / 'B/state.json')
        assert not c.ls(type='bifrostGraphShape')
        assert base.geometry(c) == before['geometry']
        state = {node: {attr: c.getAttr(node + '.' + attr) for attr in before[node]}
                 for node in render.NODE_KEYS}
        state = json.loads(json.dumps(state))
        for node in render.NODE_KEYS:
            assert state[node] == before[node], node
        assert c.listConnections('waterSG.surfaceShader') == ['water']
        assert not c.getAttr('baked_raw.visibility') and not c.getAttr('baked_no_mcf.visibility')
        c.setAttr('foamLayerShape.dso', (folder / 'foam.ass').as_posix(), type='string')
        c.setAttr('opticalEXR.prefix', (folder / 'optics').as_posix(), type='string')
        c.setAttr('defaultArnoldRenderOptions.log_filename', (folder / 'arnold.log').as_posix(), type='string')
        state.update(geometry=base.geometry(c), foam_ass_sha256=base.sha(folder / 'foam.ass'),
                     foam_arm='C', foam_layer_sha256=base.sha(OUT / 'layer.npz'),
                     mist=False, view=before['view'])
        changed = [k for k in before if before[k] != state[k]]
        assert set(changed) == {'foam_arm', 'foam_ass_sha256', 'foam_layer_sha256'}
        base.save(folder / 'state.json', state)
        c.file(rename=str(folder / 'candidate.mb'))
        c.file(save=True, type='mayaBinary', force=True)
        c.arnoldExportAss(f=str(folder / 'scene_check.ass'), mask=24)
        assert layer.ass_units(folder / 'scene_check.ass') == layer.ass_units(folder / 'foam.ass') == .01
        base.save(folder / 'export_check.json', dict(status='passed',
                  scene_sha256=base.sha(folder / 'candidate.mb'), ass_sha256=base.sha(folder / 'scene_check.ass')))
        base.save(OUT / 'prepare_worker.json', dict(status='passed', changed_state_keys=changed))
        base.exit_worker(0)
    except BaseException:
        base.save(OUT / 'prepare_worker.json', dict(status='failed', error=traceback.format_exc()))
        base.exit_worker(1)


def prepare():
    assert not base.git('status', '--porcelain')
    assert read(OUT / 'layer_check.json')['status'] == 'passed'
    env = os.environ.copy()
    env['FOAM_MIX_WORKER'] = 'prepare'
    with (OUT / 'prepare.log').open('w') as log:
        result = subprocess.run([render.MAYAPY, str(Path(__file__).resolve())], cwd=base.ROOT,
                                env=env, stdout=log, stderr=subprocess.STDOUT, timeout=180)
    assert result.returncode == 0, (OUT / 'prepare.log').read_text()[-3000:]
    assert read(OUT / 'prepare_worker.json')['status'] == 'passed'
    alignment = layer.check_imported(OUT / 'C')
    base.save(OUT / 'prepare.json', dict(status='passed', alignment=alignment))
    print(json.dumps(alignment, indent=2))


def launch():
    assert not base.git('status', '--porcelain')
    assert read(OUT / 'review.json')['status'] == 'passed'
    assert read(OUT / 'prepare.json')['status'] == 'passed'
    assert controls() == read(OUT / 'layer_check.json')['source_hashes']
    folder = OUT / 'C'
    assert base.sha(folder / 'candidate.mb') == read(folder / 'export_check.json')['scene_sha256']
    assert base.sha(folder / 'foam.ass') == read(folder / 'alignment_check.json')['foam_ass_sha256']
    assert base.sha(OUT / 'layer.npz') == read(OUT / 'layer_check.json')['layer_sha256']
    hashes = {p.name: base.sha(p) for p in base.CACHE.iterdir() if p.is_file()}
    expected = read(render.CONTROL / 'receipt.json')['cache_sha256']
    assert hashes == expected
    path = folder / 'receipt.json'
    assert not path.exists(), 'Only one new render authorized'
    launcher = (VIEW / 'B/launch.ps1').read_text().replace('FOAM_VIEW_SCRIPT', 'FOAM_MIX_SCRIPT')
    (folder / 'launch.ps1').write_text(launcher)
    env = {k: v for k, v in os.environ.items() if not k.startswith('SURF_')}
    env.update(FOAM_MIX_WORKER='render', FOAM_MIX_SCRIPT=str(Path(__file__).resolve()),
               MAYA_DISABLE_CER='1', MAYA_DISABLE_CIP='1', MAYA_DISABLE_CLIC_IPM='1')
    receipt = dict(status='running', git_sha=base.git('rev-parse', 'HEAD'), alpha_B=ALPHA,
                   cache_sha256=hashes, started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    base.save(path, receipt)
    started = time.perf_counter()
    try:
        with (folder / 'launch.log').open('w') as log:
            result = subprocess.run(['powershell.exe', '-NoProfile', '-ExecutionPolicy', 'Bypass',
                                     '-File', str(folder / 'launch.ps1')], cwd=base.ROOT, env=env,
                                    stdout=log, stderr=subprocess.STDOUT, timeout=930)
        assert result.returncode == 0, 'Read C/launch.log and C/worker_result.json'
        worker = read(folder / 'worker_result.json')
        assert worker['status'] == 'passed'
        assert hashes == {p.name: base.sha(p) for p in base.CACHE.iterdir() if p.is_file()}
        receipt.update(status='passed', exit_code=result.returncode, worker=worker,
                       png_sha256=base.sha(folder / 'beauty.png'), exr_sha256=base.sha(folder / worker['exr']))
    except BaseException:
        receipt.update(status='failed', error=traceback.format_exc())
        raise
    finally:
        receipt['process_s'] = time.perf_counter() - started
        base.save(path, receipt)
    render.check(folder)
    base.save(OUT / 'analysis.json', dict(status='passed', new_render_count=1,
              image=read(folder / 'image_check.json'), layer=read(OUT / 'layer_check.json'),
              preparation=read(OUT / 'prepare_worker.json'), alignment=read(folder / 'alignment_check.json'),
              simulation=False, remeshing=False, mist=False, source_controls_reused=str(VIEW),
              temporal_validation=False, paper_method_claim=False))


if __name__ == '__main__':
    mode = os.environ.get('FOAM_MIX_WORKER')
    if mode:
        prepare_worker() if mode == 'prepare' else render.render(OUT / 'C')
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('action', choices=['build', 'prepare', 'render'])
        {'build': build, 'prepare': prepare, 'render': launch}[parser.parse_args().action]()
