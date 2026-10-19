"""Six authorized cached-wave views with fixed 20% uniform / 80% clustered foam."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time
import ctypes
import shutil
import traceback

import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parent))
import breaker_deep_water as base
import breaker_surface_correspondence as surface
import breaker_interface_mask as masklib
import breaker_foam_layer as foam
import breaker_foam_render as render
from breaker_foam_views import mesh_flags

OUT = base.ROOT / 'out/runs/breaker_six_views_20261019'
CACHE = base.ROOT / 'out/runs/wavesea4/cache/overview'
OLD = base.ROOT / 'out/runs/breaker_foam_views_20261019'
PROTOTYPE = base.ROOT / 'out/runs/breaker_foam_mix_20261019/C'
FRAMES = {
    'pre': dict(source=CACHE, frame=1, time_s=23.025, target=[129, 13, 68], rotation=[-18, 140, 0], width=24,
                xmin=120, wide_target=[140, 12.5, 44]),
    'break': dict(source=base.ROOT/'out/runs/wavesea4/cache/breaking_overview', frame=1,
                  time_s=24.9, target=[145.5, 13.5, 70], rotation=[-15, 140, 0], width=22,
                  xmin=132, wide_target=[151, 12.5, 44]),
    'shore': dict(source=CACHE, frame=5, time_s=27.025, target=[157, 12.75, 65], rotation=[-30, 90, 0], width=24,
                  xmin=142, wide_target=[159, 13, 44]),
}
ALPHA = .8


def read(path):
    return json.loads(path.read_text())


def control():
    receipt = read(OLD/'frames/control/mesh_receipt.json')
    assert receipt['status'] == 'passed' and receipt['normal_transfer_max_error'] < 1e-5
    assert base.sha(OLD/'frames/control/mesh.npz') == receipt['mesh_sha256']
    z = np.load(OLD/'frames/control/mesh.npz')
    original = np.load(surface.FROZEN)
    assert np.array_equal(z['points'], original['points'])
    def canonical(t):
        t = np.sort(t, axis=1)
        return t[np.lexsort(t.T[::-1])]
    assert np.array_equal(canonical(z['triangles']), canonical(original['triangles']))
    original_normals = np.empty_like(z['points'])
    original_normals[original['face_vertices']] = original['normals'][original['normal_ids']]
    assert np.max(np.abs(z['normals'] - original_normals)) < 1e-5
    assert read(PROTOTYPE/'image_check.json')['status'] == 'passed'
    assert base.sha(PROTOTYPE/'candidate.mb') == read(PROTOTYPE/'export_check.json')['scene_sha256']
    base.save(OUT/'control.json', dict(status='passed', frozen_mesh_reproduced=True,
              normal_control_max_error=float(np.max(np.abs(z['normals']-original_normals))),
              stock_settings=surface.KNOBS, prototype_sha256=base.sha(PROTOTYPE/'candidate.mb')))


class MemoryStatus(ctypes.Structure):
    _fields_ = [('length', ctypes.c_ulong), ('load', ctypes.c_ulong)] + [
        (name, ctypes.c_ulonglong) for name in
        ['total_physical', 'available_physical', 'total_page', 'available_page',
         'total_virtual', 'available_virtual', 'available_extended']]


def free_gb():
    state = MemoryStatus()
    state.length = ctypes.sizeof(state)
    assert ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(state))
    return state.available_physical / 2**30


def mesh(name):
    sha = surface.clean_sha()
    assert read(OUT/'control.json')['status'] == 'passed'
    spec = FRAMES[name]
    folder = OUT/'frames'/name
    assert not folder.exists()
    folder.mkdir(parents=True)
    (folder/'cache').mkdir()
    p = np.load(spec['source']/f"liquid.{spec['frame']:04}.npy", mmap_mode='r')
    keep = (p[:, 0] >= spec['xmin']*4) & (p[:, 2] >= 10*4) & (p[:, 2] < 78*4)
    hashes = {}
    for kind in ['liquid', 'radius', 'speed']:
        source = spec['source']/f"{kind}.{spec['frame']:04}.npy"
        hashes[kind] = base.sha(source)
        np.save(folder/'cache'/f'{kind}.0001.npy', np.load(source, mmap_mode='r')[keep])
    count = int(keep.sum())
    del p, keep
    base.save(folder/'graph.json', surface.graph(1))
    np.save(folder/'probes.npy', np.array([[582, 54, 280]], np.float32))
    inputs = dict(surface.KNOBS, liquid_pattern=(folder/'cache/liquid.####').as_posix(),
                  radius_pattern=(folder/'cache/radius.####').as_posix(), frame=1,
                  path_probes=(folder/'probes.npy').as_posix())
    for key in ['positions', 'indices', 'offsets', 'sdf']:
        inputs['path_'+key] = (folder/(key+'.npy')).as_posix()
    command = [str(surface.BIF/'bin/bifcmd.exe'), str(folder/'graph.json')]
    for key, value in inputs.items():
        command += ['--set-port', key, str(value)]
    env = os.environ.copy()
    env['PATH'] = str(surface.BIF/'bin')+';'+str(surface.BIF/'thirdparty/bin')+';'+env['PATH']
    receipt = dict(status='running', git_sha=sha, time_s=spec['time_s'], source=str(spec['source']),
                   source_frame=spec['frame'], source_hashes=hashes, inputs=inputs,
                   particles=count, crop_x_min_m=spec['xmin'], crop_z_m=[10, 78], simulation=False)
    base.save(folder/'mesh_receipt.json', receipt)
    started = time.perf_counter()
    minimum = free_gb()
    assert minimum > 38, f'Insufficient initial RAM for wider reconstruction: {minimum:.1f} GiB'
    with (folder/'bifcmd.log').open('w') as log:
        process = subprocess.Popen(command, cwd=base.ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        while process.poll() is None:
            elapsed = time.perf_counter()-started
            available = free_gb()
            minimum = min(minimum, available)
            base.save(folder/'progress.json', dict(pid=process.pid, elapsed_s=elapsed, free_gib=available,
                                                  minimum_free_gib=minimum))
            if available < 5 or elapsed > 1200:
                process.kill()
                process.wait()
                receipt.update(status='failed', reason='Owned process stopped at RAM/time guard',
                               free_gib=available, elapsed_s=elapsed)
                base.save(folder/'mesh_receipt.json', receipt)
                raise RuntimeError(receipt)
            time.sleep(2)
    receipt.update(exit_code=process.returncode, process_s=time.perf_counter()-started,
                   minimum_free_gib=minimum)
    base.save(folder/'mesh_receipt.json', receipt)
    assert process.returncode == 0
    text = (folder/'bifcmd.log').read_text()
    for key in ['positions', 'indices', 'offsets', 'sdf']:
        assert f"'ok_{key}': (bool) true" in text
    positions = np.load(folder/'positions.npy')
    triangles = np.load(folder/'indices.npy').reshape(-1, 3).astype(np.int64)
    offsets = np.load(folder/'offsets.npy')
    assert np.all(np.diff(offsets) == 3) and np.isfinite(positions).all()
    assert triangles.min() >= 0 and triangles.max() < len(positions)
    env['SIX_WORKER'] = 'normals'
    env['SIX_FRAME'] = name
    with (folder/'normal.log').open('w') as log:
        proc = subprocess.run([render.MAYAPY, str(Path(__file__).resolve())], env=env,
                              stdout=log, stderr=subprocess.STDOUT, timeout=180)
    assert proc.returncode == 0
    normals = np.load(folder/'normals.npy')
    assert normals.shape == positions.shape and np.isfinite(normals).all()
    assert np.max(np.abs(np.linalg.norm(normals, axis=1)-1)) < 1e-4
    np.savez_compressed(folder/'mesh.npz', points=positions, triangles=triangles, normals=normals)
    for kind, expected in hashes.items():
        assert base.sha(spec['source']/f"{kind}.{spec['frame']:04}.npy") == expected
    receipt.update(status='passed', vertices=len(positions), triangles=len(triangles),
                   mesh_sha256=base.sha(folder/'mesh.npz'))
    base.save(folder/'mesh_receipt.json', receipt)
    print(json.dumps(receipt, indent=2))


def normal_worker(name):
    import maya.standalone
    maya.standalone.initialize(name='python')
    import maya.api.OpenMaya as om
    folder = OUT/'frames'/name
    p = np.load(folder/'positions.npy')
    t = np.load(folder/'indices.npy').reshape(-1, 3)
    fn = om.MFnMesh()
    fn.create(om.MPointArray(p.tolist()), [3]*len(t), t.ravel().tolist())
    _, ids = fn.getNormalIds()
    corners = np.asarray(fn.getNormals())[np.asarray(ids)]
    normals = np.empty_like(p)
    normals[t.ravel()] = corners
    assert np.max(np.abs(corners-normals[t.ravel()])) < 1e-5
    np.save(folder/'normals.npy', normals)
    base.exit_worker(0)


def mask(name):
    surface.clean_sha()
    source = OUT/'frames'/name
    folder = source/'mask'
    assert read(source/'mesh_receipt.json')['status'] == 'passed'
    assert not folder.exists()
    folder.mkdir()
    z = np.load(source/'mesh.npz')
    p, t = z['points'].astype(float)*.25, z['triangles']
    h, area, n, edges, counts, ta = masklib.curvature(p, t)
    assert np.all(counts == 2)
    filtered = masklib.filter_field(h, p, edges)
    lo, hi = p.min(0), p.max(0)
    cut = np.minimum.reduce([p[:, 0]-lo[0], hi[0]-p[:, 0], p[:, 2]-lo[2], hi[2]-p[:, 2]])
    m = masklib.smoothstep(2, 4, np.abs(filtered))*masklib.smoothstep(.6, 1, cut)*masklib.smoothstep(12.25, 13, p[:, 1])
    tm = m[t].mean(1)
    assert tm[cut[t].min(1) < .5].max() == 0
    np.savez_compressed(folder/'mask.npz', points_m=p, triangles=t, normals=n, curvature=filtered,
                        raw_curvature=h, mask=m, cut_distance_m=cut)
    masklib.diagrams(p, t, m, folder, center=FRAMES[name]['target'], width_m=FRAMES[name]['width'])
    base.save(folder/'mask_receipt.json', dict(status='numeric_passed_visual_review_pending',
              selected_area_m2=float(ta[tm > .05].sum()), strong_area_m2=float(ta[tm > .5].sum()),
              cut_max=0, height_gate_m=[12.25, 13], time_s=FRAMES[name]['time_s'],
              mesh_sha256=base.sha(source/'mesh.npz')))
    print(json.dumps(read(folder/'mask_receipt.json'), indent=2))


def build_layer(name):
    source = OUT/'frames'/name
    foam.MASK = source/'mask'
    foam.VIEW_CENTER = FRAMES[name]['target']
    foam.VIEW_WIDTH = FRAMES[name]['width']
    folder = source/'foam'
    foam.build(folder)
    z = np.load(folder/'layer.npz')
    arrays = {key: z[key] for key in z.files}
    arrays['density_C_per_m'] = .2*arrays['density_A_per_m'] + .8*arrays['density_B_per_m']
    assert np.array_equal(arrays['density_C_per_m'] > 0, arrays['density_A_per_m'] > 0)
    budgets = [arrays['density_'+arm+'_per_m'].sum() for arm in ['A', 'B', 'C']]
    assert max(abs(value/budgets[0]-1) for value in budgets) < 1e-12
    np.savez_compressed(folder/'layer.npz', **arrays)
    (folder/'C').mkdir()
    exported = foam.export_points(folder/'C/foam.ass', arrays['points_m'], arrays['density_C_per_m'])
    record = read(folder/'layer_check.json')
    record.update(layer_sha256=base.sha(folder/'layer.npz'), alpha_B=ALPHA)
    record['exports']['C'] = exported
    base.save(folder/'layer_check.json', record)
    print(json.dumps(exported, indent=2))


def view_spec(name, view):
    spec = FRAMES[name]
    return dict(time_s=spec['time_s'], width_m=56 if view == 'wide' else spec['width'],
                rotation=[-65, 90, 0] if view == 'wide' else spec['rotation'],
                target_m=spec['wide_target'] if view == 'wide' else spec['target'])


def camera_plot(name):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    z = np.load(OUT/'frames'/name/'mesh.npz')
    p = z['points']*.25
    fig, axes = plt.subplots(1, 2, figsize=(16, 5), constrained_layout=True)
    records = {}
    for ax, view in zip(axes, ['wide', 'detail']):
        spec = view_spec(name, view)
        pitch, yaw, _ = np.radians(spec['rotation'])
        back = np.array([np.sin(yaw)*np.cos(pitch), -np.sin(pitch), np.cos(yaw)*np.cos(pitch)])
        right = np.array([np.cos(yaw), 0, -np.sin(yaw)])
        up = np.cross(back, right)
        xy = (p-spec['target_m']) @ np.stack([right, up], axis=1)
        half = spec['width_m']/2
        visible = (np.abs(xy[:, 0]) < half) & (np.abs(xy[:, 1]) < half*540/960)
        edge = (p[:, 0] < FRAMES[name]['xmin']+.4) | (p[:, 2] < 10.4) | (p[:, 2] > 77.6)
        above_bed = p[:, 1] > (p[:, 0]-36)*.1+.3
        stride = max(1, len(p)//120000)
        ax.scatter(xy[::stride, 0], xy[::stride, 1], s=.1, c=p[::stride, 1], cmap='Blues', vmin=10, vmax=17)
        marked = edge & above_bed & visible
        ax.scatter(xy[marked, 0], xy[marked, 1], s=.5, c='red')
        ax.set(xlim=(-half, half), ylim=(-half*540/960, half*540/960), aspect='equal',
               title=f'{name} / {view}: projected mesh; red = boundary candidates', xlabel='camera horizontal / m', ylabel='camera vertical / m')
        records[view] = dict(projected_boundary_candidates=int(marked.sum()),
                             note='Projection only, not occlusion or image coverage')
    fig.savefig(OUT/'frames'/name/'camera_check.png', dpi=120)
    plt.close(fig)
    base.save(OUT/'frames'/name/'camera_check.json', records)


def prepare_worker(name):
    import maya.standalone
    maya.standalone.initialize(name='python')
    import maya.cmds as c
    import maya.api.OpenMaya as om
    source = OUT/'frames'/name
    before = read(PROTOTYPE/'state.json')
    z = np.load(source/'mesh.npz')
    points, triangles, normals = z['points'], z['triangles'], z['normals']
    for view in ['wide', 'detail']:
        parent_folder = OUT/'shots'/(name+'_'+view)
        folder = parent_folder/'C'
        spec = view_spec(name, view)
        c.file(str(PROTOTYPE/'candidate.mb'), open=True, force=True)
        assert base.geometry(c) == before['geometry'] and not c.ls(type='bifrostGraphShape')
        old_shape = c.listRelatives('baked_control', shapes=True, fullPath=True)[0]
        flags = mesh_flags(c, old_shape)
        c.delete('baked_control')
        transform = c.createNode('transform', name='baked_control')
        sel = om.MSelectionList()
        sel.add(transform)
        fn = om.MFnMesh()
        fn.create(om.MPointArray(points.tolist()), [3]*len(triangles), triangles.ravel().tolist(), parent=sel.getDependNode(0))
        fn.setVertexNormals(om.MVectorArray(normals.tolist()), list(range(len(points))))
        shape = c.listRelatives(transform, shapes=True, fullPath=True)[0]
        for attr, item in flags.items():
            if item['type'] == 'string':
                c.setAttr(shape+'.'+attr, item['value'], type='string')
            else:
                c.setAttr(shape+'.'+attr, item['value'])
        assert mesh_flags(c, shape) == flags
        c.sets(transform, e=True, forceElement='waterSG')
        _, ids = fn.getNormalIds()
        corner_error = float(np.max(np.abs(np.asarray(fn.getNormals())[np.asarray(ids)]-normals[triangles.ravel()])))
        assert corner_error < 1e-5
        assert np.array_equal(np.asarray(fn.getPoints())[:, :3].astype(np.float32), points)
        # Extend only the existing beach's transverse span; retain slope and shader.
        collider_shapes = c.listRelatives('colliders', allDescendents=True, type='mesh', fullPath=True)
        assert len(collider_shapes) == 1
        selection = om.MSelectionList()
        selection.add(collider_shapes[0])
        beach = om.MFnMesh(selection.getDagPath(0))
        old_points = beach.getPoints(om.MSpace.kWorld)
        xyz = np.asarray(old_points)[:, :3].copy()
        assert np.allclose([xyz[:, 2].min(), xyz[:, 2].max()], [240, 312])
        for point in old_points:
            point.z = (point.z-240)*320/72
        beach.setPoints(old_points, om.MSpace.kWorld)
        expanded = np.asarray(beach.getPoints(om.MSpace.kWorld))[:, :3]
        assert np.array_equal(xyz[:, :2], expanded[:, :2])
        assert np.allclose([expanded[:, 2].min(), expanded[:, 2].max()], [0, 320])
        pitch, yaw, _ = np.radians(spec['rotation'])
        direction = np.array([np.sin(yaw)*np.cos(pitch), -np.sin(pitch), np.cos(yaw)*np.cos(pitch)])
        eye = np.asarray(spec['target_m'])*4 + 500*direction
        c.setAttr('persp.translate', *eye, type='double3')
        c.setAttr('persp.rotate', *spec['rotation'], type='double3')
        c.setAttr('perspShape.orthographic', True)
        c.setAttr('perspShape.orthographicWidth', spec['width_m']*4)
        c.setAttr('foamLayerShape.dso', (folder/'foam.ass').as_posix(), type='string')
        c.setAttr('opticalEXR.prefix', (folder/'optics').as_posix(), type='string')
        c.setAttr('defaultArnoldRenderOptions.log_filename', (folder/'arnold.log').as_posix(), type='string')
        state = {node: {a: c.getAttr(node+'.'+a) for a in before[node]} for node in render.NODE_KEYS}
        state = json.loads(json.dumps(state))
        for node in render.NODE_KEYS:
            if node not in ['persp', 'perspShape']:
                assert state[node] == before[node], node
        state.update(geometry=base.geometry(c), foam_ass_sha256=base.sha(folder/'foam.ass'),
                     foam_layer_sha256=base.sha(parent_folder/'layer.npz'), foam_arm='C', alpha_B=ALPHA,
                     mist=False, view=spec, beach_bounds=c.exactWorldBoundingBox('colliders'))
        base.save(folder/'state.json', state)
        c.file(rename=str(folder/'candidate.mb'))
        c.file(save=True, type='mayaBinary', force=True)
        c.arnoldExportAss(f=str(folder/'scene_check.ass'), mask=24)
        base.save(folder/'export_check.json', dict(status='passed', scene_sha256=base.sha(folder/'candidate.mb'),
                  ass_sha256=base.sha(folder/'scene_check.ass'), normal_corner_max_error=corner_error,
                  mesh_render_flags_preserved=True, beach_change='z extent 60..78 m to 0..80 m; x/y unchanged'))
    base.save(source/'prepare_worker.json', dict(status='passed'))
    base.exit_worker(0)


def prepare(name):
    surface.clean_sha()
    source = OUT/'frames'/name
    assert read(source/'foam/review.json')['status'] == 'passed'
    for view in ['wide', 'detail']:
        folder = OUT/'shots'/(name+'_'+view)
        assert not folder.exists()
        (folder/'C').mkdir(parents=True)
        shutil.copy2(source/'foam/layer.npz', folder/'layer.npz')
        shutil.copy2(source/'foam/C/foam.ass', folder/'C/foam.ass')
    env = os.environ.copy()
    env.update(SIX_WORKER='prepare', SIX_FRAME=name)
    with (source/'prepare.log').open('w') as log:
        result = subprocess.run([render.MAYAPY, str(Path(__file__).resolve())], env=env,
                                stdout=log, stderr=subprocess.STDOUT, timeout=240)
    assert result.returncode == 0, (source/'prepare.log').read_text()[-3000:]
    assert read(source/'prepare_worker.json')['status'] == 'passed'
    for view in ['wide', 'detail']:
        foam.check_imported(OUT/'shots'/(name+'_'+view)/'C')
    base.save(source/'prepare.json', dict(status='passed'))


def launch(name, view):
    sha = surface.clean_sha()
    source = OUT/'frames'/name
    folder = OUT/'shots'/(name+'_'+view)/'C'
    assert read(source/'prepare.json')['status'] == 'passed'
    assert read(source/'camera_review.json')['status'] == 'passed'
    assert base.sha(folder/'foam.ass') == read(folder/'alignment_check.json')['foam_ass_sha256']
    assert base.sha(folder/'candidate.mb') == read(folder/'export_check.json')['scene_sha256']
    path = folder/'receipt.json'
    assert not path.exists(), 'One call per authorized image'
    finished = list((OUT/'shots').glob('*/C/receipt.json'))
    assert len(finished) < 6
    for previous in finished:
        assert read(previous)['status'] == 'passed'
        assert read(previous.parent/'image_check.json')['status'] == 'passed'
    spec = FRAMES[name]
    hashes = read(source/'mesh_receipt.json')['source_hashes']
    for kind, expected in hashes.items():
        assert base.sha(spec['source']/f"{kind}.{spec['frame']:04}.npy") == expected
    launcher = (PROTOTYPE/'launch.ps1').read_text().replace('FOAM_MIX_SCRIPT', 'SIX_SCRIPT')
    (folder/'launch.ps1').write_text(launcher)
    env = {k: v for k, v in os.environ.items() if not k.startswith('SURF_')}
    env.update(SIX_WORKER='render', SIX_FRAME=name, SIX_VIEW=view, SIX_SCRIPT=str(Path(__file__).resolve()),
               MAYA_DISABLE_CER='1', MAYA_DISABLE_CIP='1', MAYA_DISABLE_CLIC_IPM='1')
    receipt = dict(status='running', git_sha=sha, alpha_B=ALPHA, view=view_spec(name, view),
                   source_hashes=hashes, started_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime()))
    base.save(path, receipt)
    started = time.perf_counter()
    try:
        with (folder/'launch.log').open('w') as log:
            result = subprocess.run(['powershell.exe', '-NoProfile', '-ExecutionPolicy', 'Bypass',
                                     '-File', str(folder/'launch.ps1')], cwd=base.ROOT, env=env,
                                    stdout=log, stderr=subprocess.STDOUT, timeout=930)
        assert result.returncode == 0
        worker = read(folder/'worker_result.json')
        assert worker['status'] == 'passed'
        receipt.update(status='passed', exit_code=0, worker=worker,
                       png_sha256=base.sha(folder/'beauty.png'), exr_sha256=base.sha(folder/worker['exr']))
    except BaseException:
        receipt.update(status='failed', error=traceback.format_exc())
        raise
    finally:
        receipt['process_s'] = time.perf_counter()-started
        base.save(path, receipt)
    render.check(folder)


def report():
    from PIL import Image, ImageDraw, ImageFont
    font = ImageFont.truetype('C:/Windows/Fonts/arial.ttf', 17)
    sheet = Image.new('RGB', (960, 930), (20, 23, 28))
    draw = ImageDraw.Draw(sheet)
    results = {}
    for row, name in enumerate(FRAMES):
        for col, view in enumerate(['wide', 'detail']):
            key = name+'_'+view
            folder = OUT/'shots'/key/'C'
            check = read(folder/'image_check.json')
            assert check['status'] == 'passed'
            results[key] = dict(view=view_spec(name, view), image=check,
                                state=read(folder/'state.json'), alignment=read(folder/'alignment_check.json'))
            photo = Image.open(folder/'beauty_view.png').convert('RGB')
            sheet.paste(photo.resize((480, 270), Image.Resampling.LANCZOS), (col*480, row*310+40))
            draw.text((col*480+8, row*310+8), f"{FRAMES[name]['time_s']:g} s | {view} | 20% A + 80% B", fill='white', font=font)
    sheet.save(OUT/'comparison.png')
    base.save(OUT/'analysis.json', dict(status='passed', new_render_count=6, simulation=False, mist=False,
              alpha_B=ALPHA, shots=results, frames={name: read(OUT/'frames'/name/'mesh_receipt.json') for name in FRAMES},
              limitation='Static curvature/height mask, no transported foam history; no paper-fidelity claim for procedural density.'))


def scout():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    OUT.mkdir(parents=True, exist_ok=True)
    meta = json.loads((CACHE / 'cache.json').read_text())
    fig, axes = plt.subplots(4, 2, figsize=(16, 13), constrained_layout=True)
    records = []
    for index, (ax, when) in enumerate(zip(axes.ravel(), meta['times_s']), 1):
        p = np.load(CACHE / f'liquid.{index:04}.npy', mmap_mode='r')
        r = np.load(CACHE / f'radius.{index:04}.npy', mmap_mode='r')
        sel = (p[:, 0] >= 480) & (p[:, 2] >= 264) & (p[:, 2] < 288) & (r == 0)
        xyz = p[sel] * .25
        hist, x, y = np.histogram2d(xyz[:, 0], xyz[:, 1], bins=[560, 180], range=[[120, 172], [7, 20]])
        ax.imshow(np.log1p(hist.T), origin='lower', extent=[120, 172, 7, 20], cmap='Blues', aspect='equal')
        ax.plot([120, 172], [8.4, 13.6], color='sienna', linewidth=1)
        ax.axvline(156, color='sienna', linestyle=':', linewidth=1)
        ax.set(title=f'{when:g} s | particle section z = 66..72 m', xlabel='x / m', ylabel='y / m')
        records.append(dict(time_s=when, frame=index, total_particles=len(p),
                            full_width_x120_particles=int((p[:, 0] >= 480).sum()), section_particles=len(xyz)))
    fig.savefig(OUT / 'particle_scout.png', dpi=130)
    plt.close(fig)
    base.save(OUT / 'scout.json', records)
    print(json.dumps(records, indent=2))


if __name__ == '__main__':
    mode = os.environ.get('SIX_WORKER')
    if mode:
        try:
            name = os.environ['SIX_FRAME']
            if mode == 'normals':
                normal_worker(name)
            elif mode == 'prepare':
                prepare_worker(name)
            else:
                render.render(OUT/'shots'/(name+'_'+os.environ['SIX_VIEW'])/'C')
        except BaseException:
            traceback.print_exc()
            base.exit_worker(1)
    else:
        parser = argparse.ArgumentParser(description=__doc__)
        parser.add_argument('action', choices=['scout', 'control', 'mesh', 'mask', 'layer', 'camera', 'prepare', 'render', 'report'])
        parser.add_argument('name', choices=list(FRAMES), nargs='?')
        parser.add_argument('view', choices=['wide', 'detail'], nargs='?')
        args = parser.parse_args()
        if args.action in ['scout', 'control', 'report']:
            globals()[args.action]()
        elif args.action == 'render':
            assert args.name and args.view
            launch(args.name, args.view)
        else:
            assert args.name
            {'mesh': mesh, 'mask': mask, 'layer': build_layer, 'camera': camera_plot, 'prepare': prepare}[args.action](args.name)
