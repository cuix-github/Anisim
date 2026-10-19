"""Memory-bounded stock reconstruction with overlap and a watertight seam gate.

No change to particle radii, voxel spacing or smoothing. Each side has six metres
of context across the join. Only matching existing vertices are welded; no new
surface or interpolation geometry is invented.
"""
import itertools
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
import breaker_six_views as study
import breaker_deep_water as base
import breaker_surface_correspondence as surface

TILES = [(10, 50), (38, 62), (50, 78)]
JOINS = [44, 56]


def stitch(p0, t0, p1, t1, plane):
    a = t0[p0[t0, 2].mean(1) < plane]
    b = t1[p1[t1, 2].mean(1) >= plane]
    used0, inv0 = np.unique(a, return_inverse=True)
    used1, inv1 = np.unique(b, return_inverse=True)
    p0, a = p0[used0], inv0.reshape(-1, 3)
    p1, b = p1[used1], inv1.reshape(-1, 3)
    width, tolerance = .001, .0001  # scene cells: max 25 micrometres in simulation metres
    lookup = {}
    for index in np.flatnonzero(np.abs(p0[:, 2]-plane) < 2):
        key = tuple(np.floor(p0[index]/width).astype(np.int64))
        lookup.setdefault(key, []).append(index)
    remap = np.full(len(p1), -1, np.int64)
    maximum = 0.
    matches = 0
    for index in np.flatnonzero(np.abs(p1[:, 2]-plane) < 2):
        key = np.floor(p1[index]/width).astype(np.int64)
        choices = []
        for offset in itertools.product([-1, 0, 1], repeat=3):
            for other in lookup.get(tuple(key+offset), []):
                error = np.linalg.norm(p0[other].astype(float)-p1[index])
                if error <= tolerance:
                    choices.append((other, error))
        assert len(choices) <= 1, 'Ambiguous seam correspondence'
        if choices:
            remap[index], error = choices[0]
            maximum = max(maximum, float(error))
            matches += 1
    missing = remap < 0
    remap[missing] = np.arange(missing.sum())+len(p0)
    p = np.concatenate([p0, p1[missing]])
    t = np.concatenate([a, remap[b]])
    assert np.all(t[:, 0] != t[:, 1]) and np.all(t[:, 1] != t[:, 2]) and np.all(t[:, 2] != t[:, 0])
    edges = np.concatenate([t[:, [0, 1]], t[:, [1, 2]], t[:, [2, 0]]])
    edges.sort(axis=1)
    _, counts = np.unique(edges, axis=0, return_counts=True)
    assert np.all(counts == 2), f'Seam is not watertight: {int((counts != 2).sum())} invalid edges'
    cross = np.cross(p[t[:, 1]]-p[t[:, 0]], p[t[:, 2]]-p[t[:, 0]])
    assert np.all(np.linalg.norm(cross, axis=1) > 1e-12)
    return p, t, dict(status='passed', welded_vertices=matches, maximum_weld_error_cells=maximum,
                       maximum_weld_error_m=maximum*.25, all_edge_incidence_two=True,
                       join_z_m=plane*.25, no_new_vertices=True)


def control():
    z = np.load(study.OLD/'frames/control/mesh.npz')
    p, t = z['points'], z['triangles']
    plane = 69*4
    center = p[t, 2].mean(1)
    q, f, receipt = stitch(p, t[center < plane+12], p, t[center > plane-12], plane)
    def triangle_positions(p, t):
        vertices = np.sort(t, axis=1)
        return {tuple(sorted(tuple(float(x) for x in p[i]) for i in face)) for face in vertices}
    assert triangle_positions(p, t) == triangle_positions(q, f)
    assert len(q) == len(p) and len(f) == len(t)
    receipt.update(frozen_triangles_exact=True, vertices=len(q), triangles=len(f))
    base.save(study.OUT/'tile_control.json', receipt)
    print(json.dumps(receipt, indent=2))


def run_tile(name, index):
    commit = surface.clean_sha()
    assert study.read(study.OUT/'tile_control.json')['frozen_triangles_exact']
    spec = study.FRAMES[name]
    # The high-lip side needs narrower tiles; each join retains 6 m context.
    zlo, zhi = TILES[index]
    root = study.OUT/'frames'/name
    folder = root/'tiles'/str(index)
    assert not folder.exists()
    (folder/'cache').mkdir(parents=True)
    p = np.load(spec['source']/f"liquid.{spec['frame']:04}.npy", mmap_mode='r')
    keep = (p[:, 0] >= spec['xmin']*4) & (p[:, 2] >= zlo*4) & (p[:, 2] < zhi*4)
    hashes = {}
    for kind in ['liquid', 'radius']:
        source = spec['source']/f"{kind}.{spec['frame']:04}.npy"
        hashes[kind] = base.sha(source)
        np.save(folder/'cache'/f'{kind}.0001.npy', np.load(source, mmap_mode='r')[keep])
    count = int(keep.sum())
    del p, keep
    base.save(folder/'graph.json', surface.graph(1))
    np.save(folder/'probes.npy', np.array([[582, 54, 176]], np.float32))
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
    receipt = dict(status='running', git_sha=commit, source_hashes=hashes, inputs=inputs,
                   particles=count, z_bounds_m=[zlo, zhi], time_s=spec['time_s'])
    base.save(folder/'receipt.json', receipt)
    started = time.perf_counter()
    minimum = study.free_gb()
    assert minimum > 32
    with (folder/'bifcmd.log').open('w') as log:
        process = subprocess.Popen(command, cwd=base.ROOT, env=env, stdout=log, stderr=subprocess.STDOUT)
        while process.poll() is None:
            available = study.free_gb()
            minimum = min(minimum, available)
            elapsed = time.perf_counter()-started
            base.save(folder/'progress.json', dict(pid=process.pid, elapsed_s=elapsed,
                                                  free_gib=available, minimum_free_gib=minimum))
            if available < 5 or elapsed > 1200:
                process.kill()
                process.wait()
                receipt.update(status='failed', reason='RAM/time guard', elapsed_s=elapsed, free_gib=available)
                base.save(folder/'receipt.json', receipt)
                raise RuntimeError(receipt)
            time.sleep(2)
    receipt.update(exit_code=process.returncode, process_s=time.perf_counter()-started, minimum_free_gib=minimum)
    base.save(folder/'receipt.json', receipt)
    assert process.returncode == 0
    text = (folder/'bifcmd.log').read_text()
    for key in ['positions', 'indices', 'offsets', 'sdf']:
        assert f"'ok_{key}': (bool) true" in text
    for kind, expected in hashes.items():
        assert base.sha(spec['source']/f"{kind}.{spec['frame']:04}.npy") == expected
    receipt['status'] = 'passed'
    base.save(folder/'receipt.json', receipt)


def assemble(name):
    surface.clean_sha()
    root = study.OUT/'frames'/name
    tiles = [root/'tiles'/str(i) for i in range(len(TILES))]
    receipts = [study.read(folder/'receipt.json') for folder in tiles]
    assert all(r['status'] == 'passed' for r in receipts)
    meshes = [(np.load(folder/'positions.npy'), np.load(folder/'indices.npy').reshape(-1, 3)) for folder in tiles]
    p, t = meshes[0]
    seams = []
    for (other_p, other_t), plane in zip(meshes[1:], JOINS):
        p, t, seam = stitch(p, t, other_p, other_t, plane=plane*4)
        seams.append(seam)
    base.save(root/'seam_check.json', dict(status='passed', joins=seams))
    np.save(root/'positions.npy', p)
    np.save(root/'indices.npy', t.ravel())
    env = os.environ.copy()
    env.update(SIX_WORKER='normals', SIX_FRAME=name)
    with (root/'normal.log').open('w') as log:
        proc = subprocess.run([study.render.MAYAPY, str(Path(study.__file__).resolve())], env=env,
                              stdout=log, stderr=subprocess.STDOUT, timeout=180)
    assert proc.returncode == 0, (root/'normal.log').read_text()[-3000:]
    normals = np.load(root/'normals.npy')
    assert normals.shape == p.shape and np.isfinite(normals).all()
    assert np.max(np.abs(np.linalg.norm(normals, axis=1)-1)) < 1e-4
    np.savez_compressed(root/'mesh.npz', points=p, triangles=t, normals=normals)
    spec = study.FRAMES[name]
    base.save(root/'mesh_receipt.json', dict(status='passed', time_s=spec['time_s'], source=str(spec['source']),
              source_frame=spec['frame'], source_hashes=receipts[0]['source_hashes'],
              inputs=surface.KNOBS, crop_x_min_m=spec['xmin'], crop_z_m=[10, 78],
              vertices=len(p), triangles=len(t), mesh_sha256=base.sha(root/'mesh.npz'),
              process_s=sum(r['process_s'] for r in receipts), seams=seams, tiles=receipts,
              simulation=False, interpolation_geometry=False))
    print(json.dumps(study.read(root/'mesh_receipt.json'), indent=2))


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['control', 'tile', 'assemble'])
    parser.add_argument('name', choices=list(study.FRAMES), nargs='?')
    parser.add_argument('index', type=int, choices=list(range(len(TILES))), nargs='?')
    args = parser.parse_args()
    if args.action == 'control':
        control()
    elif args.action == 'tile':
        run_tile(args.name, args.index)
    else:
        assemble(args.name)
