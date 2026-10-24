"""One frozen C++ plunge, extended without added damping; no production inputs.

The stopping criteria are an observational near-rest definition, not a proof of
paper fidelity. Every continuation carries the complete state and step counter.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np

from phase_repair_2d import ROOT, NX, NY, D, E, ESCAPE, make, mask, ports
from loop3d import run_bifrost3, load_final_state

DT = float(np.float32(.0025))
CHUNK_STEPS = 8000
EVERY = 80
CRITERIA = dict(window_seconds=10, all_water_rms_cells_s=.05,
                water_speed_p99_cells_s=.2, kinetic_fraction_of_peak=1e-4,
                surface_temporal_range_cells=.25)


def write(path, value):
    path.write_text(json.dumps(value, indent=2), encoding='utf-8')


def config(s):
    p = dict(ESCAPE, **ports(E, D), phase_model=2, phase_alpha=1,
                rho0_face=s.rho0_face, tolerance=s.pressure_tolerance,
                transport_model=0, alpha_liquid=s.alpha[1], alpha_air=s.alpha[0],
                st=0, adapt=0, pressure_solver=0)
    p.update(escape=2, droplet_model=2, emit_model=2)
    return p


def run(out, s, steps, offset=0, state=None, every=EVERY):
    if out.exists():
        raise FileExistsError(out)
    p3 = np.c_[s.pos, np.full(len(s.pos), .5)].astype(np.float32)
    start = time.perf_counter()
    _, _, diags = run_bifrost3('settling_phase', NX, NY, 1,
        lambda x, y, z: mask(x, y), steps, DT,
        seed_fn=lambda *a: (p3, s.typ.copy()), fused=True,
        extra_ports=config(s), init_arrays=state, step_offset=offset,
        final_state=True, dumps=('ph', 'vel', 'rad', 'spr', 'esc'),
        dump_every=every, out_dir=str(out), log_prefix=str(out/'bifcmd'))
    assert np.isfinite(diags).all()
    events = [json.loads(line[len('SPRAY_EVENT '):]) for line in
        (out/'bifcmd.stderr.log').read_text().splitlines()
        if line.startswith('SPRAY_EVENT ')]
    assert len(events) == steps, (len(events), steps)
    write(out/'ledger.json', events)
    receipt = dict(offset=offset, steps=steps, seconds=time.perf_counter()-start,
                   max_pressure_residual=float(diags[:, 1].max()),
                   max_particle_cfl=float(diags[:, 3].max()*DT))
    write(out/'runtime.json', receipt)
    return receipt


def frame(d, j):
    return dict(pos=np.load(d/f'frame.{j:04d}.npy')[:, :2],
                vel=np.load(d/f'vel.{j:04d}.npy')[:, :2],
                typ=np.load(d/f'ph.{j:04d}.npy') > .5,
                spr=np.load(d/f'spr.{j:04d}.npy') > .5,
                esc=np.load(d/f'esc.{j:04d}.npy') > .5,
                rad=np.load(d/f'rad.{j:04d}.npy'))


def measure(a, t):
    primary = ~a['spr']
    water = primary & a['typ']
    bulk = water & ~a['esc']
    speed2 = np.sum(a['vel'].astype(float)**2, axis=1)
    p = a['pos'].astype(float)
    mass = np.where(a['typ'], 1000., 1.)
    surface = []
    # A robust display/measurement proxy, not a reconstructed free surface.
    for left in np.arange(0, NX, 2):
        q = p[bulk & (p[:, 0] >= left) & (p[:, 0] < left+2), 1]
        surface.append(float(np.quantile(q, .95)) if len(q) else None)
    return dict(time=t, water_count=int(water.sum()),
                primary_count=int(primary.sum()), escaped_water=int((water & a['esc']).sum()),
                spray=int(a['spr'].sum()),
                water_rms=float(np.sqrt(speed2[water].mean())),
                bulk_water_rms=float(np.sqrt(speed2[bulk].mean())) if bulk.any() else None,
                water_speed_p99=float(np.quantile(np.sqrt(speed2[water]), .99)),
                water_max_speed=float(np.sqrt(speed2[water].max())),
                air_rms=float(np.sqrt(speed2[primary & ~a['typ']].mean())),
                water_kinetic=float(.5*1000*speed2[water].sum()),
                primary_energy=float(np.sum(mass[primary]*(.5*speed2[primary]+9.800000190734863*p[primary, 1]))),
                water_com=p[water].mean(0).tolist(),
                height_p99=float(np.quantile(p[water, 1], .99)), surface_proxy=surface)


def rest_check(rows):
    tail = [r for r in rows if rows[-1]['time']-r['time'] <= CRITERIA['window_seconds']+1e-6]
    peak = max(r['water_kinetic'] for r in rows)
    surf = np.array([r['surface_proxy'] for r in tail], dtype=float)
    missing = not bool(np.isfinite(surf).all())
    variation = float(np.max(np.ptp(surf, axis=0))) if not missing else None
    checks = dict(long_enough=len(tail) >= 50 and rows[-1]['time'] >= 20-1e-5,
        all_water_rms=max(r['water_rms'] for r in tail) <= CRITERIA['all_water_rms_cells_s'],
        speed_p99=max(r['water_speed_p99'] for r in tail) <= CRITERIA['water_speed_p99_cells_s'],
        kinetic_fraction=max(r['water_kinetic'] for r in tail) <= peak*CRITERIA['kinetic_fraction_of_peak'],
        surface_stationary=not missing and variation <= CRITERIA['surface_temporal_range_cells'])
    return dict(pass_near_rest=all(checks.values()), checks=checks,
                window_start=tail[0]['time'], window_end=tail[-1]['time'],
                max_water_rms=max(r['water_rms'] for r in tail),
                max_speed_p99=max(r['water_speed_p99'] for r in tail),
                max_kinetic_fraction=max(r['water_kinetic'] for r in tail)/peak,
                surface_max_temporal_range=variation)


def compare_final(a, b):
    names = ['final.npy']+[f'final_{n}.npy' for n in
        ('velocities', 'mass', 'momx', 'momy', 'momz', 'tau', 'wt', 'wtph', 'scale', 'phase', 'radius', 'spray')]
    result = {name: bool(np.array_equal(np.load(a/name), np.load(b/name))) for name in names}
    assert all(result.values()), result
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--out', required=True, type=Path)
    parser.add_argument('--control', required=True, type=Path)
    parser.add_argument('--max-seconds', type=int, default=240)
    a = parser.parse_args()
    assert a.max_seconds % 20 == 0 and 40 <= a.max_seconds <= 240
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip()
    if dirty:
        raise RuntimeError('Commit experiment before launch: '+dirty)
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    control = a.control.resolve()
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['PFFLIP_SPRAY_AUDIT'] = '1'
    s = make('full', escape_rule=2)
    dll = ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    job = dict(kind='Frozen C++ long-decay plunge', commit=subprocess.check_output(
        ['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),
        control=str(control), nx=NX, ny=NY, dt=DT, chunk_steps=CHUNK_STEPS,
        dump_every=EVERY, max_seconds=a.max_seconds, ports=config(s),
        initial_seed=1, emitter_seed=137, omp_threads=1, rest_criteria=CRITERIA,
        scope='Uniform 2D experimental phase=2, escape=2; no additional damping, no solver edits',
        units='Coordinates in cells; this fixture uses 1 metre per cell. Not production-ocean scale.',
        surface_proxy='Water excluding escaped: height P95 in 2-cell-wide columns; not a mesh')
    assert job['binary_sha256'] == json.loads((control/'job.json').read_text())['binary_sha256']
    write(out/'job.json', job)
    initial = dict(pos=s.pos.astype(np.float32), vel=s.vel.astype(np.float32),
                   typ=s.typ.astype(bool), spr=s.spray, esc=s.escaped, rad=s.radius)
    np.savez_compressed(out/'initial.npz', **initial)
    rows = [measure(initial, 0.)]; chunks = []; state = None
    start = time.perf_counter()
    for k in range(a.max_seconds//20):
        d = out/f'chunk_{k:02d}'
        receipt = run(d, s, CHUNK_STEPS, k*CHUNK_STEPS, state)
        if k == 0:
            same = compare_final(d, control/'bifrost')
            for j in range(0, CHUNK_STEPS, EVERY):
                for nm in ('frame', 'vel', 'ph', 'rad', 'spr', 'esc'):
                    assert np.array_equal(np.load(d/f'{nm}.{j:04d}.npy'),
                                          np.load(control/'bifrost'/f'{nm}.{j:04d}.npy'))
            write(out/'original_20s_control.json', dict(pass_control=True,
                  final_arrays=same, sampled_arrays_compared=600))
            # Check complete-state restart in the actual active spray state.
            gates = out/'resume_gate'; gates.mkdir()
            run(gates/'continuous', s, 80, CHUNK_STEPS, load_final_state(str(d)), every=40)
            run(gates/'first', s, 40, CHUNK_STEPS, load_final_state(str(d)), every=40)
            run(gates/'second', s, 40, CHUNK_STEPS+40,
                load_final_state(str(gates/'first')), every=40)
            write(gates/'result.json', dict(pass_resume=True,
                final_arrays=compare_final(gates/'continuous', gates/'second')))
        for j in range(0, CHUNK_STEPS, EVERY):
            rows.append(measure(frame(d, j), (k*CHUNK_STEPS+j+1)*DT))
        rest = rest_check(rows)
        chunks.append(dict(directory=d.name, **receipt, near_rest=rest))
        state = load_final_state(str(d))
        write(out/'measurements.json', rows)
        write(out/'progress.json', dict(chunks=chunks, elapsed_seconds=time.perf_counter()-start,
              simulated_seconds=(k+1)*CHUNK_STEPS*DT,
              status='near_rest' if rest['pass_near_rest'] else 'continuing', near_rest=rest))
        print(json.dumps(dict(chunk=k, simulation_seconds=(k+1)*CHUNK_STEPS*DT,
              elapsed_seconds=time.perf_counter()-start, near_rest=rest)), flush=True)
        if rest['pass_near_rest']:
            break
    progress = json.loads((out/'progress.json').read_text())
    progress['status'] = 'near_rest' if rest['pass_near_rest'] else 'limit_reached_not_settled'
    write(out/'progress.json', progress)
    print('COMPLETE '+progress['status'], flush=True)


if __name__ == '__main__':
    main()
