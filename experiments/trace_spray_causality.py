"""Causal audit of the original 20-second spray movie, not a new solver.

Keep the legacy phase law and all escape/emission/flight settings. Match public
float32 ports in an independent Python implementation, then optionally replay
the actual graph's grid velocities. Shared-stage runs are NOT independent
acceptance. Never read production caches or overwrite an existing arm.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / 'experiments'), str(ROOT / 'gates'),
               str(ROOT / 'reference2d'), str(ROOT / 'graphs')]
from spray_events_2d import NX, NY, make_sim, run_cpp, save_py
from phase_contract import PhaseContractSim
from trace_phase_sequence import read_record, delta
from loop3d import FINAL_STATE


class LegacyPortSim(PhaseContractSim):
    """Diagnostic replica of the existing fused tag-mean/float32 boundary.

    Does not replace the paper's reference or production defaults. Momentum is
    rounded independently from unrounded velocity, after deletions/reordering.
    """
    def p2g(self):
        m = self._round(self.masses())
        valid = (~self.spray).astype(float)
        faces = []
        for axis, name, shape, frame in zip(range(2), ('u', 'v'),
                (self.u.shape, self.v.shape), self._face_frames()):
            mom = (self._port_momentum[:, axis] if self._port_momentum is not None
                   else self._round(m * self.vel[:, axis]))
            mass, momentum, support, tag = self._splat(shape, frame,
                [m, mom, valid, valid * self.typ])
            setattr(self, 'm' + name, mass)
            denominator = support + float(np.float32(1e-6))
            mm = self._round(mass / denominator)
            mo = self._round(momentum / denominator)
            setattr(self, name, np.divide(mo, mm, out=np.zeros_like(mm), where=mm > 0))
            faces.append(self._round(np.clip(tag / denominator, 0, 1)))
        self._contract_faces = tuple(faces)
        self._enforce_walls()


def make_matched(tol=1e-6):
    s = make_sim(2, 'tag')
    s.__class__ = LegacyPortSim
    s.port_rounding = True
    s.pressure_tolerance = float(np.float32(tol))
    s._port_momentum = None
    s._contract_faces = None
    return s


def provenance(kind, **kw):
    dirty = subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT, text=True).strip()
    if dirty:
        raise RuntimeError('Commit intended experiment first: ' + dirty)
    dll = ROOT / 'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    return dict(kind=kind, commit=subprocess.check_output(['git', 'rev-parse', 'HEAD'],
        cwd=ROOT, text=True).strip(), binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),
        original_scope='legacy phase 0, original spray/escape ON, uniform 2D', **kw)


def verify_graph(actual, control, steps, every):
    checked = 0
    for j in range(0, steps, every):
        for channel in ('frame', 'vel', 'ph', 'rad', 'spr', 'esc'):
            name = f'{channel}.{j:04d}.npy'
            if not np.array_equal(np.load(actual / name), np.load(control / name)):
                raise AssertionError('Original movie not reproduced: ' + name)
            checked += 1
    for name in ['final.npy'] + ['final_' + n + '.npy' for n in FINAL_STATE]:
        if not np.array_equal(np.load(actual / name), np.load(control / name)):
            raise AssertionError('Original final state not reproduced: ' + name)
        checked += 1
    return dict(passed=True, bit_identical_arrays=checked)


def comparison(s, d, j):
    b = {name: np.load(d / f'{prefix}.{j:04d}.npy') for name, prefix in
         [('pos', 'frame'), ('vel', 'vel'), ('typ', 'ph'), ('spr', 'spr'),
          ('esc', 'esc'), ('rad', 'rad')]}
    b['pos'], b['vel'] = b['pos'][:, :2], b['vel'][:, :2]
    result = dict(step=j, python_count=len(s.pos), cpp_count=len(b['pos']))
    for name in ('pos', 'vel', 'typ', 'spr', 'esc', 'rad'):
        a = getattr(s, dict(spr='spray', esc='escaped', rad='radius').get(name, name))
        result[name] = delta(np.asarray(a, float), np.asarray(b[name], float)) if a.shape == b[name].shape else None
    pa = s.pos[(s.typ == 1) & ~s.spray]
    pb = b['pos'][(b['typ'] > .5) & (b['spr'] < .5)].astype(float)
    hist = [np.histogram2d(p[:, 0], p[:, 1], bins=(16, 20),
            range=((0, NX), (0, NY)))[0] for p in (pa, pb)]
    result.update(water_com_difference=float(np.linalg.norm(pa.mean(0) - pb.mean(0))),
        water_height_p99_difference=float(abs(np.quantile(pa[:, 1], .99)-np.quantile(pb[:, 1], .99))),
        water_tv=float(.5*np.abs(hist[0]/hist[0].sum()-hist[1]/hist[1].sum()).sum()))
    return result


def run_reference(out, graph, mode, steps, dt, every, tol):
    out.mkdir(parents=True, exist_ok=False)
    job = provenance('Original spray causal reference; shared stages are NOT validation',
                     graph=str(graph), mode=mode, steps=steps, dt=dt, every=every, tolerance=tol)
    (out / 'job.json').write_text(json.dumps(job, indent=2))
    s = make_matched(tol)
    p2g, project = s.p2g, s.project
    rows, events = [], []
    cpp_events = json.loads((graph / 'bifrost/ledger.json').read_text())
    first_grid = None
    first_event = None
    maxima = dict(pre=0., phase=0., post=0.)
    start = time.perf_counter()
    with open(graph / 'stages.bin', 'rb') as f:
        for j in range(steps):
            fields = read_record(f, NX, NY)
            stages = {}
            def traced_p2g():
                p2g()
                stages['pre'] = max(float(abs(s.u-fields[0]).max()), float(abs(s.v-fields[1]).max()))
                stages['phase'] = max(float(abs(s._face_phase()[0]-fields[4]).max()),
                                      float(abs(s._face_phase()[1]-fields[5]).max()))
                if mode == 'shared_both':
                    s.u, s.v = [a.copy() for a in fields[:2]]
            def traced_project(h):
                if mode == 'independent':
                    project(h)
                    stages['post'] = max(float(abs(s.u-fields[2]).max()), float(abs(s.v-fields[3]).max()))
                else:
                    s.u, s.v = [a.copy() for a in fields[2:4]]
                    stages['post'] = 0.
            s.p2g, s.project = traced_p2g, traced_project
            s.step(dt)
            for name in ('pos', 'vel', 'radius', 'scale', 'tau'):
                setattr(s, name, getattr(s, name).astype(np.float32).astype(float))
            assert np.isfinite(s.pos).all() and np.isfinite(s.vel).all()
            for name in maxima:
                maxima[name] = max(maxima[name], stages[name])
            if first_grid is None and (stages['pre'] or stages['phase']):
                first_grid = dict(step=j, **stages)
            ev = dict(s.emission_stats)
            ev.update(step=j, time=(j+1)*dt, primary=int((~s.spray).sum()))
            events.append(ev)
            keys = ('candidates', 'events', 'accepted', 'evaporated', 'reentered', 'alive')
            if first_event is None and any(ev[k] != cpp_events[j][k] for k in keys):
                first_event = dict(step=j, python=ev, cpp=cpp_events[j])
            if j % every == 0:
                rows.append(comparison(s, graph / 'bifrost', j))
                save_py(out / f'frame.{j:04d}.npz', s)
            if (j+1) % 400 == 0:
                print(f'{mode} {j+1}/{steps}, spray={s.spray.sum()}, elapsed={time.perf_counter()-start:.1f}s', flush=True)
    save_py(out / 'final.npz', s)
    result = dict(seconds=time.perf_counter()-start, first_grid_difference=first_grid,
        first_event_difference=first_event, grid_maxima=maxima, rows=rows,
        accepted=sum(e['accepted'] for e in events), final_spray=int(s.spray.sum()),
        max_water_tv=max(r['water_tv'] for r in rows),
        max_water_com_difference=max(r['water_com_difference'] for r in rows),
        max_water_height_p99_difference=max(r['water_height_p99_difference'] for r in rows))
    (out / 'ledger.json').write_text(json.dumps(events))
    (out / 'receipt.json').write_text(json.dumps(result, indent=2))
    print(json.dumps({k: v for k, v in result.items() if k != 'rows'}, indent=2), flush=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=('verify', 'graph', 'python'))
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--control', type=Path)
    p.add_argument('--graph', type=Path)
    p.add_argument('--mode', choices=('independent', 'shared_pressure', 'shared_both'), default='shared_pressure')
    p.add_argument('--steps', type=int, default=8000)
    p.add_argument('--every', type=int, default=40)
    p.add_argument('--tolerance', type=float, default=1e-6)
    a = p.parse_args(); dt = float(np.float32(.0025)); out = a.out.resolve()
    if a.action == 'verify':
        result = verify_graph(out / 'bifrost', a.control.resolve() / 'bifrost', a.steps, a.every)
        (out / 'control_receipt.json').write_text(json.dumps(result, indent=2))
        print(json.dumps(result))
    elif a.action == 'graph':
        out.mkdir(parents=True, exist_ok=False)
        (out / 'job.json').write_text(json.dumps(provenance('Legacy 20s trace; arithmetic unchanged',
            control=str(a.control.resolve()), steps=a.steps, dt=dt, every=a.every), indent=2))
        os.environ['PFFLIP_FLAT_STAGE_TRACE'] = str(out / 'stages.bin')
        run_cpp(out / 'bifrost', a.steps, dt, a.every, 2)
        del os.environ['PFFLIP_FLAT_STAGE_TRACE']
        assert (out / 'stages.bin').stat().st_size == a.steps * (12 + 3*((NX+1)*NY+NX*(NY+1))*8)
        result = verify_graph(out / 'bifrost', a.control.resolve() / 'bifrost', a.steps, a.every)
        (out / 'receipt.json').write_text(json.dumps(result, indent=2)); print(json.dumps(result))
    else:
        run_reference(out, a.graph.resolve(), a.mode, a.steps, dt, a.every, a.tolerance)


if __name__ == '__main__':
    main()
