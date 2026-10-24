"""Single common-state classification intervention; not a solver repair.

Hold particle state and the original C++ projected field fixed. Change only
one air-count sample across the measured escape threshold, to measure the
discontinuity between the fluid and ballistic-droplet motion rules.
"""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'experiments'), str(ROOT/'gates')]
from trace_spray_causality import make_matched, provenance


def main():
    p = argparse.ArgumentParser(); p.add_argument('--source', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True); a = p.parse_args()
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    report = json.loads((a.source/'receipt.json').read_text())['first_branch']
    assert len(report['skip_different']) == 1
    index = report['skip_different'][0]
    measured = report['details'][0]['values']
    count_delta = measured[1]['air'] - measured[0]['air']
    (out/'job.json').write_text(json.dumps(provenance('Common-state escape discontinuity intervention',
        source=str(a.source.resolve()), index=index, count_delta=count_delta), indent=2))
    results = []
    for flipped in (False, True):
        with np.load(a.source/'first_pre_0.npz') as z:
            s = make_matched()
            for name in ('pos','vel','typ','spray','radius','scale','tau','_port_momentum'):
                setattr(s,name,z[name].copy())
            s.xi = np.ones(len(s.pos)); s.escaped = np.zeros(len(s.pos), bool)
            s.nstep = report['step']; s.t = s.nstep*float(np.float32(.0025))
            projected = [z['post_u'].copy(),z['post_v'].copy()]
            input_count = len(s.pos)
        original_splat, original_sample = s._splat, s._sample_faces
        air_fields = []; injections = []
        def splat(shape, frame, values):
            arrays = original_splat(shape,frame,values)
            if len(values) == 2:
                air_fields.append(arrays[0])
            return arrays
        def sample(grid, offset, pts):
            values = original_sample(grid,offset,pts)
            if flipped and air_fields and grid is air_fields[0]:
                assert len(pts) == len(s.pos)
                values[index] += 2*count_delta  # air fraction averages the two face samples
                injections.append(index)
            return values
        s._splat, s._sample_faces = splat, sample
        s.project = lambda dt: (setattr(s,'u',projected[0].copy()),setattr(s,'v',projected[1].copy()))
        s.step(float(np.float32(.0025)))
        for name in ('pos','vel','radius','scale','tau'):
            setattr(s,name,getattr(s,name).astype(np.float32).astype(float))
        assert injections == ([index] if flipped else [])
        assert len(s.pos) == input_count + s.emission_stats['accepted'] - s.emission_stats['reentered'] - s.emission_stats['evaporated']
        results.append(s)
        np.savez_compressed(out/('flipped.npz' if flipped else 'control.npz'),pos=s.pos,vel=s.vel,
                            escaped=s.escaped,radius=s.radius,spray=s.spray)
    before, after = results
    assert np.array_equal(before.vel[index],np.asarray(measured[0]['post_velocity']))
    assert len(before.pos) == len(after.pos)
    changed = np.flatnonzero(np.any(before.vel != after.vel,axis=1))
    result = dict(index=index, threshold=measured[0]['threshold'],
        measured_pre_position_distance=float(np.linalg.norm(np.asarray(measured[0]['position'])-measured[1]['position'])),
        control_escaped=bool(before.escaped[index]), flipped_escaped=bool(after.escaped[index]),
        velocity_control=before.vel[index].tolist(),velocity_flipped=after.vel[index].tolist(),
        velocity_jump=float(np.linalg.norm(before.vel[index]-after.vel[index])),
        position_jump=float(np.linalg.norm(before.pos[index]-after.pos[index])),
        changed_velocity_indices=changed.tolist(),
        radius_control=float(before.radius[index]),radius_flipped=float(after.radius[index]),
        event_counts_unchanged=all(before.emission_stats[k]==after.emission_stats[k] for k in
            ('events','accepted','reentered','evaporated','alive')))
    assert changed.tolist() == [index]
    assert result['control_escaped'] and not result['flipped_escaped']
    (out/'receipt.json').write_text(json.dumps(result,indent=2)); print(json.dumps(result,indent=2))


if __name__ == '__main__':
    main()
