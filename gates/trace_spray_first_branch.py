"""First discrete branch in the original spray scene, with independent pressure.

One Python arm replays the bit-verified original C++ grid. The other computes
its own pressure with matched phase and port semantics. No physical fix here.
"""
import argparse
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'experiments'), str(ROOT/'gates')]
from trace_spray_causality import make_matched, provenance, read_record


def snapshot(s):
    return {k: getattr(s, k).copy() for k in
            ('pos', 'vel', 'typ', 'spray', 'radius', 'scale', 'tau', '_port_momentum')
            if getattr(s, k, None) is not None}


def fractions(s, pre, index):
    # Rebuild counts at the pre-step locations, using unchanged kernel math.
    from copy import copy
    q = copy(s)
    for k, v in pre.items():
        setattr(q, k, v)
    fu, fv = q._face_frames()
    au, lu = q._splat(q.u.shape, fu, [(q.typ == 0).astype(float),
                      ((q.typ == 1) & ~q.spray).astype(float)])
    av, lv = q._splat(q.v.shape, fv, [(q.typ == 0).astype(float),
                      ((q.typ == 1) & ~q.spray).astype(float)])
    pos = q.pos[[index]]
    sample = lambda a,b: float(.5*(q._sample_faces(a, np.array([0.,.5]), pos)[0]+
                                   q._sample_faces(b, np.array([.5,0.]), pos)[0]))
    return dict(air=sample(au,av), liquid=sample(lu,lv),
                threshold=(1-q.esc_phi)*q.rho0_face, position=pos[0].tolist(),
                velocity=q.vel[index].tolist(), phase=int(q.typ[index]))


def main():
    p = argparse.ArgumentParser(); p.add_argument('--out', type=Path, required=True)
    p.add_argument('--graph', type=Path, required=True); p.add_argument('--steps', type=int, default=2400)
    a = p.parse_args(); out = a.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    (out/'job.json').write_text(json.dumps(provenance('First branch in original full-spray scene',
        graph=str(a.graph.resolve()), steps=a.steps), indent=2))
    sims = [make_matched(), make_matched()]
    dt = float(np.float32(.0025)); start = time.perf_counter()
    history = []
    with open(a.graph/'stages.bin', 'rb') as f:
        for j in range(a.steps):
            fields = read_record(f,64,80)
            pre = [snapshot(s) for s in sims]
            stages = []
            for arm, s in enumerate(sims):
                data = {}
                old_g2p = s.g2p
                old_emit = s._emit_clouds
                def g2p(u,v,skip=None):
                    data['skip'] = skip.copy()
                    return old_g2p(u,v,skip)
                def emit(h, idx, cells, gp):
                    data['candidates'] = idx.copy()
                    data['gp'] = gp.copy()
                    data['lf'] = s._audit_liquid_fraction.copy()
                    return old_emit(h,idx,cells,gp)
                s.g2p, s._emit_clouds = g2p, emit
                if arm == 0:
                    s.project = lambda h: (setattr(s,'u',fields[2].copy()), setattr(s,'v',fields[3].copy()))
                s.step(dt)
                s.g2p, s._emit_clouds = old_g2p, old_emit
                for name in ('pos','vel','radius','scale','tau'):
                    setattr(s,name,getattr(s,name).astype(np.float32).astype(float))
                data['post_u'], data['post_v'] = s.u.copy(), s.v.copy()
                stages.append(data)
            assert pre[0]['pos'].shape == pre[1]['pos'].shape, 'Prior topology mismatch should have stopped trace'
            skip_diff = np.flatnonzero(stages[0]['skip'] != stages[1]['skip'])
            candidate_diff = np.setxor1d(stages[0]['candidates'],stages[1]['candidates'])
            events_equal = all(sims[0].emission_stats[k] == sims[1].emission_stats[k] for k in
                               ('events','accepted','evaporated','reentered','alive'))
            row = dict(step=j,time=(j+1)*dt,
                pre_position_max=float(abs(pre[0]['pos']-pre[1]['pos']).max()),
                pre_velocity_max=float(abs(pre[0]['vel']-pre[1]['vel']).max()),
                projected_velocity_max=max(float(abs(stages[0]['post_u']-stages[1]['post_u']).max()),
                    float(abs(stages[0]['post_v']-stages[1]['post_v']).max())),
                skip_different=skip_diff.tolist(), candidate_different=candidate_diff.tolist(),
                events_equal=events_equal)
            if j%40 == 0:
                history.append(row)
            if len(skip_diff) or len(candidate_diff) or not events_equal:
                inds = np.union1d(skip_diff,candidate_diff)
                row['details'] = []
                for i in inds[:20]:
                    row['details'].append(dict(index=int(i), values=[dict(
                        **fractions(s,pre[k],i), escaped=bool(stages[k]['skip'][i]),
                        gradient=float(stages[k]['gp'][i]), liquid_fraction=float(stages[k]['lf'][i]),
                        post_velocity=s.vel[i].tolist() if len(s.vel)>i else None)
                        for k,s in enumerate(sims)]))
                row['emission'] = [s.emission_stats for s in sims]
                for k in range(2):
                    np.savez_compressed(out/f'first_pre_{k}.npz',**pre[k],
                        post_u=stages[k]['post_u'],post_v=stages[k]['post_v'],
                        cpp_pre_u=fields[0],cpp_pre_v=fields[1])
                result = dict(first_branch=row,history=history,seconds=time.perf_counter()-start)
                (out/'receipt.json').write_text(json.dumps(result,indent=2))
                print(json.dumps(row,indent=2),flush=True)
                return
            if (j+1)%400==0:
                print(f'{j+1} steps, no discrete difference yet',flush=True)
    raise RuntimeError('No branch found in requested range; do not infer its absence later')


if __name__ == '__main__':
    main()
