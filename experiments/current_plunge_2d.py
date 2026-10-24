"""User-requested visual replay of the cumulative uniform-2D candidate.

Not a causal ablation or the outstanding physical acceptance certificate.
Reuses the original seed and measured near-rest criteria, without added damping.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

import numpy as np
from settling_phase_2d import (ROOT, NX, NY, DT, CHUNK_STEPS, EVERY, CRITERIA,
    make, mask, config as old_config, frame, measure, rest_check, compare_final, write)
from loop3d import run_bifrost3, load_final_state
from validate_sparse_density import controls


def config(s):
    p=old_config(s)
    p.update(transport_model=3, droplet_transfer_model=1,
             volume_tau=.5, volume_model=1, droplet_evap_k=0.)
    return p


def run(out,s,steps,offset=0,state=None,every=EVERY):
    assert not out.exists(),out
    p=np.c_[s.pos,np.full(len(s.pos),.5)].astype(np.float32)
    start=time.perf_counter()
    _,_,diag=run_bifrost3('current_plunge',NX,NY,1,lambda x,y,z:mask(x,y),steps,DT,
        seed_fn=lambda *a:(p,s.typ.copy()),fused=True,extra_ports=config(s),
        init_arrays=state,step_offset=offset,final_state=True,
        dumps=('ph','vel','rad','spr','esc'),dump_every=every,
        out_dir=str(out),log_prefix=str(out/'bifcmd'))
    assert np.isfinite(diag).all() and diag[:,1].max()<=s.pressure_tolerance
    ev=[json.loads(l[len('SPRAY_EVENT '):]) for l in
        (out/'bifcmd.stderr.log').read_text().splitlines() if l.startswith('SPRAY_EVENT ')]
    assert len(ev)==steps,(len(ev),steps)
    previous=int((state['spray']>.5).sum()) if state is not None else 0
    for e in ev:
        assert e['alive']==previous+e['accepted']-e['reentered']-e['evaporated'],e
        assert e['evaporated']==0,e
        previous=e['alive']
    write(out/'ledger.json',ev)
    receipt=dict(offset=offset,steps=steps,seconds=time.perf_counter()-start,
        max_pressure_residual=float(diag[:,1].max()),
        max_particle_cfl=float(diag[:,3].max()*DT),
        births=sum(e['accepted'] for e in ev),returned=sum(e['reentered'] for e in ev),
        last_spray=previous,lifecycle_closes_every_step=True)
    write(out/'runtime.json',receipt)
    return receipt


def final_frame(d):
    a=load_final_state(str(d));typ=a['particle_phase']>.5
    return dict(pos=a['positions'][:,:2],vel=a['velocities'][:,:2],typ=typ,
        spr=a['spray']>.5,esc=typ & (a['radius']>0),rad=a['radius'])


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--max-seconds',type=int,default=240);a=p.parse_args()
    assert 40<=a.max_seconds<=240 and a.max_seconds%20==0
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    os.environ['OMP_NUM_THREADS']='1';os.environ['PFFLIP_SPRAY_AUDIT']='1'
    s=make('full',escape_rule=2)
    old=ROOT/'out/runs/settling_phase_20261023_r2'
    initial=dict(pos=s.pos.astype(np.float32),vel=s.vel.astype(np.float32),
        typ=s.typ.astype(bool),spr=s.spray,esc=s.escaped,rad=s.radius)
    with np.load(old/'initial.npz') as previous:
        assert all(np.array_equal(initial[k],previous[k]) for k in initial)
    np.savez_compressed(out/'initial.npz',**initial)
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    job=dict(kind='Cumulative candidate visual diagnostic; NOT physical certification',
        source_sha=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),
        nx=NX,ny=NY,dt=DT,chunk_steps=CHUNK_STEPS,dump_every=EVERY,
        max_seconds=a.max_seconds,ports=config(s),initial_seed=1,emitter_seed=137,
        omp_threads=1,rest_criteria=CRITERIA,initial_exactly_matches_old=True,
        old_long_run=str(old),original_20s=str(ROOT/'out/runs/spray_events_2d_20261023/plunge_20s_dt0025'),
        caveats=['Multiple documented model changes; not single-variable attribution.',
          'Detached primary contract requires evaporation OFF; shared port disables auxiliary evaporation too.',
          'Particle count/height proxy is not geometric volume certification.',
          'Uniform 2D, ST/adapt OFF. Production defaults/caches untouched.',
          'One cell = one metre in this diagnostic; marker size is not physical radius.'])
    write(out/'job.json',job)
    start=time.perf_counter();rows=[measure(initial,0.)];chunks=[];state=None
    progress=dict(status='control',chunks=chunks,simulated_seconds=0.)
    write(out/'progress.json',progress)
    try:
        (out/'control').mkdir()
        write(out/'control'/'result.json',controls(out/'control'))
        print('PASS existing off/log controls: 13 arrays each identical',flush=True)
        for k in range(a.max_seconds//20):
            d=out/f'chunk_{k:02d}'
            receipt=run(d,s,CHUNK_STEPS,k*CHUNK_STEPS,state)
            for j in range(0,CHUNK_STEPS,EVERY):
                rows.append(measure(frame(d,j),(k*CHUNK_STEPS+j+1)*DT))
            exact=measure(final_frame(d),(k+1)*CHUNK_STEPS*DT)
            rows.append(exact)
            assert all(r['water_count']==rows[0]['water_count'] for r in rows)
            state=load_final_state(str(d));rest=rest_check(rows)
            chunks.append(dict(directory=d.name,**receipt,near_rest=rest))
            write(out/'measurements.json',rows)
            progress=dict(status='continuing',chunks=chunks,
                elapsed_seconds=time.perf_counter()-start,simulated_seconds=exact['time'],
                near_rest=rest,primary_water_inventory_unchanged=True)
            write(out/'progress.json',progress)
            print(json.dumps(dict(simulation_seconds=exact['time'],runtime=receipt['seconds'],
                water_rms=exact['water_rms'],spray=exact['spray'],near_rest=rest)),flush=True)
            if k==0:
                gate=out/'resume_gate';gate.mkdir()
                run(gate/'continuous',s,80,CHUNK_STEPS,state,every=40)
                run(gate/'first',s,40,CHUNK_STEPS,state,every=40)
                run(gate/'second',s,40,CHUNK_STEPS+40,
                    load_final_state(str(gate/'first')),every=40)
                write(gate/'result.json',dict(pass_resume=True,
                    final_arrays=compare_final(gate/'continuous',gate/'second')))
            if rest['pass_near_rest']:break
        progress['status']='near_rest' if rest['pass_near_rest'] else 'limit_reached_not_settled'
    except BaseException as e:
        progress.update(status='failed',failure=repr(e));raise
    finally:
        progress['elapsed_seconds']=time.perf_counter()-start
        write(out/'progress.json',progress)
    print('COMPLETE '+progress['status'],flush=True)


if __name__=='__main__':main()
