"""Read-only analysis of saved pressure replays and the original 20 s caches."""
import argparse
import json
from pathlib import Path
import sys

import numpy as np

from audit_pressure_replay import (STEPS, prepare, read_dump, flux,
                                   independent_pcg, E, D)
from spray_events_2d import make_sim, load_frame
import spray_events


def source_audit(source, out, step):
    s, _, _ = prepare(source, out, step)
    d = out/f'state_{step}'/'mg'
    dt = json.loads((d/'ports.json').read_text())['dt']
    coeff, rhs, p, before = read_dump(d/'pressure_double.bin')
    after = [v-dt*f for v,f in zip(before,flux(p,coeff))]
    s.u, s.v = after[0][:,:,0], after[1][:,:,0]
    air, water = (s.typ==0).astype(float), ((s.typ==1)&~s.spray).astype(float)
    face = [s._splat(g.shape, f, [air,water]) for g,f in zip((s.u,s.v),s._face_frames())]
    counts = [sum(s._sample_faces(face[i][ph],np.array(offset),s.pos)
                  for i,offset in enumerate(((0,.5),(.5,0))))*.5/s.rho0_face for ph in range(2)]
    fa,fl = counts
    esc = np.where(s.typ==1,fa,fl) > 1-s.esc_phi
    mid = s.pos + .5*dt*s.sample_velocity(s.pos)
    moved = s.pos + dt*s.sample_velocity(np.maximum(mid,.51))
    moved = np.clip(moved,.51,[s.nx-.51,s.ny-.51])
    cell = np.clip(np.floor(moved).astype(int),0,[s.nx-1,s.ny-1])
    grad = s._grad_norm()[cell[:,0],cell[:,1]]
    cand = (s.typ==1)&~s.spray&~esc&(fl<1-s.esc_phi)&(grad>E['grad_threshold'])
    idx = np.flatnonzero(cand)
    events = spray_events.requested_events(E['rate']*(grad[idx]-E['grad_threshold'])*dt,
                                          spray_events.keys(E['seed'],step+1,idx))
    near = np.min(np.c_[moved[idx], np.array([s.nx,s.ny])-moved[idx]],axis=1)<1.5
    noair = fa[idx] < 1e-12
    accept, reject, born_noair, born_near, births = 0,0,0,0,[]
    # These samples have no budget rejection in the actual graph. Check that
    # assumption before this simplified geometry-only reconstruction is used.
    actual=json.loads((d/'stats.json').read_text())['spray_events'][0]
    assert actual['rejected_cap']==0
    for j,(pi,ne) in enumerate(zip(idx,events)):
        for ev in range(int(ne)):
            _,dp,_,_=spray_events.cloud(E['seed'],step+1,int(pi),ev,int(E['cloud_size']),
                                       D['r_sigma'],E['event_radius'],D['r_min'],E['radius_max'])
            pos=moved[pi]+E['spread']*dp
            keep=((pos>=.51)&(pos<=np.array([s.nx,s.ny])-.51)).all(1)
            n=int(keep.sum());accept+=n;reject+=int((~keep).sum())
            born_noair+=n*int(noair[j]);born_near+=n*int(near[j]);births.extend(pos[keep])
    reconstructed=dict(candidates=len(idx),events=int(events.sum()),accepted=accept,rejected_geometry=reject)
    assert all(reconstructed[k]==actual[k] for k in reconstructed), (step,reconstructed,actual)
    # Isolate the phase law, keeping exactly the same incoming velocity and RHS.
    s.p2g()
    phase=s._face_phase()
    bmass=[(1/(1+999*phi))[:,:,None] for phi in phase]+[np.zeros_like(coeff[2])]
    bmass[0][[0,-1],:]=0;bmass[1][:,[0,-1]]=0
    pmass,_=independent_pcg(rhs,bmass)
    mass_after=[v-dt*f for v,f in zip(before,flux(pmass,bmass))]
    phase_speed_delta=max(float(abs(a-b).max()) for a,b in zip(after,mass_after))
    rec=dict(time=(step+1)*dt,matched_actual_graph_counts=reconstructed,
             candidates_without_sampled_air=int(noair.sum()),candidates_near_wall=int(near.sum()),
             noair_also_near_wall=int((noair&near).sum()),
             births_without_sampled_air=born_noair,births_near_wall=born_near,
             candidate_grad_median=float(np.median(grad[idx])) if len(idx) else 0,
             candidate_grad_max=float(grad[idx].max()) if len(idx) else 0,
             phase_law_only_velocity_max_difference=phase_speed_delta)
    np.savez_compressed(d/'source_diagnosis.npz',positions=moved[idx],noair=noair,
                        near_wall=near,births=np.asarray(births).reshape(-1,2),
                        grad=grad[idx],air_fraction=fa[idx],water_fraction=fl[idx])
    return rec


def history(source):
    d=source/'bifrost'
    diag=np.array([np.load(f).ravel() for f in sorted(d.glob('diag.*.npy'))])
    ledger=json.loads((d/'ledger.json').read_text())
    dt=json.loads((source/'job.json').read_text())['dt']
    rows=[]
    for step in range(0,len(ledger),40):
        a=load_frame(source,'bifrost',step)
        row=dict(time=(step+1)*dt)
        for name,mask in [('air',a['typ']==0),('water',(a['typ']==1)&~a['spr']),('spray',a['spr'])]:
            v=a['vel'][mask].astype(float)
            row[name]=dict(count=int(mask.sum()),max_component=float(abs(v).max()) if len(v) else 0,
                           rms=float(np.sqrt(np.mean(np.sum(v*v,axis=1)))) if len(v) else 0)
        rows.append(row)
    windows=[]
    for start in range(0,len(diag),800):
        sub=ledger[start:start+800]
        windows.append(dict(time_start=start*dt,
                            max_iterations=int(diag[start:start+800,0].max()),
                            max_residual=float(diag[start:start+800,1].max()),
                            max_divergence=float(diag[start:start+800,2].max()),
                            births=sum(r['accepted'] for r in sub),
                            candidates=sum(r['candidates'] for r in sub)))
    return dict(finite=bool(np.isfinite(diag).all()),steps=len(diag),dt=dt,
                max_iterations=int(diag[:,0].max()),max_residual=float(diag[:,1].max()),
                max_divergence=float(diag[:,2].max()),windows=windows,particle_samples=rows),diag


def figure(out, h, diag):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(3,1,figsize=(11,8),sharex=True,layout='constrained')
    t=(np.arange(len(diag))+1)*h['dt']
    ax[0].semilogy(t,diag[:,1],lw=.45,color='#466a9f',label='reported pressure residual')
    ax[0].axhline(1e-6,color='#ad4459',ls='--',label='tolerance')
    ax[0].set_ylabel('relative residual');ax[0].legend(loc='lower right')
    samples=h['particle_samples'];st=[r['time'] for r in samples]
    for name,color in [('air','#738499'),('water','#2472ba'),('spray','#bc7b1f')]:
        ax[1].plot(st,[r[name]['rms'] for r in samples],label=name,color=color)
    ax[1].set_ylabel('particle RMS speed\n(cells/s)');ax[1].legend(ncol=3)
    w=h['windows']
    ax[2].bar([r['time_start'] for r in w],[r['births'] for r in w],width=1.85,align='edge',color='#bc7b1f')
    ax[2].set_ylabel('spray births per 2 s');ax[2].set_xlabel('simulation time (s)')
    fig.suptitle('Existing Bifrost 2D cache: late spray increases without growing pressure residual\nThis does not establish stability of the full fluid / classification pipeline',fontsize=12)
    fig.savefig(out/'diagnostics.png',dpi=160)
    plt.close(fig)


def main():
    p=argparse.ArgumentParser();p.add_argument('--source',required=True,type=Path)
    p.add_argument('--out',required=True,type=Path);a=p.parse_args()
    source,out=a.source.resolve(),a.out.resolve()
    rec={str(step):source_audit(source,out,step) for step in STEPS}
    print(json.dumps(rec,indent=2))
    h,diag=history(source)
    result=dict(source_audits=rec,history=h)
    (out/'flow_diagnosis.json').write_text(json.dumps(result,indent=2))
    figure(out,h,diag)


if __name__=='__main__':main()
