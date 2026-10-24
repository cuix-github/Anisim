"""Read-only interface attribution and independent-trajectory comparison."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/'experiments')]
from phase_contract_2d import make


def main():
    p=argparse.ArgumentParser();p.add_argument('--state',type=Path,required=True)
    p.add_argument('--sequence',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    (out/'job.json').write_text(json.dumps(dict(kind='Read-only pressure-work and threshold attribution',
        state=str(a.state.resolve()),sequence=str(a.sequence.resolve())),indent=2))
    z=np.load(a.state/'energy_fields.npz');s=make('plunge',2,1e-10)
    s.pos=z['pos'];s.typ=z['typ'];s.vel=z['old_velocity'];dt=float(np.float32(.0025))
    pre=[z['pre_u'],z['pre_v']];post=[z['post_u'],z['post_v']];phi=[z['phi_u'],z['phi_v']]
    before=[v.copy() for v in pre];before[1][:,1:-1]+=s.g*dt
    sample=lambda fields:np.stack([s._sample_faces(v,np.array(off),s.pos) for v,off in zip(fields,((0,.5),(.5,0)))],1)
    impulse=[v-u for u,v in zip(before,post)];dv=sample(impulse)
    m=s.masses();alpha=np.where(s.typ==1,s.alpha[1],s.alpha[0])[:,None]
    vg=alpha*s.vel+(1-alpha)*sample(pre)+sample(before)-sample(pre)
    # Exact partition of kinetic work, including all cross terms equally:
    # dK = m (v_before + delta_v/2) dot delta_v; delta_v=sum(delta_v_group).
    works={};impulses={}
    for label,pred in [('pure_air_faces',lambda p:p==0),('mixed_faces',lambda p:(p>0)&(p<1)),('pure_water_faces',lambda p:p==1)]:
        impulses[label]=sample([v*pred(p) for v,p in zip(impulse,phi)])
        works[label]=m*np.sum((vg+.5*dv)*impulses[label],axis=1)
    wsum=sum(works.values());expected=.5*m*(np.sum((vg+dv)**2,axis=1)-np.sum(vg**2,axis=1))
    assert np.allclose(wsum,expected,rtol=1e-10,atol=1e-7)
    result=dict(work_reconstruction_max_error=float(abs(wsum-expected).max()),groups={})
    near_wall=(s.pos[:,0]<3)|(s.pos[:,0]>s.nx-3)
    for label,mask in [('water',s.typ==1),('air',s.typ==0),('water_phi_05_09',(s.typ==1)&(z['phi']>=.5)&(z['phi']<.9)),
        ('water_within_3_cells_of_side_wall',(s.typ==1)&near_wall),('water_away_from_side_wall',(s.typ==1)&~near_wall)]:
        result['groups'][label]=dict(count=int(mask.sum()),total_pressure_work=float(expected[mask].sum()),
            net_by_face_group={n:float(v[mask].sum()) for n,v in works.items()},
            max_impulse_by_face_group={n:float(np.linalg.norm(v[mask],axis=1).max()) for n,v in impulses.items()})
    rows=[]
    py=make('plunge',2,1e-10);cpp=make('plunge',2,1e-10)
    for path in sorted((a.sequence/'python').glob('frame.*.npz')):
        j=int(path.name.split('.')[1]);q=np.load(path)
        py.pos=q['pos'].astype(float);py.vel=q['vel'].astype(float);py.typ=q['typ'];py._port_momentum=None;py.p2g()
        cpp.pos=np.load(a.sequence/'bifrost'/f'frame.{j:04d}.npy')[:,:2].astype(float)
        cpp.vel=np.load(a.sequence/'bifrost'/f'vel.{j:04d}.npy')[:,:2].astype(float)
        cpp.typ=(np.load(a.sequence/'bifrost'/f'ph.{j:04d}.npy')>.5).astype(int);cpp._port_momentum=None;cpp.p2g()
        mismatch=(py.cell_phase>=.5)!=(cpp.cell_phase>=.5)
        row=dict(step=j,time=(j+1)*dt,position_max=float(abs(py.pos-cpp.pos).max()),
            velocity_max=float(abs(py.vel-cpp.vel).max()),cell_tags_different=int(mismatch.sum()),
            face_phi_max=max(float(abs(a-b).max()) for a,b in zip(py._face_phase(),cpp._face_phase())))
        if mismatch.any():
            ii=np.argwhere(mismatch);row['cells']=[dict(index=ij.tolist(),python=float(py.cell_phase[tuple(ij)]),
                cpp=float(cpp.cell_phase[tuple(ij)])) for ij in ii[:12]]
        rows.append(row)
    result['independent_sequence_sampled']=rows
    result['first_sampled_cell_tag_mismatch']=next((r for r in rows if r['cell_tags_different']),None)
    (out/'receipt.json').write_text(json.dumps(result,indent=2))
    print(json.dumps({k:v for k,v in result.items() if k!='independent_sequence_sampled'},indent=2))
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axs=plt.subplots(1,2,figsize=(12,5))
    sel=s.typ==1;axs[0].scatter(s.pos[~sel,0],s.pos[~sel,1],s=.15,c='#cbd5e1')
    c=axs[0].scatter(s.pos[sel,0],s.pos[sel,1],s=2,c=works['pure_air_faces'][sel],cmap='coolwarm',vmin=-15000,vmax=15000)
    axs[0].set(xlim=(0,64),ylim=(0,80),aspect='equal',title='3.2025 s: pressure work from pure-air faces\nWater particles; clipped color scale')
    fig.colorbar(c,ax=axs[0],label='Particle kinetic work (Vp=1 units)')
    axs[1].plot([r['time'] for r in rows],[r['position_max'] for r in rows],label='Max position difference (cells)')
    axs[1].plot([r['time'] for r in rows],[r['cell_tags_different'] for r in rows],label='Different cell tags (count)')
    axs[1].set_yscale('symlog',linthresh=.001)
    axs[1].set(xlabel='Time (s)',title='Independent trajectories: sampled every 0.1 s')
    axs[1].legend();fig.tight_layout();fig.savefig(out/'causality.png',dpi=160);plt.close(fig)


if __name__=='__main__':main()
