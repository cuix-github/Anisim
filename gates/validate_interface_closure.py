"""Finite independent interface gates, frozen in independent-interface-closure.md.

No production inputs, rendering, or deletion. Default paths are checked before
new choices. Old coupling is a contrast, never judged against the one-way oracle.
"""
import argparse,hashlib,json,os,subprocess,sys,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(ROOT/p) for p in ('reference2d','gates','graphs','experiments')]
from volume_contract import VolumeContractSim,target_divergence
from validate_volume_contract import run as legacy_run
from validate_3d import run_node
from loop3d import run_bifrost3,load_final_state
from phase_area import area,particle_field
DT=float(np.float32(.0025))
LIMITS=dict(position=8e-6,velocity=5e-5,ode_relative=.01,rest_rms=.1,
            rest_com=.1,rest_energy_ratio=1.005,relative_area=.02)

def write(p,v):p.write_text(json.dumps(v,indent=2),encoding='utf-8')

def make(nx=16,ny=12,mode=1,tau=0.,g=0.,transport=0,dtmax=.001):
    s=VolumeContractSim(nx,ny,g=float(np.float32(g)),escape=2,
        esc_phi=float(np.float32(.3)),droplet_transfer_model=mode,
        droplet=dict(integrator=2,r_median=float(np.float32(.0001)),r_sigma=0.,
                     evap_k=0.,dt_max=float(np.float32(dtmax))),spray_cfl=.5,
        volume_tau=tau,transport_model=transport,
        alpha_l=float(np.float32(.97)),alpha_g=float(np.float32(.9)),
        pressure_tolerance=float(np.float32(1e-6)))
    s.seed(lambda x,y:y<ny/2,jitter=0,seed=1);s.calibrate()
    s.rho0_face=float(np.float32(s.rho0_face))
    s.pos=np.maximum(.51,np.minimum([nx-.51,ny-.51],s.pos)).astype(np.float32).astype(float)
    return s

def pack(s):
    n=len(s.pos);z=np.zeros(n,np.float32);one=np.ones(n,np.float32)
    m=s.masses().astype(np.float32)
    mom=(m[:,None]*s.vel).astype(np.float32) if s._port_momentum is None else s._port_momentum.astype(np.float32)
    return dict(positions=np.c_[s.pos,np.full(n,.5)].astype(np.float32),
        velocities=np.c_[s.vel,z].astype(np.float32),particle_phase=s.typ.astype(np.float32),
        mass=m,momx=mom[:,0],momy=mom[:,1],momz=z,tau=z,wt=one,wtph=s.typ.astype(np.float32),
        phase_st=s.typ.astype(np.float32),scale=s.scale.astype(np.float32),
        radius=s.radius.astype(np.float32),spray=s.spray.astype(np.float32))

def config(s):
    from spray_events_2d import ports,E
    p=ports(s.emit or E,s.droplet)
    p.update(gravity=s.g,phase_model=2,phase_alpha=1.,rho0_face=s.rho0_face,
        escape=2,esc_phi=s.esc_phi,droplet_model=2,emit_model=2 if s.emit else 0,
        droplet_transfer_model=s.droplet_transfer_model,volume_tau=s.volume_tau,volume_model=s.volume_model,
        transport_model=s.transport_model,spray_cfl=s.spray_cfl,st=0,adapt=0,
        alpha_liquid=s.alpha[1],alpha_air=s.alpha[0],tolerance=s.pressure_tolerance)
    return p

def cpp(s,out,steps,state=None,offset=0,every=1):
    assert not out.exists()
    _,_,diag=run_bifrost3('interface_closure',s.nx,s.ny,1,lambda x,y,z:y<s.ny/2,
        steps,DT,init_arrays=pack(s) if state is None else state,fused=True,
        extra_ports=config(s),final_state=True,step_offset=offset,
        dumps=('ph','vel','rad','spr','esc'),dump_every=every,
        out_dir=str(out),log_prefix=str(out/'bifcmd'))
    assert np.isfinite(diag).all() and diag[:,1].max()<=s.pressure_tolerance
    return load_final_state(str(out)),float(diag[:,1].max())

def py_step(s):
    s.step(DT)
    for n in ('pos','vel','radius','scale','tau'):
        setattr(s,n,getattr(s,n).astype(np.float32).astype(float))

def compare(s,a):
    assert len(a['positions'])==len(s.pos)
    dp=float(abs(a['positions'][:,:2]-s.pos).max());dv=float(abs(a['velocities'][:,:2]-s.vel).max())
    assert dp<LIMITS['position'] and dv<LIMITS['velocity'],(dp,dv)
    assert np.array_equal(a['particle_phase'],s.typ)
    assert np.array_equal(a['mass'],s.masses().astype(np.float32))
    assert np.array_equal(a['radius']>0,s.radius>0)
    return dict(position=dp,velocity=dv,count=len(s.pos))

def math_gate():
    n=np.array([.8,1.,1.2,1.4]);ref=np.array([1.,1.2,1.1,.9]);tau=.5
    d=target_divergence(n,ref,tau);f=np.log(n/ref)/tau
    assert abs(d.sum())<1e-14
    for scale in (.1,2.,17.):
        assert np.allclose(target_divergence(n,ref*scale,tau),d,rtol=0,atol=2e-15)
    rng=np.random.default_rng(721);base=np.sum((d-f)**2)
    for _ in range(50):
        e=rng.normal(size=4);e-=e.mean()
        # Orthogonality/Pythagoras, independently of pressure code.
        assert abs(np.sum((d+e-f)**2)-base-np.sum(e*e))<1e-13
    for bad in ([0.,1.],[np.nan,1.],[-1.,1.]):
        try:target_divergence(bad,np.ones(2),tau)
        except ValueError:pass
        else:raise AssertionError('invalid density accepted')
    return dict(target=d.tolist(),sum=float(d.sum()),projection_squared_cost=float(base),
                source_free_log_difference_ratio=float(np.exp(-1/tau)))

def prescribed(out):
    rows=[]
    for phi,r0,escaped in [(.2,0.,True),(.5,0.,False),(.5,.0001,True),(.8,.0001,False)]:
        d=out/f'phi{phi}_r{r0}';ind=d/'input';ind.mkdir(parents=True)
        for axis,shape in [('u',(17,12)),('v',(16,13)),('w',(16,12,2))]:
            for suffix,value in [('mass',1.),('mom',0.),('phase',phi)]:
                np.save(ind/f'{axis}_{suffix}.npy',np.full(shape,value,np.float32).ravel())
        for name,val in dict(positions=np.array([[8,6,.5],[12,6,.5]],np.float32),
             velocities=np.zeros((2,3),np.float32),particle_phase=np.ones(2,np.float32),
             radius_in=np.full(2,r0,np.float32)).items():np.save(ind/f'{name}.npy',val)
        stats=run_node(str(ind),str(d/'cpp'),nx=16,ny=12,nz=1,dt=DT,gravity=0,
            phase_model=2,escape=2,esc_phi=.3,rho0_face=3.1015625,droplet_model=1,
            droplet_transfer_model=1,droplet_r_sigma=0.,droplet_evap_k=0.,
            extra_arrays=('radius_in',),extra_outputs=('radius_out','out_mass'))
        m=np.load(d/'cpp/out_mass.npy');r=np.load(d/'cpp/radius_out.npy')
        assert np.all(m==(0. if escaped else 1000.)) and np.all((r>0)==escaped)
        assert stats['final_residual']==0
        rows.append(dict(phi=phi,previous_radius=r0,escaped=escaped,transfer=m.tolist()))
    return rows

def airborne(mode,dtmax=.001):
    s=make(mode=mode,dtmax=dtmax);s.typ[:]=0
    i=int(np.argmin(np.sum((s.pos-[8.25,6.25])**2,axis=1)))
    s.typ[i]=1;s.radius[i]=s.droplet['r_median'];s.vel[i,0]=3.
    return s,i

def ode_velocity(t,r):
    # Independent scalar RK4 for stationary carrier, no code from flight kernel.
    h=t/20000;v=3.
    def f(v):
        re=2*r*abs(v)/1.5e-5
        cd=24/max(re,1e-30)*(1+.15*re**.687) if re<=1000 else .44
        return -3*.001*cd*abs(v)*v/(8*r)
    for _ in range(20000):
        a=f(v);b=f(v+h*a/2);c=f(v+h*b/2);d=f(v+h*c)
        v+=h*(a+2*b+2*c+d)/6
    return v

def flight(out):
    rows=[];errors=[];steps=20
    for mode,dtmax in [(0,.001),(1,.002),(1,.001),(1,.0005)]:
        s,i=airborne(mode,dtmax);a,res=cpp(s,out/f'm{mode}_h{dtmax}',steps)
        for _ in range(steps):py_step(s)
        row=dict(mode=mode,dtmax=dtmax,**compare(s,a),residual=res,final_vx=float(a['velocities'][i,0]))
        if mode:
            assert np.all(a['velocities'][a['particle_phase']==0]==0),'one-way drop changed stationary carrier'
            assert a['mass'][i]==0 and a['radius'][i]>0
            exact=ode_velocity(steps*DT,s.droplet['r_median']);err=abs(row['final_vx']/exact-1)
            row.update(reference_vx=exact,relative_error=err);errors.append(err)
        rows.append(row)
    assert errors[0]>errors[1]>errors[2] and errors[-1]<LIMITS['ode_relative'],errors
    return rows

def restart(out):
    s,i=airborne(1);a,_=cpp(s,out/'full',80)
    b,_=cpp(s,out/'first',40);c,_=cpp(s,out/'second',40,b,40)
    checks={k:bool(np.array_equal(a[k],c[k])) for k in a}
    assert all(checks.values()),checks
    bad=pack(s);bad['mass'][i]=1000.
    try:cpp(s,out/'invalid_checkpoint',1,bad)
    except ValueError as e:assert 'inconsistent transfer' in str(e)
    else:raise AssertionError('old full-transfer state silently accepted')
    return checks

def common_states(out):
    rows=[]
    for step in (1280,3200,7200):
        d=ROOT/f'out/runs/acceptance_closure_20261024/impact_states_after/step_{step}/bifrost'
        s=make(64,80,tau=.5,transport=3,g=-9.8)
        for attr,file in [('pos','positions'),('vel','velocities'),('typ','particle_phase'),
                          ('radius','radius'),('spray','spray'),('scale','scale')]:
            a=np.load(d/f'{file}.npy')
            if attr in ('pos','vel'):a=a[:,:2]
            if attr=='spray':a=a>.5
            if attr=='typ':a=a.astype(int)
            setattr(s,attr,a.astype(float) if attr in ('pos','vel','radius','scale') else a)
        s.tau=np.zeros(len(s.pos));s.xi=np.zeros((len(s.pos),2));s.escaped=np.zeros(len(s.pos),bool)
        initial_water=int(((s.typ==1)&~s.spray).sum())
        actual,res=cpp(s,out/f'step_{step}',1)
        py_step(s);row=dict(step=step,**compare(s,actual),residual=res)
        assert int(((s.typ==1)&~s.spray).sum())==initial_water
        row['primary_water_unchanged']=initial_water;rows.append(row)
    return rows

def measure(s):
    water=(s.typ==1)&~s.spray;primary=~s.spray
    mass=np.where(s.typ==1,1000.,1.)
    return dict(water=int(water.sum()),air=int(((s.typ==0)&primary).sum()),
        rms=float(np.sqrt(np.mean(np.sum(s.vel[water]**2,axis=1)))),com=s.pos[water].mean(0).tolist(),
        energy=float(np.sum(mass[primary]*(.5*np.sum(s.vel[primary]**2,axis=1)-s.g*s.pos[primary,1]))),
        area=area(particle_field(s))['area'])

def rest(out,model=0,taus=(0.,.5)):
    rows=[]
    for tau in taus:
        s=make(48,40,tau=tau,g=-9.8,transport=3)
        s.volume_model=model
        s.seed(lambda x,y:y<12,jitter=.35,seed=1)
        s.pos=np.maximum(.51,np.minimum([47.49,39.49],s.pos)).astype(np.float32).astype(float)
        initial=measure(s);d=out/f'tau{tau}';a,res=cpp(s,d,2400,every=40)
        histories={}
        for backend in ('bifrost','python'):
            q=make(48,40,tau=tau,g=-9.8,transport=3)
            q.volume_model=model
            q.seed(lambda x,y:y<12,jitter=.35,seed=1);q.pos=s.pos.copy()
            h=[initial]
            for j in range(2400):
                if backend=='python':py_step(q)
                if j%40==0 or j==2399:
                    if backend=='bifrost':
                        end=j==2399
                        for attr,file in [('pos','positions'),('vel','velocities'),('typ','phase'),('spray','spray'),('radius','radius')]:
                            prefix={'positions':'frame','velocities':'vel','phase':'ph','spray':'spr','radius':'rad'}[file]
                            fn=('final.npy' if file=='positions' else f'final_{file}.npy') if end else f'{prefix}.{j:04d}.npy'
                            x=np.load(d/fn)
                            if attr in ('pos','vel'):x=x[:,:2].astype(float)
                            if attr=='spray':x=x>.5
                            setattr(q,attr,x)
                        q.scale=np.ones(len(q.pos))
                    h.append(measure(q))
            histories[backend]=h
            write(out/f'rest_tau{tau}_{backend}.json',h)
            assert all(x['water']==initial['water'] for x in h)
            assert h[-1]['rms']<=.1 and np.linalg.norm(np.array(h[-1]['com'])-initial['com'])<=.1
            assert max(x['energy'] for x in h)/initial['energy']<=1.005
            assert max(abs(x['area']/initial['area']-1) for x in h)<=.02
        rows.append(dict(tau=tau,volume_model=model,residual=res,history=histories))
        write(out/f'rest_tau{tau}.json',rows[-1])
    return rows

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True)
    p.add_argument('--stage',choices=('local','rest'),default='local');a=p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False);os.environ['OMP_NUM_THREADS']='1'
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    write(out/'job.json',dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),command=sys.argv,thresholds=LIMITS,scope=__doc__))
    r=dict(passed=False);start=time.perf_counter()
    try:
        if a.stage=='local':
            legacy_run(out/'legacy',0.,'density',4)
            old=ROOT/'out/runs/interface_closure_20261024/control/bifrost'
            names=['final.npy']+[p.name for p in old.glob('final_*.npy')]
            r['legacy']={n:bool(np.array_equal(np.load(old/n),np.load(out/'legacy/bifrost'/n))) for n in names}
            assert len(names)==13 and all(r['legacy'].values())
            for name,fn in [('math',lambda p:math_gate()),('lifecycle',prescribed),('flight',flight),
                            ('restart',restart),('common',common_states)]:
                r[name]=fn(out/name);write(out/'receipt.json',r);print('PASS '+name,flush=True)
        else:r['rest']=rest(out)
        r['passed']=True
    except BaseException as e:r['failure']=repr(e);raise
    finally:
        r['seconds']=time.perf_counter()-start;write(out/'receipt.json',r)
        print(json.dumps(dict(passed=r['passed'],failure=r.get('failure'),seconds=r['seconds'])),flush=True)

if __name__=='__main__':main()
