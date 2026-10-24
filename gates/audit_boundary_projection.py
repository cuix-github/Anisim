"""Finite separation of sampling gaps, pressure projection and density feedback.

This prescribes tiny single-phase grid channels, deliberately bypassing P2G.
It is NOT an end-to-end boundary or long-flow certificate. Five fixed arms,
one unchanged legacy control, no solver changes and no parameter search.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/p) for p in ('gates', 'reference2d', 'graphs')]
from validate_interface_closure import airborne, cpp, make, LIMITS
from loop3d import load_final_state
from validate_3d import run_node

NX, NY, DT = 16, 12, 1/256


def write(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n', encoding='utf-8')


def divergence(u, v):
    return np.diff(u, axis=0)+np.diff(v, axis=1)


def gradient(p):
    u = np.zeros((NX+1, NY)); v = np.zeros((NX, NY+1))
    u[1:-1] = np.diff(p, axis=0); v[:, 1:-1] = np.diff(p, axis=1)
    return u, v


def independent_projection(u, v):
    # Dense graph Laplacian assembled independently of the solver's stencil.
    a = np.zeros((NX*NY, NX*NY))
    for i in range(NX):
        for j in range(NY):
            c = i*NY+j
            for ni, nj in ((i+1,j), (i,j+1)):
                if ni<NX and nj<NY:
                    n = ni*NY+nj
                    a[c,c] += 1; a[n,n] += 1
                    a[c,n] -= 1; a[n,c] -= 1
    rhs = -divergence(u,v).ravel()/DT
    assert abs(rhs.sum()) < 1e-10
    p = np.zeros(NX*NY)
    p[1:] = np.linalg.solve(a[1:,1:], rhs[1:])  # fix pressure gauge only
    assert np.max(abs(a@p-rhs)) < 1e-9
    gu, gv = gradient(p.reshape(NX,NY))
    return u-DT*gu, v-DT*gv


def main():
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out = Path(sys.argv[1]).resolve(); out.mkdir(parents=True,exist_ok=False)
    os.environ['OMP_NUM_THREADS'] = '1'
    binary = ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    r = dict(passed=False,scope=__doc__,
        source_sha=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(),
        limits=dict(velocity=LIMITS['velocity'],divergence=5e-5,pressure_relative=1e-6),
        arms=[])
    write(out/'job.json',r)
    start = time.perf_counter()
    try:
        s, _ = airborne(0)
        actual, _ = cpp(s,out/'legacy_control',20)
        old = load_final_state(str(ROOT/'out/runs/interface_closure_20261024/local_after/flight/m0_h0.001'))
        assert actual.keys()==old.keys() and all(np.array_equal(actual[k],old[k]) for k in old)
        r['legacy_control_bit_identical'] = True

        # Exactly representable dyadic streamfunction: D curl(psi) = 0 and
        # psi=0 at walls. Add a gradient to require a nontrivial projection.
        x,y = np.meshgrid(np.arange(NX+1),np.arange(NY+1),indexing='ij')
        psi = x*(NX-x)*y*(NY-y)/2048
        truth_u = np.diff(psi,axis=1); truth_v = -np.diff(psi,axis=0)
        assert np.all(divergence(truth_u,truth_v)==0)
        i,j = np.meshgrid(np.arange(NX),np.arange(NY),indexing='ij')
        gu,gv = gradient((i-7)**2/16+(j-5)**2/32)
        initial_u = truth_u+DT*gu; initial_v = truth_v+DT*gv
        assert np.max(abs(divergence(initial_u,initial_v)))>0
        exact_u,exact_v = independent_projection(initial_u,initial_v)
        assert max(np.max(abs(exact_u-truth_u)),np.max(abs(exact_v-truth_v)))<1e-12
        seed = make(mode=0); seed.typ[:]=0
        regular = seed.pos.copy()
        gap = regular.copy()
        moved = np.linalg.norm(gap-[.5,6.5],axis=1)<1
        gap[moved,0] += 3  # retain every particle; remove only local sampling
        seed.pos=gap
        assert seed.number_density()[0,6]==0
        r['fixture'] = dict(nx=NX,ny=NY,dt=DT,particles=len(regular),
            relocated_particles=int(moved.sum()),centre_gap=[0,6],missing_v_face=[0,6],
            prescribed_channels_bypass_p2g=True)
        outputs = {}
        for name,pos,tau,missing_face in [
            ('regular_tau0',regular,0.,False),
            ('regular_tau05',regular,.5,False),
            ('centre_gap_tau0',gap,0.,False),
            ('centre_gap_tau05',gap,.5,False),
            ('missing_face_tau0',regular,0.,True)]:
            d=out/name; inp=d/'input'; inp.mkdir(parents=True)
            mu=np.ones_like(initial_u); mv=np.ones_like(initial_v)
            ou=initial_u.copy(); ov=initial_v.copy()
            if missing_face: mv[0,6]=0.; ov[0,6]=0.
            for axis,mass,mom in [('u',mu,ou),('v',mv,ov),
                                 ('w',np.ones((NX,NY,2)),np.zeros((NX,NY,2)))]:
                for suffix,values in [('mass',mass),('mom',mom),('phase',np.zeros_like(mass))]:
                    np.save(inp/f'{axis}_{suffix}.npy',values.astype(np.float32).ravel())
            np.save(inp/'positions.npy',np.c_[pos,np.full(len(pos),.5)].astype(np.float32))
            np.save(inp/'velocities.npy',np.zeros((len(pos),3),np.float32))
            np.save(inp/'particle_phase.npy',np.zeros(len(pos),np.float32))
            try:
                stats=run_node(str(inp),str(d/'cpp'),nx=NX,ny=NY,nz=1,dt=DT,gravity=0,
                    phase_model=2,escape=0,droplet_model=0,emit_model=0,volume_tau=tau,
                    rho0_face=3.1015625,alpha_air=0.,alpha_liquid=0.,tolerance=1e-6)
            except SystemExit as e:
                assert name=='centre_gap_tau05' and 'refuses empty/nonfinite particle support' in str(e),str(e)
                r['arms'].append(dict(name=name,expected_density_refusal=str(e)))
                continue
            assert name!='centre_gap_tau05','missing density unexpectedly accepted'
            u=np.load(d/'cpp/out_u.npy').reshape(NX+1,NY).astype(float)
            v=np.load(d/'cpp/out_v.npy').reshape(NX,NY+1).astype(float)
            eu,ev=independent_projection(ou,ov)
            error=max(np.max(abs(u-eu)),np.max(abs(v-ev)))
            truth_error=max(np.max(abs(u-truth_u)),np.max(abs(v-truth_v)))
            div=float(np.max(abs(divergence(u,v))))
            assert np.isfinite(u).all() and np.isfinite(v).all()
            assert error<LIMITS['velocity'] and div<5e-5
            assert stats['final_residual']<=1e-6
            assert np.all(u[[0,-1]]==0) and np.all(v[:,[0,-1]]==0)
            if missing_face:
                # A converged projection need not recover lost solenoidal data.
                assert truth_error>1e-3 and abs(v[0,6])>1e-3
            else: assert truth_error<LIMITS['velocity']
            outputs[name]=(u,v)
            r['arms'].append(dict(name=name,pressure_relative=stats['final_residual'],
                independent_projection_max_error=float(error),max_divergence=div,
                manufactured_flow_max_error=float(truth_error),selected_v_face=float(v[0,6]),
                wall_normal_max=0.))
        for name in ('regular_tau05','centre_gap_tau0'):
            assert all(np.array_equal(a,b) for a,b in zip(outputs[name],outputs['regular_tau0']))
        r['same_grid_with_changed_centre_sampling']=True
        r['passed']=True
    finally:
        r['seconds']=time.perf_counter()-start
        write(out/'receipt.json',r)
        print(json.dumps(r,indent=2))


if __name__=='__main__':
    main()
