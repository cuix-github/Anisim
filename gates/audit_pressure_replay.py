"""Audit pressure from saved 2D states; never advance a long simulation.

First capture --stage baseline with the unmodified binary. After adding the
opt-in double-precision dump, --stage audit compares actual Bifrost solves to an
independent NumPy face-flux operator / Jacobi PCG. Each arm is one main step.
This deliberately shares coefficients to isolate the linear solve; it does not
claim that the default Python and graph phase laws agree.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT/'gates'), str(ROOT/'experiments')]
from spray_events_2d import make_sim, load_frame, D, E
from validate_spray_events import ports
from validate_droplet_node import ESCAPE
from validate_level_set_node import flat_inputs
from validate_3d import run_node

STEPS = (0, 1280, 4000, 5600, 7360, 7960)


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def prepare(source, out, step):
    state = load_frame(source, 'bifrost', step)
    s = make_sim()
    for key, name in [('pos','pos'),('vel','vel'),('typ','typ'),('spr','spray'),
                      ('esc','escaped'),('rad','radius')]:
        setattr(s, name, state[key].copy())
    s.pos, s.vel = s.pos.astype(float), s.vel.astype(float)
    s.scale, s.tau, s.xi = np.ones(len(s.pos)), np.zeros(len(s.pos)), np.ones(len(s.pos))
    m, valid = s.masses(), (~s.spray).astype(float)
    channels = {}
    for axis, name, shape, frame in zip(range(2), ('u','v'), (s.u.shape,s.v.shape), s._face_frames()):
        sums = s._splat(shape, frame, [m, m*s.vel[:,axis], valid*s.typ, valid])
        for label, arr in zip(('mass','mom','phase'), sums[:3]):
            channels[f'{name}_{label}'] = (arr/(sums[3]+float(np.float32(1e-6)))).astype(np.float32)
    arrays = flat_inputs(s, channels)
    arrays.update(radius_in=s.radius.astype(np.float32), spray_in=s.spray.astype(np.float32))
    dest = out/f'state_{step}'/'input'
    dest.mkdir(parents=True, exist_ok=True)
    for name, arr in arrays.items():
        path = dest/(name+'.npy')
        if path.exists():
            assert np.array_equal(np.load(path), arr), path
        else:
            np.save(path, arr)
    return s, arrays, dest


def call(source, out, step, label, preconditioner, tolerance, dump, max_iterations=4000):
    s, arrays, indir = prepare(source, out, step)
    d = out/f'state_{step}'/label
    if d.exists():
        raise FileExistsError(f'Preserve prior arm: {d}')
    d.mkdir()
    dt = json.loads((source/'job.json').read_text())['dt']
    if dump:
        os.environ['PFFLIP_PRESSURE_DUMP'] = str(d/'pressure_double.bin')
    else:
        os.environ.pop('PFFLIP_PRESSURE_DUMP', None)
    os.environ['PFFLIP_SPRAY_AUDIT'] = '1'
    kw = dict(ESCAPE, **ports(E,D), nx=s.nx, ny=s.ny, nz=1, dt=dt,
              gravity=float(np.float32(-9.8)), alpha_liquid=float(np.float32(.97)),
              alpha_air=float(np.float32(.9)), rho_liquid=1000, rho_air=1,
              rho0_face=s.rho0_face, step_index=step+1, pressure_solver=0,
              preconditioner=preconditioner, tolerance=tolerance, max_iterations=max_iterations)
    stats = run_node(str(indir), str(d), extra_arrays=('radius_in','spray_in'),
                     extra_outputs=('radius_out','spray_out'), **kw)
    (d/'stats.json').write_text(json.dumps(stats, indent=2))
    (d/'ports.json').write_text(json.dumps(kw, indent=2))
    return s, arrays, d, stats, dt


def divergence(faces):
    return sum(np.diff(f, axis=axis) for axis, f in enumerate(faces))


def flux(p, coeff):
    fields = [np.zeros_like(b) for b in coeff]
    for axis, f in enumerate(fields):
        sl = [slice(None)]*3
        sl[axis] = slice(1,-1)
        f[tuple(sl)] = coeff[axis][tuple(sl)]*np.diff(p, axis=axis)
    return fields


def operator(p, coeff):
    return -divergence(flux(p, coeff))


def independent_pcg(rhs, coeff, tol=1e-11):
    """Conventional Jacobi PCG; check actual b-Ap before accepting convergence."""
    diagonal = np.zeros_like(rhs)
    for axis, b in enumerate(coeff):
        sl0, sl1 = [slice(None)]*3, [slice(None)]*3
        sl0[axis], sl1[axis] = slice(None,-1), slice(1,None)
        diagonal += b[tuple(sl0)] + b[tuple(sl1)]
    p, r = np.zeros_like(rhs), rhs.copy()
    norm = max(float(np.linalg.norm(rhs)), 1e-300)
    z = r/diagonal
    d, rz = z.copy(), float(np.sum(r*z))
    for k in range(4000):
        q = operator(d, coeff)
        dq = float(np.sum(d*q))
        assert dq > 0 and np.isfinite(dq), (k,dq)
        alpha = rz/dq
        p += alpha*d
        r -= alpha*q
        if np.linalg.norm(r)/norm < tol:
            real = rhs-operator(p,coeff)
            if np.linalg.norm(real)/norm < tol*2:
                return p, k+1
            r = real
            z = r/diagonal
            d, rz = z.copy(), float(np.sum(r*z))
            continue
        z = r/diagonal
        rz_new = float(np.sum(r*z))
        d = z + (rz_new/rz)*d
        rz = rz_new
    raise AssertionError('independent PCG failed true-residual convergence')


def read_dump(path):
    with path.open('rb') as f:
        dims = tuple(int(x) for x in np.fromfile(f, dtype='<i4', count=3))
        shapes = [tuple(n+(axis==i) for i,n in enumerate(dims)) for axis in range(3)]
        def arr(shape):
            return np.fromfile(f, dtype='<f8', count=int(np.prod(shape))).reshape(shape)
        coeff = [arr(sh) for sh in shapes]
        rhs, p = arr(dims), arr(dims)
        before = [arr(sh) for sh in shapes]
        assert f.read() == b''
    return coeff, rhs, p, before


def inspect(arrays, d, stats, dt, suffix='audit'):
    coeff, rhs, p, before = read_dump(d/'pressure_double.bin')
    after = [v-dt*f for v,f in zip(before,flux(p,coeff))]
    r = rhs-operator(p,coeff)
    rnorm = float(np.linalg.norm(r)/np.linalg.norm(rhs))
    pref, iterations = independent_pcg(rhs,coeff)
    ref_vel = [v-dt*f for v,f in zip(before,flux(pref,coeff))]
    rng = np.random.default_rng(137)
    x,y = rng.normal(size=p.shape), rng.normal(size=p.shape)
    ax,ay = operator(x,coeff),operator(y,coeff)
    sym = abs(np.sum(x*ay)-np.sum(ax*y))/max(abs(np.sum(x*ay)),abs(np.sum(ax*y)),1)
    def energy(fields):
        return sum(float(np.sum(v[b>0]**2/b[b>0])) for b,v in zip(coeff,fields))/2
    bnorm = np.linalg.norm(rhs)
    coeff_err, velocity_err = [], []
    for axis,name in enumerate(('u','v')):
        shape = before[axis].shape
        phase = arrays[name+'_phase'].reshape(shape).astype(float)
        expected = 1/(1+999*np.clip(phase,0,1))
        sl = [slice(None)]*3
        sl[axis] = [0,-1]
        expected[tuple(sl)] = 0
        coeff_err.append(float(np.max(abs(expected-coeff[axis]))))
        mass = arrays[name+'_mass'].reshape(shape).astype(float)
        mom = arrays[name+'_mom'].reshape(shape).astype(float)
        v = np.divide(mom, mass, out=np.zeros_like(mass), where=mass>0)
        if axis==1:
            # step_3d promotes dt_in to double before multiplying by gravity.
            v += float(np.float32(-9.8))*float(np.float32(dt))
        v[tuple(sl)] = 0
        velocity_err.append(float(np.max(abs(v-before[axis]))))
    exported = [np.load(d/f'out_{n}.npy').reshape(v.shape) for n,v in zip(('u','v','w'),after)]
    result = dict(iterations=stats['iterations_used'], reported_residual=stats['final_residual'],
        true_residual=rnorm, independent_iterations=iterations,
        rhs_assembly_relative=float(np.linalg.norm(rhs+divergence(before)/dt)/bnorm),
        compatibility_relative=float(abs(rhs.sum())/(np.sqrt(rhs.size)*bnorm)),
        coefficient_max_error=max(coeff_err), pre_velocity_max_error=max(velocity_err),
        symmetry_relative=float(sym), constant_nullspace_max=float(abs(operator(np.ones_like(p),coeff)).max()),
        positive_energy=float(np.sum(x*ax)),
        divergence_max=float(abs(divergence(after)).max()),
        divergence_residual_identity_max=float(abs(divergence(after)+dt*r).max()),
        velocity_max_error_vs_independent=max(float(abs(a-b).max()) for a,b in zip(after,ref_vel)),
        exported_velocity_roundtrip_exact=all(np.array_equal(a.astype(np.float32),b) for a,b in zip(after,exported)),
        weighted_energy_before=energy(before), weighted_energy_after=energy(after))
    result['passed'] = bool(rnorm < float(json.loads((d/'ports.json').read_text())['tolerance'])*1.02
        and result['coefficient_max_error']==0 and result['pre_velocity_max_error']<1e-12
        and result['rhs_assembly_relative']<1e-14 and result['compatibility_relative']<1e-14
        and sym<1e-12 and result['positive_energy']>0 and result['constant_nullspace_max']==0
        and result['divergence_residual_identity_max']<1e-10
        and result['velocity_max_error_vs_independent']<2e-4
        and result['weighted_energy_after']<=result['weighted_energy_before']*(1+1e-10)
        and result['exported_velocity_roundtrip_exact'])
    (d/(suffix+'.json')).write_text(json.dumps(result,indent=2))
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True)
    ap.add_argument('--out',type=Path,required=True);ap.add_argument('--stage',choices=('baseline','audit','recheck'),required=True)
    a=ap.parse_args();source,out=a.source.resolve(),a.out.resolve();out.mkdir(parents=True,exist_ok=True)
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    job=dict(kind='Existing 2D states: one-step pressure audits', source=str(source), command=sys.argv,
             commit=subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=ROOT).strip(),
             binary_sha256=digest(dll), production_cache='untouched', stage=a.stage)
    job_path=out/f'job_{a.stage}.json'
    assert not job_path.exists(), job_path
    job_path.write_text(json.dumps(job,indent=2))
    if a.stage=='baseline':
        call(source,out,5600,'control_unmodified',1,1e-6,False)
        return
    receipt={}
    for step in STEPS:
        rows={}
        for name,pre,tol in [('mg',1,1e-6),('jacobi',0,1e-6),('mg_tight',1,1e-10)]:
            if a.stage=='recheck':
                _,arrays,_=prepare(source,out,step)
                d=out/f'state_{step}'/name
                stats=json.loads((d/'stats.json').read_text())
                dt=json.loads((d/'ports.json').read_text())['dt']
            else:
                _,arrays,d,stats,dt=call(source,out,step,name,pre,tol,True)
            rows[name]=inspect(arrays,d,stats,dt,'audit_rechecked' if a.stage=='recheck' else 'audit')
            print(step,name,json.dumps(rows[name]),flush=True)
        if step==5600:
            control=out/f'state_{step}'/'control_unmodified'
            files=list(control.glob('*.npy'))
            assert files and {f.name for f in files} == {f.name for f in (out/f'state_{step}'/'mg').glob('*.npy')}
            rows['dump_does_not_change_outputs']=all(np.array_equal(np.load(f),np.load(out/f'state_{step}'/'mg'/f.name)) for f in files)
        receipt[str(step)]=rows
        (out/('receipt_rechecked.json' if a.stage=='recheck' else 'receipt.json')).write_text(json.dumps(receipt,indent=2))
    assert all(r['passed'] for rows in receipt.values() for r in rows.values() if isinstance(r,dict))
    assert receipt['5600']['dump_does_not_change_outputs']


if __name__=='__main__':
    main()
