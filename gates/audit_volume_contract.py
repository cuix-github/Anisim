"""Bounded missing-volume-term witness; no solver patch or production run.

The zero-motion witness is expected to PASS as a diagnostic while the overall
paper-fidelity verdict remains NOT_ACCEPTED. A diagnostic pass is not a fix.
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
sys.path[:0] = [str(ROOT/'reference2d'), str(ROOT/'graphs')]
from phase_contract import PhaseContractSim
from loop3d import run_bifrost3


def math_contract():
    # PF Eq.8 uses +T; hence div(u_after)=-T. Positive over-density
    # must expand, not compress. tau carries seconds, independently of dt.
    density = np.array([.8, 1., 1.2])
    tau = .5
    T = np.log(1/density)/tau
    target_div = -T
    assert target_div[0] < 0 and target_div[1] == 0 and target_div[2] > 0
    assert np.allclose(np.log(1/density)/(2*tau), T/2, rtol=0, atol=0)
    nx, ny = 12, 10
    def negative_lap(n):
        diagonal = np.full(n, 2.); diagonal[[0, -1]] = 1
        return np.diag(diagonal)-np.diag(np.ones(n-1),-1)-np.diag(np.ones(n-1),1)
    A = np.kron(negative_lap(nx), np.eye(ny)) + np.kron(np.eye(nx), negative_lap(ny))
    assert np.max(abs(A@np.ones(nx*ny))) == 0
    # Equal-volume +/-20% density deviations conserve mean density, but the
    # log target does NOT have zero mean: this is the concrete compatibility trap.
    n = np.ones((nx, ny)); n[:nx//2] = .8; n[nx//2:] = 1.2
    desired = np.log(n)/tau
    assert abs(n.mean()-1) < 1e-14 and abs(desired.mean()) > .01
    dt = .01
    rhs = desired.ravel()/dt  # A=-div grad, beta=1, u*=0
    # For every p, norm(Ap-rhs)/norm(rhs) >= norm(constant part)/norm(rhs).
    lower_bound = abs(rhs.mean())*np.sqrt(len(rhs))/np.linalg.norm(rhs)
    assert lower_bound > 1e-6
    compatible = rhs-rhs.mean()
    p = np.zeros(len(rhs))
    p[1:] = np.linalg.solve(A[1:, 1:], compatible[1:])
    residual = np.linalg.norm(A@p-compatible)/np.linalg.norm(compatible)
    assert residual < 1e-11
    pp = p.reshape(nx, ny)
    u = np.zeros((nx+1, ny)); v = np.zeros((nx, ny+1))
    u[1:-1] = -dt*np.diff(pp, axis=0); v[:,1:-1] = -dt*np.diff(pp, axis=1)
    measured = np.diff(u, axis=0)+np.diff(v, axis=1)
    error = np.max(abs(measured-(desired-desired.mean())))
    assert error < 1e-11
    return dict(tau_seconds=tau, density_ratio=density.tolist(), c_div=T.tolist(),
        target_divergence=target_div.tolist(), original_target_mean=float(desired.mean()),
        incompatible_rhs_min_relative_residual=float(lower_bound),
        compatible_projection_residual=float(residual), target_divergence_max_error=float(error),
        qualification='Mean removal demonstrates compatibility only; not an adopted PF closure')


def make(deformed):
    s = PhaseContractSim(24, 24, g=0., escape=0, pressure_tolerance=float(np.float32(1e-6)))
    s.seed(lambda x,y:y<12, jitter=0, seed=1); s.calibrate()
    # Start admissible, including the same boundary locations in both arms.
    s.pos = np.clip(s.pos, .51, 23.49)
    if deformed:
        # Smooth displacement supported strictly inside the tank, preserving
        # particle count, phase, and all boundary particles. No momentum.
        x,y=s.pos.T
        sx=np.clip((x-4)/16,0,1); sy=np.clip((y-4)/16,0,1)
        s.pos[:,0] += 1.25*np.sin(2*np.pi*sx)*np.sin(np.pi*sy)**2
    s.pos=s.pos.astype(np.float32).astype(float)
    s.rho0_face=float(np.float32(s.rho0_face))
    return s


def count_field(s):
    q=s.pos-.5; b=np.floor(q+.5).astype(np.int64)
    count,=s._splat((s.nx,s.ny),(b,q-b),[np.ones(len(s.pos))])
    return count


def actual_graph(out, deformed):
    s=make(deformed); initial=s.pos.copy(); counts=count_field(s)
    p3=np.c_[initial,np.full(len(initial),.5)].astype(np.float32)
    d=out/('density_error' if deformed else 'control')
    _,_,diag=run_bifrost3('volume_witness',s.nx,s.ny,1,lambda x,y,z:y<12,
        4,float(np.float32(.0025)),seed_fn=lambda *a:(p3,s.typ.copy()),fused=True,
        final_state=True,dumps=('vel','ph'),out_dir=str(d),log_prefix=str(d/'bifcmd'),
        extra_ports=dict(gravity=0.,phase_model=2,phase_alpha=1.,rho0_face=s.rho0_face,
            escape=0,emit_model=0,droplet_model=0,transport_model=0,st=0,adapt=0))
    for _ in range(4):
        s.step(float(np.float32(.0025)))
        s.pos=s.pos.astype(np.float32).astype(float)
        s.vel=s.vel.astype(np.float32).astype(float)
    pos=np.load(d/'final.npy')[:,:2]; vel=np.load(d/'final_velocities.npy')[:,:2]
    assert np.array_equal(pos,initial) and np.array_equal(s.pos,initial)
    assert np.max(abs(vel))==0 and np.max(abs(s.vel))==0
    assert np.array_equal(count_field(s),counts)
    assert np.array_equal(pos,s.pos) and np.array_equal(vel,s.vel)
    assert np.isfinite(diag).all() and np.max(diag[:,1])==0
    return dict(particles=len(pos),steps=4,python_cpp_positions_bit_identical=True,
        max_velocity=float(np.max(abs(vel))),pressure_max_residual=float(diag[:,1].max()),
        unchanged_number_density=True),counts


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out=a.out.resolve();out.mkdir(parents=True,exist_ok=False)
    os.environ['OMP_NUM_THREADS']='1'
    dll=ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    job=dict(commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256(dll.read_bytes()).hexdigest(),scope=__doc__,
        command=sys.argv,expected_outcome='missing correction witnessed, not repaired')
    (out/'job.json').write_text(json.dumps(job,indent=2))
    report=dict(diagnostic_checks_passed=False,paper_fidelity_accepted=False)
    try:
        report['math']=math_contract()
        report['control'],base=actual_graph(out,False)
        report['density_error'],shifted=actual_graph(out,True)
        relative=(shifted[4:-4,4:-4]-base[4:-4,4:-4])/base[4:-4,4:-4]
        assert np.max(relative)>.1 and np.min(relative)<-.1
        report['density_error'].update(relative_min=float(relative.min()),relative_max=float(relative.max()),
            relative_rms=float(np.sqrt(np.mean(relative**2))),
            interpretation='Existing nonzero count-density error persists at zero velocity in both implementations')
        np.savez_compressed(out/'count_fields.npz',control=base,deformed=shifted)
        report['diagnostic_checks_passed']=True
    finally:
        (out/'receipt.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
