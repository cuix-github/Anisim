"""Separate droplet/carrier feedback from cell-centred density support.

Static mode only reads retained 2D states. Flight mode runs three 0.05-second
actual-graph fixtures, first reproducing the retained legacy control exactly.
No solver edits, production runs, rendering, epsilon floors or new admission
thresholds. Numerical tolerances are those of validate_interface_closure.
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
sys.path.insert(0, str(ROOT / 'gates'))
from diagnose_interface_support import state
from validate_interface_closure import airborne, cpp, py_step, compare, DT
from loop3d import load_final_state


def write(path, obj):
    path.write_text(json.dumps(obj, indent=2) + '\n', encoding='utf-8')


def support_audit():
    retained = json.loads((ROOT / 'out/runs/interface_closure_20261024/support_after/receipt.json').read_text())
    rows = []
    for expected in retained['support']:
        step = expected['step']
        s = state(step)
        active = s.masses() > 0
        density = s.number_density()
        zero = density == 0
        cells = np.argwhere(zero)
        assert cells.tolist() == expected['zero_coordinates'], 'retained support control differs'
        # Independent nearest-neighbour/kernel evaluation, not the splat stencil.
        active_points = s.pos[active]
        active_masses = s.masses()[active]
        pts = cells + .5
        distance = np.array([np.sqrt(np.min(np.sum((active_points-p)**2, axis=1))) for p in pts])
        assert np.all(distance >= 1.), 'zero support despite a point inside the kernel'
        s.p2g()
        mu, mv = s.mu.copy(), s.mv.copy()
        holes = []
        for (i, j), d in zip(cells, distance):
            faces = [((i, j+.5), mu[i,j]), ((i+1, j+.5), mu[i+1,j]),
                     ((i+.5, j), mv[i,j]), ((i+.5, j+1), mv[i,j+1])]
            checked = []
            for point, deposited in faces:
                delta = active_points - point
                weights = np.maximum(1.-np.sum(delta*delta, axis=1), 0.)**3
                independent = float(weights @ active_masses)
                assert np.isclose(deposited, independent, rtol=2e-12, atol=2e-12)
                checked.append(float(deposited))
            holes.append(dict(cell=[int(i),int(j)], nearest_active_distance=float(d),
                boundary_ring=int(min(i,j,s.nx-1-i,s.ny-1-j)),
                face_transfer_mass=checked,
                positive_faces=int(np.count_nonzero(checked)),
                phase=float(s.cell_phase[i,j])))
        # Velocity samples use staggered face fields, NOT cell-centred counts.
        # A missing P2G face is diagnostic only: projection/boundaries may still
        # give it a finite velocity. Do not equate this with an invalid final flow.
        umissing = (mu == 0).astype(float); umissing[[0,-1],:] = 0
        vmissing = (mv == 0).astype(float); vmissing[:,[0,-1]] = 0
        exposure = (s._sample_faces(umissing, np.array([0.,.5]), s.pos) +
                    s._sample_faces(vmissing, np.array([.5,0.]), s.pos))
        airborne_mask = (s.typ == 1) & (s.radius > 0)
        rows.append(dict(step=step, retained_control=True,
            active_count=int(active.sum()), auxiliary_count=int(s.spray.sum()),
            centre_holes=len(cells), holes_with_all_four_faces_supported=sum(h['positive_faces']==4 for h in holes),
            holes_with_any_face_supported=sum(h['positive_faces']>0 for h in holes),
            holes_in_boundary_cell=sum(h['boundary_ring']==0 for h in holes),
            holes_within_two_cell_boundary_band=sum(h['boundary_ring']<=1 for h in holes),
            unsupported_interior_u_faces=int(umissing.sum()),
            unsupported_interior_v_faces=int(vmissing.sum()),
            airborne_particles=int(airborne_mask.sum()),
            airborne_sampling_any_unsupported_preprojection_face=int(np.count_nonzero(exposure[airborne_mask]>0)),
            holes=holes))
    return rows


def flight_audit(out):
    rows = []
    actual = {}
    for name, mode, auxiliary in [('legacy_primary',0,False),
                                  ('auxiliary',0,True),('detached_primary',1,False)]:
        s, i = airborne(mode, dtmax=.001)
        if auxiliary:
            s.spray[i] = True
        s.p2g()
        initial_sample = s.sample_velocity(s.pos[[i]])[0].tolist()
        a, residual = cpp(s, out/name, 20)
        if name == 'legacy_primary':
            old = load_final_state(str(ROOT / 'out/runs/interface_closure_20261024/local_after/flight/m0_h0.001'))
            assert a.keys() == old.keys()
            assert all(np.array_equal(a[k], old[k]) for k in a), 'legacy control changed'
        for _ in range(20):
            py_step(s)
        row = dict(arm=name, initial_preprojection_carrier_sample=initial_sample,
            final_drop_vx=float(a['velocities'][i,0]), pressure_residual=residual,
            max_air_speed=float(np.linalg.norm(a['velocities'][a['particle_phase']==0],axis=1).max()),
            execution=compare(s,a))
        if name != 'legacy_primary':
            assert row['max_air_speed'] == 0
            assert a['mass'][i] == 0 and a['radius'][i] > 0
        actual[name] = a
        rows.append(row)
    b, c = actual['auxiliary'], actual['detached_primary']
    checks = {k:bool(np.array_equal(b[k],c[k])) for k in b if k != 'spray'}
    assert all(checks.values()), checks
    return dict(legacy_control_bit_identical=True, steps=20, simulated_seconds=20*DT,
                arms=rows, auxiliary_vs_detached_physical_arrays_identical=checks)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out',type=Path,required=True)
    parser.add_argument('--flight',action='store_true',help='also run the three tiny actual-node fixtures')
    args = parser.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'],cwd=ROOT,text=True).strip()
    out = args.out.resolve(); out.mkdir(parents=True,exist_ok=False)
    os.environ['OMP_NUM_THREADS'] = '1'
    binary = ROOT/'nodes/build/PFFlipSpike-1.0.0/lib/PFFlipSpikeOps.dll'
    receipt = dict(scope=__doc__, source_sha=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        binary_sha256=hashlib.sha256(binary.read_bytes()).hexdigest(), command=sys.argv,
        passed=False)
    start = time.perf_counter()
    try:
        receipt['support'] = support_audit()
        write(out/'receipt.json',receipt)
        if args.flight:
            receipt['flight'] = flight_audit(out)
        receipt['passed'] = True
    finally:
        receipt['seconds'] = time.perf_counter()-start
        write(out/'receipt.json',receipt)
        print(json.dumps({k:v for k,v in receipt.items() if k!='support'},indent=2))


if __name__ == '__main__':
    main()
