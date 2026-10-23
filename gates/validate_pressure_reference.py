"""Independent 2D failure checks and adaptive-reference before/after controls."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'reference2d'))
from pfflip2d import Sim
from adaptive_solve import build_cells, Operator, Multigrid


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--out', type=Path, required=True)
    p.add_argument('--cpp', type=Path, required=True)
    p.add_argument('--before', type=Path)
    a = p.parse_args()
    assert not subprocess.check_output(['git','status','--porcelain'], cwd=ROOT, text=True).strip()
    out = a.out.resolve(); out.mkdir(parents=True, exist_ok=False)
    (out/'job.json').write_text(json.dumps(dict(commit=subprocess.check_output(
        ['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(), cpp=str(a.cpp.resolve())), indent=2))
    source = a.cpp / 'flow2d_jacobi'; inp = source/'input'
    def sim():
        s = Sim(24,20)
        phases = []
        for c, shape in [('u',(25,20)), ('v',(24,21))]:
            mass = np.load(inp/(c+'_mass.npy')).reshape(shape).astype(float)
            mom = np.load(inp/(c+'_mom.npy')).reshape(shape).astype(float)
            setattr(s,c,mom/mass)
            phases.append(np.load(inp/(c+'_phase.npy')).reshape(shape).astype(float))
        s._face_phase = lambda: tuple(phases)
        s.v += float(np.float32(-9.8))*float(np.float32(.01))
        s._enforce_walls()
        return s
    dt = float(np.float32(.01)); tol = float(np.float32(1e-6))
    s = sim(); s.project(dt,tol,4000)
    errors = {c:float(abs(getattr(s,c)-np.load(source/'converged'/('out_'+c+'.npy')).reshape(getattr(s,c).shape)).max()) for c in ('u','v')}
    assert max(errors.values()) < 5e-6, errors
    s = sim(); old = [s.u.copy(), s.v.copy()]
    try: s.project(dt,tol,1)
    except RuntimeError: pass
    else: raise AssertionError('uniform reference accepted capped pressure')
    assert all(np.array_equal(x,y) for x,y in zip(old,[s.u,s.v]))
    uniform_failure_residual = s.pressure_residual
    z = Sim(12,12); z._face_phase = lambda:(np.ones_like(z.u), np.ones_like(z.v))
    z.project(dt,tol,0); assert z.pressure_residual == 0
    rows = dict(uniform_cpp_max_face_error=errors, uniform_capped_residual=uniform_failure_residual,
                uniform_capped_faces_unchanged=True, uniform_zero_budget_zero_rhs=True)
    i,j = np.indices((24,24)); target = np.where((i-12)**2+(j-12)**2<36,0,2)
    cells = build_cells(24,24,2,target)
    beta = lambda x,y,axis,span: 1/(1+999*(y<11))
    for name,solver in [('jacobi',Operator(cells,beta)),('mg',Multigrid(cells,beta))]:
        op = solver if name=='jacobi' else solver.ops[0]
        rng = np.random.default_rng(37); b = rng.normal(size=len(cells)); b -= b.mean()
        x,it = solver.solve(b,tol=1e-8,max_iter=4000)
        residual = float(np.linalg.norm(b-op.matvec(x))/np.linalg.norm(b))
        assert residual < 1e-8
        np.save(out/(name+'.npy'),x)
        if a.before: assert np.array_equal(x,np.load(a.before/(name+'.npy')))
        refused = False
        try: solver.solve(b,tol=1e-8,max_iter=1)
        except RuntimeError: refused = True
        assert refused == bool(a.before), (name,refused)
        zx,_ = solver.solve(np.zeros_like(b),tol=1e-8,max_iter=0)
        assert np.array_equal(zx,np.zeros_like(b))
        rows[name] = dict(residual=residual, iterations=it, capped_refused=refused,
                          bit_identical=bool(a.before), zero_budget_zero_rhs=True)
    (out/'receipt.json').write_text(json.dumps(rows,indent=2)); print(json.dumps(rows,indent=2))


if __name__ == '__main__': main()
