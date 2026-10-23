"""Independent laws and actual Bifrost graph checks for opt-in spray events.

Run with the main checkout's venv. All arrays/logs go to the supplied --out;
no production cache is read. Exit nonzero if any gate fails.
"""
import argparse
import json
import os
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE), str(HERE.parent / "reference2d")]
import spray_events as sampler
from pfflip2d import Sim
import validate_spray_node as legacy
import validate_droplet_node as drops
from validate_level_set_node import load


def f32(d):
    return {k: float(np.float32(v)) for k, v in d.items()}


D = f32(dict(Sim.DROPLET_DEFAULTS, integrator=2, dt_max=.005, metres_per_cell=.25))
E = f32(dict(Sim.EMIT_DEFAULTS, model=2, grad_threshold=3, rate=50,
             cloud_size=16, max_per_cell=512, max_active=100000, seed=137))


def ports(e=E, d=D):
    p = dict(drops.PORTS, droplet_model=int(d["integrator"]),
             spray_dt_max=drops.pf(d["dt_max"]), spray_cfl=drops.pf(.5))
    for name in ("r_median", "r_sigma", "evap_k", "r_min"):
        p["droplet_" + name] = drops.pf(d[name])
    for name in ("metres_per_cell", "air_nu"):
        p[name] = drops.pf(d[name])
    ints = ("model", "cloud_size", "max_per_cell", "max_active", "seed")
    p.update({"emit_" + k: int(v) if k in ints else drops.pf(v) for k, v in e.items()})
    return p


def laws():
    assert Sim(8,8,escape=True,droplet=D,emit=dict(model=0),spray_cfl=.5).emit is None
    n = 300000
    k = sampler.keys(137, 12, np.arange(n, dtype=np.uint64))
    counts = {}
    for lam in (.05, .5, 1.2, 5.75):
        draws = sampler.requested_events(lam, k)
        mean = float(draws.mean())
        se = np.sqrt((lam % 1) * (1-lam % 1) / n)
        assert abs(mean-lam) < 6*se, (lam, mean, se)
        counts[str(lam)] = dict(mean=mean, target=lam, standard_error=float(se))
    # Expected total event count is unchanged when an interval is subdivided.
    temporal = {}
    for steps in (1, 4, 16):
        total = sum(sampler.requested_events(2.8/steps,
                    sampler.keys(137, step, np.arange(30000))) for step in range(steps))
        temporal[str(steps)] = float(total.mean())
        assert abs(total.mean()-2.8) < .04
    radius = {}
    for N in (1, 8, 64):
        sigma, R = .4, 3e-4
        key = sampler.keys(137, 0, np.arange(n), 0, N)
        z = np.sqrt(-2*np.log(1-sampler.uniform(key, 1))) * np.cos(2*np.pi*sampler.uniform(key, 2))
        r = R*N**(-1/3)*np.exp(sigma*z-1.5*sigma*sigma)
        ratio = float(np.mean(N*r**3/R**3))
        se = np.sqrt((np.exp(9*sigma*sigma)-1)/n)
        assert abs(ratio-1) < 6*se
        assert abs(z.mean()) < .015 and abs(z.var()-1) < .02
        radius[str(N)] = dict(expected_volume_ratio=ratio, standard_error=float(se))
    geometry = {}
    for dim in (2, 3):
        p = sampler.ball(k, 3, dim)
        mean_square = float((p*p).sum(1).mean())
        assert np.linalg.norm(p.mean(0)) < .01
        assert (np.linalg.norm(p, axis=1) <= 1+1e-14).all()
        assert abs(mean_square-dim/(dim+2)) < .005
        geometry[str(dim)] = dict(mean_square_radius=mean_square, target=dim/(dim+2))
    bads = [dict(E, rate=-1), dict(E, cloud_size=0), dict(E, seed=.5),
            dict(E, max_active=np.nan), dict(E, radius_max=0)]
    for bad in bads:
        try:
            sampler.validate(bad, D)
        except ValueError:
            continue
        raise AssertionError("invalid configuration accepted")
    return dict(explicit_model_zero_disabled=True, event_count=counts, temporal_expectation=temporal, radius=radius,
                geometry=geometry, invalid_rejected=len(bads))


def flight():
    # A smooth, spatially varying air field; the oracle integrates the continuous
    # particle ODE with classical RK4, not the solver's exponential discretization.
    radius, nu, ratio, gravity, T = 2e-4, 1.5e-5, .001, -.05, .4
    def rhs(t, y):
        vg = np.array([.3 + .04*y[0], .02*y[1]])
        speed = np.linalg.norm(vg-y[2:])
        re = 2*radius*speed/nu
        cdv = 12*nu/radius*(1+.15*re**.687) if re <= 1000 else .44*speed
        rate = .375*ratio*cdv/radius
        return np.r_[y[2:], rate*(vg-y[2:])+[0, gravity]]
    initial = np.array([2., 4., .08, .01])
    def rk4(steps):
        y, h = initial.copy(), T/steps
        for j in range(steps):
            a = rhs(j*h, y); b = rhs((j+.5)*h, y+.5*h*a)
            c = rhs((j+.5)*h, y+.5*h*b); d = rhs((j+1)*h, y+h*c)
            y += h/6*(a+2*b+2*c+d)
        return y
    truth = rk4(8000)
    oracle_error = float(np.linalg.norm(truth-rk4(4000)))
    assert oracle_error < 1e-10
    errors = []
    for h in (.04, .02, .01):
        s = Sim(8, 8, escape=True, g=gravity,
                droplet=dict(D, air_nu=nu, metres_per_cell=1, dt_max=h), spray_cfl=.5)
        s.pos, s.vel = initial[None, :2].copy(), initial[None, 2:].copy()
        s.radius = np.array([radius])
        x = np.arange(9)[:, None]
        y = np.arange(9)[None, :]
        s.u = np.broadcast_to(.3+.04*x, (9, 8)).copy()
        s.v = np.broadcast_to(.02*y, (8, 9)).copy()
        s._droplet_fly(np.array([True]), np.array([T]))
        errors.append(float(np.linalg.norm(np.r_[s.pos[0], s.vel[0]]-truth)))
    ratios = [errors[i]/errors[i+1] for i in range(2)]
    assert min(ratios) > 3.2 and max(ratios) < 4.8, (errors, ratios)
    assert errors[-1] < 1e-4
    good, anchors = drops.check_anchors()
    assert good, anchors
    return dict(nonuniform_errors=errors, refinement_ratios=ratios,
                oracle_refinement_error=oracle_error, legacy_physics=anchors)


def graph_cases(out):
    legacy.OUT = str(out / "graph")
    os.environ["PFFLIP_SPRAY_AUDIT"] = "1"
    a, off = legacy.check_off()
    b, lock = legacy.check_lockstep()
    assert a and b
    baseline = out.parent / "control"
    if baseline.exists():
        for tag in ("off_zero_out", "off_unset_out", "lock_out"):
            for p in (baseline / tag).glob("*.npy"):
                assert np.array_equal(np.load(p), np.load(out/"graph"/tag/p.name)), str(p)
    results = dict(legacy_off=off, legacy_lockstep=lock)
    sim = legacy.tank(wind=30, shear=10)
    idx = legacy.plant(sim, 60, 26, 36)
    radii = np.zeros(len(sim.pos)); flags = np.zeros(len(sim.pos))
    radii[idx[:40]] = np.geomspace(2e-5, 2e-3, 40)
    radii[idx[40:]] = 1e-4; flags[idx[40:]] = 1
    cases = {"clouds": E, "budget": dict(E, max_active=27, max_per_cell=3),
             "boundary": dict(E, spread=30), "no_room": dict(E, max_active=0)}
    for name, config in cases.items():
        dest, ch, stats = legacy.node_step(name, sim, radii, flags, step_index=7, **ports(config))
        ref, gone, emitted = legacy.twin(sim, ch, radii, flags, 7, config, .5, droplet=D)
        good, rec = legacy.compare(name, dest, ref)
        assert good and rec["escape_mismatches"] == 0
        node = stats["spray_events"][0]
        for key, value in ref.emission_stats.items():
            if key.endswith("r3"):
                assert np.isclose(node[key], value, rtol=1e-10, atol=1e-30), (key,node[key],value)
            else:
                assert node[key] == value, (name, key,node[key],value)
        assert node["requested"] == sum(node[k] for k in ("accepted", "rejected_cap", "rejected_geometry"))
        assert node["alive"] == rec["spray"]
        if name == "clouds":
            assert emitted > 50
            assert (load(dest, "radius_out")[-emitted:] > 0).all()
        results[name] = dict(comparison=rec, accounting=node)
    # Invalid values must cause the actual graph to refuse, not quietly default.
    try:
        legacy.node_step("invalid", sim, radii, flags, **ports(dict(E, event_radius=-1)))
    except SystemExit as exc:
        assert "refused" in str(exc)
        results["invalid_graph_refused"] = str(exc)
    else:
        raise AssertionError("actual graph accepted negative event radius")
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--check", choices=("all", "laws", "flight", "graph"), default="all")
    args = parser.parse_args()
    args.out = args.out.resolve()
    args.out.mkdir(parents=True, exist_ok=True)
    receipt = dict(check=args.check, passed=False)
    start = time.perf_counter()
    try:
        for name, fn in (("laws", laws), ("flight", flight), ("graph", lambda: graph_cases(args.out))):
            if args.check in ("all", name):
                receipt[name] = fn()
                print(name, "PASS", flush=True)
        receipt["passed"] = True
    finally:
        receipt["seconds"] = time.perf_counter()-start
        (args.out / ("receipt_"+args.check+".json")).write_text(json.dumps(receipt, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
