"""Milestone 6.5a stage-3a gate: the adaptive projection inside the 2D loop.

Both junction treatments that survived stage 3a (adaptive_projection.py):
naive, the first-order conservative stencil the PF-FLIP paper adopts, and
cleanup, the second-order solve plus a conservative projection of its
leftover. Each runs against the uniform projection on the same scenes, and
every bar is set by the uniform solver's own particle-seed spread: two runs
of the SAME uniform solver from different random seedings differ by some
amount, and an arm passes when it sits closer to the uniform run than twice
that. Three checks plus a cost record:

  books   one projection of a real dam-break state: every uniform face,
          inside coarse cells included, must leave with max |divergence|
          within 10x of the uniform projection's own leftover at the same
          tolerance (1e-6).
  tank    a half-full tank at rest, three seeds, 40 frames: rms particle
          speed per phase over frames 30-40, paired by seed against uniform.
          cleanup is gated; naive's parasitic currents (the known first-order
          junction error, dt g / 3 in hydrostatics) are measured and
          reported, the input to the 3D decision, not gated. Seed 1 runs on
          to frame 200 and reports frames 190-200 too, because naive's
          stirring had not levelled off by frame 40. (The liquid's own rms
          grows over those 200 frames in every arm, uniform included: a
          property of the 2D reference at rest, not of the adaptive solve.)
  dam     a 160x80 dam break, 41 frames: front, centroid, occupied area and
          per-phase kinetic energy against the uniform run, both arms
          gated on the seed-spread rule; no blow-up.
  cost    unknowns saved and pressure iterations per projection (main solve
          and cleanup solve), recorded, not gated.

Run from reference2d:  ..\\.venv\\Scripts\\python.exe validate_adaptive_projection.py
"""
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from adaptive_projection import AdaptiveProjection  # noqa: E402
from pfflip2d import Sim  # noqa: E402
from validate_adaptive_solve import dam_state  # noqa: E402

RES = os.path.join(HERE, "..", "results", "2d", "adaptive")
ARMS = ("naive", "cleanup")
SEEDS = (1, 2, 3)


def div_max(s):
    d = (s.u[1:, :] - s.u[:-1, :]) + (s.v[:, 1:] - s.v[:, :-1])
    return float(np.abs(d).max())


def costs(proj):
    lg = proj.log
    return {"projections": len(lg),
            "unknowns_saved": float(1 - np.mean([g["unknowns"] / g["uniform"] for g in lg])),
            "iters_main": float(np.mean([g["iters"] for g in lg])),
            "iters_cleanup": float(np.mean([g["iters_cleanup"] for g in lg]))}


def check_books():
    print("=== books: one projection of a real dam-break state, every uniform face")
    dt = 0.06
    base, _, _ = dam_state()

    def clone():
        c = Sim(base.nx, base.ny)
        c.__dict__.update({k: (v.copy() if isinstance(v, np.ndarray) else v)
                           for k, v in base.__dict__.items()})
        return c
    ref = clone()
    ref.project(dt)
    floor = div_max(ref)
    rows = {"uniform": floor}
    ok = True
    for arm in ARMS:
        c = clone()
        c.project = AdaptiveProjection(c, junction=arm)
        c.project(dt)
        rows[arm] = div_max(c)
        print(f"    {arm}: max |div| {rows[arm]:.2e} (uniform {floor:.2e})")
        ok = ok and rows[arm] < 10 * floor
    print("    PASS" if ok else "    FAIL")
    return ok, rows


def check_tank():
    print("=== tank: half-full tank at rest, rms speed per phase, paired by seed")
    out = {}
    for arm in ("uniform",) + ARMS:
        for seed in SEEDS:
            s = Sim(64, 64)
            s.seed(lambda x, y: y < 32, seed=seed)
            s.calibrate()
            proj = None
            if arm != "uniform":
                proj = AdaptiveProjection(s, junction=arm)
                s.project = proj
            air, liq = [], []

            def cb(f, q):
                sp2 = (q.vel ** 2).sum(1)
                air.append(float(np.sqrt(sp2[q.typ == 0].mean())))
                liq.append(float(np.sqrt(sp2[q.typ == 1].mean())))
            s.run(200 if seed == 1 else 40, 0.05, on_frame=cb)
            out[(arm, seed)] = {"air": float(np.mean(air[30:40])),
                                "liquid": float(np.mean(liq[30:40])),
                                "air_curve": air, "liquid_curve": liq}
            if seed == 1:
                out[(arm, seed)].update(air_190_200=float(np.mean(air[190:200])),
                                        liquid_190_200=float(np.mean(liq[190:200])))
            if proj:
                out[(arm, seed)].update(costs(proj))
    spread = {ph: max(out[("uniform", k)][ph] for k in SEEDS) /
              min(out[("uniform", k)][ph] for k in SEEDS) - 1 for ph in ("air", "liquid")}
    rows = {"uniform_seed_spread": spread, "arms": {}}
    print(f"    uniform seed spread: air {spread['air']:.1%}, liquid {spread['liquid']:.1%}")
    ok = True
    for arm in ARMS:
        ratio = {ph: [out[(arm, k)][ph] / out[("uniform", k)][ph] for k in SEEDS]
                 for ph in ("air", "liquid")}
        gated = arm == "cleanup"
        inside = all(max(ratio[ph]) <= 1 + 2 * spread[ph] for ph in ratio)
        rows["arms"][arm] = {"ratio": ratio, "inside_floor": inside,
                             "per_seed": {k: out[(arm, k)] for k in SEEDS}}
        print(f"    {arm}: air rms x{min(ratio['air']):.2f}-{max(ratio['air']):.2f} of uniform, "
              f"liquid x{min(ratio['liquid']):.2f}-{max(ratio['liquid']):.2f} | "
              f"iters {out[(arm, 1)]['iters_main']:.1f} + cleanup "
              f"{out[(arm, 1)]['iters_cleanup']:.1f} | "
              + ("gated: " + ("inside" if inside else "OUTSIDE") if gated
                 else "reported: " + ("inside" if inside else "outside")) + " the floor")
        if gated:
            ok = ok and inside
    for arm in ("uniform",) + ARMS:
        r = out[(arm, 1)]
        print(f"    seed 1 at frames 190-200: {arm} air rms {r['air_190_200']:.3f}, "
              f"liquid {r['liquid_190_200']:.3f}")
    rows["uniform"] = {k: out[("uniform", k)] for k in SEEDS}
    print("    PASS" if ok else "    FAIL")
    return ok, rows


def check_dam():
    print("=== dam: 160x80 dam break, arms vs uniform inside the seed spread")
    nx, ny = 160, 80

    def run(arm, seed):
        s = Sim(nx, ny)
        s.seed(lambda x, y: (x < 40) & (y < 56), seed=seed)
        s.calibrate()
        proj = None
        if arm != "uniform":
            proj = AdaptiveProjection(s, junction=arm)
            s.project = proj
        rec = {"front": [], "centroid": [], "occupied": [], "ke_air": [], "ke_liquid": []}

        def cb(f, q):
            liq = q.pos[q.typ == 1]
            sp2 = (q.vel ** 2).sum(1)
            if sp2[q.typ == 0].max() > 200.0 ** 2:
                raise RuntimeError(f"air speed blow-up at frame {f}")
            rec["front"].append(float(np.percentile(liq[:, 0], 99.5)))
            rec["centroid"].append(liq.mean(0).tolist())
            h = np.zeros((nx, ny))
            np.add.at(h, (np.clip(liq[:, 0].astype(int), 0, nx - 1),
                          np.clip(liq[:, 1].astype(int), 0, ny - 1)), 1.0)
            rec["occupied"].append(int((h >= 2).sum()))   # cells at least half full
            rec["ke_air"].append(float(sp2[q.typ == 0].mean()))
            rec["ke_liquid"].append(float(sp2[q.typ == 1].mean()))
        try:
            s.run(41, 0.06, on_frame=cb)
        except RuntimeError as e:
            rec["blowup"] = str(e)
        if proj:
            rec.update(costs(proj))
        return rec

    def dist(a, b):
        return {"front": float(np.abs(np.subtract(a["front"], b["front"])).max()),
                "centroid": float(np.linalg.norm(np.subtract(a["centroid"], b["centroid"]),
                                                 axis=1).max()),
                "occupied": float(np.abs(np.divide(a["occupied"], b["occupied"]) - 1).max()),
                "ke_air": abs(np.mean(a["ke_air"]) / np.mean(b["ke_air"]) - 1),
                "ke_liquid": abs(np.mean(a["ke_liquid"]) / np.mean(b["ke_liquid"]) - 1)}

    uni = {k: run("uniform", k) for k in SEEDS}
    floor = {}
    for k in SEEDS[1:]:
        for m, v in dist(uni[k], uni[1]).items():
            floor[m] = max(floor.get(m, 0.0), v)
    print("    uniform seed floor: " + ", ".join(f"{m} {v:.3g}" for m, v in floor.items()))
    rows = {"floor": floor, "arms": {}}
    ok = True
    for arm in ARMS:
        r = run(arm, 1)
        if "blowup" in r:
            print(f"    {arm}: {r['blowup']}")
            rows["arms"][arm] = {"blowup": r["blowup"]}
            ok = False
            continue
        d = dist(r, uni[1])
        inside = all(d[m] <= 2 * floor[m] for m in d)
        rows["arms"][arm] = {"dist": d, "inside_floor": inside,
                             **{k: r[k] for k in ("unknowns_saved", "iters_main",
                                                  "iters_cleanup")}}
        print(f"    {arm}: " + ", ".join(f"{m} {v:.3g}" for m, v in d.items())
              + f" | saved {r['unknowns_saved']:.0%}, iters {r['iters_main']:.1f}"
              f" + cleanup {r['iters_cleanup']:.1f} | "
              + ("inside" if inside else "OUTSIDE") + " 2x the floor")
        ok = ok and inside
    print("    PASS" if ok else "    FAIL")
    return ok, rows


def main():
    os.makedirs(RES, exist_ok=True)
    ok_b, books = check_books()
    ok_t, tank = check_tank()
    ok_d, dam = check_dam()
    for arm in ARMS:
        t, d = tank["arms"][arm]["per_seed"][1], dam["arms"][arm]
        if "blowup" not in d:
            print(f"cost {arm}: tank {t['iters_main'] + t['iters_cleanup']:.1f} iters "
                  f"per projection, dam {d['iters_main'] + d['iters_cleanup']:.1f}; "
                  f"unknowns saved tank {t['unknowns_saved']:.0%}, dam {d['unknowns_saved']:.0%}")
    with open(os.path.join(RES, "stage3a_gate.json"), "w", newline="\n") as f:
        json.dump({"books": books, "tank": tank, "dam": dam}, f, indent=1)
    ok = ok_b and ok_t and ok_d
    print("ALL PASS" if ok else "FAILURES PRESENT")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
