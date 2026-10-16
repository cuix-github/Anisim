"""A solitary wave in the 2D reference: the control for the beach runs.

A tsunami near shore is one long hump of water, and the hump has an exact enough
description: the solitary wave of shallow-water theory (Boussinesq). Over still water
of depth h, a hump of height H has the profile

    eta(x) = H sech^2(k (x - x0)),   k = sqrt(3 H / (4 h^3)),

travels without changing shape at c = sqrt(g (h + H)), and carries a depth-averaged
horizontal velocity u = c eta / (h + eta); continuity then asks for a vertical
velocity w = -y du/dx (y from the bed). Seed the particles with eta as the free
surface and (u, w) as their velocity and the solver should move the hump at c and
keep its height; seed the same hump at rest and it must instead split into two half
height waves running apart, the textbook result, which is how the check can fail.

Three scenes, cells of one metre (g = 9.8 cells/s^2):

  flat    400 x 48 cells, h = 20, H = 6 (H/h = 0.3), the hump at x0 = 100 with its
          velocity field, run 11 s. Checks: the crest speed (a straight-line fit of
          crest position against time after the first second) within 6 percent of c;
          the crest height at the end at least 0.75 H; nothing above 0.15 H left
          behind at mid-run.
  rest    the same hump with no velocity, 6 s. Checks at mid-run: a crest runs left as
          well, at least 0.3 H high; the right-going crest is below 0.7 H. The
          negative control.
  beach   520 x 56 cells, the same hump at x0 = 80 over 20 m of water, a plane beach
          of slope 1:15 from x = 160 to the shoreline at x = 460, as a cut-cell level
          set solid (pfflip2d's solid_phi, solid_cut). Run until the wave breaks or
          reaches the shoreline. Two marks on the way, both read off the front face as
          the liquid's front (the x with a cell's worth of particles beyond it) in one-metre bands above
          still water around the crest: the face is "steep" when some three bands rise
          more than 2.5 m per metre of advance (about 70 degrees, Grilli's vertical
          face within the noise of a particle surface) under a crest that has shoaled
          to at least 1.1 H, and the wave has "turned over"
          when a band's front exceeds a lower band's by more than 1.5 m, the lip
          overhanging. Each mark's first frame is recorded with the crest's position,
          the still-water depth under it and its height. Checks: the flat scenes never
          overhang; the beach's front turns steep between the toe and the shoreline.
          Reported beside it, not gated: Grilli, Svendsen and Subramanya (1997) put a
          solitary wave's breaking (the front face vertical) at H_b / h_b = 0.841 exp(6.421 S0) with
          S0 = 1.521 s / sqrt(H/h), the slope parameter; S0 = 0.185 here, a plunging
          breaker, H_b / h_b about 2.8.

  fine    `solitary_wave.py fine` runs the beach alone at half-metre cells (every
          length doubled in cells, g = 19.6 cells/s^2, the same seconds): does the
          breaking depth move with the resolution? Its own receipt_fine.json.

Receipts: results/2d/solitary_wave/receipt.json; frame strips flat.png, rest.png,
beach.png and a close-up of the last two seconds before the break, beach_zoom.png.
Run from reference2d:  ..\\.venv\\Scripts\\python.exe solitary_wave.py  (three minutes,
one core; `fine` about fifteen). Exits nonzero on failure.
"""
import collections
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pfflip2d import Sim, draw, save_strip  # noqa: E402

RESULTS = os.path.join(HERE, "..", "results", "2d", "solitary_wave")
G = 9.8                                   # m/s^2; in cells/s^2 it is G * cells per metre


def profile(x, h, H, x0):
    """The solitary wave's free surface above still water, and its wavenumber k."""
    k = np.sqrt(3.0 * H / (4.0 * h ** 3))
    return H / np.cosh(k * (x - x0)) ** 2, k


def velocity(x, y, h, H, x0, g=G):
    """The wave's velocity at (x, y): depth-averaged u, and the w continuity asks for."""
    c = np.sqrt(g * (h + H))
    eta, k = profile(x, h, H, x0)
    u = c * eta / (h + eta)
    deta = -2.0 * k * eta * np.tanh(k * (x - x0))
    du = c * h * deta / (h + eta) ** 2
    return u, -y * du, c


def beach_bed(x, toe, slope):
    """The bed's height above the flat bottom: zero to the toe, then the slope."""
    return np.maximum(0.0, (x - toe) * slope)


def beach_phi(nx, ny, toe, slope):
    """Signed distance to the beach at the grid nodes, negative inside the sand."""
    i, j = np.meshgrid(np.arange(nx + 1, dtype=float), np.arange(ny + 1, dtype=float), indexing="ij")
    return (j - beach_bed(i, toe, slope)) * np.cos(np.arctan(slope))


def make(nx, ny, h, H, x0, moving=True, beach=None, g=G, **kw):
    """A seeded, calibrated simulation: the hump as the free surface, the wave's
    velocity on the liquid particles (or none), a beach as a cut-cell solid if given."""
    phi = beach_phi(nx, ny, *beach) if beach else None
    sim = Sim(nx, ny, g=-g, solid_phi=phi, solid_cut=phi is not None, **kw)
    eta_of = lambda x: profile(x, h, H, x0)[0]   # noqa: E731
    sim.seed(lambda x, y: y < h + eta_of(x))
    sim.calibrate()
    if moving:
        liquid = sim.typ == 1
        u, w, _ = velocity(sim.pos[:, 0], sim.pos[:, 1], h, H, x0, g)
        sim.vel[liquid, 0] = u[liquid]
        sim.vel[liquid, 1] = w[liquid]
    return sim


def surface(sim):
    """The free surface per one-cell column: the highest liquid particle, raw and
    smoothed over five columns; -inf where a column holds no liquid."""
    liq = sim.pos[sim.typ == 1]
    col = np.clip(liq[:, 0].astype(np.int64), 0, sim.nx - 1)
    top = np.full(sim.nx, -np.inf)
    np.maximum.at(top, col, liq[:, 1])
    wet = np.isfinite(top)
    sm = top.copy()
    k = np.ones(5) / 5.0
    filled = np.where(wet, top, 0.0)
    sm[wet] = (np.convolve(filled, k, "same") / np.maximum(np.convolve(wet.astype(float), k, "same"), 1e-9))[wet]
    return top, sm


def crest(top, lo=0, hi=None):
    """The crest between columns lo and hi: its position (parabola through the three
    highest columns) and its height."""
    hi = len(top) if hi is None else hi
    seg = top[lo:hi]
    i = int(np.argmax(seg))
    x = lo + i + 0.5
    if 0 < i < len(seg) - 1 and np.all(np.isfinite(seg[i - 1:i + 2])):
        a, b, c = seg[i - 1], seg[i], seg[i + 1]
        den = a - 2 * b + c
        if den < 0:
            x += 0.5 * (a - c) / den
    return float(x), float(seg[i])


def face(pos, liquid, h, crest_x, s=1):
    """The wave's front face, from the liquid's front in one-cell bands above still
    water, 10 m behind the crest to 40 m ahead of it. A band's front is the x with one
    cell's worth of particles (4 s) beyond it, so a stray particle ahead does not move
    it and a long band is not read short; a band with fewer than five cells' worth
    (the tip) ends the count. Returns the face's steepness,
    the largest rise over advance across any three bands (cells per cell, 99 when
    vertical or beyond), and the overhang, the largest excess of a band's front over
    any lower band's (zero for a face that leans forward without turning over)."""
    p = pos[liquid]
    p = p[(p[:, 0] > crest_x - 10 * s) & (p[:, 0] < crest_x + 40 * s) & (p[:, 1] >= h)]
    if len(p) < 20:
        return 0.0, 0.0
    band = ((p[:, 1] - h) / s).astype(np.int64)
    k = int(4 * s)                      # the particles of one cell of a band
    fronts = []
    for b in range(int(band.max()) + 1):
        xs = np.sort(p[band == b, 0])
        if len(xs) < 5 * k:
            break
        fronts.append(float(xs[-k]))    # a cell's worth of particles stands behind it
    steep = 0.0
    for j in range(len(fronts) - 3):
        dx = fronts[j] - fronts[j + 3]
        steep = max(steep, 99.0 if dx <= 1e-6 else min(99.0, 3.0 * s / dx))
    over = max([fronts[j] - min(fronts[:j]) for j in range(1, len(fronts))] + [0.0])
    return steep, over


def close_up(sim, cx, half=40, scale=6):
    """The particles around column cx, `half` cells either side, the whole height."""
    img = draw(sim, scale=scale)
    x0 = int(max(0, (cx - half) * scale))
    x1 = int(min(img.size[0], (cx + half) * scale))
    return img.crop((x0, 0, x1, img.size[1]))


def run_flat(moving, seconds, label):
    nx, ny, h, H, x0 = 400, 48, 20.0, 6.0, 100.0
    sim = make(nx, ny, h, H, x0, moving=moving)
    _, _, c = velocity(np.zeros(1), np.zeros(1), h, H, x0)
    dt_frame, frames = 0.1, int(round(seconds / 0.1))
    ts, xs, heights, left, shots, over = [], [], [], [], [], 0.0
    t0 = time.time()

    def cb(f, s):
        nonlocal over
        t = (f + 1) * dt_frame
        _, top = surface(s)
        x, y = crest(top, lo=int(x0) - 20)          # the right-going crest
        xl, yl = crest(top, hi=int(x0) - 20)        # whatever runs left
        ts.append(t)
        xs.append(x)
        heights.append(y - h)
        left.append(yl - h)
        over = max(over, face(s.pos, s.typ == 1, h, x)[1])
        if f % max(1, frames // 5) == 0 or f == frames - 1:
            shots.append(draw(s, scale=2))
    sim.run(frames, dt_frame, on_frame=cb)
    wall = time.time() - t0
    os.makedirs(RESULTS, exist_ok=True)
    save_strip(shots, os.path.join(RESULTS, f"{label}.png"), cols=1)
    ts, xs, heights, left = map(np.array, (ts, xs, heights, left))
    fit = ts >= 1.0
    speed = float(np.polyfit(ts[fit], xs[fit], 1)[0]) if fit.sum() > 2 else float("nan")
    mid = len(ts) // 2
    return {"label": label, "grid": [nx, ny], "h": h, "H": H, "x0": x0, "moving": moving,
            "seconds": seconds, "wall_s": round(wall, 1), "particles": int(len(sim.pos)),
            "c_analytic": round(float(c), 3), "crest_speed": round(speed, 3),
            "speed_ratio": round(speed / float(c), 4),
            "crest_height": {"first": round(float(heights[0]), 2), "mid": round(float(heights[mid]), 2),
                             "last": round(float(heights[-1]), 2)},
            "left_height_mid": round(float(left[mid]), 2), "left_height_max": round(float(left.max()), 2),
            "crest_x": {"first": round(float(xs[0]), 1), "last": round(float(xs[-1]), 1)},
            "overhang_max_m": round(float(over), 2)}


PRODUCTION = dict(st=True, escape=True, adapt=True)   # the 3D runs' config, in the reference


def run_beach(cells_per_metre=1, seconds=26.0, label="beach", config=None, dt_frame=0.1, cfl=0.5):
    """The beach scene; config adds Sim options (PRODUCTION: ST sampling, escaped
    particles, adaptive air), cfl and dt_frame set the substep (plain: CFL 0.5 under
    dt_frame / 2; production: dt_frame 0.2 and a CFL cap of 10 make every substep 0.1 s,
    the 3D runs' step)."""
    s = cells_per_metre
    nx, ny, h, H, x0 = 520 * s, 56 * s, 20.0 * s, 6.0 * s, 80.0 * s
    toe, slope = 160.0 * s, 1.0 / 15.0
    shoreline = toe + h / slope
    g = G * s
    sim = make(nx, ny, h, H, x0, moving=True, beach=(toe, slope), g=g, **(config or {}))
    bed_of = lambda x: beach_bed(x, toe, slope)   # noqa: E731
    frames = int(round(seconds / dt_frame))
    shots, log, hit, steep = [], [], None, None
    recent = collections.deque(maxlen=20)
    t0 = time.time()

    def metres(v):
        return round(float(v) / s, 2)

    def cb(f, sim_):
        nonlocal hit, steep
        t = (f + 1) * dt_frame
        _, top = surface(sim_)
        x, y = crest(top, lo=20 * s)
        steepness, over = face(sim_.pos, sim_.typ == 1, h, x, s)
        entry = {"t": round(t, 2), "crest_x": metres(x), "crest_height": metres(y - h),
                 "depth_under_crest": metres(h - float(bed_of(x))), "face_steepness": round(steepness, 2),
                 "overhang_m": metres(over)}
        log.append(entry)
        if steep is None and steepness >= 2.5 and (y - h) >= 1.1 * H and toe < x < shoreline:
            steep = dict(entry)             # a shoaled crest (1.1 H) with a steep face, not noise
        if hit is None and over >= 1.5 * s and toe < x < shoreline:
            hit = dict(entry)
        recent.append((t, close_up(sim_, x, half=40 * s, scale=max(1, 6 // s))))
        if (f % 20 == 0 and t < 10) or (f % 6 == 0 and t >= 10) or f == frames - 1:
            shots.append(draw(sim_, scale=max(1, 2 // s)))
        done = hit is not None and steep is not None and t > max(hit["t"], steep["t"]) + 1.0
        if done or x > shoreline - 10 * s:          # both marks set, or the wave ran up
            raise StopIteration
    try:
        sim.run(frames, dt_frame, cfl=cfl, on_frame=cb)
    except StopIteration:
        pass
    wall = time.time() - t0
    os.makedirs(RESULTS, exist_ok=True)
    save_strip(shots[-8:], os.path.join(RESULTS, f"{label}.png"), cols=1)
    zoom = list(recent)[::2][-8:]
    save_strip([im for _, im in zoom], os.path.join(RESULTS, f"{label}_zoom.png"), cols=2)
    S0 = 1.521 * slope / np.sqrt(H / h)
    ratio = 0.841 * np.exp(6.421 * S0)
    return {"label": label, "cells_per_metre": s, "grid": [nx, ny], "h_m": 20.0, "H_m": 6.0, "x0_m": 80.0,
            "toe_m": 160.0, "slope": "1:15", "shoreline_m": 460.0, "wall_s": round(wall, 1),
            "config": {k: bool(v) for k, v in (config or {}).items()}, "dt_frame": dt_frame, "cfl": cfl,
            "particles": int(len(sim.pos)), "frames_run": len(log),
            "steep_front": steep, "turned_over": hit,
            "crest_at_toe_passage": next((e for e in log if e["crest_x"] >= 160.0), None),
            "grilli": {"S0": round(float(S0), 3), "Hb_over_hb": round(float(ratio), 2),
                       "breaker": "plunging" if 0.025 < S0 < 0.3 else ("spilling" if S0 <= 0.025 else "surging"),
                       "measured_Hb_over_hb_at_steep": (round(steep["crest_height"] / max(steep["depth_under_crest"], 1e-6), 2)
                                                        if steep else None),
                       "measured_Hb_over_hb_turned_over": (round(hit["crest_height"] / max(hit["depth_under_crest"], 1e-6), 2)
                                                           if hit else None)},
            "zoom_times": [round(t, 1) for t, _ in zoom],
            "log_every_2s": [e for i, e in enumerate(log) if i % 20 == 0] + log[-1:]}


def report_beach(beach):
    b, st = beach["turned_over"], beach["steep_front"]
    if st:
        print(f"    front face steep (>= 2.5 m rise per m) at t {st['t']} s, crest x {st['crest_x']} m over "
              f"{st['depth_under_crest']} m of water, crest {st['crest_height']} m "
              f"(H_b/h_b {beach['grilli']['measured_Hb_over_hb_at_steep']})")
    else:
        print(f"    the front never turned steep in {beach['frames_run']} frames")
    if b:
        print(f"    lip overhangs by {b['overhang_m']} m at t {b['t']} s, crest x {b['crest_x']} m over "
              f"{b['depth_under_crest']} m of water, crest {b['crest_height']} m "
              f"(H_b/h_b {beach['grilli']['measured_Hb_over_hb_turned_over']})")
    else:
        print(f"    no overhang above 1.5 m in {beach['frames_run']} frames")
    print(f"    Grilli: {beach['grilli']['breaker']}, vertical face at H_b/h_b {beach['grilli']['Hb_over_hb']}; "
          f"{beach['particles']:,} particles, {beach['wall_s']} s")
    ok = st is not None and 160.0 < st["crest_x"] < 460.0
    print("    " + ("OK" if ok else "FAIL"))
    return ok


def main(argv):
    os.makedirs(RESULTS, exist_ok=True)
    if argv[:1] == ["fine"]:
        print("=== beach at half-metre cells")
        beach = run_beach(cells_per_metre=2, label="beach_fine")
        ok = report_beach(beach)
        with open(os.path.join(RESULTS, "receipt_fine.json"), "w") as f:
            json.dump({"beach_fine": beach}, f, indent=1)
        print("SOLITARY WAVE (fine):", "PASS" if ok else "FAIL")
        sys.exit(0 if ok else 1)
    if argv[:1] == ["production"]:
        print("=== beach in the production config (ST, escape, adaptive air) at substeps of 0.1 s")
        beach = run_beach(label="beach_production", config=PRODUCTION, dt_frame=0.2, cfl=10.0)
        ok = report_beach(beach)
        with open(os.path.join(RESULTS, "receipt_production.json"), "w") as f:
            json.dump({"beach_production": beach}, f, indent=1)
        print("SOLITARY WAVE (production):", "PASS" if ok else "FAIL")
        sys.exit(0 if ok else 1)
    print("=== flat: the hump with its velocity field")
    flat = run_flat(True, 11.0, "flat")
    print(f"    crest speed {flat['crest_speed']} against c {flat['c_analytic']} ({flat['speed_ratio']:.3f}); "
          f"height {flat['crest_height']['first']} -> {flat['crest_height']['last']} of H 6; left behind at mid-run "
          f"{flat['left_height_mid']}; largest overhang {flat['overhang_max_m']} m; {flat['wall_s']} s")
    ok_flat = (abs(flat["speed_ratio"] - 1.0) <= 0.06 and flat["crest_height"]["last"] >= 0.75 * 6
               and abs(flat["left_height_mid"]) <= 0.15 * 6 and flat["overhang_max_m"] < 1.5)
    print("    " + ("OK" if ok_flat else "FAIL"))
    print("=== rest: the same hump at rest (must split)")
    rest = run_flat(False, 6.0, "rest")
    print(f"    at mid-run: right crest {rest['crest_height']['mid']} of H 6, left crest {rest['left_height_mid']}; "
          f"{rest['wall_s']} s")
    ok_rest = rest["left_height_mid"] >= 0.3 * 6 and rest["crest_height"]["mid"] <= 0.7 * 6 \
        and rest["overhang_max_m"] < 1.5
    print("    " + ("OK" if ok_rest else "FAIL"))
    print("=== beach: the wave onto a 1:15 slope")
    beach = run_beach()
    ok_beach = report_beach(beach)
    with open(os.path.join(RESULTS, "receipt.json"), "w") as f:
        json.dump({"flat": flat, "rest": rest, "beach": beach}, f, indent=1)
    ok = ok_flat and ok_rest and ok_beach
    print("SOLITARY WAVE:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
