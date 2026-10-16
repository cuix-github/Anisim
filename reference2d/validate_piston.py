"""The piston wavemaker in the 2D reference: pfflip2d.Sim(piston=(S, T, x0, width, top))
puts the paddle of a wave flume in the tank, a plate on the floor from x0 to x0 + width
reaching up to y = top, and slides it as x(t) = x0 + S/2 (1 - cos 2 pi t / T). The plate
stays below the ceiling so the air it displaces flows over it: a whole wall sweeping into
the sealed tank of incompressible fluid has no pressure solution (the divergence sums to
the swept rate), which is how the first version of this gate spun for hours on its first
frame. Behind the plate is air down to the floor (a dry-back paddle); in front, the water.
Three checks, cells of one metre, g = 9.8:

  rest        a plate with no stroke is the static solid of the same shape: the dam break
              of pfflip2d.py run 100 steps with piston=(0, 1, 8, 2, 56.5) and with that
              plate as a fixed cut-cell level set, every particle position and velocity
              and every face velocity identical to the bit.
  dispersion  a flat channel 600 x 40, h = 20 m, S = 2 m, T = 8 s, 50 s (six periods).
              Linear theory (Airy) ties the period to the wavelength through the
              dispersion relation w^2 = g k tanh(k h): k = 0.0708 per m, L = 88.7 m
              here. The surface at 50 s between 100 and 350 m from the plate's rest
              face, where the train is established and nothing has come back from the
              far wall, is measured for its wavelength (the mean spacing of its
              up-crossings) and its height (mean crest to trough). Checks: the
              wavelength within 10 percent of L; the height within 30 percent of the
              piston transfer function of Dean and Dalrymple, H / S = 2 (cosh 2kh - 1) /
              (sinh 2kh + 2kh) = 1.33, so 2.66 m; at least two waves in the window.
  fast        the hard case: a 4 m stroke at a period of 4 s (the plate up to 3.1 m/s,
              sweeping water four times as hard as the dispersion check) in a 200 x 60
              channel 20 m deep, the plate's top at 50 m, 12 s. After every frame: no
              liquid particle inside the plate, none behind it below its top (water can
              only get behind by going over), the fastest liquid particle under 30 m/s (the
              air squeezed over the plate's top runs faster and is reported), and the
              liquid count unchanged.

The surface is the highest liquid particle per one-metre column, smoothed over five
(solitary_wave.surface). At the default FLIP blend the reference's resting water is not
at rest: its surface layer thins into a mixing layer and that surface sinks by metres
over a minute (docs/open-questions.md), so the profile is read against its own line
over the window, smoothed again over nine columns, with crossings counted only once
the profile has swung a quarter of its amplitude either way; the window's level every
10 s is recorded beside the result. A count of liquid particles per column was tried
as a drift-proof surface and is noisier than the top (the mixing layer shuffles
particles sideways too). A third check, the train onto a 1:15 beach until its leading
crest breaks, was built and dropped: by the time the train reaches the slope the
surface is that mixing layer and no mark read off it means anything; the periodic
breaker is the node's test, whose tank stays quiet.

Receipts in results/2d/piston/ (receipt.json, the strip dispersion.png). Run from
reference2d:  ..\\.venv\\Scripts\\python.exe validate_piston.py [rest|dispersion|fast]
[--alpha a]; --alpha runs the dispersion check at another FLIP blend (both phases) with
the suffix _alpha<a> on the files. About nine minutes on one core, a line every five
simulated seconds. Exits nonzero on failure.
"""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from pfflip2d import Sim, draw, save_strip  # noqa: E402
from solitary_wave import surface  # noqa: E402

RESULTS = os.path.join(HERE, "..", "results", "2d", "piston")
G = 9.8
PLATE_X0, PLATE_W = 8.0, 2.0          # the plate's rest position and thickness, metres
SIM_KW = {}           # extra Sim arguments for the wave check (--alpha sets the FLIP blend)
SUFFIX = ""           # on every output file name


def wavenumber(T, h):
    """k from the dispersion relation w^2 = g k tanh(k h), by fixed-point iteration."""
    w = 2.0 * np.pi / T
    k = w * w / G
    for _ in range(100):
        k = w * w / (G * np.tanh(k * h))
    return float(k)


def transfer(k, h):
    """Dean and Dalrymple's piston wavemaker: wave height over stroke."""
    return 2.0 * (np.cosh(2 * k * h) - 1.0) / (np.sinh(2 * k * h) + 2 * k * h)


def progress(sim, f, dt_frame, t0):
    t = (f + 1) * dt_frame
    if f % 50 == 49:
        vmax = float(np.abs(sim.vel).max())
        print(f"    t {t:5.1f} s  plate x {sim.piston_x:6.2f}  solve iters {sim.iters:4d}  "
              f"max speed {vmax:5.2f} m/s  {time.time() - t0:5.0f} s", flush=True)


def check_rest():
    print("=== rest: a plate with no stroke is the static solid, bit for bit", flush=True)
    top = 56.5
    sims = [Sim(160, 80, piston=(0.0, 1.0, PLATE_X0, PLATE_W, top)),
            Sim(160, 80, solid_phi=Sim.plate_phi(160, 80, PLATE_X0, PLATE_W, top), solid_cut=True)]
    for sim in sims:
        sim.seed(lambda x, y: (x < 40) & (y < 56))
        sim.calibrate()
    a, b = sims
    worst = 0.0
    for k in range(100):
        for sim in sims:
            sim.step(0.03)
        worst = max(worst, float(np.abs(a.pos - b.pos).max()), float(np.abs(a.vel - b.vel).max()),
                    float(np.abs(a.u - b.u).max()), float(np.abs(a.v - b.v).max()))
    same = worst == 0.0 and len(a.pos) == len(b.pos)
    print(f"    100 steps of the dam break around the plate: max |diff| {worst:.1e} over positions, "
          f"velocities, faces; {len(a.pos):,} particles  " + ("OK" if same else "FAIL"), flush=True)
    return same, {"steps": 100, "max_diff": worst, "particles": len(a.pos)}


def wavelength_and_height(eta, smooth=9):
    """Mean up-crossing spacing and mean crest-to-trough height of a profile. The
    profile is smoothed over `smooth` columns first (a 1.7 percent loss on an 89 m
    wave at nine) and a crossing counts only once the profile has been below minus a
    quarter of the amplitude and rises above plus a quarter, so the surface's noise
    near the mean does not make waves of its own."""
    m = eta - eta.mean()
    if smooth > 1:
        half = smooth // 2
        m = np.convolve(np.pad(m, (half, half), mode="edge"), np.ones(smooth) / smooth, "valid")
    delta = 0.25 * 0.5 * (m.max() - m.min())
    state, ups = 0, []
    for i in range(1, len(m)):
        if m[i] < -delta:
            state = -1
        elif m[i] > delta and state == -1:
            j = i
            while j > 0 and m[j] >= 0.0:
                j -= 1
            ups.append(j + (-m[j]) / (m[j + 1] - m[j]))    # the crossing, refined linearly
            state = 1
    if len(ups) < 2:
        return float("nan"), float("nan"), 0
    L = float(np.mean(np.diff(ups)))
    heights = []
    for a, b in zip(ups[:-1], ups[1:]):
        seg = m[int(a):int(b) + 2]
        heights.append(float(seg.max() - seg.min()))
    return L, float(np.mean(heights)), len(heights)


def check_dispersion():
    print("=== dispersion: a flat channel, S = 2 m, T = 8 s, 50 s", flush=True)
    nx, ny, h, S, T = 600, 40, 20.0, 2.0, 8.0
    front = PLATE_X0 + PLATE_W                    # the plate's rest front face
    k = wavenumber(T, h)
    L_lin, H_lin = 2 * np.pi / k, S * transfer(k, h)
    sim = Sim(nx, ny, piston=(S, T, PLATE_X0, PLATE_W, h + 8.5), **SIM_KW)
    sim.seed(lambda x, y: (x > front) & (y < h))
    sim.calibrate()
    dt_frame, frames = 0.1, 500
    shots, series = [], []
    t0 = time.time()
    lo, hi = int(front + 100), int(front + 350)

    def cb(f, s):
        progress(s, f, dt_frame, t0)
        if f % 50 == 49:
            shots.append(draw(s, scale=2))
        if f % 100 == 99:
            _, top = surface(s)
            series.append(((f + 1) * dt_frame, top - h))
    sim.run(frames, dt_frame, on_frame=cb)
    wall = time.time() - t0
    os.makedirs(RESULTS, exist_ok=True)
    save_strip(shots, os.path.join(RESULTS, f"dispersion{SUFFIX}.png"), cols=1)
    # the window's profile against its own line: the level drifts (the docstring)
    eta = series[-1][1][lo:hi]
    xw = np.arange(lo, hi) + 0.5
    fin = np.isfinite(eta)
    trend = np.polyfit(xw[fin], eta[fin], 1)
    detrended = np.where(fin, eta - np.polyval(trend, xw), 0.0)
    L, H, n = wavelength_and_height(detrended)
    drift = [round(float(np.median(e[lo:hi][np.isfinite(e[lo:hi])])), 3) for _, e in series]
    ok = abs(L / L_lin - 1.0) <= 0.10 and abs(H / H_lin - 1.0) <= 0.30 and n >= 2
    print(f"    linear theory: k {k:.4f}/m, L {L_lin:.1f} m, H/S {transfer(k, h):.2f} -> H {H_lin:.2f} m")
    print(f"    measured at {series[-1][0]:.0f} s over {lo}..{hi} m: L {L:.1f} m ({L / L_lin:.3f}), "
          f"H {H:.2f} m ({H / H_lin:.3f}) over {n} waves; the window's level {drift[-1]:+.2f} m "
          f"(every 10 s: {drift}); plate at {sim.piston_x:.2f} m; {len(sim.pos):,} particles, {wall:.0f} s  "
          + ("OK" if ok else "FAIL"), flush=True)
    return ok, {"grid": [nx, ny], "h": h, "S": S, "T": T, "plate": [PLATE_X0, PLATE_W, h + 8.5],
                "sim_kw": SIM_KW, "k_linear": round(k, 5), "L_linear": round(L_lin, 2),
                "H_over_S_linear": round(float(transfer(k, h)), 3), "H_linear": round(float(H_lin), 3),
                "window_m": [lo, hi], "L_measured": round(L, 2), "H_measured": round(H, 3),
                "waves_measured": n, "L_ratio": round(L / L_lin, 4), "H_ratio": round(H / H_lin, 4),
                "window_level_every_10s_m": drift, "wall_s": round(wall, 1),
                "surface_every_5m": {f"{t:.0f}s": [round(float(v), 3) for v in e[::5]] for t, e in series}}


def check_fast():
    print("=== fast: a 4 m stroke at 4 s, nothing leaks through the plate", flush=True)
    nx, ny, h, S, T, top = 200, 60, 20.0, 4.0, 4.0, 50.0
    front = PLATE_X0 + PLATE_W
    sim = Sim(nx, ny, piston=(S, T, PLATE_X0, PLATE_W, top))
    sim.seed(lambda x, y: (x > front) & (y < h))
    sim.calibrate()
    n_liq = int((sim.typ == 1).sum())
    worst = {"inside": 0, "behind": 0, "speed": 0.0, "air_speed": 0.0, "lost": 0}
    t0 = time.time()

    def cb(f, s):
        progress(s, f, 0.1, t0)
        liq = s.typ == 1
        p = s.pos[liq]
        inside = int((s._phi_at(p, s._plate_phi)[0] < 0.0).sum())
        behind = int(((p[:, 0] < s.piston_x) & (p[:, 1] < top)).sum())
        worst["inside"] = max(worst["inside"], inside)
        worst["behind"] = max(worst["behind"], behind)
        worst["speed"] = max(worst["speed"], float(np.linalg.norm(s.vel[liq], axis=1).max()))
        worst["air_speed"] = max(worst["air_speed"], float(np.linalg.norm(s.vel[~liq], axis=1).max()))
        worst["lost"] = max(worst["lost"], abs(int(liq.sum()) - n_liq))
    sim.run(120, 0.1, on_frame=cb)
    ok = worst["inside"] == 0 and worst["behind"] == 0 and worst["speed"] < 30.0 and worst["lost"] == 0
    print(f"    12 s, worst frame: liquid inside the plate {worst['inside']}, behind it {worst['behind']}, "
          f"fastest liquid {worst['speed']:.1f} m/s (air {worst['air_speed']:.1f}), liquid count change {worst['lost']}; "
          f"{n_liq:,} liquid particles, {time.time() - t0:.0f} s  " + ("OK" if ok else "FAIL"), flush=True)
    return ok, dict(worst, grid=[nx, ny], S=S, T=T, top=top, liquid=n_liq)


def main(argv):
    global SIM_KW, SUFFIX
    args = list(argv)
    if "--alpha" in args:
        i = args.index("--alpha")
        a = float(args[i + 1])
        del args[i:i + 2]
        SIM_KW = dict(alpha_l=a, alpha_g=a)
        SUFFIX = f"_alpha{a:g}"
    which = args[0] if args else "all"
    os.makedirs(RESULTS, exist_ok=True)
    results, ok = {}, True
    for name, fn in (("rest", check_rest), ("dispersion", check_dispersion), ("fast", check_fast)):
        if which in ("all", name):
            good, rec = fn()
            ok &= good
            results[name] = rec
    path = os.path.join(RESULTS, f"receipt{SUFFIX}.json" if which == "all" else f"receipt_{which}{SUFFIX}.json")
    with open(path, "w") as f:
        json.dump(results, f, indent=1)
    print("PISTON:", "PASS" if ok else "FAIL", flush=True)
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main(sys.argv[1:])
