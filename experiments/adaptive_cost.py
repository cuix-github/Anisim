"""Milestone 6.5a stage 3c: what the adaptive pressure solve costs at finale
size, uniform against naive against cleanup.

INDICATIVE, not the canonical receipt: the canonical profiling state (the
mid-bore dump in out/loop3_D_finale_c1/) is the bore with two pillars, and
obstacles still send the adaptive solve back to the uniform one. This is a
synthetic state at the same grid, 256x64x128 (2.1M cells): a dam with a
shallow bed, 1000:1, a divergence-free swirl in it, liquid particles seeded
8 per cell (the tree reads only liquid counts), one step_3d step through
solve_step3_test.json. PFFLIP_PROFILE=1 times the node's phases from inside
it, so Bifrost's startup is not in the numbers; the pressure phase is
assembly + solve + correction (the adaptive branch's tree, hierarchy, both
solves and local fills all fall in it). Each configuration runs REPEATS
times, serially, nothing else running; the median is reported with the
spread.

From the repo root:  .venv\\Scripts\\python.exe experiments\\adaptive_cost.py
"""
import json
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
sys.path[:0] = [os.path.join(ROOT, "gates"), os.path.join(ROOT, "graphs")]
from loop3d import seed3d  # noqa: E402

NX, NY, NZ, DT = 256, 64, 128, 0.02
RHO_L, RHO_G = 1000.0, 1.0
REPEATS = 3
OUT = os.path.join(ROOT, "out", "adaptive_cost")
RESULTS = os.path.join(ROOT, "results", "3d", "adaptive")
NAMES = ("u_mass", "u_mom", "u_phase", "v_mass", "v_mom", "v_phase", "w_mass", "w_mom",
         "w_phase", "particle_phase", "positions", "velocities")
OUTS = ("out_u", "out_v", "out_w", "pressure", "out_positions", "out_velocities", "out_escaped")
PROF = re.compile(r"PFPROF .*iters=(\d+) asm=([\d.]+) solve=([\d.]+) corr=([\d.]+) "
                  r"part=([\d.]+) total=([\d.]+)")
SUB = re.compile(r"PFPROF_ADAPTIVE cells=(\d+) of \d+ levels=\d+ iters=(\d+)\+(\d+) tree=([\d.]+) "
                 r"hier=([\d.]+) solve=([\d.]+) grad=([\d.]+) cleanup=([\d.]+) fill=([\d.]+)")


def mask(x, y, z):
    return (y < 8) | ((x < 64) & (y < 48))


def swirl(x, y, z):
    """Divergence-free, no normal flow through the walls, a few cells per second."""
    A, B = 100.0, 60.0
    u = A * np.pi / NY * np.sin(np.pi * x / NX) * np.cos(np.pi * y / NY)
    v = (-A * np.pi / NX * np.cos(np.pi * x / NX) * np.sin(np.pi * y / NY)
         + B * np.pi / NZ * np.sin(np.pi * y / NY) * np.cos(np.pi * z / NZ))
    w = -B * np.pi / NY * np.cos(np.pi * y / NY) * np.sin(np.pi * z / NZ)
    return u, v, w


def package(d):
    os.makedirs(d, exist_ok=True)
    if os.path.exists(os.path.join(d, "positions.npy")):
        return
    offs = {"u": (0.0, 0.5, 0.5), "v": (0.5, 0.0, 0.5), "w": (0.5, 0.5, 0.0)}
    dims = {"u": (NX + 1, NY, NZ), "v": (NX, NY + 1, NZ), "w": (NX, NY, NZ + 1)}
    for c in "uvw":
        x, y, z = np.meshgrid(*[np.arange(n) + o for n, o in zip(dims[c], offs[c])], indexing="ij")
        phi = mask(x, y, z).astype(float)
        mass = RHO_G + (RHO_L - RHO_G) * phi
        vel = swirl(x, y, z)["uvw".index(c)]
        np.save(os.path.join(d, f"{c}_mass.npy"), mass.ravel().astype(np.float32))
        np.save(os.path.join(d, f"{c}_mom.npy"), (mass * vel).ravel().astype(np.float32))
        np.save(os.path.join(d, f"{c}_phase.npy"), phi.ravel().astype(np.float32))
    pos, typ = seed3d(NX, NY, NZ, mask)
    keep = typ == 1
    pos = pos[keep]
    np.save(os.path.join(d, "positions.npy"), pos.astype(np.float32))
    np.save(os.path.join(d, "velocities.npy"),
            np.stack(swirl(pos[:, 0], pos[:, 1], pos[:, 2]), 1).astype(np.float32))
    np.save(os.path.join(d, "particle_phase.npy"), np.ones(len(pos), np.float32))


def run(indir, outdir, **ports):
    os.makedirs(outdir, exist_ok=True)
    args = [os.path.join(ROOT, "graphs", "run_solve_step3.bat")]
    for k, v in ports.items():
        args += ["--set-port", k, str(v)]
    for n in NAMES:
        args += ["--set-port", "path_" + n, os.path.join(indir, n + ".npy").replace("\\", "/")]
    for n in OUTS:
        args += ["--set-port", "path_" + n, os.path.join(outdir, n + ".npy").replace("\\", "/")]
    env = dict(os.environ, PFFLIP_PROFILE="1")
    res = subprocess.run(args, capture_output=True, text=True, env=env)
    m = PROF.search(res.stdout + res.stderr)
    if res.returncode != 0 or m is None:
        sys.stderr.write(res.stdout[-2500:] + res.stderr[-1500:])
        raise SystemExit("bifcmd step_3d failed")
    it, asm, solve, corr, part, total = (float(g) for g in m.groups())
    row = {"iterations": int(it), "pressure_ms": asm + solve + corr, "asm_ms": asm,
           "solve_ms": solve, "corr_ms": corr, "particles_ms": part, "node_ms": total}
    a = SUB.search(res.stdout + res.stderr)
    if a:   # the adaptive solve's own parts, inside solve_ms
        g = [float(x) for x in a.groups()]
        row["adaptive"] = {"cells": int(g[0]), "iterations": [int(g[1]), int(g[2])],
                           "tree_ms": g[3], "hierarchy_ms": g[4], "solve_ms": g[5],
                           "gradient_ms": g[6], "cleanup_ms": g[7], "fill_ms": g[8]}
    return row


def main():
    d = os.path.join(OUT, "state_in")
    package(d)
    n = len(np.load(os.path.join(d, "positions.npy"), mmap_mode="r"))
    print(f"{NX}x{NY}x{NZ} ({NX * NY * NZ / 1e6:.1f}M cells), {n / 1e6:.1f}M liquid particles",
          flush=True)
    rows = {}
    for tol in (1e-6, 1e-4):
        for name, solver in (("uniform", 0), ("naive", 1), ("cleanup", 2)):
            reps = [run(d, os.path.join(OUT, f"{name}_tol{tol:g}_r{r}"), nx=NX, ny=NY, nz=NZ, dt=DT,
                        pressure_solver=solver, preconditioner=1, tolerance=tol)
                    for r in range(REPEATS)]
            p = sorted(r["pressure_ms"] for r in reps)
            row = {"iterations": reps[0]["iterations"],
                   "pressure_ms_median": p[len(p) // 2], "pressure_ms_min": p[0],
                   "pressure_ms_max": p[-1],
                   "node_ms_median": sorted(r["node_ms"] for r in reps)[len(p) // 2]}
            if "adaptive" in reps[0]:
                mid = sorted(reps, key=lambda r: r["pressure_ms"])[len(p) // 2]
                row["adaptive_parts_of_median_run"] = mid["adaptive"]
            rows[f"{name}_tol{tol:g}"] = row
            print(f"  tol {tol:g} {name}: {row['iterations']} iterations, pressure phase "
                  f"{row['pressure_ms_median']:.0f} ms (min {p[0]:.0f}, max {p[-1]:.0f}); "
                  f"whole node {row['node_ms_median']:.0f} ms", flush=True)
            if "adaptive_parts_of_median_run" in row:
                a = row["adaptive_parts_of_median_run"]
                print(f"      {a['cells']} adaptive cells; tree {a['tree_ms']:.0f}, hierarchy "
                      f"{a['hierarchy_ms']:.0f}, solve {a['solve_ms']:.0f}, gradient "
                      f"{a['gradient_ms']:.0f}, cleanup {a['cleanup_ms']:.0f}, fill "
                      f"{a['fill_ms']:.0f} ms", flush=True)
    os.makedirs(RESULTS, exist_ok=True)
    with open(os.path.join(RESULTS, "cost_indicative.json"), "w", newline="\n") as f:
        json.dump({"state": "synthetic dam + bed + swirl, no obstacles (indicative; the canonical "
                            "state has pillars the adaptive solve cannot take yet)",
                   "grid": [NX, NY, NZ], "liquid_particles": n, "repeats": REPEATS,
                   "timer": "PFFLIP_PROFILE inside step_3d; pressure = asm + solve + corr; "
                            "the adaptive parts (PFPROF_ADAPTIVE) sit inside solve",
                   "rows": rows}, f, indent=1)


if __name__ == "__main__":
    main()
