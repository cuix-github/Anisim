"""What a level-set solid costs step_3d: one step on the same state with no solid, the
box path, and the same box as a field (binary faces and cut cells), timed per phase
through PFFLIP_PROFILE (assembly, solve, correction, particle pass). Each arm repeats
and reports its median; the field arms pay for the open fractions in the assembly and
for one distance lookup per particle in the particle pass.

Run from the repo root:
  .venv\\Scripts\\python.exe experiments\\level_set_cost.py [nx ny nz] [repeats]
Set OMP_NUM_THREADS to fix the worker count; every arm runs with the same one.
"""
import os
import re
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "gates"))
from validate_level_set_node import (OUT, flow_inputs, node_grid3, save,  # noqa: E402
                                     sdf_box3)

ARGS = [int(a) for a in sys.argv[1:4]] if len(sys.argv) >= 4 else [96, 48, 64]
REPEATS = int(sys.argv[4]) if len(sys.argv) >= 5 else 3


def main():
    nx, ny, nz = ARGS
    box = (0.40 * nx, 0.0, 0.40 * nz, 0.55 * nx, 0.7 * ny, 0.60 * nz)
    box = tuple(float(round(v)) for v in box)
    phi = sdf_box3(*node_grid3(nx, ny, nz), box)
    base = flow_inputs(nx, ny, nz, phi, level=0.4 * ny)
    print(f"{nx}x{ny}x{nz}, {len(base['positions'])} particles, a pillar {box};"
          f" OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', 'unset')}, median of {REPEATS}")
    arms = {"no solid": ({}, {}),
            "box path": ({"obstacles_min": np.array([box[:3]], np.float32),
                          "obstacles_max": np.array([box[3:]], np.float32)}, {}),
            "field, binary": ({"solid_phi": phi.astype(np.float32).ravel()},
                              {"solid_mode": 1}),
            "field, cut cells": ({"solid_phi": phi.astype(np.float32).ravel()},
                                 {"solid_mode": 2})}
    bat = os.path.join(HERE, "..", "..", "graphs", "run_solve_step3.bat")
    env = dict(os.environ, PFFLIP_PROFILE="1")
    print(f"  {'arm':18s} {'assembly':>9s} {'solve':>8s} {'correct':>8s} {'particles':>10s}"
          f" {'total':>8s}  iterations (ms)")
    for name, (extra, ports) in arms.items():
        d_in = os.path.join(OUT, "cost_in_" + name.replace(" ", "_").replace(",", ""))
        d_out = d_in.replace("cost_in_", "cost_out_")
        for k, v in dict(base, **extra).items():
            save(d_in, k, v)
        names = ("u_mass", "u_mom", "u_phase", "v_mass", "v_mom", "v_phase", "w_mass",
                 "w_mom", "w_phase", "particle_phase", "positions", "velocities") + tuple(extra)
        args = [bat]
        for k, v in dict({"nx": nx, "ny": ny, "nz": nz, "dt": 0.02, "preconditioner": 1},
                         **ports).items():
            args += ["--set-port", k, str(v)]
        for n in names:
            args += ["--set-port", "path_" + n, os.path.join(d_in, n + ".npy").replace("\\", "/")]
        args += ["--set-port", "path_out_u", os.path.join(d_out, "out_u.npy").replace("\\", "/")]
        rows = []
        for _ in range(REPEATS):
            res = subprocess.run(args, capture_output=True, text=True, env=env)
            m = re.search(r"iters=(\d+) asm=([\d.]+) solve=([\d.]+) corr=([\d.]+) "
                          r"part=([\d.]+) total=([\d.]+)", res.stderr + res.stdout)
            if res.returncode != 0 or not m:
                sys.stderr.write(res.stdout[-2000:] + res.stderr[-2000:])
                raise SystemExit(f"step failed for {name}")
            rows.append([float(g) for g in m.groups()])
        med = np.median(np.array(rows), 0)
        print(f"  {name:18s} {med[1]:9.1f} {med[2]:8.1f} {med[3]:8.1f} {med[4]:10.1f}"
              f" {med[5]:8.1f}  {int(med[0])}")


if __name__ == "__main__":
    main()
