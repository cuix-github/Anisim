"""The 3D loop runner: seeds the particle state, saves it as .npy, and runs the
whole simulation headlessly through run_sim_3d.bat (or run_sim_3d_fused.bat) via
bifcmd. Shared by the loop gates and the 3D experiments; experiments/demo3d.py
re-exports these names.
"""
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))


def seed3d(nx, ny, nz, mask, jitter=0.35, seed=1):
    rng = np.random.default_rng(seed)
    cx, cy, cz = np.meshgrid(np.arange(nx), np.arange(ny), np.arange(nz), indexing="ij")
    ps, ts = [], []
    for ox in (0.25, 0.75):
        for oy in (0.25, 0.75):
            for oz in (0.25, 0.75):
                p = np.stack([cx + ox, cy + oy, cz + oz], -1).reshape(-1, 3).astype(np.float64)
                p += rng.uniform(-jitter, jitter, p.shape) * 0.5
                ps.append(p)
                ts.append(mask(p[:, 0], p[:, 1], p[:, 2]).astype(np.int32))
    return np.concatenate(ps), np.concatenate(ts)


def lattice(ni, nj, nk):
    i, j, k = np.meshgrid(np.arange(ni) + 0.5, np.arange(nj) + 0.5, np.arange(nk) + 0.5,
                          indexing="ij")
    return np.stack([i.ravel(), j.ravel(), k.ravel()], 1).astype(np.float32)


FINAL_STATE = ("velocities", "mass", "momx", "momy", "momz", "tau", "wt", "wtph",
               "scale", "phase")
# per-substep dumps beside frame.#### (always written): escaped flags, particle
# scales, phases; "vel" (particle velocities) is opt-in
DUMPS = ("esc", "scl", "ph")


def load_final_state(d):
    """The full particle state a run left in d (final_state=True), as init_arrays."""
    st = {"positions": np.load(os.path.join(d, "final.npy"))}
    for nm in FINAL_STATE:
        st["phase_st" if nm == "phase" else nm] = np.load(os.path.join(d, f"final_{nm}.npy"))
    st["particle_phase"] = st["phase_st"]
    return st


def run_bifrost3(tag, nx, ny, nz, mask, substeps, dt, rho_l=1000.0, rho_g=1.0,
                 extra_ports=None, seed_fn=None, obstacles=None, fused=False,
                 dump_every=1, init_arrays=None, final_state=False, step_offset=0,
                 dumps=DUMPS, solid=None):
    """init_arrays resumes from a saved state instead of seeding; final_state also
    writes the whole state at the end (final_<name>.npy beside final.npy); a
    resumed run passes step_offset = substeps already run, so the ST jitter
    continues its stream rather than restarting it. dumps names the
    per-substep dumps written beside frame.#### (<name>.####). solid is a
    level-set solid, (distance at the (nx+1, ny+1, nz+1) grid nodes, negative
    inside; solid_mode 1 binary faces or 2 cut cells); seeding stays the
    caller's (seed_fn), so particles inside it are the caller's to drop."""
    d = os.path.join(HERE, "..", "out", "loop3_" + tag)
    os.makedirs(d, exist_ok=True)
    if init_arrays is not None:
        pos, typ = init_arrays["positions"], init_arrays["particle_phase"].astype(np.int32)
    else:
        pos, typ = (seed_fn or seed3d)(nx, ny, nz, mask)
    mass = np.where(typ == 1, rho_l, rho_g).astype(np.float32)
    zeros = np.zeros(len(pos), np.float32)
    arrs = {"positions": pos.astype(np.float32), "velocities": np.zeros_like(pos, dtype=np.float32),
            "particle_phase": typ.astype(np.float32), "mass": mass,
            "momx": zeros, "momy": zeros, "momz": zeros,
            "tau": zeros, "wt": np.ones(len(pos), np.float32),
            "wtph": typ.astype(np.float32),
            "scale": np.ones(len(pos), np.float32),
            "phase_st": typ.astype(np.float32),
            "obs_min": (obstacles[0].astype(np.float32) if obstacles is not None and len(obstacles[0])
                        else np.array([[-100.0, -100.0, -100.0]], np.float32)),
            "obs_max": (obstacles[1].astype(np.float32) if obstacles is not None and len(obstacles[1])
                        else np.array([[-99.0, -99.0, -99.0]], np.float32)),
            "probes_u": lattice(nx + 1, ny, nz), "probes_v": lattice(nx, ny + 1, nz),
            "probes_w": lattice(nx, ny, nz + 1),
            # no solid: one value, which the node ignores at solid_mode 0
            "solid_phi": (np.asarray(solid[0], np.float32).ravel() if solid is not None
                          else np.ones(1, np.float32))}
    if solid is not None:
        assert np.shape(solid[0]) == (nx + 1, ny + 1, nz + 1), np.shape(solid[0])
        extra_ports = dict(extra_ports or {}, solid_mode=int(solid[1]))
    if init_arrays is not None:
        for k2 in ("velocities", "mass", "momx", "momy", "momz", "tau", "wt", "wtph",
                   "scale", "phase_st"):
            arrs[k2] = init_arrays[k2].astype(np.float32)
    for k, v in arrs.items():
        np.save(os.path.join(d, k + ".npy"), v)
    args = [os.path.join(HERE, "..", "graphs", "run_sim_3d_fused.bat" if fused else "run_sim_3d.bat"),
            "--set-port", "nx", str(nx), "--set-port", "ny", str(ny),
            "--set-port", "nz", str(nz), "--set-port", "dt", str(dt),
            "--set-port", "substeps", str(substeps), "--set-port", "preconditioner", "1",
            "--set-port", "dump_every", str(dump_every),
            "--set-port", "rho_liquid", str(rho_l), "--set-port", "rho_air", str(rho_g),
            "--set-port", "pos_pattern", (d + "/frame.####").replace("\\", "/"),
            "--set-port", "diag_pattern", (d + "/diag.####").replace("\\", "/"),
            "--set-port", "path_final_positions", (d + "/final.npy").replace("\\", "/"),
            *[x for nm in (FINAL_STATE if final_state else ())
              for x in ("--set-port", "path_final_" + nm,
                        (d + "/final_" + nm + ".npy").replace("\\", "/"))],
            *(("--set-port", "step_offset", str(step_offset)) if step_offset else ()),
            *[x for k2, v2 in (extra_ports or {}).items()
              for x in ("--set-port", k2, str(v2))],
            *[x for nm in dumps
              for x in ("--set-port", nm + "_pattern", (d + f"/{nm}.####").replace("\\", "/"))]]
    for k in arrs:
        args += ["--set-port", "path_" + k, os.path.join(d, k + ".npy").replace("\\", "/")]
    res = subprocess.run(args, capture_output=True, text=True)
    last = ((substeps - 1) // dump_every) * dump_every
    if res.returncode != 0 or not os.path.exists(os.path.join(d, f"frame.{last:04d}.npy")):
        sys.stderr.write(res.stdout[-3000:] + res.stderr[-2000:])
        raise SystemExit(f"bifcmd 3D loop failed for {tag}")
    diags = np.stack([np.load(os.path.join(d, f"diag.{k:04d}.npy")) for k in range(substeps)])
    # the node refuses bad inputs (a solid_phi of the wrong size or dtype) with
    # iterations_used = -1, empty outputs and a "step_3d: error:" line on stderr
    if (diags[:, 0] < 0).any():
        errs = [ln for ln in res.stderr.splitlines() if ln.startswith("step_3d: error:")]
        raise SystemExit(f"step_3d refused its inputs for {tag}: "
                         + (errs[0] if errs else "see diag"))
    return d, typ, diags
