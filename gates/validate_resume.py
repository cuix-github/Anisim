"""Chunked resume: a run split into two launches, the second resumed from the
first's full-state dump, must reproduce the uninterrupted run bit for bit. This
is what lets a long 3D stage run as a chain of launches under the background
runner's time cap (experiments/run_bore_suite.py's stage/finalize mode).

Production config (ST, escape, adaptive particles) on a small wet-bed bore with
one obstacle, through both 3D graphs: 2K substeps in one launch against K + K
resumed with step_offset = K. Control: the same resume with step_offset = 0
replays the first launch's ST jitter stream, so it must NOT match - proof the
check can fail and that the offset is what makes the resume exact.

Run from gates:  ..\\.venv\\Scripts\\python.exe validate_resume.py
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
sys.path.insert(0, os.path.join(HERE, "..", "reference2d"))
from pfflip2d import Sim  # noqa: E402
from loop3d import FINAL_STATE, load_final_state, run_bifrost3  # noqa: E402

NX, NY, NZ, DT, K = 40, 20, 16, 0.12, 4
PMIN = np.array([[22.0, 0.0, 6.0]], np.float32)
PMAX = np.array([[25.0, 10.0, 10.0]], np.float32)


def mask(x, y, z):
    return (y < 3) | ((x < 10) & (y < 15))


def seed(nx, ny, nz, m):
    from loop3d import seed3d
    pos, typ = seed3d(nx, ny, nz, m)
    a, b = PMIN[0], PMAX[0]
    keep = ~((pos[:, 0] > a[0]) & (pos[:, 0] < b[0]) & (pos[:, 1] > a[1]) &
             (pos[:, 1] < b[1]) & (pos[:, 2] > a[2]) & (pos[:, 2] < b[2]))
    return pos[keep], typ[keep]


def state(d):
    s = {"final": np.load(os.path.join(d, "final.npy"))}
    for nm in FINAL_STATE:
        s[nm] = np.load(os.path.join(d, f"final_{nm}.npy"))
    return s


def maxdiff(a, b):
    if a.shape != b.shape:
        return float("inf")
    return float(np.abs(a.astype(np.float64) - b.astype(np.float64)).max())


def main():
    cal = Sim(64, 64, adapt=True, st=True)
    cal.seed(lambda x, y: np.ones_like(x, dtype=bool))
    cal.calibrate()
    ports = {"st": 1, "escape": 1, "adapt": 1, "rho0_face": cal.rho0_face,
             "coarse_gain": cal.coarse_gain, "ws_epsilon": 16}
    common = dict(extra_ports=ports, seed_fn=seed, obstacles=(PMIN, PMAX),
                  final_state=True)
    ok = True
    for fused in (False, True):
        g = "fused" if fused else "stock"
        d_full, _, diag_full = run_bifrost3(f"resume_{g}_full", NX, NY, NZ, mask, 2 * K,
                                            DT, fused=fused, **common)
        d_c0, _, _ = run_bifrost3(f"resume_{g}_c0", NX, NY, NZ, mask, K, DT,
                                  fused=fused, **common)
        init = load_final_state(d_c0)
        d_c1, _, diag_c1 = run_bifrost3(f"resume_{g}_c1", NX, NY, NZ, mask, K, DT,
                                        fused=fused, init_arrays=init, step_offset=K,
                                        **common)
        d_ctl, _, _ = run_bifrost3(f"resume_{g}_c1_nooffset", NX, NY, NZ, mask, K, DT,
                                   fused=fused, init_arrays=init, **common)

        full, c1, ctl = state(d_full), state(d_c1), state(d_ctl)
        worst = max(maxdiff(full[k], c1[k]) for k in full)
        frames = max(
            maxdiff(np.load(os.path.join(d_full, f"{p}.{K + k:04d}.npy")),
                    np.load(os.path.join(d_c1, f"{p}.{k:04d}.npy")))
            for k in range(K) for p in ("frame", "ph", "scl"))
        diags = maxdiff(diag_full[K:], diag_c1)
        ctl_d = max(maxdiff(full[k], ctl[k]) for k in full)
        n0 = len(np.load(os.path.join(d_full, "positions.npy")))
        print(f"{g}: {n0} particles seeded, {len(full['final'])} after {2 * K} substeps")
        print(f"  resumed vs uninterrupted: final state max |diff| {worst:.3e}, "
              f"per-substep dumps {frames:.3e}, diagnostics {diags:.3e}")
        print(f"  control (no offset): final state max |diff| {ctl_d:.3e}")
        g_ok = worst == 0.0 and frames == 0.0 and diags == 0.0 and ctl_d > 0.0
        print(f"  {g}: {'PASS' if g_ok else 'FAIL'}")
        ok &= g_ok
    print("CHUNKED RESUME:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
