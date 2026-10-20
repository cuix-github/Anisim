"""Bake an independent spectral body-acceleration source, then run Bifrost.

No target velocities or prescribed surface. The graph reads acceleration
fields per absolute substep and projects them with the actual particle phase.
The Gaussian source amplitude uses a shallow-water linear transfer estimate;
measured response, not that estimate, determines the produced wave height.
"""
import argparse
import json
import math
import os
from pathlib import Path
import sys
import time
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "graphs"), str(ROOT / "experiments")]
from loop3d import run_bifrost3, load_final_state, DUMPS
import wave_basin as basin
from solitary_beach import calibration
from wave_features import atomic_json, spray_ports


def recipe(single=False, flat=False):
    if single:
        return np.array([[1., 6., 0., 0.]])
    rng = np.random.default_rng(20)
    leads = np.array([[4., 8., 15., 0.], [3., 9.5, -18., 0.], [1.5, 6.5, 5., 0.]])
    freq = np.linspace(1 / 7.5, 1 / 3.0, 997)
    weight = freq ** -2.5 * rng.rayleigh(1, len(freq))
    heights = weight * (1.5 / 4) / np.sqrt(np.sum(weight ** 2 / 8))
    background = np.column_stack((heights, 1 / freq, rng.normal(0, 18, len(freq)),
                                  rng.uniform(0, 2 * np.pi, len(freq))))
    background[:, 2] = np.clip(background[:, 2], -35, 35)
    rows = np.vstack((leads, background))
    if flat:
        rows[:, 2] = 0
    for row in rows[:3]:
        k = basin._k(row[1], 12)
        theta = np.radians(row[2])
        row[3] = 2 * np.pi / row[1] * 18 - k * (62 * np.cos(theta) + 36 * np.sin(theta))
    return rows


def bake(root, grid, dx, dt, steps, rows, gain):
    d = root / "inputs/forces"
    d.mkdir(parents=True)
    nx, ny, nz = grid
    omega = 2 * np.pi / rows[:, 1]
    k = np.array([basin._k(T, 12) for T in rows[:, 1]])
    theta = np.radians(rows[:, 2])
    sigma = np.clip(.12 * 2 * np.pi / k, 2, 6)
    kx, kz = k * np.cos(theta), k * np.sin(theta)
    amplitude = gain * 9.8 * rows[:, 0] * np.cos(theta) / (
        np.sqrt(2 * np.pi) * sigma * np.exp(-.5 * (kx * sigma) ** 2))
    times = (np.arange(steps) + .5) * dt
    ramp = np.minimum(times / 4, 1)
    ramp = ramp ** 2 * (3 - 2 * ramp)
    tc, ts = np.cos(times[:, None] * omega), np.sin(times[:, None] * omega)
    for axis in "uw":
        if axis == "w" and nz == 1:
            continue
        x = (np.arange(nx + (axis == "u")) + (0 if axis == "u" else .5)) * dx
        z = (np.arange(nz + (axis == "w")) + (0 if axis == "w" else .5)) * dx
        active = np.flatnonzero(np.abs(x - 38) <= 24)
        shape = (len(x), ny, len(z))
        xx, zz = np.meshgrid(x[active], z, indexing="ij")
        envelope = amplitude[:, None] * np.exp(-.5 * ((xx.ravel()[None, :] - 38) / sigma[:, None]) ** 2)
        envelope *= (np.cos(theta) if axis == "u" else np.sin(theta))[:, None]
        phase = kz[:, None] * zz.ravel()[None, :] + rows[:, 3, None]
        # Expensive spectrum sum is baked once, outside the solver. No water state
        # is read here and there is no feedback toward a prescribed velocity.
        field = (tc @ (envelope * np.cos(phase)) + ts @ (envelope * np.sin(phase))) / dx
        for step in range(steps):
            a = np.zeros(shape, np.float32)
            a[active, :, :] = (field[step].reshape(len(active), len(z)) * ramp[step])[:, None, :]
            np.save(d / f"a_{axis}.{step:04d}.npy", a.ravel())
    np.save(d / "zero.npy", np.zeros(1, np.float32))
    ports = {"external_forces": 1}
    for axis in "uvw":
        dynamic = axis == "u" or (axis == "w" and nz > 1)
        ports["path_liquid_acceleration_" + axis] = str(d / (f"a_{axis}.####" if dynamic else "zero.npy")).replace("\\", "/")
        ni, nj, nk = nx + (axis == "u"), ny + (axis == "v"), nz + (axis == "w")
        x = (np.arange(ni) + (0 if axis == "u" else .5)) * dx
        z = (np.arange(nk) + (0 if axis == "w" else .5)) * dx
        rate = np.maximum(0, (20 - x) / 20) ** 2
        side = np.zeros(nk)
        if nz > 1:
            side = .6 * (np.maximum(0, (6 - z) / 6) ** 2 +
                         np.maximum(0, (z - (nz * dx - 6)) / 6) ** 2)
        drag = np.broadcast_to(rate[:, None, None] + side[None, None, :], (ni, nj, nk)).astype(np.float32)
        np.save(d / f"drag_{axis}.npy", drag.ravel())
        ports["path_linear_drag_" + axis] = str(d / f"drag_{axis}.npy").replace("\\", "/")
    return ports


def export_frame(d, step, dest, frame, dx, record):
    names = dict(pos="frame", phase="ph", escape="esc", spray="spr", radius="rad", vel="vel")
    a = {key: np.load(d / f"{prefix}.{step:04d}.npy", mmap_mode="r") for key, prefix in names.items()}
    if any(not np.isfinite(v).all() for v in a.values()):
        raise ValueError("Nonfinite particle output")
    water, esc, spr = a["phase"] >= .5, a["escape"] > .5, a["spray"] > .5
    masks = dict(water=water & ~esc & ~spr, droplets=water & esc & ~spr, spray=spr, air=~water)
    counts = {}
    for layer, mask in masks.items():
        ix = np.flatnonzero(mask)
        counts[layer] = len(ix)
        if layer in ("water", "air"):
            cap = 120_000 if layer == "water" else 10_000
            ix = ix[::max(1, math.ceil(len(ix) / cap))]
        ld = dest / layer
        ld.mkdir(parents=True, exist_ok=True)
        np.save(ld / f"liquid.{frame:04d}.npy", np.asarray(a["pos"][ix] * dx, np.float32))
        np.save(ld / f"speed.{frame:04d}.npy", np.asarray(np.linalg.norm(a["vel"][ix], axis=1) * dx, np.float32))
    interior = water & (a["pos"][:, 0] * dx > 68)
    pos = a["pos"][interior] * dx
    record.update(frame=frame, counts=counts, escaped_beyond_source=int(np.count_nonzero(interior & esc)),
                  spray_beyond_source=int(np.count_nonzero(interior & spr)),
                  liquid_max_y_m=float(pos[:, 1].max()) if len(pos) else 0)


def run(args):
    root = (ROOT / "out/runs" / args.tag).resolve()
    root.mkdir(parents=True, exist_ok=False)
    (root / "receipts").mkdir()
    (root / "logs").mkdir()
    start = time.monotonic()
    dx, dt = args.dx, 1 / 24
    length, height, width = (260, 32, dx) if args.single else (180, 32, dx if args.kind == "2d" else 72)
    grid = tuple(round(v / dx) for v in (length, height, width))
    rows = recipe(args.single, grid[2] == 1)
    steps = round(args.seconds / dt)
    force_ports = bake(root, grid, dx, dt, steps, rows, args.gain)
    run_spec = basin.Run(1 / dx)
    run_spec.grid = grid
    sc = dict(m=1 / dx, h=12 / dx, toe=(300 if args.single else 70) / dx, slope=.2)
    ports = dict(calibration(), **spray_ports(1 / dx), gravity=-9.8 / dx,
                 zone_length=0, **force_ports)
    if args.single:
        ports.update(escape=0, st=0, adapt=0, droplet_model=0, emit_model=0)
    plan = dict(tag=args.tag, grid=grid, domain_m=[n * dx for n in grid], dx=dx,
                dt=dt, seconds=args.seconds, ports=ports, force_paths=force_ports,
                source_center_m=38, source_support_m=[14, 62], beach_toe_m=sc["toe"] * dx,
                beach_slope=.2, still_depth_m=12, rows=rows.tolist(), gain=args.gain,
                interpretation="Phase-weighted body acceleration before projection; no velocity target",
                display_fps=12, adaptive_pressure="Uniform because the solid beach is present")
    atomic_json(root / "inputs/plan.json", plan)
    init = basin.seed(run_spec, sc)
    if grid[2] == 1:
        init["positions"][:, 2] = .5
    solid = (basin.beach_phi(run_spec, sc), 2)
    records, chunks, frame = [], [], 0
    full = root / "cache/full"
    chunk_steps = 120
    for chunk, offset in enumerate(range(0, steps, chunk_steps)):
        n = min(chunk_steps, steps - offset)
        d = root / "dumps" / f"c{chunk}"
        t0 = time.monotonic()
        _, _, diag = run_bifrost3(args.tag, *grid, None, n, dt, extra_ports=ports,
            init_arrays=init, release_init=True, final_state=True, step_offset=offset,
            fused=True, dump_every=2, dumps=DUMPS + ("vel",), solid=solid,
            out_dir=str(d), log_prefix=str(root / "logs" / f"c{chunk}"))
        elapsed = time.monotonic() - t0
        if not np.isfinite(diag).all() or np.max(diag[:, 1]) > 1.05e-6:
            raise RuntimeError("Pressure convergence/finite check failed")
        for s in range(0, n, 2):
            frame += 1
            rec = dict(time_s=(offset + s + 1) * dt)
            export_frame(d, s, full, frame, dx, rec)
            records.append(rec)
        chunks.append(dict(chunk=chunk, end_s=(offset + n) * dt, wall_s=elapsed,
                           iterations_mean=float(diag[:, 0].mean()), residual_max=float(diag[:, 1].max()),
                           max_speed_m_s=float(diag[:, 3].max()) * dx))
        atomic_json(root / "receipts/progress.json", dict(chunks=chunks, frames=records, elapsed_s=time.monotonic() - start))
        print(json.dumps(chunks[-1]), flush=True)
        init = load_final_state(str(d))
    atomic_json(full / "cache.json", dict(frames=frame, fps=12, units="metres",
        times=[r["time_s"] for r in records], grid=grid, dt=dt, step=2,
        bulk_display_cap=120000, escaped_and_spray_not_decimated=True))
    best = max(records, key=lambda r: r["escaped_beyond_source"] + r["spray_beyond_source"])
    atomic_json(root / "receipts/particle_audit.json", dict(best_simultaneous_frame=best,
        peak_escaped_beyond_source=max(r["escaped_beyond_source"] for r in records),
        peak_spray_beyond_source=max(r["spray_beyond_source"] for r in records),
        elapsed_s=time.monotonic() - start, chunks=chunks))
    print("COMPLETE", root, time.monotonic() - start, flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--kind", choices=("2d", "3d"), default="2d")
    ap.add_argument("--seconds", type=float, default=30)
    ap.add_argument("--dx", type=float, default=1)
    ap.add_argument("--gain", type=float, default=1)
    ap.add_argument("--single", action="store_true")
    args = ap.parse_args()
    if not args.tag.replace("_", "").replace("-", "").isalnum():
        ap.error("Simple run tag required")
    if not all(np.isfinite(x) and x > 0 for x in (args.seconds, args.dx, args.gain)):
        ap.error("Positive finite parameters required")
    run(args)
