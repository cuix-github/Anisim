"""Shared, explicit spray settings and exact-frame exports for wave experiments.

No time interpolation: particle birth/death invalidates index correspondence.
Positions in solver dumps are cells; radii are SI metres, including zero on
newly emitted particles until their first droplet integration step.
"""
import json
from pathlib import Path

import numpy as np


def spray_ports(cells_per_metre):
    return dict(droplet_model=1, metres_per_cell=1.0 / cells_per_metre,
                droplet_r_median=0.0001, droplet_r_sigma=0.7,
                droplet_evap_k=1e-9, droplet_r_min=1e-6, air_nu=1.5e-5,
                spray_cfl=0.5, emit_model=1, emit_grad_threshold=10.0,
                emit_rate=0.05, emit_max_per_cell=8, pressure_solver=0)


def atomic_json(path, data):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    tmp.replace(path)


def exact_frame(source, step, destination, frame, cells_per_metre=1.0, layers=False):
    """Export one actual dump; keep radius and spray identity, never invent drops."""
    source, destination = Path(source), Path(destination)
    a = {key: np.load(source / f"{prefix}.{step:04d}.npy", mmap_mode="r")
         for key, prefix in (("pos", "frame"), ("vel", "vel"), ("phase", "ph"),
                             ("radius", "rad"), ("spray", "spr"), ("escape", "esc"))}
    n = len(a["pos"])
    if any(len(v) != n for v in a.values()):
        raise ValueError("Mismatched dump channels")
    liquid, spray, escape = a["phase"] >= 0.5, a["spray"] > 0.5, a["escape"] > 0.5
    selections = {"full": liquid}
    if layers:
        selections.update(water=liquid & ~escape & ~spray,
                          droplets=liquid & escape & ~spray, spray=spray, air=~liquid)
    counts = {}
    for name, select in selections.items():
        d = destination if name == "full" else destination / name
        d.mkdir(parents=True, exist_ok=True)
        # Chunk the selection so a full basin never needs all channels resident.
        indices = np.flatnonzero(select)
        count = len(indices)
        arrays = {
            "liquid": np.lib.format.open_memmap(d / f"liquid.{frame:04d}.npy", mode="w+",
                                               dtype=np.float32, shape=(count, 3)),
            **{key: np.lib.format.open_memmap(d / f"{key}.{frame:04d}.npy", mode="w+",
                                             dtype=np.float32, shape=(count,))
               for key in ("speed", "radius", "spray", "escape")}}
        for start in range(0, count, 250_000):
            ix = indices[start:start + 250_000]
            stop = start + len(ix)
            arrays["liquid"][start:stop] = a["pos"][ix] / cells_per_metre
            arrays["speed"][start:stop] = np.linalg.norm(a["vel"][ix], axis=1) / cells_per_metre
            for key in ("radius", "spray", "escape"):
                arrays[key][start:stop] = a[key][ix]
        for value in arrays.values():
            value.flush()
        counts[name] = count
    counts["fresh_spray_without_radius"] = int(np.count_nonzero(spray & (a["radius"] == 0)))
    return counts
