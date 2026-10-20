"""A discrete curl is an independent 3D, pre-projection acceleration anchor."""
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "gates"))
from validate_3d import run_node
from validate_external_forces import FIELDS


def main():
    out = ROOT / "out/runs/force_wave_validation_20261020/curl3d"
    source, result = out / "in", out / "out"
    source.mkdir(parents=True, exist_ok=True)
    nx, ny, nz = 16, 12, 8
    x = np.arange(nx + 1)[:, None, None]
    y = (np.arange(ny) + .5)[None, :, None]
    z = np.arange(nz + 1)[None, None, :]
    psi = np.sin(np.pi * x / nx) * np.sin(np.pi * z / nz) * (1 + .2 * np.cos(np.pi * y / ny))
    au, aw = np.diff(psi, axis=2), -np.diff(psi, axis=0)
    av = np.zeros((nx, ny + 1, nz))
    dt = float(np.float32(.02))
    arrays = {}
    for axis, acc in zip("uvw", (au, av, aw)):
        arrays.update({axis + "_mass": np.ones(acc.size), axis + "_phase": np.ones(acc.size),
                       axis + "_mom": np.zeros(acc.size), "liquid_acceleration_" + axis: acc.ravel(),
                       "linear_drag_" + axis: np.zeros(1)})
    arrays["positions"] = np.stack(np.meshgrid(np.arange(nx) + .5, np.arange(ny) + .5,
                                               np.arange(nz) + .5, indexing="ij"), -1).reshape(-1, 3)
    arrays["velocities"] = np.zeros_like(arrays["positions"])
    arrays["particle_phase"] = np.ones(nx * ny * nz)
    for name, values in arrays.items():
        np.save(source / (name + ".npy"), values.astype(np.float32))
    stats = run_node(str(source), str(result), extra_arrays=FIELDS, nx=nx, ny=ny, nz=nz,
                     dt=dt, gravity=0, external_forces=1, preconditioner=1)
    errors = {axis: float(np.max(np.abs(np.load(result / f"out_{axis}.npy") -
                                       dt * acc.astype(np.float32).ravel())))
              for axis, acc in zip("uvw", (au, av, aw))}
    assert max(errors.values()) < 1e-7, errors
    assert np.max(np.abs(np.load(result / "out_w.npy"))) > .001
    # The field already has zero discrete divergence and zero normal wall flow.
    # Projection must preserve it, independently checking w storage and indexing.
    receipt = dict(passed=True, grid=[nx, ny, nz], max_velocity_errors=errors, statistics=stats)
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
