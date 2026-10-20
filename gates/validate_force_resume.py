"""Check absolute field sampling across a real Bifrost checkpoint boundary."""
import json
from pathlib import Path
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "graphs"))
from loop3d import run_bifrost3, load_final_state, FINAL_STATE


def main():
    source = ROOT / "out/runs/force_single_calibrated_20261020"
    plan = json.loads((source / "inputs/plan.json").read_text())
    first = source / "dumps/c0"
    init = {k: np.load(first / f"{k}.npy") for k in
            ("positions", "velocities", "particle_phase", "mass", "momx", "momy", "momz",
             "tau", "wt", "wtph", "scale", "phase_st", "radius", "spray")}
    solid = (np.load(first / "solid_phi.npy").reshape(tuple(n + 1 for n in plan["grid"])), 2)
    out = ROOT / "out/runs/force_wave_validation_20261020/resume"
    def run(name, state, steps, offset):
        d = out / name
        run_bifrost3(name, *plan["grid"], None, steps, plan["dt"], extra_ports=plan["ports"],
                     init_arrays=state, final_state=True, step_offset=offset, fused=True,
                     solid=solid, out_dir=str(d))
        return d
    control = run("full", init, 16, 0)
    a = run("split_a", init, 8, 0)
    b = run("split_b", load_final_state(str(a)), 8, 8)
    fields = ["final"] + ["final_" + n for n in FINAL_STATE]
    identical = all(np.array_equal(np.load(control / f"{n}.npy"), np.load(b / f"{n}.npy")) for n in fields)
    assert identical
    bad = run("wrong_clock", load_final_state(str(a)), 8, 0)
    delta = float(np.max(np.abs(np.load(b / "final_velocities.npy") - np.load(bad / "final_velocities.npy"))))
    assert delta > 1e-6, delta
    receipt = dict(passed=True, final_state_bit_identical=identical,
                   wrong_clock_velocity_difference=delta, compared_channels=fields)
    (out / "receipt.json").write_text(json.dumps(receipt, indent=2) + "\n")
    print(json.dumps(receipt, indent=2))


if __name__ == "__main__":
    main()
