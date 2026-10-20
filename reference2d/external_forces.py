"""Independent numpy specification for generic pre-projection external fields.

Acceleration is in cells/s^2, phase is the normalized liquid volume fraction,
and nonnegative linear drag is in 1/s. Gravity has already been applied.
"""
import numpy as np


def apply_external(velocity, liquid_fraction, acceleration, drag, dt):
    a, rate = np.asarray(acceleration, float), np.asarray(drag, float)
    if not np.isfinite(a).all() or not np.isfinite(rate).all() or (rate < 0).any():
        raise ValueError("Invalid external field")
    return (velocity + dt * np.clip(liquid_fraction, 0, 1) * a) / (1 + dt * rate)
