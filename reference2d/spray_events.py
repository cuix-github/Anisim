"""Independent NumPy specification of the opt-in spray event sampler.

See docs/spray-events-2d.md. No solver or C++ calls are used here.
"""
import numpy as np

DEFAULTS = dict(model=1, cloud_size=16, event_radius=3e-4, spread=0.25,
                velocity_spread=0.1, radius_max=3e-3, max_active=200000, seed=1)
GOLD = np.uint64(0x9E3779B97F4A7C15)


def mix(s):
    with np.errstate(over="ignore"):
        z = np.asarray(s, np.uint64) + GOLD
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))


def keys(seed, step, source, event=0, child=0):
    # Ordered composition: exchanging source/child must not produce the same key.
    key = mix(seed)
    for value in (step, source, event, child):
        key = mix(key ^ mix(np.asarray(value, np.uint64)))
    return key


def uniform(key, channel):
    with np.errstate(over="ignore"):
        s = key + np.uint64(channel) * GOLD
    return (mix(s) >> np.uint64(11)).astype(np.float64) * 2.0 ** -53


def requested_events(lam, key):
    lam = np.asarray(lam, float)
    if not np.isfinite(lam).all() or (lam < 0).any() or (lam > 1e6).any():
        raise ValueError("spray expected events must be finite and in [0, 1e6]")
    base = np.floor(lam).astype(np.int64)
    return base + (uniform(key, 0) < lam - base)


def ball(key, channel, dim):
    a, b = uniform(key, channel), uniform(key, channel + 1)
    theta = 2.0 * np.pi * b
    if dim == 2:
        r = np.sqrt(a)
        return np.stack([r * np.cos(theta), r * np.sin(theta)], axis=-1)
    z = 2.0 * a - 1.0
    r = uniform(key, channel + 2) ** (1.0 / 3.0)
    xy = np.sqrt(np.maximum(1.0 - z*z, 0.0))
    return r[..., None] * np.stack([xy*np.cos(theta), xy*np.sin(theta), z], axis=-1)


def cloud(seed, step, source, event, n, sigma, radius, rmin, rmax, dim=2):
    key = keys(seed, step, source, event, np.arange(n, dtype=np.uint64))
    u1 = 1.0 - uniform(key, 1)
    z = np.sqrt(-2.0*np.log(u1)) * np.cos(2.0*np.pi*uniform(key, 2))
    raw = radius * n**(-1.0/3.0) * np.exp(sigma*z - 1.5*sigma*sigma)
    return np.clip(raw, rmin, rmax), ball(key, 3, dim), ball(key, 6, dim), raw


def validate(e, d, dt=None):
    vals = [e[k] for k in ("grad_threshold", "rate", "event_radius", "spread",
                           "velocity_spread", "radius_max")]
    vals += [d[k] for k in ("r_sigma", "r_min", "metres_per_cell", "air_nu", "evap_k")]
    vals += [e[k] for k in ("cloud_size", "max_active", "max_per_cell", "seed")]
    if not np.isfinite(vals).all():
        raise ValueError("nonfinite spray configuration")
    if any(int(e[k]) != e[k] for k in ("cloud_size", "max_active", "max_per_cell", "seed")):
        raise ValueError("spray count, budget and seed parameters must be integers")
    if not (0 <= e["grad_threshold"] and 0 <= e["rate"] <= 1e3
            and 1 <= e["cloud_size"] <= 256 and int(e["cloud_size"]) == e["cloud_size"]
            and e["event_radius"] > 0 and 0 <= d["r_sigma"] <= 2
            and 0 < d["r_min"] <= e["radius_max"] and d["metres_per_cell"] > 0
            and d["air_nu"] > 0 and d["evap_k"] >= 0
            and e["spread"] >= 0 and e["velocity_spread"] >= 0
            and 0 <= e["max_active"] <= 10000000 and e["max_per_cell"] >= 0
            and 0 <= e["seed"] <= 2147483647):
        raise ValueError("invalid spray configuration")
    if dt is not None and (not np.isfinite(dt) or dt <= 0):
        raise ValueError("spray dt must be positive and finite")
