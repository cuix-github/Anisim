"""A minimal 2D reference implementation of the core PF-FLIP loop, following Section 3 of
"Adaptive Phase-Field-FLIP for Very Large Scale Two-Phase Fluid Simulation" (Braun, Bender,
Thuerey, SIGGRAPH 2025). This is the executable spec for the Bifrost build: no adaptivity, no
escaped-particle droplets, no divergence correction — just the two-phase heart of the method.

Both phases are particles. Each step: splat mass and momentum to MAC faces with the paper's
smooth kernel (Eq. 6), read the phase field off the accumulated face masses (Eq. 7), add
gravity, solve the variable-coefficient pressure Poisson equation (Eq. 8) with face
coefficients 1/rho(phi), correct velocities, FLIP-update the particles with a per-phase blend
(Eq. 12, which doubles as viscosity via Eq. 13), and advect with RK2 through the grid field.

Scenes (run this file):  python pfflip2d.py [hydrostatic|dam_break|rayleigh_taylor|all]
Outputs go to ../results/2d/: a frame-strip PNG per scene plus printed pass/fail style numbers.
Needs numpy and Pillow only.
"""
import os
import sys
import time

import numpy as np
from PIL import Image
import spray_events

# ------------------------------------------------------------------ simulation core


class Sim:
    # Cut cells close a face less open than this: the solve weights a face by its open
    # part and barely sees a sliver, so the velocity it leaves there is not pinned down
    # (found by the node in 3D, where the water surface meets a curved wall: faces
    # under 1% open came out at up to 3 cells/s from water at rest). The node's
    # Solids::kMinOpen is the same constant.
    SOLID_MIN_OPEN = 0.01
    # Droplet radii (PF-FLIP Section 7), the defaults in SI: the log-normal's median
    # radius and the spread of ln r (Babinsky and Sojka 2002 review the form; a
    # 100 um median with sigma 0.7 spans the 30 to 300 um of spume drops over
    # breaking waves), the d^2-law evaporation constant of a water drop in air
    # (Sirignano 2010; about 1e-9 m^2/s at room temperature in dry air, slower in
    # humid marine air), the radius below which a drop is deleted, the metres one
    # cell spans, and air's kinematic viscosity for the drag's Reynolds number.
    DROPLET_DEFAULTS = dict(r_median=1e-4, r_sigma=0.7, evap_k=1e-9, r_min=1e-6,
                            metres_per_cell=1.0, air_nu=1.5e-5, integrator=1, dt_max=0.01)
    # spray emission: the gradient (1/s) above which an interface particle may emit,
    # the chance per unit excess gradient per second, the cap per cell per step
    EMIT_DEFAULTS = dict(grad_threshold=10.0, rate=0.05, max_per_cell=8, **spray_events.DEFAULTS)

    def __init__(self, nx, ny, rho_l=1000.0, rho_g=1.0, alpha_l=0.97, alpha_g=0.9, g=-9.8,
                 escape=False, esc_phi=0.3, drag_droplet=1.0, drag_bubble=8.0, buoyancy=2.0,
                 st=False, sub_advect=False, adapt=False, obstacles=None,
                 solid_phi=None, solid_cut=False, solid_support=True, piston=None,
                 droplet=None, emit=None, spray_cfl=0.0, zone=None):
        self.nx, self.ny, self.dx = nx, ny, 1.0
        self.rho_l, self.rho_g = rho_l, rho_g
        self.alpha = {1: alpha_l, 0: alpha_g}
        self.g = g
        # Escaped-particle treatment (the paper's droplets and bubbles), OFF by default
        # so every validated behaviour is unchanged. A liquid particle where the phase
        # field says air is a droplet (ballistic, gravity + drag); an air particle where
        # it says liquid is a bubble (buoyancy + drag). Both skip the FLIP blend and the
        # grid advection until they rejoin their own phase.
        self.escape, self.esc_phi = escape, esc_phi
        if escape==2 and not callable(getattr(self,'escape_phase_samples',None)):
            raise ValueError('phase-based escape requires an explicit phase contract')
        self.drag_droplet, self.drag_bubble, self.buoyancy = drag_droplet, drag_bubble, buoyancy
        self.escaped = np.zeros(0, dtype=bool)
        # Droplet radii (PF-FLIP Section 7), OFF by default: with droplet=dict(...)
        # (keys of DROPLET_DEFAULTS, any subset) an escaped liquid particle becomes a
        # droplet with its own radius in metres, drawn once from the log-normal by a
        # hash of (step, index) so the node draws the same number (droplet_normal),
        # shrunk every step by the d^2 law (r^2 falls by evap_k / 4 per second) and
        # deleted below r_min. Its drag then follows its size: a sphere's relaxation
        # rate (3/8)(rho_a / rho_l)(C_d / r)|v_rel| with the Schiller-Naumann C_d(Re),
        # applied implicitly (a 10 um drop relaxes thousands of times a second), in
        # place of the constant drag_droplet. Rejoining the water resets the radius to
        # 0, a FLIP particle again. Bubbles are untouched. Needs escape=True.
        self.droplet = None if droplet is None else dict(self.DROPLET_DEFAULTS, **droplet)
        if self.droplet is not None:
            assert escape, "droplet radii need escape=True"
            assert set(self.droplet) == set(self.DROPLET_DEFAULTS), self.droplet
        self.radius = np.zeros(0)
        self.nstep = 0
        # Extra spray emission (PF-FLIP Section 7, after Ihmsen et al. 2012), OFF by
        # default: with emit=dict(...) (keys of EMIT_DEFAULTS, any subset) a liquid
        # particle at the interface (its liquid fraction under 1 - esc_phi, the escape
        # test's own line seen from the water side) in a cell whose velocity gradient
        # (the Frobenius norm of the grid velocity's Jacobian, 1/s) exceeds
        # grad_threshold emits, with probability min(1, rate (G - G0) dt) from a hash of
        # (step, index), a spray particle at its own position and velocity: liquid,
        # massless (it adds nothing to any splat, so the flow never feels it), always
        # a droplet (a radius drawn at its first step, the drag, the evaporation),
        # deleted when it falls into the water (its liquid fraction above 1 - esc_phi)
        # or evaporates; at most max_per_cell a cell a step, in index order. The paper
        # names the method and no thresholds; these are documented choices (the
        # breaker's interface gradients at 24 s: median 2, p90 9, p99 20 per second,
        # 10 selecting the lip). Needs droplet (so escape).
        self.emit = None if emit is None else dict(self.EMIT_DEFAULTS, **emit)
        if self.emit is not None and self.emit["model"] == 0:
            self.emit = None
        if self.emit is not None:
            if self.emit["model"] not in (1, 2):
                raise ValueError("spray emission model must be 0, 1 or 2")
            assert self.droplet is not None, "spray emission needs the droplet model"
            assert set(self.emit) == set(self.EMIT_DEFAULTS), self.emit
        self.spray = np.zeros(0, dtype=bool)
        self.n_emitted = 0
        self.emission_stats = {}
        if self.emit is not None and self.emit["model"] == 2:
            spray_events.validate(self.emit, self.droplet)
        if self.droplet is not None and self.droplet["integrator"] == 2:
            if not (np.isfinite(self.droplet["dt_max"]) and self.droplet["dt_max"] > 0
                    and np.isfinite(spray_cfl) and spray_cfl > 0):
                raise ValueError("midpoint droplet flight requires positive dt_max and spray_cfl")
        # The droplets' own integrator (PF-FLIP Section 7: RK2 with locally adaptive
        # time stepping), OFF at 0: with spray_cfl > 0 a droplet's step is cut into
        # sub-steps that move it at most spray_cfl cells each; over a sub-step the
        # drag relaxation toward the air velocity sampled at its start is taken
        # exactly, with gravity, at the rate of the sub-step's midpoint (the rate
        # depends on the slip: a half step at the starting rate finds it), and the
        # position advances by the trapezoid of the sub-step's two velocities:
        # second order in the sub-step. An explicit RK2 on the drag itself would
        # need thousands of sub-steps for a 10 um drop, which relaxes in
        # microseconds; the relaxation taken exactly needs none.
        self.spray_cfl = float(spray_cfl)
        # ST-FLIP (Braun et al. 2026) spatiotemporal sampling, OFF by default. Each
        # particle carries a time jitter tau in (-1/2, 1/2] slab units; splats are
        # weighted by the paper's one-sided temporal poly6 kernel (peak at +1/2), and
        # advection folds the jitter change into one per-particle step of length
        # dt*(1 + tau_new - tau_old), sub-stepped at local CFL <= 1. Spatial kernel
        # deviation from the paper (ours is radial, theirs separable) noted; the
        # milestone-4 finding says kernel shape barely moves the fields.
        self.st = st
        self.sub_advect = sub_advect
        # Adaptive particles (paper Section 5, 2D spec): air particles in deep air
        # merge 4-to-1 into scale-2 particles (mass x4, splat radius x2); coarse
        # particles split back near the interface or walls. The per-size blend
        # correction equalizes viscosity: (1-alpha_coarse) = (1-alpha_fine)/scale^2.
        self.adapt = adapt
        self.scale = np.zeros(0)
        self.coarse_gain = 1.0
        # Solid obstacles as axis-aligned boxes [(xmin, ymin, xmax, ymax), ...]:
        # faces inside a box carry beta = 0 and zero velocity (a wall, anywhere);
        # particles that end up inside are pushed out along the axis of least
        # penetration, losing that velocity component, exactly like the outer walls.
        self.obstacles = [tuple(map(float, o)) for o in (obstacles or [])]
        # Solids from a signed distance field, a mesh collider voxelized to a level
        # set: solid_phi holds the distance at the grid nodes (cell corners, shape
        # (nx+1, ny+1)) in cells, negative inside. Binary faces close when the solid
        # covers their midpoint, the boxes' rule (phi taken linear along a face, so the
        # midpoint value is the mean of its two ends). With solid_cut a face is weighted
        # by its open fraction instead, the variational cut-cell treatment of Batty,
        # Bertails and Bridson (2007): a slope is a slope, not a staircase. Particles
        # that end up inside move out along the field's gradient and lose their normal
        # velocity, as the boxes' least-penetration push does. solid_support reads each
        # face's phase against the share of its splat reach outside the solid (below).
        self.solid_phi = None if solid_phi is None else np.asarray(solid_phi, np.float64)
        if self.solid_phi is not None:
            assert self.solid_phi.shape == (nx + 1, ny + 1), self.solid_phi.shape
        self.solid_cut = bool(solid_cut)
        self._open = None
        self.solid_support = bool(solid_support)
        self._support = None
        # A piston wavemaker (2026-10-16): the paddle of a wave flume, a plate standing
        # on the floor that slides along x. piston = (S, T, x0, width, top): the plate's
        # left face at x0 + S/2 (1 - cos 2 pi t / T), its thickness width, its top at
        # y = top, below the ceiling so the air it displaces flows over it. (A whole
        # wall sweeping into the sealed tank of incompressible fluid has no pressure
        # solution: the divergence sums to the swept rate, the solve spins to its cap
        # and the velocities blow up. A plate sweeps nothing net.) Every step rebuilds
        # the level set, folded with the static solid_phi given here if any, as a
        # cut-cell solid moving at the step's mean speed: the faces it closes carry
        # that speed, its part of a cut face carries the speed into the divergence
        # (the variational form), a pushed-out particle keeps the plate's speed, and
        # the support share is re-probed within three cells of the plate only. At rest
        # (S = 0) every path is the static solid's, bit for bit (validate_piston.py).
        self.piston = None
        if piston is not None:
            given = tuple(map(float, piston))
            self.piston = given + (8.0, 2.0, 0.7 * ny)[len(given) - 2:]
        self.t = 0.0
        self.static_phi = self.solid_phi
        self.piston_x = 0.0
        self.solid_vel = np.zeros(2)
        self.support_band = None
        self._support_static = None
        self._open_static = None
        self._plate_phi = None
        if self.piston is not None:
            self.solid_cut, self.solid_support = True, True
            self._move_piston(0.0, 0.0)
        # the wave zone (below): None leaves every path untouched
        self.zone = None
        self._zone_t = None          # the end of the step under way, for the inlet
        if zone is not None:
            self.zone = dict(self.ZONE_DEFAULTS, **zone)
            self.zone["waves"] = [tuple(map(float, wv)) for wv in self.zone["waves"]]
        self.st_rng = np.random.default_rng(7)
        self.tau = np.zeros(0)
        self.u = np.zeros((nx + 1, ny))   # x-velocity on vertical faces
        self.v = np.zeros((nx, ny + 1))   # y-velocity on horizontal faces
        self.pos = np.zeros((0, 2))
        self.vel = np.zeros((0, 2))
        self.typ = np.zeros(0, dtype=np.int32)   # 1 liquid, 0 air
        self.rho0_face = None                    # calibrated splat mass of a full face

    # --- seeding: 2x2 particles per cell, jittered, everywhere; type from a mask function
    def seed(self, liquid_mask, jitter=0.35, seed=1):
        rng = np.random.default_rng(seed)
        cx, cy = np.meshgrid(np.arange(self.nx), np.arange(self.ny), indexing="ij")
        offs = [(0.25, 0.25), (0.75, 0.25), (0.25, 0.75), (0.75, 0.75)]
        ps, ts = [], []
        for ox, oy in offs:
            p = np.stack([cx + ox, cy + oy], -1).reshape(-1, 2).astype(np.float64)
            p += rng.uniform(-jitter, jitter, p.shape) * 0.5
            ps.append(p)
            ts.append(liquid_mask(p[:, 0], p[:, 1]).astype(np.int32))
        self.pos = np.concatenate(ps)
        self.typ = np.concatenate(ts)
        if self.has_solids:
            keep = ~self._inside_solid(self.pos)
            self.pos, self.typ = self.pos[keep], self.typ[keep]
        self.vel = np.zeros_like(self.pos)
        self.scale = np.ones(len(self.pos))
        self.tau = np.zeros(len(self.pos))      # particles start synchronized
        self.xi = np.ones(len(self.pos))        # per-particle jitter strength
        self.radius = np.zeros(len(self.pos))   # droplet radius in metres, 0 = none
        self.spray = np.zeros(len(self.pos), dtype=bool)   # emitted spray, massless

    @staticmethod
    def _splitmix(s):
        """splitmix64's output for the uint64 state array s."""
        z = s + np.uint64(0x9E3779B97F4A7C15)
        z = (z ^ (z >> np.uint64(30))) * np.uint64(0xBF58476D1CE4E5B9)
        z = (z ^ (z >> np.uint64(27))) * np.uint64(0x94D049BB133111EB)
        return z ^ (z >> np.uint64(31))

    @staticmethod
    def droplet_normal(step_index, idx):
        """A standard normal per particle from a hash of (step, index), the same in the
        node: splitmix64 at (step + 1) * golden + index and at the next state, the two
        outputs' top 53 bits as uniforms (the first on (0, 1], the second on [0, 1)),
        Box-Muller."""
        with np.errstate(over="ignore"):
            s = (np.uint64(step_index + 1) * np.uint64(0x9E3779B97F4A7C15)
                 + np.asarray(idx, np.uint64))
            z1 = Sim._splitmix(s)
            z2 = Sim._splitmix(s + np.uint64(1))
        u1 = ((z1 >> np.uint64(11)).astype(np.float64) + 1.0) * 2.0 ** -53
        u2 = (z2 >> np.uint64(11)).astype(np.float64) * 2.0 ** -53
        return np.sqrt(-2.0 * np.log(u1)) * np.cos(2.0 * np.pi * u2)

    @staticmethod
    def droplet_uniform(step_index, idx):
        """A uniform on [0, 1) per particle from the hash of (step, index), the state
        two past the normal's first (splitmix64 at (step + 1) * golden + index + 2),
        the same in the node."""
        with np.errstate(over="ignore"):
            s = (np.uint64(step_index + 1) * np.uint64(0x9E3779B97F4A7C15)
                 + np.asarray(idx, np.uint64) + np.uint64(2))
            z = Sim._splitmix(s)
        return (z >> np.uint64(11)).astype(np.float64) * 2.0 ** -53

    def _droplet_radii(self, drop, dt):
        """The radius bookkeeping of one step: rejoined droplets back to 0, new ones
        drawn, every droplet evaporated; returns the mask of those gone below r_min."""
        D, r = self.droplet, self.radius
        new = drop & (r <= 0.0)
        r[~drop] = 0.0
        if new.any():
            idx = np.flatnonzero(new)
            r[idx] = D["r_median"] * np.exp(D["r_sigma"] * self.droplet_normal(self.nstep, idx))
        if drop.any():
            r[drop] = np.sqrt(np.maximum(r[drop] ** 2 - 0.25 * D["evap_k"] * dt, 0.0))
        return drop & (r < D["r_min"])

    def _droplet_drag(self, d, vg, dt):
        """Size-dependent drag on the droplets d toward the air velocity vg (cells per
        second), implicit in the relaxation rate."""
        D = self.droplet
        r = np.maximum(self.radius[d], D["r_min"])
        rel = (vg[d] - self.vel[d]) * D["metres_per_cell"]
        speed = np.sqrt((rel ** 2).sum(1))
        re = 2.0 * r * speed / D["air_nu"]
        # C_d |v_rel|: Schiller-Naumann below Re 1000 (12 nu / r (1 + 0.15 Re^0.687),
        # Stokes' 12 nu / r at the limit), the constant 0.44 above
        cdv = np.where(re <= 1000.0, 12.0 * D["air_nu"] / r * (1.0 + 0.15 * re ** 0.687),
                       0.44 * speed)
        k = (0.375 * (self.rho_g / self.rho_l) * cdv / r * dt)[:, None]
        self.vel[d] = (self.vel[d] + k * vg[d]) / (1.0 + k)

    def _droplet_fly(self, d, times):
        """The droplets d over their own times (seconds each): sub-steps that move a
        droplet at most spray_cfl cells, the drag relaxation toward the air velocity
        (sampled at the sub-step's start) taken exactly over each with gravity, the
        position by the trapezoid of the sub-step's two velocities."""
        D = self.droplet
        idx = np.flatnonzero(d)
        rem = np.asarray(times, float).copy()
        pos, vel = self.pos[idx].copy(), self.vel[idx].copy()
        r = np.maximum(self.radius[idx], D["r_min"])
        g = np.array([0.0, self.g])
        hi = np.array([self.nx - 0.51, self.ny - 0.51])
        def rate(v, vg):
            """The relaxation rate k (1/s) at slip vg - v: (3/8)(rho_a / rho_l)(C_d / r)|v_rel|."""
            rel = (vg - v) * D["metres_per_cell"]
            sp = np.sqrt((rel ** 2).sum(1))
            re = 2.0 * r * sp / D["air_nu"]
            cdv = np.where(re <= 1000.0, 12.0 * D["air_nu"] / r * (1.0 + 0.15 * re ** 0.687),
                           0.44 * sp)
            return 0.375 * (self.rho_g / self.rho_l) * cdv / r

        def relax(v, vg, k, h):
            """v after h seconds of relaxation at rate k toward vg under gravity."""
            e = np.exp(-k * h)
            gk = np.where(k * h > 1e-6, (1.0 - e) / np.maximum(k, 1e-300), h * (1.0 - 0.5 * k * h))
            return vg + (v - vg) * e[:, None] + g[None, :] * gk[:, None]

        for _ in range(1000000):
            act = rem > 1e-12
            if not act.any():
                break
            vg = self.sample_velocity(np.clip(pos, 0.51, hi))
            mid_carrier = D["integrator"] == 2
            speed = np.maximum(np.abs(vel).max(1), 1e-9)
            if mid_carrier:
                speed = np.maximum(speed, np.abs(vg).max(1))
            h = np.where(act, np.minimum(rem, self.spray_cfl / speed), 0.0)
            if mid_carrier:
                if not np.isfinite(D["dt_max"]) or D["dt_max"] <= 0:
                    raise ValueError("droplet dt_max must be positive and finite")
                h = np.minimum(h, D["dt_max"])
            # the rate depends on the slip, so it is taken at the sub-step's midpoint
            # (a half step at the starting rate finds it), then the whole step at it
            vmid = relax(vel, vg, rate(vel, vg), 0.5 * h)
            if mid_carrier:
                pmid = pos + 0.25 * h[:, None] * (vel + vmid)
                vg = self.sample_velocity(np.clip(pmid, 0.51, hi))
            vnew = relax(vel, vg, rate(vmid, vg), h)
            pos = pos + 0.5 * h[:, None] * (vel + vnew)
            vel, rem = vnew, rem - h
        self.pos[idx], self.vel[idx] = pos, vel

    def _grad_norm(self):
        """The Frobenius norm of the grid velocity's Jacobian per cell (1/s): the
        cell-centred velocity from the faces, central differences inside and one-sided
        at the walls, a cell as the unit of length."""
        uc = 0.5 * (self.u[:-1, :] + self.u[1:, :])
        vc = 0.5 * (self.v[:, :-1] + self.v[:, 1:])
        g2 = np.zeros((self.nx, self.ny))
        for f in (uc, vc):
            gx, gy = np.gradient(f)
            g2 += gx ** 2 + gy ** 2
        return np.sqrt(g2)

    def _emit(self, dt, lf):
        """Spray emission for this step from the interface particles (liquid fraction
        lf per particle), appended to the arrays; returns how many."""
        E, n = self.emit, len(self.pos)
        c = np.clip(np.floor(self.pos).astype(np.int64), 0, [self.nx - 1, self.ny - 1])
        gp = self._grad_norm()[c[:, 0], c[:, 1]]
        cand = ((self.typ == 1) & ~self.spray & ~self.escaped & (lf < 1.0 - self.esc_phi)
                & (gp > E["grad_threshold"]))
        idx = np.flatnonzero(cand)
        if E["model"] == 2:
            return self._emit_clouds(dt, idx, c, gp)
        if len(idx) == 0:
            return 0
        prob = np.clip(E["rate"] * (gp[idx] - E["grad_threshold"]) * dt, 0.0, 1.0)
        hit = idx[self.droplet_uniform(self.nstep, idx) < prob]
        if len(hit) == 0:
            return 0
        # at most max_per_cell a cell, the lowest indices first
        key = c[hit, 0] * self.ny + c[hit, 1]
        order = np.argsort(key, kind="stable")
        ks = key[order]
        starts = np.r_[0, np.flatnonzero(np.diff(ks)) + 1]
        rank = np.empty(len(ks), np.int64)
        rank[order] = np.arange(len(ks)) - np.repeat(starts, np.diff(np.r_[starts, len(ks)]))
        hit = hit[rank < E["max_per_cell"]]
        m = len(hit)
        self.pos = np.concatenate([self.pos, self.pos[hit]])
        self.vel = np.concatenate([self.vel, self.vel[hit]])
        self.typ = np.concatenate([self.typ, np.ones(m, dtype=self.typ.dtype)])
        self.scale = np.concatenate([self.scale, np.ones(m)])
        self.tau = np.concatenate([self.tau, np.zeros(m)])
        if len(self.xi) == n:
            self.xi = np.concatenate([self.xi, np.ones(m)])
        self.radius = np.concatenate([self.radius, np.zeros(m)])
        self.spray = np.concatenate([self.spray, np.ones(m, dtype=bool)])
        self.escaped = np.concatenate([self.escaped, np.ones(m, dtype=bool)])
        return m

    def _emit_clouds(self, dt, idx, cells, gp):
        E, D = self.emit, self.droplet
        spray_events.validate(E, D, dt)
        N, seed = int(E["cloud_size"]), int(E["seed"])
        lam = E["rate"] * (gp[idx] - E["grad_threshold"]) * dt
        events = spray_events.requested_events(lam, spray_events.keys(seed, self.nstep, idx))
        stats = dict(candidates=len(idx), events=int(events.sum()), requested=int(events.sum())*N,
                     accepted=0, rejected_cap=0, rejected_geometry=0, raw_r3=0.0, clipped_r3=0.0)
        per_cell, points, velocities, radii = {}, [], [], []
        # Include particles scheduled for removal this step in the resource budget;
        # this conservative rule is identical in the node and reference.
        room = max(0, int(E["max_active"]) - int(self.spray.sum()))
        hi = np.array([self.nx-0.51, self.ny-0.51])
        for pi, count in zip(idx, events):
            cell = tuple(cells[pi])
            for ev in range(int(count)):
                available = min(room, int(E["max_per_cell"]) - per_cell.get(cell, 0))
                if available <= 0:
                    stats["rejected_cap"] += (int(count)-ev)*N
                    break
                r, dp, dv, raw = spray_events.cloud(seed, self.nstep, int(pi), ev, N,
                    D["r_sigma"], E["event_radius"], D["r_min"], E["radius_max"])
                p = self.pos[pi] + E["spread"]*dp
                v = self.vel[pi] + E["velocity_spread"]*np.linalg.norm(self.vel[pi])*dv
                good = ((p >= 0.51) & (p <= hi)).all(1)
                if self.has_solids:
                    good &= ~self._inside_solid(p)
                # Rejection precedence is geometry first, budget second.
                stats["rejected_geometry"] += int((~good).sum())
                keep = np.flatnonzero(good)[:available]
                m = len(keep)
                stats["rejected_cap"] += int(good.sum()) - m
                points.extend(p[keep]); velocities.extend(v[keep]); radii.extend(r[keep])
                stats["raw_r3"] += float((raw[keep]**3).sum())
                stats["clipped_r3"] += float((r[keep]**3).sum())
                stats["accepted"] += m
                per_cell[cell] = per_cell.get(cell, 0) + m
                room -= m
        m, n = len(radii), len(self.pos)
        self.emission_stats = stats
        if m:
            self.pos = np.concatenate([self.pos, np.asarray(points)])
            self.vel = np.concatenate([self.vel, np.asarray(velocities)])
            self.typ = np.concatenate([self.typ, np.ones(m, dtype=self.typ.dtype)])
            self.scale = np.concatenate([self.scale, np.ones(m)])
            self.tau = np.concatenate([self.tau, np.zeros(m)])
            if len(self.xi) == n:
                self.xi = np.concatenate([self.xi, np.ones(m)])
            self.radius = np.concatenate([self.radius, np.asarray(radii)])
            self.spray = np.concatenate([self.spray, np.ones(m, bool)])
            self.escaped = np.concatenate([self.escaped, np.ones(m, bool)])
        return m

    def _delete(self, keep):
        n = len(keep)
        for name in ("pos", "vel", "typ", "scale", "tau", "xi", "radius", "escaped", "spray"):
            a = getattr(self, name, None)       # xi exists only once seeded
            if a is not None and len(a) == n:
                setattr(self, name, a[keep])

    @staticmethod
    def plate_phi(nx, ny, x_left, width, top):
        """The node distances of a plate on the floor: x_left to x_left + width, up to top."""
        X, Y = np.meshgrid(np.arange(nx + 1, dtype=float), np.arange(ny + 1, dtype=float),
                           indexing="ij")
        return np.maximum(np.maximum(x_left - X, X - (x_left + width)), Y - top)

    # --- the wave zone (2026-10-18): the wavemaker of a numerical wave basin without
    # a paddle, the relaxation zone of ship and coastal CFD (Jacobsen et al. 2012).
    # zone = dict(length, depth, ramp, waves=[(H, T, phase), ...]) in cells, seconds and
    # radians: after every projection the face velocities within `length` of the
    # x = 0 wall are blended toward the sum of the linear trains' velocity fields, and
    # so are the liquid particles' velocities after their FLIP update, with the weight
    # 1 - (exp(s^3.5) - 1) / (e - 1), s = x / length, a face's times its phase (the
    # water is forced where it is, the air left alone); the x = 0 wall faces are the
    # trains' inlet, carrying their velocity (times the face's phase) instead of 0
    # whenever a step is under way, so the flux the trains need is not forced into the
    # wall's cells. (The first zone faded the weight over the cell above the trains'
    # summed surface instead, which pumped the water one way: lifted on the crests, the
    # water above the trains' troughs was left where it was, and the level rose.) Each train's wavenumber
    # solves w^2 = g k tanh(k h) at the still depth h; its height rises from 0 over
    # `ramp` seconds (a smoothstep). The field is taken at the end of the step. The 3D
    # node's zone_* ports and `waves` array are this (its trains may also have a
    # direction; here they run along x).
    ZONE_DEFAULTS = dict(length=24.0, depth=12.0, ramp=0.0, waves=())

    @staticmethod
    def wavenumber(T, h, g):
        """k from w^2 = g k tanh(k h), Newton from the deep-water k (as the node)."""
        om = 2.0 * np.pi / T
        k = om * om / g
        for _ in range(50):
            th, ch = np.tanh(k * h), np.cosh(k * h)
            k -= (g * k * th - om * om) / (g * (th + k * h / (ch * ch)))
        return k

    def _zone_trains(self, t):
        z = self.zone
        g = abs(self.g)
        ramp = 1.0
        if z["ramp"] > 0.0:
            s = min(1.0, t / z["ramp"])
            ramp = s * s * (3.0 - 2.0 * s)
        out = []
        for H, T, ph in z["waves"]:
            if not (T > 0.0) or not (H > 0.0):
                continue
            out.append((0.5 * H * ramp, self.wavenumber(T, z["depth"], g), 2.0 * np.pi / T, ph))
        return out

    def _zone_at(self, x, y, t):
        """The blend weight at the points (x, y) and the trains' velocity there."""
        z = self.zone
        h = z["depth"]
        x, y = np.asarray(x, float), np.asarray(y, float)
        s = np.clip(x / z["length"], 0.0, 1.0)
        wz = np.where(x < z["length"], 1.0 - (np.exp(s ** 3.5) - 1.0) / 1.718281828459045, 0.0)
        ye = np.minimum(y, h)
        ut, vt = np.zeros_like(x), np.zeros_like(x)
        for a, k, om, ph in self._zone_trains(t):
            psi = k * x - om * t + ph
            sh = np.sinh(k * h)
            ut += a * om * np.cosh(k * ye) / sh * np.cos(psi)
            vt += a * om * np.sinh(k * ye) / sh * np.sin(psi)
        return wz, ut, vt

    def _zone_faces(self, t):
        nx, ny = self.nx, self.ny
        pu, pv = self._face_phase()
        i, j = np.meshgrid(np.arange(nx + 1, dtype=float), np.arange(ny, dtype=float) + 0.5, indexing="ij")
        w, ut, _ = self._zone_at(i, j, t)
        w = w * pu
        self.u = (1.0 - w) * self.u + w * ut
        i, j = np.meshgrid(np.arange(nx, dtype=float) + 0.5, np.arange(ny + 1, dtype=float), indexing="ij")
        w, _, vt = self._zone_at(i, j, t)
        w = w * pv
        self.v = (1.0 - w) * self.v + w * vt

    def _zone_particles(self, t, skip=None):
        liquid = self.typ == 1
        if skip is not None:
            liquid = liquid & ~skip
        if not liquid.any():
            return
        w, ut, vt = self._zone_at(self.pos[liquid, 0], self.pos[liquid, 1], t)
        v = self.vel[liquid]
        v[:, 0] = (1.0 - w) * v[:, 0] + w * ut
        v[:, 1] = (1.0 - w) * v[:, 1] + w * vt
        self.vel[liquid] = v

    def _move_piston(self, t, dt):
        """The plate at time t, moving at the step's mean speed towards t + dt."""
        S, T, x0, width, top = self.piston
        x_left = lambda tt: x0 + 0.5 * S * (1.0 - np.cos(2.0 * np.pi * tt / T))   # noqa: E731
        xl = float(x_left(t))
        self.piston_x = xl
        self.solid_vel = np.array([(float(x_left(t + dt)) - xl) / dt if dt > 0 else 0.0, 0.0])
        phi = self.plate_phi(self.nx, self.ny, xl, width, top)
        self._plate_phi = phi
        self.solid_phi = phi if self.static_phi is None else np.minimum(phi, self.static_phi)
        self._open = self._support = None
        self.support_band = (xl - 3.0, xl + width + 3.0)

    def masses(self):
        m = np.where(self.typ == 1, self.rho_l, self.rho_g)
        if self.adapt:
            m = m * self.scale ** 2 * np.where(self.scale > 1.5, self.coarse_gain, 1.0)
        if len(self.spray) == len(m):
            m = np.where(self.spray, 0.0, m)      # emitted spray adds nothing to the grid
        return m

    # --- Eq. 6 kernel weights against a set of face centres, vectorised over a 3x3 stencil
    def _splat(self, centres_shape, face_of, values_list):
        """Splat each entry of values_list (per-particle scalars) onto faces. face_of maps a
        particle to its base face index. Uniform runs use the r = dx kernel on a 3x3
        stencil; adaptive runs splat each size group with its own radius."""
        out = [np.zeros(centres_shape) for _ in values_list]
        if self.adapt and len(self.scale) == len(face_of[0]):
            for sc in np.unique(self.scale):
                sel = self.scale == sc
                self._splat_group(out, centres_shape,
                                  (face_of[0][sel], face_of[1][sel]),
                                  [np.asarray(v)[sel] if np.ndim(v) else v
                                   for v in values_list], float(sc))
        else:
            self._splat_group(out, centres_shape, face_of, values_list, 1.0)
        return out

    @staticmethod
    def _splat_group(out, centres_shape, face_of, values_list, r):
        base, frac = face_of
        span = int(np.ceil(r))
        for di in range(-span, span + 1):
            for dj in range(-span, span + 1):
                ii = base[:, 0] + di
                jj = base[:, 1] + dj
                ok = (ii >= 0) & (ii < centres_shape[0]) & (jj >= 0) & (jj < centres_shape[1])
                d2 = ((frac[:, 0] - di) ** 2 + (frac[:, 1] - dj) ** 2) / r ** 2
                w = np.maximum(1.0 - d2, 0.0) ** 3
                w = np.where(ok, w, 0.0)
                iis = np.clip(ii, 0, centres_shape[0] - 1)
                jjs = np.clip(jj, 0, centres_shape[1] - 1)
                for o, val in zip(out, values_list):
                    np.add.at(o, (iis, jjs), w * val)

    def _face_frames(self):
        """Base indices and fractional offsets of every particle against u-faces and v-faces."""
        pu = self.pos - np.array([0.0, 0.5])   # u-face centres are at (i, j+0.5)
        pv = self.pos - np.array([0.5, 0.0])   # v-face centres are at (i+0.5, j)
        bu = np.floor(pu + 0.5).astype(np.int64)
        bv = np.floor(pv + 0.5).astype(np.int64)
        return (bu, pu - bu), (bv, pv - bv)

    @staticmethod
    def _wt(tau):
        """ST-FLIP temporal kernel: one-sided poly6 peaking at the slab end (+1/2)."""
        return (35.0 / 16.0) * np.maximum(1.0 - (tau - 0.5) ** 2, 0.0) ** 3

    _GL_NODES = 0.5 * np.polynomial.legendre.leggauss(8)[0]
    _GL_WEIGHTS = 0.5 * np.polynomial.legendre.leggauss(8)[1]

    def _wt_norm(self, xi):
        """E[wT(xi*T)], T~U(-1/2,1/2): the weight renormalization for attenuated
        jitter, so calm (xi<1) and fast (xi=1) regions deposit consistent mass."""
        taus = xi[:, None] * self._GL_NODES[None, :]
        return (self._wt(taus) * self._GL_WEIGHTS[None, :]).sum(1)

    def p2g(self):
        m = self.masses()
        if self.st:
            m = m * self._wt(self.tau) / np.maximum(self._wt_norm(self.xi), 1e-9)
        fu, fv = self._face_frames()
        mu, pu = self._splat(self.u.shape, fu, [m, m * self.vel[:, 0]])
        mv, pv = self._splat(self.v.shape, fv, [m, m * self.vel[:, 1]])
        self.mu, self.mv = mu, mv
        with np.errstate(divide="ignore", invalid="ignore"):
            self.u = np.where(mu > 0, pu / np.maximum(mu, 1e-300), 0.0)
            self.v = np.where(mv > 0, pv / np.maximum(mv, 1e-300), 0.0)
        self._enforce_walls()

    @property
    def has_solids(self):
        return bool(self.obstacles) or self.solid_phi is not None

    def _inside_solid(self, pts):
        inside = np.zeros(len(pts), dtype=bool)
        for x0, y0, x1, y1 in self.obstacles:
            inside |= ((pts[:, 0] > x0) & (pts[:, 0] < x1) &
                       (pts[:, 1] > y0) & (pts[:, 1] < y1))
        if self.solid_phi is not None:
            inside |= self._phi_at(pts)[0] < 0.0
        return inside

    def _phi_at(self, pts, P=None):
        """The solid's signed distance at points, bilinear over the nodes, and its gradient."""
        P = self.solid_phi if P is None else P
        q = np.clip(pts, 0.0, [self.nx, self.ny])
        b = np.minimum(np.floor(q).astype(np.int64), [self.nx - 1, self.ny - 1])
        fx, fy = (q - b)[:, 0], (q - b)[:, 1]
        i, j = b[:, 0], b[:, 1]
        p00, p10, p01, p11 = P[i, j], P[i + 1, j], P[i, j + 1], P[i + 1, j + 1]
        phi = ((1 - fx) * (1 - fy) * p00 + fx * (1 - fy) * p10 +
               (1 - fx) * fy * p01 + fx * fy * p11)
        gx = (1 - fy) * (p10 - p00) + fy * (p11 - p01)
        gy = (1 - fx) * (p01 - p00) + fx * (p11 - p10)
        return phi, np.stack([gx, gy], 1)

    @staticmethod
    def _open_fraction(a, b):
        """The part of a face where the solid's distance is positive, the distance
        linear between the face's two ends a and b."""
        with np.errstate(divide="ignore", invalid="ignore"):
            t = a / (a - b)                  # the zero crossing, measured from a
        return np.where((a > 0) & (b > 0), 1.0,
                        np.where((a <= 0) & (b <= 0), 0.0, np.where(a > 0, t, 1.0 - t)))

    def _solid_support_share(self):
        """Per-face share of the splat kernel's reach that lies outside the level-set
        solid (1 away from it): a fine all-liquid probe of the domain, splatted with and
        without the solid's interior, as calibrate() measures a full face. The phase
        reads a face's mass against this share, so a face beside the solid is not taken
        for air because the solid holds no particles."""
        if self._support is None:
            if self.piston is None:
                self._support = self._probe_share(None, self.solid_phi)
            else:
                # the static solids once over the whole domain; the plate's
                # neighbourhood every step, written over the faces whose kernel reach
                # (one cell) lies inside the probed band
                if self._support_static is None:
                    self._support_static = self._probe_share(None, self.static_phi)
                lo, hi = self.support_band
                band = self._probe_share((lo, hi), self.solid_phi)
                su, sv = (a.copy() for a in self._support_static)
                xu = np.arange(self.nx + 1, dtype=float)
                xv = np.arange(self.nx, dtype=float) + 0.5
                cu = (xu >= lo + 1.0) & (xu <= hi - 1.0)
                cv = (xv >= lo + 1.0) & (xv <= hi - 1.0)
                su[cu, :] = band[0][cu, :]
                sv[cv, :] = band[1][cv, :]
                self._support = (su, sv)
        return self._support

    def _probe_share(self, band, P):
        """The share for every face from a fine all-liquid probe of the domain (of the
        columns between band's bounds when given), splatted with and without the
        interior of the solid whose node distances are P; all ones when P is None."""
        if P is None:
            return np.ones(self.u.shape), np.ones(self.v.shape)
        k = 6
        pts = np.stack(np.meshgrid((np.arange(self.nx * k) + 0.5) / k,
                                   (np.arange(self.ny * k) + 0.5) / k,
                                   indexing="ij"), -1).reshape(-1, 2)
        if band is not None:
            pts = pts[(pts[:, 0] >= band[0]) & (pts[:, 0] <= band[1])]
        keep = (self._phi_at(pts, P)[0] > 0.0).astype(float)
        out = []
        for shape, off in ((self.u.shape, (0.0, 0.5)), (self.v.shape, (0.5, 0.0))):
            q = pts - np.array(off)
            base = np.floor(q + 0.5).astype(np.int64)
            acc = [np.zeros(shape), np.zeros(shape)]
            self._splat_group(acc, shape, (base, q - base), [np.ones(len(pts)), keep], 1.0)
            out.append(np.where(acc[0] > 0, acc[1] / np.maximum(acc[0], 1e-300), 1.0))
        return tuple(out)

    def _solid_open(self):
        """Open fraction of every u and v face against the level-set solid (1 = open)."""
        if self._open is None:
            P = self.solid_phi
            a_u, b_u = P[:, :-1], P[:, 1:]    # u face (i, j): nodes (i, j) and (i, j+1)
            a_v, b_v = P[:-1, :], P[1:, :]    # v face (i, j): nodes (i, j) and (i+1, j)
            if self.solid_cut:
                self._open = tuple(np.where(w < self.SOLID_MIN_OPEN, 0.0, w) for w in
                                   (self._open_fraction(a_u, b_u),
                                    self._open_fraction(a_v, b_v)))
            else:
                self._open = ((0.5 * (a_u + b_u) > 0).astype(float),
                              (0.5 * (a_v + b_v) > 0).astype(float))
        return self._open

    def _static_open(self):
        """Open fraction of every face against the static solids alone (all open without
        any): what the moving plate covers of a face is this less the combined fraction."""
        if self._open_static is None:
            if self.static_phi is None:
                self._open_static = (np.ones(self.u.shape), np.ones(self.v.shape))
            else:
                P = self.static_phi
                self._open_static = tuple(np.where(w < self.SOLID_MIN_OPEN, 0.0, w) for w in
                                          (self._open_fraction(P[:, :-1], P[:, 1:]),
                                           self._open_fraction(P[:-1, :], P[1:, :])))
        return self._open_static

    def _solid_face_masks(self):
        iu, ju = np.meshgrid(np.arange(self.nx + 1), np.arange(self.ny) + 0.5,
                             indexing="ij")
        iv, jv = np.meshgrid(np.arange(self.nx) + 0.5, np.arange(self.ny + 1),
                             indexing="ij")
        mu = np.zeros(self.u.shape, dtype=bool)
        mv = np.zeros(self.v.shape, dtype=bool)
        for x0, y0, x1, y1 in self.obstacles:
            mu |= (iu >= x0) & (iu <= x1) & (ju >= y0) & (ju <= y1)
            mv |= (iv >= x0) & (iv <= x1) & (jv >= y0) & (jv <= y1)
        if self.solid_phi is not None:
            wu, wv = self._solid_open()
            mu |= wu <= 0.0
            mv |= wv <= 0.0
        return mu, mv

    def _push_out(self):
        for x0, y0, x1, y1 in self.obstacles:
            inside = ((self.pos[:, 0] > x0) & (self.pos[:, 0] < x1) &
                      (self.pos[:, 1] > y0) & (self.pos[:, 1] < y1))
            if not inside.any():
                continue
            px, py = self.pos[inside, 0], self.pos[inside, 1]
            dxm = np.stack([px - x0, x1 - px, py - y0, y1 - py], 1)
            side = np.argmin(dxm, 1)
            tx = np.where(side == 0, x0 - 0.01, np.where(side == 1, x1 + 0.01, px))
            ty = np.where(side == 2, y0 - 0.01, np.where(side == 3, y1 + 0.01, py))
            self.pos[inside, 0], self.pos[inside, 1] = tx, ty
            vz = self.vel[inside]
            vz[side <= 1, 0] = 0.0
            vz[side >= 2, 1] = 0.0
            self.vel[inside] = vz
        if self.solid_phi is not None:
            for _ in range(3):
                phi, grad = self._phi_at(self.pos)
                inside = phi < 0.0
                if not inside.any():
                    break
                g = grad[inside]
                g2 = np.maximum((g * g).sum(1), 1e-12)
                n = g / np.sqrt(g2)[:, None]
                # a Newton step to phi = +0.01: one unit step on a true distance field,
                # longer where the bilinear field is no distance (in a corner cell its
                # slope can be a tenth of one), at most 1.5 cells
                step = -(phi[inside] - 0.01)[:, None] * g / g2[:, None]
                ln = np.linalg.norm(step, axis=1)
                vs = self.solid_vel
                if self.piston is not None and self.static_phi is not None:
                    # the plate's speed for a particle inside the plate, none inside a
                    # static solid; read before the move (after it the particle is out)
                    in_plate = self._phi_at(self.pos[inside], self._plate_phi)[0] < 0.0
                    vs = self.solid_vel[None, :] * in_plate[:, None]
                self.pos[inside] += step * np.minimum(1.0, 1.5 / np.maximum(ln, 1e-12))[:, None]
                vz = self.vel[inside]
                vz -= ((vz - vs) * n).sum(1)[:, None] * n
                self.vel[inside] = vz

    def _enforce_walls(self):
        self.u[0, :] = self.u[-1, :] = 0.0
        self.v[:, 0] = self.v[:, -1] = 0.0
        if self.zone is not None and self._zone_t is not None:
            # the wave inlet: the x = 0 faces carry the trains' velocity less its
            # flux-weighted mean (a sealed box has no pressure solution under a net
            # inflow; the wall's own column breathes by the mean instead)
            w, ut, _ = self._zone_at(np.zeros(self.ny), np.arange(self.ny, dtype=float) + 0.5, self._zone_t)
            w = w * self._face_phase()[0][0, :]
            c = float(np.sum(w * ut) / np.sum(w)) if np.sum(w) > 0.0 else 0.0
            self.u[0, :] = w * (ut - c)
        if self.has_solids:
            mu, mv = self._solid_face_masks()
            self.u[mu] = 0.0
            self.v[mv] = 0.0
            if self.piston is not None:          # the faces the plate closes move with it
                wu, wv = self._solid_open()
                ou, ov = self._static_open()
                self.u[(wu <= 0.0) & (ou > 0.0)] = self.solid_vel[0]
                self.v[(wv <= 0.0) & (ov > 0.0)] = self.solid_vel[1]

    # --- Eq. 7: the phase field from raw splatted face mass
    def phase(self, raw):
        if self.rho0_face is None:
            raise RuntimeError("call calibrate() after seeding")
        eta = np.log(self.rho_l / self.rho_g)
        rho_min = eta * self.rho_g * self.rho0_face
        phi = np.sqrt(np.maximum(raw - rho_min, 0.0) / (self.rho0_face * self.rho_l))
        return np.minimum(phi, 1.0)

    def calibrate(self):
        """rho0_face: the raw mass a face accumulates from a uniform all-liquid seeding,
        with the same temporal weighting the simulation will use (averaged over a few
        jitter realizations when ST sampling is on, per the ST-FLIP m0 recipe)."""
        probe = Sim(8, 8, self.rho_l, self.rho_g, st=self.st)
        probe.seed(lambda x, y: np.ones_like(x, dtype=bool), jitter=0.0)
        fu, _ = probe._face_frames()
        reps = 5 if self.st else 1
        acc = []
        for _ in range(reps):
            m = probe.masses()
            if self.st:
                probe.tau = probe.st_rng.uniform(-0.5, 0.5, len(probe.pos))
                m = m * self._wt(probe.tau)   # full-jitter reference scale (xi = 1)
            (mu,) = probe._splat(probe.u.shape, fu, [m])
            acc.append(np.median(mu[2:-2, 2:-2]))
        self.rho0_face = float(np.mean(acc)) / self.rho_l
        if self.adapt:
            # measure the raw face mass an all-coarse uniform seeding produces, and
            # gain-correct coarse masses so fine and coarse fields read identically.
            probe2 = Sim(8, 8, self.rho_l, self.rho_g)
            probe2.seed(lambda x, y: np.ones_like(x, dtype=bool), jitter=0.0)
            keep = (probe2.pos[:, 0] % 2 < 1) & (probe2.pos[:, 1] % 2 < 1)
            probe2.pos = probe2.pos[keep] + 0.5
            probe2.typ = probe2.typ[keep]
            probe2.adapt = True
            probe2.scale = np.full(len(probe2.pos), 2.0)
            probe2.coarse_gain = 1.0
            m2 = probe2.masses() * 4.0 / 4.0
            fu2, _ = probe2._face_frames()
            (mu2,) = probe2._splat(probe2.u.shape, fu2, [probe2.masses()])
            raw2 = np.median(mu2[2:-2, 2:-2]) / self.rho_l
            self.coarse_gain = float(np.mean(acc)) / self.rho_l / max(raw2, 1e-9)

    def _face_phase(self):
        """The faces' phase (Eq. 7) as the projection takes it."""
        if self.solid_phi is not None and self.solid_support:
            su, sv = self._solid_support_share()
            return self.phase(self.mu / np.maximum(su, 0.1)), self.phase(self.mv / np.maximum(sv, 0.1))
        return self.phase(self.mu), self.phase(self.mv)

    def project(self, dt, tol=1e-6, max_iter=4000):
        if not np.isfinite([dt,tol]).all() or dt<=0 or tol<=0 or max_iter<0:
            raise ValueError('invalid pressure timestep, tolerance or iteration limit')
        phi_u, phi_v = self._face_phase()
        rho_u = self.rho_g + (self.rho_l - self.rho_g) * phi_u
        rho_v = self.rho_g + (self.rho_l - self.rho_g) * phi_v
        bu = 1.0 / rho_u
        bv = 1.0 / rho_v
        bu[0, :] = bu[-1, :] = 0.0     # closed walls
        bv[:, 0] = bv[:, -1] = 0.0
        if self.has_solids:
            mu_s, mv_s = self._solid_face_masks()
            bu[mu_s] = 0.0
            bv[mv_s] = 0.0

        div = (self.u[1:, :] - self.u[:-1, :]) + (self.v[:, 1:] - self.v[:, :-1])
        ba_u, ba_v = bu, bv           # face coefficients of the operator
        if self.solid_phi is not None and self.solid_cut:
            # cut cells: flux passes through the open part of a face only, so the
            # operator's coefficients and the divergence take the open fraction; the
            # velocity update below does not. A moving solid's part of a face carries
            # the solid's speed into the divergence (Batty, Bertails, Bridson 2007).
            wu, wv = self._solid_open()
            uw, vw = wu * self.u, wv * self.v
            if self.piston is not None:
                # the part of a face the plate covers (open against the static solids,
                # closed against all) moves at the plate's speed; a static solid's part
                # does not
                ou, ov = self._static_open()
                uw = uw + (ou - wu) * self.solid_vel[0]
                vw = vw + (ov - wv) * self.solid_vel[1]
            div = (uw[1:, :] - uw[:-1, :]) + (vw[:, 1:] - vw[:, :-1])
            ba_u, ba_v = bu * wu, bv * wv
        b = -div / dt
        # Default None leaves every legacy operation unchanged. The opt-in
        # volume subclass supplies a compatible target, div(u_after)=target.
        target_div = getattr(self, '_target_divergence', None)
        if target_div is not None:
            b = b + target_div / dt

        def A(p):
            gpu = np.zeros_like(self.u)
            gpv = np.zeros_like(self.v)
            gpu[1:-1, :] = (p[1:, :] - p[:-1, :]) * ba_u[1:-1, :]
            gpv[:, 1:-1] = (p[:, 1:] - p[:, :-1]) * ba_v[:, 1:-1]
            return -((gpu[1:, :] - gpu[:-1, :]) + (gpv[:, 1:] - gpv[:, :-1]))

        diag = ba_u[1:, :] + ba_u[:-1, :] + ba_v[:, 1:] + ba_v[:, :-1]
        inv_diag = np.where(diag > 0, 1.0 / np.maximum(diag, 1e-300), 0.0)

        p = np.zeros((self.nx, self.ny))
        r = b - A(p)
        z = inv_diag * r
        d = z.copy()
        rz = float((r * z).sum())
        r0 = np.sqrt(float((r * r).sum())) or 1.0
        it = 0
        self.pressure_restarts = 0
        for it in range(max_iter):
            if np.sqrt(float((r * r).sum())) / r0 < tol:
                actual = b - A(p)
                certified = float(np.linalg.norm(actual) / r0)
                if not np.isfinite(certified) or certified <= tol:
                    break
                # Recursive residuals accumulate rounding error. Continue from
                # this iterate with b-Ap, within the original iteration budget.
                r = actual
                z = inv_diag * r
                d = z.copy()
                rz = float((r * z).sum())
                self.pressure_restarts += 1
            q = A(d)
            a = rz / float((d * q).sum())
            p += a * d
            r -= a * q
            z = inv_diag * r
            rz_new = float((r * z).sum())
            d = z + (rz_new / rz) * d
            rz = rz_new
        self.iters = it

        self.pressure_residual = float(np.linalg.norm(b-A(p))/r0)
        if not np.isfinite(self.pressure_residual) or self.pressure_residual>tol:
            raise RuntimeError(f'pressure did not converge: true residual {self.pressure_residual:.17g}, tolerance {tol:.17g}')

        self.u[1:-1, :] -= dt * bu[1:-1, :] * (p[1:, :] - p[:-1, :])
        self.v[:, 1:-1] -= dt * bv[:, 1:-1] * (p[:, 1:] - p[:, :-1])
        self._enforce_walls()

    def _sample_faces(self, grid, offset, pts):
        q = pts - offset
        b = np.floor(q).astype(np.int64)
        f = q - b
        out = np.zeros(len(pts))
        for di in (0, 1):
            for dj in (0, 1):
                w = (f[:, 0] if di else 1 - f[:, 0]) * (f[:, 1] if dj else 1 - f[:, 1])
                ii = np.clip(b[:, 0] + di, 0, grid.shape[0] - 1)
                jj = np.clip(b[:, 1] + dj, 0, grid.shape[1] - 1)
                out += w * grid[ii, jj]
        return out

    def sample_velocity(self, pts):
        return np.stack([self._sample_faces(self.u, np.array([0.0, 0.5]), pts),
                         self._sample_faces(self.v, np.array([0.5, 0.0]), pts)], -1)

    def g2p(self, u_old, v_old, skip=None):
        du, dv = self.u - u_old, self.v - v_old
        delta = np.stack([self._sample_faces(du, np.array([0.0, 0.5]), self.pos),
                          self._sample_faces(dv, np.array([0.5, 0.0]), self.pos)], -1)
        pic = self.sample_velocity(self.pos)
        a1 = np.where(self.typ == 1, self.alpha[1], self.alpha[0])
        if self.adapt:
            a1 = 1.0 - (1.0 - a1) / self.scale ** 2
        a = a1[:, None]
        blended = a * (self.vel + delta) + (1 - a) * pic
        if skip is not None:
            blended[skip] = self.vel[skip]
        self.vel = blended

    def advect(self, dt, ballistic=None, fly=None):
        if fly is not None and fly.any():
            self._droplet_fly(fly, np.full(int(fly.sum()), float(dt)))
        old = self.pos
        mid = self.pos + 0.5 * dt * self.sample_velocity(self.pos)
        self.pos = self.pos + dt * self.sample_velocity(np.clip(mid, 0.51, None))
        if ballistic is not None and ballistic.any():
            move = ballistic if fly is None else (ballistic & ~fly)
            self.pos[move] = old[move] + dt * self.vel[move]
        if fly is not None and fly.any():
            self.pos[fly] = old[fly]              # flown already, above
        lo = 0.51
        hix, hiy = self.nx - 0.51, self.ny - 0.51
        for k, hi in ((0, hix), (1, hiy)):
            below, above = self.pos[:, k] < lo, self.pos[:, k] > hi
            self.pos[:, k] = np.clip(self.pos[:, k], lo, hi)
            self.vel[below | above, k] = 0.0

    def step(self, dt):
        if self.piston is not None:
            self._move_piston(self.t, dt)
        if self.zone is not None:
            self._zone_t = self.t + dt
        self.p2g()
        u_star, v_star = self.u.copy(), self.v.copy()
        self.v += self.g * dt          # gravity on y-faces
        self._enforce_walls()
        self.project(dt)
        if self.zone is not None:
            self._zone_faces(self.t + dt)
            self._enforce_walls()
        esc = None
        if self.escape:
            # Escape detection from the OTHER phase's local number density, which is
            # self-excluding: a particle's own splat cannot mask it (phi-based
            # detection fails for droplets because a lone liquid particle's own mass,
            # through Eq. 7's square root, still reads as phi ~ 0.45).
            w_rep = self.scale ** 2 if self.adapt else np.ones(len(self.typ))
            ones_air = (self.typ == 0).astype(float) * w_rep
            ones_liq = ((self.typ == 1) & ~self.spray).astype(float) * w_rep
            fu, fv = self._face_frames()
            au, lu = self._splat(self.u.shape, fu, [ones_air, ones_liq])
            av, lv = self._splat(self.v.shape, fv, [ones_air, ones_liq])
            thresh = (1.0 - self.esc_phi) * self.rho0_face
            frac_air = 0.5 * (self._sample_faces(au, np.array([0.0, 0.5]), self.pos) +
                              self._sample_faces(av, np.array([0.5, 0.0]), self.pos))
            frac_liq = 0.5 * (self._sample_faces(lu, np.array([0.0, 0.5]), self.pos) +
                              self._sample_faces(lv, np.array([0.5, 0.0]), self.pos))
            drop = (self.typ == 1) & (frac_air > thresh)
            # Opt-in classification audit. esc_phi=.3 is the lower threshold
            # (paper delta_phi=.2), NOT the paper's delta itself. Keep the
            # legacy count rule byte-for-byte for escape=True / 1.
            if self.escape == 2:
                phase_at_particle = self.escape_phase_samples()
                drop = (self.typ == 1) & (phase_at_particle < self.esc_phi)
            if self.zone is not None:
                drop &= ~(self.pos[:, 0] < self.zone["length"])   # the zone's liquid is the zone's to move
            bub = (self.typ == 0) & (frac_liq > thresh)
            if self.escape == 2:
                bub = (self.typ == 0) & (phase_at_particle > 1.0-self.esc_phi)
            # emitted spray is a droplet whatever the test says, and falls back into
            # the water (deleted below) where the test would call it liquid
            drop = drop | self.spray
            inside = self.spray & (frac_liq > thresh)
            if self.escape == 2:
                inside = self.spray & (phase_at_particle > 1.0-self.esc_phi)
            esc = drop | bub
            self.escaped = esc
            lf = frac_liq / self.rho0_face
        self.g2p(u_star, v_star, skip=esc)
        if self.zone is not None:
            self._zone_particles(self.t + dt, esc)
        gone = None
        if self.droplet is not None and esc is not None:
            gone = self._droplet_radii(drop, dt)
            deaths = dict(evaporated=int((gone & self.spray & ~inside).sum()),
                          reentered=int(inside.sum()))
            if inside.any():
                gone = gone | inside
        if self.escape == 2:
            # PF section 3.4 deletes escaped air inclusions. This is independent
            # of whether the optional water-droplet flight model is enabled.
            gone = bub.copy() if gone is None else gone | bub
        fly = None
        if esc is not None and esc.any():
            vg = self.sample_velocity(self.pos)
            d = drop
            if self.droplet is not None and self.spray_cfl > 0.0:
                fly = d                      # the droplets' own integrator, in the advection
            else:
                self.vel[d, 1] += self.g * dt
                if self.droplet is None:
                    self.vel[d] += self.drag_droplet * dt * (vg[d] - self.vel[d])
                elif d.any():
                    self._droplet_drag(d, vg, dt)
            b = bub
            self.vel[b, 1] += self.buoyancy * abs(self.g) * dt
            self.vel[b] += self.drag_bubble * dt * (vg[b] - self.vel[b])
        if self.st or self.sub_advect:
            self._advect_st(dt, ballistic=esc, fly=fly)
        else:
            self.advect(dt, ballistic=esc, fly=fly)
        if self.has_solids:
            self._push_out()
        if self.emit is not None and esc is not None:
            m = self._emit(dt, lf)
            self.n_emitted = m
            if self.emit["model"] == 2:
                self.emission_stats.update(deaths, step=self.nstep)
            if m and gone is not None:
                gone = np.concatenate([gone, np.zeros(m, bool)])
        if gone is not None and gone.any():
            self._delete(~gone)
        if self.adapt:
            self._adapt_particles()
        if self.emit is not None and self.spray.any():
            # the spray kept behind every other particle, so the others keep the
            # order (and the ST jitter stream) they would have without it
            order = np.concatenate([np.flatnonzero(~self.spray), np.flatnonzero(self.spray)])
            for name in ("pos", "vel", "typ", "scale", "tau", "xi", "radius", "escaped", "spray"):
                a = getattr(self, name, None)
                if a is not None and len(a) == len(order):
                    setattr(self, name, a[order])
        self.t += dt
        if self.emit is not None and self.emit["model"] == 2:
            self.emission_stats["alive"] = int(self.spray.sum())
        self.nstep += 1

    def _liq_frac(self):
        ones_liq = ((self.typ == 1) & ~self.spray).astype(float) * self.scale ** 2
        fu, fv = self._face_frames()
        (lu,) = self._splat(self.u.shape, fu, [ones_liq])
        (lv,) = self._splat(self.v.shape, fv, [ones_liq])
        return 0.5 * (self._sample_faces(lu, np.array([0.0, 0.5]), self.pos) +
                      self._sample_faces(lv, np.array([0.5, 0.0]), self.pos)) / self.rho0_face

    def _adapt_particles(self, merge_thresh=0.02, split_thresh=0.05, wall_margin=2.0):
        lf = self._liq_frac()
        # ---- split coarse near the interface or walls
        near_wall = ((self.pos[:, 0] < wall_margin) | (self.pos[:, 0] > self.nx - wall_margin) |
                     (self.pos[:, 1] < wall_margin) | (self.pos[:, 1] > self.ny - wall_margin))
        to_split = (self.scale > 1.5) & ((lf > split_thresh) | near_wall)
        if to_split.any():
            base = self.pos[to_split]
            offs = np.array([[-0.5, -0.5], [0.5, -0.5], [-0.5, 0.5], [0.5, 0.5]])
            new_pos = (base[:, None, :] + offs[None, :, :]).reshape(-1, 2)
            new_vel = np.repeat(self.vel[to_split], 4, axis=0)
            keep = ~to_split
            self.pos = np.concatenate([self.pos[keep], new_pos])
            self.vel = np.concatenate([self.vel[keep], new_vel])
            self.typ = np.concatenate([self.typ[keep], np.repeat(self.typ[to_split], 4)])
            self.scale = np.concatenate([self.scale[keep], np.ones(len(new_pos))])
            self.tau = np.concatenate([self.tau[keep], np.repeat(self.tau[to_split], 4)])
            self.radius = np.concatenate([self.radius[keep], np.zeros(len(new_pos))])
            self.spray = np.concatenate([self.spray[keep], np.zeros(len(new_pos), bool)])
            if len(self.escaped) == len(keep):
                self.escaped = np.concatenate([self.escaped[keep], np.zeros(len(new_pos), bool)])
            self.xi = np.concatenate([self.xi[keep] if len(self.xi) == len(keep) else
                                      np.ones(keep.sum()), np.ones(len(new_pos))])
            lf = np.concatenate([lf[keep], np.full(len(new_pos), 1.0)])  # keep them fine this step
        # ---- merge deep-air fine 4-to-1 per 2x2 block
        cand = (self.scale < 1.5) & (self.typ == 0) & (lf < merge_thresh)
        if cand.sum() >= 4:
            idx = np.where(cand)[0]
            blocks = (self.pos[idx, 0] // 2).astype(np.int64) * 100000 + \
                     (self.pos[idx, 1] // 2).astype(np.int64)
            order = np.argsort(blocks, kind="stable")
            idx, blocks = idx[order], blocks[order]
            _, starts, counts = np.unique(blocks, return_index=True, return_counts=True)
            merged_members, new_p, new_v, new_tau = [], [], [], []
            for st_i, ct in zip(starts, counts):
                take = (ct // 4) * 4
                for g in range(0, take, 4):
                    mem = idx[st_i + g: st_i + g + 4]
                    merged_members.append(mem)
                    new_p.append(self.pos[mem].mean(0))
                    new_v.append(self.vel[mem].mean(0))
                    new_tau.append(self.tau[mem].mean())
            if merged_members:
                drop = np.concatenate(merged_members)
                keep = np.ones(len(self.pos), bool)
                keep[drop] = False
                self.pos = np.concatenate([self.pos[keep], np.array(new_p)])
                self.vel = np.concatenate([self.vel[keep], np.array(new_v)])
                self.typ = np.concatenate([self.typ[keep],
                                           np.zeros(len(new_p), dtype=np.int32)])
                self.scale = np.concatenate([self.scale[keep], np.full(len(new_p), 2.0)])
                self.tau = np.concatenate([self.tau[keep], np.array(new_tau)])
                self.radius = np.concatenate([self.radius[keep], np.zeros(len(new_p))])
                self.spray = np.concatenate([self.spray[keep], np.zeros(len(new_p), bool)])
                if len(self.escaped) == len(keep):
                    self.escaped = np.concatenate([self.escaped[keep], np.zeros(len(new_p), bool)])
                self.xi = np.concatenate([self.xi[keep] if len(self.xi) == len(keep) else
                                          np.ones(int(keep.sum())), np.ones(len(new_p))])

    def _advect_st(self, dt, ballistic=None, cfl_local=1.0, max_rounds=64, fly=None):
        """ST-FLIP advection: each particle advances by dt*(1 + tau_new - tau_old),
        folding the jitter change into one step, sub-stepped at local CFL <= 1."""
        if self.st:
            # Velocity-adaptive jitter attenuation (the paper's calm-water fix):
            # full jitter where a particle crosses a cell per step, none where the
            # flow is still, via a smoothstep of the local per-particle CFL.
            cfl_p = np.abs(self.vel).max(1) * dt / self.dx
            t = np.clip(cfl_p, 0.0, 1.0)
            self.xi = t * t * (3.0 - 2.0 * t)
            # the jitter drawn for every particle but the emitted spray, which flies
            # unjittered (tau 0, massless): the stream the others consume is then the
            # one they would consume without it
            tau_new = np.zeros(len(self.pos))
            prim = ~self.spray if len(self.spray) == len(self.pos) else np.ones(len(self.pos), bool)
            tau_new[prim] = self.xi[prim] * self.st_rng.uniform(-0.5, 0.5, int(prim.sum()))
            remaining = dt * (1.0 + tau_new - self.tau)
            self.tau = tau_new
        else:
            remaining = np.full(len(self.pos), float(dt))
        if fly is not None and fly.any():
            # the droplets' own integrator over their full jittered step
            self._droplet_fly(fly, remaining[fly])
            remaining = np.where(fly, 0.0, remaining)
        if ballistic is not None and ballistic.any():
            # escaped particles fly on their own velocity for their full jittered step
            move = ballistic if fly is None else (ballistic & ~fly)
            self.pos[move] += remaining[move, None] * self.vel[move]
            remaining = np.where(ballistic, 0.0, remaining)
        for _ in range(max_rounds):
            active = remaining > 1e-12
            if not active.any():
                break
            v = self.sample_velocity(self.pos)
            speed = np.maximum(np.abs(v).max(1), 1e-9)
            dt_sub = np.where(active, np.minimum(remaining, cfl_local * self.dx / speed), 0.0)
            mid = self.pos + 0.5 * dt_sub[:, None] * v
            v2 = self.sample_velocity(np.clip(mid, 0.51, None))
            self.pos = self.pos + dt_sub[:, None] * v2
            remaining = remaining - dt_sub
        lo = 0.51
        hix, hiy = self.nx - 0.51, self.ny - 0.51
        for k, hi in ((0, hix), (1, hiy)):
            below, above = self.pos[:, k] < lo, self.pos[:, k] > hi
            self.pos[:, k] = np.clip(self.pos[:, k], lo, hi)
            self.vel[below | above, k] = 0.0

    def run(self, frames, dt_frame, cfl=0.5, on_frame=None, blowup=1e3):
        """Frames of dt_frame, each in CFL substeps. A particle faster than blowup cells
        per frame aborts the run (a solve gone wrong makes picosecond substeps that
        never end); None disables the check."""
        for f in range(frames):
            remaining = dt_frame
            while remaining > 1e-9:
                vmax = max(1e-6, float(np.abs(self.vel).max()))
                if blowup is not None and vmax > blowup * self.dx / dt_frame:
                    raise RuntimeError(f"velocity blew up: {vmax:.3g} cells/s at t {self.t:.3f} s, frame {f}")
                dt = min(remaining, cfl * self.dx / vmax, dt_frame / 2)
                self.step(dt)
                remaining -= dt
            if on_frame:
                on_frame(f, self)


# ------------------------------------------------------------------ rendering helpers

def draw(sim, scale=3):
    img = np.full((sim.nx * scale, sim.ny * scale, 3), 245, np.uint8)
    pix = (sim.pos * scale).astype(np.int64)
    pix[:, 0] = np.clip(pix[:, 0], 0, sim.nx * scale - 1)
    pix[:, 1] = np.clip(pix[:, 1], 0, sim.ny * scale - 1)
    air = sim.typ == 0
    img[pix[air, 0], pix[air, 1]] = (210, 210, 215)
    liq = sim.typ == 1
    img[pix[liq, 0], pix[liq, 1]] = (30, 90, 200)
    return Image.fromarray(np.rot90(img))


def save_strip(frames_imgs, path, cols):
    w, h = frames_imgs[0].size
    rows = (len(frames_imgs) + cols - 1) // cols
    sheet = Image.new("RGB", (w * cols, h * rows), "white")
    for i, im in enumerate(frames_imgs):
        sheet.paste(im, ((i % cols) * w, (i // cols) * h))
    sheet.save(path)
    print("wrote", path)


# ------------------------------------------------------------------ scenes

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "2d")


def hydrostatic():
    sim = Sim(64, 64)
    sim.seed(lambda x, y: y < 32)
    sim.calibrate()
    vmax_log = []
    sim.run(30, 0.05, on_frame=lambda f, s: vmax_log.append(float(np.abs(s.vel).max())))
    print("hydrostatic: max |v| over last 5 frames = %.4g (want: small, no blow-up); iters/frame ~%d"
          % (max(vmax_log[-5:]), sim.iters))
    return max(vmax_log[-5:]) < 1.0


def dam_break():
    sim = Sim(160, 80)
    sim.seed(lambda x, y: (x < 40) & (y < 56))
    sim.calibrate()
    shots, fronts = [], []

    def cb(f, s):
        liq_x = s.pos[s.typ == 1, 0]
        fronts.append(float(np.percentile(liq_x, 99.5)))
        if f % 4 == 0:
            shots.append(draw(s))

    sim.run(41, 0.06, on_frame=cb)
    save_strip(shots, os.path.join(RESULTS, "dam_break.png"), cols=3)
    mono = all(b >= a - 1.0 for a, b in zip(fronts, fronts[1:]))
    print("dam_break: front x from %.1f to %.1f of 160 (want: advancing); monotonic-ish: %s"
          % (fronts[0], max(fronts), mono))
    return mono and max(fronts) > 80  # mildly damped by the FLIP blend; tune alpha when matching Bifrost


def rayleigh_taylor():
    sim = Sim(64, 192, rho_l=19.0, rho_g=1.0, alpha_l=0.97, alpha_g=0.97)  # Atwood 0.9, as the paper's RT setup
    mid, amp = 96, 3.0
    sim.seed(lambda x, y: y > mid + amp * np.cos(2 * np.pi * x / 64))       # heavy phase on top
    sim.calibrate()
    shots, tips = [], []

    def cb(f, s):
        heavy_y = s.pos[s.typ == 1, 1]
        tips.append(float(np.percentile(heavy_y, 0.5)))
        if f % 5 == 0:
            shots.append(draw(s, scale=2))

    sim.run(46, 0.12, on_frame=cb)
    save_strip(shots, os.path.join(RESULTS, "rayleigh_taylor.png"), cols=5)
    print("rayleigh_taylor: heavy-phase tip fell from y=%.1f to y=%.1f (want: fingers growing downward)"
          % (tips[0], min(tips)))
    return min(tips) < tips[0] - 30


if __name__ == "__main__":
    os.makedirs(RESULTS, exist_ok=True)
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    ok = True
    for name, fn in (("hydrostatic", hydrostatic), ("dam_break", dam_break),
                     ("rayleigh_taylor", rayleigh_taylor)):
        if which in (name, "all"):
            print("===", name)
            ok = fn() and ok
    print("total %.1f s, overall: %s" % (time.time() - t0, "PASS" if ok else "CHECK OUTPUTS"))
