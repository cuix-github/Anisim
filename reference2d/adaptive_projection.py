"""Milestone 6.5a, stage 3a: the adaptive pressure solve inside the simulation loop.

Stages 1 and 2 (adaptive_solve.py) validated one solve on frozen states and
checked corrected velocities only on the always-fine interface band. A running
simulation needs more: the particles read velocities from the UNIFORM face
grid everywhere, so every uniform face must leave the projection corrected,
including the faces inside a coarse cell, where the adaptive solve has a single
pressure and says nothing. This module supplies that last mile and plugs the
adaptive solve into pfflip2d.Sim in place of its uniform projection.

Faces between two adaptive cells take the pressure gradient of the adaptive
face they lie on (Delta / d, the stage-1 interpolated difference over the
center distance: s for same-level faces, 1.5 s for junction sub-faces), each
fine face with its own beta, so the fine fluxes sum to the adaptive face's
flux. Faces strictly inside a coarse cell get a small fine-level Neumann solve
per coarse cell, boundary faces fixed: the least-divergence interior field
compatible with the adaptive boundary fluxes. Cheaper fills (no correction,
or gradients interpolated linearly across the cell) left divergence of 0.2 to
0.8 per step inside coarse cells on the gate's dam state, against 1e-6 here.

Two junction treatments survive validate_adaptive_projection.py, and which one
production uses is a 3D cost-versus-quality question, left open on purpose:

  naive    the conservative two-cell stencil over the 1.5 s center distance
           (Losasso et al. 2004, the treatment the PF-FLIP paper adopts): one
           solve, every face flux leaves one cell and enters the other, but
           first order at junctions, so a fluid at rest picks up small
           parasitic currents there (its error is dt g / 3 in hydrostatics)
  cleanup  the second-order variational solve, whose junction fluxes balance
           only in a weighted sense (the uniform faces keep a divergence of
           ~1.8 per step next to junctions), followed by a naive-stencil
           projection of just that leftover: exact books, and the first-order
           error now acts on the small leftover instead of on gravity, at the
           cost of a second solve

Rejected on the way, numbers in the README's stage-3a record: the variational
solve alone (books do not balance), and moving the leftover junction flux
across the face the coarse cell shares with its lateral neighbor (balances
the books without touching the pressure, but is not a projection, so it can
add energy; the dam break blew up by frame 23).
"""
import numpy as np

from adaptive_solve import Multigrid, Operator, build_cells


def band_target(sim, band=5, levels=2):
    """Required finest level per fine cell: fine within `band` cells of the
    liquid interface and next to the walls, one level coarser out to 2*band,
    `levels` beyond. The stage-1 criterion, whose dam-break curve put the
    operating point at band >= 5."""
    nx, ny = sim.nx, sim.ny
    liq = np.zeros((nx, ny))
    np.add.at(liq, (np.clip(sim.pos[:, 0].astype(int), 0, nx - 1),
                    np.clip(sim.pos[:, 1].astype(int), 0, ny - 1)),
              (sim.typ == 1).astype(float))
    frac = liq / 4.0
    interface = (frac > 0.05) & (frac < 0.95)
    for ax, sh in ((0, 1), (0, -1), (1, 1), (1, -1)):
        interface |= (frac > 0.5) != (np.roll(frac, sh, axis=ax) > 0.5)

    def grown(m, k):
        for _ in range(k):
            g = m.copy()
            for ax, sh in ((0, 1), (0, -1), (1, 1), (1, -1)):
                g |= np.roll(m, sh, axis=ax)
            m = g
        return m
    near = grown(interface, band)
    mid = grown(near, 2 * band)
    target = np.full((nx, ny), levels, int)
    target[mid] = 1
    target[near] = 0
    target[:2, :] = target[-2:, :] = 0
    target[:, :2] = target[:, -2:] = 0
    return target


class AdaptiveProjection:
    """Callable replacement for Sim.project: sim.project = AdaptiveProjection(sim)."""

    def __init__(self, sim, junction="cleanup", band=5, levels=2, tol=1e-6,
                 multigrid=True):
        assert junction in ("naive", "cleanup")
        assert not sim.obstacles, "obstacles: not in the stage-3a scope yet"
        self.sim, self.junction = sim, junction
        self.band, self.levels, self.tol, self.multigrid = band, levels, tol, multigrid
        self.log = []     # per-call: unknowns, iterations, divergence stats

    # the same face betas as Sim.project
    def _betas(self):
        s = self.sim
        rho_u = s.rho_g + (s.rho_l - s.rho_g) * s.phase(s.mu)
        rho_v = s.rho_g + (s.rho_l - s.rho_g) * s.phase(s.mv)
        bu, bv = 1.0 / rho_u, 1.0 / rho_v
        bu[0, :] = bu[-1, :] = 0.0
        bv[:, 0] = bv[:, -1] = 0.0
        return bu, bv

    def __call__(self, dt, tol=None, max_iter=None):
        s = self.sim
        nx, ny = s.nx, s.ny
        bu, bv = self._betas()

        def beta(x, y, axis, span):
            if axis == 0:
                ys = (y - span / 2.0 + np.arange(span)).astype(int)
                return float(bu[int(round(x)), ys].mean())
            xs = (x - span / 2.0 + np.arange(span)).astype(int)
            return float(bv[xs, int(round(y))].mean())

        cells = build_cells(nx, ny, self.levels, band_target(s, self.band, self.levels))
        stencil = "naive" if self.junction == "naive" else "variational"
        op, p, it, bnorm = self._solve(cells, beta, stencil, dt, self.tol)
        gu, gv = self._gradient(op, p, s)
        it2 = 0

        # faces inside coarse cells (NaN gradients) stay at u* until the local
        # solves below fill them
        ustar, vstar = s.u.copy(), s.v.copy()
        s.u -= dt * bu * np.nan_to_num(gu)
        s.v -= dt * bv * np.nan_to_num(gv)
        if self.junction == "cleanup":
            # the variational fluxes balance only in the weighted sense; a
            # second, conservative (naive-stencil) projection of the corrected
            # field removes the leftover junction dipoles. It is a true
            # projection, so it cannot add energy, and its first-order error
            # only acts on those dipoles, not on gravity. Its tolerance is set
            # so its absolute residual matches the main solve's.
            op2 = Operator(cells, beta, "naive")
            b2n = float(np.linalg.norm(op2.rhs(s.u, s.v, dt)))
            if b2n > 0:
                tol2 = min(1.0, self.tol * bnorm / b2n)
                op2, q, it2, _ = self._solve(cells, beta, "naive", dt, tol2)
                g2u, g2v = self._gradient(op2, q, s)
                s.u -= dt * bu * np.nan_to_num(g2u)
                s.v -= dt * bv * np.nan_to_num(g2v)
        # boundary faces of every coarse cell are final now; solve the interior
        # faces so each fine cell's divergence is as small as the cell's net
        # boundary flux allows (spread evenly)
        cell_net = [self._local(l, i, j, s, bu, bv, ustar, vstar, dt)
                    for (l, i, j) in cells if l > 0]
        s._enforce_walls()

        div = (s.u[1:, :] - s.u[:-1, :]) + (s.v[:, 1:] - s.v[:, :-1])
        lev = np.zeros((nx, ny), int)
        for (l, i, j) in cells:
            lev[i << l:(i + 1) << l, j << l:(j + 1) << l] = l
        fine = lev == 0
        self.log.append({
            "unknowns": len(cells), "uniform": nx * ny, "iters": int(it),
            "iters_cleanup": int(it2),
            "div_fine": float(np.abs(div[fine]).max()) if fine.any() else 0.0,
            "div_coarse": float(np.abs(div[~fine]).max()) if (~fine).any() else 0.0,
            "net_coarse": float(max(cell_net)) if cell_net else 0.0})

    def _solve(self, cells, beta, stencil, dt, tol):
        """Assemble on cells with the given junction stencil and solve for the
        current sim face fields; returns the fine operator, pressure,
        iterations and the norm of the right-hand side."""
        s = self.sim
        nx, ny = s.nx, s.ny
        if self.multigrid:
            # every hierarchy level lifts each cell one level, so after k
            # coarsenings the coarsest cells sit at self.levels + k; a cell
            # at level L fits the domain only if 2^L divides both grid sizes
            depth = 0
            while nx % (2 << depth) == 0 and ny % (2 << depth) == 0:
                depth += 1
            mg = Multigrid(cells, beta, max_levels=depth - self.levels + 1,
                           stencil=stencil)
            op = mg.ops[0]
            b = op.rhs(s.u, s.v, dt)
            p, it = mg.solve(b, tol=tol)
        else:
            op = Operator(cells, beta, stencil)
            b = op.rhs(s.u, s.v, dt)
            p, it = op.solve(b, tol=tol)
        return op, p, it, float(np.linalg.norm(b - b.mean()))

    @staticmethod
    def _gradient(op, p, s):
        """Pressure gradient (Delta / d) on every uniform face that lies on a
        face between two adaptive cells; NaN elsewhere, 0 on the walls."""
        gu = np.full(s.u.shape, np.nan)
        gv = np.full(s.v.shape, np.nan)
        for ent, _c, axis, area, kind, fx, fy in op.faces:
            delta = sum(w * p[idx] for idx, w in ent)
            g = delta / (area if kind == "same" else 1.5 * area)
            if axis == 0:
                gu[int(round(fx)), (fy - area / 2.0 + np.arange(area)).astype(int)] = g
            else:
                gv[(fx - area / 2.0 + np.arange(area)).astype(int), int(round(fy))] = g
        gu[0, :] = gu[-1, :] = 0.0
        gv[:, 0] = gv[:, -1] = 0.0
        return gu, gv

    @staticmethod
    def _local(l, i, j, s, bu, bv, ustar, vstar, dt):
        """Fine Neumann solve inside one coarse cell, boundary faces fixed."""
        sz = 1 << l
        x0, y0 = i * sz, j * sz
        n = sz * sz
        idx = lambda a, b: a * sz + b          # noqa: E731
        L = np.zeros((n, n))
        rhs = np.zeros(n)
        # divergence of the boundary-corrected field with interior faces at u*
        uu = s.u[x0:x0 + sz + 1, y0:y0 + sz].copy()
        vv = s.v[x0:x0 + sz, y0:y0 + sz + 1].copy()
        uu[1:-1, :] = ustar[x0 + 1:x0 + sz, y0:y0 + sz]
        vv[:, 1:-1] = vstar[x0:x0 + sz, y0 + 1:y0 + sz]
        d0 = (uu[1:, :] - uu[:-1, :]) + (vv[:, 1:] - vv[:, :-1])
        net = float(d0.sum())
        rhs[:] = -(d0 - net / n).ravel() / dt
        for a in range(sz):
            for b in range(sz):
                if a + 1 < sz:
                    c = bu[x0 + a + 1, y0 + b]
                    p_, q_ = idx(a, b), idx(a + 1, b)
                    L[p_, p_] += c
                    L[q_, q_] += c
                    L[p_, q_] -= c
                    L[q_, p_] -= c
                if b + 1 < sz:
                    c = bv[x0 + a, y0 + b + 1]
                    p_, q_ = idx(a, b), idx(a, b + 1)
                    L[p_, p_] += c
                    L[q_, q_] += c
                    L[p_, q_] -= c
                    L[q_, p_] -= c
        # an interior face moved by -dt beta (q_b - q_a) raises cell a's
        # divergence by dt beta (q_a - q_b) and lowers b's by the same, so with
        # L the graph Laplacian the divergence after is d0 + dt L q; solving
        # L q = -(d0 - net/n) / dt leaves every fine cell at net/n.
        q = np.linalg.lstsq(L, rhs, rcond=None)[0]
        q = q.reshape(sz, sz)
        s.u[x0 + 1:x0 + sz, y0:y0 + sz] = (ustar[x0 + 1:x0 + sz, y0:y0 + sz] -
                                           dt * bu[x0 + 1:x0 + sz, y0:y0 + sz] *
                                           (q[1:, :] - q[:-1, :]))
        s.v[x0:x0 + sz, y0 + 1:y0 + sz] = (vstar[x0:x0 + sz, y0 + 1:y0 + sz] -
                                           dt * bv[x0:x0 + sz, y0 + 1:y0 + sz] *
                                           (q[:, 1:] - q[:, :-1]))
        return abs(net)
