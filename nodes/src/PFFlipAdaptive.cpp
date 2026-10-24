#include "PFFlipAdaptive.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <utility>

// The structure follows the reference module by module: band_target ->
// bandTarget, build_cells/_grade -> buildCells/grade, Operator -> buildOp,
// Multigrid -> Hierarchy + solve, and AdaptiveProjection.__call__ ->
// project. Orders that decide rounding
// (cell lists, face enumeration, restriction sums) mirror the reference so
// the lockstep stays tight; the CSR matvec sums each row in column order,
// which differs from numpy's scatter only at rounding level.

namespace PFFlip {
namespace Adaptive {
namespace {

struct Cell {
    int l, i, j, k;
};

struct Dims {
    int nx, ny, nz;
    bool flat;
    size_t ncell() const { return static_cast<size_t>(nx) * ny * nz; }
    size_t ic(int i, int j, int k) const { return (static_cast<size_t>(i) * ny + j) * nz + k; }
    size_t iu(int i, int j, int k) const { return (static_cast<size_t>(i) * ny + j) * nz + k; }
    size_t iv(int i, int j, int k) const { return (static_cast<size_t>(i) * (ny + 1) + j) * nz + k; }
    size_t iw(int i, int j, int k) const { return (static_cast<size_t>(i) * ny + j) * (nz + 1) + k; }
    // a level-l cell spans 2^l fine cells in x and y, and in z unless flat
    int zs(int l) const { return flat ? 1 : (1 << l); }
    int axes() const { return flat ? 2 : 3; }
};

// ---- refinement criterion (band_target) -----------------------------------
// Neighbor comparisons and growth wrap around the domain, as np.roll does in
// the reference. Kept for lockstep; a known quirk to retire in both places
// together: a band touching one wall also refines the opposite wall.

template <class Fn>
void forNeighbors(const Dims& D, int i, int j, int k, Fn&& fn) {
    auto wrap = [](int a, int n) { return a < 0 ? a + n : (a >= n ? a - n : a); };
    fn(D.ic(wrap(i - 1, D.nx), j, k));
    fn(D.ic(wrap(i + 1, D.nx), j, k));
    fn(D.ic(i, wrap(j - 1, D.ny), k));
    fn(D.ic(i, wrap(j + 1, D.ny), k));
    if (!D.flat) {
        fn(D.ic(i, j, wrap(k - 1, D.nz)));
        fn(D.ic(i, j, wrap(k + 1, D.nz)));
    }
}

std::vector<char> grown(const Dims& D, std::vector<char> m, int steps) {
    std::vector<char> g(m.size());
    for (int s = 0; s < steps; ++s) {
        for (int i = 0; i < D.nx; ++i)
            for (int j = 0; j < D.ny; ++j)
                for (int k = 0; k < D.nz; ++k) {
                    const size_t c = D.ic(i, j, k);
                    char on = m[c];
                    forNeighbors(D, i, j, k, [&](size_t nb) { on |= m[nb]; });
                    g[c] = on;
                }
        m.swap(g);
    }
    return m;
}

std::vector<int> bandTarget(const Dims& D, const Vec& frac, int band, int levels) {
    const size_t n = D.ncell();
    std::vector<char> iface(n, 0);
    for (int i = 0; i < D.nx; ++i)
        for (int j = 0; j < D.ny; ++j)
            for (int k = 0; k < D.nz; ++k) {
                const size_t c = D.ic(i, j, k);
                const bool me = frac[c] > 0.5;
                char on = frac[c] > 0.05 && frac[c] < 0.95;
                forNeighbors(D, i, j, k, [&](size_t nb) {
                    if ((frac[nb] > 0.5) != me) on = 1;
                });
                iface[c] = on;
            }
    const std::vector<char> nearBand = grown(D, iface, band);
    const std::vector<char> midBand = grown(D, nearBand, 2 * band);
    std::vector<int> target(n, levels);
    for (int i = 0; i < D.nx; ++i)
        for (int j = 0; j < D.ny; ++j)
            for (int k = 0; k < D.nz; ++k) {
                const size_t c = D.ic(i, j, k);
                if (midBand[c]) target[c] = 1;
                if (nearBand[c]) target[c] = 0;
                if (i < 2 || i >= D.nx - 2 || j < 2 || j >= D.ny - 2 ||
                    (!D.flat && (k < 2 || k >= D.nz - 2)))
                    target[c] = 0;
            }
    return target;
}

// ---- the graded tree (build_cells, _grade) ---------------------------------

void appendChildren(const Dims& D, const Cell& c, std::vector<Cell>& out) {
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int e = 0; e < (D.flat ? 1 : 2); ++e)
                out.push_back({c.l - 1, 2 * c.i + a, 2 * c.j + b, D.flat ? c.k : 2 * c.k + e});
}

template <class Fn>
void forFootprintCells(const Dims& D, const Cell& c, Fn&& fn) {
    const int s = 1 << c.l, zs = D.zs(c.l);
    for (int i = c.i * s; i < (c.i + 1) * s; ++i)
        for (int j = c.j * s; j < (c.j + 1) * s; ++j)
            for (int k = c.k * zs; k < (c.k + 1) * zs; ++k)
                fn(i, j, k);
}

// per fine cell, the index of the adaptive cell covering it
std::vector<int> ownerMap(const Dims& D, const std::vector<Cell>& cells) {
    std::vector<int> own(D.ncell(), -1);
    for (size_t n = 0; n < cells.size(); ++n)
        forFootprintCells(D, cells[n], [&](int i, int j, int k) {
            own[D.ic(i, j, k)] = static_cast<int>(n);
        });
    return own;
}

std::vector<Cell> grade(const Dims& D, std::vector<Cell> cells) {
    // 2:1 balance: a cell with a face neighbor two or more levels finer
    // splits; repeat to a fixpoint (splitting in place keeps the order)
    while (true) {
        const std::vector<int> own = ownerMap(D, cells);
        std::vector<char> split(cells.size(), 0);
        bool any = false;
        for (size_t n = 0; n < cells.size(); ++n) {
            const Cell& c = cells[n];
            if (c.l == 0) continue;
            const int s = 1 << c.l, zs = D.zs(c.l);
            const int x0 = c.i * s, y0 = c.j * s, z0 = c.k * zs;
            auto tooFine = [&](int i, int j, int k) {
                return cells[own[D.ic(i, j, k)]].l <= c.l - 2;
            };
            bool bad = false;
            for (int fi : {x0 - 1, x0 + s}) {
                if (fi < 0 || fi >= D.nx) continue;
                for (int j = y0; j < y0 + s && !bad; ++j)
                    for (int k = z0; k < z0 + zs && !bad; ++k) bad = tooFine(fi, j, k);
            }
            for (int fj : {y0 - 1, y0 + s}) {
                if (fj < 0 || fj >= D.ny) continue;
                for (int i = x0; i < x0 + s && !bad; ++i)
                    for (int k = z0; k < z0 + zs && !bad; ++k) bad = tooFine(i, fj, k);
            }
            if (!D.flat)
                for (int fk : {z0 - 1, z0 + zs}) {
                    if (fk < 0 || fk >= D.nz) continue;
                    for (int i = x0; i < x0 + s && !bad; ++i)
                        for (int j = y0; j < y0 + s && !bad; ++j) bad = tooFine(i, j, fk);
                }
            if (bad) { split[n] = 1; any = true; }
        }
        if (!any) return cells;
        std::vector<Cell> next;
        next.reserve(cells.size() * 2);
        for (size_t n = 0; n < cells.size(); ++n) {
            if (split[n]) appendChildren(D, cells[n], next);
            else next.push_back(cells[n]);
        }
        cells.swap(next);
    }
}

std::vector<Cell> buildCells(const Dims& D, int levels, const std::vector<int>& target) {
    std::vector<Cell> cells;
    const int nk = D.flat ? 1 : D.nz >> levels;
    for (int i = 0; i < D.nx >> levels; ++i)
        for (int j = 0; j < D.ny >> levels; ++j)
            for (int k = 0; k < nk; ++k) cells.push_back({levels, i, j, k});
    // coarse-to-fine sweep: split any cell whose footprint demands finer
    for (int cur = levels; cur > 0; --cur) {
        std::vector<Cell> next;
        next.reserve(cells.size() * 2);
        for (const Cell& c : cells) {
            if (c.l != cur) { next.push_back(c); continue; }
            int lo = 1 << 30;
            forFootprintCells(D, c, [&](int i, int j, int k) {
                lo = std::min(lo, target[D.ic(i, j, k)]);
            });
            if (lo < c.l) appendChildren(D, c, next);
            else next.push_back(c);
        }
        cells.swap(next);
    }
    return grade(D, std::move(cells));
}

// ---- the operator (Operator) -----------------------------------------------

struct Face {
    int a, b;          // cells on the negative and positive side: Delta = p[b] - p[a]
    int axis;          // 0 x, 1 y, 2 z
    int plane;         // fine index of the face plane along the axis
    int t0, n0;        // first transverse axis (y for x faces, else x): origin, extent
    int t1, n1;        // second transverse axis (z, or y for z faces): origin, extent
    double d;          // center distance along the axis: s, or 1.5 s at a junction
    double coef;       // beta * area / d
    // a variational junction face (en > 0) replaces p[b] - p[a] with the
    // weighted sum over Op::entIdx/entW[e0, e0 + en), in the reference's
    // entry order: fine cell, coarse cell, lateral cells
    int e0 = 0, en = 0;
};

// visit the uniform faces a face covers, by index into its axis's array
template <class Fn>
void forFootprint(const Dims& D, const Face& f, Fn&& fn) {
    for (int p = f.t0; p < f.t0 + f.n0; ++p)
        for (int q = f.t1; q < f.t1 + f.n1; ++q) {
            if (f.axis == 0) fn(D.iu(f.plane, p, q));
            else if (f.axis == 1) fn(D.iv(p, f.plane, q));
            else fn(D.iw(p, q, f.plane));
        }
}

struct Op {
    std::vector<Cell> cells;
    std::vector<Face> faces;
    std::vector<int> entIdx;    // variational junction entries (Face::e0, en)
    Vec entW;
    std::vector<size_t> start;
    std::vector<int> col;
    Vec val, diag, dinv;
    size_t size() const { return cells.size(); }
    void matvec(const Vec& x, Vec& y) const {
        const long long n = static_cast<long long>(cells.size());
#pragma omp parallel for schedule(static)
        for (long long r = 0; r < n; ++r) {
            double acc = 0.0;
            for (size_t e = start[r]; e < start[r + 1]; ++e) acc += val[e] * x[col[e]];
            y[r] = acc;
        }
    }
};

uint64_t cellKey(const Cell& c) {
    return (static_cast<uint64_t>(c.l) << 60) | (static_cast<uint64_t>(c.i) << 40) |
           (static_cast<uint64_t>(c.j) << 20) | static_cast<uint64_t>(c.k);
}

// Delta across a face: p[b] - p[a], or a variational junction's weighted sum
double delta(const Op& op, const Face& f, const Vec& p) {
    if (f.en == 0) return p[f.b] - p[f.a];
    double s = 0.0;
    for (int e = f.e0; e < f.e0 + f.en; ++e) s += op.entW[e] * p[op.entIdx[e]];
    return s;
}

// variational = false is the reference's stencil = "naive"; true its
// "variational" stencil (at nz > 1 with both tangent axes as lateral axes)
Op buildOp(const Dims& D, std::vector<Cell> cells, const Vec* beta[3], bool variational) {
    Op op;
    op.cells = std::move(cells);
    const int N = static_cast<int>(op.cells.size());
    // neighbors are found by exact (level, index) keys, as the reference's
    // dict: the coarser multigrid levels are not partitions (a parent whose
    // children were mixed leaves and subdivided cells collapses next to the
    // collapsed subdivided child, inside its own footprint), and the key
    // lookups define which of those overlapping cells connect
    std::unordered_map<uint64_t, int> index;
    index.reserve(op.cells.size() * 2);
    for (int n = 0; n < N; ++n) index.emplace(cellKey(op.cells[n]), n);
    auto find = [&](int l, const int at[3]) {
        if (at[0] < 0 || at[1] < 0 || at[2] < 0) return -1;
        const auto it = index.find(cellKey({l, at[0], at[1], at[2]}));
        return it == index.end() ? -1 : it->second;
    };
    // Operator's lateral(): the coarse cell K's neighbor in direction latDir
    // along latAxis - the same level, one coarser, or the children of the
    // same-level neighbor that face K (a pair at nz = 1, four in 3D: both
    // halves along every other subdivided axis). Fills (cell, fraction)
    // entries summing to 1 and the center distance along latAxis in fine
    // units; false at a wall (no lateral at all).
    auto lateral = [&](const Cell& K, int latDir, int latAxis, int idx[4], double frac[4],
                       int& n, double& dist) {
        const int kc[3] = {K.i, K.j, K.k};
        int at[3] = {kc[0], kc[1], kc[2]};
        at[latAxis] += latDir;
        const int same = find(K.l, at);
        if (same >= 0) {
            idx[0] = same, frac[0] = 1.0, n = 1, dist = static_cast<double>(1 << K.l);
            return true;
        }
        // one coarser, not K's own parent; the shifts keep the reference's
        // floor semantics at -1 (find rejects it); z stays put when flat
        int up[3], mine[3];
        for (int a = 0; a < 3; ++a) {
            const bool keep = a == 2 && D.flat;
            up[a] = keep ? at[a] : at[a] >> 1;
            mine[a] = keep ? kc[a] : kc[a] >> 1;
        }
        const int coarser = find(K.l + 1, up);
        if (coarser >= 0 && !(up[0] == mine[0] && up[1] == mine[1] && up[2] == mine[2])) {
            idx[0] = coarser, frac[0] = 1.0, n = 1, dist = 1.5 * (1 << K.l);
            return true;
        }
        const int fl = K.l - 1;
        if (fl < 0) return false;
        int base[3], vary[2], nv = 0;
        for (int a = 0; a < 3; ++a) {
            if (a == latAxis)
                base[a] = latDir > 0 ? 2 * (kc[a] + latDir) : 2 * kc[a] - 1;
            else if (a == 2 && D.flat)
                base[a] = kc[a];
            else
                base[a] = 2 * kc[a], vary[nv++] = a;
        }
        const int count = 1 << nv;
        for (int c = 0; c < count; ++c) {
            int ch[3] = {base[0], base[1], base[2]};
            for (int v = 0; v < nv; ++v) ch[vary[v]] += (c >> (nv - 1 - v)) & 1;
            idx[c] = find(fl, ch);
            if (idx[c] < 0) return false;
            frac[c] = 1.0 / count;
        }
        n = count;
        dist = 1.5 * (1 << fl);
        return true;
    };
    for (int me = 0; me < N; ++me) {
        const Cell& C = op.cells[me];
        const int s = 1 << C.l, zs = D.zs(C.l);
        const int idx[3] = {C.i, C.j, C.k};
        const int org[3] = {C.i * s, C.j * s, C.k * zs};
        const int ext[3] = {s, s, zs};
        for (int axis = 0; axis < D.axes(); ++axis) {
            const int ta = axis == 0 ? 1 : 0, tb = axis == 2 ? 1 : 2;
            const int lo = org[axis], hi = org[axis] + ext[axis];
            auto push = [&](int a, int b, int plane, double d) {
                op.faces.push_back({a, b, axis, plane, org[ta], ext[ta], org[tb], ext[tb], d, 0.0});
            };
            // the same-level neighbor on the positive side
            int at[3] = {idx[0], idx[1], idx[2]};
            at[axis] += 1;
            const int same = find(C.l, at);
            if (same >= 0) push(me, same, hi, s);
            // a one-level-coarser cell holding the fine cell just across either
            // side at my origin row: the junction sub-face is mine
            for (int side : {+1, -1}) {
                const int fi = side > 0 ? hi : lo - 1;
                int up[3];
                for (int a = 0; a < 3; ++a) {
                    const int fine = a == axis ? fi : org[a];
                    up[a] = (a == 2 && D.flat) ? fine : (fine < 0 ? -1 : fine >> (C.l + 1));
                }
                const int o = find(C.l + 1, up);
                if (o < 0) continue;
                if (side > 0) push(me, o, hi, 1.5 * s);
                else push(o, me, lo, 1.5 * s);
                if (!variational) continue;
                // the variational sample: the coarse side is read where the
                // fine cell's center sits, s/2 off the coarse axis along each
                // lateral axis (one at nz = 1, both tangents in 3D), toward
                // the lateral neighbor on the fine cell's side:
                //   p_K + sum_t w_t (p_lateral_t - p_K),  w_t = (s/2) / D_t
                // (3/4 coarse + 1/4 a same-level lateral in 2D; 1/2 + 1/4 +
                // 1/4 in 3D), exact for linear fields without the edge
                // neighbor. An axis with no lateral (a wall) drops its term;
                // with none at all the face stays the plain two-cell one, as
                // the reference
                int li[2][4], ln[2] = {0, 0}, nt = 0;
                double lf[2][4], wt[2];
                for (int t : {ta, tb}) {
                    if (t == 2 && D.flat) continue;
                    const int ofs = org[t];
                    const int latDir = ((ofs >> (C.l + 1)) << (C.l + 1)) != ofs ? 1 : -1;
                    double dist = 0.0;
                    if (!lateral(op.cells[o], latDir, t, li[nt], lf[nt], ln[nt], dist)) continue;
                    wt[nt++] = (static_cast<double>(s) / 2.0) / dist;
                }
                if (nt == 0) continue;
                double wsum = 0.0;
                for (int q = 0; q < nt; ++q) wsum += wt[q];
                const double sgn = static_cast<double>(side);
                Face& f = op.faces.back();
                f.e0 = static_cast<int>(op.entIdx.size());
                f.en = 2;
                op.entIdx.push_back(me), op.entW.push_back(-sgn);
                op.entIdx.push_back(o), op.entW.push_back((1.0 - wsum) * sgn);
                for (int q = 0; q < nt; ++q)
                    for (int e = 0; e < ln[q]; ++e, ++f.en)
                        op.entIdx.push_back(li[q][e]), op.entW.push_back(wt[q] * lf[q][e] * sgn);
            }
        }
    }
    // coefficients: the face-averaged uniform beta times area over distance
    for (Face& f : op.faces) {
        double sum = 0.0;
        int n = 0;
        forFootprint(D, f, [&](size_t q) { sum += (*beta[f.axis])[q]; ++n; });
        const double area = static_cast<double>(f.n0) * f.n1;
        f.coef = (sum / n) * area / f.d;
    }
    // CSR, rows sorted by column, duplicates merged
    // (a variational face adds coef * w_x * w_y for every entry pair: the
    // energy coef * Delta^2, so the operator stays symmetric with zero row sums)
    std::vector<std::vector<std::pair<int, double>>> rows(N);
    for (const Face& f : op.faces) {
        if (f.en > 0) {
            for (int x = f.e0; x < f.e0 + f.en; ++x)
                for (int y = f.e0; y < f.e0 + f.en; ++y)
                    rows[op.entIdx[x]].push_back({op.entIdx[y], f.coef * op.entW[x] * op.entW[y]});
            continue;
        }
        rows[f.a].push_back({f.a, f.coef});
        rows[f.a].push_back({f.b, -f.coef});
        rows[f.b].push_back({f.b, f.coef});
        rows[f.b].push_back({f.a, -f.coef});
    }
    op.start.assign(N + 1, 0);
    op.diag.assign(N, 0.0);
    for (int r = 0; r < N; ++r) {
        auto& row = rows[r];
        std::stable_sort(row.begin(), row.end(),
                         [](const auto& x, const auto& y) { return x.first < y.first; });
        for (size_t e = 0; e < row.size(); ++e) {
            if (!op.col.empty() && op.start[r] < op.col.size() && op.col.back() == row[e].first)
                op.val.back() += row[e].second;
            else {
                op.col.push_back(row[e].first);
                op.val.push_back(row[e].second);
            }
            if (row[e].first == r) op.diag[r] += row[e].second;
        }
        op.start[r + 1] = op.col.size();
    }
    op.dinv.assign(N, 0.0);
    for (int r = 0; r < N; ++r)
        op.dinv[r] = op.diag[r] > 0.0 ? 1.0 / std::max(op.diag[r], 1e-300) : 0.0;
    return op;
}

// b = G^T (A u*): each face's flux, area times the mean uniform velocity it
// covers, leaves its negative-side cell and enters its positive-side cell (a
// variational face spreads it over its entries with the gradient's weights)
Vec rhsOf(const Dims& D, const Op& op, const Vec* vel[3], double dt) {
    Vec b(op.size(), 0.0);
    for (const Face& f : op.faces) {
        double sum = 0.0;
        int n = 0;
        forFootprint(D, f, [&](size_t q) { sum += (*vel[f.axis])[q]; ++n; });
        const double ustar = sum / n;
        const double area = static_cast<double>(f.n0) * f.n1;
        if (f.en > 0) {
            for (int e = f.e0; e < f.e0 + f.en; ++e)
                b[op.entIdx[e]] += area * ustar * op.entW[e] / dt;
            continue;
        }
        b[f.a] += area * ustar * -1.0 / dt;
        b[f.b] += area * ustar * 1.0 / dt;
    }
    return b;
}

// ---- multigrid on the collapsed tree (Multigrid) ----------------------------

struct Hierarchy {
    std::vector<Op> ops;
    std::vector<std::vector<int>> parent;   // fine index -> coarse index, per level
};

// every level is assembled with the same stencil, as the reference
Hierarchy buildHierarchy(const Dims& D, std::vector<Cell> cells, const Vec* beta[3],
                         int maxLevels, size_t coarsest, bool variational) {
    Hierarchy h;
    h.ops.push_back(buildOp(D, std::move(cells), beta, variational));
    while (static_cast<int>(h.ops.size()) < maxLevels && h.ops.back().size() > coarsest) {
        // every cell collapses into its parent, first occurrence keeps the order
        const std::vector<Cell>& fine = h.ops.back().cells;
        std::unordered_map<uint64_t, int> index;
        std::vector<Cell> coarse;
        std::vector<int> parent(fine.size());
        for (size_t n = 0; n < fine.size(); ++n) {
            const Cell& c = fine[n];
            const Cell pc{c.l + 1, c.i >> 1, c.j >> 1, D.flat ? c.k : c.k >> 1};
            const uint64_t key = cellKey(pc);
            auto it = index.find(key);
            if (it == index.end()) {
                it = index.emplace(key, static_cast<int>(coarse.size())).first;
                coarse.push_back(pc);
            }
            parent[n] = it->second;
        }
        h.parent.push_back(std::move(parent));
        h.ops.push_back(buildOp(D, std::move(coarse), beta, variational));
    }
    return h;
}

void smooth(const Op& op, Vec& x, const Vec& b, int sweeps, Vec& tmp) {
    const long long n = static_cast<long long>(op.size());
    for (int s = 0; s < sweeps; ++s) {
        op.matvec(x, tmp);
#pragma omp parallel for schedule(static)
        for (long long c = 0; c < n; ++c) x[c] += 0.8 * op.dinv[c] * (b[c] - tmp[c]);
    }
}

Vec vcycle(const Hierarchy& h, size_t k, const Vec& b) {
    const Op& op = h.ops[k];
    Vec x(op.size(), 0.0), tmp(op.size());
    if (k + 1 == h.ops.size()) {
        smooth(op, x, b, 40, tmp);
        return x;
    }
    smooth(op, x, b, 2, tmp);
    op.matvec(x, tmp);
    Vec bc(h.ops[k + 1].size(), 0.0);
    const std::vector<int>& par = h.parent[k];
    for (size_t c = 0; c < op.size(); ++c) bc[par[c]] += b[c] - tmp[c];
    const Vec xc = vcycle(h, k + 1, bc);
    for (size_t c = 0; c < op.size(); ++c) x[c] += xc[par[c]];
    smooth(op, x, b, 2, tmp);
    return x;
}

double dot(const Vec& a, const Vec& b) {
    double s = 0.0;   // serial: the summation order is part of the result
    for (size_t c = 0; c < a.size(); ++c) s += a[c] * b[c];
    return s;
}

// flexible (Polak-Ribiere) PCG with the V-cycle; mean-free right-hand side
// and solution (all-Neumann walls), as Multigrid.solve
Vec solve(const Hierarchy& h, Vec b, double tol, int maxIter, int& iterations, double& relRes) {
    const Op& op = h.ops[0];
    const size_t n = op.size();
    double mean = 0.0;
    for (double x : b) mean += x;
    mean /= static_cast<double>(n);
    for (double& x : b) x -= mean;
    Vec x(n, 0.0), r = b, q(n), rPrev;
    Vec z = vcycle(h, 0, r), d = z;
    double rz = dot(r, z);
    const double r0 = std::max(std::sqrt(dot(r, r)), 1e-300);
    int it = 0;
    for (it = 0; it < maxIter; ++it) {
        if (std::sqrt(dot(r, r)) / r0 < tol) {
            // Certification includes the gauge normalization performed on return.
            Vec candidate = x;
            double xm = 0.0;
            for (double v : candidate) xm += v;
            xm /= static_cast<double>(n);
            for (double& v : candidate) v -= xm;
            op.matvec(candidate, q);
            Vec actual(n);
            for (size_t c = 0; c < n; ++c) actual[c] = b[c] - q[c];
            const double certified = std::sqrt(dot(actual, actual)) / r0;
            if (!std::isfinite(certified) || certified <= tol) break;
            x.swap(candidate); r.swap(actual);
            z = vcycle(h, 0, r); d = z; rz = dot(r, z);
            if (std::getenv("PFFLIP_RESIDUAL_RESTART_TRACE"))
                std::fprintf(stderr,"residual-restart: adaptive iteration=%d true=%.17g tolerance=%.17g\n",it,certified,tol);
        }
        op.matvec(d, q);
        const double dq = dot(d, q);
        if (dq == 0.0) break;
        const double a = rz / dq;
        for (size_t c = 0; c < n; ++c) x[c] += a * d[c];
        rPrev = r;
        for (size_t c = 0; c < n; ++c) r[c] -= a * q[c];
        z = vcycle(h, 0, r);
        double rzNew = 0.0;
        for (size_t c = 0; c < n; ++c) rzNew += (r[c] - rPrev[c]) * z[c];
        rzNew = std::max(0.0, rzNew);
        const double betaCG = rzNew / std::max(rz, 1e-300);
        for (size_t c = 0; c < n; ++c) d[c] = z[c] + betaCG * d[c];
        rz = dot(r, z);
        if (rz == 0.0) break;
    }
    if (it == maxIter) it = std::max(0, maxIter - 1); // reference count; zero budget is not a failure code
    iterations = it;
    mean = 0.0;
    for (double v : x) mean += v;
    mean /= static_cast<double>(n);
    for (double& v : x) v -= mean;
    op.matvec(x, q);
    for (size_t c = 0; c < n; ++c) r[c] = b[c] - q[c];
    relRes = std::sqrt(dot(r, r)) / r0;
    return x;
}

// ---- inside a coarse cell (AdaptiveProjection._local) -----------------------

// Dense Cholesky solve of the SPD system M y = r (n small); M row-major.
void choleskySolve(std::vector<double>& M, std::vector<double>& r, int n) {
    for (int c = 0; c < n; ++c) {
        double s = M[c * n + c];
        for (int k = 0; k < c; ++k) s -= M[c * n + k] * M[c * n + k];
        const double dg = std::sqrt(std::max(s, 1e-300));
        M[c * n + c] = dg;
        for (int r2 = c + 1; r2 < n; ++r2) {
            double t = M[r2 * n + c];
            for (int k = 0; k < c; ++k) t -= M[r2 * n + k] * M[c * n + k];
            M[r2 * n + c] = t / dg;
        }
    }
    for (int c = 0; c < n; ++c) {
        double t = r[c];
        for (int k = 0; k < c; ++k) t -= M[c * n + k] * r[k];
        r[c] = t / M[c * n + c];
    }
    for (int c = n - 1; c >= 0; --c) {
        double t = r[c];
        for (int k = c + 1; k < n; ++k) t -= M[k * n + c] * r[k];
        r[c] = t / M[c * n + c];
    }
}

// Boundary faces of the coarse cell are final; solve the fine Neumann problem
// for its interior faces (starting from u*) so every fine cell's divergence
// equals the cell's net boundary flux spread evenly. The pressure is defined
// up to a constant; pinning the first cell gives the reference's minimum-norm
// lstsq solution plus a constant, so every face difference agrees.
void localFill(const Dims& D, const Cell& C, double dt, Vec* vel[3], const Vec* star[3],
               const Vec* beta[3]) {
    const int s = 1 << C.l, zs = D.zs(C.l);
    const int x0 = C.i * s, y0 = C.j * s, z0 = C.k * zs;
    const int n = s * s * zs;
    auto idx = [&](int a, int b, int c) { return (a * s + b) * zs + c; };
    Vec& u = *vel[0];
    Vec& v = *vel[1];
    Vec& w = *vel[2];
    // divergence with boundary faces corrected and interior faces at u*
    auto uAt = [&](int a, int b, int c) {
        const size_t q = D.iu(x0 + a, y0 + b, z0 + c);
        return (a == 0 || a == s) ? u[q] : (*star[0])[q];
    };
    auto vAt = [&](int a, int b, int c) {
        const size_t q = D.iv(x0 + a, y0 + b, z0 + c);
        return (b == 0 || b == s) ? v[q] : (*star[1])[q];
    };
    auto wAt = [&](int a, int b, int c) {
        const size_t q = D.iw(x0 + a, y0 + b, z0 + c);
        return (c == 0 || c == zs) ? w[q] : (*star[2])[q];
    };
    Vec d0(n);
    double net = 0.0;
    for (int a = 0; a < s; ++a)
        for (int b = 0; b < s; ++b)
            for (int c = 0; c < zs; ++c) {
                double dv = (uAt(a + 1, b, c) - uAt(a, b, c)) + (vAt(a, b + 1, c) - vAt(a, b, c));
                if (!D.flat) dv += wAt(a, b, c + 1) - wAt(a, b, c);
                d0[idx(a, b, c)] = dv;
            }
    for (int c = 0; c < n; ++c) net += d0[c];
    // graph Laplacian over the interior faces; divergence after = d0 + dt L q
    std::vector<double> L(static_cast<size_t>(n) * n, 0.0);
    auto link = [&](int p, int q, double coef) {
        L[p * n + p] += coef;
        L[q * n + q] += coef;
        L[p * n + q] -= coef;
        L[q * n + p] -= coef;
    };
    for (int a = 0; a < s; ++a)
        for (int b = 0; b < s; ++b)
            for (int c = 0; c < zs; ++c) {
                if (a + 1 < s) link(idx(a, b, c), idx(a + 1, b, c), (*beta[0])[D.iu(x0 + a + 1, y0 + b, z0 + c)]);
                if (b + 1 < s) link(idx(a, b, c), idx(a, b + 1, c), (*beta[1])[D.iv(x0 + a, y0 + b + 1, z0 + c)]);
                if (c + 1 < zs) link(idx(a, b, c), idx(a, b, c + 1), (*beta[2])[D.iw(x0 + a, y0 + b, z0 + c + 1)]);
            }
    // pin q[0] = 0 and solve the reduced SPD system
    const int m = n - 1;
    std::vector<double> M(static_cast<size_t>(m) * m), r(m);
    for (int p = 1; p < n; ++p) {
        r[p - 1] = -(d0[p] - net / n) / dt;
        for (int q = 1; q < n; ++q) M[(p - 1) * m + (q - 1)] = L[p * n + q];
    }
    choleskySolve(M, r, m);
    Vec q(n, 0.0);
    for (int p = 1; p < n; ++p) q[p] = r[p - 1];
    for (int a = 1; a < s; ++a)
        for (int b = 0; b < s; ++b)
            for (int c = 0; c < zs; ++c) {
                const size_t f = D.iu(x0 + a, y0 + b, z0 + c);
                u[f] = (*star[0])[f] - dt * (*beta[0])[f] * (q[idx(a, b, c)] - q[idx(a - 1, b, c)]);
            }
    for (int a = 0; a < s; ++a)
        for (int b = 1; b < s; ++b)
            for (int c = 0; c < zs; ++c) {
                const size_t f = D.iv(x0 + a, y0 + b, z0 + c);
                v[f] = (*star[1])[f] - dt * (*beta[1])[f] * (q[idx(a, b, c)] - q[idx(a, b - 1, c)]);
            }
    for (int a = 0; a < s; ++a)
        for (int b = 0; b < s; ++b)
            for (int c = 1; c < zs; ++c) {
                const size_t f = D.iw(x0 + a, y0 + b, z0 + c);
                w[f] = (*star[2])[f] - dt * (*beta[2])[f] * (q[idx(a, b, c)] - q[idx(a, b, c - 1)]);
            }
}

// faces between adaptive cells take their face's gradient, each uniform face
// with its own beta (AdaptiveProjection._gradient and the update after it)
void applyGradient(const Dims& D, const Op& op, const Vec& p, double dt, Vec* out[3],
                   const Vec* beta[3]) {
    for (const Face& f : op.faces) {
        const double g = delta(op, f, p) / f.d;
        Vec& field = *out[f.axis];
        const Vec& b = *beta[f.axis];
        forFootprint(D, f, [&](size_t q) { field[q] -= dt * b[q] * g; });
    }
}

} // namespace

Result project(int nx, int ny, int nz, double dt, const Settings& settings, Junction junction,
               const Vec& liquidFraction, Vec& u, Vec& v, Vec& w,
               const Vec& bu, const Vec& bv, const Vec& bw) {
    const Dims D{nx, ny, nz, nz == 1};
    // how many times every coarsened axis halves evenly; leaf levels and the
    // hierarchy depth must fit inside it (the reference asserts divisibility)
    int depth = 0;
    while (nx % (2 << depth) == 0 && ny % (2 << depth) == 0 &&
           (D.flat || nz % (2 << depth) == 0))
        ++depth;
    const int levels = std::clamp(settings.levels, 0, depth);
    // PFFLIP_PROFILE=1 splits the step node's solve phase (spike 1) into the
    // adaptive solve's own parts, printed to stderr as a PFPROF_ADAPTIVE line
    static const bool prof = [] {
        const char* e = std::getenv("PFFLIP_PROFILE");
        return e && e[0] == '1';
    }();
    using Clock = std::chrono::steady_clock;
    Clock::time_point mark = Clock::now();
    double msTree = 0.0, msHier = 0.0, msSolve = 0.0, msGrad = 0.0, msClean = 0.0, msFill = 0.0;
    auto lap = [&](double& acc) {
        if (!prof) return;
        const Clock::time_point now = Clock::now();
        acc += std::chrono::duration<double, std::milli>(now - mark).count();
        mark = now;
    };

    const std::vector<int> target = bandTarget(D, liquidFraction, settings.band, levels);
    std::vector<Cell> cells = buildCells(D, levels, target);
    lap(msTree);

    const bool cleanup = junction == Junction::Cleanup;
    const Vec* beta[3] = {&bu, &bv, &bw};
    const int maxLevels = depth - levels + 1;
    const Hierarchy h = buildHierarchy(D, std::move(cells), beta, maxLevels, 64, cleanup);
    const Op& op = h.ops[0];
    lap(msHier);
    // PFFLIP_ADAPTIVE_DUMP=<path> writes the leaf operator for inspection: one
    // "C l i j k" line per cell, then one "F axis a b d en [cell weight]..."
    // line per face (en = 0: Delta = p[b] - p[a]); gates/validate_adaptive_3d.py
    // checks every face's exactness for linear pressure from it
    if (const char* path = std::getenv("PFFLIP_ADAPTIVE_DUMP")) {
        if (FILE* fp = std::fopen(path, "w")) {
            for (const Cell& c : op.cells) std::fprintf(fp, "C %d %d %d %d\n", c.l, c.i, c.j, c.k);
            for (const Face& f : op.faces) {
                std::fprintf(fp, "F %d %d %d %.17g %d", f.axis, f.a, f.b, f.d, f.en);
                for (int e = f.e0; e < f.e0 + f.en; ++e)
                    std::fprintf(fp, " %d %.17g", op.entIdx[e], op.entW[e]);
                std::fprintf(fp, "\n");
            }
            std::fclose(fp);
        }
    }

    const Vec* vel[3] = {&u, &v, &w};
    Result res;
    res.levels = levels;
    res.unknowns = op.size();
    const Vec b = rhsOf(D, op, vel, dt);
    const Vec p = solve(h, b, settings.tolerance, settings.maxIterations, res.iterations,
                        res.residual);
    if (!std::isfinite(res.residual) || res.residual > settings.tolerance) {
        res.converged = false;
        res.failureStage = "adaptive main";
        res.failureResidual = res.residual;
        res.failureTolerance = settings.tolerance;
        return res;
    }
    lap(msSolve);

    // faces inside coarse cells keep u* until the local fills below
    const Vec uStar = u, vStar = v, wStar = w;
    Vec* out[3] = {&u, &v, &w};
    applyGradient(D, op, p, dt, out, beta);
    lap(msGrad);
    Vec q;
    if (cleanup) {
        // the variational fluxes balance only in the weighted sense; a second,
        // conservative (naive-stencil) projection of the corrected field removes
        // the leftover junction dipoles. Its tolerance makes its absolute
        // residual match the main solve's.
        double mean = 0.0;
        for (double x : b) mean += x;
        mean /= static_cast<double>(b.size());
        double bnorm = 0.0;
        for (double x : b) bnorm += (x - mean) * (x - mean);
        bnorm = std::sqrt(bnorm);
        const Hierarchy h2 = buildHierarchy(D, op.cells, beta, maxLevels, 64, false);
        const Vec b2 = rhsOf(D, h2.ops[0], vel, dt);
        const double b2n = std::sqrt(dot(b2, b2));
        if (b2n > 0.0) {
            const double tol2 = std::min(1.0, settings.tolerance * bnorm / b2n);
            q = solve(h2, b2, tol2, settings.maxIterations, res.cleanupIterations, res.cleanupResidual);
            if (!std::isfinite(res.cleanupResidual) || res.cleanupResidual > tol2) {
                res.converged = false;
                res.failureStage = "adaptive cleanup";
                res.failureResidual = res.cleanupResidual;
                res.failureTolerance = tol2;
                return res;
            }
            applyGradient(D, h2.ops[0], q, dt, out, beta);
        }
    }
    lap(msClean);
    const Vec* star[3] = {&uStar, &vStar, &wStar};
    const long long nc = static_cast<long long>(op.size());
#pragma omp parallel for schedule(dynamic, 16)
    for (long long c = 0; c < nc; ++c)
        if (op.cells[c].l > 0) localFill(D, op.cells[c], dt, out, star, beta);

    res.pressure.assign(D.ncell(), 0.0);
    for (size_t n = 0; n < op.size(); ++n) {
        const double pn = q.empty() ? p[n] : p[n] + q[n];
        forFootprintCells(D, op.cells[n], [&](int i, int j, int k) {
            res.pressure[D.ic(i, j, k)] = pn;
        });
    }
    lap(msFill);
    if (prof)
        std::fprintf(stderr,
                     "PFPROF_ADAPTIVE cells=%zu of %zu levels=%d iters=%d+%d tree=%.1f hier=%.1f "
                     "solve=%.1f grad=%.1f cleanup=%.1f fill=%.1f\n",
                     op.size(), D.ncell(), levels, res.iterations, res.cleanupIterations, msTree,
                     msHier, msSolve, msGrad, msClean, msFill);
    return res;
}

} // namespace Adaptive
} // namespace PFFlip
