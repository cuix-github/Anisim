//-
// Level-set solids (mesh colliders) for step_3d. The solid arrives as its
// signed distance at the grid nodes, in cells, negative inside: what the
// shipped mesh_to_level_set plus one sample_volume at the node positions
// hands the graph. Two jobs live here, both inside the step: how much of
// each face the solid leaves open (the pressure solve's weights), and pushing
// a particle that ended inside back out. reference2d/pfflip2d.py is the spec:
// on a flat grid (nz = 1) every function below repeats its 2D arithmetic
// operation for operation, which is what the nz = 1 lockstep relies on.
//+

#ifndef PF_FLIP_SOLIDS_H
#define PF_FLIP_SOLIDS_H

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace PFFlip {
namespace Solids {

// Cut cells close a face less open than this. The solve weights a face by its
// open part, so it barely sees a sliver, and the velocity it leaves there is
// not pinned down: where the water surface meets a curved wall, faces under 1%
// open came out at up to 3 cells/s from water at rest, and the particles beside
// them picked it up. The reference's Sim.SOLID_MIN_OPEN is the same constant.
constexpr double kMinOpen = 0.01;

// The part of a face edge where the distance is positive, linear between the
// edge's ends a and b (the reference's _open_fraction).
inline double edgeOpen(double a, double b) {
    if (a > 0.0 && b > 0.0) return 1.0;
    if (a <= 0.0 && b <= 0.0) return 0.0;
    const double t = a / (a - b);
    return a > 0.0 ? t : 1.0 - t;
}

// The part of a triangle where a distance linear over it is positive, from
// its three corner values.
inline double triangleOpen(double a, double b, double c) {
    const int pos = (a > 0.0) + (b > 0.0) + (c > 0.0);
    if (pos == 3) return 1.0;
    if (pos == 0) return 0.0;
    if (pos == 1) {
        // the open part is the corner triangle at the one positive vertex
        const double p = a > 0.0 ? a : (b > 0.0 ? b : c);
        const double n1 = a > 0.0 ? b : a;
        const double n2 = c > 0.0 ? b : c;
        return (p / (p - n1)) * (p / (p - n2));
    }
    // two positive: the closed part is the corner triangle at the third vertex
    const double m = a <= 0.0 ? a : (b <= 0.0 ? b : c);
    const double p1 = a <= 0.0 ? b : a;
    const double p2 = c <= 0.0 ? b : c;
    return 1.0 - (m / (m - p1)) * (m / (m - p2));
}

// The open area of a unit square face from its corner values, (s, t) =
// (0,0) a, (1,0) b, (0,1) c, (1,1) d: four triangles around the centre, the
// centre taking the corners' mean. Exact whenever the distance is linear
// across the face (a planar surface), and equal to edgeOpen for a distance
// that varies along one side only.
inline double faceOpen(double a, double b, double c, double d) {
    const double m = 0.25 * ((a + b) + (c + d));
    return 0.25 * (triangleOpen(a, b, m) + triangleOpen(b, d, m) +
                   triangleOpen(d, c, m) + triangleOpen(c, a, m));
}

// The distance and its gradient at a point. Flat grids interpolate the k = 0
// node layer bilinearly, as the reference's _phi_at; otherwise trilinear.
// Node layout: (nx+1) x (ny+1) x (nz+1), row-major, index (i*(ny+1)+j)*(nz+1)+k.
template <class Arr>
inline double phiAt(const Arr& P, int nx, int ny, int nz, bool flat, double px,
                    double py, double pz, double& gx, double& gy, double& gz) {
    const size_t sj = static_cast<size_t>(nz) + 1;
    const size_t si = (static_cast<size_t>(ny) + 1) * sj;
    const double qx = std::clamp(px, 0.0, static_cast<double>(nx));
    const double qy = std::clamp(py, 0.0, static_cast<double>(ny));
    const int i = std::min(static_cast<int>(std::floor(qx)), nx - 1);
    const int j = std::min(static_cast<int>(std::floor(qy)), ny - 1);
    const double fx = qx - i, fy = qy - j;
    if (flat) {
        const size_t o = i * si + j * sj;
        const double p00 = P[o], p10 = P[o + si], p01 = P[o + sj], p11 = P[o + si + sj];
        gx = (1 - fy) * (p10 - p00) + fy * (p11 - p01);
        gy = (1 - fx) * (p01 - p00) + fx * (p11 - p10);
        gz = 0.0;
        return ((1 - fx) * (1 - fy) * p00 + fx * (1 - fy) * p10 +
                (1 - fx) * fy * p01 + fx * fy * p11);
    }
    const double qz = std::clamp(pz, 0.0, static_cast<double>(nz));
    const int k = std::min(static_cast<int>(std::floor(qz)), nz - 1);
    const double fz = qz - k;
    const size_t o = i * si + j * sj + k;
    double c[2][2][2];
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int e = 0; e < 2; ++e) c[a][b][e] = P[o + a * si + b * sj + e];
    const double wx[2] = {1 - fx, fx}, wy[2] = {1 - fy, fy}, wz[2] = {1 - fz, fz};
    double phi = 0.0;
    gx = gy = gz = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            gx += wy[a] * wz[b] * (c[1][a][b] - c[0][a][b]);
            gy += wx[a] * wz[b] * (c[a][1][b] - c[a][0][b]);
            gz += wx[a] * wy[b] * (c[a][b][1] - c[a][b][0]);
            for (int e = 0; e < 2; ++e) phi += wx[a] * wy[b] * wz[e] * c[a][b][e];
        }
    return phi;
}

// The distance alone, for the inside test that every particle takes: the same
// arithmetic as phiAt's value, without the gradient.
template <class Arr>
inline double phiValue(const Arr& P, int nx, int ny, int nz, bool flat, double px,
                       double py, double pz) {
    const size_t sj = static_cast<size_t>(nz) + 1;
    const size_t si = (static_cast<size_t>(ny) + 1) * sj;
    const double qx = std::clamp(px, 0.0, static_cast<double>(nx));
    const double qy = std::clamp(py, 0.0, static_cast<double>(ny));
    const int i = std::min(static_cast<int>(std::floor(qx)), nx - 1);
    const int j = std::min(static_cast<int>(std::floor(qy)), ny - 1);
    const double fx = qx - i, fy = qy - j;
    if (flat) {
        const size_t o = i * si + j * sj;
        const double p00 = P[o], p10 = P[o + si], p01 = P[o + sj], p11 = P[o + si + sj];
        return ((1 - fx) * (1 - fy) * p00 + fx * (1 - fy) * p10 +
                (1 - fx) * fy * p01 + fx * fy * p11);
    }
    const double qz = std::clamp(pz, 0.0, static_cast<double>(nz));
    const int k = std::min(static_cast<int>(std::floor(qz)), nz - 1);
    const double fz = qz - k;
    const size_t o = i * si + j * sj + k;
    const double wx[2] = {1 - fx, fx}, wy[2] = {1 - fy, fy}, wz[2] = {1 - fz, fz};
    double phi = 0.0;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b)
            for (int e = 0; e < 2; ++e)
                phi += wx[a] * wy[b] * wz[e] * static_cast<double>(P[o + a * si + b * sj + e]);
    return phi;
}

// Push a point that is inside the solid back out: up to three Newton steps to
// a distance of +0.01 along the gradient, each at most 1.5 cells (one unit
// step on a true distance; longer where the interpolated field is no distance,
// as in a box's corner cell where its slope can be a tenth of one). With vel
// given, each step also removes the velocity's component along the normal.
// Returns whether the point started inside. Mirrors the reference's
// _push_out (whose passes skip particles already outside).
template <class Arr>
inline bool pushOut(const Arr& P, int nx, int ny, int nz, bool flat, double& px,
                    double& py, double& pz, double* vx, double* vy, double* vz) {
    // most points are outside: settle that without the gradient
    if (!(phiValue(P, nx, ny, nz, flat, px, py, pz) < 0.0)) return false;
    bool moved = false;
    for (int pass = 0; pass < 3; ++pass) {
        double gx, gy, gz;
        const double phi = phiAt(P, nx, ny, nz, flat, px, py, pz, gx, gy, gz);
        if (!(phi < 0.0)) break;
        moved = true;
        const double g2 = std::max(flat ? gx * gx + gy * gy : gx * gx + gy * gy + gz * gz,
                                   1e-12);
        const double rg = std::sqrt(g2);
        const double nxv = gx / rg, nyv = gy / rg, nzv = gz / rg;
        const double c = -(phi - 0.01);
        const double sx = c * gx / g2, sy = c * gy / g2, sz = flat ? 0.0 : c * gz / g2;
        const double ln = std::sqrt(flat ? sx * sx + sy * sy : sx * sx + sy * sy + sz * sz);
        const double f = std::min(1.0, 1.5 / std::max(ln, 1e-12));
        px += sx * f;
        py += sy * f;
        if (!flat) pz += sz * f;
        if (vx) {
            const double vn = flat ? *vx * nxv + *vy * nyv
                                   : *vx * nxv + *vy * nyv + *vz * nzv;
            *vx -= vn * nxv;
            *vy -= vn * nyv;
            if (!flat) *vz -= vn * nzv;
        }
    }
    return moved;
}

// pushOut beside a solid that moves at (svx, 0, 0): the velocity loses its
// normal component relative to the solid's. With plate given (a static solid
// folded in), the plate's speed applies only where the point sits inside the
// plate before the pass moves it; inside the static part the solid is at rest.
// The reference's _push_out with a piston, operation for operation at nz = 1.
template <class Arr, class Arr2>
inline bool pushOutMoving(const Arr& P, const Arr2* plate, double svx, int nx, int ny,
                          int nz, bool flat, double& px, double& py, double& pz,
                          double* vx, double* vy, double* vz) {
    if (!(phiValue(P, nx, ny, nz, flat, px, py, pz) < 0.0)) return false;
    bool moved = false;
    for (int pass = 0; pass < 3; ++pass) {
        double gx, gy, gz;
        const double phi = phiAt(P, nx, ny, nz, flat, px, py, pz, gx, gy, gz);
        if (!(phi < 0.0)) break;
        moved = true;
        const double g2 = std::max(flat ? gx * gx + gy * gy : gx * gx + gy * gy + gz * gz,
                                   1e-12);
        const double rg = std::sqrt(g2);
        const double nxv = gx / rg, nyv = gy / rg, nzv = gz / rg;
        const double c = -(phi - 0.01);
        const double sx = c * gx / g2, sy = c * gy / g2, sz = flat ? 0.0 : c * gz / g2;
        const double ln = std::sqrt(flat ? sx * sx + sy * sy : sx * sx + sy * sy + sz * sz);
        const double f = std::min(1.0, 1.5 / std::max(ln, 1e-12));
        const double vsx = (!plate || phiValue(*plate, nx, ny, nz, flat, px, py, pz) < 0.0)
                               ? svx : 0.0;
        px += sx * f;
        py += sy * f;
        if (!flat) pz += sz * f;
        if (vx) {
            const double vn = flat ? (*vx - vsx) * nxv + *vy * nyv
                                   : (*vx - vsx) * nxv + *vy * nyv + *vz * nzv;
            *vx -= vn * nxv;
            *vy -= vn * nyv;
            if (!flat) *vz -= vn * nzv;
        }
    }
    return moved;
}

} // namespace Solids
} // namespace PFFlip

#endif // PF_FLIP_SOLIDS_H
