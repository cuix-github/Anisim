#include "PFFlipSolve.h"
#include "PFFlipAdaptive.h"
#include "PFFlipScatter.h"
#include "PFFlipSolids.h"
#include "PFFlipSprayEvents.h"

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <random>
#include <cmath>
#include <unordered_map>
#include <vector>

// Milestone 6: the PF-FLIP step in 3D with a multigrid-preconditioned solve.
// Same structure as the 2D node (which mirrors reference2d/pfflip2d.py);
// at nz = 1 the w faces are all walls and this reduces exactly to the 2D
// method, which validate_3d.py uses as the lockstep bridge. Double precision
// internally. The multigrid is a geometric V-cycle on the variable-beta
// Poisson operator: damped-Jacobi smoothing, cell-average restriction,
// constant prolongation, face-averaged coarse coefficients, explicit 1/h^2
// per level; used inside flexible PCG so mild asymmetry cannot stall CG.

namespace {

using Vec = std::vector<double>;

struct Level {
    int nx, ny, nz;
    double inv_h2;
    Vec bu, bv, bw;        // face coefficients, walls zero
    Vec diagInv;           // inverse diagonal of A
    Vec x, b, r;           // work vectors (cells)
    size_t nc() const { return static_cast<size_t>(nx) * ny * nz; }
    size_t iu(int i, int j, int k) const { return (static_cast<size_t>(i) * ny + j) * nz + k; }
    size_t iv(int i, int j, int k) const { return (static_cast<size_t>(i) * (ny + 1) + j) * nz + k; }
    size_t iw(int i, int j, int k) const { return (static_cast<size_t>(i) * ny + j) * (nz + 1) + k; }
    size_t ic(int i, int j, int k) const { return (static_cast<size_t>(i) * ny + j) * nz + k; }
};

void applyA(const Level& L, const Vec& x, Vec& out) {
    const int NX = L.nx, NY = L.ny, NZ = L.nz;
#pragma omp parallel for schedule(static)
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NY; ++j)
            for (int k = 0; k < NZ; ++k) {
                const size_t c = L.ic(i, j, k);
                double acc = 0.0;
                if (i + 1 < NX) acc += L.bu[L.iu(i + 1, j, k)] * (x[L.ic(i + 1, j, k)] - x[c]);
                if (i > 0)      acc -= L.bu[L.iu(i, j, k)] * (x[c] - x[L.ic(i - 1, j, k)]);
                if (j + 1 < NY) acc += L.bv[L.iv(i, j + 1, k)] * (x[L.ic(i, j + 1, k)] - x[c]);
                if (j > 0)      acc -= L.bv[L.iv(i, j, k)] * (x[c] - x[L.ic(i, j - 1, k)]);
                if (k + 1 < NZ) acc += L.bw[L.iw(i, j, k + 1)] * (x[L.ic(i, j, k + 1)] - x[c]);
                if (k > 0)      acc -= L.bw[L.iw(i, j, k)] * (x[c] - x[L.ic(i, j, k - 1)]);
                out[c] = -acc * L.inv_h2;
            }
}

void buildDiag(Level& L) {
    L.diagInv.assign(L.nc(), 0.0);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L.nx; ++i)
        for (int j = 0; j < L.ny; ++j)
            for (int k = 0; k < L.nz; ++k) {
                const double dg = (L.bu[L.iu(i, j, k)] + L.bu[L.iu(i + 1, j, k)] +
                                   L.bv[L.iv(i, j, k)] + L.bv[L.iv(i, j + 1, k)] +
                                   L.bw[L.iw(i, j, k)] + L.bw[L.iw(i, j, k + 1)]) * L.inv_h2;
                L.diagInv[L.ic(i, j, k)] = dg > 0.0 ? 1.0 / dg : 0.0;
            }
}

void smooth(Level& L, Vec& x, const Vec& b, int sweeps, Vec& tmp) {
    const double omega = 0.8;
    for (int s = 0; s < sweeps; ++s) {
        applyA(L, x, tmp);
        const long long n = static_cast<long long>(x.size());
#pragma omp parallel for schedule(static)
        for (long long c = 0; c < n; ++c)
            x[c] += omega * L.diagInv[c] * (b[c] - tmp[c]);
    }
}

int half(int n) { return n > 1 ? (n + 1) / 2 : 1; }

// Average the <=4 fine faces a coarse face covers; clipped at odd boundaries.
double faceAvgU(const Level& f, int I, int J, int K, int cny, int cnz) {
    (void)cny; (void)cnz;
    const int i = std::min(2 * I, f.nx);
    double s = 0.0; int n = 0;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            const int j = 2 * J + a, k = 2 * K + b;
            if (j < f.ny && k < f.nz) { s += f.bu[f.iu(i, j, k)]; ++n; }
        }
    return n ? s / n : 0.0;
}
double faceAvgV(const Level& f, int I, int J, int K) {
    const int j = std::min(2 * J, f.ny);
    double s = 0.0; int n = 0;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            const int i = 2 * I + a, k = 2 * K + b;
            if (i < f.nx && k < f.nz) { s += f.bv[f.iv(i, j, k)]; ++n; }
        }
    return n ? s / n : 0.0;
}
double faceAvgW(const Level& f, int I, int J, int K) {
    const int k = std::min(2 * K, f.nz);
    double s = 0.0; int n = 0;
    for (int a = 0; a < 2; ++a)
        for (int b = 0; b < 2; ++b) {
            const int i = 2 * I + a, j = 2 * J + b;
            if (i < f.nx && j < f.ny) { s += f.bw[f.iw(i, j, k)]; ++n; }
        }
    return n ? s / n : 0.0;
}

std::vector<Level> buildHierarchy(Level&& fine) {
    std::vector<Level> ls;
    ls.push_back(std::move(fine));
    while (ls.size() < 12) {
        const Level& f = ls.back();
        if (f.nc() <= 64) break;
        Level c;
        c.nx = half(f.nx); c.ny = half(f.ny); c.nz = half(f.nz);
        c.inv_h2 = f.inv_h2 * 0.25;
        c.bu.assign(static_cast<size_t>(c.nx + 1) * c.ny * c.nz, 0.0);
        c.bv.assign(static_cast<size_t>(c.nx) * (c.ny + 1) * c.nz, 0.0);
        c.bw.assign(static_cast<size_t>(c.nx) * c.ny * (c.nz + 1), 0.0);
#pragma omp parallel for schedule(static)
        for (int I = 0; I <= c.nx; ++I)
            for (int J = 0; J < c.ny; ++J)
                for (int K = 0; K < c.nz; ++K)
                    c.bu[c.iu(I, J, K)] = (I == 0 || I == c.nx) ? 0.0 : faceAvgU(f, I, J, K, c.ny, c.nz);
#pragma omp parallel for schedule(static)
        for (int I = 0; I < c.nx; ++I)
            for (int J = 0; J <= c.ny; ++J)
                for (int K = 0; K < c.nz; ++K)
                    c.bv[c.iv(I, J, K)] = (J == 0 || J == c.ny) ? 0.0 : faceAvgV(f, I, J, K);
#pragma omp parallel for schedule(static)
        for (int I = 0; I < c.nx; ++I)
            for (int J = 0; J < c.ny; ++J)
                for (int K = 0; K <= c.nz; ++K)
                    c.bw[c.iw(I, J, K)] = (K == 0 || K == c.nz) ? 0.0 : faceAvgW(f, I, J, K);
        buildDiag(c);
        c.x.assign(c.nc(), 0.0); c.b.assign(c.nc(), 0.0); c.r.assign(c.nc(), 0.0);
        ls.push_back(std::move(c));
    }
    return ls;
}

void vcycle(std::vector<Level>& ls, size_t l, Vec& tmp) {
    Level& L = ls[l];
    if (l + 1 == ls.size()) {
        smooth(L, L.x, L.b, 40, tmp);
        return;
    }
    smooth(L, L.x, L.b, 2, tmp);
    applyA(L, L.x, tmp);
    const long long ncl = static_cast<long long>(L.nc());
#pragma omp parallel for schedule(static)
    for (long long c = 0; c < ncl; ++c) L.r[c] = L.b[c] - tmp[c];
    Level& C = ls[l + 1];
    std::fill(C.x.begin(), C.x.end(), 0.0);
    // Gather form of the cell-average restriction: each coarse cell sums its
    // fine children in fine-lexicographic order - the exact addition order
    // the old fine-loop scatter produced - so coarse cells become independent
    // writers and the result stays bit-identical to the serial version.
#pragma omp parallel for schedule(static)
    for (int I = 0; I < C.nx; ++I)
        for (int J = 0; J < C.ny; ++J)
            for (int K = 0; K < C.nz; ++K) {
                double s = 0.0;
                int n = 0;
                for (int a = 0; a < 2; ++a) {
                    const int i = 2 * I + a;
                    if (i >= L.nx) continue;
                    for (int b2 = 0; b2 < 2; ++b2) {
                        const int j = 2 * J + b2;
                        if (j >= L.ny) continue;
                        for (int d = 0; d < 2; ++d) {
                            const int k = 2 * K + d;
                            if (k >= L.nz) continue;
                            s += L.r[L.ic(i, j, k)];
                            ++n;
                        }
                    }
                }
                C.b[C.ic(I, J, K)] = n ? s / n : 0.0;
            }
    vcycle(ls, l + 1, tmp);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < L.nx; ++i)
        for (int j = 0; j < L.ny; ++j)
            for (int k = 0; k < L.nz; ++k)
                L.x[L.ic(i, j, k)] +=
                    C.x[C.ic(i / 2, std::min(j / 2, C.ny - 1), std::min(k / 2, C.nz - 1))];
    smooth(L, L.x, L.b, 2, tmp);
}

// Trilinear sampling of a staggered grid with per-axis index clipping,
// mirroring the 2D _sample_faces (and collapsing in any dimension of size 1).
double sampleGrid(const Vec& g, int gx, int gy, int gz,
                  double offx, double offy, double offz,
                  double px, double py, double pz) {
    const double qx = px - offx, qy = py - offy, qz = pz - offz;
    const double bx = std::floor(qx), by = std::floor(qy), bz = std::floor(qz);
    const double fx = qx - bx, fy = qy - by, fz = qz - bz;
    double out = 0.0;
    for (int di = 0; di <= 1; ++di)
        for (int dj = 0; dj <= 1; ++dj)
            for (int dk = 0; dk <= 1; ++dk) {
                const double w = (di ? fx : 1.0 - fx) * (dj ? fy : 1.0 - fy) *
                                 (dk ? fz : 1.0 - fz);
                const int ii = std::clamp(static_cast<int>(bx) + di, 0, gx - 1);
                const int jj = std::clamp(static_cast<int>(by) + dj, 0, gy - 1);
                const int kk = std::clamp(static_cast<int>(bz) + dk, 0, gz - 1);
                out += w * g[(static_cast<size_t>(ii) * gy + jj) * gz + kk];
            }
    return out;
}

// ---- ST-FLIP helpers: temporal kernel, its Gauss-Legendre normalization, smoothstep.
double wtKernel(double tau) {
    const double d = 1.0 - (tau - 0.5) * (tau - 0.5);
    return d > 0.0 ? (35.0 / 16.0) * d * d * d : 0.0;
}
// numpy.polynomial.legendre.leggauss(8), both scaled by 0.5 (as the 2D reference).
const double GL_N[8] = {-0.4801449282487681, -0.3983332387068134, -0.2627662049581645,
                        -0.0917173212478249, 0.0917173212478249, 0.2627662049581645,
                        0.3983332387068134, 0.4801449282487681};
const double GL_W[8] = {0.0506142681451881, 0.1111905172266872, 0.1568533229389436,
                        0.1813418916891810, 0.1813418916891810, 0.1568533229389436,
                        0.1111905172266872, 0.0506142681451881};
double wtNorm(double xi) {
    double z = 0.0;
    for (int i = 0; i < 8; ++i) z += GL_W[i] * wtKernel(xi * GL_N[i]);
    return std::max(z, 1e-9);
}

// Droplet radii (PF-FLIP Section 7): a standard normal from a hash of (step,
// index), the reference's Sim.droplet_normal up to the libm calls: splitmix64 at
// (step + 1) * golden + index and at the next state, the two outputs' top 53 bits
// as uniforms (the first on (0, 1], the second on [0, 1)), Box-Muller.
inline uint64_t splitmix64(uint64_t s) {
    uint64_t z = s + 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

inline double dropletNormal(int step_index, uint64_t index) {
    const uint64_t s = static_cast<uint64_t>(step_index + 1) * 0x9E3779B97F4A7C15ULL + index;
    const double u1 = (static_cast<double>(splitmix64(s) >> 11) + 1.0) * 0x1.0p-53;
    const double u2 = static_cast<double>(splitmix64(s + 1) >> 11) * 0x1.0p-53;
    return std::sqrt(-2.0 * std::log(u1)) * std::cos(6.283185307179586 * u2);
}

// the reference's Sim.droplet_uniform: a uniform on [0, 1) from the state two past
// the normal's first (splitmix64 at (step + 1) * golden + index + 2)
inline double dropletUniform(int step_index, uint64_t index) {
    const uint64_t s = static_cast<uint64_t>(step_index + 1) * 0x9E3779B97F4A7C15ULL + index + 2;
    return static_cast<double>(splitmix64(s) >> 11) * 0x1.0p-53;
}

} // namespace

namespace PFFlip {
namespace Solve {

void step_3d(int nx, int ny, int nz, float dt_in, float gravity, float rho_liquid,
             float rho_air, float alpha_liquid, float alpha_air, int max_iterations,
             float tolerance, int preconditioner, int pressure_solver,
             int adaptive_band, int adaptive_levels, int escape, float esc_phi,
             float drag_droplet, float drag_bubble, float buoyancy, float rho0_face,
             int st, int st_seed, int step_index,
             int adapt, float coarse_gain, float ws_epsilon,
             const Amino::Array<float>& scale_in,
             const Amino::Array<Bifrost::Math::float3>& obstacles_min,
             const Amino::Array<Bifrost::Math::float3>& obstacles_max,
             int solid_mode, const Amino::Array<float>& solid_phi,
             float piston_stroke, float piston_period, float piston_x0,
             float piston_width, float piston_top,
             int droplet_model, float droplet_r_median, float droplet_r_sigma,
             float droplet_evap_k, float droplet_r_min, float metres_per_cell,
             float air_nu, const Amino::Array<float>& radius_in,
             int emit_model, float emit_grad_threshold, float emit_rate, int emit_max_per_cell,
             int emit_cloud_size, float emit_event_radius, float emit_spread,
             float emit_velocity_spread, float emit_radius_max, int emit_max_active,
             int emit_seed, float spray_dt_max,
             float spray_cfl, const Amino::Array<float>& spray_in,
             float zone_length, float zone_depth, float zone_ramp,
             const Amino::Array<float>& waves,
             int external_forces,
             const Amino::Array<float>& liquid_acceleration_u,
             const Amino::Array<float>& liquid_acceleration_v,
             const Amino::Array<float>& liquid_acceleration_w,
             const Amino::Array<float>& linear_drag_u,
             const Amino::Array<float>& linear_drag_v,
             const Amino::Array<float>& linear_drag_w,
             const Amino::Array<float>& tau_in,
             const Amino::Array<float>& u_wt, const Amino::Array<float>& v_wt,
             const Amino::Array<float>& w_wt,
             const Amino::Array<float>& u_ws, const Amino::Array<float>& v_ws,
             const Amino::Array<float>& w_ws,
             const Amino::Array<float>& u2_mass, const Amino::Array<float>& u2_mom,
             const Amino::Array<float>& u2_phase, const Amino::Array<float>& u2_wt,
             const Amino::Array<float>& u2_ws,
             const Amino::Array<float>& v2_mass, const Amino::Array<float>& v2_mom,
             const Amino::Array<float>& v2_phase, const Amino::Array<float>& v2_wt,
             const Amino::Array<float>& v2_ws,
             const Amino::Array<float>& w2_mass, const Amino::Array<float>& w2_mom,
             const Amino::Array<float>& w2_phase, const Amino::Array<float>& w2_wt,
             const Amino::Array<float>& w2_ws,
             const Amino::Array<float>& u_mass, const Amino::Array<float>& u_mom,
             const Amino::Array<float>& u_phase, const Amino::Array<float>& v_mass,
             const Amino::Array<float>& v_mom, const Amino::Array<float>& v_phase,
             const Amino::Array<float>& w_mass, const Amino::Array<float>& w_mom,
             const Amino::Array<float>& w_phase,
             const Amino::Array<Bifrost::Math::float3>& positions,
             const Amino::Array<Bifrost::Math::float3>& velocities,
             const Amino::Array<float>& particle_phase,
             Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& out_positions,
             Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& out_velocities,
             Amino::Ptr<Amino::Array<float>>& out_mass,
             Amino::Ptr<Amino::Array<float>>& out_mom_x,
             Amino::Ptr<Amino::Array<float>>& out_mom_y,
             Amino::Ptr<Amino::Array<float>>& out_mom_z,
             Amino::Ptr<Amino::Array<float>>& out_u,
             Amino::Ptr<Amino::Array<float>>& out_v,
             Amino::Ptr<Amino::Array<float>>& out_w,
             Amino::Ptr<Amino::Array<float>>& pressure,
             Amino::Ptr<Amino::Array<float>>& out_escaped,
             Amino::Ptr<Amino::Array<float>>& tau_out,
             Amino::Ptr<Amino::Array<float>>& scale_out,
             Amino::Ptr<Amino::Array<float>>& radius_out,
             Amino::Ptr<Amino::Array<float>>& spray_out,
             Amino::Ptr<Amino::Array<bool>>& mask_fine,
             Amino::Ptr<Amino::Array<bool>>& mask_coarse,
             Amino::Ptr<Amino::Array<float>>& out_phase_state,
             Amino::Ptr<Amino::Array<float>>& out_wt,
             Amino::Ptr<Amino::Array<float>>& out_wtph,
             Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& out_pos_synced,
             int& iterations_used, float& final_residual, float& max_divergence_after,
             float& max_speed, const Amino::RuntimeServices& services) {
    const int NX = std::max(2, nx), NY = std::max(2, ny), NZ = std::max(1, nz);
    const double dt = dt_in;
    const size_t np = positions.size();
    auto refuseSpray = [&]() {
        using namespace Amino::MessageLiterals;
        services.logError("step_3d: invalid spray event configuration or intensity"_amt);
        std::fprintf(stderr, "step_3d: error: invalid spray event configuration or intensity\n");
        out_positions = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>>();
        out_velocities = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>>();
        out_pos_synced = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>>();
        for (auto* o : {&out_mass,&out_mom_x,&out_mom_y,&out_mom_z,&out_u,&out_v,&out_w,
                       &pressure,&out_escaped,&tau_out,&scale_out,&radius_out,&spray_out,
                       &out_phase_state,&out_wt,&out_wtph})
            *o = Amino::newClassPtr<Amino::Array<float>>();
        mask_fine = Amino::newClassPtr<Amino::Array<bool>>();
        mask_coarse = Amino::newClassPtr<Amino::Array<bool>>();
        iterations_used=-1; final_residual=max_divergence_after=max_speed=0;
    };
    if (emit_model==2 || droplet_model==2) {
        bool valid=escape!=0 && droplet_model!=0 && std::isfinite(dt) && dt>0 &&
            std::isfinite(droplet_r_sigma) && droplet_r_sigma>=0 && droplet_r_sigma<=2 &&
            std::isfinite(droplet_r_min) && droplet_r_min>0 &&
            std::isfinite(metres_per_cell) && metres_per_cell>0 &&
            std::isfinite(air_nu) && air_nu>0 &&
            std::isfinite(droplet_evap_k) && droplet_evap_k>=0;
        if (droplet_model==2) valid &= std::isfinite(spray_dt_max) && spray_dt_max>0 &&
                                      std::isfinite(spray_cfl) && spray_cfl>0;
        if (emit_model==2) valid &= emit_cloud_size>=1 && emit_cloud_size<=256 &&
            std::isfinite(emit_grad_threshold) && emit_grad_threshold>=0 &&
            std::isfinite(emit_rate) && emit_rate>=0 && emit_rate<=1000 &&
            std::isfinite(emit_event_radius) && emit_event_radius>0 &&
            std::isfinite(emit_radius_max) && emit_radius_max>=droplet_r_min &&
            std::isfinite(emit_spread) && emit_spread>=0 &&
            std::isfinite(emit_velocity_spread) && emit_velocity_spread>=0 &&
            emit_max_active>=0 && emit_max_active<=10000000 && emit_max_per_cell>=0 && emit_seed>=0;
        if (!valid) { refuseSpray(); return; }
    }

    // PFFLIP_PROFILE=1 prints per-phase wall times to stderr (spike 1); the
    // phase boundaries are the ---- section markers below.
    static const bool prof = []() {
        const char* e = std::getenv("PFFLIP_PROFILE");
        return e && e[0] == '1';
    }();
    using Clock = std::chrono::steady_clock;
    Clock::time_point tMark = Clock::now();
    const Clock::time_point tBegin = tMark;
    double msAsm = 0.0, msSolve = 0.0, msCorr = 0.0, msPart = 0.0;
    auto lap = [&](double& acc) {
        if (!prof) return;
        const Clock::time_point now = Clock::now();
        acc = std::chrono::duration<double, std::milli>(now - tMark).count();
        tMark = now;
    };

    Level L;
    L.nx = NX; L.ny = NY; L.nz = NZ; L.inv_h2 = 1.0;
    const size_t nu = static_cast<size_t>(NX + 1) * NY * NZ;
    const size_t nv = static_cast<size_t>(NX) * (NY + 1) * NZ;
    const size_t nw = static_cast<size_t>(NX) * NY * (NZ + 1);

    // Generic external fields, independent of the scene recipe. Accelerations
    // are cells/s^2, multiplied by the measured liquid fraction. Drag is 1/s,
    // acting on both phases. A scalar broadcasts; enabled empty/malformed or
    // nonfinite fields are errors (a missing cache must not silently stop forcing).
    if (external_forces) {
        auto validField = [](const Amino::Array<float>& a, size_t n, bool drag) {
            if (a.size() != 1 && a.size() != n) return false;
            for (size_t i = 0; i < a.size(); ++i)
                if (!std::isfinite(a[i]) || (drag && a[i] < 0.0f)) return false;
            return true;
        };
        if (!validField(liquid_acceleration_u, nu, false) ||
            !validField(liquid_acceleration_v, nv, false) ||
            !validField(liquid_acceleration_w, nw, false) ||
            !validField(linear_drag_u, nu, true) ||
            !validField(linear_drag_v, nv, true) ||
            !validField(linear_drag_w, nw, true)) {
            using namespace Amino::MessageLiterals;
            services.logError("step_3d: invalid external acceleration/drag field; nothing simulated"_amt);
            std::fprintf(stderr, "step_3d: error: invalid external acceleration/drag field; nothing simulated\n");
            out_positions = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>>();
            out_velocities = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>>();
            out_pos_synced = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>>();
            for (auto* o : {&out_mass, &out_mom_x, &out_mom_y, &out_mom_z, &out_u, &out_v,
                           &out_w, &pressure, &out_escaped, &tau_out, &scale_out, &radius_out,
                           &spray_out, &out_phase_state, &out_wt, &out_wtph})
                *o = Amino::newClassPtr<Amino::Array<float>>();
            mask_fine = Amino::newClassPtr<Amino::Array<bool>>();
            mask_coarse = Amino::newClassPtr<Amino::Array<bool>>();
            iterations_used = -1;
            final_residual = max_divergence_after = max_speed = 0.0f;
            return;
        }
    }

    // ---- channels -> face velocities (u = mom/mass where mass > 0).
    Vec u(nu, 0.0), v(nv, 0.0), w(nw, 0.0);
    // Two-tier combination: the stock splat returns weighted MEANS per tier; the
    // weight SUMS are recovered from the epsilon-inversion channel ws, where the
    // graph splats the constant 1 with add_to_denominator = eps, so
    // mean = S/(S+eps)  =>  S = eps*mean/(1-mean). Combined mean of q is then
    // (S_A q_A + S_B q_B)/(S_A + S_B). With adapt off (or no tier-B channels) the
    // tier-A channels pass straight through.
    const bool twoTier = adapt != 0 && u2_mass.size() > 0;
    const double epsw = ws_epsilon > 0.0f ? ws_epsilon : 16.0;
    auto wsum = [&](const Amino::Array<float>& ws, size_t k) {
        const double a = k < ws.size()
            ? std::clamp(static_cast<double>(ws[k]), 0.0, 0.999999) : 0.0;
        return epsw * a / (1.0 - a);
    };
    auto at = [](const Amino::Array<float>& arr, size_t k) {
        return k < arr.size() ? static_cast<double>(arr[k]) : 0.0;
    };
    Vec chM_u(nu), chO_u(nu), chP_u(nu), chW_u(nu);
    Vec chM_v(nv), chO_v(nv), chP_v(nv), chW_v(nv);
    Vec chM_w(nw), chO_w(nw), chP_w(nw), chW_w(nw);
    const bool hasWt = u_wt.size() > 0;
    auto combine = [&](size_t n, const Amino::Array<float>& mA, const Amino::Array<float>& oA,
                       const Amino::Array<float>& pA, const Amino::Array<float>& wtA,
                       const Amino::Array<float>& wsA, const Amino::Array<float>& mB,
                       const Amino::Array<float>& oB, const Amino::Array<float>& pB,
                       const Amino::Array<float>& wtB, const Amino::Array<float>& wsB,
                       Vec& M, Vec& O, Vec& P, Vec& W) {
        for (size_t k = 0; k < n; ++k) {
            if (!twoTier) {
                M[k] = at(mA, k); O[k] = at(oA, k); P[k] = at(pA, k);
                W[k] = hasWt ? at(wtA, k) : 1.0;
                continue;
            }
            const double sA = wsum(wsA, k), sB = wsum(wsB, k);
            const double den = sA + sB;
            if (den < 1e-12) { M[k] = O[k] = P[k] = 0.0; W[k] = hasWt ? 0.0 : 1.0; continue; }
            M[k] = (sA * at(mA, k) + sB * at(mB, k)) / den;
            O[k] = (sA * at(oA, k) + sB * at(oB, k)) / den;
            P[k] = (sA * at(pA, k) + sB * at(pB, k)) / den;
            W[k] = hasWt ? (sA * at(wtA, k) + sB * at(wtB, k)) / den : 1.0;
        }
    };
    combine(nu, u_mass, u_mom, u_phase, u_wt, u_ws,
            u2_mass, u2_mom, u2_phase, u2_wt, u2_ws, chM_u, chO_u, chP_u, chW_u);
    combine(nv, v_mass, v_mom, v_phase, v_wt, v_ws,
            v2_mass, v2_mom, v2_phase, v2_wt, v2_ws, chM_v, chO_v, chP_v, chW_v);
    combine(nw, w_mass, w_mom, w_phase, w_wt, w_ws,
            w2_mass, w2_mom, w2_phase, w2_wt, w2_ws, chM_w, chO_w, chP_w, chW_w);
    auto normalize = [](Vec& g, const Vec& mass, const Vec& mom) {
        for (size_t c = 0; c < g.size(); ++c)
            g[c] = mass[c] > 0.0 ? mom[c] / std::max(mass[c], 1e-30) : 0.0;
    };
    normalize(u, chM_u, chO_u);
    normalize(v, chM_v, chO_v);
    normalize(w, chM_w, chO_w);

    const size_t nObs = std::min(obstacles_min.size(), obstacles_max.size());
    // apply fn(face-kind, i, j, k) over every face whose centre lies inside a box.
    auto forSolidFaces = [&](auto&& fn) {
        for (size_t b = 0; b < nObs; ++b) {
            const auto& m0 = obstacles_min[b];
            const auto& m1 = obstacles_max[b];
            for (int i = std::max(0, (int)std::ceil(m0.x));
                 i <= std::min(NX, (int)std::floor(m1.x)); ++i)
                for (int j = std::max(0, (int)std::ceil(m0.y - 0.5));
                     j <= std::min(NY - 1, (int)std::floor(m1.y - 0.5)); ++j)
                    for (int k = std::max(0, (int)std::ceil(m0.z - 0.5));
                         k <= std::min(NZ - 1, (int)std::floor(m1.z - 0.5)); ++k)
                        fn(0, i, j, k);
            for (int i = std::max(0, (int)std::ceil(m0.x - 0.5));
                 i <= std::min(NX - 1, (int)std::floor(m1.x - 0.5)); ++i)
                for (int j = std::max(0, (int)std::ceil(m0.y));
                     j <= std::min(NY, (int)std::floor(m1.y)); ++j)
                    for (int k = std::max(0, (int)std::ceil(m0.z - 0.5));
                         k <= std::min(NZ - 1, (int)std::floor(m1.z - 0.5)); ++k)
                        fn(1, i, j, k);
            for (int i = std::max(0, (int)std::ceil(m0.x - 0.5));
                 i <= std::min(NX - 1, (int)std::floor(m1.x - 0.5)); ++i)
                for (int j = std::max(0, (int)std::ceil(m0.y - 0.5));
                     j <= std::min(NY - 1, (int)std::floor(m1.y - 0.5)); ++j)
                    for (int k = std::max(0, (int)std::ceil(m0.z));
                         k <= std::min(NZ, (int)std::floor(m1.z)); ++k)
                        fn(2, i, j, k);
        }
    };
    // ---- level-set solids: every face's open fraction, 1 open to 0 closed.
    // solid_mode 1 rounds it by the face centre (closed when the solid covers
    // the centre, the boxes' rule); 2 keeps the fraction (cut cells).
    const bool flatGrid = NZ == 1;
    const size_t nNode = static_cast<size_t>(NX + 1) * (NY + 1) * (NZ + 1);
    // the piston wavemaker (PFFlipSolve.h): a plate that moves, folded with the
    // static solid; it makes the solve cut-cell, as the reference's piston does
    const bool pistonOn = piston_period > 0.0f && piston_width > 0.0f;
    const bool staticOn = solid_mode != 0;
    const bool lsOn = staticOn || pistonOn;
    if (staticOn && solid_phi.size() != nNode) {
        // Refuse rather than simulate without the solid: a .npy of the wrong
        // dtype, for one, arrives here as an empty array. The error goes to
        // the graph's own channel (the runtime services) and to stderr, the
        // outputs stay empty and iterations_used reads -1.
        using namespace Amino::MessageLiterals;
        services.logError("step_3d: solid_phi holds {} values, not the {} grid nodes "
                          "of {}x{}x{}: one float32 per node, (nx+1)(ny+1)(nz+1); "
                          "nothing simulated"_amt,
                          static_cast<size_t>(solid_phi.size()), nNode, NX, NY, NZ);
        std::fprintf(stderr, "step_3d: error: solid_phi holds %zu values, not the %zu grid "
                             "nodes of %dx%dx%d: one float32 per node, (nx+1)(ny+1)(nz+1); "
                             "nothing simulated\n",
                     static_cast<size_t>(solid_phi.size()), nNode, NX, NY, NZ);
        std::fflush(stderr);
        out_positions = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>>();
        out_velocities = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>>();
        out_pos_synced = Amino::newClassPtr<Amino::Array<Bifrost::Math::float3>>();
        for (auto* o : {&out_mass, &out_mom_x, &out_mom_y, &out_mom_z, &out_u, &out_v,
                        &out_w, &pressure, &out_escaped, &tau_out, &scale_out,
                        &out_phase_state, &out_wt, &out_wtph, &radius_out})
            *o = Amino::newClassPtr<Amino::Array<float>>();
        mask_fine = Amino::newClassPtr<Amino::Array<bool>>();
        mask_coarse = Amino::newClassPtr<Amino::Array<bool>>();
        iterations_used = -1;
        final_residual = max_divergence_after = max_speed = 0.0f;
        return;
    }
    const bool lsCut = pistonOn || (staticOn && solid_mode == 2);
    // the plate at this step: its left face's x at t = step_index * dt, the
    // step's mean speed towards t + dt, its distance at the nodes, and the
    // distance of everything solid (the plate folded with the static solid)
    double pistonX = 0.0, pistonVx = 0.0;
    std::vector<double> phiPlate, phiAll;
    if (pistonOn) {
        const double S = piston_stroke, T = piston_period, x0 = piston_x0;
        const double wd = piston_width, top = piston_top;
        const double twoPi = 2.0 * 3.14159265358979323846;
        auto xLeft = [&](double tt) { return x0 + 0.5 * S * (1.0 - std::cos(twoPi * tt / T)); };
        const double t = static_cast<double>(step_index) * dt;
        pistonX = xLeft(t);
        pistonVx = dt > 0.0 ? (xLeft(t + dt) - pistonX) / dt : 0.0;
        phiPlate.resize(nNode);
        phiAll.resize(nNode);
        const size_t sj = static_cast<size_t>(NZ) + 1;
        const size_t si = (static_cast<size_t>(NY) + 1) * sj;
        for (int i = 0; i <= NX; ++i)
            for (int j = 0; j <= NY; ++j)
                for (int k = 0; k <= NZ; ++k) {
                    const size_t o = i * si + j * sj + k;
                    const double X = i, Y = j;
                    const double pp = std::max(std::max(pistonX - X, X - (pistonX + wd)), Y - top);
                    phiPlate[o] = pp;
                    phiAll[o] = staticOn ? std::min(pp, static_cast<double>(solid_phi[o])) : pp;
                }
    }
    Vec openU, openV, openW;
    Vec openSU, openSV, openSW;   // open against the static solid alone (piston only)
    bool lsTouches = false;
    if (lsOn) {
        const size_t sj = static_cast<size_t>(NZ) + 1;
        const size_t si = (static_cast<size_t>(NY) + 1) * sj;
        auto P = [&](int i, int j, int k) {
            return pistonOn ? phiAll[i * si + j * sj + k]
                            : static_cast<double>(solid_phi[i * si + j * sj + k]);
        };
        auto PS = [&](int i, int j, int k) {
            return static_cast<double>(solid_phi[i * si + j * sj + k]);
        };
        // one face from its corners (a, b along its first tangent, c, d the
        // same one step along its second): a flat grid takes the 2D
        // reference's rule on the k = 0 ends, a 3D grid the four-corner area
        auto open = [&](double a, double b, double c, double d) {
            if (!lsCut)
                return flatGrid ? (0.5 * (a + b) > 0 ? 1.0 : 0.0)
                                : (0.25 * ((a + b) + (c + d)) > 0 ? 1.0 : 0.0);
            const double f = flatGrid ? Solids::edgeOpen(a, b) : Solids::faceOpen(a, b, c, d);
            return f < Solids::kMinOpen ? 0.0 : f;
        };
        openU.assign(nu, 1.0); openV.assign(nv, 1.0); openW.assign(nw, 1.0);
#pragma omp parallel for schedule(static)
        for (int i = 0; i <= NX; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k < NZ; ++k)
                    openU[L.iu(i, j, k)] = open(P(i, j, k), P(i, j + 1, k),
                                                P(i, j, k + 1), P(i, j + 1, k + 1));
#pragma omp parallel for schedule(static)
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j <= NY; ++j)
                for (int k = 0; k < NZ; ++k)
                    openV[L.iv(i, j, k)] = open(P(i, j, k), P(i + 1, j, k),
                                                P(i, j, k + 1), P(i + 1, j, k + 1));
#pragma omp parallel for schedule(static)
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k <= NZ; ++k)
                    openW[L.iw(i, j, k)] = open(P(i, j, k), P(i + 1, j, k),
                                                P(i, j + 1, k), P(i + 1, j + 1, k));
        if (pistonOn) {
            openSU.assign(nu, 1.0); openSV.assign(nv, 1.0); openSW.assign(nw, 1.0);
            if (staticOn) {
#pragma omp parallel for schedule(static)
                for (int i = 0; i <= NX; ++i)
                    for (int j = 0; j < NY; ++j)
                        for (int k = 0; k < NZ; ++k)
                            openSU[L.iu(i, j, k)] = open(PS(i, j, k), PS(i, j + 1, k),
                                                         PS(i, j, k + 1), PS(i, j + 1, k + 1));
#pragma omp parallel for schedule(static)
                for (int i = 0; i < NX; ++i)
                    for (int j = 0; j <= NY; ++j)
                        for (int k = 0; k < NZ; ++k)
                            openSV[L.iv(i, j, k)] = open(PS(i, j, k), PS(i + 1, j, k),
                                                         PS(i, j, k + 1), PS(i + 1, j, k + 1));
#pragma omp parallel for schedule(static)
                for (int i = 0; i < NX; ++i)
                    for (int j = 0; j < NY; ++j)
                        for (int k = 0; k <= NZ; ++k)
                            openSW[L.iw(i, j, k)] = open(PS(i, j, k), PS(i + 1, j, k),
                                                         PS(i, j + 1, k), PS(i + 1, j + 1, k));
            }
        }
        for (size_t c = 0; c < nu && !lsTouches; ++c) lsTouches = openU[c] < 1.0;
        for (size_t c = 0; c < nv && !lsTouches; ++c) lsTouches = openV[c] < 1.0;
        for (size_t c = 0; c < nw && !lsTouches && !flatGrid; ++c) lsTouches = openW[c] < 1.0;
        // PFFLIP_SOLID_DUMP=<path> writes the three fractions (u, v, w, as
        // doubles) for gates/validate_level_set_node.py
        if (const char* path = std::getenv("PFFLIP_SOLID_DUMP")) {
            if (FILE* fp = std::fopen(path, "wb")) {
                std::fwrite(openU.data(), sizeof(double), nu, fp);
                std::fwrite(openV.data(), sizeof(double), nv, fp);
                std::fwrite(openW.data(), sizeof(double), nw, fp);
                std::fclose(fp);
            }
        }
    }
    // vel: velocities (the plate's faces take its speed), not the coefficients
    // ---- the wave zone (PFFlipSolve.h): after the projection the face velocities
    // within zone_length of the x = 0 wall are blended toward the sum of linear wave
    // trains at the end of the step, the weight the zone's exp weight times the face's
    // liquid fraction, and the liquid particles get the same blend after
    // their FLIP update; the wall faces the zone touches (x = 0, and the side walls
    // along it) are the trains' inlet, carrying their velocity in `walls` below
    // instead of 0, so the flux the trains need is not forced into the wall's cells.
    // The reference's _zone_faces / _zone_particles / _enforce_walls.
    struct ZoneTrain { double a, k, w, dx, dz, ph; };
    std::vector<ZoneTrain> zoneTrains;
    const bool zoneOn = zone_length > 0.0f && waves.size() >= 4 && zone_depth > 0.0f;
    const double zoneLen = zone_length, zoneH = zone_depth;
    double zoneT = 0.0;
    if (zoneOn) {
        const double g = std::fabs(static_cast<double>(gravity));
        zoneT = (static_cast<double>(step_index) + 1.0) * dt;
        double ramp = 1.0;
        if (zone_ramp > 0.0f) {
            const double s = std::min(1.0, zoneT / static_cast<double>(zone_ramp));
            ramp = s * s * (3.0 - 2.0 * s);
        }
        const double deg = 3.14159265358979323846 / 180.0;
        for (size_t q = 0; q + 3 < waves.size(); q += 4) {
            const double H = waves[q], T = waves[q + 1];
            const double th = static_cast<double>(waves[q + 2]) * deg, ph = waves[q + 3];
            if (!(T > 0.0) || !(H > 0.0) || !(g > 0.0)) continue;
            const double om = 2.0 * 3.14159265358979323846 / T;
            double k = om * om / g;
            for (int it2 = 0; it2 < 50; ++it2) {
                const double th2 = std::tanh(k * zoneH), ch = std::cosh(k * zoneH);
                const double f = g * k * th2 - om * om;
                const double df = g * (th2 + k * zoneH / (ch * ch));
                k -= f / df;
            }
            zoneTrains.push_back({0.5 * H * ramp, k, om, std::cos(th),
                                  flatGrid ? 0.0 : std::sin(th), static_cast<double>(ph)});
        }
    }
    // the blend weight at a point and the trains' velocity there
    auto zoneAt = [&](double x, double y, double z, double& ut, double& vt, double& wt) {
        ut = vt = wt = 0.0;
        if (!zoneOn || x >= zoneLen) return 0.0;
        const double s = x / zoneLen;
        double wz = 1.0 - (std::exp(std::pow(s, 3.5)) - 1.0) / 1.718281828459045;
        if (wz <= 0.0) return 0.0;
        const double ye = std::min(y, zoneH);
        for (const ZoneTrain& tr : zoneTrains) {
            const double psi = tr.k * (x * tr.dx + z * tr.dz) - tr.w * zoneT + tr.ph;
            const double sh = std::sinh(tr.k * zoneH);
            const double uh = tr.a * tr.w * std::cosh(tr.k * ye) / sh * std::cos(psi);
            ut += uh * tr.dx;
            wt += uh * tr.dz;
            vt += tr.a * tr.w * std::sinh(tr.k * ye) / sh * std::sin(psi);
        }
        return wz;
    };
    // the zone forces the water where the water is, not where the trains put their
    // surface: a face's weight is times its liquid fraction, the air is left alone.
    // Forcing the water below the trains' surface only (the first zone's fade) pumped
    // it one way: lifted on the crests, the water above the trains' troughs was left
    // where it was, and the zone's level rose a metre in 8 s under the mixed sea, the
    // unforced layer on top sheared to 10 m/s and the air over it to 18 (2026-10-18)
    auto liquidOf = [&](const Vec& P, const Vec& Wt, size_t c) {
        double p = P[c];
        if (hasWt) p = Wt[c] > 1e-12 ? p / Wt[c] : 0.0;
        return std::clamp(p, 0.0, 1.0);
    };
    auto walls = [&](Vec& uu, Vec& vv, Vec& ww, bool vel) {
        for (int j = 0; j < NY; ++j)
            for (int k = 0; k < NZ; ++k) { uu[L.iu(0, j, k)] = 0.0; uu[L.iu(NX, j, k)] = 0.0; }
        for (int i = 0; i < NX; ++i)
            for (int k = 0; k < NZ; ++k) { vv[L.iv(i, 0, k)] = 0.0; vv[L.iv(i, NY, k)] = 0.0; }
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NY; ++j) { ww[L.iw(i, j, 0)] = 0.0; ww[L.iw(i, j, NZ)] = 0.0; }
        forSolidFaces([&](int kind, int i, int j, int k) {
            if (kind == 0) uu[L.iu(i, j, k)] = 0.0;
            else if (kind == 1) vv[L.iv(i, j, k)] = 0.0;
            else ww[L.iw(i, j, k)] = 0.0;
        });
        if (lsOn) {
            for (size_t c = 0; c < nu; ++c) if (openU[c] <= 0.0) uu[c] = 0.0;
            for (size_t c = 0; c < nv; ++c) if (openV[c] <= 0.0) vv[c] = 0.0;
            for (size_t c = 0; c < nw; ++c) if (openW[c] <= 0.0) ww[c] = 0.0;
            // the faces the plate closes move with it (y and z: its speed is 0)
            if (vel && pistonOn)
                for (size_t c = 0; c < nu; ++c)
                    if (openU[c] <= 0.0 && openSU[c] > 0.0) uu[c] = pistonVx;
        }
        if (vel && zoneOn) {
            // the wave inlet: the x = 0 faces carry the trains' velocity (times the
            // face's liquid fraction), and so do the side-wall faces along the zone,
            // less the flux-weighted mean over every inlet face: a sealed box of
            // incompressible fluid has no pressure solution under a net inflow, so the
            // inlet keeps the box's volume and the wall's own column breathes by the mean
            struct Inlet { size_t c; double w, t, sign; };
            std::vector<Inlet> faces;
            double Q = 0.0, W = 0.0;
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k < NZ; ++k) {
                    double ut, vt, wt;
                    const double wgt = zoneAt(0.0, j + 0.5, k + 0.5, ut, vt, wt) * liquidOf(chP_u, chW_u, L.iu(0, j, k));
                    faces.push_back({L.iu(0, j, k), wgt, ut, 1.0});
                    Q += wgt * ut; W += wgt;
                }
            if (!flatGrid) {
                const int iMax = std::min(NX, static_cast<int>(std::ceil(zoneLen)) + 1);
                for (int i = 0; i < iMax; ++i)
                    for (int j = 0; j < NY; ++j) {
                        double ut, vt, wt;
                        double wgt = zoneAt(i + 0.5, j + 0.5, 0.0, ut, vt, wt) * liquidOf(chP_w, chW_w, L.iw(i, j, 0));
                        faces.push_back({L.iw(i, j, 0), wgt, wt, 1.0});
                        Q += wgt * wt; W += wgt;
                        wgt = zoneAt(i + 0.5, j + 0.5, static_cast<double>(NZ), ut, vt, wt) * liquidOf(chP_w, chW_w, L.iw(i, j, NZ));
                        faces.push_back({L.iw(i, j, NZ), wgt, wt, -1.0});
                        Q -= wgt * wt; W += wgt;
                    }
            }
            const double cmean = W > 0.0 ? Q / W : 0.0;
            // the x = 0 faces come first, then the side faces in pairs
            size_t q = 0;
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k < NZ; ++k, ++q)
                    uu[faces[q].c] = faces[q].w * (faces[q].t - cmean);
            for (; q < faces.size(); ++q)
                ww[faces[q].c] = faces[q].w * (faces[q].t - faces[q].sign * cmean);
        }
    };
    walls(u, v, w, true);
    Vec uStar = u, vStar = v, wStar = w;

    // ---- gravity on y faces, then walls.
    for (size_t c = 0; c < nv; ++c) v[c] += gravity * dt;
    if (external_forces) {
        auto force = [&](Vec& vel, const Vec& phase, const Vec& wt,
                         const Amino::Array<float>& a, const Amino::Array<float>& drag) {
            const long long n = static_cast<long long>(vel.size());
#pragma omp parallel for schedule(static)
            for (long long c = 0; c < n; ++c) {
                const double acc = a[a.size() == 1 ? 0 : c];
                const double rate = drag[drag.size() == 1 ? 0 : c];
                vel[c] = (vel[c] + dt * liquidOf(phase, wt, c) * acc) / (1.0 + dt * rate);
            }
        };
        force(u, chP_u, chW_u, liquid_acceleration_u, linear_drag_u);
        force(v, chP_v, chW_v, liquid_acceleration_v, linear_drag_v);
        force(w, chP_w, chW_w, liquid_acceleration_w, linear_drag_w);
    }
    walls(u, v, w, true);

    // ---- face coefficients beta = 1/rho(phi), walls zero.
    L.bu.assign(nu, 0.0); L.bv.assign(nv, 0.0); L.bw.assign(nw, 0.0);
    const double rl = rho_liquid, rg = rho_air;
    // With ST graphs the phase channel arrives premultiplied by the temporal weight
    // and must be divided by the splatted weight channel; an empty wt channel means
    // the plain scheme (wt = 1), keeping every existing graph and driver valid.
    auto beta = [&](Vec& bb, const Vec& ph, const Vec& wt) {
        for (size_t c = 0; c < bb.size(); ++c) {
            double p = ph[c];
            if (hasWt) p = wt[c] > 1e-12 ? p / wt[c] : 0.0;
            bb[c] = 1.0 / (rg + (rl - rg) * std::clamp(p, 0.0, 1.0));
        }
    };
    beta(L.bu, chP_u, chW_u); beta(L.bv, chP_v, chW_v); beta(L.bw, chP_w, chW_w);
    walls(L.bu, L.bv, L.bw, false);
    // cut cells: flux passes through a face's open part only, so the
    // operator's coefficients and the divergence take the open fraction; the
    // velocity update keeps the unweighted beta (the solid is static)
    Vec betaU, betaV, betaW;
    if (lsCut) {
        betaU = L.bu; betaV = L.bv; betaW = L.bw;
        for (size_t c = 0; c < nu; ++c) L.bu[c] *= openU[c];
        for (size_t c = 0; c < nv; ++c) L.bv[c] *= openV[c];
        for (size_t c = 0; c < nw; ++c) L.bw[c] *= openW[c];
    }
    buildDiag(L);
    L.x.assign(L.nc(), 0.0); L.b.assign(L.nc(), 0.0); L.r.assign(L.nc(), 0.0);

    // ---- right-hand side b = -div(u)/dt (cut cells: of the open-weighted u).
    Vec rhs(L.nc());
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NY; ++j)
            for (int k = 0; k < NZ; ++k)
                rhs[L.ic(i, j, k)] = pistonOn
                    // the part of a face the plate covers moves at its speed
                    // (open against the static solid, closed against all)
                    ? -(((openU[L.iu(i + 1, j, k)] * u[L.iu(i + 1, j, k)] +
                          (openSU[L.iu(i + 1, j, k)] - openU[L.iu(i + 1, j, k)]) * pistonVx) -
                         (openU[L.iu(i, j, k)] * u[L.iu(i, j, k)] +
                          (openSU[L.iu(i, j, k)] - openU[L.iu(i, j, k)]) * pistonVx)) +
                        (openV[L.iv(i, j + 1, k)] * v[L.iv(i, j + 1, k)] -
                         openV[L.iv(i, j, k)] * v[L.iv(i, j, k)]) +
                        (openW[L.iw(i, j, k + 1)] * w[L.iw(i, j, k + 1)] -
                         openW[L.iw(i, j, k)] * w[L.iw(i, j, k)])) / dt
                    : lsCut
                    ? -((openU[L.iu(i + 1, j, k)] * u[L.iu(i + 1, j, k)] -
                         openU[L.iu(i, j, k)] * u[L.iu(i, j, k)]) +
                        (openV[L.iv(i, j + 1, k)] * v[L.iv(i, j + 1, k)] -
                         openV[L.iv(i, j, k)] * v[L.iv(i, j, k)]) +
                        (openW[L.iw(i, j, k + 1)] * w[L.iw(i, j, k + 1)] -
                         openW[L.iw(i, j, k)] * w[L.iw(i, j, k)])) / dt
                    : -((u[L.iu(i + 1, j, k)] - u[L.iu(i, j, k)]) +
                        (v[L.iv(i, j + 1, k)] - v[L.iv(i, j, k)]) +
                        (w[L.iw(i, j, k + 1)] - w[L.iw(i, j, k)])) / dt;

    lap(msAsm);
    const size_t ncell = L.nc();
    Vec p;
    int it = 0;
    double finalResidual = 0.0;
    // solids are outside the adaptive solve's validated scope so far: the
    // uniform solve stays in charge when any are present. Only boxes reaching
    // the domain count: one wholly outside it (the loop drivers' far
    // placeholder) marks no face and holds no particle. A level-set solid
    // counts when it closes or cuts any face.
    size_t obsInDomain = lsTouches ? 1 : 0;
    for (size_t b = 0; b < nObs; ++b) {
        const auto& m0 = obstacles_min[b];
        const auto& m1 = obstacles_max[b];
        if (m1.x >= 0.0f && m0.x <= NX && m1.y >= 0.0f && m0.y <= NY && m1.z >= 0.0f && m0.z <= NZ)
            ++obsInDomain;
    }
    const bool adaptiveAsked = pressure_solver == 1 || pressure_solver == 2;
    const bool adaptiveSolve = adaptiveAsked && obsInDomain == 0;
    if (adaptiveAsked && obsInDomain > 0) {
        static bool warned = false;
        if (!warned) {
            std::fprintf(stderr, "step_3d: the adaptive solve does not handle solids "
                                 "yet; using the uniform solve\n");
            warned = true;
        }
    }
    if (adaptiveSolve) {
        // ---- the adaptive solve (milestone 6.5a): its refinement input is the
        // liquid particle count per cell over the seeding count (4 per cell in
        // 2D, 8 in 3D), as the reference's band_target.
        Vec liquid(ncell, 0.0);
        for (size_t pi = 0; pi < np; ++pi) {
            if (!(pi < particle_phase.size() && particle_phase[pi] > 0.5f)) continue;
            const int ci = std::clamp(static_cast<int>(std::floor(positions[pi].x)), 0, NX - 1);
            const int cj = std::clamp(static_cast<int>(std::floor(positions[pi].y)), 0, NY - 1);
            const int ck = std::clamp(static_cast<int>(std::floor(positions[pi].z)), 0, NZ - 1);
            liquid[L.ic(ci, cj, ck)] += 1.0;
        }
        const double perCell = NZ == 1 ? 4.0 : 8.0;
        for (double& c : liquid) c /= perCell;
        Adaptive::Settings as;
        as.band = adaptive_band > 0 ? adaptive_band : 5;
        as.levels = adaptive_levels > 0 ? adaptive_levels : 2;
        as.tolerance = tolerance;
        as.maxIterations = max_iterations;
        const Adaptive::Junction junction =
            pressure_solver == 2 ? Adaptive::Junction::Cleanup : Adaptive::Junction::Naive;
        Adaptive::Result ar = Adaptive::project(NX, NY, NZ, dt, as, junction, liquid, u, v, w,
                                                L.bu, L.bv, L.bw);
        p = std::move(ar.pressure);
        // both solves' iterations: the cost the port reports
        it = ar.iterations + ar.cleanupIterations;
        finalResidual = ar.residual;
        lap(msSolve);
    } else {
        // ---- flexible PCG with Jacobi (0) or multigrid V-cycle (1) preconditioner.
        std::vector<Level> ls;
        Vec tmp(L.nc());
        if (preconditioner == 1) ls = buildHierarchy(std::move(L));
        Level& F = preconditioner == 1 ? ls[0] : L;

        p.assign(ncell, 0.0);
        Vec r = rhs, z(ncell), d(ncell), q(ncell), rPrev(ncell, 0.0);
        auto precond = [&](const Vec& rr, Vec& zz) {
            if (preconditioner == 1) {
                ls[0].b = rr;
                std::fill(ls[0].x.begin(), ls[0].x.end(), 0.0);
                vcycle(ls, 0, tmp);
                zz = ls[0].x;
            } else {
                const long long nl = static_cast<long long>(ncell);
#pragma omp parallel for schedule(static)
                for (long long c = 0; c < nl; ++c) zz[c] = F.diagInv[c] * rr[c];
            }
        };
        precond(r, z);
        d = z;
        double rz = 0.0, rr2 = 0.0;
        for (size_t c = 0; c < ncell; ++c) { rz += r[c] * z[c]; rr2 += r[c] * r[c]; }
        double r0 = std::sqrt(rr2);
        if (r0 == 0.0) r0 = 1.0;
        for (it = 0; it < max_iterations; ++it) {
            rr2 = 0.0;
            for (size_t c = 0; c < ncell; ++c) rr2 += r[c] * r[c];
            if (std::sqrt(rr2) / r0 < tolerance) break;
            applyA(F, d, q);
            double dq = 0.0;
            for (size_t c = 0; c < ncell; ++c) dq += d[c] * q[c];
            if (dq == 0.0) break;
            const double a = rz / dq;
            rPrev = r;
            const long long ncl = static_cast<long long>(ncell);
#pragma omp parallel for schedule(static)
            for (long long c = 0; c < ncl; ++c) { p[c] += a * d[c]; r[c] -= a * q[c]; }
            precond(r, z);
            double rzNew = 0.0;   // flexible (Polak-Ribiere) update tolerates the V-cycle
            for (size_t c = 0; c < ncell; ++c) rzNew += (r[c] - rPrev[c]) * z[c];
            const double betaCG = std::max(0.0, rzNew / rz);
#pragma omp parallel for schedule(static)
            for (long long c = 0; c < ncl; ++c) d[c] = z[c] + betaCG * d[c];
            rz = 0.0;
            for (size_t c = 0; c < ncell; ++c) rz += r[c] * z[c];
            if (rz == 0.0) break;
        }

        finalResidual = std::sqrt(rr2) / r0;
        // Opt-in, single-step pressure audit. Preserve the actual doubles before
        // float32 output quantization: int32 dimensions, then bu/bv/bw, rhs, p,
        // and pre-projection u/v/w. The independent reader reassembles -D beta G
        // and checks b-Ap; this does not change the solve or its stopping rule.
        // Caller must provide a distinct path per invocation (not a loop path).
        if (const char* path = std::getenv("PFFLIP_PRESSURE_DUMP")) {
            if (FILE* fp = std::fopen(path, "wb")) {
                const int32_t dims[3] = {NX, NY, NZ};
                std::fwrite(dims, sizeof(int32_t), 3, fp);
                for (const Vec* a : {&F.bu, &F.bv, &F.bw, &rhs, &p, &u, &v, &w})
                    std::fwrite(a->data(), sizeof(double), a->size(), fp);
                std::fclose(fp);
            } else {
                std::fprintf(stderr, "step_3d: pressure audit dump failed: %s\n", path);
            }
        }
        lap(msSolve);
        // ---- velocity correction on interior faces, then walls.
        const Vec& cbu = lsCut ? betaU : F.bu;
        const Vec& cbv = lsCut ? betaV : F.bv;
        const Vec& cbw = lsCut ? betaW : F.bw;
        for (int i = 1; i < NX; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k < NZ; ++k)
                    u[F.iu(i, j, k)] -= dt * cbu[F.iu(i, j, k)] *
                                        (p[F.ic(i, j, k)] - p[F.ic(i - 1, j, k)]);
        for (int i = 0; i < NX; ++i)
            for (int j = 1; j < NY; ++j)
                for (int k = 0; k < NZ; ++k)
                    v[F.iv(i, j, k)] -= dt * cbv[F.iv(i, j, k)] *
                                        (p[F.ic(i, j, k)] - p[F.ic(i, j - 1, k)]);
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 1; k < NZ; ++k)
                    w[F.iw(i, j, k)] -= dt * cbw[F.iw(i, j, k)] *
                                        (p[F.ic(i, j, k)] - p[F.ic(i, j, k - 1)]);
    }
    walls(u, v, w, true);

    if (zoneOn) {
        const int iMax = std::min(NX, static_cast<int>(std::ceil(zoneLen)) + 1);
#pragma omp parallel for schedule(static)
        for (int i = 0; i <= iMax; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k < NZ; ++k) {
                    double ut, vt, wt;
                    const double wgt = zoneAt(i, j + 0.5, k + 0.5, ut, vt, wt) * liquidOf(chP_u, chW_u, L.iu(i, j, k));
                    if (wgt > 0.0) {
                        double& c = u[L.iu(i, j, k)];
                        c = (1.0 - wgt) * c + wgt * ut;
                    }
                }
#pragma omp parallel for schedule(static)
        for (int i = 0; i < iMax; ++i)
            for (int j = 0; j <= NY; ++j)
                for (int k = 0; k < NZ; ++k) {
                    double ut, vt, wt;
                    const double wgt = zoneAt(i + 0.5, j, k + 0.5, ut, vt, wt) * liquidOf(chP_v, chW_v, L.iv(i, j, k));
                    if (wgt > 0.0) {
                        double& c = v[L.iv(i, j, k)];
                        c = (1.0 - wgt) * c + wgt * vt;
                    }
                }
        if (!flatGrid) {
#pragma omp parallel for schedule(static)
            for (int i = 0; i < iMax; ++i)
                for (int j = 0; j < NY; ++j)
                    for (int k = 0; k <= NZ; ++k) {
                        double ut, vt, wt;
                        const double wgt = zoneAt(i + 0.5, j + 0.5, k, ut, vt, wt) * liquidOf(chP_w, chW_w, L.iw(i, j, k));
                        if (wgt > 0.0) {
                            double& c = w[L.iw(i, j, k)];
                            c = (1.0 - wgt) * c + wgt * wt;
                        }
                    }
        }
        walls(u, v, w, true);
    }

    // the divergence the solve removed (cut cells: of the open-weighted u)
    auto fu = [&](size_t c) {
        return pistonOn ? openU[c] * u[c] + (openSU[c] - openU[c]) * pistonVx
                        : (lsCut ? openU[c] * u[c] : u[c]);
    };
    auto fv = [&](size_t c) { return lsCut ? openV[c] * v[c] : v[c]; };
    auto fw = [&](size_t c) { return lsCut ? openW[c] * w[c] : w[c]; };
    double maxDiv = 0.0;
    for (int i = 0; i < NX; ++i)
        for (int j = 0; j < NY; ++j)
            for (int k = 0; k < NZ; ++k)
                maxDiv = std::max(maxDiv, std::fabs(
                    (fu(L.iu(i + 1, j, k)) - fu(L.iu(i, j, k))) +
                    (fv(L.iv(i, j + 1, k)) - fv(L.iv(i, j, k))) +
                    (fw(L.iw(i, j, k + 1)) - fw(L.iw(i, j, k)))));

    lap(msCorr);
    // ---- particle machinery. Sizes: with adapt on, particles carry a scale
    // (1 fine, 2 coarse); masses scale as scale^d (d = 2 flat, 3 otherwise) with a
    // calibrated coarse gain; the FLIP blend uses (1-alpha)/scale^2; deep-air fine
    // air merges 4-to-1 (8-to-1 in 3D) and coarse splits near the interface or
    // walls, mirroring the 2D reference's rules.
    const bool flat = NZ == 1;
    // a particle inside a solid back out: the static field as before, or with the
    // piston the folded field and the plate's speed
    const std::vector<double>* plateOnly = staticOn ? &phiPlate : nullptr;
    auto pushP = [&](double& x, double& y, double& z, double* vx, double* vy, double* vz) {
        if (pistonOn)
            return Solids::pushOutMoving(phiAll, plateOnly, pistonVx, NX, NY, NZ, flat, x, y, z,
                                         vx, vy, vz);
        return Solids::pushOut(solid_phi, NX, NY, NZ, flat, x, y, z, vx, vy, vz);
    };
    const bool stOn = st != 0;
    const bool adOn = adapt != 0;
    const double dExp = flat ? 2.0 : 3.0;
    auto scaleOf = [&](size_t pi) {
        return adOn && pi < scale_in.size() ? static_cast<double>(scale_in[pi]) : 1.0;
    };

    // count splats (Eq. 6 kernel) at given positions; value = scale^d per particle.
    auto countSplat = [&](const std::vector<double>& qx, const std::vector<double>& qy,
                          const std::vector<double>& qz, const std::vector<char>& isLiq,
                          const std::vector<char>& skipP,
                          const std::vector<double>& scl, Vec cntU[2], Vec cntV[2],
                          Vec cntW[2]) {
        const size_t n = qx.size();
        const Detail::OrderedScatterSlabs slabs(n, NX,
            [&](size_t i) { return qx[i]; }, [](size_t) { return 1.0; });
#pragma omp parallel for schedule(dynamic, 1) if(slabs.count() > 1)
        for (int slab = 0; slab < slabs.count(); ++slab) {
            for (size_t slot = slabs.begin(slab); slot < slabs.end(slab); ++slot) {
                const size_t pi = slabs.particle(slot);
                if (skipP[pi]) continue;           // emitted spray counts in no phase
                const int ph = isLiq[pi] ? 1 : 0;
                const double val = std::pow(scl[pi], dExp);
                const int ngrids = NZ == 1 ? 2 : 3;
                for (int grid = 0; grid < ngrids; ++grid) {
                    const double ox = grid == 0 ? 0.0 : 0.5;
                    const double oy = grid == 1 ? 0.0 : 0.5;
                    const double oz = grid == 2 ? 0.0 : 0.5;
                    const int gx = grid == 0 ? NX + 1 : NX;
                    const int gy = grid == 1 ? NY + 1 : NY;
                    const int gzd = grid == 2 ? NZ + 1 : NZ;
                    Vec* cnt = grid == 0 ? cntU : (grid == 1 ? cntV : cntW);
                    const double ax = qx[pi] - ox, ay = qy[pi] - oy, az = qz[pi] - oz;
                    const int bx = static_cast<int>(std::floor(ax + 0.5));
                    const int by = static_cast<int>(std::floor(ay + 0.5));
                    const int bz = static_cast<int>(std::floor(az + 0.5));
                    const double fx = ax - bx, fy = ay - by, fz = az - bz;
                    for (int di = -1; di <= 1; ++di)
                        for (int dj = -1; dj <= 1; ++dj)
                            for (int dk = -1; dk <= 1; ++dk) {
                                const int i = bx + di, j = by + dj, k = bz + dk;
                                if (i < slabs.lo(slab) || i >= slabs.hi(slab) || i >= gx ||
                                    j < 0 || j >= gy || k < 0 || k >= gzd)
                                    continue;
                                const double d2 = (fx - di) * (fx - di) + (fy - dj) * (fy - dj) +
                                                  (fz - dk) * (fz - dk);
                                const double wgt = std::max(1.0 - d2, 0.0);
                                if (wgt > 0.0)
                                    cnt[ph][(static_cast<size_t>(i) * gy + j) * gzd + k] +=
                                        wgt * wgt * wgt * val;
                            }
                }
            }
        }
    };
    auto sampleCnt = [&](const Vec& cu, const Vec& cv, const Vec& cw,
                         double px, double py, double pz) {
        const double su2 = sampleGrid(cu, NX + 1, NY, NZ, 0.0, 0.5, 0.5, px, py, pz);
        const double sv2 = sampleGrid(cv, NX, NY + 1, NZ, 0.5, 0.0, 0.5, px, py, pz);
        if (NZ == 1) return 0.5 * (su2 + sv2);
        const double sw2 = sampleGrid(cw, NX, NY, NZ + 1, 0.5, 0.5, 0.0, px, py, pz);
        return (su2 + sv2 + sw2) / 3.0;
    };

    std::vector<char> escFlag(np, 0);
    // emitted spray (the reference's spray flag): a massless liquid particle flagged
    // in spray_in, always a droplet, counted in no splat, deleted when it falls into
    // the water (its liquid fraction above 1 - esc_phi); the liquid fraction of the
    // other liquid particles kept for the emission below
    const bool sprayOn = droplet_model != 0 && escape != 0;
    const bool emitOn = sprayOn && emit_model != 0;
    std::vector<char> sprayF(np, 0), insideF(np, 0);
    std::vector<double> lfP(emitOn ? np : 0, 1.0);
    if ((escape != 0 || adOn) && np > 0) {
        std::vector<double> qx(np), qy(np), qz(np), scl(np);
        std::vector<char> isLiq(np);
#pragma omp parallel for schedule(static)
        for (long long pi = 0; pi < static_cast<long long>(np); ++pi) {
            qx[pi] = positions[pi].x; qy[pi] = positions[pi].y; qz[pi] = positions[pi].z;
            scl[pi] = scaleOf(pi);
            isLiq[pi] = (pi < particle_phase.size() ? particle_phase[pi] : 0.0f) >= 0.5f;
            sprayF[pi] = (sprayOn && static_cast<size_t>(pi) < spray_in.size()
                          && spray_in[pi] >= 0.5f && isLiq[pi]) ? 1 : 0;
        }
        Vec cntU[2] = {Vec(nu, 0.0), Vec(nu, 0.0)};
        Vec cntV[2] = {Vec(nv, 0.0), Vec(nv, 0.0)};
        Vec cntW[2] = {Vec(nw, 0.0), Vec(nw, 0.0)};
        countSplat(qx, qy, qz, isLiq, sprayF, scl, cntU, cntV, cntW);
        if (escape != 0) {
            const double thresh = (1.0 - esc_phi) * rho0_face;
#pragma omp parallel for schedule(static)
            for (long long pi = 0; pi < static_cast<long long>(np); ++pi) {
                const int other = isLiq[pi] ? 0 : 1;
                const double frac = sampleCnt(cntU[other], cntV[other], cntW[other],
                                              qx[pi], qy[pi], qz[pi]);
                escFlag[pi] = frac > thresh ? 1 : 0;
                // the wave zone's liquid is the zone's to move, never a droplet
                if (zoneOn && isLiq[pi] && qx[pi] < zoneLen) escFlag[pi] = 0;
                if (sprayF[pi]) {
                    escFlag[pi] = 1;
                    const double fl = sampleCnt(cntU[1], cntV[1], cntW[1], qx[pi], qy[pi], qz[pi]);
                    insideF[pi] = fl > thresh ? 1 : 0;
                } else if (emitOn && isLiq[pi]) {
                    lfP[pi] = sampleCnt(cntU[1], cntV[1], cntW[1], qx[pi], qy[pi], qz[pi]) /
                              static_cast<double>(rho0_face);
                }
            }
        }
    }

    Vec du(nu), dv(nv), dw(nw);
    const long long lnu = static_cast<long long>(nu), lnv = static_cast<long long>(nv);
    const long long lnw = static_cast<long long>(nw);
#pragma omp parallel for schedule(static)
    for (long long c = 0; c < lnu; ++c) du[c] = u[c] - uStar[c];
#pragma omp parallel for schedule(static)
    for (long long c = 0; c < lnv; ++c) dv[c] = v[c] - vStar[c];
#pragma omp parallel for schedule(static)
    for (long long c = 0; c < lnw; ++c) dw[c] = w[c] - wStar[c];

    auto sU = [&](const Vec& g, double x, double y, double zc) {
        return sampleGrid(g, NX + 1, NY, NZ, 0.0, 0.5, 0.5, x, y, zc); };
    auto sV = [&](const Vec& g, double x, double y, double zc) {
        return sampleGrid(g, NX, NY + 1, NZ, 0.5, 0.0, 0.5, x, y, zc); };
    auto sW = [&](const Vec& g, double x, double y, double zc) {
        return sampleGrid(g, NX, NY, NZ + 1, 0.5, 0.5, 0.0, x, y, zc); };

    const double lo = 0.51, hx = NX - 0.51, hy = NY - 0.51, hz = NZ - 0.51;
    const long long lnp = static_cast<long long>(np);

    // droplet radii (PF-FLIP Section 7; the reference's _droplet_radii): a liquid
    // particle the escape flag marks is a droplet; new, it draws its radius; every
    // step it evaporates by the d^2 law; below droplet_r_min it is gone (dropped
    // after the advection below). Rejoined, its radius returns to 0.
    const bool dropOn = droplet_model != 0;
    std::vector<double> radNew(np, 0.0);
    std::vector<char> gone(np, 0);
    if (dropOn) {
        const double rmed = droplet_r_median, rsig = droplet_r_sigma;
        const double evap = 0.25 * static_cast<double>(droplet_evap_k) * dt;
#pragma omp parallel for schedule(static)
        for (long long pi = 0; pi < lnp; ++pi) {
            const bool liquid = (static_cast<size_t>(pi) < particle_phase.size()
                                 ? particle_phase[pi] : 0.0f) >= 0.5f;
            if (!escFlag[pi] || !liquid) continue;
            double r = static_cast<size_t>(pi) < radius_in.size()
                ? static_cast<double>(radius_in[pi]) : 0.0;
            if (r <= 0.0)
                r = rmed * std::exp(rsig * dropletNormal(step_index, static_cast<uint64_t>(pi)));
            r = std::sqrt(std::max(r * r - evap, 0.0));
            radNew[pi] = r;
            gone[pi] = r < static_cast<double>(droplet_r_min) ? 1 : 0;
            if (sprayF[pi] && insideF[pi]) gone[pi] = 1;
        }
    }

    // the droplets' own integrator (spray_cfl > 0): gravity and drag applied in the
    // advection below, sub-stepped; here they only keep their velocity
    const bool flyOn = dropOn && spray_cfl > 0.0f;
    const double sprayCfl = spray_cfl;
    std::vector<char> flyP(np, 0);
    std::vector<double> nvx(np), nvy(np), nvz(np);
#pragma omp parallel for schedule(static)
    for (long long pi = 0; pi < lnp; ++pi) {
        const double px = positions[pi].x, py = positions[pi].y, pz = positions[pi].z;
        double vx = pi < static_cast<long long>(velocities.size()) ? velocities[pi].x : 0.0;
        double vy = pi < static_cast<long long>(velocities.size()) ? velocities[pi].y : 0.0;
        double vz = pi < static_cast<long long>(velocities.size()) ? velocities[pi].z : 0.0;
        const bool liquid = (static_cast<size_t>(pi) < particle_phase.size()
                             ? particle_phase[pi] : 0.0f) >= 0.5f;
        double a = liquid ? alpha_liquid : alpha_air;
        if (adOn) {
            const double sc = scaleOf(pi);
            a = 1.0 - (1.0 - a) / (sc * sc);
        }
        if (!escFlag[pi]) {
            vx = a * (vx + sU(du, px, py, pz)) + (1.0 - a) * sU(u, px, py, pz);
            vy = a * (vy + sV(dv, px, py, pz)) + (1.0 - a) * sV(v, px, py, pz);
            vz = flat ? 0.0 : a * (vz + sW(dw, px, py, pz)) + (1.0 - a) * sW(w, px, py, pz);
            if (zoneOn && liquid) {
                // the wave zone: the water's velocity blended to the trains' (above)
                double ut, vt, wt;
                const double wgt = zoneAt(px, py, pz, ut, vt, wt);
                if (wgt > 0.0) {
                    vx = (1.0 - wgt) * vx + wgt * ut;
                    vy = (1.0 - wgt) * vy + wgt * vt;
                    vz = flat ? 0.0 : (1.0 - wgt) * vz + wgt * wt;
                }
            }
        } else if (flyOn && liquid) {
            flyP[pi] = 1;
        } else {
            const double gx = sU(u, px, py, pz), gy = sV(v, px, py, pz);
            const double gz = flat ? 0.0 : sW(w, px, py, pz);
            vy += (liquid ? gravity : buoyancy * std::fabs(gravity)) * dt;
            if (dropOn && liquid) {
                // the droplet's own drag (the reference's _droplet_drag): the sphere's
                // relaxation rate (3/8)(rho_air/rho_liquid)(C_d/r)|v_rel|, Schiller-
                // Naumann C_d(Re) below Re 1000 and 0.44 above, implicit in the rate
                const double r = std::max(radNew[pi], static_cast<double>(droplet_r_min));
                const double mpc = metres_per_cell, nu = air_nu;
                const double rx = (gx - vx) * mpc, ry = (gy - vy) * mpc;
                const double rz = flat ? 0.0 : (gz - vz) * mpc;
                const double speed = std::sqrt(rx * rx + ry * ry + rz * rz);
                const double re = 2.0 * r * speed / nu;
                const double cdv = re <= 1000.0
                    ? 12.0 * nu / r * (1.0 + 0.15 * std::pow(re, 0.687)) : 0.44 * speed;
                const double k = 0.375 * (static_cast<double>(rho_air) / rho_liquid) * cdv / r * dt;
                vx = (vx + k * gx) / (1.0 + k);
                vy = (vy + k * gy) / (1.0 + k);
                vz = flat ? 0.0 : (vz + k * gz) / (1.0 + k);
            } else {
                const double drag = liquid ? drag_droplet : drag_bubble;
                vx += drag * dt * (gx - vx);
                vy += drag * dt * (gy - vy);
                vz = flat ? 0.0 : vz + drag * dt * (gz - vz);
            }
        }
        nvx[pi] = vx; nvy[pi] = vy; nvz[pi] = vz;
    }

    std::vector<double> tauNew(np, 0.0), wtNewV(np, 1.0);
    if (stOn) {
        std::mt19937 rng(static_cast<unsigned>(st_seed) * 2654435761u ^
                         static_cast<unsigned>(step_index) * 2246822519u);
        std::uniform_real_distribution<double> uni(-0.5, 0.5);
        // Keep the original RNG stream in input order. Only the independent
        // velocity attenuation and eight-point normalization run in parallel;
        // per-thread RNGs would change trajectories with the worker count.
        for (size_t pi = 0; pi < np; ++pi) tauNew[pi] = uni(rng);
#pragma omp parallel for schedule(static)
        for (long long pi = 0; pi < lnp; ++pi) {
            const double sp = std::max({std::fabs(nvx[pi]), std::fabs(nvy[pi]),
                                        std::fabs(nvz[pi])});
            const double c = std::clamp(sp * dt, 0.0, 1.0);
            const double xi = c * c * (3.0 - 2.0 * c);
            const double tp = xi * tauNew[pi];
            tauNew[pi] = tp;
            wtNewV[pi] = wtKernel(tp) / wtNorm(xi);
        }
    }

    // a droplet's flight over `total` seconds (the reference's _droplet_fly): sub-steps
    // of at most sprayCfl cells; over each the drag relaxation toward the air
    // velocity sampled at its start is taken exactly with gravity (rate k =
    // (3/8)(rho_air/rho_liquid)(C_d/r)|v_rel|, Schiller-Naumann, at the sub-step's
    // midpoint slip), the position by the trapezoid of the sub-step's two velocities
    auto flyDroplet = [&](double& px, double& py, double& pz, double& vx, double& vy,
                          double& vz, double total, double r) {
        const double mpc = metres_per_cell, nu = air_nu;
        const double rr = static_cast<double>(rho_air) / rho_liquid;
        double rem = total;
        for (int it2 = 0; rem > 1e-12 && it2 < 1000000; ++it2) {
            const double cx = std::clamp(px, lo, hx), cy = std::clamp(py, lo, hy);
            const double cz = flat ? pz : std::clamp(pz, lo, hz);
            double gx = sU(u, cx, cy, cz), gy = sV(v, cx, cy, cz);
            double gz = flat ? 0.0 : sW(w, cx, cy, cz);
            double speed = std::max({std::fabs(vx), std::fabs(vy), std::fabs(vz), 1e-9});
            if (droplet_model==2) speed=std::max({speed,std::fabs(gx),std::fabs(gy),std::fabs(gz)});
            double h = std::min(rem, sprayCfl / speed);
            if (droplet_model==2) h=std::min(h,static_cast<double>(spray_dt_max));
            // the relaxation rate at a slip (the reference's rate)
            auto rate = [&](double ax, double ay, double az) {
                const double rx = (gx - ax) * mpc, ry = (gy - ay) * mpc;
                const double rz = flat ? 0.0 : (gz - az) * mpc;
                const double sp = std::sqrt(rx * rx + ry * ry + rz * rz);
                const double re = 2.0 * r * sp / nu;
                const double cdv = re <= 1000.0
                    ? 12.0 * nu / r * (1.0 + 0.15 * std::pow(re, 0.687)) : 0.44 * sp;
                return 0.375 * rr * cdv / r;
            };
            // the velocity after hh seconds at rate k toward the air under gravity
            auto relax = [&](double k, double hh, double& ox, double& oy, double& oz) {
                const double e = std::exp(-k * hh);
                const double gk = k * hh > 1e-6 ? (1.0 - e) / std::max(k, 1e-300)
                                                : hh * (1.0 - 0.5 * k * hh);
                ox = gx + (vx - gx) * e;
                oy = gy + (vy - gy) * e + gravity * gk;
                oz = flat ? 0.0 : gz + (vz - gz) * e;
            };
            // the rate taken at the sub-step's midpoint: a half step at the starting
            // rate finds it, then the whole step at that rate
            double mx, my, mz, nx2, ny2, nz2;
            relax(rate(vx, vy, vz), 0.5 * h, mx, my, mz);
            if (droplet_model==2) {
                const double qx=std::clamp(px+0.25*h*(vx+mx),lo,hx);
                const double qy=std::clamp(py+0.25*h*(vy+my),lo,hy);
                const double qz=flat ? pz : std::clamp(pz+0.25*h*(vz+mz),lo,hz);
                gx=sU(u,qx,qy,qz); gy=sV(v,qx,qy,qz); gz=flat ? 0 : sW(w,qx,qy,qz);
            }
            relax(rate(mx, my, mz), h, nx2, ny2, nz2);
            px += 0.5 * h * (vx + nx2);
            py += 0.5 * h * (vy + ny2);
            if (!flat) pz += 0.5 * h * (vz + nz2);
            vx = nx2; vy = ny2; vz = nz2;
            rem -= h;
        }
    };

    // working vectors sized np, then merge/split edits them.
    std::vector<double> Px(np), Py(np), Pz(np), Vx(np), Vy(np), Vz(np);
    std::vector<double> Tau(np), Wt(np), Scl(np), Ph01(np);
    std::vector<char> Esc(np), Spray(np);
    std::vector<double> Rad(np);
    // PFFLIP_SOLID_DEBUG=1 reports how many particles the level set pushed out
    static const bool lsDebug = []() {
        const char* e = std::getenv("PFFLIP_SOLID_DEBUG");
        return e && e[0] == '1';
    }();
    long long lsPushed = 0;
#pragma omp parallel for schedule(dynamic, 1024) reduction(+ : lsPushed)
    for (long long pi = 0; pi < lnp; ++pi) {
        double px = positions[pi].x, py = positions[pi].y, pz = positions[pi].z;
        double vx = nvx[pi], vy = nvy[pi], vz = nvz[pi];
        const bool liquid = (static_cast<size_t>(pi) < particle_phase.size()
                             ? particle_phase[pi] : 0.0f) >= 0.5f;
        const double tauOld = static_cast<size_t>(pi) < tau_in.size()
            ? static_cast<double>(tau_in[pi]) : 0.0;
        const double dtAct = stOn
            ? std::clamp(dt * (1.0 + tauNew[pi] - tauOld), 0.0, 2.0 * dt) : dt;
        if (escFlag[pi] && flyP[pi]) {
            flyDroplet(px, py, pz, vx, vy, vz, dtAct,
                       std::max(radNew[pi], static_cast<double>(droplet_r_min)));
        } else if (escFlag[pi]) {
            px += dtAct * vx; py += dtAct * vy;
            if (!flat) pz += dtAct * vz;
        } else if (stOn) {
            double rem = dtAct;
            for (int r = 0; rem > 1e-12 && r < 64; ++r) {
                const double v1x = sU(u, px, py, pz), v1y = sV(v, px, py, pz);
                const double v1z = flat ? 0.0 : sW(w, px, py, pz);
                const double sp = std::max({std::fabs(v1x), std::fabs(v1y),
                                            std::fabs(v1z), 1e-9});
                const double dts = std::min(rem, 1.0 / sp);
                const double mx = px + 0.5 * dts * v1x, my = py + 0.5 * dts * v1y;
                const double mz = flat ? pz : pz + 0.5 * dts * v1z;
                const double cx = std::max(mx, 0.51), cy = std::max(my, 0.51);
                const double cz = flat ? pz : std::max(mz, 0.51);
                px += dts * sU(u, cx, cy, cz);
                py += dts * sV(v, cx, cy, cz);
                if (!flat) pz += dts * sW(w, cx, cy, cz);
                rem -= dts;
            }
        } else {
            const double mx = px + 0.5 * dt * sU(u, px, py, pz);
            const double my = py + 0.5 * dt * sV(v, px, py, pz);
            const double mz = flat ? pz : pz + 0.5 * dt * sW(w, px, py, pz);
            const double cx = std::max(mx, 0.51), cy = std::max(my, 0.51);
            const double cz = flat ? pz : std::max(mz, 0.51);
            px += dt * sU(u, cx, cy, cz);
            py += dt * sV(v, cx, cy, cz);
            if (!flat) pz += dt * sW(w, cx, cy, cz);
        }
        if (px < lo || px > hx) { px = std::clamp(px, lo, hx); vx = 0.0; }
        if (py < lo || py > hy) { py = std::clamp(py, lo, hy); vy = 0.0; }
        if (!flat && (pz < lo || pz > hz)) { pz = std::clamp(pz, lo, hz); vz = 0.0; }
        for (size_t b = 0; b < nObs; ++b) {
            const auto& m0 = obstacles_min[b];
            const auto& m1 = obstacles_max[b];
            const bool inZ = flat || (pz > m0.z && pz < m1.z);
            if (px > m0.x && px < m1.x && py > m0.y && py < m1.y && inZ) {
                double pen[6] = {px - m0.x, m1.x - px, py - m0.y, m1.y - py,
                                 flat ? 1e30 : pz - m0.z, flat ? 1e30 : m1.z - pz};
                int side = 0;
                for (int q2 = 1; q2 < 6; ++q2) if (pen[q2] < pen[side]) side = q2;
                if (side == 0) { px = m0.x - 0.01; vx = 0.0; }
                else if (side == 1) { px = m1.x + 0.01; vx = 0.0; }
                else if (side == 2) { py = m0.y - 0.01; vy = 0.0; }
                else if (side == 3) { py = m1.y + 0.01; vy = 0.0; }
                else if (side == 4) { pz = m0.z - 0.01; vz = 0.0; }
                else { pz = m1.z + 0.01; vz = 0.0; }
            }
        }
        if (lsOn && pushP(px, py, pz, &vx, &vy, &vz))
            ++lsPushed;
        Px[pi] = px; Py[pi] = py; Pz[pi] = pz;
        Vx[pi] = vx; Vy[pi] = vy; Vz[pi] = vz;
        Tau[pi] = tauNew[pi]; Wt[pi] = wtNewV[pi];
        Scl[pi] = scaleOf(pi); Ph01[pi] = liquid ? 1.0 : 0.0;
        Esc[pi] = escFlag[pi]; Rad[pi] = radNew[pi]; Spray[pi] = sprayF[pi];
    }

    if (lsDebug)
        std::fprintf(stderr, "step_3d: level set %s; pushed %lld of %zu particles out\n",
                     lsOn ? "on" : "off", lsPushed, np);

    // extra spray emission (PF-FLIP Section 7 after Ihmsen et al. 2012; the
    // reference's _emit): a liquid particle at the interface (its liquid fraction
    // under 1 - esc_phi) in a cell whose velocity gradient exceeds
    // emit_grad_threshold emits, with probability min(1, emit_rate (G - G0) dt) from
    // the hash of (step_index, index), a massless spray particle at its position and
    // velocity, at most emit_max_per_cell a cell, in index order
    long long nEmitted = 0;
    long long nCandidates=0,nEvents=0,nRequested=0,nCap=0,nGeometry=0;
    double rawR3=0,clippedR3=0;
    long long oldSpray=0;
    if (emit_model==2) for (size_t i=0;i<np;++i) oldSpray+=Spray[i]!=0;
    if (emitOn && np > 0) {
        // the Frobenius norm of the Jacobian of the cell-centred velocity, per cell:
        // central differences inside, one-sided at the walls, a cell the unit
        Vec uc(ncell), vc(ncell), wc(ncell);
#pragma omp parallel for schedule(static)
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k < NZ; ++k) {
                    const size_t c = L.ic(i, j, k);
                    uc[c] = 0.5 * (u[L.iu(i, j, k)] + u[L.iu(i + 1, j, k)]);
                    vc[c] = 0.5 * (v[L.iv(i, j, k)] + v[L.iv(i, j + 1, k)]);
                    wc[c] = flat ? 0.0 : 0.5 * (w[L.iw(i, j, k)] + w[L.iw(i, j, k + 1)]);
                }
        auto dAxis = [&](const Vec& f, int i, int j, int k, int ax) {
            const int n = ax == 0 ? NX : (ax == 1 ? NY : NZ);
            const int q = ax == 0 ? i : (ax == 1 ? j : k);
            if (n < 2) return 0.0;
            auto at = [&](int d) {
                return f[L.ic(ax == 0 ? i + d : i, ax == 1 ? j + d : j, ax == 2 ? k + d : k)];
            };
            if (q == 0) return at(1) - at(0);
            if (q == n - 1) return at(0) - at(-1);
            return 0.5 * (at(1) - at(-1));
        };
        Vec gradN(ncell, 0.0);
#pragma omp parallel for schedule(static)
        for (int i = 0; i < NX; ++i)
            for (int j = 0; j < NY; ++j)
                for (int k = 0; k < NZ; ++k) {
                    double g2 = 0.0;
                    for (int ax = 0; ax < 3; ++ax)
                        for (const Vec* f : {&uc, &vc, &wc}) {
                            const double d = dAxis(*f, i, j, k, ax);
                            g2 += d * d;
                        }
                    gradN[L.ic(i, j, k)] = std::sqrt(g2);
                }
        const double thr = emit_grad_threshold, rate = emit_rate;
        const double lfLine = 1.0 - static_cast<double>(esc_phi);
        std::unordered_map<size_t, int> perCell;
        std::vector<size_t> hits;
        for (size_t pi = 0; pi < np; ++pi) {
            if (Ph01[pi] < 0.5 || Spray[pi] || Esc[pi] || lfP[pi] >= lfLine) continue;
            const int ci = std::clamp(static_cast<int>(std::floor(Px[pi])), 0, NX - 1);
            const int cj = std::clamp(static_cast<int>(std::floor(Py[pi])), 0, NY - 1);
            const int ck = flat ? 0 : std::clamp(static_cast<int>(std::floor(Pz[pi])), 0, NZ - 1);
            const size_t cell = L.ic(ci, cj, ck);
            const double g = gradN[cell];
            if (g <= thr) continue;
            if (emit_model==2) {
                ++nCandidates;
                const double lam=rate*(g-thr)*dt;
                if (!std::isfinite(lam) || lam>1e6) { refuseSpray(); return; }
                const auto key=SprayEvents::key(emit_seed,step_index,pi);
                const int ne=static_cast<int>(std::floor(lam)) +
                    (SprayEvents::uniform(key,0)<lam-std::floor(lam));
                nEvents+=ne; nRequested+=static_cast<long long>(ne)*emit_cloud_size;
                // Live budget includes deaths awaiting this step's compaction.
                for (int ev=0; ev<ne; ++ev) {
                    if (perCell[cell]>=emit_max_per_cell || oldSpray+nEmitted>=emit_max_active) {
                        nCap+=static_cast<long long>(ne-ev)*emit_cloud_size; break;
                    }
                    const double speed=std::sqrt(Vx[pi]*Vx[pi]+Vy[pi]*Vy[pi]+Vz[pi]*Vz[pi]);
                    for (int child=0; child<emit_cloud_size; ++child) {
                        const auto k=SprayEvents::key(emit_seed,step_index,pi,ev,child);
                        const auto dp=SprayEvents::ball(k,3,flat), dv=SprayEvents::ball(k,6,flat);
                        const double x=Px[pi]+emit_spread*dp[0], y=Py[pi]+emit_spread*dp[1];
                        const double z=flat ? Pz[pi] : Pz[pi]+emit_spread*dp[2];
                        bool bad=x<lo||x>hx||y<lo||y>hy||(!flat&&(z<lo||z>hz));
                        for (size_t b=0;b<nObs&&!bad;++b) {
                            const auto& a=obstacles_min[b]; const auto& c=obstacles_max[b];
                            bad=x>a.x&&x<c.x&&y>a.y&&y<c.y&&(flat||(z>a.z&&z<c.z));
                        }
                        if (lsOn && !bad) bad=(pistonOn ?
                            Solids::phiValue(phiAll,NX,NY,NZ,flat,x,y,z) :
                            Solids::phiValue(solid_phi,NX,NY,NZ,flat,x,y,z))<0;
                        if (bad) { ++nGeometry; continue; }
                        if (perCell[cell]>=emit_max_per_cell || oldSpray+nEmitted>=emit_max_active) {
                            ++nCap; continue;
                        }
                        const double raw=SprayEvents::rawRadius(k,emit_cloud_size,droplet_r_sigma,emit_event_radius);
                        const double r=std::clamp(raw,static_cast<double>(droplet_r_min),static_cast<double>(emit_radius_max));
                        Px.push_back(x); Py.push_back(y); Pz.push_back(z);
                        Vx.push_back(Vx[pi]+emit_velocity_spread*speed*dv[0]);
                        Vy.push_back(Vy[pi]+emit_velocity_spread*speed*dv[1]);
                        Vz.push_back(flat ? 0 : Vz[pi]+emit_velocity_spread*speed*dv[2]);
                        Tau.push_back(0); Wt.push_back(stOn ? wtKernel(0)/wtNorm(1):1);
                        Scl.push_back(1); Ph01.push_back(1); Esc.push_back(1);
                        Rad.push_back(r); Spray.push_back(1); gone.push_back(0);
                        ++perCell[cell]; ++nEmitted; rawR3+=raw*raw*raw; clippedR3+=r*r*r;
                    }
                }
                continue;
            }
            const double prob = std::clamp(rate * (g - thr) * dt, 0.0, 1.0);
            if (dropletUniform(step_index, static_cast<uint64_t>(pi)) >= prob) continue;
            int& cnt = perCell[cell];
            if (cnt >= emit_max_per_cell) continue;
            ++cnt;
            hits.push_back(pi);
        }
        for (size_t pi : hits) {
            const double px = Px[pi], py = Py[pi], pz = Pz[pi];
            const double vx = Vx[pi], vy = Vy[pi], vz = Vz[pi];
            Px.push_back(px); Py.push_back(py); Pz.push_back(pz);
            Vx.push_back(vx); Vy.push_back(vy); Vz.push_back(vz);
            Tau.push_back(0.0); Wt.push_back(stOn ? wtKernel(0.0) / wtNorm(1.0) : 1.0);
            Scl.push_back(1.0); Ph01.push_back(1.0); Esc.push_back(1);
            Rad.push_back(0.0); Spray.push_back(1);
            gone.push_back(0);
        }
        if (emit_model!=2) nEmitted = static_cast<long long>(hits.size());
    }
    if (emit_model==2 && std::getenv("PFFLIP_SPRAY_AUDIT")) {
        long long evap=0,reentry=0,alive=0;
        for (size_t i=0;i<np;++i) if (Spray[i]) {
            if (insideF[i]) ++reentry;
            else if (gone[i]) ++evap;
            else ++alive;
        }
        std::fprintf(stderr,"SPRAY_EVENT {\"step\":%d,\"candidates\":%lld,\"events\":%lld,\"requested\":%lld,\"accepted\":%lld,\"rejected_cap\":%lld,\"rejected_geometry\":%lld,\"raw_r3\":%.17g,\"clipped_r3\":%.17g,\"evaporated\":%lld,\"reentered\":%lld,\"alive\":%lld}\n",
            step_index,nCandidates,nEvents,nRequested,nEmitted,nCap,nGeometry,rawR3,clippedR3,evap,reentry,alive+nEmitted);
    }
    if (lsDebug && emitOn)
        std::fprintf(stderr, "step_3d: emitted %lld spray particles\n", nEmitted);

    if (dropOn) {
        // droplets evaporated below droplet_r_min leave the simulation here, before
        // the adaptive merge and split see the particle set
        bool any = false;
        for (size_t i = 0; i < gone.size(); ++i) if (gone[i]) { any = true; break; }
        if (any) {
            auto dropGone = [&](auto& vec) {
                size_t wr = 0;
                for (size_t i = 0; i < gone.size(); ++i)
                    if (!gone[i]) vec[wr++] = vec[i];
                vec.resize(wr);
            };
            dropGone(Px); dropGone(Py); dropGone(Pz);
            dropGone(Vx); dropGone(Vy); dropGone(Vz);
            dropGone(Tau); dropGone(Wt); dropGone(Scl); dropGone(Ph01);
            dropGone(Rad); dropGone(Esc); dropGone(Spray);
        }
    }

    if (adOn && np > 0) {
        // fresh liquid-count field at the new positions, as the reference does.
        std::vector<char> isLiq(Px.size());
#pragma omp parallel for schedule(static)
        for (long long i = 0; i < static_cast<long long>(Px.size()); ++i)
            isLiq[i] = Ph01[i] >= 0.5 ? 1 : 0;
        Vec cntU[2] = {Vec(nu, 0.0), Vec(nu, 0.0)};
        Vec cntV[2] = {Vec(nv, 0.0), Vec(nv, 0.0)};
        Vec cntW[2] = {Vec(nw, 0.0), Vec(nw, 0.0)};
        countSplat(Px, Py, Pz, isLiq, Spray, Scl, cntU, cntV, cntW);
        std::vector<double> lf(Px.size());
#pragma omp parallel for schedule(static)
        for (long long i = 0; i < static_cast<long long>(Px.size()); ++i)
            lf[i] = sampleCnt(cntU[1], cntV[1], cntW[1], Px[i], Py[i], Pz[i]) /
                    std::max(static_cast<double>(rho0_face), 1e-9);

        const int groupN = flat ? 4 : 8;
        const double wallM = 2.0;
        // ---- split coarse near interface or walls
        std::vector<size_t> splitIdx;
        for (size_t i = 0; i < Px.size(); ++i) {
            if (Scl[i] < 1.5) continue;
            const bool nearWall = Px[i] < wallM || Px[i] > NX - wallM ||
                                  Py[i] < wallM || Py[i] > NY - wallM ||
                                  (!flat && (Pz[i] < wallM || Pz[i] > NZ - wallM));
            if (lf[i] > 0.05 || nearWall) splitIdx.push_back(i);
        }
        if (!splitIdx.empty()) {
            std::vector<char> drop(Px.size(), 0);
            for (size_t i : splitIdx) drop[i] = 1;
            auto keepFilter = [&](auto& vec) {
                size_t wr = 0;
                for (size_t i = 0; i < drop.size(); ++i)
                    if (!drop[i]) vec[wr++] = vec[i];
                vec.resize(wr);
            };
            std::vector<double> sp_[10];
            // stash split parents before filtering
            std::vector<double> pxs, pys, pzs, vxs, vys, vzs, tas, phs;
            for (size_t i : splitIdx) {
                pxs.push_back(Px[i]); pys.push_back(Py[i]); pzs.push_back(Pz[i]);
                vxs.push_back(Vx[i]); vys.push_back(Vy[i]); vzs.push_back(Vz[i]);
                tas.push_back(Tau[i]); phs.push_back(Ph01[i]);
            }
            keepFilter(Px); keepFilter(Py); keepFilter(Pz);
            keepFilter(Vx); keepFilter(Vy); keepFilter(Vz);
            keepFilter(Tau); keepFilter(Wt); keepFilter(Scl); keepFilter(Ph01);
            keepFilter(Rad); keepFilter(Spray);
            // The merge criterion belongs to the same particles as the state.
            // Keeping old indices here let newly emitted spray shift which air
            // particles merged, despite spray having zero transfer weight.
            keepFilter(lf);
            {
                size_t wr = 0;
                for (size_t i = 0; i < drop.size(); ++i)
                    if (!drop[i]) Esc[wr++] = Esc[i];
                Esc.resize(wr);
            }
            for (size_t m2 = 0; m2 < pxs.size(); ++m2) {
                for (int ci = 0; ci < groupN; ++ci) {
                    const double ox = (ci & 1) ? 0.5 : -0.5;
                    const double oy = (ci & 2) ? 0.5 : -0.5;
                    const double oz = flat ? 0.0 : ((ci & 4) ? 0.5 : -0.5);
                    Px.push_back(std::clamp(pxs[m2] + ox, lo, hx));
                    Py.push_back(std::clamp(pys[m2] + oy, lo, hy));
                    Pz.push_back(flat ? pzs[m2] : std::clamp(pzs[m2] + oz, lo, hz));
                    Vx.push_back(vxs[m2]); Vy.push_back(vys[m2]); Vz.push_back(vzs[m2]);
                    Tau.push_back(tas[m2]);
                    Wt.push_back(wtKernel(tas[m2]) / wtNorm(1.0));
                    Scl.push_back(1.0); Ph01.push_back(phs[m2]); Esc.push_back(0);
                    Rad.push_back(0.0); Spray.push_back(0);
                    lf.push_back(1.0); // split children stay fine for this step
                }
            }
        }
        // ---- merge deep-air fine, groupN-to-1 per 2-cell block
        std::vector<size_t> cand;
        for (size_t i = 0; i < Px.size(); ++i)
            if (Scl[i] < 1.5 && Ph01[i] < 0.5 &&
                (i < lf.size() ? lf[i] : 1.0) < 0.02)
                cand.push_back(i);
        if (static_cast<int>(cand.size()) >= groupN) {
            std::stable_sort(cand.begin(), cand.end(), [&](size_t a2, size_t b2) {
                const long long ka =
                    (static_cast<long long>(Px[a2] / 2) * 100000 +
                     static_cast<long long>(Py[a2] / 2)) * 100000 +
                    (flat ? 0 : static_cast<long long>(Pz[a2] / 2));
                const long long kb =
                    (static_cast<long long>(Px[b2] / 2) * 100000 +
                     static_cast<long long>(Py[b2] / 2)) * 100000 +
                    (flat ? 0 : static_cast<long long>(Pz[b2] / 2));
                return ka < kb;
            });
            auto key = [&](size_t i) {
                return (static_cast<long long>(Px[i] / 2) * 100000 +
                        static_cast<long long>(Py[i] / 2)) * 100000 +
                       (flat ? 0 : static_cast<long long>(Pz[i] / 2));
            };
            std::vector<char> drop(Px.size(), 0);
            std::vector<double> mpx, mpy, mpz, mvx, mvy, mvz, mta;
            size_t i = 0;
            while (i < cand.size()) {
                size_t j = i;
                while (j < cand.size() && key(cand[j]) == key(cand[i])) ++j;
                const size_t take = ((j - i) / groupN) * groupN;
                for (size_t g2 = 0; g2 < take; g2 += groupN) {
                    double ax = 0, ay = 0, az = 0, bx = 0, by = 0, bz = 0, tt = 0;
                    for (int m2 = 0; m2 < groupN; ++m2) {
                        const size_t id = cand[i + g2 + m2];
                        drop[id] = 1;
                        ax += Px[id]; ay += Py[id]; az += Pz[id];
                        bx += Vx[id]; by += Vy[id]; bz += Vz[id];
                        tt += Tau[id];
                    }
                    mpx.push_back(ax / groupN); mpy.push_back(ay / groupN);
                    mpz.push_back(az / groupN);
                    mvx.push_back(bx / groupN); mvy.push_back(by / groupN);
                    mvz.push_back(bz / groupN);
                    mta.push_back(tt / groupN);
                }
                i = j;
            }
            if (!mpx.empty()) {
                auto keepFilter = [&](auto& vec) {
                    size_t wr = 0;
                    for (size_t k2 = 0; k2 < drop.size(); ++k2)
                        if (!drop[k2]) vec[wr++] = vec[k2];
                    vec.resize(wr);
                };
                keepFilter(Px); keepFilter(Py); keepFilter(Pz);
                keepFilter(Vx); keepFilter(Vy); keepFilter(Vz);
                keepFilter(Tau); keepFilter(Wt); keepFilter(Scl); keepFilter(Ph01);
                keepFilter(Rad); keepFilter(Spray);
                {
                    size_t wr = 0;
                    for (size_t k2 = 0; k2 < drop.size(); ++k2)
                        if (!drop[k2]) Esc[wr++] = Esc[k2];
                    Esc.resize(wr);
                }
                for (size_t m2 = 0; m2 < mpx.size(); ++m2) {
                    Px.push_back(mpx[m2]); Py.push_back(mpy[m2]); Pz.push_back(mpz[m2]);
                    Vx.push_back(mvx[m2]); Vy.push_back(mvy[m2]); Vz.push_back(mvz[m2]);
                    Tau.push_back(mta[m2]);
                    Wt.push_back(wtKernel(mta[m2]) / wtNorm(1.0));
                    Scl.push_back(2.0); Ph01.push_back(0.0); Esc.push_back(0);
                    Rad.push_back(0.0); Spray.push_back(0);
                }
            }
        }
    }

    // ---- final obstacle sweep: merge centroids and split children may land
    // inside a box; push every particle out along least penetration.
    if (nObs > 0) {
#pragma omp parallel for schedule(static)
        for (long long i = 0; i < static_cast<long long>(Px.size()); ++i) {
            for (size_t b = 0; b < nObs; ++b) {
                const auto& m0 = obstacles_min[b];
                const auto& m1 = obstacles_max[b];
                const bool inZ = flat || (Pz[i] > m0.z && Pz[i] < m1.z);
                if (Px[i] > m0.x && Px[i] < m1.x && Py[i] > m0.y && Py[i] < m1.y && inZ) {
                    double pen[6] = {Px[i] - m0.x, m1.x - Px[i], Py[i] - m0.y,
                                     m1.y - Py[i], flat ? 1e30 : Pz[i] - m0.z,
                                     flat ? 1e30 : m1.z - Pz[i]};
                    int side = 0;
                    for (int q2 = 1; q2 < 6; ++q2) if (pen[q2] < pen[side]) side = q2;
                    if (side == 0) { Px[i] = m0.x - 0.01; Vx[i] = 0.0; }
                    else if (side == 1) { Px[i] = m1.x + 0.01; Vx[i] = 0.0; }
                    else if (side == 2) { Py[i] = m0.y - 0.01; Vy[i] = 0.0; }
                    else if (side == 3) { Py[i] = m1.y + 0.01; Vy[i] = 0.0; }
                    else if (side == 4) { Pz[i] = m0.z - 0.01; Vz[i] = 0.0; }
                    else { Pz[i] = m1.z + 0.01; Vz[i] = 0.0; }
                }
            }
        }
    }
    // the level-set solid's sweep runs only when particles merged or split:
    // every other particle was pushed out right after its advection above
    if (lsOn && adOn) {
#pragma omp parallel for schedule(static)
        for (long long i = 0; i < static_cast<long long>(Px.size()); ++i)
            pushP(Px[i], Py[i], Pz[i], &Vx[i], &Vy[i], &Vz[i]);
    }

    // ---- pack outputs at the final particle count.
    // the emitted spray kept behind every other particle, so the others keep the
    // order (and the ST jitter stream) they would have without it
    if (sprayOn) {
        bool any = false;
        for (char s2 : Spray) if (s2) { any = true; break; }
        if (any) {
            std::vector<size_t> order;
            order.reserve(Px.size());
            for (size_t i = 0; i < Px.size(); ++i) if (!Spray[i]) order.push_back(i);
            for (size_t i = 0; i < Px.size(); ++i) if (Spray[i]) order.push_back(i);
            auto permute = [&](auto& vec) {
                auto copy = vec;
                for (size_t i = 0; i < order.size(); ++i) vec[i] = copy[order[i]];
            };
            permute(Px); permute(Py); permute(Pz); permute(Vx); permute(Vy); permute(Vz);
            permute(Tau); permute(Wt); permute(Scl); permute(Ph01); permute(Rad);
            permute(Esc); permute(Spray);
        }
    }
    const size_t nf = Px.size();
    auto outPos = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>(nf);
    auto outVel = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>(nf);
    auto outSync = Amino::newMutablePtr<Amino::Array<Bifrost::Math::float3>>(nf);
    auto outM = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outMx = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outMy = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outMz = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outTau = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outScl = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outPh = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outWt = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outWtPh = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outEsc = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outRad = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outSpray = Amino::newMutablePtr<Amino::Array<float>>(nf);
    auto outMF = Amino::newMutablePtr<Amino::Array<bool>>(nf);
    auto outMC = Amino::newMutablePtr<Amino::Array<bool>>(nf);
    double maxSpeed = 0.0;
    for (size_t i = 0; i < nf; ++i) {
        (*outPos)[i] = {static_cast<float>(Px[i]), static_cast<float>(Py[i]),
                        static_cast<float>(Pz[i])};
        (*outVel)[i] = {static_cast<float>(Vx[i]), static_cast<float>(Vy[i]),
                        static_cast<float>(Vz[i])};
        const double sdt = Tau[i] * dt;
        // The time-resync can carry a display position back inside a box the
        // sweep above just pushed the actual position out of; clamp the
        // display position too (position only - no state or velocity change).
        double sx = Px[i] - sdt * Vx[i];
        double sy = Py[i] - sdt * Vy[i];
        double sz = flat ? Pz[i] : Pz[i] - sdt * Vz[i];
        for (size_t b = 0; b < nObs; ++b) {
            const auto& m0 = obstacles_min[b];
            const auto& m1 = obstacles_max[b];
            const bool inZ = flat || (sz > m0.z && sz < m1.z);
            if (sx > m0.x && sx < m1.x && sy > m0.y && sy < m1.y && inZ) {
                double pen[6] = {sx - m0.x, m1.x - sx, sy - m0.y, m1.y - sy,
                                 flat ? 1e30 : sz - m0.z, flat ? 1e30 : m1.z - sz};
                int side = 0;
                for (int q2 = 1; q2 < 6; ++q2) if (pen[q2] < pen[side]) side = q2;
                if (side == 0) sx = m0.x - 0.01;
                else if (side == 1) sx = m1.x + 0.01;
                else if (side == 2) sy = m0.y - 0.01;
                else if (side == 3) sy = m1.y + 0.01;
                else if (side == 4) sz = m0.z - 0.01;
                else sz = m1.z + 0.01;
            }
        }
        // (without a time jitter the display position is the actual one, already out)
        if (lsOn && sdt != 0.0)
            pushP(sx, sy, sz, nullptr, nullptr, nullptr);
        (*outSync)[i] = {static_cast<float>(sx), static_cast<float>(sy),
                         static_cast<float>(sz)};
        const double rhoP = Ph01[i] >= 0.5 ? rl : rg;
        const double mEff = Spray[i] ? 0.0 : rhoP * std::pow(Scl[i], dExp) *
                            (Scl[i] > 1.5 ? coarse_gain : 1.0);
        const double wp = Wt[i];
        (*outM)[i] = static_cast<float>(wp * mEff);
        (*outMx)[i] = static_cast<float>(wp * mEff * Vx[i]);
        (*outMy)[i] = static_cast<float>(wp * mEff * Vy[i]);
        (*outMz)[i] = static_cast<float>(wp * mEff * Vz[i]);
        (*outTau)[i] = static_cast<float>(Tau[i]);
        (*outScl)[i] = static_cast<float>(Scl[i]);
        (*outPh)[i] = static_cast<float>(Ph01[i]);
        (*outWt)[i] = static_cast<float>(wp);
        (*outWtPh)[i] = static_cast<float>(wp * Ph01[i]);
        (*outEsc)[i] = Esc[i] ? 1.0f : 0.0f;
        (*outRad)[i] = static_cast<float>(Rad[i]);
        (*outSpray)[i] = Spray[i] ? 1.0f : 0.0f;
        (*outMF)[i] = Scl[i] < 1.5;
        (*outMC)[i] = Scl[i] >= 1.5;
        maxSpeed = std::max({maxSpeed, std::fabs(Vx[i]), std::fabs(Vy[i]),
                             std::fabs(Vz[i])});
    }

    auto outU = Amino::newMutablePtr<Amino::Array<float>>(nu);
    auto outV = Amino::newMutablePtr<Amino::Array<float>>(nv);
    auto outW = Amino::newMutablePtr<Amino::Array<float>>(nw);
    auto outP = Amino::newMutablePtr<Amino::Array<float>>(ncell);
    for (size_t c = 0; c < nu; ++c) (*outU)[c] = static_cast<float>(u[c]);
    for (size_t c = 0; c < nv; ++c) (*outV)[c] = static_cast<float>(v[c]);
    for (size_t c = 0; c < nw; ++c) (*outW)[c] = static_cast<float>(w[c]);
    for (size_t c = 0; c < ncell; ++c) (*outP)[c] = static_cast<float>(p[c]);

    out_positions = outPos.toImmutable();
    out_velocities = outVel.toImmutable();
    out_mass = outM.toImmutable();
    out_mom_x = outMx.toImmutable();
    out_mom_y = outMy.toImmutable();
    out_mom_z = outMz.toImmutable();
    out_u = outU.toImmutable();
    out_v = outV.toImmutable();
    out_w = outW.toImmutable();
    pressure = outP.toImmutable();
    out_escaped = outEsc.toImmutable();
    radius_out = outRad.toImmutable();
    spray_out = outSpray.toImmutable();
    tau_out = outTau.toImmutable();
    scale_out = outScl.toImmutable();
    mask_fine = outMF.toImmutable();
    mask_coarse = outMC.toImmutable();
    out_phase_state = outPh.toImmutable();
    out_wt = outWt.toImmutable();
    out_wtph = outWtPh.toImmutable();
    out_pos_synced = outSync.toImmutable();
    iterations_used = it;
    final_residual = static_cast<float>(finalResidual);
    max_divergence_after = static_cast<float>(maxDiv);
    max_speed = static_cast<float>(maxSpeed);
    lap(msPart);
    if (prof) {
        std::fprintf(stderr,
                     "PFPROF step=%d grid=%dx%dx%d np=%zu iters=%d asm=%.1f "
                     "solve=%.1f corr=%.1f part=%.1f total=%.1f\n",
                     step_index, NX, NY, NZ, np, it, msAsm, msSolve, msCorr,
                     msPart,
                     std::chrono::duration<double, std::milli>(Clock::now() -
                                                               tBegin).count());
        std::fflush(stderr);
    }
}

} // namespace Solve
} // namespace PFFlip
