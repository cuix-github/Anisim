#pragma once

#include <cstddef>
#include <vector>

// Milestone 6.5a: the adaptive pressure projection behind step_3d's
// pressure_solver = 1 and 2, a port of reference2d/adaptive_projection.py
// with its two junction treatments. Dimension-generic: at nz = 1 the tree is
// a quadtree in x-y whose cells stay one fine cell thick, and the method is
// exactly the 2D reference's, which gates/validate_adaptive_step.py checks in
// lockstep. Internal to the pack: not parsed for nodes.

namespace PFFlip {
namespace Adaptive {

using Vec = std::vector<double>;

// The reference's junction argument.
//   Naive    the conservative two-cell junction stencil the PF-FLIP paper
//            adopts: one solve, first order at junctions (stage 3b)
//   Cleanup  the second-order variational solve, then a naive-stencil
//            projection of what its junctions leave (stage 3b2; in 3D,
//            stage 3c: the junction sample corrects along both tangents).
enum class Junction { Naive, Cleanup };

struct Settings {
    int band = 5;              // fine cells kept fine around the liquid interface
    int levels = 2;            // coarsest leaf level: cells up to 2^levels wide
    double tolerance = 1e-6;   // relative residual of the flexible PCG
    int maxIterations = 2000;
};

struct Result {
    int iterations = 0;        // the main solve's
    int cleanupIterations = 0; // Cleanup's second solve (0 for Naive)
    double residual = 0.0;     // the main solve's relative residual at exit
    size_t unknowns = 0;       // adaptive cells
    int levels = 0;            // leaf levels actually used (clamped to fit the grid)
    Vec pressure;              // per uniform cell: the covering adaptive cell's value
                               // (Cleanup: the two solves' pressures summed)
};

// u, v, w: the uniform MAC face velocities after gravity and walls (in/out),
// in the step node's layouts; bu, bv, bw: face coefficients 1/rho, zero on
// walls; liquidFraction: per uniform cell, liquid particles over the seeding
// count per cell (4 in 2D, 8 in 3D), the reference's refinement input.
Result project(int nx, int ny, int nz, double dt, const Settings& settings, Junction junction,
               const Vec& liquidFraction, Vec& u, Vec& v, Vec& w,
               const Vec& bu, const Vec& bv, const Vec& bw);

} // namespace Adaptive
} // namespace PFFlip
