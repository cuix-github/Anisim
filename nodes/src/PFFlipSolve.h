//-
// Milestone 5: the PF-FLIP step for Bifrost, 2D on a MAC grid.
//
// The graph assembles the staggered channels with stock nodes (two shifted
// splat_points_into_volume passes produce, per face grid, the weighted mean
// particle mass, mean mass-weighted momentum, and mean phase), samples them
// into flat arrays, and this node does the rest of one simulation step:
// normalize to face velocities, add gravity, solve the variable-coefficient
// pressure Poisson equation (Jacobi-preconditioned CG), correct, grid-to-
// particle with the per-phase FLIP blend, and RK2 advection with wall clamps.
//
// The executable spec is reference2d/pfflip2d.py: this file mirrors its
// project(), g2p() and advect() line by line (double precision internally),
// so the two can be compared in lockstep on identical channel inputs.
//
// Flat layouts, row-major matching numpy's ravel():
//   u faces: (nx+1) x ny,  index i*ny + j
//   v faces: nx x (ny+1),  index i*(ny+1) + j
//   pressure: nx x ny,     index i*ny + j
// Particles are float3 with z carried through untouched.
//+

#ifndef PF_FLIP_SOLVE_H
#define PF_FLIP_SOLVE_H

#include "PFFlipSpikeExport.h"

#include <Amino/Core/Array.h>
#include <Amino/Core/Ptr.h>
#include <Amino/Core/RuntimeServices.h>
#include <Amino/Core/String.h>
#include <Amino/Cpp/Annotate.h>
#include <Bifrost/Math/Types.h>

namespace PFFlip {
namespace Solve {

PF_FLIP_SPIKE_DECL
void step_2d(int   nx,
             int   ny,
             float dt,
             float gravity,
             float rho_liquid,
             float rho_air,
             float alpha_liquid,
             float alpha_air,
             int   max_iterations,
             float tolerance,
             const Amino::Array<float>& u_mass,
             const Amino::Array<float>& u_mom,
             const Amino::Array<float>& u_phase,
             const Amino::Array<float>& v_mass,
             const Amino::Array<float>& v_mom,
             const Amino::Array<float>& v_phase,
             const Amino::Array<Bifrost::Math::float3>& positions,
             const Amino::Array<Bifrost::Math::float3>& velocities,
             const Amino::Array<float>& particle_phase,
             Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& out_positions,
             Amino::Ptr<Amino::Array<Bifrost::Math::float3>>& out_velocities,
             Amino::Ptr<Amino::Array<float>>& out_mass,
             Amino::Ptr<Amino::Array<float>>& out_mom_x,
             Amino::Ptr<Amino::Array<float>>& out_mom_y,
             Amino::Ptr<Amino::Array<float>>& out_u,
             Amino::Ptr<Amino::Array<float>>& out_v,
             Amino::Ptr<Amino::Array<float>>& pressure,
             int&   iterations_used,
             float& final_residual,
             float& max_divergence_after,
             float& max_speed)
    AMINO_ANNOTATE("Amino::Node");

// Milestone 6: the same step in 3D on a uniform MAC grid, with a geometric
// multigrid-preconditioned CG solve (preconditioner 0 = Jacobi, 1 = multigrid
// V-cycle). At nz = 1 the w faces are all walls and the step reduces exactly
// to the 2D method, which is the lockstep bridge used by validate_3d.py.
// pressure_solver 1 and 2 replace the uniform solve with the adaptive one
// (milestone 6.5a, PFFlipAdaptive.h): a graded tree fine within
// adaptive_band cells of the liquid interface (0 = 5) and up to
// adaptive_levels coarser elsewhere (0 = 2), its own multigrid PCG, and an
// exact local fill of the uniform faces inside coarse cells. 1 uses the
// PF-FLIP paper's naive junction stencil (one solve); 2 the cleanup junction
// (a second-order variational solve, then a naive-stencil projection of its
// leftover), and iterations_used then counts both solves. 0 (the default) is
// the uniform solve, untouched; the preconditioner port then applies, and the
// uniform solve also stands in whenever obstacles reach the domain.
// Solids: obstacles_min/max are axis-aligned boxes. solid_phi is any solid as
// its signed distance at the grid nodes, in cells, negative inside (a mesh
// collider, through mesh_to_level_set and one sample at the nodes);
// solid_mode 0 ignores it, 1 closes a face whose centre the solid covers (the
// boxes' rule), 2 weights each face by its open area (cut cells, Batty,
// Bertails and Bridson 2007). Particles inside are pushed out along the
// field's gradient (PFFlipSolids.h). Both kinds keep the uniform solve.
// With solid_mode on, a field that is not one float per grid node (the wrong
// grid, or a .npy of the wrong dtype, which read_NumPy delivers empty) is an
// error: the node reports it through the runtime services (the graph's own
// error channel, and stderr under bifcmd), returns empty arrays and sets
// iterations_used to -1, rather than simulating without the solid.
// Flat layouts, row-major like numpy ravel() on (nx[+1], ny[+1], nz[+1]):
//   u: (nx+1)*ny*nz   index (i*ny + j)*nz + k
//   v: nx*(ny+1)*nz   index (i*(ny+1) + j)*nz + k
//   w: nx*ny*(nz+1)   index (i*ny + j)*(nz+1) + k
//   p: nx*ny*nz       index (i*ny + j)*nz + k
//   solid_phi: (nx+1)*(ny+1)*(nz+1)   index (i*(ny+1) + j)*(nz+1) + k
// A piston wavemaker (2026-10-16): with piston_period > 0 and piston_width > 0 a
// plate stands on the floor, the whole depth in z, from x to x + piston_width up
// to y = piston_top (cells), its x moving as piston_x0 + S/2 (1 - cos 2 pi t / T)
// with S = piston_stroke (cells) and T = piston_period (seconds), t = step_index
// * dt. The plate is a cut-cell solid folded with solid_phi (if any) and moving
// at the step's mean speed: on the faces it closes, in the divergence of the cut
// faces it covers, in the push-out (the normal velocity relative to the plate's).
// The reference's Sim(piston=...) is the spec (validate_piston.py) and
// gates/validate_piston_node.py the lockstep. Period 0, every port's default,
// leaves the step as it was. Seed nothing inside the plate (as the reference's seed
// drops what lies inside a solid): a particle left there is pushed out of its front
// face into the water on the first step.
// Droplet radii (2026-10-16, PF-FLIP Section 7): with droplet_model 1 (and
// escape on) a liquid particle the escape flag marks is a droplet carrying its
// own radius in metres (radius_in / radius_out, 0 = none): drawn once from a
// log-normal of median droplet_r_median and ln-spread droplet_r_sigma by a hash
// of (step_index, index), the same number the reference draws; shrunk every step
// by the d^2 law of evaporation, r^2 falling by droplet_evap_k / 4 per second;
// deleted below droplet_r_min (the particle leaves, every output shorter by it).
// Its drag is then the sphere's own: the relaxation rate (3/8)(rho_air /
// rho_liquid)(C_d / r)|v_rel|, Schiller-Naumann C_d(Re) below Re 1000 and 0.44
// above, Re from air_nu (air's kinematic viscosity) and metres_per_cell, applied
// implicitly; drag_droplet is then unused. Rejoining the water resets the radius
// to 0. Bubbles are untouched. droplet_model 0, the default, leaves the step as
// it was (radius_out = 0 for every particle). The reference's Sim(droplet=...)
// is the spec and gates/validate_droplet_node.py the lockstep.
// Extra spray emission (2026-10-17, PF-FLIP Section 7 after Ihmsen et al. 2012):
// with emit_model 1 (and the droplet model and escape on) a liquid particle at the
// interface (its liquid fraction under 1 - esc_phi, the escape test's own line seen
// from the water) in a cell whose velocity gradient (the Frobenius norm of the
// Jacobian of the cell-centred velocity, 1/s) exceeds emit_grad_threshold emits,
// with probability min(1, emit_rate (G - G0) dt) from a hash of (step_index,
// index), a spray particle at its own position and velocity, at most
// emit_max_per_cell a cell a step in index order: liquid, massless (mass 0 in
// out_mass, counted in no splat, so the flow never feels it), flagged in spray_in /
// spray_out (1 = spray), always a droplet (a radius drawn at its first step, the
// drag, the evaporation), deleted when it falls into the water (its liquid
// fraction above 1 - esc_phi) or evaporates, and kept behind every other particle
// so the others keep their order. The paper names the method and no thresholds.
// spray_cfl > 0 gives every droplet its own integrator (the paper's RK2 with
// locally adaptive time stepping): sub-steps of at most spray_cfl cells, the drag
// relaxation toward the air velocity (sampled at the sub-step's start) taken
// exactly with gravity at the rate of the sub-step's midpoint slip, the position by
// the trapezoid of the sub-step's two velocities; 0, the default, keeps the one
// implicit step. The reference's
// Sim(emit=..., spray_cfl=...) is the spec and gates/validate_spray_node.py the
// lockstep.
// A wave zone (2026-10-18): the wavemaker of a numerical wave basin without a
// paddle, the relaxation zone of ship and coastal CFD (Jacobsen et al. 2012). With
// zone_length > 0 (cells) and waves holding at least one train, after every
// projection the face velocities within zone_length of the x = 0 wall are blended
// toward the velocity field of a sum of linear wave trains, u = (1 - w) u + w u_t,
// and so are the liquid particles' velocities after their FLIP update. waves holds
// four floats a train: its height (cells), period (seconds), direction (degrees from
// +x toward +z) and phase (radians); zone_depth is the still water depth (cells)
// the linear theory is evaluated at and zone_ramp the seconds over which every
// height rises from 0 (a smoothstep; 0 = full from the first step). The weight w is
// the zone's exp weight, 1 - (exp(s^3.5) - 1) / (e - 1) with s = x / zone_length,
// times the face's liquid fraction (the water forced where it is, the air left
// alone; a fade over the cell above the trains' summed surface, the first zone's,
// pumped the water up one way, see docs/history.md). Each train: k from
// w^2 = g k tanh(k h) (Newton), the horizontal velocity a w cosh(k y) / sinh(k h)
// cos(psi) along its direction, the vertical a w sinh(k y) / sinh(k h) sin(psi), psi = k (x cos + z sin) - w t + phase, y
// clamped at the still depth, t the end of the step, (step_index + 1) dt. (The
// trains' Stokes drift is not subtracted: tried as a closed basin's return flow, it
// cut the radiated height to 0.68 of the train's; the particles the zone carries
// drift forward by it, an open question.) The wall
// faces the zone touches are its inlet: the x = 0 faces, and in 3D the side-wall
// faces along the zone, carry the trains' velocity (times the liquid fraction) less its
// flux-weighted mean over every inlet face, so no net flux is forced into the
// sealed box (a plain inlet has no pressure solution there; returning each column's
// flux through the air faces above it, tried, cut the radiated height to 0.67); a
// particle carried into a wall is held at it, as everywhere, and the escape test
// leaves the zone's liquid alone, so the wall's own column, which breathes by the
// mean and piles at the wall, never spawns droplets. Below
// zone_length 0, every port's default, the step is untouched bit for bit. The
// reference's Sim(zone=...) is the spec (angles are 0 at nz = 1) and
// gates/validate_wave_zone.py the lockstep, the train against linear theory and the
// oblique train's angle.
PF_FLIP_SPIKE_DECL
void step_3d(int   nx,
             int   ny,
             int   nz,
             float dt,
             float gravity,
             float rho_liquid,
             float rho_air,
             float alpha_liquid,
             float alpha_air,
             int   max_iterations,
             float tolerance,
             int   preconditioner,
             int   pressure_solver,
             int   adaptive_band,
             int   adaptive_levels,
             int   escape,
             float esc_phi,
             float drag_droplet,
             float drag_bubble,
             float buoyancy,
             float rho0_face,
             int phase_model,
             int transport_model,
             float volume_tau,
             int   st,
             int   st_seed,
             int   step_index,
             int   adapt,
             float coarse_gain,
             float ws_epsilon,
             const Amino::Array<float>& scale_in,
             const Amino::Array<Bifrost::Math::float3>& obstacles_min,
             const Amino::Array<Bifrost::Math::float3>& obstacles_max,
             int   solid_mode,
             const Amino::Array<float>& solid_phi,
             float piston_stroke,
             float piston_period,
             float piston_x0,
             float piston_width,
             float piston_top,
             int   droplet_model,
             float droplet_r_median,
             float droplet_r_sigma,
             float droplet_evap_k,
             float droplet_r_min,
             float metres_per_cell,
             float air_nu,
             const Amino::Array<float>& radius_in,
             int   emit_model,
             float emit_grad_threshold,
             float emit_rate,
             int   emit_max_per_cell,
             int   emit_cloud_size, float emit_event_radius, float emit_spread,
             float emit_velocity_spread, float emit_radius_max, int emit_max_active,
             int emit_seed, float spray_dt_max,
             float spray_cfl,
             const Amino::Array<float>& spray_in,
             float zone_length,
             float zone_depth,
             float zone_ramp,
             const Amino::Array<float>& waves,
             int external_forces,
             const Amino::Array<float>& liquid_acceleration_u,
             const Amino::Array<float>& liquid_acceleration_v,
             const Amino::Array<float>& liquid_acceleration_w,
             const Amino::Array<float>& linear_drag_u,
             const Amino::Array<float>& linear_drag_v,
             const Amino::Array<float>& linear_drag_w,
             const Amino::Array<float>& tau_in,
             const Amino::Array<float>& u_wt,
             const Amino::Array<float>& v_wt,
             const Amino::Array<float>& w_wt,
             const Amino::Array<float>& u_ws,
             const Amino::Array<float>& v_ws,
             const Amino::Array<float>& w_ws,
             const Amino::Array<float>& u2_mass,
             const Amino::Array<float>& u2_mom,
             const Amino::Array<float>& u2_phase,
             const Amino::Array<float>& u2_wt,
             const Amino::Array<float>& u2_ws,
             const Amino::Array<float>& v2_mass,
             const Amino::Array<float>& v2_mom,
             const Amino::Array<float>& v2_phase,
             const Amino::Array<float>& v2_wt,
             const Amino::Array<float>& v2_ws,
             const Amino::Array<float>& w2_mass,
             const Amino::Array<float>& w2_mom,
             const Amino::Array<float>& w2_phase,
             const Amino::Array<float>& w2_wt,
             const Amino::Array<float>& w2_ws,
             const Amino::Array<float>& u_mass,
             const Amino::Array<float>& u_mom,
             const Amino::Array<float>& u_phase,
             const Amino::Array<float>& v_mass,
             const Amino::Array<float>& v_mom,
             const Amino::Array<float>& v_phase,
             const Amino::Array<float>& w_mass,
             const Amino::Array<float>& w_mom,
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
             int&   iterations_used,
             float& final_residual,
             float& max_divergence_after,
             float& max_speed,
             const Amino::RuntimeServices& services)
    AMINO_ANNOTATE("Amino::Node");

// Fused P2G (milestone 7 prep): one pass over the particles scatters every
// channel the step node consumes - mass, per-axis momentum, phase, and the
// ST temporal weight - onto all three MAC face grids and both size tiers,
// directly into flat arrays, replacing the per-channel stock splat + sample
// round-trips. The math is the reference's (pfflip2d.py _splat_group): Eq. 6
// kernel w = max(1 - d^2/r^2, 0)^3 on a (2 ceil r + 1)^d stencil, r = 1 fine
// and r = 2 coarse, face frames base = floor(p' + 0.5) with p' = pos minus
// the face offset. Channel semantics stay what the stock path produced, so
// step_3d is untouched: weighted mean with eps_mean added to the
// denominator, and ws = K/(K + ws_epsilon) so the step node's eps-inversion
// recovers the exact kernel sum K. st == 0 leaves the wt channels empty
// (plain-phi path); adapt == 0 leaves ws and every tier-2 channel empty and
// splats all particles into tier 1, matching the old graphs exactly.
PF_FLIP_SPIKE_DECL
void p2g_3d(int   nx,
            int   ny,
            int   nz,
            int   adapt,
            int   st,
            float eps_mean,
            float ws_epsilon,
            int phase_model, float phase_alpha, float rho0_face,
            float rho_liquid, float rho_air,
            const Amino::Array<Bifrost::Math::float3>& positions,
            const Amino::Array<float>& scale,
            const Amino::Array<float>& mass,
            const Amino::Array<float>& momx,
            const Amino::Array<float>& momy,
            const Amino::Array<float>& momz,
            const Amino::Array<float>& wt,
            const Amino::Array<float>& phase,
            Amino::Ptr<Amino::Array<float>>& u_mass,
            Amino::Ptr<Amino::Array<float>>& u_mom,
            Amino::Ptr<Amino::Array<float>>& u_phase,
            Amino::Ptr<Amino::Array<float>>& u_wt,
            Amino::Ptr<Amino::Array<float>>& u_ws,
            Amino::Ptr<Amino::Array<float>>& v_mass,
            Amino::Ptr<Amino::Array<float>>& v_mom,
            Amino::Ptr<Amino::Array<float>>& v_phase,
            Amino::Ptr<Amino::Array<float>>& v_wt,
            Amino::Ptr<Amino::Array<float>>& v_ws,
            Amino::Ptr<Amino::Array<float>>& w_mass,
            Amino::Ptr<Amino::Array<float>>& w_mom,
            Amino::Ptr<Amino::Array<float>>& w_phase,
            Amino::Ptr<Amino::Array<float>>& w_wt,
            Amino::Ptr<Amino::Array<float>>& w_ws,
            Amino::Ptr<Amino::Array<float>>& u2_mass,
            Amino::Ptr<Amino::Array<float>>& u2_mom,
            Amino::Ptr<Amino::Array<float>>& u2_phase,
            Amino::Ptr<Amino::Array<float>>& u2_wt,
            Amino::Ptr<Amino::Array<float>>& u2_ws,
            Amino::Ptr<Amino::Array<float>>& v2_mass,
            Amino::Ptr<Amino::Array<float>>& v2_mom,
            Amino::Ptr<Amino::Array<float>>& v2_phase,
            Amino::Ptr<Amino::Array<float>>& v2_wt,
            Amino::Ptr<Amino::Array<float>>& v2_ws,
            Amino::Ptr<Amino::Array<float>>& w2_mass,
            Amino::Ptr<Amino::Array<float>>& w2_mom,
            Amino::Ptr<Amino::Array<float>>& w2_phase,
            Amino::Ptr<Amino::Array<float>>& w2_wt,
            Amino::Ptr<Amino::Array<float>>& w2_ws,
            const Amino::RuntimeServices& services)
    AMINO_ANNOTATE("Amino::Node");

// Dump policy: pass the write path through only every Nth substep (step 0
// always dumps; every <= 1 means every substep). An empty path makes the
// downstream write node skip the dump, which keeps per-substep file I/O out
// of production runs without touching the write nodes themselves.
PF_FLIP_SPIKE_DECL
void dump_gate(int step_index,
               int every,
               const Amino::String& path,
               Amino::String& path_out)
    AMINO_ANNOTATE("Amino::Node");

} // namespace Solve
} // namespace PFFlip

#endif // PF_FLIP_SOLVE_H
