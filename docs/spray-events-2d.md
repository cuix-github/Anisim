# Spray event implementation and 2D validation

2026-10-23. This stage covers the emission specification, Python/C++ mechanism
checks, and a particle-only 2D comparison. It does not authorize coarse 3D,
production simulation, rendering, or deletion of the existing 40-second cache.

## Scope and sources

PF-FLIP Section 7 describes gradient-triggered spray, a lognormal droplet size
distribution, and drag/evaporation with the simulated air as carrier. Ihmsen
et al. separate source eligibility from an emission rate integrated over time,
and distribute births in space. These public ideas motivate the model below;
its exact coefficients and random sampler are independent implementation
choices, not a reproduction of the authors' unreleased solver.

- https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf
- https://cg.informatik.uni-freiburg.de/publications/2012_CGI_sprayFoamBubbles.pdf
- https://people.maths.bris.ac.uk/~majge/rpp8_3_036601.pdf

The original escaped FLIP particles retain their existing classification,
numerical mass and motion. This revision adds clouds to **additional** spray
emission. Splitting an escaped particle's represented liquid volume is a
separate unresolved contract; it is not silently replaced by massless copies.
No private email or private source code is included in this public specification.

## Opt-in event model

`emit_model=0` disables emission; `1` preserves the existing single-birth path;
`2` selects this model. Source eligibility and gradient G (in inverse seconds)
remain the existing interface test. A source's expected event count in a step
is lambda = emit_rate * max(G - emit_grad_threshold, 0) * dt. Its requested
event count is floor(lambda) plus a Bernoulli draw with probability frac(lambda).
Unlike a probability clamped to one, this has the correct expected count even
when lambda exceeds one. It is stochastic rounding, **not a Poisson process**;
variance and individual trajectories need not be timestep invariant.

Each event requests N=`emit_cloud_size` independent droplets. An artist can
change event rate (more represented spray water per second) independently of N
(finer subdivision at the same expected event water volume). The chosen radius
law is

    r_i = R * N^(-1/3) * exp(sigma * Z_i - 1.5 * sigma^2), Z_i ~ Normal(0,1).

R=`emit_event_radius` is a spherical-volume-equivalent radius in metres;
sigma=`droplet_r_sigma` is dimensionless. Before clipping,
E[sum(r_i^3)] = R^3 because E[exp(3 sigma Z)] = exp(4.5 sigma^2).
This is expected volume, not exact per-event conservation. Radii are clipped
to [droplet_r_min, emit_radius_max]; the ledger reports both unclipped and
clipped sums of r cubed. Rejected births lose their represented volume rather
than increasing the weights/radii of survivors. The main water loses no mass
in this auxiliary model; it is not a conservative atomization solver.

Births have a uniform disk (2D) or ball (3D) offset of at most `emit_spread`
cells about the source's end-of-step position. Velocity offsets are independently
sampled in a disk/ball with radius `emit_velocity_spread * |source velocity|`.
They have zero mean in expectation; no backwards direction or artificial wind
is prescribed. Positions outside the domain or inside a solid are rejected,
not clamped into piles. These are explicit bounded sampling choices rather
than claims about the authors' exact emission kernel. Newborns start at the
end of the step with a positive radius and begin flight/evaporation next step.

Random samples are keyed by seed, main-step index, source index, event ordinal,
child ordinal and channel. C++ and NumPy use the same integer hash; normal
transforms can differ by floating-point roundoff. These are birth keys, not
persistent particle IDs across adaptive reordering. Validation compares common
states strictly; independent full runs use phase/trajectory distributions.

`emit_max_per_cell` remains a per-step resource cap, not a physical emission
rate. `emit_max_active` limits live auxiliary droplets; the ledger records
requests, accepted births and each rejection cause. Source traversal is ordered
and deterministic, so saturated budgets can bias which sources survive. A
visual preset is accepted only with this saturation visible and quantified.

## Flight and dimensional consistency

Positions and velocities use cells and cells/second internally. Physical radii,
air viscosity and evaporation use SI; `metres_per_cell` converts relative speed.
Both implementations retain the same Schiller–Naumann sphere drag and diameter
squared evaporation law. `droplet_model=1` preserves the old integrator.
The opt-in `droplet_model=2` samples carrier velocity at a predicted midpoint
and limits flight substeps by both CFL and `spray_dt_max`. Stiff drag uses
exponential relaxation, not an unstable explicit RK2 kick. Convergence in a
spatially varying carrier is measured separately from uniform-wind relaxation.
This does not imply second-order accuracy for the complete split flow solver.

2D droplets still use spherical 3D drag/radius physics, moving in a planar
carrier. Their visual surface density is not a quantitative prediction for 3D.

## Acceptance sequence

First reproduce the original off and single-birth lockstep controls. Then check
event-count expectation above and below one event per step, random moments,
radius bounds and expected volume, budget accounting, birth spread, and finite
state/error handling. Compare the new path through the actual Bifrost graph
against NumPy on identical float32-rounded inputs. Include steady-air drag,
evaporation and a nonuniform-carrier step-refinement test.

The final independent 2D runs use identical seeding, domain, physical constants,
step schedule and source settings, with both phases enabled. Preserve a common
baseline and change one feature at a time. Record actual runtime configurations,
code and binary hashes, birth/death statistics, bulk motion and spray geometry.
Show visible particles and the full impact/recovery sequence, not just an
attractive still. A long-run disagreement is investigated, not excused as
chaos without a common-state check. Known Python/Bifrost long-run differences
in the repository remain relevant; this test must not conceal them.
