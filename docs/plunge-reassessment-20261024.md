# Reassessing the cumulative 2D plunge candidate

2026-10-24. **The cumulative candidate is not accepted for physical fidelity or
coarse 3D.** Its lower terminal speed does not establish that it reproduces the
papers. This read-only reassessment identifies a changed spray lifecycle and a
defect in using globally centered, all-phase sampling feedback as liquid-volume
control after deleting air samples. It does not establish a single cause for
every visually sticky feature. No solver edit, build, simulation, or render was
performed for this reassessment; the original artifacts remain evidence.

The frozen experiment and numerical checks are in
[current-plunge-visual.md](current-plunge-visual.md) and
[current_plunge_20261024.json](../results/current_plunge_20261024.json).
The local raw jobs are `out/runs/current_plunge_20261024/` and
`out/runs/settling_phase_20261023_r2/`. Both use the same initial arrays,
64 x 80 x 1 grid, approximately 0.0025 s timestep, 1 m/cell, and water/air
FLIP fractions 0.97/0.90. Several other choices changed together, so their
cumulative visual comparison is not a one-variable causal test.

## Shared zero evaporation changed the spray experiment

The current job sets `droplet_evap_k=0`; the previous long job used `1e-9`.
The new detached-primary guard rejects any nonzero coefficient, although the
coefficient is shared with auxiliary spray. Thus preserving numerical primary
material also disabled auxiliary evaporation. The code still removes auxiliary
droplets on water entry or at the radius cutoff, but it has no age timeout.
This was an experiment-design regression relative to the intended auxiliary
lifecycle, not evidence that longer-lived white points validate the paper.
See the guard, shared radius update, and deletion logic in
[PFFlipSolve3.cpp](../nodes/src/PFFlipSolve3.cpp).

The current ledger closes as **21,737 births = 21,383 returns + 0 evaporations +
354 survivors**. Last birth occurred at 30.172499 s. Read-only inspection of
the 240 s snapshot found survivor radius min/median/max of
6.82/26.08/42.93 micrometres and vertical velocity median of -0.0952 m/s.
None was within 0.01 cell of the allowed sidewall or ceiling clamp. These are
small, mostly descending survivors, not persistent late emission. The white
two-pixel display markers do not represent physical droplet diameter.
The ledger is in the [committed result](../results/current_plunge_20261024.json);
radius and velocity measurements use `chunk_11/final_radius.npy`,
`final_velocities.npy`, `final.npy`, and `final_spray.npy` in the current run.

PF-FLIP section 7 includes evaporation and section 3.4 describes returning
escaped primary droplets to FLIP. Neither promises that every small droplet
falls back within a fixed duration. Ihmsen's foam-age rule is also not a
universal spray timeout. These mechanisms must not be conflated.
([PF-FLIP](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf#page=15),
[Ihmsen et al. 2012, section 3.2](https://cg.informatik.uni-freiburg.de/publications/2012_CGI_sprayFoamBubbles.pdf#page=5))

For the implementation's law `r_new^2 = max(r_old^2 - k*dt/4, 0)`, applying
the old coefficient only to the saved survivor radii gives cutoff times of
0.18/2.72/7.37 s (min/median/max). This is a scalar counterfactual, not a rerun
or a prediction that positions would follow the same trajectory. Restoring
auxiliary evaporation requires independent lifecycle control; removing the
primary guard blindly would reintroduce a separate material-loss risk.

## Global count compatibility does not conserve each phase's volume

The finite candidate deposits all active primary air and water samples against
a fixed initial reference. It excludes auxiliary spray and detached primary
water. With `q_i = n_i / reference_i`, the implemented target is

```text
s_i = (q_i - 1) / tau
d_i = s_i - mean(s) = (q_i - mean(q)) / tau
```

Consequently, **any spatially constant q is a zero-correction state**; the
operator does not require q=1. Subtracting a global mean makes the closed-domain
pressure target compatible, but it does not preserve the volume of each phase.
This follows from [targetFinite](../nodes/src/PFFlipVolume.h) and its assembly
in [PFFlipSolve3.cpp](../nodes/src/PFFlipSolve3.cpp), independent of pressure
tolerance or floating-point differences.

A two-region counterexample fixes the problem precisely. Let water occupy
fraction f of the domain with correct sampling q_water=1. Remove some air
samples, leaving q_air=1-delta. Then q_mean=1-(1-f)*delta and
`d_water=(1-f)*delta/tau > 0`. The correction asks already-correct water to
expand while air contracts. The whole-domain divergence integral can vanish
while the liquid-region integral does not. For an advected liquid region,
`dV_water/dt = integral_water div(u) dV`; whole-domain compatibility is
insufficient. This is a general defect in the candidate's intended use with
changing phase sample inventories, not a scene-specific tuning issue.

The current snapshot supports that mechanism. Initial primary counts were
6,008 water and 14,472 air. Final counts are 6,008 water and 12,781 air:
**1,691 air samples were deleted**, or 8.26% of the initial all-phase inventory.
Twelve primary waters are detached, leaving 18,777 active carriers. Independent
read-only reconstruction of the implemented density/reference field gives:

| Final count-density statistic | Value |
| --- | ---: |
| Whole-domain q mean | 0.92068261 |
| Whole-domain q standard deviation | 0.00148610 |
| Interior deep-water q mean | 0.92073595 |
| Interior deep-air q mean | 0.92075383 |
| Centered target RMS at tau=0.5, 1/s | 0.00297220 |

The two phases have nearly equal sampling density about 8% below their initial
reference. This is the equilibrium requested by the operator, not recovery to
q=1. The measurement deposits current active particles using
`max(1-dx^2-dy^2,0)^3` and the exact initial wall-clamped four-sample reference
from [PFFlipVolume.h](../nodes/src/PFFlipVolume.h). Interior slices are
x=3:61, y=3:17 for water and x=3:61, y=35:75 for air. Positions and reference
coordinates are first stored as float32, then evaluated in float64 for this
independent diagnostic. No solver or integration routine is called.

Final water-height P99 is 25.2576 cells and COM_y is 12.8271, versus
21.5423 and 10.9598 in the prior run. These are **particle-distribution
statistics, not reconstructed geometric volumes**. The prescribed initial
water area is `64*18 + 14*25 = 1502 cells^2`, with ideal flat depth
23.46875 cells before boundary and reconstruction calibration. A crude common
density estimate `80*6008/18789 = 25.58 cells` is consistent with the higher
layer; it is not a measured volume error or proof of every trajectory's cause.
The older, lower waterline is not a geometric ground truth either.
([Scene mask](../experiments/spray_events_2d.py),
[retained measurements](../results/current_plunge_20261024.json))

PF-FLIP section 3.5 motivates a particle-density correction but does not specify
this complete discretization. Losasso's density treatment distinguishes the
physical meaning of the sampled population and warns about inappropriate
fixed-density forcing in diffuse material. The papers do not establish this
finite, globally centered operator as a per-phase volume guarantee under
changing sample inventories. With sample deletion, the
continuity accounting needs its source term:
`partial_t n + div(n*u) = S_n`. Replacing finite feedback with a centered
logarithm or a small epsilon does not resolve the lost distinction between
sample inventory and material volume.
([PF-FLIP section 3.5](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf#page=8),
[Losasso et al. 2008, sections 3-4](https://jerrytalton.net/research/ltkf-tcspl-08/paper.pdf#page=3))

Auxiliary spray already transfers zero P2G mass and disappears on return
without adding a primary water particle. Its birth/return ledger therefore
does not explain increased primary-water count; that count is unchanged.
The concern above is main-fluid sampling feedback after air deletion.

## Stickiness remains a separate causal question

No new material-viscosity parameter or explicit extra damping was selected.
The velocity FLIP/PIC blend is unchanged; transport mode 3 changes particle
position advection through water-side extension and local RK3. Extension
stencils, detached-primary transfer, density feedback, and repaired air
deletion all differ from the older run. Their effects cannot be separated by
watching the cumulative movie alone.
([Transport implementation](../nodes/src/PFFlipTransport.h),
[velocity and position updates](../nodes/src/PFFlipSolve3.cpp))

PF-FLIP equation 13 estimates mixing-related numerical viscosity as
`nu=(1-alpha_FLIP)*h^2/(6*dt)`. The diagnostic's h=1 m, dt=0.0025 s and
alpha=0.97/0.90 formally give about **2/6.67 m^2/s** for water/air.
These are scale estimates under the paper's approximation, not measured
constitutive viscosities. They show that no *new* damping does not imply
little numerical dissipation. Because these parameters match the older run,
they cannot alone explain the new visual regression.
([PF-FLIP equation 13](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf#page=8))

ST-FLIP measures remaining volume loss and dissipation in particular tests;
that is not permission to accept this candidate's phase-volume error.
Its Appendix B describes render-time surface reconstruction and protected
curvature smoothing, which were not used in this particle video and cannot
repair particle trajectories. A smaller late RMS can mean that a wrong target
state was reached efficiently. The existing stationarity test also still
fails its surface-height condition.
([ST-FLIP sections 4.3-4.4 and Appendix B](https://ge.in.tum.de/download/ST-FLIP.pdf#page=12),
[original experiment result](current-plunge-visual.md))

## Keep verified fixes separate from candidate models

Matching the phase formula, fixing float32 momentum-port comparison, preserving
complete restart state, and preventing unintended primary deletion under the
declared no-primary-evaporation contract remain supported corrections. Opt-in
bubble deletion has its own fixture. These do not validate the combined
density, transfer, and transport model. Zero-P2G airborne primary transfer,
hysteresis details, extension stencils, and the finite density target remain
documented implementation choices with different evidence requirements.
([Phase ladder](phase-repair-ladder.md),
[interface checks](independent-interface-closure.md),
[density candidate contract](sparse-density-contract.md))

The current multi-droplet emitter is also a candidate, not a full reproduction
of every cited diffuse-material mechanism. It uses gradient-excess event rate,
16 children per event, and disk/ball offsets about the source endpoint.
Ihmsen equation 8 instead includes a kinetic-energy potential with wave-crest
and trapped-air potentials, and distributes births through a swept cylinder.
The current expected-r-cubed event budget and clipping policy are explicit
choices. More particles must not hide a defective main-fluid state.
([Emitter contract](spray-events-2d.md),
[implementation](../nodes/src/PFFlipSprayEvents.h),
[Ihmsen equation 8](https://cg.informatik.uni-freiburg.de/publications/2012_CGI_sprayFoamBubbles.pdf#page=4))

## Proposed acceptance route stays finite

**No following experiment has been started.** First retain the current result
as the failed integrated control and calibrate the geometric-area measure
against known shapes; the existing coarse-circle absolute-area failure must
not be waived. Freeze a phase/material ledger and numeric acceptance bounds
before another run. Arithmetic validity of the density unit tests remains
useful, but integrated liquid-volume acceptance is withheld.

The first fluid comparison should replay the same short impact, roughly
20-40 s, with only density feedback disabled (`volume_tau=0`). Hold code,
initial state, seeds, dt, mixing, transfer, transport and evaporation fixed.
Measure calibrated liquid area, phase counts/deletions, morphology and named
energy measures. The predeclared question is whether removing this feedback
removes the air-deletion-driven deformation. Only if sticky motion remains
should a second arm change transport 3 to 0, with every other setting fixed.
Do not simultaneously tune emission or add surface smoothing.

Separately, an explicit primary/auxiliary evaporation split needs a small
fixture checking the radius-squared law, drag and water return against actual
Bifrost C++ and independent NumPy. Retain the primary-material witness and
verify auxiliary on/off leaves the one-way carrier unchanged. Only after these
two issues close should one fixed integrated 2D replay test area, inventories,
lifecycles and visible recovery together. An emitter redesign, coarse 3D,
adaptivity, or rendering is outside this acceptance step.

The conclusion required for promotion is a coherent model and measured phase
behavior, not merely matching Python/C++ arithmetic, pressure convergence,
particle-count preservation, a lower waterline, or a quieter final frame.
