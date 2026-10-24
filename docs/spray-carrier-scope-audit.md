# Spray carrier feedback and density support are separate contracts

2026-10-24. This bounded follow-up separates two issues previously combined in
the interface-closure admission gate. It uses the unchanged solver binary from
that experiment, retained impact inputs, and three tiny actual Bifrost graph
fixtures. It changes neither production defaults nor the C++ solver. No long
simulation, render, production-cache change, or deletion was performed.

The executable audit is [audit_spray_carrier_scope.py](../gates/audit_spray_carrier_scope.py).
The full numerical [receipt](../results/spray_carrier_scope_20261024.json)
includes the source commit, binary hash, all hole coordinates and tolerances
inherited from the existing interface test. Raw results and the failed harness
attempt remain under `out/runs/spray_carrier_scope_20261024/`.

## What is being tested

PF-FLIP section 7 describes droplets driven by gravity and air-velocity drag,
integrated as point trajectories, and separately describes higher-resolution
rasterization for rendering. This does not require a second spray pressure
solve. The word Lagrangian alone does not specify feedback; the one-way tracker
is the explicit contract tested here. Full water-air coupling in the carrier
is distinct from adding droplet feedback to that carrier.

The current auxiliary spray already has zero transfer mass/momentum and is
skipped by fused P2G, including its support channels. Its nonzero aerodynamic
radius still controls inertia/drag: zero *transfer weight* is not a claim that
the physical droplet is massless. Legacy escaped primary water retains its
original grid weight. Experimental `droplet_transfer_model=1` excludes airborne
primaries too; it remains default-off and limited to uniform 2D.

Sources: [PF-FLIP sections 3.4, 3.5 and 7](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf)
and [Losasso et al. 2008 sections 3-5](https://jerrytalton.net/research/ltkf-tcspl-08/paper.pdf).
Private correspondence and private code excerpts are not reproduced here.

## Same droplet, three transfer paths

All arms start from the existing stationary-air fixture: 16 x 12 cells, one
water droplet at 3 cells/s, radius 0.0001 m, zero gravity/evaporation, and the
same drag integrator. Twenty steps advance 0.05 seconds. Density correction
and new emission are disabled so this tests flight/feedback only.

The legacy arm first reproduces every retained final array bit for bit. The
other arms alter only the auxiliary flag or the primary-transfer mode. C++
stably places auxiliary particles last even with emission disabled; Python
does not reorder in this fixture. The known input-identity permutation is
inverted before comparing arrays, without spatial matching or fitting.

| Path | Air velocity sampled at the drop immediately after initial P2G, before projection (cells/s) | Final drop speed (cells/s) | Maximum final air-particle speed (cells/s) |
| --- | ---: | ---: | ---: |
| Legacy primary, still transfers to grid | 2.744910 | 2.683784 | 1.707543 |
| Auxiliary spray, zero transfer | 0 | 1.323624 | 0 |
| Experimental detached primary, zero transfer | 0 | 1.323624 | 0 |

The auxiliary and detached arms have all 13 compared physical arrays exactly
equal; only their spray flags differ. Both also match the NumPy implementation
exactly in position and velocity. The legacy C++/NumPy differences remain below
the existing limits (position 9.54e-7 cells, velocity 2.38e-7 cells/s). The earlier
independent drag ODE/refinement evidence remains in the interface-closure report.

This establishes a concrete mechanism: a full-weight airborne primary changes
the carrier it subsequently samples, while auxiliary spray does not. It is not
evidence that every two-way physical model is wrong, or that this single mechanism
explains the full old 20-second movie. It does not validate conservative breakup,
reentry impulses, evaporation of primary represented volume, or 3D.

## Empty counts are local sampling geometry, not a spray pressure solve

Three retained impact inputs are inspected without advancing them. Their zero
cell-centre counts exactly reproduce the previous support diagnostic. A separate
direct distance/kernel evaluation verifies the splat results; it does not reuse
the splat index stencil. Every zero has its nearest active particle at least one
cell away, outside the compact kernel support. No cancellation or underflow is
needed to explain the zero.

| Saved input step | Empty cell centres | In outermost cell ring | Within outer two-cell band | Centres with any of four faces still supported | Airborne particles touching an unsupported interior P2G velocity face at their current sampling position |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1280 (about 3.2 s) | 2 | 2 | 2 | 2 | 114 / 1715 |
| 3200 (about 8 s) | 20 | 15 | 20 | 9 | 0 / 3133 |
| 7200 (about 18 s) | 14 | 14 | 14 | 10 | 84 / 14286 |

These are 36 cell-time observations, not 36 distinct cells over the full run.
All classify as air. Thirty-five observations already existed before excluding
airborne primaries; only one additional zero is introduced by that exclusion.
All are within two cells of a wall, and 21 retain some adjacent face support.
Nearest active-particle distances span 1.00145 to 2.03419 cells.

Cell-centre count and staggered velocity-face support are different quantities.
A zero at a cell centre is not proof that every velocity sample nearby is absent.
At the 8-second input none of the 3133 current airborne sampling stencils touches
an unsupported interior face, although the global positive-count density formula
still refuses that input. This directly separates that refusal from an immediate
spray interpolation failure at the inspected positions.

This does NOT prove that boundary support is harmless: the other two inputs do
have exposed sampling stencils; pressure projection and boundaries may provide
finite velocities on unsupported P2G faces; and pressure is globally coupled.
This audit does not measure the accuracy of those post-projection velocities or
attribute the earlier motion that created the gaps. It identifies the present
geometric reason for zeros and narrows further investigation to wall sampling.

## Acceptance scope after this audit

Keep two explicit outcomes. The existing one-way auxiliary flight contract and
the opt-in detached-primary flight contract pass this controlled carrier test.
The centred logarithmic density candidate still fails on empty air sampling
support. Its failure remains a base-fluid result; it must not be relabelled as
evidence that spray requires density restoration, P2G feedback or another pressure
solver. Conversely, flight passing cannot certify the carrier's volume behavior.

Prescribed-carrier flight, emission and lifetime tests can proceed independently
of that density candidate. A full coupled long-impact/holdout certificate still
requires a selected base-flow configuration and the previously defined physical
checks. None is waived here, and no coarse 3D or production run is admitted.
The next base-flow question is specifically a justified boundary/empty-support
contract, not arbitrary spray feedback or an unbounded search over all cells.
No epsilon floor, reseeding rule, alternative kernel or relaxed test threshold
has been introduced to turn the failed candidate into a pass.

## Execution record

The first launch found SciPy absent from the pinned environment before any run;
the diagnostic was changed to direct NumPy distances without adding dependencies.
The first graph attempt reproduced the legacy control, then exposed the auxiliary
ordering difference in the comparison harness. That failed attempt is retained.
After applying the known identity permutation in the harness, the complete audit
passed in 12.64 seconds including graph calls and Python comparison. No solver
source or binary changed during these corrections.
