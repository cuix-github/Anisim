# Independent PF interface closure

2026-10-24. Xue authorized deriving and testing an explicit implementation rather
than waiting for unpublished author details. This supplements the bounded route;
it does not erase its failed geometry calibration or 240-second near-rest gate.
No production simulation, rendering, data deletion or automatic PR merge is part
of this work.

## Candidate frozen before its experiments

`droplet_transfer_model=1` is an opt-in **one-way** approximation, restricted to
uniform two-dimensional `phase_model=2`, `escape=2`, with the droplet integrator
enabled and `droplet_evap_k=0`. Model 0 retains the existing implementation.
Primary evaporation is deliberately excluded from this conservative parcel
experiment; its physical-radius loss cannot silently erase a larger represented
parcel. Existing evaporation law tests remain separate evidence.

The previous radius is the persistent airborne flag. Ordinary primary water
enters the droplet state at phi < .3 and remains there until phi > .7 (using the
existing esc_phi port for the lower bound). The hysteresis is our explicit
choice. While airborne, a primary particle exports zero transfer mass and
momentum, and contributes neither phase nor kernel support to the next P2G.
Its phase, scale, position, velocity and radius survive. On reentry its original
sampling weight is reconstructed from phase/scale; it does not create a second
particle. The existing FLIP/PIC update applies on reentry. This is not a claim
of exact two-way drag or impact-momentum conservation.

Classification follows the current pressure solve. Thus newly escaping primary
particles leave P2G on the **next** step. Restart state must include all existing
arrays and transfer weights consistent with the radius flag. Old caches are
not silently converted. Auxiliary emitted spray remains zero-transfer and is
removed on return; it never adds primary water volume.

For an initial uniform two-dimensional four-per-cell sample, each primary has
book area h²/4 (or volume h²/4 times a declared slab thickness). This is a
separate material inventory, not the integral of reconstructed phi and not
4*pi*r³/3. Splitting or variable weights are outside this candidate. With
evaporation off, primary water inventory must remain exactly constant; deleting
classified air inclusions changes the air sampling inventory only.

The density term is explicitly a **relative sampling regularizer**, counting
only active continuum primaries. For positive n and reference n0,
`d = (log(n/n0) - mean(log(n/n0))) / tau`. The fixed reference and tau=.5 s are
retained; the new option removes only the fixed-population restriction. Empty
or nonfinite support is still refused, never floored to conceal a hole.
This permits a declared redistribution response to localized particle removal;
it does not prove liquid geometric volume conservation. A global multiplier of
n0 cancels exactly and therefore cannot repair local conversion errors.

On equal-volume closed cells, subtracting the mean is the unique minimum-L2
change satisfying zero total requested divergence. This proves pressure
compatibility, not zero energy injection. In source-free continuum transport
with a constant reference, log-density differences decay exponentially with
time constant tau. It does not prove monotonic Eulerian RMS near our walls or
after particle sources/sinks. Those require the finite experiments below.

PF sections 3.4/7 do not explicitly specify the airborne P2G mask. Losasso's
ballistic/slip construction includes P2G, so our earlier implication that
"purely Lagrangian" necessarily requires zero P2G was too strong. The reason
to test this candidate is isolation of the droplet from its own sampled carrier
velocity, not an assertion about unpublished author code.

## Fixed finite gates

1. Reproduce retained baseline state arrays before and after the change; model 0
   must remain bit-identical. Check formula scaling invariance, zero sum,
   projection optimality, and positive-support refusal independently.
2. Prescribed constant phase .2/.5/.8: test entry, memory through the mixed band,
   reentry, restored transfer mass, and unchanged represented water inventory.
   Positions/radius flags and exact zero/1000 transfer weights have exact oracles.
3. Stationary-air primary flight: test the actual graph against Python and an
   independent drag ODE at successively smaller droplet steps. The new mode
   must leave carrier air exactly at rest. The old mode is retained as a
   coupling contrast, not judged against the new mode's one-way oracle.
4. Local deletion/return and three existing common impact states: Python versus
   graph max position <8e-6 cells, velocity <5e-5 cells/s; no loss of primary
   water with evaporation off. Continuous and checkpointed graph output must
   be bit-identical across all 13 state arrays. Nonfinite/empty support refuses.
5. Combined six-second 48x40 rest: existing RMS<=.1, COM displacement<=.1,
   energy peak ratio<=1.005 and relative reconstructed-area drift<=.02 remain.
   Area is explicitly the old reconstruction measurement, not absolute volume.
6. Only after these gates pass, run the fixed complete impact recovery and an
   independent offset falling-water geometry. Preserve the old near-rest gate
   (10-second window, RMS<=.05, p99<=.2, kinetic fraction<=1e-4, surface range
   <=.25). No long independent trajectory pointwise-match gate is introduced.

The suite records each failure and does not launch downstream expensive cases
after a prerequisite fails. One additional targeted diagnostic is permitted per
failed interface; it cannot adjust a threshold after seeing the result.

Sources: [PF-FLIP](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf),
[Losasso et al. 2008](https://jerrytalton.net/research/ltkf-tcspl-08/paper.pdf).
Detailed independent notes and derivations remain in this run's research folder.

## First result and bounded diagnostic

The first local run at `54c5d6c` passes the thirteen-array legacy control,
algebra, prescribed lifecycle, stationary-air flight refinement and exact restart.
It refuses the first common impact input because sampling density is zero in
some cells. The combined candidate has therefore **failed** its support gate.
Before considering a different model, one diagnostic counts holes before/after
primary exclusion at all three retained inputs, then checks the same candidate
with only `volume_tau=0` at those inputs. The latter isolates lifecycle/transport;
passing it cannot qualify the integrated density candidate. Empty cells are not
floored, reseeded or given a new reference in this round. The long recovery and
holdout are withheld because the prerequisite failed.

The diagnostic finds 2/20/13 pre-existing empty cells at the three impact inputs,
and 2/20/14 after detached-primary exclusion. All are classified air. Thus the
first refusal is not caused by the new mask; the fixed positive-support law
already has an unsupported-domain problem. All three one-step common-state
checks pass with density feedback off. No zero-support policy is added here.

A separate lifecycle witness exposes a concrete material-loss bug in the new
candidate: both a newly drawn radius and a carried radius below `r_min` delete
a primary even when evaporation is zero. The correction clamps only these
new-mode primary aerodynamic radii to `r_min`; the represented parcel remains.
It is an explicit cutoff choice, not an inference of physical droplet volume.
Auxiliary droplets and model 0 retain their deletion behavior. Finite positive
radius parameters and finite nonnegative carried radii are required. The same
prescribed fixture is retained before/after, including the old-mode controls.

After this correction, the local gates are repeated to qualify the changed
binary. The six-second static case may run as an independent supported-domain
control, not as admission past the failed impact-support gate. It cannot enable
the long impact or holdout. All original numerical limits remain fixed.

## Executed outcome: integrated candidate refused, no production promotion

The corrected implementation is `0d0f9fc`; DLL SHA256 is
`fb1eacf574587b391187887150a717506bf5d7f5fc0c6ed92a3d33efd61dc07b`.
The compact [evidence receipt](../results/interface_closure_20261024.json)
includes source/binary provenance, thresholds, successful subgates, both radius
witnesses, and the failing integrated run. All raw arrays and logs remain under
`out/runs/interface_closure_20261024/`.

| Question | Actual result | Scope of conclusion |
| --- | --- | --- |
| Does default mode change? | Thirteen retained final-state arrays are bit-identical; nine old transport cases, analytic RK3/extension and three configuration refusals pass. | The checked old paths remain unchanged; no new production default. |
| Does one-way flight implement its specified drag? | At .05 s, 3 cells/s becomes 1.32362437 at max droplet step .001 s; independent RK4 gives 1.32361601. Halving .002/.001/.0005 s reduces relative errors to 2.469e-5/6.321e-6/1.818e-6. Python/C++ states coincide in these new-mode arms; carrier air stays exactly still. | Verifies the prescribed stationary-air one-way limit. |
| Is the old coupling observably different? | Same radius, initial state and .001 s droplet step give 2.68378377 cells/s with full transfer. The primary remains airborne at radius 1e-4 throughout both arms. | Full transfer changes its own sampled carrier. This is not proof all two-way coupling is wrong or that this explains the old long movie. |
| Is conversion/reentry internally consistent? | Phase .2/.5/.8 fixtures pass; positive radius persists in the mixed band; reentry restores transfer weight. 80 continuous steps equal 40+40 restart across every state array. An inconsistent old checkpoint refuses. | Supports the declared persistent state and transfer-weight contract. |
| Can a primary disappear with zero evaporation? | Before correction, both carried and newly drawn sub-cutoff primary radii delete the water parcel. After correction both survive at `r_min`; old-mode controls still delete it. | Closes this concrete primary inventory defect only in the new mode. |
| Do common impact states match? | With density feedback OFF, all three actual graph/Python cases pass: max position 7.6294e-6 cells, max velocity 2.4796e-5 cells/s. Each retains its input primary-water count. | Common-input execution agreement; these fixtures construct candidate transfer weights/momentum from the same state, not an exact continuation of old checkpoints. |
| Does the combined density candidate work on these states? | FAIL: first state refuses zero support. Counts before/after exclusion are 2/2, 20/20, 13/14; all holes classify as air. | The fixed log law lacks an admissible empty-support model. This is an interface/domain failure, not a divergent pressure iteration. |
| Is the short supported-domain static control stable? | Both tau=0 and .5 pass the six-second fixed limits in both backends. C++ final water RMS is .006722/.009992, COM shift .02143/.01253 cells, energy ratio 1/1.0000985, area drift .0782%/.1077%; all sampled water counts remain 2304. | Static supported states pass. This does not override the impact failure or certify absolute geometric volume. |

The local run took 56.45 s, support diagnostic 12.66 s, inventory-after witness
16.17 s, and static control 139.74 s including graph calls and independent Python.
These are separate wall-clock receipts, not per-solver performance benchmarks.

The prescribed local deletion/return density-response fixtures, new-mode
emission-on integrated test, complete impact recovery and independent holdout
have **not** been executed. Downstream admission stops at the failed support
prerequisite. The old 240-second near-rest failure and coarse geometric
calibration failure remain unchanged. No 3D, rendering, production-cache change
or new author correspondence occurred.

The next research question is now finite: how to define density feedback in
under-sampled air without pretending that zero samples mean zero physical air,
while keeping the closed-domain divergence constraint and a bounded energy
response. Liquid-only feedback, a supported-domain correction, or conservative
sampling maintenance are different models, not interchangeable tolerance fixes.
Select and derive one before a new preregistered test. An epsilon floor, broader
kernel or refreshed reference added just to make this case run would invalidate
the current candidate's interpretation. An author reply is not a prerequisite.
