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
