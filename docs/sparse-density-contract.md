# Sparse-support density candidate: fixed experiment

2026-10-24. This is one new candidate after the positive-count log law failed on
retained impact states. No production defaults, phase law, P2G kernel, pressure
tolerance, escape thresholds, emission parameters or particle positions are
changed to make the candidate pass. Long impact/holdout and 3D remain deferred.

## Formula and scope, frozen before execution

Losasso et al. 2008 Eq. (5) discretizes the density term using the target density
in the denominator: `(n_target - n_current)/(n_target * time_interval)`.
The following paragraph replaces this with a log-averaged rate; the previous
candidate used that log form. This candidate instead retains the finite
difference, with the existing relaxation time tau=.5 seconds rather than the
simulation substep. This is an explicit PF implementation choice, not a claim
to reproduce unpublished author code or to use Losasso's entire diffuse-SPH
method. Source: https://jerrytalton.net/research/ltkf-tcspl-08/paper.pdf, §3,
Eq. (5) and following paragraph. Sparse primary air is not ballistic spray.

Let q_i=n_i/r_i, where r_i>0 is the same calibrated wall-aware reference used
by the old candidate. The raw target velocity divergence is

    s_i = (q_i - 1) / tau

In the closed uniform grid, require zero total divergence and choose the
least-squares correction to s:

    d_i = s_i - mean(s) = (q_i - mean(q)) / tau.

The actual pressure equation is unchanged: A=-D beta G,
`A p=(d-D u*)/dt`, `u_new=u*-dt beta G p`; hence `D u_new=d`.
The global mean projection is valid only for the connected closed domain in
scope. It does not replace per-component compatibility for future obstacles.

Zero local n is admissible and yields finite negative divergence. There is no
epsilon floor or logarithm of zero. The interpretation is a finite request to
redistribute carrier samples into an undersampled region; this does not create
particles, fill velocity faces directly, or guarantee the region will refill.
Negative/nonfinite counts, nonpositive reference, all-zero carrier support,
nonpositive tau and nonfinite resulting targets are refused.

Near q=1, q-1 and log(q) agree to first order. They differ for large deviations;
in particular this candidate is not invariant to uniformly scaling all counts
or references. Fixed r retains the initial four-samples-per-cell calibration.
Uniform q still gives d=0 after compatibility projection, so this law also cannot
repair global loss. It is a sampling regularizer, not a volume-conservation proof.

For constant r and a continuum n satisfying continuity with closed walls, the
domain mean qbar is constant. Substituting this d gives a useful independent
check: the entropy-like quantity integral(q log q) has derivative
`-integral((q-qbar)^2)/tau <= 0`. This does not prove the discrete particle
implementation stable: the wall reference is nonconstant, kernel sums are not
partition-of-unity material measures, and conversion/deletion add sources.
Nonzero d can inject kinetic energy; the existing energy and geometric checks
remain mandatory. No stability theorem for the whole solver is claimed.

Only active primary carrier samples contribute, exactly as in the current
detached-primary contract. Airborne primaries and auxiliary spray do not.
Conversion/deletion affects next-step support; primary represented water stays
in its separate material inventory, and reentry restores its previous weight.
The new `volume_model=1` is opt-in; 0 preserves the old logarithmic calculation
and its refusal behavior. `volume_tau=0` keeps density feedback off. The new
mode remains restricted to the existing uniform-2D phase=2 scope, with the same
no-primary-evaporation restriction when detached primaries are used.

## Finite tests and stop conditions

Before rebuilding, retain actual-graph default-off and supported log controls
and the original DLL. After rebuilding, reproduce all saved state arrays bit
for bit with model=0. Test the exact density vector `[0,.5,1,1.5,2]` with unit
reference and tau=.5, obtaining `[-2,-1,0,1,2]`; check compatibility, invalid
inputs, and near-reference first-order agreement separately.

Run the new model on the three frozen impact inputs at steps 1280/3200/7200,
and compare actual Bifrost against independently computed NumPy with identical
float32 ports: position <8e-6 cells, velocity <5e-5 cells/s, requested pressure
relative residual <=1e-6, finite outputs and unchanged primary water inventory.
The old model must still refuse its zero-count input. A separate prescribed-grid
single step checks actual divergence against the finite target, independent of
the particle integrator and of whether two implementations share a formula.

Only if local checks pass: the existing one-second density perturbation must
reduce centred log-density RMS by at least 20%; the existing six-second 48×40
rest case must retain RMS<=.1 cells/s, COM drift<=.1 cells and sampled energy
ratio<=1.005. Continuous versus split restart must match every state array.
Use a generic displaced sampling patch to ensure support acceptance is not
specific to the old impact demo. Retain rejected inputs and every failed arm.

This round does not redefine the geometric area gate: coarse-circle calibration
5.118% still exceeds 5%, while refined values are retained evidence. No absolute
volume certificate or long-suite admission follows solely from these tests.
The area error budget must be resolved explicitly before the integrated suite.

Stop downstream execution at a failed physical gate. At most two bounded
diagnostics may clarify a failure; no second candidate, epsilon tuning, changed
thresholds or automatic long run is introduced in this experiment.
