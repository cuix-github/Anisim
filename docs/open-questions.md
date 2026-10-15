# Open questions

Kept current; the answered ones stay at the bottom with their answer and date, because the answers are findings.

## Open

- **A CFL-bound substep.** Every long run steps at a fixed dt. The 512 run's air reached 900 to 1450 cells per second and crossed 76 to 123 cells in one substep, past anything the push-out reaches, which is how a handful of air particles ended up inside a boulder for one frame; the 2D long runs' splashes made particles at 120 to 460 cells per second in every arm. Fine for an A/B between solvers, not for a quality claim.
- **Moving colliders.** Solids are static: one signed-distance field per run, velocity zero. The piston wavemaker of the breaking-wave ocean needs a solid that moves every substep and carries its velocity into the pressure solve and the push-out. The plan: the 2D reference first, gated against linear wave theory, then the node.
- **Rendering droplets and bubbles.** How escaped particles are drawn, decided now that the shot exists: in progress on the `surface-preview` branch as a meshed surface plus a fog volume for the spray.
- **A fourth push-out pass.** The push-out stops once a pass lands outside, so a settled particle sits at or just outside the surface, and a third pass that lands a rounding error short gets no fourth; cheap to add, not a leak.
- **The reference's refinement criterion wraps around the domain** (`band_target` uses `np.roll`), so a band touching one wall also refines the opposite wall; harmless, wasteful. The node mirrors it for lockstep; retire both together and re-record the stage-3a gate.
- **The refinement criterion reads part-filled cells as interface**, which assumes every deep cell keeps its four seeded particles; once motion shuffles them it halves the tree's savings (44 to 22 percent in the resting tank under naive junctions). Candidates: a smoothed fraction or the grid's own phase field, changed in the reference and the node together with the wrap above.
- **The 2D reference's resting tank is not silent.** Its liquid's rms speed grows from 0.10 to 0.54 cells per second over 200 frames under the plain uniform projection. The Bifrost loop does not show it (the same tank's rms decays from 0.013 to 0.003 over 100 s), so it is a trait of the Python stepping, not of the method; worth a look before any rest test leans on the reference.
- **Volume-drift correction** (the paper's divergence correction term): add if a run visibly loses volume; none has so far.
- **Pressure tolerance.** The authors solve to 1e-4 where we run 1e-6; a gated lever on the solve's 40 percent of a substep.

## Answered

- ~~Why does the stock-graph ST plunge produce 56 peak bubbles against the reference's 259?~~ Answered 2026-10-11 (Aniso PR #48): not a solver bug. The validator forwarded only `rho0_face` while the graph's escape parameters defaulted to zero, so the two sides had run different physics since the test was born; with parameters matched the plunge scores droplets 99 against 94 and bubbles 237 against 259. The story is under "Both gaps closed" in the [history](history.md).
- ~~The 20-second bore has zero actual-state obstacle penetrations but 1,225 display-position penetrations.~~ Answered 2026-10-11 (Aniso PRs #43 and #44): the time resynchronization ran after the final obstacle sweep; the synced display position now gets the same push-out, position only, and the re-rendered demo audits at zero.
- ~~Hero scene.~~ Answered 2026-10-13: a tsunami into a coastline corner of a licensed captured scene (milestone 7).
- ~~Preconditioner for the 3D solve: incomplete Cholesky or a multigrid V-cycle.~~ Answered 2026-10-06 (Aniso PR #25): a geometric multigrid V-cycle inside flexible PCG solves the dam-break system in 9 iterations instead of 533, same solution; incomplete Cholesky was never needed.
- ~~Is the stock splat's kernel smooth enough for a quiet phase field?~~ Answered at milestone 4: the stock linear kernel's field matches the reference, and swapping the paper's Eq. 6 kernel for a tent kernel moves the interface midpoint by under 0.01 voxels.
- ~~Junction treatment for the adaptive solve.~~ Settled 2026-10-13: both kept as options, uniform the default; see [decisions](decisions.md).
- ~~Arbitrary solids.~~ Answered 2026-10-14 by the mesh-collider stages 1 to 3 and the showcases; see the [history](history.md).
