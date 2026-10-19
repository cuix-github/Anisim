# Decisions: what was settled, and what measurement reversed

Four decisions were settled at the start. Two were later reversed by measurement, and the positions that moved are recorded here with their dates, because the reasons matter as much as the outcome.

## Settled at the start

- **Uniform-grid pressure solve, adaptive particles.** The paper's grid adaptivity (its Sections 4 and 6) exists to reach ocean scale; a room-scale shot does not need it, and the particle side (Section 5, including the per-size FLIP-blend correction of 5.3) is where the cost explosion actually is, because air fills the whole domain with particles. Consequence accepted: domain size is capped by the uniform solve, a snug box around the action.
- **One custom C++ node, everything else stock.** The two-phase pressure projection does not exist in the product (the shipped liquid and gas projections each assume a single density) and is our one SDK node. Stock nodes were verified for the rest: `splat_points_into_volume`, `sample_property`, `transport_particles`, standard field math and the public iteration machinery the shipped solvers are built from.
- **Public surface only.** Shipped files, public SDK, public docs. The custom node reads and writes plain arrays attached as custom properties, deliberately independent of Bifrost's internal volume layout.
- **Reference first.** A small pure-Python 2D implementation of the core loop (`reference2d/`) is the executable spec; every Bifrost stage is validated against it and against the papers' own test cases.

## Reversed by measurement

- **The transfer moved into a custom node (milestone 6.75, 2026-10-09).** Profiling at finale scale put the stock splat layer at 94 percent of a substep: 91 of 96.4 s. One fused node replaced thirty splats and thirty samples, 9.4x a substep, and brought the pipeline closer to the paper's kernel, not further. The stock graph stays as the control for comparisons.
- **The adaptive pressure solve came back as a staged milestone (6.5a).** Dropped as a non-goal, revived when the solve became the largest remaining slice after the transfer was fused and threaded, built and gated in 2D and 3D, and then kept as an option rather than the default, because at snug-box scale it does not pay. The position moved five times; the dates are below.

## Positions that moved

**Grid adaptivity.** 2026-10-05: no longer a flat non-goal but a staged spike after the uniform 3D build stands. 2026-10-09: spike 1 profiled the pressure share at 0.9 percent of wall, so the adaptive solve is craft, not production. 2026-10-10: the fused transfer shrank the denominator and its sharper kernel raised the solve's cost to 22 percent of a substep, the share the authors' own table reports, so milestone 6.5a was back; with every stage threaded the solve carries 40 percent of a 3.66 s substep. 2026-10-13: stage 3c's cost receipt and a re-read of the paper settled it: the paper calls uniform grids sufficient for several hundred cells per dimension and motivates adaptivity by domains far beyond that, where it is feasibility, not speed; its grid keeps an adaptive cell nearly as cheap as a uniform one, ours is a tree rebuilt every step. Break-even, derived from the receipt: the adaptive solve wins only where the kept share of cells is small, a wide domain that is mostly calm water and open air, and even a free solve speeds a substep at most 1.67x. The decisive test is the hero scene's starting state. The lesson kept from the 0.9-to-22-percent journey: a share is set by the whole pie, so never quote one from measurements taken in different pipeline configurations.

**Junction treatment for the adaptive solve.** Settled 2026-10-13 as the user's choice rather than ours: both junctions stay as gated options (`pressure_solver` 1 naive, the paper's choice, and 2 cleanup) and the uniform solve (0) stays the default. Naive is cheaper and stirs air at rest; cleanup keeps air at rest clean and lands further from the uniform solve in moving liquid. The stirring shows only where something rides the air, and the paper's showcase effects, mist and spray, are exactly that.

**SDK stability posture** (2026-10-04). The node exchanges flat arrays as custom properties and lets stock nodes touch volumes at the edges, which keeps us on the SDK's compact, property-centric core. Deeper native volume traversal exists for extension authors without that core's stability expectations, so if profiling ever demands it, it goes behind a one-file adapter rather than into the solver.

**Colliders** (2026-10-14). The artist's input is the shipped `Simulation::Common::collider` node, chosen over a compound of our own for compatibility with existing liquid setups; the adapter does what the shipped liquid solver's own graph does, through public nodes. Static colliders only, so far; a moving collider needs the field per frame plus its velocity, which the wavemaker of milestone 7 will require.

## Rules of the repository

**Rendering scope, Xue's explicit direction on 2026-10-19.** First reproduce the
PF-FLIP and ST-FLIP demo appearance as closely as the Bifrost implementation and
available cache permit; incremental improvements come after a satisfactory
reproduction. Do not develop the newly proposed hybrid particle classifier or
sphere/ellipsoid replacement scheme as the next step. Retain the valid clustered
foam comparison as experimental evidence, but do not extend its procedural noise
as though it were the paper's method: both papers explicitly exclude procedural
textures. Prioritize conformance of particle-derived surface reconstruction,
feature-preserving smoothing, interface scattering and the paper's spray path.
Arnold-specific numerical approximations must be labelled and validated; the
papers do not supply a complete material/lighting preset. This direction does
did not initially authorize another run. Later on October 19 Xue authorized
offline reconstruction and scattering checks together, still deferring spray
and solver changes. After seeing the valid A/B foam images on desktop, he judged
both promising and requested additional paired angles, distances and later cached
frames before selecting either. That observation study should preserve both
existing looks; it does not turn procedural B into a claim of paper fidelity.

- Every change lands through its own pull request; the repository squash-merges, so a branch whose pull request merged is never pushed onto again.
- The repository stores code and small numeric evidence (JSON receipts, charts, frame sheets). Videos, rendered frames, caches and Maya scenes stay in the ignored `out/` (rule set 2026-10-13; the videos committed before then were untracked the same day).
- Every solver feature ships with its gate: lockstep against the 2D reference, a system gate through the real graph, then physics anchors. Compatibility paths (nz=1, st=0, adapt=0, dump_every=1) stay bit-identical to what the old gates measured.
- Measurement discipline: profile before optimizing, difference the walls of two run lengths so startup cancels, never combine measurements from different pipeline configurations into one ratio, serialize timing runs, keep every arm's artifacts until the conclusion lands.
- Licensed assets (the captured scene) and anything derived from them never enter the repository.
