# Wave features and the October 20 continuation

On October 19, 2026, Xue authorized a separate continuation of the existing
wave basin from 18 to 30 seconds with gradient spray emission and the droplet
model enabled. The original simulation remains the control and is never
overwritten. He also requested a small Bifrost-driven 2D particle demonstration,
then asked for larger waves and 60 seconds after seeing that the initial
30-second demonstration produced very little additional spray.

## What the wave zone does

The existing implementation blends velocities toward a prescribed wave field
inside the first 24 metres of the basin. It acts on the grid after projection
and on liquid particles after their FLIP update, with spatial blending and an
eight-second ramp. It is an external wavemaker, based on the relaxation-zone
approach already documented in [history](history.md). It is not an initial
displacement alone, and it does not prescribe the later breaking shape.
Outside the zone, propagation, interaction with the beach, breaking and
particle escape are computed by the solver.

The larger 2D diagnostic doubles the component heights, including the two lead
waves from 3/1.5 metres to 6/3 metres, without changing their periods. These
component heights are not measurements of the resulting breaking crest.
All headings become unidirectional, with phases matched at x=100 m, z=40 m.
This is a diagnostic wave flume, not a slice of the 3D solution or a calibrated
coastal engineering prediction.

## Particle playback

The actual `sim_3d_fused` Bifrost compound runs our custom transfer and solver
nodes with nz=1. ST, particle adaptivity, escape, per-droplet radii, the spray
integrator and gradient emission are enabled. Adaptive pressure remains off:
the existing adaptive pressure path does not support the solid beach.

The original 344×48×1 diagnostic completed 30 seconds and 720 exact frames in
115.85 seconds including cache export. It reached 147 escaped liquid particles
and only two additional spray particles. Doubling wave heights alone reached
512 escaped particles and six spray particles in the first 30 seconds.
Increasing only the 2D emission rate from 0.05 to 1.0, while retaining the
10 s⁻¹ gradient threshold and eight-particle cell cap, made the additional
spray visible. This diagnostic rate does not change the scheduled 3D settings.

The first stronger 60-second run exposed another limitation: spray approached
the original 24-metre ceiling and the short dry beach ended at a nearby wall.
The final diagnostic uses 260×40 metres at 0.5-metre cells, or 520×80×1, while
retaining the original water depth, beach toe/slope and forcing. The extra
space is air above and dry beach after the original shoreline. All runs remain
closed finite domains; moving the boundaries is not an open-ocean boundary
condition.

The final local run is
`out/runs/wave2d_big_open_60s_20261020/`. Its `open_in_maya.bat` opens
`scenes/wave_2d.ma`. The timeline has 1,440 frames at 24 fps, sampled directly
every two 1/48-second solver steps. No index-based interpolation is used:
spray births and deaths change particle counts. Separate cached objects hold
water, escaped liquid, emitted spray and hidden air. Marker sizes are for
inspection and are not physical droplet diameters. In a wireframe viewport,
switch to shaded display to inspect the configured point colors.

The scene also contains the fully configured solver compound named
`wave_solver_OPEN_IN_BIFROST`. Its outputs are disconnected so ordinary
playback only reads cache. Manually evaluating that compound can start a new
simulation into `manual_solver_output`; it is not a real-time solver.
The scene builder checks metre/centimetre conversion against actual Maya
world-space bounds at representative frames for all nonempty liquid layers.

`experiments/wave_2d_audit.py` checks the saved particle positions and measures
both populations outside the forcing zone and away from the outer boundaries.
It creates a diagnostic point plot, not an Arnold image. Per-run receipts
record the measured counts and chosen inspection frame.

The final 60-second simulation and cache export took 294.17 seconds; the
twelve Bifrost launches account for 213.66 seconds of that total, including
their process startup. It reached 745 escaped liquid particles and 134
emitted spray particles. At frame 939 (39.104 s), the interior region contains
363 escaped particles and 133 spray particles simultaneously. There are 571
sampled frames meeting the inspection criterion of at least 200 escaped and
50 emitted particles in that region. Neither population enters the top or
far-wall two-metre bands in any saved frame. This checks those boundary
interactions, not absence of all finite-domain effects.

Maya standalone evaluated six representative frames and verified all nonempty
liquid-layer bounds. After the first evaluation (81 ms), those sampled
evaluations took 14–27 ms including cache reads and the bounds checks.
These are not measured GUI frame rates; the replacement's GUI playback was
not tested after the Computer Use interruption. The scene opens at frame 939.

Xue explicitly requested removal of the original 30-second 2D run.
Automatic approval blocked the recursive removal before execution, reporting
only “blocked by policy”; `wave2d_spray_20261020` therefore remains in place.
The stronger pilot runs are retained as experimental evidence. Computer Use
was stopped by a physical Escape event; the final scene is prepared through
Maya standalone, without further UI actions or a claim of GUI playback QA
for this replacement.

## Scheduled 3D continuation

`out/runs/wavesea4_spray_18_30_20261020/inputs/plan.json` identifies the
unchanged 18-second source, `wavesea4/dumps/c2`, containing 61,086,383 particles.
All thirteen full-state channels were checked for shape, finite values,
spray mass/radius consistency and SHA-256 identity.
The immutable runtime snapshot records commit `a6a7cc2`; later 2D changes do
not alter it. The DLL and graph hashes match the preflight validation.

The 3D plan retains 688×96×320 cells and dt=0.025, emission threshold 10 s⁻¹,
rate 0.05 and the original waves. It enables the existing droplet and spray
paths only from 18 seconds onward; it cannot recreate spray missing from the
earlier history. Twenty-four half-second segments retain full checkpoints and
exact raw dumps every 0.125 seconds. Old data and failed attempts are retained.
No rendering is chained to the simulation.

An app heartbeat named **Anisim 18–30秒喷雾续跑** starts the supervisor at
07:00 Toronto on October 20 and checks hourly between 07:00 and 23:00.
The supervisor checks progress/resources every 15 seconds without model calls.
Normal progress stays quiet. Failure triggers bounded investigation and at
most two repair attempts per segment, without disabling features, relaxing
pressure convergence or lowering resolution. A process lock and orphan
check prevent concurrent workers from this frozen runtime. The local
supervisor can continue after 23:00; the next model inspection is at 07:00.
Completion or an unresolved failure pauses the heartbeat.

The previous 18–30-second interval took about 7.8 hours; additional spray
overhead is unmeasured. Retaining all new inputs, checkpoints and raw dumps is
estimated at roughly 500 GiB, with 250 GiB disk reserve. The supervisor
checks commit memory and halts only its own worker tree on exhausted reserves
or thirty minutes without a new completed step. This estimate is not a fixed
storage cap. The computer and Codex app must remain running.

Preflight reproduced stock/fused resume controls bit for bit, including
per-step positions, phase, scale and diagnostics; resetting the step offset
failed the negative control as expected. Spray off/lockstep/birth-death tests
passed, as did an all-feature beach/zone restart comparison across all thirteen
final arrays. File gates reject partial/nonfinite states and modified runtime
files, and preserve newborn spray with a zero radius until its first
integration step. An early launch was rejected, and no production segment
was started while preparing the job. Local evidence is under
`out/runs/wave_preflight_20261020/`.
