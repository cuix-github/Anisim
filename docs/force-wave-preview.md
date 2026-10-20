# Force-driven wave previews, October 20, 2026

A one-hour budget covers building and inspecting a two-dimensional preview
and a coarse three-dimensional preview before choosing a production resolution
or start time. These are new runs from still water. The previous 18–30 s
continuation was cancelled; old caches are retained. Nothing in this study
launches a production run or an Arnold render.

The new source applies acceleration before the pressure projection. It does
not blend toward target velocities or reset particle velocities after the
projection. `zone_length=0` disables the previous wave zone. The generic solver
interface accepts liquid acceleration on the three MAC face grids, in cells/s²,
and nonnegative linear damping rates in 1/s. The measured liquid fraction
weights the acceleration; damping acts on both phases. Each component follows
`u_new = (u_old + dt * liquid_fraction * acceleration) / (1 + dt * damping)`.
Pressure projection, boundary conditions, transfer and advection then proceed
through the existing solver.

The scene recipe stays outside the C++ kernel. `experiments/force_wave_preview.py`
bakes acceleration fields without reading the water state. The Bifrost graph
reads them at each absolute substep, including across checkpoints. This is
an offline spectral source feeding the actual custom solver graph, not yet a
live editable spectrum generator compound inside Maya. The saved Maya scene
includes that solver graph, disconnected from the playback outputs so playing
the timeline only reads caches.

## The matched preview recipe

| Property | Two-dimensional preview | Coarse three-dimensional preview |
| --- | --- | --- |
| Domain | 180 × 32 × 1 m | 180 × 32 × 72 m |
| Grid | 180 × 32 × 1 | 180 × 32 × 72 |
| Cell size / timestep | 1 m / 1/24 s | Same |
| Simulated duration | 30 s | Same |
| Still-water depth | 12 m | Same |
| Beach | Toe at x=70 m, slope 1:5, shoreline x=130 m | Same |
| Source | Gaussian envelopes around x=38 m, sampled over x=14–62 m, ramp over 4 s | Same |
| Spectrum | 3 lead waves + 997 background components, seed 20 | Same amplitudes and periods |
| Lead components | Target heights 4, 3, 1.5 m; periods 8, 9.5, 6.5 s | Same, headings +15°, −18°, +5° |
| Background | Target significant height 1.5 m; periods 3–7.5 s | Same, headings clipped to ±35° |
| Dimensional reduction | Headings set to zero; lead phases refocused | Directional phases retained |
| Playback | 360 actual saved frames at 12 fps, no interpolation | Same |

The height targets specify the source's linear transfer estimate, not guaranteed
nonlinear crest heights. Three lead phases are chosen to meet near x=100 m,
z=36 m at t=18 s in that estimate. The two-dimensional run preserves component
amplitudes, not a literal central slice of the three-dimensional interference.
The 997-component background is an explicit experimental distribution, not a
claimed reproduction of a named ocean spectrum from the papers.

Both previews enable ST-FLIP, escaped liquid particles, adaptive air particles,
the droplet model and gradient-triggered spray emission. They use the same
10/s emission threshold, 0.05 emission rate and cap of eight new particles per
eligible cell; there is no diagnostic emission boost. The solid beach requires
the existing uniform pressure path. Adaptive pressure is not active.

The source emits waves both ways. A back damping region at x<20 m and, in 3D,
six-metre side damping strips reduce reflected motion, but these are finite
closed-domain boundaries. Their reflection performance has not been validated.
The 2D preview already shows rear-wall pile-up; force-driven generation does
not by itself remove boundary problems.

## Validation before the mixed sea

The original no-wave-zone control ran on the old DLL before modification.
After rebuilding, all compared grid and particle arrays matched that control
bit for bit. The new disabled-field path also matches the ordinary path.

An independent NumPy acceleration/damping specification, followed by the
existing reference projection and particle update, matches the real Bifrost
node: maximum grid error 4.27e-7, position error 3.82e-6 and velocity error
3.93e-7. The pressure residual is 4.92e-7. A nonfinite enabled field is rejected.
These are small two-dimensional checks of the generic mechanism, not a
convergence study of the ocean scene.

A separate 16 × 12 × 8 single-phase anchor applies a discrete curl field
with zero divergence and zero normal wall velocity. Projection must preserve
`dt * acceleration`, including its varying z component. Maximum velocity
error is 4.66e-10. This tests three-dimensional field indexing independently
of the ocean recipe.

The checkpoint gate compares a continuous 16-step run with two eight-step
runs using the same absolute force clock. All 13 saved state channels match
bit for bit. Resetting the force clock in the second segment is a negative
control and changes velocity by 0.00533 cells/s.

A single 1 m, 6 s target wave on a flat 260 m basin calibrates the source
gain. Gain 1 produced approximately 1.47 m in the first inspection; gain 0.7
was selected and run for 35 s. A harmonic fit over 20–34.96 s gives heights
1.127 m at x=75 m and 1.007 m at x=95 m. Their phase difference gives a 52.65 m
wavelength, versus 50.69 m from finite-depth linear dispersion: 3.85% error.
The probes use the 97.5th particle-height percentile in a 5 m band, not a
reconstructed surface. The original and calibrated runs remain separate.

## Local artifacts and inspection

Each run is under `out/runs/force_sea_2d_20261020/` or
`out/runs/force_sea_3d_20261020/`. `inputs/plan.json` records the complete recipe,
`inputs/forces/` stores acceleration fields, `dumps/c*/` stores actual solver
outputs and checkpoint states, and `receipts/` stores measured timings and
particle checks. `particle_diagnostic.png` shows particle-only side views.

`maya/build_wave_2d_scene.py` builds either dimensionality. The saved
`scenes/wave_2d.ma` or `scenes/wave_3d.ma` and `open_in_maya.bat` load the
four playback layers. Intended colors are blue for bulk water, cyan for escaped
water and orange for additional spray; grey air is hidden by default. These
colors and full GUI playback have not been accepted. Bulk-water display is limited
to 120,000 points per frame, and air to 10,000, solely for viewport speed.
Escaped water and emitted spray are not thinned. Original dumps retain all
particles. Viewport point sizes are visualization aids, not physical radii.

The Maya builder checks world-space bounds against saved metre-scale points
at six representative frames. It records evaluation times separately from
simulation timings; those headless times are not a measured GUI playback rate.
The 3D cameras use perspective. The source band, beach and still water are
useful scale references, not a rendered water surface.

Both 30 s runs completed. The 2D run took 42.7 s including field baking and
playback export (39.8 s across six graph calls). The 3D run took 594.5 s
including baking/export (579.9 s across six graph calls), using 32 OpenMP
threads. Graph-call times include Bifrost startup and disk output; these are
workflow timings, not isolated solver performance comparisons.

Beyond x=68 m, the maximum simultaneously present escaped-water counts in
saved frames are 95 in 2D and 2,747 in 3D. Extra spray peaks at only 1 and 7,
respectively. These are sampled populations, not cumulative emitted counts.
The 3D plot projects escaped particles from all 72 m of width against a
central 12 m slice of bulk water; use Maya perspective views to judge the
actual spatial distribution.

All simulation steps meet the configured pressure tolerance of 1e-6. At
five-second checkpoint inspections, no liquid lies within one metre of the
ceiling. The largest inspected liquid displacement per step is 0.65 cells in
2D and 0.72 cells in 3D; these checkpoints do not bound every intervening
substep. The rear wall reaches liquid particle heights of 20.0 m in 2D and
23.4 m in 3D at 20 s, above the original 12 m surface. This is a material
boundary defect for scene selection, not something to hide in the camera.
The available evidence does not yet establish the intended open-sea
interactions or enough spray for the final shot.

Both Maya scenes pass world-bounds checks at six frames. Headless evaluation
of those samples takes 2–28 ms in 2D and 20–44 ms in 3D. The original 3D
scene-build attempt failed because Maya appended a suffix to camera names;
the builder now uses the actual returned names. No simulation was repeated
for that scene-only correction. Numeric receipts are in
[`results/3d/force_wave_preview/`](../results/3d/force_wave_preview/).

Production resolution, duration beyond the agreed 30 s, boundary treatment,
source strength and eventual start time remain subject to inspection. Sparse
spray at metre-sized cells cannot establish its production density. A successful
pressure solve also does not establish timestep or spatial convergence.

The first GUI inspection failed: almost no particles were visibly drawn despite
correct headless bounds. The numeric Maya receipts above are not a playback
acceptance. The user rejected the visible result and asked to prioritize a
denser, clearly readable 2D preview before any further coarse 3D work.

The replacement `out/runs/force_sea_2d_half_20261020/` has now completed:
360×64×1 cells at 0.5 m, dt=1/48 s, 30 simulated seconds and 720 playback frames
at 24 fps. It took 103.3 s including baking/export, with 95.2 s across the graph
calls, using eight threads. Initial bulk water is about 38,400 particles.
Beyond x=68 m, sampled escaped-water and extra-spray populations peak at 306
and 3. The rear boundary still piles up: the highest inspected liquid particle
in x<10 m reaches 23.73 m at 20 s against the initial 12 m surface.

The latest scene candidate is `scenes/wave_2d_points.ma` in that run, and its
`open_in_maya.bat` points there. The earlier `scenes/wave_2d.ma` remains a failed
display variant. The builder copies positions to a centimetre display cache,
uses identity parent transforms, and outputs points directly in wireframe.
The six-frame bounds check passes, but the latest candidate has **not** passed
GUI visibility or playback acceptance. Earlier scope/shaded variants failed
visually despite passing bounds. No production launcher or replacement
scheduled run is ready; original caches remain intact.

The [October 20 handoff](handoffs/2026-10-20-ocean-rendering-and-force-waves.md)
records the display diagnosis, exact candidate paths, rendering configurations,
remaining scene defects and the next inspection steps.
