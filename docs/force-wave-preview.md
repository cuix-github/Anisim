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
four playback layers. Blue is bulk water, cyan is escaped water, orange is
additional spray; grey air is hidden by default. Bulk-water display is limited
to 120,000 points per frame, and air to 10,000, solely for viewport speed.
Escaped water and emitted spray are not thinned. Original dumps retain all
particles. Viewport point sizes are visualization aids, not physical radii.

The Maya builder checks world-space bounds against saved metre-scale points
at six representative frames. It records evaluation times separately from
simulation timings; those headless times are not a measured GUI playback rate.
The 3D cameras use perspective. The source band, beach and still water are
useful scale references, not a rendered water surface.

Production resolution, duration beyond the agreed 30 s, boundary treatment,
source strength and eventual start time remain subject to inspection. Sparse
spray at metre-sized cells cannot establish its production density. A successful
pressure solve also does not establish timestep or spatial convergence.
