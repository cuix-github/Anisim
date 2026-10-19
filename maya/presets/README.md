# Render presets: the knobs and what each one controls

A preset is a PowerShell file here that sets the launcher's environment knobs, each only when the caller has not set it, so `render_surface_preview.ps1 -Preset paper` gives the paper look and an explicit variable on top of it changes one thing. The point is a fixed instrument: when two simulations are compared, the render is the same and only the simulation differs.

`paper.ps1` is the look settled on 2026-10-16/17: the breaking-wave stills of PF-FLIP and ST-FLIP on our particles. Two knobs belong to the scene and are not in the preset: `SURF_METRES_PER_CELL` (the simulation's cell in metres, 0.25 for the breaker) and `SURF_CAMERA`.

`breaker_inspect.ps1` is the water-only diagnostic candidate from the [24.9 s comparison](../../docs/breaker-render-study.md), 2026-10-19. It bypasses curvature foam, hides the existing mist, uses 10 MCF iterations and doubles the original ocean environment RGB. It exposes sheets and holes but still looks glassy; it is not a replacement for `paper_clean` or an accepted realistic ocean look. Hiding mist also hides its escaped particles; it never changes the cache or returns those particles to the surface.

`paper_clean.ps1` is the same look, clean on calm water (2026-10-18). On the breaker every frame was a break, so the paper look's rule that any water the smoothed surface leaves out is a droplet cloud read as white water. On the wave basin's calm 2.8 m swell at 18 s, far from breaking and with no wind, the same rule painted the surface as snow: a real swell like it is glassy. Two knobs remove the dots, one render and one cache:

| Knob | paper | paper_clean | What it removes |
| --- | --- | --- | --- |
| `SURF_MIST_THRESHOLD` (render) | 0 | 0.75 | water less than 0.75 cells outside the smoothed surface, left out only because the smoothing irons single-particle bumps away |
| `--spray-above` (cache: `experiments/wave_basin.py cache` or `section`) | 0 | 1 | escaped particles within one cell of the bulk water: a surface particle poking above its neighbours sees mostly air and trips the solver's escape test without flying |

The test, one knob at a time on the basin's 18 s frame with the paper look's camera framing (`out/runs/wavebasin4/renders/full_height/S0_t18`, `S1_t18`, `S3_t18`): the threshold alone took out most dots; the 502 flagged particles left in the slab all sat within a cell of the water, and the cache knob took them out, leaving none. Not yet tested on a break: real spray leaves the crest fast and should be more than a cell above the water within a frame or two, which the basin's break frames are to show (arms S0 and S3 of the night of 2026-10-18). The cache knob is not an environment variable because the play cache's `radius.####.npy` carries it; a cache made without it renders as the paper look's spray whatever the preset.

## The surface (`maya/particle_surface`, ST-FLIP Appendix B)

| Knob | What it controls | Where it acts | Paper |
| --- | --- | --- | --- |
| `SURF_POINT_SIZE` | the sphere each particle becomes, in cells; the union of those spheres is the raw surface | `points_to_volume` | 0.5 dx |
| `SURF_DETAIL` | the surface grid's voxel in cells; 0.5 is twice the simulation's resolution | `points_to_volume` | dx / 2 |
| `SURF_SMOOTH_ITERS` | Gaussian passes over the whole level set before the flow; each erases what is thinner than about two voxels | `smooth_level_set_property` (Gaussian) | none |
| `SURF_MCF_ITERS` | iterations of mean-curvature flow under the mask; 0 is off, more rounds bulk edges more (the block's top edge: 0.34 cells at 30); the paper's one free parameter | `smooth_level_set_property` (MeanCurvature) | k_psi = 30 |
| `SURF_MASK_SIGMA` | the blur of the mask's density, in cells; wider protects larger features | `smooth_voxel_property` | 2 dx |
| `SURF_MASK_THETA` | where protection starts: a feature whose density is theta times its blur begins to be held | the field math | 2 |
| `SURF_MASK_ZETA` | where protection is complete | the field math | 5 |
| `SURF_MASK_FLOOR` | the mask's floor; 1 runs the flow everywhere (the unmasked arm) | the field math | 0 |
| `SURF_MASK_DILATE` | the mask's numerator pushed out by this many voxels so both sides of a kept feature's band are held; our departure, measured in the gate | the field math | 0 |
| `SURF_MESH_SMOOTH` | the mesher's own smoothing of the triangles | `volume_to_mesh` | n/a |
| `SURF_OFFSET` | pushes the surface out (+) or in (-), in cells | `points_to_volume` | 0 |

## The spray (`maya/particle_mist`, PF-FLIP Section 7)

| Knob | What it controls | Paper |
| --- | --- | --- |
| `SURF_MIST_THRESHOLD` | the depth above which a particle the surface leaves out is mist; every escaped droplet (a radius) is mist regardless | the spray particles |
| `SURF_MIST_DETAIL` | the spray grid's voxel in cells | 2 to 4 times finer than the pressure grid |
| `SURF_MIST_RADIUS` | the cloud each mist particle is spread over, in cells; wider is softer and fainter | n/a |
| `SURF_MIST_DENSITY` | the fraction of a parcel's water drawn as droplets of its radius; 1 is the whole parcel and opaque, 0.001 reads as spray; a choice, recorded in `docs/open-questions.md` | not stated |
| `SURF_DEFAULT_RADIUS` | the droplet radius in metres for a mist particle without one (a cache without the droplet model) | 100 um median |
| `SURF_METRES_PER_CELL` | the cell in metres; the droplet clouds' extinction is computed from it | the scene's |

## The look (`maya/build_surface_scene.py`)

| Knob | What it controls |
| --- | --- |
| `SURF_LOOK` | `ocean`: dark absorbing water, a near-black sky the camera does not see, a dark floor; anything else is the pale default |
| `SURF_WATER_DEPTH` | how far light travels in the water before the transmission colour takes over, in cells |
| `SURF_WATER_COLOR` | the transmission colour |
| `SURF_WATER_ROUGHNESS` | water specular roughness, 0 to 1; default 0.03 preserves the original; the proposed 0.12 arm is not visually validated |
| `SURF_DIAGNOSTIC_MATTE` | 1 assigns a neutral diffuse material to the existing surface geometry for inspection; default 0 preserves normal material assignment; mist visibility remains separate |
| `SURF_SUN_ROTATE` | the sun's direction; the preset puts it low behind the spray so the mist is lit from behind |
| `SURF_SUN_INTENSITY`, `SURF_SKY_INTENSITY` | the two lights |
| `SURF_SKY_VISIBLE` | whether the camera sees the sky dome |
| `SURF_SKY_COLOR` | optional linear RGB override for the dome; its intensity is a separate multiplier |
| `SURF_FOAM_STRENGTH` | 0 bypasses curvature foam; 0 to 1 scales its curvature multiplier; default 1 preserves the original |
| `SURF_MIST_VISIBLE` | 0 hides the existing mist shape; default 1 shows it; particle classification is unchanged |
| `SURF_CAMERA` | "x y z tx ty tz" in cells |

The white matte foam material is mixed into the glass by convex curvature (`aiCurvature`, radius 1.2, threshold 0.15), the render-time rule the papers use. The strength control defaults to the original multiplier of 1.6; zero assigns the water shader directly. Bubbles are not drawn (`docs/open-questions.md`).

## Arnold (`maya/render_surface_preview.py`)

`SURF_WIDTH`, `SURF_HEIGHT`, `SURF_SAMPLES` (AA samples), `SURF_THREADS` (0 = all; 8 is the daytime cap on the shared machine), `SURF_GPU`, `SURF_GI_VOLUME_DEPTH` (scattering bounces in the mist), `SURF_VOLUME_SAMPLES`.
