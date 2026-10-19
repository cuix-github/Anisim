# A preset-derived breaker candidate, 2026-10-19

The two authorized 960 by 540 previews completed using the existing 24.9 s
cache and frozen control mesh. This candidate is **not an accepted replacement
look**: the original view is more saturated blue and darker, and the broad
glassy sheets remain. The side view exposes the cropped slab boundary and
provides a silhouette check, not the reference's kilometer-scale composition.
No solver, surface reconstruction, spray, or curvature-foam pass was run.

## Candidate and provenance

The existing control was validated in [PR #39](https://github.com/cuix-github/Anisim/pull/39).
Its `B/optical.mb` was reused, with its SHA-256 pinned in the launcher. The
original-camera candidate changes a coherent material/lighting package; it is
not a one-variable causal comparison. The second image changes only the camera.
The user's agreed two-image scope reuses the existing control instead of spending
another render reproducing it. All input hashes are checked before and after.

| Setting | Value | Basis |
| --- | --- | --- |
| Standard Surface base / metalness | 0 / 0 | Water dielectric |
| Specular / transmission weights | 1 / 1 | Shipped `aiStandardSurface/Deep_Water.mel` |
| IOR / roughness | 1.333 / 0.12 | Same preset |
| Transmission color | (0.744, 0.9400832057, 1) | Same preset |
| Transmission scatter | (0.0675, 0.1340401471, 0.5) | Same preset |
| Scatter anisotropy | 0.75 | Same preset |
| Transmission depth | 10 solver-cell units = 2.5 m | Agreed candidate optical scale, not a universal water constant |
| Physical Sky elevation / azimuth | 30 / 160 degrees | Candidate direction |
| Turbidity / sun size | 6 / 5 degrees | Agreed candidate; enlarged sun is an artistic softening |
| Sky / dome intensity, dome exposure | 1 / 1, 0 | Fixed for both views |
| Dome resolution | 2048 | Importance sampling map, not spray density resolution |
| Separate old directional sun | Disabled | One sky/sun system |
| AA / specular depth / transmission depth | 4 / 3 / 8 | Validated optical baseline, ray counts not meters |
| Volume depth / samples | 3 / 2 | Preserved |
| PNG / EXR | ACES view / raw ACEScg | Independently compared offline |

The installed preset does not declare the color space in which its RGB constants
were authored. This candidate interprets them directly in the working ACEScg
space; they are reference-derived shader settings, not measured ocean optical
coefficients. See Arnold's [Transmission](https://help.autodesk.com/cloudhelp/ENU/AR-Core/files/ac-shading/ac-surface-shaders/ac-standard-surface/arnold_user_guide_ac_standard_surface_ac_transmission_html.html)
and [Physical Sky](https://help.autodesk.com/cloudhelp/ENU/AR-Core/files/ac-shading/ac-texture-shaders/arnold_core_ac_physical_sky_html.html)
documentation for the mechanisms and scale dependence.

The side camera is orthographic, with zero pitch and roll and yaw 168 degrees.
Its horizontal sightline is perpendicular to the dominant +12 degree wave train
(verified dot product 2.78e-17). Other wave trains have different headings. Its
width is 183.02936 cell units, about 45.76 m, to fit the existing mesh bounds.
It does not create an ocean extension or move the scene. Both scenes preserve
the original black camera background and the beach/ground geometry.

## Results and checks

| Measurement | Existing camera | Side camera |
| --- | --- | --- |
| Arnold render call | 78.20 s | 17.60 s |
| Maya process including startup/setup/collection | 89.57 s | 27.35 s |
| PNG versus independent OCIO view, mean 8-bit difference | 0.179 | 0.144 |
| Maximum display difference | 1 level | 1 level |
| Raw component sum relative L1 error | 9.09e-8 | 9.33e-8 |

Preparation, inspection and documentation are additional to those process times.
Both Maya sessions exited zero, wrote licensed 54-channel float EXRs and complete
logs, and passed output collection without recovery or a repeat render. Arnold
warned that the original view traces 2365 rays per pixel. Sampling was not
silently reduced or increased to obtain another image.

The visible mesh is still 146,326 vertices and 292,516 triangles. World points,
stored normals and topology hashes are identical before/after the candidate
changes and between cameras. The original view's P AOV is bit-identical to the
baseline. Its N AOV is not bit-identical: maximum absolute difference is
1.79e-7, consistent with floating-point roundoff rather than an edited normal
field. The raw, undenoised beauty retains the same large sheets and blue cast.

The saturated blue is concentrated in transmission paths. Whole-image mean
linear RGB for the candidate's specular AOV is (0.01425, 0.01432, 0.01402), while
transmission is (0.00513, 0.01008, 0.04744). Baseline transmission was
(0.03481, 0.04915, 0.05555). These are image statistics including background and
non-water geometry, not isolated material reflectance. AOV labels classify the
first scattering event; transmission can contain later reflections and internal
scattering. A zero separate `volume` AOV does not mean the water's internal
scattering is disabled.

The package did not solve the target appearance. Its combined changes do not
establish which individual setting caused every difference, nor prove that all
uses of Deep Water fail. No further parameter sweep was launched. The corrected
previous baseline remains the better visual reference.

## Relation to PF-FLIP Section 7

The [paper, Section 7 and Table 7](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf)
describes a general spray/rendering pipeline, not a camera mask that makes wave
crest pixels higher resolution. Escaped liquid particles become independently
integrated droplets; additional droplets are emitted near interfaces with steep
velocity gradients. They are represented at a finer spray scale and rasterized
into a sparse density field, at two to four times the linear pressure-grid
resolution. A volumetric renderer handles scattering; liquid particles separately
provide a reconstructed surface. High-curvature interfaces receive foam-like
scattering. This is not the same as a persistent simulated surface-foam layer.

The breaking-wave table lists pressure 3072 by 1536 by 1024 and spray density
12288 by 6144 by 4096. Across the stated 1 km width, those spacings are about
0.326 m and 0.0814 m. Four times the linear density resolution would cost 64
times as many voxels for a dense fixed volume, which explains the importance
of sparse storage. This does not require a pressure solve on that finer grid.

Our implementation borrows this separation, but is not a reproduction of the
authors' renderer:

- `maya/generate_particle_surface.py` reconstructs the liquid surface and creates
  a separate `particle_mist` volume. Mist selection uses cached droplet radii
  and particles outside the smoothed surface. It splats weighted clouds into
  `voxel_fog_density`, with an independently configurable `spray_detail`.
- `maya/build_surface_scene.py` assigns that field an Arnold Standard Volume
  shader, while a separate Standard Surface shades the mesh. Curvature foam is
  currently a white surface-shader blend, not the paper's specified interface
  scattering implementation. Density scaling and kernel/cloud size are our
  documented choices rather than parameters recovered from the paper.
- The prior configured mist voxel was 0.5 cell = 0.125 m: finer than our 0.25 m
  pressure cell, but coarser than the reference spray density. These two new
  previews contain no active mist volume. The cache's 1,192 positive-radius
  particles are excluded from the surface, so hiding mist also hides their
  intended representation. Other escaped particles can remain in the surface
  after cache classification; 1,192 is not the number of all escaped particles.
- `experiments/wave_basin.py` explicitly runs without the per-droplet model or
  additional spray emission. The implementation exists elsewhere, but its
  presence in the solver is not evidence that this cache exercised it. The cache
  also applies the `spray_above_cells=1` filter.

Arnold combines surfaces and volumes in one scene: rays sample the fine density
where it exists and interact with the water mesh where it is hit. Its volume
[step size](https://help.autodesk.com/cloudhelp/ENU/AR-Core/files/ac-shapes/arnold_user_guide_ac_shapes_ac_volumes_html.html)
must resolve the field; it is distinct from density voxel size, pixel resolution,
AA samples, and light samples. The renderer need not recognize a crest visually.

A finer volume can be rasterized from existing particle positions without
rerunning the pressure solve. That alone does not add spray particles, reconstruct
missing motion, or shrink the splat kernels: sampling a broad cloud more finely
still gives a broad cloud. Capturing extra droplet dynamics/emission is a separate
question from rendering the existing cache and is not authorized by this study.

## Local artifacts and reproduction

`out/runs/breaker_deep_water_20261019/` contains `original/` and `side/`, each
with its scene, PNG/EXR, state, input hashes, receipts and logs. `comparison.png`
includes the corrected prior baseline; `aov_panels.png` shows the first view's
components. Original PNGs remain untouched; `*_view.png` are verified offline
conversions. Raw imagery and scenes stay outside Git.

The compact committed evidence is `results/rendering/breaker_deep_water_20261019.json`.
The render instrument was committed at `f571d24` before either measurement.

Offline validation, with local artifacts present:

```
.venv\Scripts\python.exe experiments\analyze_breaker_deep_water.py
```

Separately authorized reruns use `experiments/breaker_deep_water.py original`
and then `side`, both with the same fresh `--out` directory. The launcher checks
a clean checkout and existing control, refuses overwrites and another Maya
process, and bounds each owned hidden GUI Render Sequence process to six minutes.
