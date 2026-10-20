# Six cached-wave views with 20/80 foam

On October 19, 2026, Xue reacted enthusiastically to the 25/75 intermediate
preview and authorized six more images using 20 percent uniform and 80 percent
clustered foam. These views sample 23.025 s (lip before impact), 24.9 s
(breaking) and 27.025 s (broken water advancing shoreward), with a 56 m wide
view and a 22–24 m detail at each time. They use existing liquid/radius caches;
no solver, mist, additional emission or procedural ocean was enabled.

The mixture is applied to the spatial extinction density, C = 0.2 A + 0.8 B,
with the same integrated primitive extinction budget as A and B for each frame.
It is not an image blend, a coverage fraction or a physically calibrated foam
model. Corrected-baseline water, lighting, 960 x 540 sampling and ACES display
remain fixed. The beach collider's x/y profile is unchanged; its render-only
z extent expands from 60..78 m to 0..80 m to match the wider water.

## Geometry preparation and validation

Five reconstruction attempts were stopped by the 5 GiB free-memory guard and
retained. Overlapping tiles demonstrated a watertight seam with a maximum
0.00763 mm weld discrepancy, but smaller tiles still showed unexpectedly high
memory use. The successful images therefore do not use stitched geometry.

A public-node bounds audit found that an approximately 8.4 m-high particle
region produced volume_bounds spanning 108 m vertically. Our graph fed those
loose bounds into field_to_volume for the dense density mask. Restricting only
that field to particle bounds plus 8 m of empty padding, aligned to the same
0.5-cell voxel spacing, reproduced the completed 155,546-vertex control exactly:
all positions and triangle indices match. Raw particle spheres, level set,
10 masked-MCF iterations and stock meshing remain unchanged. The paired control
fell from about 80 s to 12.02 s; this is one measured preparation case, not a
general solver speed claim. [Public volume_bounds reference](https://help.autodesk.com/cloudhelp/2023/ENU/Bifrost-Common/files/reference/geometries/Bifrost_Common_reference_geometries_volume_bounds_html.html).

The final full-width reconstructions take 20–26 s each with more than 36 GiB
free RAM remaining. A disconnected six-vertex/eight-face component at 27.025 s
had collapsed to a 4.3-micrometre extent, with zero-area faces and undefined
normals. It alone was removed; original arrays remain preserved. No visible
fragments were rounded or replaced. Boundary-crossing foam triangles are
excluded so vertex-mask interpolation cannot spill across the cut exclusion.

Primary mesh-intersection diagnostics found no visible cut faces in any wide
view or the shore detail. Pre detail has 0.30 percent and break detail 0.16
percent of sampled water hits near a remote boundary silhouette. These are
coarse primary-ray samples, not refracted-ray coverage or proof that every
boundary artifact is gone. Expanded-procedural checks verify actual foam
positions, radii, transforms and centimetre units before rendering.

## Composition and limits

The first 65-degree high view flattened the unbroken wave and devoted too much
of the image to shallow water. It is retained as an overhead inspection view.
The two remaining, unrendered wide cameras were lowered to 40 degrees and given
more upstream water before their first calls. Exactly six renders are allowed;
no repeat of the first image is part of this study.

The foam mask is instantaneous curvature/height selection and has no transported
history. These stills do not establish temporal consistency or a finished ocean
look. The physical domain is finite; 56 m framing is wider context, not a
kilometre-scale ocean or the entire 172 m simulation domain.

## Render results and visual review

All six licensed GUI-Maya renders passed on their first calls. The final image
checks completed at 00:28 UTC on October 20 (October 19 locally). Render calls
total 864.59 s, with 919.42 s including process startup. Preparation and checking
occupied the larger part of this roughly seventy-minute task; the revised
80–100 minute estimate was not exceeded. No seventh render was started.

Every image has nonzero volume and zero SSS. The maximum component-closure
relative L1 error is 8.77e-8. Independent ACES display conversion agrees within
0.264 mean 8-bit level, and every source-cache hash remains unchanged.

| Cached time | Wide view / render seconds | Detail view / render seconds |
| --- | --- | --- |
| 23.025 s | 56 m, 65-degree overhead / 19.73 | 24 m oblique / 47.66 |
| 24.9 s | 56 m, 40-degree elevated / 90.14 | 22 m oblique / 360.18 |
| 27.025 s | 56 m, 40-degree elevated / 106.72 | 24 m shore-facing / 240.15 |

Independent visual review favors the 24.9 s detail: a longer breaking front
retains granular foam and open gaps. The pre-break detail shows the projecting
lip with much less whitewater, exposing the still-glossy bulk water. The first
high view is an inspection image, not a recommended presentation camera. The
shore wide view also exposes the finite beach model's straight end. The shore
detail shows a particulate white front with little lingering foam behind it;
these stills cannot diagnose transport, but the current mask has no foam history.
The batch does not justify declaring the full ocean appearance solved, and it
awaits Xue's review. The earlier enthusiastic reaction was to the 25/75 image.

Outputs are under `out/runs/breaker_six_views_20261019/`: `comparison.png` is the
compact contact sheet and `comparison_full.png` preserves each 960 x 540 image.
Individual images are `shots/{pre,break,shore}_{wide,detail}/C/beauty_view.png`.
`final_plan.json` records the actual final cameras, geometry and source hashes;
`plan.json` preserves the initial proposal. Each shot keeps its Maya scene,
expanded-geometry checks, EXR/AOVs, logs, image checks and receipt. Failed
reconstructions and superseded unrendered preparation folders remain intact.
Numeric evidence is `results/rendering/breaker_six_views_20261019.json`.

The relevant instruments are `breaker_six_views.py`, `breaker_bounded_mesh.py`
and `breaker_camera_probe.py`. `breaker_tiled_mesh.py` and
`breaker_tile_bounds.py` retain the memory investigation. The tiled and unbounded
paths are diagnostic history, not the path used by these six images. New
rendering requires a new user-approved scope; existing receipts prevent spending
another call on any of these six shots.


## Fixed-camera 60/40 follow-up

A subsequent request authorized one 60 percent uniform / 40 percent clustered
preview of the existing 24.9 s, 22 m detail. `breaker_foam_mix60.py` reuses its
frozen scene and samples, changing only extinction density to 0.6 A + 0.4 B.
Camera, geometry, normals, water, lighting, foam radius and total primitive
extinction are unchanged. This is a density mixture, not an image blend.

The offline control first stopped because its exact comparison used `1 - .8`
where the source used literal `.2`. Matching the source arithmetic reproduced
the control exactly; no density export or render preceded that correction.
The note remains in the run folder. Expanded-import checks again verify all
984,932 active samples, identity transforms and a maximum 3.001 cm extent
from the source surface. The 64-bin density export has 3.06 percent relative
L1 quantization error and 1.71e-9 relative total-budget error.

The single licensed render passed in 404.22 s (413.73 s with process startup),
with nonzero volume, zero SSS, component closure 8.60e-8 and mean display
difference 0.2402 of an 8-bit level. Cache hashes remain unchanged. There was
no simulation, mist, remeshing or repeat render.

Visual inspection finds more connected and softer whitewater, with some fuller,
more matte patches. The earlier 20/80 retains stronger granular breakup. The
change does not remove the glossy bulk water or add transported foam history.
The choice between these appearances remains open for review.

The new image is `out/runs/breaker_foam_mix60_20261019/C/beauty_view.png`.
`comparison.png` pairs the original on the left with the new image on the right;
`comparison_vertical.png` preserves the same order from top to bottom.
Numeric evidence is `results/rendering/breaker_foam_mix60_20261019.json`.
