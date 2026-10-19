# Surface foam representation and angular fragments, 2026-10-19

Xue Cui approved two close side previews after inspecting the foam distribution,
thickness and boundaries: a uniform thin layer and a clustered layer, using the
24.9 s cache, corrected baseline, mist off and no simulation. Work began at
19:17:15 UTC, with advance notice required if it was expected to exceed 70 minutes.

Xue then clarified the interpretation of the previous SSS image. Its white blocks
are useful evidence that the appearance is moving in the right direction, not a
reason to reject whitening itself. A wider shot may read those blocks as fine
whitewater. Near-view material quality and distant readability should therefore
be evaluated separately. Downsampling a completed image can test readability,
but it is not a new wider camera shot and does not add the missing surroundings.

## What is being compared

The macro support comes from the already checked curvature mask. Equal-area
surface sampling produces 324,601 quadrature points. Each represents a small
render volume, not a physical bubble or a newly simulated liquid particle.
The underlying water mesh, camera, optics and lighting remain frozen. The
previous SSS replacement is disabled: the water uses corrected baseline again.

Point spheres have radius 2.25 cm, with centres shifted 0.75 cm along the local
surface normal. Their nominal extent is 1.5 cm inward and 3 cm outward. Every
sphere is within 3 cm of its generating surface point; this is a geometric bound,
not a claim of a measured foam depth in the simulation. Sampling the interpolated
mask alone initially allowed some support into the crop-exclusion band. The
analytic exclusion is now reapplied at each sample. The nearest active sphere
is at least 58.57 cm from a crop plane. Side, top and point-level plots were
inspected before Maya preparation. Broad calm water remains outside the support.

A distributes its density uniformly subject to the macro mask and quadrature
weights. B redistributes it using deterministic fields with characteristic
lengths 24 cm and 6.5 cm; about 62.3 percent of its logical points have zero
density. At the 15 m / 960 px view, one pixel is 1.5625 cm. These field lengths
describe foam aggregates, not individual real bubble diameters. They are
appearance choices rather than values supplied by the paper or cache.

The normalization matches the sum of each point sphere's extinction coefficient
times its volume, 16.81893 square metres. This is a primitive-sum budget. It is
not rendered energy, screen coverage, or proof of the same actual optical depth
where Arnold overlaps volume primitives. Both use 5 mm ray steps, white scatter,
scatter weight 0.98 and anisotropy 0.2. The target normal optical-depth parameter
0.8 is a construction choice, not a measured value. No new airborne mist is added.

## The first A did not exercise foam scattering

The first native point cloud stored density as uniform per-point user data and
linked a user_data_float node to Standard Volume density. Native Arnold scene
export/readback preserved every position, radius and density value. Maya's
stand-in export referenced that file and retained its embedded shaders.

Nevertheless, the first licensed GUI render at commit 7d988b1 had an exactly
zero volume AOV. It took 12.51 s inside the render call and 22.36 s including the
process; mean absolute raw-beauty difference from the existing baseline was only
5.88e-7. The image check failed and B was not launched. The process receipt and
failed image check are both preserved. The layer contributed no visible volume;
density evaluation was not isolated by that check. This does not prove Arnold
generally cannot read per-point data for volumes, which its documentation says
is supported. Export/readback alone did not validate live shading.

A replacement setup uses 64 explicit constant-density shader groups instead
of per-point shader lookup. It preserves the logical positions and radii, omits
zero-density spheres, and replaces each group's density by its mean. Density
quantization has relative L1 error 1.17 percent in A and 1.88 percent in B, while
the primitive-sum budgets remain matched to better than 3e-9 relative error
after native readback. Shader assignments and positive constant densities were
also inspected directly. This is a prepared alternative, not a validated look.

The two replacement scenes passed geometry, optics, camera and export checks.
Because producing both valid previews after the failed A would spend three
render calls instead of two, approval was requested for the extra call. At this
point, replacement A and B had not rendered. Xue subsequently approved replacement
A, with B conditional on a nonzero scattering check and an explicit stop if A
failed. Replacement A ran at 19:50:49 UTC on October 19 at commit ee3a070. Its
licensed render completed in 12.00 s (21.09 s including the process), but its
volume AOV was again exactly zero. The image check failed, and B was not launched.
There have been two render calls in this thin-foam study, both ineffective A
attempts, and no valid foam comparison. No simulation ran.

A subsequent read-only audit identified a concrete integration error. The native
foam ASS header declares meters_per_unit=1; the Maya parent ASS declares 0.01.
The points were already expressed in the parent's coordinate convention, so
importing the asset introduces an unintended factor of 100. Native maximum
bounds (635.233215, 66.927734, 310.057190) predict imported maxima
(63523.3215, 6692.77344, 31005.7190); the actual render log reports
(63523.3203, 6692.77344, 31005.7188). The layer is displaced and scaled away from
the wave. The first failed A log shows the same oversized bounds.

This accounts for the absent in-frame foam without establishing a defect in
per-point density support. The earlier export/readback and identity-transform
checks missed the parent/child unit conversion. Changing density representation
did not address that error. The [Arnold units documentation](https://help.autodesk.com/cloudhelp/ENU/AR-Core/files/ac-render-settings/arnold_user_guide_ac_render_settings_ac_units_html.html)
describes this automatic procedural scaling. The next prerequisite is to align
the asset's unit metadata with the frozen parent scene and check imported bounds
and camera projection before any new render. Neither failed artifact was altered;
no unit-corrected scene or further render was produced after the conditional stop.

## Unit-aligned comparison completed

At 19:58:22 UTC on October 19, Xue approved correcting units, checking actual
imported alignment, then rendering A and conditionally B. Advance notice was
required if the work was expected to exceed ten minutes. The native ASS now
declares 0.01 meters per scene unit, matching the frozen Maya parent. This
matches the parent DCC convention; it does not change the inherited mapping of
four scene-coordinate units per simulation metre, or rescale the baseline.
The layer data hash, positions, radii, density groups and shaders remain unchanged.

An initial preparation attempt in _v5 stopped without rendering because Maya
standalone returned an empty display bounding box for the unloaded stand-in.
The replacement check uses public AiProceduralExpand on the exported parent ASS,
with no ray tracing. It verifies every active point and radius, the expanded
transform chain and distance to its original surface anchor. Both arms passed
in _v6 at commit fdf0e38. Maximum centre-to-anchor distance is 7.50843 mm;
including sphere radius, the bound is 30.00843 mm (float32 rounding included).
The former unit-mismatched _v4 asset is rejected by the new preflight gate.

A completed first and passed its nonzero volume check before B was launched.
Both licensed 960 by 540 images pass finite EXR, component closure and independent
ACES display checks. Volume contributes 10.0027 percent of summed raw RGB in A
and 7.4918 percent in B; these are whole-image radiance shares, not foam coverage
or physically measured foam fractions. SSS remains zero. A took 185.39 s inside
rendering (194.63 s process); B took 107.53 s (116.35 s process). The earlier
ineffective A images were much faster because the layer missed the wave. At
20:04:55 UTC, before the ten-minute limit, the agent revised the total estimate
to 12–15 minutes based on the now-active volume work.

Visual review finds a smoother continuous white coating in A and visibly broken,
granular patches in B. B changes the coverage texture, while retaining the same
large underlying fragments and glassy main water. It is a valid comparison and
a working rendering path, not an accepted ocean preset or a claim of reference
quality. A/B scene states differ only in foam asset and arm identifiers. Cache,
water geometry, corrected-baseline optics, close side camera and mist-off state
remain fixed. No simulation or remeshing ran. This turn spent exactly two render
calls; the whole thin-foam study has four calls, including the two earlier failed
A attempts. The images, smaller-display comparison and receipts are in _v6.

## The circled hard fragments are a separate geometry finding

Xue supplied a marked image highlighting angular free-surface fragments. An
offline audit uses the already saved raw, no-MCF and control meshes, not a new
reconstruction or render. The floating diamond near original image pixel
(463, 194) has six vertices and eight triangles in both raw and current meshes.
Its current world-space bounding extent is about 17.1 by 11.3 by 11.2 cm;
before smoothing it was about 22.3 by 13.5 by 14.2 cm. Two cached particle
centres lie within the expanded local bounding box. That is proximity evidence,
not exact particle-to-mesh provenance.

The angular silhouette is therefore in the reconstructed geometry already;
lighting or changing normals cannot round its outline. This identifies a
render-reconstruction limitation for that fragment, not a solver failure or a
diagnosis of every attached crest corner. The larger attached protrusions still
need their own checks. A later comparison could reconstruct selected fragments
more finely from the same cache, or represent isolated droplets separately. That
would be a separate geometry experiment; this foam comparison does not silently
alter the water mesh to hide the issue.

The geometry audit and close-up plot are in
out/runs/breaker_corner_audit_20261019/. The first failed support build is in
out/runs/breaker_foam_layer_20261019/. The _v2 folder retains the valid layer and
a scene-preparation tuple/list comparison failure; _v3 retains the completed
but ineffective A. The _v4 folder holds the ineffective replacement A, its failed
image_check.json and unit_audit.json, and the unrendered B. These outputs are
local, not committed; compact evidence is in results/rendering/.

[Arnold points](https://help.autodesk.com/cloudhelp/ENU/AR-Core/files/ac-shapes/arnold_user_guide_ac_shapes_ac_points_html.html)
documents volume point rendering.
[Standard Volume](https://help.autodesk.com/cloudhelp/ENU/AR-Core/files/ac-shading/ac-volume-shaders/arnold_core_ac_standard_volume_html.html)
documents density and scattering controls.
[PF-FLIP Section 7](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf#page=15)
states that no procedural textures were used. The clustered layer here is an
explicit render approximation, not a reproduction of that paper's full pipeline.
