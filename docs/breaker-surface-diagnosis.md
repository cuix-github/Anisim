# Locating the breaker's rendering gap

On 2026-10-19, the user approved following the frozen 24.9 s curl from particles
through unsmoothed and smoothed reconstruction, including normals and a reusable
mesh. He also challenged the efficiency of trying one material parameter at a
time. The resulting experiment asks where structure changes, while the optical
diagnosis separates evidence from inference. It starts no solver or beauty render.

The shader audit already gives stronger priors than another roughness sweep.
The saved control uses Standard Surface with IOR 1.33, transmission 1, depth 40
scene units (10 m at 0.25 m per cell), and transmission colour (0.3, 0.55, 0.6).
It sets no transmission scattering. The surface has no explicit fine normal map.
The sky is a uniform dark blue-grey dome, accompanied by a directional sun and
dark ground. Curvature foam and the mist shape are disabled in this diagnostic.
The liquid-only render cache has no independent underwater bubble input.

Those facts suggest a clear absorbing water body with little internal fill and
missing whitewater cues, not a material that roughness alone can turn into a
breaking ocean. This is a causal hypothesis, not a measured allocation of the
image's darkness. Under the Standard Surface absorption model, five metres of
water transmit about (0.548, 0.742, 0.775) of incident RGB light before interface
effects: each colour component raised to 5/10. Absorption by itself therefore
does not establish why short paths look black. A flat water patch at the camera's
target has about 76.1 degrees incidence; ideal unpolarized Fresnel reflection at
IOR 1.33 is about 23.6 percent, compared with 2.0 percent at normal incidence.
That explains why a low view can emphasize reflections without any metalness.
Neither calculation includes the actual curl's normals, multiple interfaces,
shadowing or lighting, and neither predicts final pixel values.

The [Arnold transmission documentation](https://help.autodesk.com/cloudhelp/ENU/AR-Core/files/ac-shading/ac-surface-shaders/ac-standard-surface/arnold_user_guide_ac_standard_surface_ac_transmission_html.html)
distinguishes absorption depth, internal scattering, and transmission roughness.
Increasing specular roughness also increases transmission roughness when extra
roughness is zero. Thus the previous 0.03-to-0.12 arm was a legitimate single
parameter change, but it did not isolate reflection from refraction or introduce
whitewater scattering. [Arnold AOVs](https://help.autodesk.com/cloudhelp/ENU/AR-Core/files/ac-output-aovs/arnold_user_guide_ac_output_aovs_ac_aovs_html.html)
can output specular, transmission, volume, world-space normals and depth in one
render. A future approved optical check should use those layers, linear EXR output
and unclipped display exposure instead of another sequence of blind preset tweaks.

The geometry instrument is `experiments/breaker_surface_diagnosis.py`. It loads
the saved approved control scene, checks its cache and scene hashes, and verifies
that the surface compounds match those used in that render. It first reproduces
the 292,516-triangle control. The next stage changes only MCF iterations from 10
to 0, and the last changes only the mesher's smoothing from 0.1 to 0. All stages
retain the full source crop, point radius, voxel size, radius-based exclusion,
camera and lighting. There are already zero Gaussian passes in the control.

The measurement region is x=134-152 m, y=13-23 m, z=63-75 m, inside the source
crop rather than on its artificial cut faces. Exact closest-point queries use all
46,455 particle centres there, including 609 radius-flagged particles. The analysis
compares mesh connectivity, geometric and transferred shading normals, distances
and actual planar sections. A closest-normal side test is local: it is not a
global winding-number test or a liquid-mass measurement. Interior liquid centres
are not expected to sit on the rendered boundary.

The first export retained the Maya converter's default `transferNormals=false`.
Its geometry and closest-point arrays were written, but a bit-exact saved-scene
normal comparison failed. The next attempt enabled transfer and locking but still
used a bit-exact comparison immediately after locking. Both failed runs are kept
in `out/runs/breaker_surface_diagnosis_20261019/` and its `_v2` sibling. A cheap
static-scene probe measured the normal changes at about 2.4e-7 per component.
These were overly strict instrumentation assertions, not evidence of a defect in
the previous direct Bifrost-to-Arnold renders. The final export, in the `_v3`
sibling, transfers and locks normals and checks face-corner correspondence with
absolute tolerance 2e-6. Geometry remains bit-exact through baking and reopening.

The static Maya scene preserves the three stages, with only the control visible,
and the source camera, lights and material. It can avoid repeatedly
reconstructing the particle surface during look development. Geometry/normal
round-trip validation is not a beauty-render equivalence test; that comparison
must precede attributing later image differences solely to shading.

The final export passed from clean commit `07da114` in 327.27 s including process
startup. A separate static-scene verifier from clean commit `d2ffb9a` passed in
16.64 s including startup; opening and checking the scene took 11.26 s within it.
All three stages have identical positions and triangles after reopening, and
their maximum face-corner normal component error is 2.384e-7. The 132.9 MB scene
contains no Bifrost graph shapes. Both processes returned zero; their logs still
record MtoA teardown warnings after the completed checks. The runner now uses an
explicit 64-bit process handle for the existing scene-builder termination pattern;
that exit operation was separately smoke-checked. No geometry run was repeated
to test process cleanup. The two failed exports consumed another 542.42 s, and
are not hidden in the final export's clock.

| Whole cached surface | Raw reconstruction | Mesher smoothing 0.1 only | Current MCF 10 + smoothing 0.1 |
| --- | ---: | ---: | ---: |
| Triangles | 381,022 | 381,022 | 292,516 |
| Connected components | 309 | 309 | 89 |
| Surface area, square metres | 2,281.21 | 2,279.10 | 1,962.68 |
| Median adjacent-face angle | 8.81 degrees | 8.62 degrees | 1.47 degrees |

The mesher's 0.1 smoothing changes much less than MCF here. In the interior curl
region, 12.45 percent of the raw surface's 19,377 vertices lie more than one
12.5 cm surface voxel from the current surface; 4.28 percent lie more than 25 cm
away. The median distance is 2.75 cm, the 99th percentile 43.97 cm. Conversely,
only 0.73 percent of the current surface's 14,877 vertices lie more than 12.5 cm
from the raw surface. This asymmetric change is consistent with trimming protruding
features. These are vertex-weighted distances, not percentages of missing water,
surface area or meaningful detail. Component counts alone cannot distinguish
removal from merging, or real spray fragments from reconstruction noise.

Of 6,950 unflagged particle centres within 12.5 cm of the raw mesh in that region,
138 end up more than 12.5 cm on the local outside side of the current surface.
The closest-surface position shifts by a median 3.4 cm and a 99th percentile
34.7 cm for that set. This local-normal test normalizes the interpolated closest
normals first and does not replace a global inside test. Radius-flagged particles
are intentionally excluded before all three stages and require their separate
render representation. The unsmoothed reconstruction itself may already omit
structures; these comparisons do not quantify that earlier loss.

The meshes are closed in their edge incidence, with no nonmanifold edges. The
current surface has no degenerate triangles or nonfinite normals. About 0.040
percent of its face-corner normals oppose their geometric face normal. That rules
out a global reversal but not local artifacts. An Arnold N AOV is still needed to
verify the final shading normals in the actual render path. The section figure at
z=68 m shows smoother crests, smaller detached contours disappearing and a rounded
lip. Particle dots use a 0.25 m slab, while the mesh curves are exact planes; they
are not identical sampling measures.

The evidence now prioritizes MCF's structural changes and the missing optical
components over a new simulation or a whole-grid refinement. It does not justify
turning smoothing off everywhere: the raw surface contains noise too. The next
useful image experiment is an approved AOV diagnostic on reusable geometry,
with a conversion control, before selecting a physically coherent water/whitewater
look. No such render was launched. The numeric receipt is
`results/rendering/breaker_surface_20261019.json`; the local section plot and scene
are `curl_sections.png` and `baked_breaker.mb` in the `_v3` output directory.
