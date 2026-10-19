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
Its geometry and closest-point arrays were written, but the saved-scene normal
comparison failed; it is retained as a failed run in
`out/runs/breaker_surface_diagnosis_20261019/`. This is an instrumentation problem,
not evidence of a problem in the previous direct Bifrost-to-Arnold renders. The
corrected export enables normal transfer, locks the baked corner normals, and
compares them by face corner after reopening, allowing normal-table reordering.
Its outputs are separate in `out/runs/breaker_surface_diagnosis_20261019_v2/`.

The static Maya scene preserves the three stages, with only the control visible,
and the source camera, lights and material. Once validated, it can avoid repeatedly
reconstructing the particle surface during look development. Geometry/normal
round-trip validation is not a beauty-render equivalence test; that comparison
must precede attributing later image differences solely to shading.
