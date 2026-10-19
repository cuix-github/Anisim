# Curvature-selected interface scattering, 2026-10-19

Xue authorized a mask check followed by two close side-view comparisons using
the corrected baseline, the existing 24.9 s cache and the frozen surface. Mist
stays off and no solver is run. He asked for advance notice if the work looked
likely to exceed 70 minutes. Work began at 18:10:59 UTC.

The paper's Section 7 distinguishes extra spray dynamics from rendering. Escaped
FLIP particles become droplets, and additional particles are emitted at steep
interface velocity gradients. High-curvature interfaces receive foam-like
scattering coefficients at render time. It does not describe a separate evolving
surface-foam simulation. Our solver has gradient-triggered emission, but wave_basin
did not enable it or the per-droplet model for this cache. Existing escaped
particles and additional emitted particles are different populations.

Xue wants to keep investigating rendering first. If this cache still cannot give
the desired fine crest spray, a future fourth- or fifth-segment continuation can
be considered. The choice must use the existing motion/crest history and the time
available for spray to accumulate before the target frame. Starting at 24 s gives
only 0.9 s before 24.9 s. This discussion is not authorization to run a continuation.

## The mask passed; that does not validate the volume

experiments/breaker_interface_mask.py computes signed cotangent mean curvature
on the frozen 146,326-vertex, 292,516-triangle mesh in metres. Four normalized
neighbor diffusion passes smooth the scalar, not the surface or normals. The
analytic sphere check gives area-weighted H = 1.000012 for a unit sphere; the
plane gives zero. All mesh edges have valence two; self-intersections were not
ruled out by that test.

The mask ramps from zero to one over absolute mean curvature 2 to 4 per metre.
It excludes artificial X/Z crop boundaries with a 0.6 to 1 m fade and uses an
explicit height gate from 12.25 to 13 m to focus on the curling interface. The
height gate is a shot-specific restriction, not a physical foam model. Selected
area is 39.78 square metres, about 2.03 percent of the entire closed mesh area
(including its large artificial faces). Strongly selected area is 19.40 square
metres. The calm-region and cut-edge selection checks both give zero. The side
and top plots were inspected before scene preparation.

The initially proposed independent VDB band was not validated. It targeted a
10 cm inner layer, nominal 3.125 cm voxels and peak extinction 20 per metre.
The first export was empty; later exports contained a grid but the roundtrip
sample test had a maximum error of 5 per cell unit and approximately 0.615
relative L1 error. Direct array equality and a nonempty VDB were insufficient
checks. Sparse volume sampling/refinement semantics remain unresolved; this is
not evidence of a Bifrost defect. Several preparation variants were retained,
including one process stopped after its memory rose to about 30 GB. No VDB
candidate was rendered and none is an accepted preset. The prototype now
requires an explicit --experimental flag. Its nominal voxel size must not be
reported as verified effective resolution.

## The alternative material and the failed first B

The agent informed Xue of a change to stock Arnold randomwalk subsurface
scattering (SSS). This is a different render-only approximation, not the promised
independent 10 cm layer and not an exact reproduction of the paper. Its 5 cm
radius scale is an initial appearance choice. The paper supplies no numeric
curvature thresholds, coefficients or layer thickness for a direct match.

The checked scalar is stored as a vertex color set, exported as varying RGBA.
A mix shader blends the original water with neutral Standard Surface SSS, with
subsurface weight one, color 0.98, radius (1,1,1), scale 0.2 scene units = 5 cm,
and anisotropy 0.2. The dielectric specular settings remain IOR 1.33, weight one
and roughness 0.03 in both components. No procedural texture or new geometry is
added. The original uniform dome, warm directional light, ACES display transform
and 15 m orthographic side camera are retained. This does not create advected,
aging surface foam or new spray.

At commit 9b5dfdd, A and B rendered successfully through licensed hidden GUI Maya
Render Sequence. A took 12.49 s inside the render call and 23.36 s including the
process; B took 12.79 s and 21.50 s respectively. Both are 960 by 540, 54-channel
float EXRs with a valid Arnold license. Their AOV sum closes to 8.89e-8 relative
L1; PNGs match the independently applied ACES display transform within 0.105
mean code values. Cache and frozen mesh hashes are unchanged, and volume AOVs
are zero as required by mist-off.

However, B's SSS AOV is exactly zero and its raw beauty equals A exactly. The
experiment therefore failed to exercise scattering and says nothing about its
visual benefit. Maya's multDoubleLinear did not survive as a shading network:
the rendered Arnold node inventory contains neither the mask reader nor its
multiply node. The successful renderer exit is not a successful comparison.

The repaired scene uses Arnold's native aiMultiply. A non-rendering ASS export
now contains interfaceMaskReader, its interfaceMask attribute, the native
multiply, and the link into interfaceMix. The repaired control scene has the
same geometry, optics, lighting and camera as the validated A; only the enabling
node's storage changes from scalar to RGB. The old A can be reused with its
original receipt and image hashes recorded explicitly. The original failed B
and its receipts are preserved. A separate image_check.json marks its actual
failure rather than rewriting the successful process receipt.

The repair was prepared and export-checked before requesting permission for one
additional B render, because that exceeds the originally agreed image count.
The requested allowance is 2–5 minutes including checks, with no simulation.
At the time of this record, that extra render has not started. The exported
connection check is not yet evidence that the SSS contribution renders correctly.

## Artifacts and sources

The mask and review are in out/runs/breaker_interface_scatter_20261019/.
The same prefix through _v9 retains failed/unvalidated VDB preparation variants.
The two first images are in out/runs/breaker_interface_sss_20261019_v3/;
the repaired, unrendered scenes and ASS checks are in the _v4 folder. Earlier
SSS folders retain preparation failures. No production defaults changed.

The scripts are breaker_interface_mask.py, the explicitly unvalidated
breaker_interface_volume.py, breaker_interface_render.py and
analyze_breaker_interface.py under experiments/. Stock public Maya, Bifrost
and Arnold interfaces are used; the Arnold renderer itself is not modified.

[PF-FLIP Section 7, page 15](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf#page=15)
describes the reference pipeline. [Arnold Standard Surface SSS](https://help.autodesk.com/cloudhelp/ENU/AR-Core/files/ac-shading/ac-surface-shaders/ac-standard-surface/arnold_user_guide_ac_standard_surface_ac_standard_subsurface_html.html)
documents randomwalk and the radius scale. [Bifrost Foam](https://help.autodesk.com/cloudhelp/2024/ENU/Maya-Bifrost/files/GUID-B8ABEF0F-4B49-4A7A-8EE5-86A596102C75.htm)
is a separate particle-simulation feature covering bubbles, foam and spray; it
is not the static scattering treatment used here, nor automatically compatible
with this custom solver's cache.
