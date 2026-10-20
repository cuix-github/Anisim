# Neutral water, finer cached mist and a closer side camera, 2026-10-19

This follow-up uses the existing 24.9 s crop and frozen surface. Xue preferred
its corrected baseline color to the previous saturated Deep Water candidate,
wanted a closer side camera, and explicitly authorized trying render-only mist
refinements. The agreed scope is two 960 x 540 previews, with no solver run.
The original view changes a combined look package; it does not isolate causation.
The second view changes only the camera. Curvature foam remains bypassed.

The water keeps the preset's dielectric weights, IOR 1.333, roughness 0.12 and
forward-scattering anisotropy 0.75, but restores baseline transmission RGB
(0.3, 0.55, 0.6) and optical depth 40 cell units = 10 m. Scattering becomes
neutral (0.15, 0.15, 0.15), an artistic candidate, not measured sea coefficients.
Physical Sky keeps elevation 30, azimuth 160, turbidity 6 and sun size 5 degrees.
Dome intensity is 0.44 / 0.171667 = 2.5631, matching the baseline dome's mean
RGB energy to the previous Physical Sky map's measured average. This does not
match total illumination: the former separate directional sun is disabled.

Mist uses existing cached positions through the existing particle_mist compound.
Selection includes positive cached radii plus particles outside the reconstructed
classification surface by its threshold; 1,192 positive-radius particles is not
the complete selected count. The classifier uses the existing MCF 10 setting,
but the visible water mesh is reused unchanged. Computing this classification
field and rasterizing mist is render preparation, not a pressure solve.
Mist voxel spacing is 0.25 cell = 6.25 cm instead of 12.5 cm; cloud radius is
0.5 cell = 12.5 cm instead of 25 cm. Density fraction falls from 0.001 to
0.00025. The existing kernel weight scales as density / radius^3, so reducing
density fourfold approximately preserves the central optical depth of a cloud
whose radius halves. It reduces total represented spray water fourfold; this is
an appearance choice, not a mass-conserving refinement or added simulation detail.
Standard Volume uses white scattering, anisotropy 0.4 and density scale 10000,
which compensates the existing density-field normalization.

The close orthographic side camera keeps zero pitch/roll and yaw 168 degrees,
normal to the dominant 12 degree wave train. It frames 15 m rather than 45.76 m,
a 3.05-fold enlargement, centered at (145.5, 13.5, 70) m around the curling lip.
All static mesh points, normals and topology must remain identical in both views.

## Fifth-chunk feasibility: read only

The fourth segment's full checkpoint at 24.0 s still exists under
out/runs/wavesea4/dumps/c3. All 13 required arrays have consistent 61,280,955
particle counts and float32 storage; combined size is 4,167,105,980 bytes.
Radius and extra-spray flags are both entirely zero. Header and field inspection
is not a large-checkpoint continuation test, and none was launched here.

The fifth segment covers 24-30 s. Its original receipt records 240 steps in
14,235.8 s (3 h 57 min), 59.32 s/step. Reaching only 24.9 s from the checkpoint
requires 36 steps; multiplying by the old average gives 2,135.52 s, about
35.6 minutes, excluding startup and added model cost. This is an extrapolation,
not a benchmark of the proposed configuration. It supplies just 0.9 s of new
spray history and is not equivalent to running the model from time zero.

The solver carries droplet radii and spray flags through full-state resume;
piston_beach already has the --from, --droplets, --integrator and --spray pattern.
wave_basin does not expose those controls or a new output tag. Its cache writer
interpolates stable particle indices and manufactures a fixed render radius from
escape flags. It needs an export path for changing counts/identities and actual
radii before new emission can be used correctly. For a 24.9 s still, an exact dump
can avoid temporal interpolation. Any future run must write a separate tag and
preserve the original c4. Enabling the droplet model also changes escaped liquid
motion/deletion, so it must not be presented as guaranteed bit-identical bulk flow.

Rendering existing mist is a visibility/reconstruction choice. Generating new
spray and changing its motion is a solver choice. Neither a fine volume nor a
later start can invent the missing earlier spray history.

## Execution record

The first preparation failed before rendering: importing build_surface_scene
following Maya initialization reset BIFROST_LIB_CONFIG_FILES and removed the
shipped libraries that Maya had prepended. Importing the helper before initialization
fixes that ordering. The failed preparation, logs and receipts remain under
out/runs/breaker_neutral_mist_20261019; successful preparation and all current
artifacts are under out/runs/breaker_neutral_mist_20261019_v2. No extra render was
spent on the failed preparation. The corrected render instrument was committed
at 1ad4251 before either image.

The original image's Arnold render completed, but PNG collection initially
failed because the images workspace rule was set only during scene preparation.
Maya does not save that runtime rule inside the binary scene; its GUI process
wrote the PNG under the default project instead. The failed exit-1 receipt is
retained. The exact path in maya.log, matching write timestamp, licensed EXR and
completed error-free Arnold log identify the finished output. It was copied into
the run directory with a separate recovery.json, without rendering again.
The second render instrument at 8d52784 sets the output workspace after opening
the scene and permits a hash-verified recovered first image. The candidate scene
and shading are unchanged by this collection-only correction. The two camera
states are compared independently rather than pretending the instruments share
one SHA.

## Results

Both images completed and passed offline verification. The original-view render
call took 318.15 s and its full Maya process 341.25 s; the close side view took
297.70 s and 314.97 s respectively. These render-call times include Bifrost mist
preparation and scene export, not only ray tracing. Preparation, inspection and
documentation are additional. A possible thirty-minute overrun was communicated
before the second image completed; no further rendering was added.

The original image is less deep blue but remains turquoise and glassy. The side
view makes the fragmented lip and isolated mist clouds much easier to inspect,
but also shows the crop's cut surface and broad transparent sheets. The finer
mist does not form the reference's continuous crest spray. This is not an
accepted replacement baseline or evidence that cached particles cannot support
any better appearance. Neither foam nor new emitted spray was added.

Both EXRs are valid licensed 960 x 540, 54-channel float outputs. Independent
ACEScg-to-sRGB display conversion matches the PNGs within one 8-bit level (mean
0.189 original, 0.133 side). Raw AOV sums close to relative L1 below 9e-8.
The separate volume AOV is nonzero: 6.68 percent and 19.77 percent of summed raw
RGB respectively. These are image-energy shares, not particle counts or coverage;
they confirm that mist was rendered, rather than hidden. Original-view mean
linear transmission RGB is (0.02861, 0.04380, 0.04303), compared with the prior
blue candidate's (0.00513, 0.01008, 0.04744), and baseline's
(0.03481, 0.04915, 0.05555). The cast changed without matching the baseline exactly.

The two scenes differ only in camera state. Visible mesh points, stored normals
and topology hashes match the frozen control (146,326 vertices, 292,516 triangles).
The original-view position AOV matches baseline exactly; normal AOV maximum
absolute difference is 1.79e-7. This confirms that the visible geometry was not
quietly altered as part of the appearance change.

The numeric receipt is results/rendering/breaker_neutral_mist_20261019.json.
The successful local run folder contains comparison.png, aov_panels.png,
analysis.json, chunk5_audit.json and each camera's scene, state, raw EXR, display
PNG, logs and receipts. Original failure receipts remain unmodified.
Offline verification can be repeated with:

```
.venv\Scripts\python.exe experiments\analyze_breaker_neutral_mist.py
```

An independently authorized future repeat uses breaker_neutral_mist.py prepare,
original and side with the same fresh --out path. The launcher refuses overwrites
and requires a clean checkout, hash-pinned cache/control scene, no other Maya
render session, and enough free memory. It bounds each owned GUI render to eight
minutes. No production look default was changed and no solver was started.
