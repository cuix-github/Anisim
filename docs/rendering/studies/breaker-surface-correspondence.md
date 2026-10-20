# Cached breaker surface correspondence

On October 19, 2026, Xue authorized investigating reconstruction and interface
scattering together, while deferring spray changes and simulation. The experiment
in `experiments/breaker_surface_correspondence.py` uses public stock Bifrost 3.1.0.8
nodes and the unchanged 24.9 s cache. It launches no Maya or Arnold rendering.

The authored baseline has **zero Gaussian passes**, ten masked mean-curvature
passes, sphere radius 0.5 cell, surface voxel size 0.5 cell and mesher smoothing
0.1. An earlier correspondence note mistakenly included a Gaussian pass by
reading the generator default instead of the saved scene override. The simulation
cell is 0.25 m. ST-FLIP's density evolution and our masked signed-distance evolution
remain different; matching iteration counts would not establish equivalence.

The instrument was committed before execution. Its initial synthetic smoke test
failed compilation because NumPy writing does not support the mesh's unsigned
topology arrays. Explicit integer conversion fixed that interface; the second
smoke test passed. Both attempts remain in the ignored run directory.

The full control then reproduced all 146,326 vertex positions exactly and all
292,516 triangles, allowing for face order/winding through Maya conversion. It
took 258.86 s including process startup. Setting the mesher's detail scale to 0.5
and enabling subdivision produced 332,491 vertices and 664,846 triangles in
237.47 s. These are job times, not a performance comparison. The sampled input
field was bit-identical, but the circled floating diamond remained six vertices
and eight triangles. This mesher change is not a demonstrated repair.

Curvature was recomputed and its sphere/plane analytic gates passed. Strong mask
coverage changed from 19.40 to 19.82 square metres; coverage above 0.05 changed
from 39.78 to 39.58 square metres. Crop-edge coverage remained zero. The legacy
four-pass graph filter has a different physical width on a refined mesh, so these
numbers are a support audit, not a perfectly matched optical comparison. The
floating fragment remains fully selected by the mask. Small water geometry can
have high curvature without supplying evidence of entrained air; curvature alone
does not resolve the white-coated-fragment appearance.

A follow-up crop retained particles within four metres of the fragment. Before
changing anything, its control reproduced the six fragment vertices exactly.
Its input-field slice differs from the full control by at most 2.38e-7 cell.
Using stock `resample_volume` at detail scale 0.5 before meshing then produced a
nearby connected component with 86 vertices and 168 triangles, versus 6/8. This
is more than tessellation: its x extent grew from 0.171 to 0.432 m and its y
extent from 0.113 to 0.251 m. Both meshes are closed/manifold, but spatial
correspondence is not proof of exact particle provenance, and resampling may
change which neighboring structures connect. The component's mask still
saturates to one. The crop jobs each took about 70 s.

Therefore neither candidate replaces the existing render geometry. The local
result identifies volume sampling/contouring as a worthwhile investigation, but
does not establish that a larger, rounder shape is more faithful. The sampled
linear-SDF zero contour and the extracted mesh also differ substantially near
this small component; their interpolation, smoothing and contour conventions
need a controlled correspondence check before using that contour as ground truth.
No solver failure or need to rerun the simulation follows from these results.

During this work Xue first saw the valid thin-foam A/B images on desktop. He finds
both promising and wants different distances, angles and later cached frames
before choosing between them. Preserve both existing appearances for that viewing
study. The procedural B remains an explicitly labelled look experiment, not a
claim that the papers use procedural texture. No new classification or droplet
replacement is being introduced.

Artifacts are under `out/runs/breaker_surface_correspondence_20261019/`: four
reconstruction arms, synthetic checks, graphs, logs, source hashes, mesh/mask
arrays, `comparison.png`, `local_comparison.png`, and numeric analysis. Small
receipts are committed in `results/rendering/breaker_surface_correspondence_20261019.json`.
All source cache hashes remained unchanged.

The method references remain [ST-FLIP Section 4 and Appendix B](https://ge.in.tum.de/download/ST-FLIP.pdf)
and [PF-FLIP Section 7](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf).
Stock-node port behavior was checked against the shipped public documentation
for `volume_to_mesh`, `resample_volume` and `sample_volume` in Bifrost 3.1.0.8.
