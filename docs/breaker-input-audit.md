# What reaches the breaker render

On 2026-10-19, the particle playback was reported to show much more structure than
the water render. The follow-up inspected the existing 24.9 s cache, source graphs,
saved Maya scenes and Arnold logs, and compared real wave photographs with the
paper reference. No Maya scene was evaluated, no simulation or render was started,
and no production defaults were changed.

`experiments/audit_breaker_render_inputs.py` reads the existing inputs and writes a
numeric audit and two side-view particle plots. The verified run is in
`out/runs/breaker_input_audit_20261019_verified/`; the numeric receipt is also in
`results/rendering/breaker_inputs_20261019.json`. It checks the five source-file
hashes against every original render receipt and again after reading. Scene hashes
also match. The audit took 0.40 s, excluding Python startup, from a clean checkout.
An earlier local output directory is retained but superseded: its visibility parser
missed Maya's `setAttr -k off` flag. The verified version asserts visibility and
material connections for all six saved scenes.

The cache contains 811,717 particles. The surface reader admits 810,525 whose radius
channel is below 1e-9. It excludes 1,192 particles. These are sent to mist and are
hidden when mist visibility is off. This is a display decision applied to existing
data, distinct from generating new spray in the solver. It explains some missing
flying water, but does not explain the entire bulk surface's glass-like appearance.
The radius channel is a classification and droplet-model quantity, not the geometric
radius of every simulated water parcel.

The particle plot separates an 18 m-deep projection from a 0.5 m-thick slice at
z = 67.75 to 68.25 m. The slice contains 22,640 particles, including 36 excluded by
radius. The curl, thin regions and disconnected points are present. Plot marker sizes
are deliberately for legibility, not physical droplet radii. Neither plot measures
which unflagged particles survive surface reconstruction. A production mesh for this
frame was not saved as explicit vertices, so that loss remains unmeasured.

The graph's comments say isolated particles can fail to produce a surface at the
default radius/voxel ratio. That is a hypothesis to test separately from filtering:
the current small gate's explicitly marked droplet is excluded by the reader before
meshing. Its zero nearby vertices do not independently prove a voxelization failure.
Other legacy probe arrays lack a complete recorded configuration and should not be
used as a controlled comparison for this frame.

There are different resolutions in the pipeline:

| Stage | This cache / render | Reference evidence |
| --- | --- | --- |
| Pressure grid | 0.25 m cells, 688 x 96 x 320 | About 0.326 m nominal cells from the displayed 1 km / 3072 |
| Water surface reconstruction | 0.125 m configured voxels; 0.125 m particle sphere radii | PF-FLIP Section 7 describes skinning and smoothing, without specifying a matching voxel size |
| Spray-density grid | 0.125 m configured voxels | Table 7, Figures 4/5/25: 12288 x 6144 x 4096, four times the pressure resolution; about 0.081 m nominal spacing |
| Screen at the camera target | About 6.22 pixels per surface voxel in the 960-wide previews; 15.78 in the 1920-wide candidate | Different framing; no matched screen-space comparison has been made |

The screen estimates assume a 35 mm lens, 36 mm horizontal aperture and horizontal
film fit: pixels = image width x 35 / 36 x voxel size / target distance. Distances are
18.748 m for A-E and 14.790 m for the candidate. These are projected lengths, not a
claim that the renderer produces square blocks of that size. Smaller pressure cells
alone do not establish finer visible detail, equivalent particle dynamics or an
equivalent droplet-density representation.

The saved water material has roughness 0.03 in every arm. The Arnold logs give
diffuse depth 1, specular depth 1, transmission depth 8 and indirect-volume depth 3.
The environment is a uniform dome plus a directional sun and dark ground. The prior
brightness arm changed the dome's magnitude, not its spatial structure. No deliberate
micro-normal texture was added. Neither roughness, a structured reflection environment,
nor specular depth was tested by the five-preview study. An IOR of 1.33 is already
assigned; appearance alone does not establish an incorrect IOR.

Two photographed close waves provide useful counterexamples to simply making water
matte: [Tim Marshall's wave](https://unsplash.com/photos/ocean-tunnel-wave-uanoYn1AmPs)
combines smooth dark wall regions with dense small highlights and broken white water;
[Jon Wright's The Curl](https://www.jonwrightphoto.com/products/the-curl-ocean-art-collection)
has a transparent, reflective wall with fine relief, thin edges and scattered drops.
Both were visually inspected. They are qualitative references, not measurements of
our solver or calibrated lighting targets. The source images were not copied into
the repository.

The [PF-FLIP paper](https://ge.in.tum.de/download/Adaptive_Phase_Field_FLIP_preprint.pdf),
Section 7 and Table 7, distinguishes surface reconstruction, curvature-dependent
foam-like scattering and finer spray-density rendering. It also includes additional
spray emission and droplet motion under air drag. Its statement that no procedural
textures were applied means it would be wrong to claim its fine detail comes from a
noise map. Our viewport uses speed-dependent particle colours; that diagnostic
contrast cannot be equated with a physically shaded water surface.

The next bounded experiment should separate geometry from shading before refining
the entire reconstruction grid. A proposed three-image comparison keeps the existing
cache, A-E camera, 960 x 540 resolution, 10 MCF iterations, and foam/mist visibility
off: (1) reproduce D, (2) assign a neutral matte surface to the identical geometry,
(3) use the original water shader with only roughness changed from 0.03 to 0.12.
The opt-in builder controls `SURF_DIAGNOSTIC_MATTE=1` and
`SURF_WATER_ROUGHNESS=0.12` prepare those two changes. Their defaults preserve the
existing assignment and roughness; they have been syntax-checked but not evaluated
in Maya or visually validated.
The matte image asks what relief and discontinuities actually reach the renderer;
the roughness arm asks how much of the coarse appearance is sharp reflection. It
cannot establish that roughness alone is the final ocean look. All three need a
fixed camera and unchanged cache hashes. Estimated end-to-end work is 20-30 minutes
based on the previous roughly five-minute full launches, with notice before extending
past 30 minutes. This proposal has not been run or approved as a new render batch.

If geometry remains suspect, measure particle-to-surface coverage and normal detail
on a bounded, padded curl patch before changing reconstruction resolution. Halving
voxel spacing can multiply dense-grid storage by eight; the previous full scene
already incurred substantial memory pressure. Any real-frame reconstruction evaluation
also needs a separately agreed time budget. No new simulation is justified by this
inspection alone.
