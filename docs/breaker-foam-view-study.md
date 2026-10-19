# Paired foam views and the water-crop boundary

On October 19, 2026, Xue found both the uniform thin surface and clustered foam
promising after seeing the images on desktop. He requested different angles,
distances and cached times after 24.9 s before choosing. Eight 960 x 540 images
were completed by 22:46 UTC, within the announced 40–60 minute preparation,
rendering and inspection allowance. A is uniform thin foam; B redistributes its
integrated extinction into the existing procedural clusters. B remains a look
experiment, not a claim that the papers use procedural textures.

All images retain corrected-baseline optics, lighting and ACES display conversion.
Mist is off. No solver, spray emission, new particle classifier or renderer-source
change was used. The experimental geometry from the
[surface correspondence audit](breaker-surface-correspondence.md) was not adopted.

| Cached time and view | Orthographic width | A render call | B render call |
| --- | --- | --- | --- |
| 24.9 s, wider side | 24 m | 75.68 s | 46.52 s |
| 24.9 s, elevated oblique | 18 m | 247.44 s | 146.10 s |
| 26.025 s, side | 18 m | 35.61 s | 22.38 s |
| 28.025 s, side | 18 m | 8.28 s | 7.29 s |

These eight calls total 589.29 s inside rendering, excluding preparation and
process startup. Each A passed scattering and image checks before its B launched.
All eight licensed GUI-Maya EXRs have nonzero volume contribution, zero SSS and
component closure below 1e-7 relative L1. Independent ACES display conversion
agrees within 0.23 mean 8-bit level. Within each pair the recorded state differs
only in foam arm and foam asset hash. Appearance was also inspected visually.

The first four images reuse the exact validated 24.9 s mesh and foam assets.
Later frames use overview cache frames 4 and 6, retaining x >= 120 m and
60 <= z < 78 m. This crop rule reproduces the original 24.9 s particle selection
exactly when applied to its corresponding cache. The existing radius labels
still govern the surface/mist split; no new classification was introduced.
The later crops contain 815,247 and 782,268 particles. Stock reconstruction uses
the frozen control's sphere radius, volume resolution, disabled Gaussian pass,
ten masked MCF iterations and meshing settings. It took 244.67 and 237.38 s
including bifcmd startup, without simulating fluid motion.

A control first reproduced all 146,326 original vertex positions and the
292,516 triangles, allowing face order/winding differences. Attempting to export
`point_normal` returned a single default zero vector because that property was
absent; the normal gate stopped before any later-frame render. Maya's normal
construction on the identical mesh then matched the frozen converted mesh's
corner normals to 1.79e-7. This validated path is used for later frames, with
corner-normal readback, exact point checks and preserved mesh render flags.
The failed property export and receipt remain beside the successful control;
the reconstruction was reused rather than rerun.

One earlier owned GUI launch also failed before scene loading because `runpy`
did not add the script import directory. The launcher was repaired, its failed
artifacts retained, and no Arnold render call was spent on that launch.

For each later frame, side and top mask plots, fine density plots, layer thickness
and crop exclusions were checked before scene preparation. Expanded Arnold
procedurals verify actual imported positions, radii and transforms: primitives
stay within 3 cm of their source surface and more than 0.58 m from a crop edge.
The same curvature and height gates were retained. Selected surface area above
mask 0.05 falls from 31.90 m2 at 26.025 s to 3.89 m2 at 28.025 s. The reduced
white coverage is an instantaneous-mask result; it is not simulated foam decay.
This pass has no transported foam history, and these stills do not validate
temporal continuity for a movie.

During inspection Xue said the white crest was approaching the target, while the
whole water body still felt like a one-metre fish tank instead of a sea tens of
metres across. The position and normal AOVs exposed a concrete presentation
problem: the camera sees the artificial side of the cropped water mesh.

`breaker_crop_audit.py` reads existing AOVs, without rendering. Among transmitted
samples inside the water-crop bounds, 54.2% in the 24.9 s wider side view and
61.8% in the 26.025 s side view lie within the front z < 60.25 m band and have
abs(Nz) > 0.7. The oblique view's share is 12.6%. The diagnostic map places these
samples on the conspicuous glassy wall. These are thresholded shading samples,
not exact object-ID coverage or percentages of the whole image. The measurements
support a substantial crop-side contribution, not an explanation of every
glassy reflection or every scale cue. Arnold documents P and N as world-space
shading position and smooth normal in its
[AOV reference](https://help.autodesk.com/cloudhelp/ENU/AR-Core/files/ac-output-aovs/arnold_user_guide_ac_output_aovs_ac_aovs_html.html).

There is also a framing mismatch. A 24 m view at 960 pixels spans approximately
2.5 cm per pixel, so a 25 cm pressure cell projects to about ten pixels before
orientation effects. Taking the user's reference screenshot's 1000 m scale and
2048-pixel width literally gives about 49 cm per pixel: its roughly 33 cm cell
projects to less than one pixel. This is a rough sampling comparison, not a claim
that pressure cells are rendered directly or that the images have identical
projections. A smaller absolute cell size does not imply finer-looking detail
when the camera magnifies it much more strongly.

The configured domain is 172 x 24 x 80 m. Existing overview particles span roughly
72 x 80 m horizontally; the current reconstruction is only an 18 m-wide strip,
and these views cover 18–24 m. A fuller cached-water view with artificial sides
outside the important sight lines should be evaluated before another material
sweep. Its geometry preparation and memory requirements need an estimate first;
no expanded reconstruction or extra render was started here. Orthographic
projection and the large visible floor may also weaken scale cues, but neither
has been isolated as a cause, and the reference also has a black background.

Large-ocean production work illustrates why surface detail and simulated particle
count are separate concerns: the public
[Moana pipeline description](https://history.siggraph.org/wp-content/uploads/2022/09/2017-Talks-Palmer_The-Ocean-and-Water-Pipeline-of-Disneys-Moana.pdf)
separates wave frequency bands and adds high-frequency detail during rendering.
That is context for diagnosis, not evidence that PF-FLIP/ST-FLIP used that method,
nor authorization to add procedural waves to this reproduction.

Artifacts live under `out/runs/breaker_foam_views_20261019/`. Each view has an A/B
pair and full-size `comparison.png`; the root contact sheet shows all eight.
`frames/` keeps cached crops, reconstruction and normal controls, masks and foam
checks. `crop_audit.png` and its JSON record the boundary evidence. Numeric
receipts are committed in `results/rendering/breaker_foam_views_20261019.json`.
All original source-cache hashes remain unchanged. No final ocean preset is
selected, and the next expanded-water render remains a proposal.
