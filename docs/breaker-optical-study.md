# The breaker's optical layers and PNG export

On 2026-10-19, two authorized 960 x 540 images reused the existing 24.9 s
cache's baked control mesh. No solver, surface reconstruction, spray, foam,
material change or new light was run. A kept specular ray depth 1; B changed
only that depth to 3, following the [official ocean tutorial](https://help.autodesk.com/cloudhelp/2025/ENU/AR-Maya/files/am-Arnold_for_Maya_User_Guide/tutorials/am-Shading/am-Displacement_Tutorials/arnold_for_maya_shading_am_Rendering_an_Ocean_with_Displacement_html.html).
The static scene came from the geometry study in
[PR #38](https://github.com/cuix-github/Anisim/pull/38); its binary and all images
remain local. Both optical renders used clean commit `7cc63d5`.

The strongest finding was a display-export error. The PNG driver retained
`Use Output Transform`, while `outputTransformEnabled` was false. Its RGB
values match the denoised linear ACEScg EXR multiplied by 255 and clipped:
mean absolute error 0.282 of an 8-bit level. They do not match the configured
ACES 1.0 SDR-video (sRGB) view (23.59 levels of difference). The PNG has no
ICC profile or transfer-function tag that would make an ordinary image viewer
interpret those values as scene-linear ACEScg. The renderer's textual
`arnold/color_space` metadata does not perform a display conversion.

The same EXR, converted using the OCIO library and configuration shipped with
Maya, reveals substantially more midtone structure. It still looks like
smooth, glassy water; correcting the display does not restore missing shapes
or supply whitewater. This finding applies to this pipeline and its reproduced
control, not every earlier image that has not been examined.

`render_surface_preview.py` now explicitly selects `Use View Transform` for
PNG output. The installed driver enum was checked without rendering:
`Raw`, `Use View Transform`, `Use Output Transform`; the selected mode read
back correctly. The offline EXR conversion was visually and numerically
inspected. A new production render through that changed PNG path was not
run. The mode-query process returned 0 but printed the installed MtoA's
shutdown stack trace after its successful query; this is not a claim of a
clean plugin shutdown. The diagnostic runner deliberately preserves the old
PNG output for comparison and writes a separate raw, full-float ACEScg EXR.
See [Arnold output color management](https://help.autodesk.com/cloudhelp/ENU/AR-Maya/files/am-Arnold_for_Maya_User_Guide/render-settings/arnold_for_maya_render_settings_am_Common_Output_html.html).

## What the controlled comparison establishes

Before rendering, conversion limits were fixed at a mean difference of one
8-bit RGB level and at most one percent of pixels changing by more than eight
levels, both for the full image and the curl rectangle (200,175)-(800,530).
A passed: full-image mean 0.0435 levels; curl mean 0.0897 levels, with 0.0441
percent of curl pixels above eight levels. This is close agreement, not
bitwise equivalence; the maximum difference is 36 levels at sparse pixels.

Recorded scene state differs only in `GISpecularDepth` after excluding the
documented per-process `lightBlindData` address. Position and normal AOVs are
bit-identical between A and B. Geometry, camera, materials, lights, cache and
baked-scene hashes are unchanged. Transmission depth stays 8; total depth
stays 10. The linear component sum reconstructs raw beauty with relative L1
error about 9.1e-8, independently of denoising.

| Curl-rectangle measurement | A: depth 1 | B: depth 3 |
| --- | ---: | ---: |
| Specular share of raw linear RGB sum | 33.99% | 32.41% |
| Transmission share | 61.93% | 63.75% |
| Diffuse share | 4.08% | 3.84% |
| Raw total RGB sum, relative to A | 100% | 106.19% |
| Arnold render wall | 12.90 s | 15.34 s |
| Render-call wall including translation/output | 13.99 s | 16.28 s |
| Process wall including startup and failed collection | 25.65 s | 27.95 s |

These shares describe a rectangular image region, including non-water pixels,
not material reflectance or perceptual luminance. The AOV category describes
the first scattering event: the transmission layer can contain later internal
reflections. Most of the additional RGB contribution appears in transmission
(5.77% of A's total), with 0.424% in indirect specular. Direct specular is
unchanged. Thus an apparently mirror-like patch need not belong to the first
reflection layer. The 1-to-3 change improves some internal paths but does not
produce the desired realism; it is not a convergence test at unlimited depth.

The raw and denoised images both retain the large smooth shapes. At this
preview resolution, denoising is not supported as their principal cause.
Normals show those broad sheets too. The remaining questions concern the
surface representation, the structured environment and scale-consistent
internal optics; the current geometry and whitewater limitations remain.
No additional candidate look or simulation was launched.

## Artifacts, failures and reproduction

Canonical outputs are `out/runs/breaker_optical_aovs_20261019_v3/`:
`A/` and `B/` hold the scenes, PNGs, 54-channel EXRs, states and logs;
`analysis.json` contains the checks. `color_comparison.png` isolates the
display transform on A. `comparison.png` compares A/B through the same ACES
view. `aov_panels.png` shows B's components at the same exposure; the normal
visualization maps [-1,1] to [0,1]. `*_view.png` are offline display conversions,
not extra renders. PFM files retain extracted floating-point channels.

Two early setup probes failed on metadata attribute queries before rendering;
an intervening inspection passed. They remain in the unsuffixed and `_v2`
folders. Both actual renders completed and wrote their PNG/EXR, but the
collector expected `tmp/*.png` while MtoA wrote
`tmp/beauty/beauty.0001.png`. Their real process exits were 1, not 0. The files
were validated and collected offline without rendering again; original
failures remain in `failed_collection_receipt.json` and `worker_result.json`.
The current receipts explicitly say `passed_after_offline_collection_recovery`
and retain exit 1. The collector now searches the nested directory.

The shipped `oiiotool` lacks OCIO support; channel extraction uses it, while
`convert_breaker_aov_views.py` uses Maya's Python OCIO binding without starting
Maya or rendering. The first display-equivalence assertion failed and led to
the PNG finding; its threshold was not loosened. The verifier instead tests
the explicit linear-export explanation against the saved pixels and settings.
Preparation/debugging and offline inspection were additional to the process
times above; a total-time overrun was reported during the work.

To audit saved outputs without rendering:

```
.venv\Scripts\python.exe experiments\analyze_breaker_optics.py
```

To run a separately authorized new comparison, provide a fresh output directory:

```
.venv\Scripts\python.exe experiments\breaker_optical_aovs.py A --out out/runs/<new-tag>
.venv\Scripts\python.exe experiments\breaker_optical_aovs.py B --out out/runs/<new-tag>
```

The runner requires a clean checkout, verifies the frozen input hashes, refuses
existing arm directories or another Maya process, bounds each owned Maya session
to six minutes, and prevents B until A passes. The small committed numeric
record is `results/rendering/breaker_optical_20261019.json`.
