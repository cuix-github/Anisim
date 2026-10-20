# The 24.9 s breaker: surface inspection, 2026-10-19

Five previews isolate the rendering choices on the existing `wavesea4` breaking
section. No simulation, cache extraction, particle editing or new spray emission
was run. The source is `out/runs/wavesea4/cache/breaking_section`: 811,717 particles,
one frame at 24.9 s, 0.25 m per cell, with the existing radius classification kept.
All five source files have matching SHA-256 hashes before and after every render.
The numeric record is `results/rendering/breaker_20261019.json`; the images and
scenes stay local under `out/runs/breaker_render_study_20261019/`.

The original scene and control A have identical saved nodes and connections
except for the camera. A was rendered and inspected before adding the controls.
Each arm started from a clean recorded Git revision, cleared inherited `SURF_*`
settings, and used a new output directory. B through E were rendered serially.
The saved Maya scenes were compared independently of the environment settings,
including their serialized Bifrost graph inputs. The verifier also checks the
actual Arnold image size, AA samples and water triangle counts.

| Preview | Comparison | Complete wall, s | Frame-render call, s | Water triangles |
| --- | --- | ---: | ---: | ---: |
| A | `paper_clean` control | 610.751 | 62.7 | 274,090 |
| B | A, curvature foam bypassed | 567.264 | 41.5 | 274,090 |
| C | B, mist shape hidden | 291.153 | 15.1 | 274,090 |
| D | C, masked MCF iterations 30 to 10 | 290.756 | 16.4 | 292,516 |
| E | C, environment RGB multiplied by four | 292.361 | 14.8 | 274,090 |

All previews are 960 x 540, CPU Arnold at AA 4 with the existing OIDN denoiser.
Their fixed camera, in cells, is `596 76 188 570 58 256` (eye then target).
The particle sphere radius and surface voxel size both remain 0.5 cells,
or 0.125 m. The sun, water absorption, IOR, roughness and camera exposure stay
fixed. These are individual elapsed times on a shared machine, not a repeated
performance benchmark. The old original's 158.7 s was only its frame-render
call; its completion marker records 648 s for scene evaluation plus rendering.
Most of the preview wall is before the frame-render timer, in scene evaluation.

Removing curvature foam takes away the hard white patches, but leaves the thick
water silhouettes. Hiding mist reveals the holes and hanging sheets; it does
not create more geometry. Reducing MCF retains more relief and also more small
surface ripples. Four times the environment light makes the shapes easier to
read but gives strong blue highlights and retains a metallic/glassy appearance.
None of these previews establishes reference-level realism.

The candidate preset, `breaker_inspect`, combines D's 10 MCF iterations with
twice the original environment RGB (`0.06 0.07 0.09`), an intermediate choice
between C and E rather than one of the five measured arms. Curvature foam and
mist remain off. Its 1920 x 1080 presentation camera is
`624 76 248 574 58 274`: eye (156, 19, 62) m, target (143.5, 14.5, 68.5) m.
The eye is inside the cache's z = 60 to 78 m band, facing the wave instead of
the near artificial cut face. This changes framing as well as the selected
settings, so the candidate is not a single-variable comparison against A.

The candidate completed in 361.589 s, including an 84.8 s frame-render call,
with the same 292,516 water triangles as D. It exposes the folded sheets and
holes more directly, but the glassy/metallic appearance remains conspicuous.
It does not achieve the desired large realism improvement and is not promoted
to the default look. The six renders total 2,413.874 s (40.23 minutes) including
their launches and scene evaluations; setup, inspection and recording are extra.

Hiding the mist also hides the escaped particles that the current pipeline
represents through that volume, including the 1,192 positive-radius particles
in this cache. It does not return them to the surface. This preset is therefore
for inspecting the bulk and sheets, not judging spray capture or final
whitewater. No lower claim about the solver follows from their absence here.
The finite cache's artificial boundary and the particle-to-surface and
particle-to-mist representations remain important limits to investigate on
this frozen state before deciding that another simulation is needed.

To reproduce, invoke `.venv/Scripts/python.exe experiments/breaker_render_study.py`
once for each of `A`, `B`, `C`, `D`, `E`, inspecting each before proceeding.
Use a fresh `--out` directory; existing arm directories are never overwritten.
The candidate command is `experiments/breaker_render_study.py D --candidate
--preset breaker_inspect --camera "624 76 248 574 58 274"` with the same `--out`.
The candidate uses the named preset without the arm's environment overrides,
so its saved scene tests the preset itself. Run
`experiments/verify_breaker_render_study.py --out <study-directory>` afterward.
These are rendering commands, not a request to launch the study again.
