# Outputs: what the runs leave on the machine

Everything a script writes lands under `out\` at the repository root, which git ignores. Nothing in it is in the repository; this page says what is there, what it is for, and which folders are canonical. `graphs\runs.py` is the code's single idea of where a run's files live; `runs.py report` prints what is on disk. Sizes are from 2026-10-15.

The 2026-10-19 [Deep Water candidate](breaker-deep-water-study.md) lives in
`out/runs/breaker_deep_water_20261019/`: `original/` and `side/` contain the two
960 by 540 previews, static scenes, raw component EXRs, hashes and successful
receipts. `comparison.png`, `aov_panels.png` and `analysis.json` record the
comparison and offline checks. The candidate was not accepted as an improvement.
It reuses the validated local scene from `breaker_optical_aovs_20261019_v3/B/`;
there is no new simulation cache or surface reconstruction.

The later [neutral water and finer mist study](breaker-neutral-mist-study.md)
lives in `out/runs/breaker_neutral_mist_20261019_v2/`, with two closer-inspection
previews and offline checks. The unversioned-suffix folder preserves the failed
first preparation. The successful folder also retains the original PNG-collection
failure plus its separate verified recovery; that image was not rerendered.
`chunk5_audit.json` records checkpoint availability and timing estimates only.

## One folder per run

The intermediate foam preview lives in `out/runs/breaker_foam_mix_20261019/`.
`C/beauty_view.png` is the single new 24.9 s, 18 m oblique image, using
25 percent uniform plus 75 percent clustered volume density. `C/candidate.mb`
uses `C/foam.ass`; retain both. `layer.npz`, `density_check.png`, alignment and
image checks, and `analysis.json` preserve the experiment. Existing A/B images
in `breaker_foam_views_20261019/oblique/` are its controls, not newly rendered.

The [thin-foam study](breaker-foam-layer-study.md) retains failed support checks
under out/runs/breaker_foam_layer_20261019/, scene-preparation checks under _v2,
the ineffective first A under _v3, and the ineffective explicit-density A under
_v4, including its failed image check and unit-mismatch audit. B never rendered;
that round stopped after the approved replacement A failed. The subsequently
approved unit repair stopped in _v5 on an unavailable Maya display bounding box
without rendering; _v6 contains the valid A/B images, expanded-geometry alignment
checks, comparison.png, comparison_small.png and analysis.json. The existing-mesh corner audit and focused
eight-triangle fragment plot are in out/runs/breaker_corner_audit_20261019/.

The [interface-scattering study](breaker-interface-study.md) retains its mask
in out/runs/breaker_interface_scatter_20261019/ and unvalidated volume attempts
through _v9. The initial SSS A and failed-effect B are in
out/runs/breaker_interface_sss_20261019_v3/. The _v4 folder holds repaired scenes,
ASS export checks, the approved replacement B, analysis.json and comparison.png.
The comparison reuses v3/A; v4/A was not rendered. B now has a validated SSS
contribution, but the white patches remain thick and the body glassy. These
artifacts are a working scattering comparison, not an accepted ocean look.

Since 2026-10-16 a run writes everything it produces under one folder named by its tag:

```
out\runs\<tag>\
    dumps\c<N>\        chunk N as the loop runner writes it: the arrays it started from,
                       the per-substep dumps (frame, vel, esc, scl, ph, diag), the resume state (final*.npy)
    cache\states\      the Maya play caches (liquid.####.npy, speed.####.npy, cache.json):
    cache\full\          one frame per chunk's end; every dumped frame; the 24 fps blend
    cache\blend24\       of the dumps
    inputs\            the colliders: solid_phi.npy, colliders.obj, the setup graph's scratch
    receipts\          c<N>.json, summary.json, field.json, cost.json
    renders\<set>\     one render set's frames, logs and mp4 (audit = the painter's frames
                       the run draws as it goes, and its preview video)
    scenes\            the Maya scenes the launcher builds (play_<mode>.ma)
    logs\
```

The dumps are the only copy of a run: losing them means simulating again. Every cache and render is derived from them and can be rebuilt: `collider_showcase.py play --tag <tag>` writes the states cache, `collider_showcase.py blend --tag <tag>` the 24 fps blend (hours for the 512 run), and the full cache is written by the run itself. A run named by its tag alone takes its shape (grid, chunking, substep, graph) from its receipts, so these commands need nothing else.

| To | Run |
| --- | --- |
| See every run, its parts and sizes, what is derived from what, and everything else under `out\` | `.venv\Scripts\python.exe graphs\runs.py report` |
| Play a run in Maya by name (writes the states cache and builds the scene if missing) | `maya\play_showcase_in_maya.bat blend boulders512` |
| See the renames that would bring the older runs into `out\runs\` (a plan; nothing moves) | `.venv\Scripts\python.exe graphs\runs.py migrate` |

## The runs made before 2026-10-16

They lie in the older layout, and `runs.py` finds them there (the new layout first, the older second), so nothing had to move for the code to change:

| Folder | What | Size |
| --- | --- | --- |
| `out\loop3_<tag>_c<N>` | One chunk of a 3D run: per-substep dumps (`frame`, `vel`, `esc`, `scl`, `ph`, `diag`), the resume state (`final*.npy`) | the 512 run: 457 GB over `boulders512_c0` to `_c5` |
| `out\collider_showcase\<tag>_play24` | The 24 fps blend cache of a run | 252 GB for `boulders512` (487 frames) |
| `out\collider_showcase\<tag>_play` | The per-dump play cache, superseded by the blend cache | 124 GB for `boulders512` |
| `out\collider_showcase\<tag>_states` | One frame per chunk end, for a quick look (the finale run's is `play_states`) | 3 GB |
| `out\collider_showcase\<tag>_c<N>.json`, `<tag>_summary.json`, `field_<tag>.json` | The run's receipts: wall per substep, iterations, peak speeds, penetration audit | small |
| `out\collider_showcase\<tag>_png`, `<tag>.mp4` | The painter's audit frames and the preview video | 0.1 GB |
| `out\collider_showcase\<tag>_render`, `<tag>_surface`, `<tag>_mist*` | Render sets: the playblast, the surface and mist tries | 0.1 GB |
| `out\collider_showcase\solid_phi_<tag>.npy`, `<tag>.obj`, `out\collider_<tag>_ab`, `_c` | The colliders as a signed-distance field, the boulders as a mesh, the setup graph's scratch | up to 0.1 GB |
| `out\play_<tag>_<mode>.ma` | The Maya scenes the launchers built | small |

`runs.py migrate` lists the renames that would bring each of these runs into `out\runs\<tag>\`; `--apply` does them. They are renames on one volume: nothing is copied, nothing deleted, and `runs\<tag>\migration.json` records every move. A Maya scene built against the old place is rebuilt by the launcher the next time it is opened. Not done; the owner's call.

## Leftovers that can go

Listed, not deleted; deleting a run's data is the owner's call.

- `out\loop3_dry512_c0` (8.6 GB) and `out\loop3_dry384_c0` (6.3 GB) with `out\collider_showcase\dry384_*`: pre-flight memory tests of 2026-10-15, inputs only.
- `out\collider_showcase\boulders512_play` (124 GB): the per-dump play cache, every frame of which the blend cache also holds, interpolated.
- `out\loop3_smoke_c0`, `loop3_smoke_c1`, `loop3_smoke_fused_c0`, `out\smoke_blend24`, `out\surface_test`: smoke-scene test runs.
- The loose logs and sweep scripts at the `out\` root.

## The gates and the studies

The paired foam viewing study is `out/runs/breaker_foam_views_20261019/`: `wide`
and `oblique` retain the 24.9 s assets; `later26` and `later28` use existing caches
at 26.025 and 28.025 s. Every view has A/B scenes, licensed renders, image checks
and a full-size pair; the root contact sheet contains eight previews. `frames/`
retains reconstruction, normal controls and pre-render mask/layer checks.
`crop_audit.png` identifies the visible artificial water-crop side using existing
position/normal AOVs, without new rendering. [Findings and limits](breaker-foam-view-study.md).

The frozen 24.9 s render study is `out/runs/breaker_render_study_20261019/`: `A` through `E` contain the five 960 x 540 previews, `candidate_D` the 1920 x 1080 inspection candidate, and `comparison.png` the labeled contact sheet. Each arm keeps its saved scene, complete launch log, Arnold/Maya logs and hashed receipt. `verified_study.json` holds the independent scene audit. The source cache at `out/runs/wavesea4/cache/breaking_section` was not rewritten. [Settings, findings and reproduction](breaker-render-study.md).

Each gate's and study's working folder stays at the `out\` root under its own name, as before: `out\loop3_<gate run>` for the loop gates (`resume_*`, `lsn_*`, `st3d`, ...), `out\step_in` and `out\step_out` and their kin for the node gates, `out\play_gate`, `out\p2g_gate`, `out\level_set_node`, `out\threading`, `out\adaptive_*`, `out\bore_suite`. The report lists them under "everything else" with a label each. `out\loop3_D_finale_c1` is the canonical profiling state, the mid-bore finale every speed receipt resumes from; it stays where CLAUDE.md names it. Ask before regenerating it.
