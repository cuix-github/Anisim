# Rendering: configurations, images and evidence

Use the [validated setup reference](../../maya/presets/README.md) to find the
current common water, lighting, foam settings and five tested density variants.
The [frozen-scene registry](../../maya/presets/breaker_validated_20261019.json)
links those variants to exact images, scenes, cameras and validation receipts.
No final preset or temporally validated sequence has been selected.

The local [gallery](../../out/gallery/index.html) browses existing images,
comparisons, videos and grouped frame sequences. It starts with current
candidates; filters expose older diagnostics and failed image checks. Its
thumbnails and manifest are generated under ignored `out/gallery/`, while every
original result stays in place. Links to a Maya file or ASS are file references,
not buttons that launch rendering. The gallery works offline when opened locally.

To regenerate the index after new results arrive, from any working directory:

```powershell
C:/Users/5d149/source/Anisim/.venv/Scripts/python.exe C:/Users/5d149/source/Anisim/experiments/rendering/build_gallery.py
```

This command scans existing media and writes only the gallery and thumbnails.
It does not run Maya, simulation, surface reconstruction or Arnold rendering.
On another checkout, substitute its absolute repository location.

| Study | What it establishes |
| --- | --- |
| [Water-only inspection](studies/breaker-render-study.md) | Geometry visibility, mist/foam controls and remaining glassy appearance |
| [Deep Water / Physical Sky](studies/breaker-deep-water-study.md) | Preset-derived candidate; not selected as a replacement |
| [Neutral water and finer mist](studies/breaker-neutral-mist-study.md) | Existing-cache mist and color observations |
| [Curvature interface](studies/breaker-interface-study.md) | Mask wiring, scattering validation and failed volume attempts |
| [Thin surface foam](studies/breaker-foam-layer-study.md) | Corrected units; working uniform and clustered volume layers |
| [Paired views and 25/75](studies/breaker-foam-view-study.md) | Angles, later frames, crop-wall evidence and intermediate density |
| [Six views and 60/40](studies/breaker-six-view-study.md) | Wider context and a fixed-camera ratio comparison |
| [Surface correspondence](studies/breaker-surface-correspondence.md) | Frozen control and fragment reconstruction limits |

The [grouped instruments](../../experiments/rendering/breaker/README.md)
retain old entry paths as compatibility shims. The former study pages similarly
redirect readers here, so historical links continue to resolve. Compact numerical
evidence remains in `results/rendering/`; media, cache and scene locations remain
under their existing `out/` and legacy `results/` paths. The legacy media migration
is deferred. No historical file was deleted or deduplicated for this organization.

The October 19 layout check verifies 25 unchanged algorithm ASTs (only root
discovery changes), 24 guarded imports through old paths, grouped paths and
package paths, and all 24 old forwarding files explicitly. The historical
`breaker_camera_probe.py` initializes Maya and evaluates a mesh at module scope;
it receives static equivalence/path checking only and was not executed.
The check compares 53,068 original file size/timestamp records, hashes 187 saved
Maya/ASS files and checks 83 extracted scene path strings. There are no changed
historical files, scene hashes or newly missing paths. This is static relocation
validation, not a fresh render or a claim that every old failed scene is usable.

The gallery indexes 2,786 real media files into 1,195 entries, including 18 videos
and 13 grouped sequences; synthetic test fixtures are excluded. Ten current
candidates appear first. One pre-existing unreadable temporary PNG is marked
as such and retained at its original path. Browser checks cover filters, search,
full-size details, asset links, video metadata, frame selection and narrow-screen
layout, with no external web requests. Full local checks and browser screenshots
are under `out/repo_organization_20261019/`.
