# Breaker rendering experiments

Implementations live together here. The former `experiments/<name>.py` paths
remain compatibility entry points for saved launchers and imports. Both routes
execute this implementation; no solver or shading algorithm was changed by
the relocation. Repository-root calculations account for the new directory.

| Area | Instruments |
| --- | --- |
| Initial water comparison | `breaker_render_study.py`, `verify_breaker_render_study.py` |
| Optical candidates | `breaker_deep_water.py`, `breaker_neutral_mist.py`, corresponding `analyze_*.py` |
| Interface mask and scattering | `breaker_interface_mask.py`, `breaker_interface_render.py`, `analyze_breaker_interface.py` |
| Thin surface foam | `breaker_foam_layer.py`, `breaker_foam_render.py`, `breaker_foam_mix.py`, `breaker_foam_mix60.py` |
| Frame and camera observations | `breaker_foam_frames.py`, `breaker_foam_views.py`, `breaker_six_views.py`, `breaker_camera_probe.py` |
| Surface and boundary audit | `breaker_surface_correspondence.py`, `breaker_corner_audit.py`, `breaker_corner_detail.py`, `breaker_crop_audit.py` |
| Bounded reconstruction | `breaker_bounded_mesh.py` |
| Preserved diagnostic attempts | `breaker_interface_volume.py`, `breaker_tiled_mesh.py`, `breaker_tile_bounds.py` |

Start with the [rendering navigation](../../../docs/rendering/README.md) and
[validated configurations](../../../maya/presets/README.md). Failed prototypes
are historical evidence, not interchangeable presets. The scripts can prepare
geometry, write artifacts or start costly rendering: organization and browsing
do not authorize running them. Saved output receipts still prevent overwrites.

`relocation.json` records previous/new paths and original source hashes.
The layout verifier checks algorithm AST equivalence, guarded imports through
old and new paths, saved-scene hashes and local references without evaluating
Maya scenes or tracing rays.
