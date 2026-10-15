# Outputs: what the runs leave on the machine

Everything a script writes lands under `out\` at the repository root, which git ignores. Nothing in it is in the repository; this page says what is there, what it is for, and which folders are canonical. Sizes are from 2026-10-15.

## The layout today

| Folder | What | Size |
| --- | --- | --- |
| `out\loop3_<tag>_c<N>` | One chunk of a 3D run: per-substep dumps (`frame`, `vel`, `esc`, `scl`, `ph`, `diag`), the resume state (`final*.npy`) | the 512 run: 457 GB over `boulders512_c0` to `_c5` |
| `out\collider_showcase\<tag>_play24` | The 24 fps blend cache of a run: `liquid.####.npy` and `speed.####.npy`, `cache.json` | 252 GB for `boulders512` (487 frames) |
| `out\collider_showcase\<tag>_play` | The per-dump play cache, superseded by the blend cache | 124 GB for `boulders512` |
| `out\collider_showcase\<tag>_states` | One frame per chunk end, for a quick look | 3 GB |
| `out\collider_showcase\<tag>_c<N>.json`, `<tag>_summary.json` | The run's receipts: wall per substep, iterations, peak speeds, penetration audit | small |
| `out\collider_showcase\<tag>_render` | The playblast frames and the mp4 of a run | 0.1 GB |
| `out\collider_showcase\solid_phi_<tag>.npy` | The colliders as a signed-distance field, computed once per run | up to 0.1 GB |
| `out\play_<tag>_<mode>.ma` | The Maya scenes the launchers build | small |
| `out\loop3_D_finale_c1` | The canonical profiling state, the mid-bore finale every speed receipt resumes from. Ask before regenerating. | 0.4 GB |
| `out\loop3_*` (small) | The gates' and studies' runs: level-set sweeps, ST and adaptive checks, smoke tests | 0.1 GB each |
| `out\threading`, `out\level_set_node`, `out\adaptive_long`, `out\play_gate`, `out\surface_gate` | Each gate's or study's working folder | under 1 GB |
| loose files at `out\` root | Logs and sweep scripts from the first weeks (`taskb_*.log`, `*_sweep.py`), stray `.npy` | small |

Videos rendered by the demos sit beside their receipts under `results\3d\` (the fused-transfer side-by-side in `p2g\`, the threading demo in `threading\`, the 3D demo in `demo3d\`), and the 2D demos' videos and frame caches under `results\2d\validation\`; all ignored.

## Leftovers that can go

Listed, not deleted; deleting a run's data is the owner's call.

- `out\loop3_dry512_c0` (8.6 GB) and `out\loop3_dry384_c0` (6.3 GB) with `out\collider_showcase\dry384_*`: pre-flight memory tests of 2026-10-15, inputs only.
- `out\collider_showcase\boulders512_play` (124 GB): the per-dump play cache, every frame of which the blend cache also holds, interpolated.
- `out\loop3_smoke_c0`, `loop3_smoke_c1`, `loop3_smoke_fused_c0`, `out\smoke_blend24`, `out\surface_test`: smoke-scene test runs.
- The loose logs and sweep scripts at the `out\` root.

## The convention going forward

One folder per run, named by the run, holding everything the run produced:

```
out\runs\<tag>\
    dumps\c<N>\        the chunks as today
    cache\             the blend cache (and the per-dump cache if asked for)
    states\
    renders\<set>\     frames, logs, the mp4
    scenes\            the Maya scenes
    receipts\          the chunk and summary JSON
out\gates\<name>\      a gate's working folder
out\studies\<name>\    an experiment's working folder
out\profiling_state\   the canonical mid-bore state
```

The launchers and the loop runner write today's layout; moving to this one is a code change with the resume and playback gates re-run, planned, not done.
