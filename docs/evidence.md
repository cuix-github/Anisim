# Evidence: every chart and image, and what it shows

All of it lives under `results/`: `2d/` for two-dimensional runs (the Python reference and the 2D Bifrost loop), `3d/` for the 3D node and graphs, including their one-cell-thick runs. In the order the work happened. Each picture names the script that made it and the gate or receipt it belongs to; the numbers quoted are the ones recorded at the time in the [history](history.md). Videos are rendered by the same scripts but stay local, see [outputs](outputs.md).

## 1. The 2D reference (milestone 2, September 2026)

`reference2d/pfflip2d.py`, three scenes.

![Dam break frame strip](../results/2d/dam_break.png)

A collapsing column running along the floor under air; the front advances monotonically, mildly damped by the FLIP blend.

![Rayleigh-Taylor](../results/2d/rayleigh_taylor.png)

Heavy fluid over light at Atwood 0.9, the paper's own validation setup. Fingers grow from a seeded ripple into mushrooms; a single-phase solver cannot produce this at all.

## 2. The validation gate (milestone 3, 2026-09-30)

`reference2d/validation.py`; verdicts in [VALIDATION.md](../reference2d/VALIDATION.md).

![Validation overview](../results/2d/validation/validation_overview.png)

The overview sheet of the four experiments.

![Air pocket curves](../results/2d/validation/air_pocket_curves.png)

Claim 1, supported decisively: a slab falls on a sealed pocket of air. The two-phase slab settles on the cushion and the gap holds at 10.6 cells; the single-phase slab falls through its vacuum, 10.6 down to 3.3 and still accelerating. The frame strips: single-phase above, two-phase below.

![Single-phase strip](../results/2d/validation/air_pocket_single_phase.png)
![Two-phase strip](../results/2d/validation/air_pocket_two_phase.png)

![Transport](../results/2d/validation/transport.png)

Claim 2, supported with a caveat: the reversing-vortex benchmark. The grid scalar loses the thin filament to numerical diffusion (L1 error 293 cell units); the particle-carried phase returns essentially exact.

![Viscosity](../results/2d/validation/viscosity.png)

Claim 3, the law holds, the constant is ours: measured viscosity is proportional to (1 - alpha) / dt, but the constant is about 1/30 rather than the paper's 1/6, which belongs to their kernels. Calibrate per implementation.

![Mixing](../results/2d/validation/mixing.png)

Claim 4, supported: the wrong-side particle fraction in the dam break peaks at 0.14 percent and does not grow.

`reference2d/compare_frames.py`: plain single-phase FLIP (top row) against PF-FLIP (bottom row), same scene, same seeds, same instants.

![Compare dam break](../results/2d/validation/compare_dam_break.png)

The dam break, where air escapes freely: broadly similar, the honest negative control.

![Compare air cushion](../results/2d/validation/compare_air_cushion.png)

The sealed air cushion, where only PF-FLIP has the physics.

![Compare plunge](../results/2d/validation/compare_plunge.png)

The plunge: trapped pockets and a cushioned crash against a free fall.

## 3. The phase field from stock nodes (milestone 4, 2026-09-30)

![Phase field, Bifrost against the reference](../results/2d/validation/phase_field_bifrost.png)

`gates/compare_phase_field.py`: the same jittered slab splatted by Bifrost's stock node and by the reference. Mean difference 0.002, maximum 0.08 confined to the one-voxel interface band, interface midpoint within 0.03 voxels.

## 4. The loop inside Bifrost (milestone 5, 2026-10-05)

`gates/validate_loop.py`: the whole simulation inside one headless graph, against the reference and against physics.

![Bifrost dam break sheet](../results/2d/validation/bifrost_dam_break.png)

Side by side with the reference; the front stays within 1.9 cells over the full run (1.0 cell once the fused transfer replaced the stock kernel).

![Bifrost dam front](../results/2d/validation/bifrost_dam_front.png)
![Bifrost cushion gap](../results/2d/validation/bifrost_cushion_gap.png)

The front-position curve, and the sealed cushion's gap holding at 10.6 cells with no reference in the loop.

![Bifrost air cushion](../results/2d/validation/bifrost_air_cushion.png)
![Bifrost Rayleigh-Taylor](../results/2d/validation/bifrost_rayleigh_taylor.png)

The physics anchors: the air cushion and 19:1 Rayleigh-Taylor growing textbook mushrooms.

![Splash demo](../results/2d/validation/bifrost_splash_demo.png)

`experiments/demo_splash.py`: the violent three-way plunge. Air-jet speeds reach 505 cells per second as the trapped film squirts out at impact, CFL 5 at that instant; the measured reason milestone 6 needed adaptive stepping.

## 5. ST-FLIP in 2D (milestone 6, October 2026)

`reference2d/st_compare.py`: each scene three ways, the trusted small dt, PF-FLIP at a large dt, and ST-FLIP at the same large dt.

![ST-FLIP dam break](../results/2d/validation/stflip_dam_break.png)
![ST-FLIP plunge](../results/2d/validation/stflip_plunge.png)

![ST energy curves](../results/2d/validation/st_energy_curves.png)

`reference2d/st_energy.py`: kinetic energy and surface roughness against time; aliasing shows as sustained agitation the reference does not have.

## 6. The fused transfer layer (milestone 6.75, 2026-10-09)

![P2G speed chart](../results/3d/p2g/p2g_speed_chart.png)

`experiments/demo_p2g_speed.py`, receipts in [p2g_speed_summary.json](../results/3d/p2g/p2g_speed_summary.json): the same 20-second bore through the stock-splat graph and the fused node. At finale scale, 96.4 s a substep to 10.22 s, 9.4x (`experiments/bench_p2g_finale.py`); 8.8x on the 2D demo. The side-by-side video with live clocks is local.

## 7. Threading (2026-10-10)

![Threading phases](../results/3d/threading/threading_phases.png)

`experiments/threading_experiment.py`, the record in [results/3d/threading/README.md](../results/3d/threading/README.md) with [receipts.json](../results/3d/threading/receipts.json) and [validation.txt](../results/3d/threading/validation.txt): destination-owned slabs, every saved array bit-identical to the serial build, 8.3 to 4.5 s a substep at 32 workers, then 3.66 s with the solve sweep threaded.

## 8. The adaptive pressure solve (milestone 6.5, 2026-10-10 to 2026-10-13)

Gates: [stage1_gate.json](../results/2d/adaptive/stage1_gate.json) and [stage3a_gate.json](../results/2d/adaptive/stage3a_gate.json) from `reference2d/validate_adaptive_solve.py` and `validate_adaptive_projection.py`.

![Long run, resting tank](../results/3d/adaptive/long_run_tank_chart.png)

`experiments/adaptive_long_run.py`, [long_run_tank.json](../results/3d/adaptive/long_run_tank.json): a half-full tank at rest for 100 s under the uniform, naive and cleanup solves, three seeds each. Cleanup is indistinguishable from uniform; naive's air stirring grows without levelling off, 1.49 cells per second at frame 2000, 430x uniform's.

![Long run, dam break](../results/3d/adaptive/long_run_dam_chart.png)

[long_run_dam.json](../results/3d/adaptive/long_run_dam.json): the dam break for 60 s. Every arm inside the uniform seed bar at every checkpoint; iterations per substep 36 uniform, 40 naive, 60 cleanup.

[cost_indicative.json](../results/3d/adaptive/cost_indicative.json), `experiments/adaptive_cost.py`: one finale-size step, no obstacles. The pressure phase costs 612 ms uniform, 2237 naive and 4336 cleanup at 1e-6; most of it is setup rebuilt every step.

## 9. The bore suite (milestone 6, 2026-10-08)

[summary.json](../results/3d/bore_suite/summary.json), `experiments/run_bore_suite.py`, judged against Stoker's analytic bore speed.

| Configuration | Grid | Bore speed / Stoker | Wall |
| --- | --- | --- | --- |
| Trusted (small dt) | 80 x 32 x 40 | 1.059 | 86 min |
| ST at 8x dt | 80 x 32 x 40 | 1.063 | an eighth of trusted |
| Production (ST + adaptive particles) | 80 x 32 x 40 | 0.916 | 16.6x cheaper than trusted |
| Finale, production | 256 x 64 x 128 | 0.983 | 107 min |

## 10. The showcases (2026-10-14 and 2026-10-15)

Numbers only, in `out/collider_showcase/*_summary.json` locally: the 256 bore over three mesh boulders (4.13 h on the stock graph, no penetrations), and the 512 run (2.40 h on the fused graph, 134 M particles seeded, 52 M at the end, 20.38 s simulated). The frames and videos are described in [outputs](outputs.md).

## 11. The solitary wave (2026-10-16)

`reference2d/solitary_wave.py`, the 2D control for the beach runs; receipts in [receipt.json](../results/2d/solitary_wave/receipt.json), [receipt_production.json](../results/2d/solitary_wave/receipt_production.json) and [receipt_fine.json](../results/2d/solitary_wave/receipt_fine.json). Cells of one metre; the strips are drawn at two pixels a cell, liquid blue, air grey.

![Solitary wave on a flat bed](../results/2d/solitary_wave/flat.png)

The hump with its velocity field, 400 x 48 cells, h = 20 m, H = 6 m: it runs at 0.966 of the analytic speed and keeps its height, 5.88 to 5.80 m over 11 s.

![The same hump at rest](../results/2d/solitary_wave/rest.png)

The negative control: seeded without its velocity it splits into two half-height waves, 2.96 m each at mid-run.

![The beach](../results/2d/solitary_wave/beach.png)

The wave onto the 1:15 beach (toe at 160 m, shoreline at 460 m), the last eight frames before it turns over.

![The beach, close up](../results/2d/solitary_wave/beach_zoom.png)

The last two seconds around the crest, six pixels a cell: the face goes steep at 4.16 m of water and the lip overhangs 1.5 m at 3.17 m. Grilli's vertical face for this slope is at about 2.6 m; the Bifrost node reaches it (the next entry). The production settings ([beach_production_zoom.png](../results/2d/solitary_wave/beach_production_zoom.png)) and half-metre cells ([beach_fine_zoom.png](../results/2d/solitary_wave/beach_fine_zoom.png)) change the reference's marks by less than a metre of depth.
