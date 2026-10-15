# Archived studies

One-off studies that settled a question and are not run any more. Each stays runnable from here (its paths point two levels up), its receipts and charts stay under `results/`, and the record of what it found is in [docs/history.md](../../docs/history.md). Archived 2026-10-15.

| Script | Question it settled | Date | Receipts |
| --- | --- | --- | --- |
| `adaptive_long_run.py` | Does the naive adaptive solve's air stirring level off? No: 430x uniform's at frame 2000; cleanup is indistinguishable from uniform | 2026-10-12 | `results/3d/adaptive/long_run_*.json`, the two charts |
| `adaptive_cost.py` | What the adaptive pressure solve costs at finale size: 612 ms uniform, 2237 naive, 4336 cleanup a step, most of it setup rebuilt every step | 2026-10-13 | `results/3d/adaptive/cost_indicative.json` |
| `level_set_cost.py` | What a level-set solid costs a step: 2 to 4 ms in the assembly, the whole step within noise of no solid | 2026-10-14 | the numbers in the history, milestone "Mesh colliders" |
| `demo_p2g_speed.py` | The fused transfer against the stock splats on the 20-second 2D bore: 8.8x, with the side-by-side video | 2026-10-09 | `results/3d/p2g/p2g_speed_summary.json`, `p2g_speed_chart.png` |
| `render_threading_demo.py` | The threaded build against the serial one on the same bore, as a video with live clocks | 2026-10-10 | `results/3d/threading/` |

`experiments/threading_experiment.py` still imports the demo scene from `demo_p2g_speed.py` here.
