# The solitary wave onto the beach, through the node

The 2D control of `reference2d/solitary_wave.py` run through the Bifrost loop by
`experiments/solitary_beach.py` on 2026-10-16: 20 m of still water, a 6 m hump started
at x = 80 m with the solitary wave's velocity field, a 1:15 beach from x = 160 m to the
shoreline at 460 m as a cut-cell level set, one-metre cells, the production
configuration (ST sampling, escaped particles, adaptive air, multigrid, the fused
transfer), dt 0.1 s, a dump every fifth substep, 26 s. The receipts here are what
`solitary_beach.py finalize --record` copies out of `out/runs/<tag>/`: the summary
with the wave's log every 2 s, and the breaker close-ups.

| Run | Grid | Particles | Wall for 26 s | Face steep, depth under the crest | Lip overhangs 1.5 m, depth | Crest at the overhang |
| --- | --- | --- | --- | --- | --- | --- |
| `beach520x8` | 520 x 56 x 8, an 8 m slice | 1.9 M (428 K liquid) | 147 s on eight threads | 2.84 m at 22.1 s | 1.73 m at 23.1 s | 8.3 m |
| `beach520x1` | 520 x 56 x 1 | 0.24 M (54 K liquid) | 30 s | 2.88 m at 22.1 s | 1.75 m at 23.1 s | 8.3 m |

The marks are the 2D control's (the face steep when three one-metre bands rise more
than 2.5 m per metre of advance under a crest shoaled to 1.1 H; turned over when a
band's front exceeds a lower band's by 1.5 m), measured per dump on the liquid
projected along z. Grilli, Svendsen and Subramanya (1997) put the vertical face for
this slope at H_b / h_b = 2.76, about 2.6 m of water. The 2D reference on the same
beach turns steep at 4.16 m and over at 3.17 m (its receipts in
`../../2d/solitary_wave/`); the one-cell-wide run agreeing with the 8 m slice makes
that a reference-against-node difference, not 2D against 3D, recorded in
`docs/open-questions.md`.

`beach520x8_breaker.png` and `beach520x1_breaker.png`: the side view around the crest
over the last dumps before and after each mark, six pixels a cell, liquid blue, air
grey. The dumps, the play cache and the audit video stay under `out/runs/` (the full
32 m slice, `--grid 520 56 32`, is about 14 GB of dumps and was not run).
