# The first breaker from the wavemaker

`experiments/piston_beach.py`, the night of 2026-10-16: the piston plate of `step_3d`
(gates/validate_piston_node.py) makes a wave train in 12 m of water that runs onto a
1:10 beach and plunges. The scene in metres: the plate 1 m thick at x = 2 m, its top at
18 m, a 4 m stroke at a period of 8 s; the toe at 36 m, the shoreline at 156 m; the
domain 172 x 24 x 16 m at four cells a metre (688 x 96 x 64, 4.2 M cells), the
production configuration (ST, escaped particles, adaptive air, multigrid, the fused
transfer), dt 0.025 s, 24 s.

| Run | Grid | Particles | Wall for 24 s | The lip over the face (overhang mark) |
| --- | --- | --- | --- | --- |
| `pistonbeach4x1`, the preview | 688 x 96 x 1 | 0.41 M (143 K liquid) | 109 s on eight threads | 21.4 s at x = 113 m, 4.3 m of water |
| `pistonbeach4` | 688 x 96 x 64 | 26 M seeded, 12 M after the air coarsens (9.1 M liquid) | 1.79 h on all threads, 5.4 to 8.1 s a substep | 21.2 s at x = 111 m, 4.5 m of water |

The leading crests of the train are small (1.4 to 1.9 m); the first full-height wave
reaches the slope at about 18 s, shoals to a 4.4 m crest and plunges. The close-ups
(`*_breaker.png`, the middle of the width, liquid blue over air grey, 0.4 s apart from
20.0 s) show the face steepen, the lip curl over, the jet land and trap a tube of air,
and the splash throw water forward; the one-cell-wide preview and the 3D run agree on
the break within 0.2 s and 2 m. The escaped spray grows from 5,300 particles at the lip
to 14,500 at 23.8 s, in the splash-up.

The marks are indicative here, the close-ups are the evidence: the "face steep" mark of
the solitary-wave scenes never fires on this beach, because the run-up water ahead of
the crest stands above still water in the same bands and is read as the front; the
overhang mark (1.52 m at 21.2 s) agrees with the close-ups but rests on the same bands.

`pistonbeach4.json` also records the chunk receipts and the settings of two renders of
this run through `maya/render_surface_preview.bat` (the ocean look, 1920 x 1080, frames
at 22.4 and 23.6 s); the images stay in `out/runs/pistonbeach4/renders/`, as footage
does. The dumps (29 GB, a dump every second to 18 s, every 0.2 s after) stay in
`out/runs/pistonbeach4/`.
