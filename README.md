# Anisim: PF-FLIP and ST-FLIP in Bifrost

Anisim implements the core of *Adaptive Phase-Field-FLIP for Very Large Scale Two-Phase Fluid Simulation* (Braun, Bender, Thuerey, SIGGRAPH 2025) and its follow-up *ST-FLIP* (SIGGRAPH 2026) inside Maya's Bifrost, from the public SDK and shipped nodes only. The authors released their grid data structure, not their solver, so this is an independent implementation from the papers' equations. Its sister repository, [Aniso](https://github.com/cuix-github/Aniso), renders Gaussian splats; the two meet at milestone 7, a tsunami in a captured scene, exchanged as per-frame particle caches.

Why these papers, in three sentences. PF-FLIP fills the air with FLIP particles too, so wind, splashes, spray and air entrainment come from physics instead of white-water heuristics. Its phase field costs nothing: it is the splatted particle mass, normalized, so there is no level set and no surface reconstruction in the solver. The only new machinery over textbook FLIP is a variable-coefficient pressure solve, where each cell face is weighted by the mass that landed on it.

## Status

| Milestone | Done | In one line |
| --- | --- | --- |
| 1. The spike | 2026-09-30 | The two-phase projection builds as a Bifrost node and holds a 1000:1 column at rest to 1e-6 |
| 2. The 2D reference | 2026-09-30 | A 300-line numpy spec: hydrostatic rest, dam break, Rayleigh-Taylor |
| 3. The validation gate | 2026-09-30 | The papers' four core claims tested against the obvious alternatives; three supported, one sharpened |
| 4. Phase field from stock nodes | 2026-09-30 | Bifrost's splat reproduces the reference field, mean difference 0.002 |
| 5. The node runs the loop | 2026-10-05 | Lockstep to 1e-6; the whole simulation inside one headless graph |
| 6. 3D at small scale | 2026-10-08 | Multigrid solve, droplets and bubbles, ST-FLIP, adaptive particles, obstacles; the overnight bore suite matches Stoker's bore speed to 1.7 percent |
| 6.5. The adaptive solve | 2026-10-12 | Built and gated in 2D and 3D, kept as an option: at this scale it does not pay |
| 6.75. The fused transfer | 2026-10-09 | One custom node replaces sixty stock splats: 9.4x a substep, threaded to 26x by 2026-10-10 |
| Mesh colliders | 2026-10-14 | Any mesh through the shipped collider node, cut cells, gated in 2D, the node and the graph |
| The 512 run | 2026-10-15 | 512 x 128 x 256 cells, 52 M particles, 20.4 s simulated in 2.4 h, every substep kept |
| 7. The shot | open | Scene chosen 2026-10-13; the next piece is a breaking-wave ocean with a wavemaker |

The full record of every milestone, with its numbers, stands in [docs/history.md](docs/history.md). The settled decisions and their reversals are in [docs/decisions.md](docs/decisions.md), what is still open in [docs/open-questions.md](docs/open-questions.md), every chart and image with what it shows in [docs/evidence.md](docs/evidence.md), and the layout of the local run outputs in [docs/outputs.md](docs/outputs.md).

The October 23 [phase and escape correspondence audit](docs/phase-repair-ladder.md)
records an important qualification to the milestones: independent 20-second
NumPy/Bifrost comparisons still fail their similarity gate after explicitly
matching the phase law. Local contracts and pressure convergence pass, but
long-time interface/transport behavior is not accepted. Experimental phase and
escape options remain off in production; existing production caches are retained.

The [paper correspondence and long-decay audit](docs/paper-correspondence-and-settling.md)
separates verified local contracts, known differences from the papers, and the
visual behavior of one frozen C++ plunge. Independent trajectories diverging is
not by itself proof of a physical defect, and settling is not proof of complete
paper reproduction.

The [bounded verification route](docs/verification-route.md) fixes the next
five gates and records an opt-in 2D density-relaxation candidate. Its local
projection, density-relaxation and still-water checks pass; integration with
escape, transport, ST/adapt and 3D remains unaccepted. Production defaults
remain unchanged.

The [independent interface closure](docs/independent-interface-closure.md) adds
an opt-in, uniform-2D one-way primary-droplet candidate. Controlled drag,
lifecycle, exact restart, common-input execution without density feedback, and
six-second static controls pass. The integrated candidate **fails** on empty
air sampling support in retained impact states; its long recovery/holdout and
3D admission remain withheld. A generic sub-cutoff primary-radius deletion is
fixed only in the new mode. See the [receipt](results/interface_closure_20261024.json).

For the current breaker images, start with [rendering navigation](docs/rendering/README.md),
the [validated setups](maya/presets/README.md) and the local
[media gallery](out/gallery/index.html). The gallery is generated from existing
local outputs and is not hosted on GitHub. Its builder and the grouped breaker
instruments live under `experiments/rendering/`; historical script paths remain
compatible. No simulation or render is needed to regenerate the gallery.

The [wave-feature diagnostic and October 20 continuation](docs/wave-feature-diagnostics.md)
record the small Maya particle demo, its stronger 60-second replacement, and
the separately scheduled 18–30-second 3D spray run. The
[distant perspective observation](docs/rendering/studies/breaker-perspective-study.md)
reuses the 24.9-second cache and 60/40 foam appearance.

The October 20 [force-driven sea plan](docs/force-wave-scene-plan.md) records
the agreed 2D/coarse-3D preview direction and a 30-second production budget
of approximately fifteen hours. The old 18–30-second continuation is cancelled;
its automation was deleted. The force generator has since been implemented for
the previews below; no replacement production run is scheduled.

The new [force-driven wave previews](docs/force-wave-preview.md) replace the
velocity-target wave zone with pre-projection acceleration. They are particle
previews for scene approval, with independent force and checkpoint gates;
production resolution and start time remain undecided.
The [October 20 handoff](docs/handoffs/2026-10-20-ocean-rendering-and-force-waves.md)
collects the validated rendering configurations and current preview state,
including the unresolved Maya display issue.

## How to run

Windows, Maya 2027 with Bifrost 3.1, Visual Studio 2022, Python 3.13.

```
py -3.13 -m venv .venv
.venv\Scripts\pip install -r requirements.txt     (pinned: every recorded result used these versions)
nodes\build.bat                                   (the custom node pack; needed by every 3D run)
```

Then, from the repository root:

| To | Run |
| --- | --- |
| Check the 2D reference against the papers' claims | `cd reference2d` then `..\.venv\Scripts\python.exe validation.py` (80 s) |
| Check the Bifrost node against the reference | `cd gates` then `..\.venv\Scripts\python.exe validate_step.py` (and the other `validate_*.py`; each exits nonzero on failure) |
| Run the boulder showcase, a chunk at a time | `.venv\Scripts\python.exe experiments\collider_showcase.py stage 0 --grid 256 64 128` |
| Play a run in Maya as particles | `maya\play_showcase_in_maya.bat blend boulders512` |
| Render that playback to a video | `maya\render_play_scene.bat out\play_boulders512_blend.ma out\collider_showcase\boulders512_render boulders512` |
| See every run on disk, its size, and what can be rebuilt from what | `.venv\Scripts\python.exe graphs\runs.py report` |

Every run writes under `out\runs\<tag>\` at the root, which is ignored (the layout is in [docs/outputs.md](docs/outputs.md)); the repository keeps code and small numeric evidence only.

## The map

One folder per job.

- `reference2d/`: the executable spec, a pure-Python 2D implementation every Bifrost stage is validated against, with its scenes and gates. Start with its [README](reference2d/README.md) and [VALIDATION.md](reference2d/VALIDATION.md).
- `nodes/`: the custom Bifrost SDK operator pack. `PFFlipSolve3.cpp` is the 3D step with its multigrid solve, `PFFlipP2G.cpp` the fused transfer, `PFFlipAdaptive.cpp` the adaptive pressure solve, `PFFlipSolids.h` the level-set solids.
- `graphs/`: the Bifrost graphs, their generators, and the `.bat` launchers that run them headlessly through bifcmd. `loop3d.py` is the 3D loop runner the gates and experiments share; `colliders.py` turns meshes into a solid field; `runs.py` says where a run's files live under `out/` and reports what is there.
- `gates/`: the validation ladder, one `validate_*.py` per feature, each exiting nonzero on failure.
- `experiments/`: the showcase (`collider_showcase.py`), the solitary wave onto the beach through the node (`solitary_beach.py`), the first breaker from the wavemaker (`piston_beach.py`), the overnight bore suite, the speed receipts and the threading experiment; `archive/` holds the one-off studies that settled a question, with an index.
- `maya/`: playback and rendering in Maya: the play cache, the Bifrost compounds, the scene builders, the one-click launchers and the render presets (`maya/presets/`, the knobs explained).
- `results/`: the committed evidence, charts and receipts: `2d/` for the two-dimensional runs (the Python reference and the 2D Bifrost loop), `3d/` for the 3D node and graphs; the videos they describe stay local.
- `docs/`: the history, the decisions, the open questions, the evidence gallery and the output layout.

## About this repository

The simulator was born inside Aniso under `sim/` and moved here on 2026-10-11 so each repository has one job. Every Aniso commit that touched the simulator came along with its author, date and message; each carries an `Aniso-commit:` trailer naming its original hash, and references like `Aniso#38` point at the pull requests where the work was reviewed. The interface between the two repositories is the per-frame particle cache: positions, velocities, phase, escape flags and particle scale as `.npy` arrays.

Papers: [PF-FLIP](https://ge.in.tum.de/?p=3691) (preprint, video on that page) and [ST-FLIP](https://ge.in.tum.de/?p=4006).

The [literature library](references/README.md) groups the core papers and the
spray, droplet and surface-reconstruction sources cited in PF-FLIP Section 7.
The [October 22 audit](docs/literature/spray-adaptivity-audit.md) separates the
paper's multiresolution spray method from implemented adaptivity, the production
settings and what the existing cache can support.
