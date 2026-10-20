# Force-driven sea: agreed direction and capacity envelope

On October 20, 2026 the owner accepted a sea driven by roughly one thousand
components with different amplitudes, periods, phases and directions, and a
two-stage particle preview: 2D for propagation and breaking through time, then
coarse 3D for directional intersections and variation along the crest. Neither
preview needs Arnold. The new force generator has not been implemented or
validated. This page is a design and budget, not an executable run configuration.

The old 18–30 s continuation prepared in Anisim PR #41 is cancelled. Its
`anisim-18-30` automation was deleted on October 20, before any production
segment ran. Its frozen plan and all previous simulation data remain local.
Do not launch that plan as a fallback for the new scene.

The owner revised the next production target to **30 simulated seconds in
approximately 15 wall-clock hours**, including preparation of an inspectable
particle cache. Use the workstation's 32 logical threads. The earlier machine
availability was October 20 at about 07:00–22:00/23:00 Toronto time; no replacement
automation has been created. The accepted sequence still requires inspecting
the previews before choosing the production scene. Development and validation
of a new driving method are separate from the 15-hour production budget.

## Physical design

The solver should accept a general external force/acceleration field before
pressure projection. The scene generator belongs in a separate Bifrost node
or compound, not in scene-specific branches of `step_3d`. In physical units
the intended sequence is `u_star = u_P2G + dt * (g + f/rho)`, followed by the
variable-density pressure projection and the existing particle update. Units
and phase weighting must be explicit in the eventual interface.

Do not use a velocity-target controller disguised as a force, and do not
overwrite grid or particle velocities after projection. The new scene uses
`zone_length=0`; it must not inherit the zone's exemption from particle escape.
Adding forces does not by itself eliminate pressure compatibility problems,
boundary reflections or high air speeds.

Two or three energetic shoreward wave groups should carry most of the sea's
energy. Oblique companions add crossing and crest variation, with a resolved
short-wave band providing surface structure. Component count and total spectral
energy are independent controls. Fix the seed and preserve the prescribed
energy when discretizing the spectrum or reducing it to 2D. A thousand
unrelated point kicks would be a different experiment.

The proposed primary wave-height target remains 4–6 m, measured trough to
crest, with periods around 8–12 s. These are calibration targets outside the
source region, not guaranteed output from guessed force amplitudes. Wave
speed follows the chosen depth, period and wavelength; do not multiply
velocities to make the shot faster.

For a 30 s production run, shorten the uninteresting travel/startup interval
through source placement and a smooth energy ramp. Aim for the first useful
breaking event around 10–15 s, with subsequent events before 30 s. These are
scene-design targets to check in the preview. The earlier recipe only began
nearshore breaking around 23 s and left little time for its aftermath.

The source region, clean observation region, beach/run-up area and boundary
treatment must fit together. Return flow, reflection and air headroom are part
of the design. A damping force outside the observation region is a possible
absorber, but its reflection and influence must be measured; there is no
validated absorbing boundary in this proposal yet. Avoid assuming that a
larger closed box is an open ocean.

Published momentum-source wavemakers establish feasibility, not a drop-in
PF-FLIP formula: [Feng and Wu, 2019](https://doi.org/10.1155/2019/1308960).
The spectral organization is inspired by
[Tessendorf's course notes](https://people.computing.clemson.edu/~jtessen/reports/papers_files/coursenotes2004.pdf);
FFT surface synthesis alone does not supply the nonlinear breaking dynamics.

## Preview and validation sequence

First reproduce a no-force control, then compare the force-enabled node with
an independent 2D reference on identical inputs. A single small wave gives
measurable wave height, period, propagation speed and mean-level drift; use
it to establish the forcing transfer before assembling the mixed sea.

The mixed-sea 2D preview uses the same depth, longitudinal bed, physical
forcing specification and spray settings as the intended 3D run. Preserve
energy in the reduction of the directional spectrum. It cannot predict the
same cross-wave events as 3D, and a coarser grid can change breaking and
emission. Report these differences explicitly. The older 2D diagnostic's
twenty-fold emission-rate increase is not the new production default.

The coarse 3D preview checks directions, crest variation, intersections,
side-boundary effects and the location/time of breaking. Compare the particle
behavior, including escaped liquid and newly emitted spray, before spending
the production budget. Check selected intervals at a finer resolution.
Keep controls and candidate outputs in separate directories.

Enable the existing ST, escaped-particle, droplet, extra-spray and adaptive-air
particle paths with consistent settings. Adaptive pressure is a separate
feature: the current implementation falls back to the uniform solve with
the solid beach. Do not claim an adaptive pressure speedup or relax pressure
convergence to meet the deadline.

## CPU budget

The machine reports a Ryzen 9 9950X3D, 16 physical cores and 32 logical
processors, with about 62 GiB of usable RAM. The old production run already
used 32 threads (`out/night/sea_night2.log`), so there is no additional
32-thread speedup to apply to its receipts.

The five receipts for `wavesea4` total **16.893 hours** for 30 s at
688×96×320 cells, dt=0.025 s, or 1,200 steps. The mean is 50.678 s/step;
the individual chunks average 31.19, 50.13, 54.40, 58.35 and 59.32 s/step.
The initial quiet chunk alone is therefore an unsuitable timing sample.

An illustrative spatial budget is approximately 180×32×72 m, with more air
headroom than the old 24 m-high domain. It is not a final choice of physical
extent, beach or absorber width. A tentative dt=1/48 s gives 1,440 steps for
30 s and simple exact sampling at 12 or 24 Hz. Its stability remains to be
tested under the new forcing.

`experiments/wave_budget.py` reads the old receipts and computes the following
screening scenarios. It scales cost linearly with cell count and step count,
then adds hypothetical overhead and one hour for cache preparation/checks.
It does not measure the new solver configuration. Changed water occupancy,
pressure iterations, source evaluation, spray growth and output cadence can
all invalidate that scaling.

| Cell width | Rounded grid | Cells | With 25% overhead + 1 h | With 50% overhead + 1 h |
| --- | --- | ---: | ---: | ---: |
| 0.30 m | 600×107×240 | 15.41 M | 19.47 h | 23.17 h |
| 1/3 m | 540×96×216 | 11.20 M | 14.42 h | 17.11 h |
| 0.35 m | 515×92×206 | 9.76 M | 12.70 h | 15.04 h |

Start capacity planning around **10–11 M cells**, using 0.35 m as the more
conservative candidate and 1/3 m if actual timing permits. These cells are
larger than the old 0.25 m cells, so this trades small-scale resolution for
runtime. The proposed gain is a better designed and more eventful scene, not
finer spatial detail.

Select the production configuration using serial measurements of startup,
propagation and representative breaking/spray states, including I/O.
At dt=1/48 s the entire solver workload must average at most **35 s/step**
to leave one of the fifteen hours for cache work. Aim below that boundary;
do not treat it as a fixed timestep limit or a convergence tolerance.
If smaller stable steps are required, recalculate the budget. Notify the
owner if the projected completion exceeds the agreed window; do not silently
alter resolution, physics or tolerance mid-run.

## Disk budget

The October 20 read-only audit found approximately **1,168.6 GiB free**.
Anisim's `out/` occupies 635.86 GiB by logical file size. The old
`wavesea4/dumps` is 439.31 GiB, while `wavesea4/cache` is only 9.47 GiB.
Deleting the latter would save little. The dumps are the original simulation
evidence, not just a redundant playback copy.

The owner permits deleting old-scene caches if necessary. No deletion is
currently needed or performed. Preserve original dumps, accepted rendering
artifacts and gallery links while the new budget fits. Do not repeat the
previous blocked deletion of the initial 2D demo as an indirect cleanup step.

For an explicit storage scenario, assume 32 million total particles, exact
12 Hz full dumps over 30 s, and two-second restart segments:

| Allocation | Estimated GiB |
| --- | ---: |
| Full dumps: positions, velocities, phase, escape, scale, radius and spray | 472.1 |
| Both segment inputs and complete restart outputs | 60.8 |
| Lightweight particle-view cache | 18.8 |
| Other inputs, previews, logs and growth allowance | 60.0 |
| Total | 611.6 |

Use **650 GiB as a planning cap**, leaving **250 GiB free disk reserve**.
Those two allocations leave approximately 269 GiB of additional margin on
the observed drive. Particle counts and actual output sizes must replace
these assumptions before launch. Neither the storage guard nor the new
exporter is implemented by this planning change.

The proposed display copy samples only bulk liquid for interactive inspection
and retains escaped liquid and emitted spray. Full dumps retain all particles.
The 18.8 GiB display estimate allows one million bulk points plus one million
escape/spray points per frame; it is not permission to drop excess spray.
Inspect the real counts and budget any excess. Do not apply index interpolation
across particle births/deaths. A 12 Hz cache plays at its real timestamps;
showing it on a 24 fps Maya timeline must not halve the duration.

Avoid duplicating the whole basin into a second full particle cache before
specific frames or regions are requested. The deliverable for this run is a
simulation and usable particle inspection cache, not an Arnold sequence.

The complete local audit and arithmetic are reproducible with:

```
.venv\Scripts\python.exe experiments\wave_budget.py --output out\wave_budget_20261020.json
```
