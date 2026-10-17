# The 2D reference implementation

`pfflip2d.py` is the executable spec for the Bifrost build: the core PF-FLIP loop from Section 3 of the paper, in about 300 lines of numpy, with none of the scale machinery. Both phases are particles; the phase field is read off the splatted face masses; the pressure solve is variable-coefficient with Jacobi-preconditioned conjugate gradient; the per-phase FLIP blend doubles as viscosity, as in the paper.

Run it (needs numpy and Pillow):

```
python pfflip2d.py            # all three scenes, ~35 s total
python pfflip2d.py dam_break  # one scene
```

## Scenes and what they prove

- **hydrostatic**: liquid under air at 1000:1 in a closed box. Proves the two-phase projection holds the column: velocities stay at FLIP-jitter level with no blow-up.
- **dam_break** (`../results/2d/dam_break.png`): a collapsing column running along the floor under air. Proves the coupled transport works; the front advances monotonically. The collapse is mildly damped by the FLIP blend, a knob to tune when matching the Bifrost build.
- **rayleigh_taylor** (`../results/2d/rayleigh_taylor.png`): heavy fluid over light at Atwood 0.9, the paper's own validation setup. Fingers grow from a seeded ripple into mushrooms. This phenomenon is two-phase-only: a single-phase solver with a free surface cannot produce it at all, which is what makes it the signature test.

## Gates beside the spec

- **solitary_wave.py**: the solitary wave of shallow-water theory, flat and onto a 1:15 beach, the control for the 3D beach runs (`../results/2d/solitary_wave/`).
- **validate_piston.py**: the piston wavemaker, a plate sliding in the tank (`Sim(piston=...)`): at rest the static solid bit for bit, a wave train against the linear dispersion relation and the piston transfer function, and the train onto the beach (`../results/2d/piston/`). Its docstring records why the wave checks read the surface from the liquid's count per column.
- **droplet radii** (`Sim(droplet=...)`, PF-FLIP Section 7): an escaped liquid particle carries its own radius in metres, drawn from a log-normal by a hash of (step, index), shrunk by the d^2 law of evaporation, deleted below a threshold, its drag the sphere's own with the Schiller-Naumann C_d(Re). Gated from `../gates/validate_droplet_node.py`, which runs the law's anchors on this spec (the d^2 law exactly, the hash normal's statistics, Stokes' relaxation, a 1 mm drop's terminal speed) and the node in lockstep against it.
- **spray emission and the droplets' integrator** (`Sim(emit=..., spray_cfl=...)`, PF-FLIP Section 7 after Ihmsen et al. 2012): a liquid particle at the interface in a cell whose velocity gradient exceeds a threshold emits, by a hash of (step, index), a massless spray particle that is always a droplet and that no splat counts, so the other particles are bit-identical with emission on and off; with `spray_cfl` every droplet flies in sub-steps of that many cells, the drag relaxation exact over each at the midpoint slip's rate, the position by the trapezoid. Gated from `../gates/validate_spray_node.py`.

## Deliberate omissions, for later milestones

Adaptive particles and the per-size blend correction, escaped-particle droplets, the divergence-drift correction, and liquid-velocity extrapolation. Each is listed in `../README.md`'s milestones.
