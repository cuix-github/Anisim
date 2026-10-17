# The wave zone's receipt

`receipt.json` is what `gates/validate_wave_zone.py` measured on 2026-10-17 (the wave
zone of `step_3d`: a relaxation zone at the x = 0 wall whose faces and liquid particles
are blended toward a sum of linear wave trains, the wall its flux-neutral inlet):

- `off`: the zone ports at 0 against the plain step, the largest difference over the
  grid, pressure, positions and velocities (0 = bit for bit).
- `lockstep`: one step of a tank the zone has driven for 45 steps, the node at nz = 1
  against the reference's twin: the largest grid, position and velocity differences.
- `train`: the fused loop at 600 x 40 x 2, 20 m deep, one train of 2 m at 8 s ramped
  over 8 s, 50 s: the wavelength and height measured 100 to 400 m beyond the zone
  against linear theory (`L_ratio`, `H_ratio`), the window's still level over time.
- `oblique`: 3D under the production settings, 160 x 20 x 96, 10 m deep, one train of
  1 m at 5 s heading 20 degrees, 15 s: the crest lines' heading from the shift between
  surface profiles 24 m apart, the highest surface.
