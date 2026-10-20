# Distant perspective observation

On October 19, 2026, Xue requested one higher-resolution distant perspective
image of the existing `mix_u60_c40` appearance at 24.9 seconds. The new image
is under `out/runs/breaker_mix60_perspective_20261020/beauty_view.png`.

The source is the saved candidate from `breaker_foam_mix60_20261019/C`.
The mesh, normals, topology, water/light settings and foam ASS were checked
against that source. Only the camera and resolution changed: 1920×1080,
perspective rather than orthographic, a 100-metre eye-to-target distance,
and a 56-metre horizontal field at the target plane. The target is
(143, 12.5, 48) metres, with the eye at approximately
(205.089, 38.382, −25.994) metres, 26.382 metres above still water.
Depth of field and mist remain off.

This is still the existing approximately 50×68-metre reconstructed patch,
not the complete 172×80-metre basin. The wider view makes the long breaking
front easier to read, but exposes the finite beach model's straight edges
and does not remove the glossy appearance of the bulk water.

One actual Arnold call completed in 370.0 seconds, with 380.3 seconds for the
launch process. A prior PowerShell execution-policy rejection occurred before
Maya or Arnold started; its receipt/log remain beside the successful run.
The corrected launcher uses the same process-local policy as the existing
render workflow. No second ray-tracing image, simulation or reconstruction
was needed.

The AOV closure relative L1 error is 8.83e-8; volume scattering contributes
about 3.0% of the summed raw RGB components and SSS remains zero. Display
conversion uses the existing ACES-to-sRGB path. These are technical checks,
not a claim that the image is a final production ocean look. Camera state,
source hashes and all check receipts are retained with the image.
