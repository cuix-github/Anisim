# The paper look, clean on calm water (2026-10-18, the basin's arm S3). Dot-sourced by
# render_surface_preview.ps1 -Preset paper_clean. It is paper.ps1 with one render knob
# changed, and it expects a cache made with one cache-side knob; README.md says why.
#
#   render:  SURF_MIST_THRESHOLD 0.75 (paper: 0). Water the smoothed surface leaves
#            less than 0.75 cells outside stays water instead of becoming a droplet cloud.
#   cache:   escaped particles within one cell of the bulk water carry no droplet radius
#            (experiments/wave_basin.py cache or section --spray-above 1). This one is not
#            an environment knob: the play cache's radius.####.npy carries it, so the cache
#            has to be written with it.
#
# Every other knob is paper.ps1's, set only where the caller has not set it.
if ([string]::IsNullOrEmpty($env:SURF_MIST_THRESHOLD)) { $env:SURF_MIST_THRESHOLD = "0.75" }
. (Join-Path $PSScriptRoot "paper.ps1")
