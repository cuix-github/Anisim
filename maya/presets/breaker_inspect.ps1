# Surface inspection candidate from the 2026-10-19 wavesea4 24.9 s study.
# D's 10 MCF passes retain more relief; the environment is halfway in stops
# between C (too dark) and E (4x, strong blue highlights). This is a diagnostic
# water-only look, not a completed spray/whitewater preset. Hiding mist also
# hides escaped particles carried by that volume; it never changes the cache.
# Camera and physical cell size remain scene inputs, as with paper.ps1.
if ([string]::IsNullOrEmpty($env:SURF_FOAM_STRENGTH)) { $env:SURF_FOAM_STRENGTH = "0" }
if ([string]::IsNullOrEmpty($env:SURF_MIST_VISIBLE)) { $env:SURF_MIST_VISIBLE = "0" }
if ([string]::IsNullOrEmpty($env:SURF_MCF_ITERS)) { $env:SURF_MCF_ITERS = "10" }
if ([string]::IsNullOrEmpty($env:SURF_SKY_COLOR)) { $env:SURF_SKY_COLOR = "0.06 0.07 0.09" }
. (Join-Path $PSScriptRoot "paper_clean.ps1")
