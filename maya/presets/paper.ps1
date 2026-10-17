# The paper look: the breaking-wave stills of PF-FLIP and ST-FLIP, as settled on
# 2026-10-16/17 (Anisim PR #33). Dot-sourced by render_surface_preview.ps1 -Preset paper
# before the scene is built. Every knob is set only when the caller has not set it, so
# an explicit environment variable wins over the preset. maya/presets/README.md says
# what each knob controls; the launcher's header says where each one is read.
#
# Two things the preset does not set, because they belong to the scene, not the look:
#   SURF_METRES_PER_CELL  the simulation's cell in metres (the breaker: 0.25); the spray's
#                         droplet clouds are sized from it
#   SURF_CAMERA           "x y z tx ty tz" in cells; the breaker's crop camera is
#                         "505 96 -52 505 65 15"
function Set-Default([string]$name, [string]$value) {
    if (-not (Test-Path "Env:$name") -or [string]::IsNullOrEmpty((Get-Item "Env:$name").Value)) {
        Set-Item "Env:$name" $value
    }
}

# --- the surface: ST-FLIP Appendix B (maya/particle_surface) ---------------------
Set-Default SURF_POINT_SIZE   "0.5"   # each particle a sphere of half a cell (the paper's 0.5 dx)
Set-Default SURF_DETAIL       "0.5"   # the surface grid's voxel, half a cell (twice the simulation's resolution)
Set-Default SURF_SMOOTH_ITERS "0"     # Gaussian passes before the flow: none, the paper has none
Set-Default SURF_MCF_ITERS    "30"    # mean-curvature flow iterations under the mask (the paper's k_psi)
Set-Default SURF_MASK_SIGMA   "2"     # the mask's blur in cells (the paper's 2 dx)
Set-Default SURF_MASK_THETA   "2"     # the mask's noise threshold (the paper's theta)
Set-Default SURF_MASK_ZETA    "5"     # the mask's slope factor (the paper's zeta)
Set-Default SURF_MASK_FLOOR   "0"     # 0: the mask as the paper's; 1: the flow unmasked
Set-Default SURF_MASK_DILATE  "1"     # our one departure: the mask's numerator pushed out a voxel
Set-Default SURF_MESH_SMOOTH  "0.1"   # the mesher's own smoothing (volume_to_mesh)

# --- the spray: PF-FLIP Section 7 (maya/particle_mist) ---------------------------
Set-Default SURF_MIST_THRESHOLD "0"      # particles the surface leaves out are mist (and every droplet is)
Set-Default SURF_MIST_DETAIL    "0.5"    # the spray grid's voxel, half a cell (the paper's 2x; it allows 4x)
Set-Default SURF_MIST_RADIUS    "1"      # each mist particle's droplet cloud, in cells
Set-Default SURF_MIST_DENSITY   "0.001"  # the fraction of a parcel's water drawn as droplets: a choice, not a measurement
Set-Default SURF_DEFAULT_RADIUS "0.0001" # the droplet radius (m) for a mist particle without one

# --- the look: dark absorbing water, black sky, the sun low behind the spray --------
Set-Default SURF_LOOK          "ocean"
Set-Default SURF_WATER_DEPTH   "40"            # transmission depth in cells (10 m at the breaker's cell)
Set-Default SURF_WATER_COLOR   "0.3 0.55 0.6"  # transmission colour
Set-Default SURF_SUN_ROTATE    "-40 160 0"     # the sun behind the spray, toward the camera
Set-Default SURF_SUN_INTENSITY "6"
Set-Default SURF_SKY_INTENSITY "12"
Set-Default SURF_SKY_VISIBLE   "0"             # the sky lights but the camera sees black

# --- Arnold -----------------------------------------------------------------------
Set-Default SURF_SAMPLES         "4"   # AA samples
Set-Default SURF_GI_VOLUME_DEPTH "3"   # scattering bounces in the mist
Set-Default SURF_VOLUME_SAMPLES  "2"
