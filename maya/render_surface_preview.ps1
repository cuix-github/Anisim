# Renders chosen frames of a play cache as Arnold water: the scene is built
# headlessly (build_surface_scene.py, mayapy: the surface and mist compounds, the
# water, foam and fog shaders, sky, sun, colliders, ground, camera), then a minimized GUI Maya session
# opens it and renders the frames (render_surface_preview.py). Arnold batch is a
# separate entitlement; the session render is what the subscription licenses. Output
# stays in out/.
#
#   render_surface_preview.bat <cache_dir> <out_dir> <name> "<frames>" [obj] ["nx ny nz"]
#     frames land in <out_dir>\<name>.####.png with a .log per frame and
#     <name>.notes.txt (times). Knobs beyond these go through the environment:
#     SURF_WIDTH/HEIGHT (1920x1080), SURF_SAMPLES (4), SURF_THREADS (0 = all),
#     SURF_GPU (1), SURF_CAMERA ("x y z tx ty tz"), SURF_POINT_SIZE (0.5, the sphere
#     radius in cells), SURF_DETAIL (0.5), SURF_SMOOTH_ITERS (1), SURF_MESH_SMOOTH (0.1),
#     SURF_MIST_RADIUS (1.5), SURF_MIST_DETAIL (1.0), SURF_MIST_DENSITY (0.08 per cell),
#     SURF_MIST_THRESHOLD (0: particles the smoothed surface leaves out are mist),
#     SURF_GI_VOLUME_DEPTH (3 scattering bounces in the mist), SURF_VOLUME_SAMPLES (2).
#     SURF_LOOK (ocean: dark absorbing water, near-black sky, the sun low behind the spray),
#     SURF_WATER_DEPTH (transmission depth, cells), SURF_SKY_INTENSITY, SURF_SUN_INTENSITY,
#     SURF_SUN_ROTATE ("x y z" degrees; the light shines down its -Z).
#     SURF_SKY_VISIBLE (ocean: 0, the sky dome hidden from the camera), SURF_WATER_COLOR ("r g b"
#     transmission colour).
# Refuses to start while another Maya is open or with less than 20 GB free (the
# machine is shared: a 35 M point surface needs the headroom).
param(
    [Parameter(Mandatory = $true)][string]$Cache,
    [Parameter(Mandatory = $true)][string]$Out,
    [string]$Name = "surface",
    [Parameter(Mandatory = $true)][string]$Frames,
    [string]$Obj = "",
    [string]$Grid = "",
    [int]$MaxMinutes = 240
)
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
if (Get-Process maya -ErrorAction SilentlyContinue) { Write-Error "a Maya is already open; close it first"; exit 2 }
$free = (Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB
if ($free -lt 20) { Write-Error ("only {0:N1} GB free; need 20" -f $free); exit 2 }
$Cache = (Resolve-Path $Cache).Path
if ($Grid -eq "") {
    $index = Get-Content (Join-Path $Cache "cache.json") | ConvertFrom-Json
    if ($index.grid) { $Grid = ($index.grid -join " ") } else { Write-Error "no grid in cache.json; pass it"; exit 2 }
}
$env:BIFROST_LIB_CONFIG_FILES = Join-Path $here "pfflip_graphs_config.json"
$env:MAYA_DISABLE_CER = "1"; $env:MAYA_DISABLE_CIP = "1"; $env:MAYA_DISABLE_CLIC_IPM = "1"
$env:SURF_CACHE = $Cache
$env:SURF_OUT = [System.IO.Path]::GetFullPath($Out)
$env:SURF_NAME = $Name; $env:SURF_FRAMES = $Frames; $env:SURF_GRID = $Grid
$env:SURF_OBJ = if ($Obj -ne "") { (Resolve-Path $Obj).Path } else { "" }
$env:SURF_SCENE = Join-Path $env:SURF_OUT ($Name + ".ma")
$script = (Join-Path $here "render_surface_preview.py") -replace "\\", "/"
$arg = '"python(\"import runpy; runpy.run_path(r''' + $script + ''')\")"'
$marker = Join-Path $env:SURF_OUT ($Name + ".done.txt")
if (-not (Test-Path $env:SURF_OUT)) { New-Item -ItemType Directory -Force $env:SURF_OUT | Out-Null }
if (Test-Path $marker) { Remove-Item $marker }
$notes = Join-Path $env:SURF_OUT ($Name + ".notes.txt")
if (Test-Path $notes) { Remove-Item $notes }

# the scene is built headlessly every run (the knobs are baked into it)
Write-Output "Building $($env:SURF_SCENE) ..."
& "C:\Program Files\Autodesk\Maya2027\bin\mayapy.exe" (Join-Path $here "build_surface_scene.py") 2>&1 |
    Where-Object { $_ -match "^saved|Error|Traceback|line \d+" } | ForEach-Object { "$_" }
if (-not (Test-Path $env:SURF_SCENE)) { Write-Error "the scene was not built"; exit 1 }

Write-Output "Rendering frames $Frames of $Cache (grid $Grid) ..."
$t0 = Get-Date
$p = Start-Process -FilePath "C:\Program Files\Autodesk\Maya2027\bin\maya.exe" -ArgumentList @("-command", $arg) -PassThru -WindowStyle Minimized
while (-not (Test-Path $marker) -and -not $p.HasExited -and ((Get-Date) - $t0).TotalMinutes -lt $MaxMinutes) { Start-Sleep -Seconds 5 }
if (-not $p.HasExited) {
    Stop-Process -Id $p.Id -Force
    if (-not $p.WaitForExit(30000)) { Write-Warning "Maya (pid $($p.Id)) did not close" }
}
Write-Output ("{0:N0} s; marker {1}" -f ((Get-Date) - $t0).TotalSeconds, (Test-Path $marker))
if (Test-Path $notes) { Get-Content $notes }
Get-ChildItem (Join-Path $env:SURF_OUT ($Name + ".*.png")) -ErrorAction SilentlyContinue | ForEach-Object { "{0}  {1:N0} B" -f $_.Name, $_.Length }
if (-not (Test-Path $marker)) { Write-Error "Maya stopped before the marker"; exit 1 }
