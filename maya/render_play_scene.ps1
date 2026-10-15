# Renders a particle play scene (build_particle_play_scene.py) to PNG frames and an
# mp4: what the viewport shows, the liquid particles coloured by speed, the colliders
# and the domain box from the scene's persp camera, 1920x1080, 8x multisampling, no
# grid or HUD. The frames are an off-screen viewport playblast from a GUI Maya
# session (playblast_play_scene.py): Maya's Hardware 2.0 batch renderer evaluates
# the Bifrost graph (its bounding box is right) but does not draw its points, and
# Arnold batch is not licensed here. Maya opens minimized, writes the frames, and is
# closed here once its done-marker appears (its own quit does not take). Output
# stays in out/ (footage never enters the repository). render_play_scene.bat wraps it:
#
#   render_play_scene.bat <scene.ma> <out_dir> <name> [first] [last] [fps]
#     frames land in <out_dir>\<name>.####.png, the video in <out_dir>\<name>.mp4;
#     first/last default to the scene's timeline, fps to 24.
param(
    [Parameter(Mandatory = $true)][string]$Scene,
    [Parameter(Mandatory = $true)][string]$Out,
    [string]$Name = "play",
    [string]$First = "",
    [string]$Last = "",
    [int]$Fps = 24
)
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
$env:BIFROST_LIB_CONFIG_FILES = Join-Path $here "pfflip_graphs_config.json"
$env:MAYA_DISABLE_CER = "1"; $env:MAYA_DISABLE_CIP = "1"; $env:MAYA_DISABLE_CLIC_IPM = "1"
$env:PLAY_SCENE = (Resolve-Path $Scene).Path
$env:PLAY_OUT = [System.IO.Path]::GetFullPath($Out)
$env:PLAY_NAME = $Name; $env:PLAY_FIRST = $First; $env:PLAY_LAST = $Last
$env:PLAY_WIDTH = ""; $env:PLAY_HEIGHT = ""
$script = (Join-Path $here "playblast_play_scene.py") -replace "\\", "/"
$arg = '"python(\"import runpy; runpy.run_path(r''' + $script + ''')\")"'
$marker = Join-Path $env:PLAY_OUT ($Name + ".playblast.txt")
if (-not (Test-Path $env:PLAY_OUT)) { New-Item -ItemType Directory -Force $env:PLAY_OUT | Out-Null }
if (Test-Path $marker) { Remove-Item $marker }

Write-Output "Playblasting $($env:PLAY_SCENE) ..."
$t0 = Get-Date
$p = Start-Process -FilePath "C:\Program Files\Autodesk\Maya2027\bin\maya.exe" -ArgumentList @("-command", $arg) -PassThru -WindowStyle Minimized
while (-not (Test-Path $marker) -and -not $p.HasExited) { Start-Sleep -Seconds 5 }
if (-not $p.HasExited) {
    Stop-Process -Id $p.Id -Force
    if (-not $p.WaitForExit(30000)) { Write-Warning "Maya (pid $($p.Id)) did not close" }
}
if (-not (Test-Path $marker)) { Write-Error "Maya exited before writing the frames"; exit 1 }
$frames = Get-ChildItem (Join-Path $env:PLAY_OUT ($Name + ".*.png")) | Sort-Object Name
Write-Output ("{0} frames in {1:N0} s" -f $frames.Count, ((Get-Date) - $t0).TotalSeconds)
$start = [int]($frames[0].BaseName -replace "^.*\.", "")
& ffmpeg -y -loglevel error -framerate $Fps -start_number $start -i (Join-Path $env:PLAY_OUT ($Name + ".%04d.png")) -c:v libx264 -pix_fmt yuv420p -crf 18 (Join-Path $env:PLAY_OUT ($Name + ".mp4"))
if ($LASTEXITCODE -ne 0) { Write-Error "ffmpeg failed"; exit 1 }
Write-Output ("Rendered " + (Join-Path $env:PLAY_OUT ($Name + ".mp4")))
