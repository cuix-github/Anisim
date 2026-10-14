@echo off
rem One click from the boulder showcase (experiments/collider_showcase.py) to Maya:
rem writes the play cache if missing, builds the scene if missing, then opens Maya on
rem it with the PF-FLIP compound library on BIFROST_LIB_CONFIG_FILES. Scrub or play:
rem each Maya frame is one cache frame, the liquid particles coloured by speed, with
rem the boulders and the domain box.
rem
rem   play_showcase_in_maya.bat         the run's kept states: 4 frames, 5 s apart
rem   play_showcase_in_maya.bat full    every dumped frame (from a run that kept them)

setlocal
set "HERE=%~dp0"
set "BIFROST_LIB_CONFIG_FILES=%HERE%pfflip_graphs_config.json"
set "SHOW=%HERE%..\out\collider_showcase"
set "CACHE=%SHOW%\play_states"
set "SCENE=%HERE%..\out\play_boulders_states.ma"
if /i "%~1"=="full" (
    set "CACHE=%SHOW%\boulders_play"
    set "SCENE=%HERE%..\out\play_boulders_full.ma"
)

if not exist "%CACHE%\cache.json" (
    if /i "%~1"=="full" (
        echo No full play cache yet: a showcase run writes it.
        exit /b 1
    )
    echo Writing the play cache from the kept states...
    "%HERE%..\.venv\Scripts\python.exe" "%HERE%..\experiments\collider_showcase.py" play || exit /b 1
)
if not exist "%SCENE%" (
    echo Building the Maya scene...
    "C:\Program Files\Autodesk\Maya2027\bin\mayapy.exe" "%HERE%build_particle_play_scene.py" "%CACHE%" "%SCENE%" --obj "%SHOW%\boulders.obj" --box 256 64 128 || exit /b 1
)
echo Opening Maya...
start "" "C:\Program Files\Autodesk\Maya2027\bin\maya.exe" "%SCENE%"
