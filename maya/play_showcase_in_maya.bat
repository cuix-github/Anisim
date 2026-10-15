@echo off
rem One click from the boulder showcase (experiments/collider_showcase.py) to Maya:
rem writes the play cache if missing, builds the scene if missing, then opens Maya on
rem it with the PF-FLIP compound library on BIFROST_LIB_CONFIG_FILES. Scrub or play:
rem each Maya frame is one cache frame, the liquid particles coloured by speed, with
rem the boulders and the domain box (the grid, read from the cache's index).
rem
rem   play_showcase_in_maya.bat [states / full / blend] [tag]
rem     states  the run's kept states, one frame per chunk's end (the default)
rem     full    every dumped frame, as the run wrote them (<tag>_play)
rem     blend   the video-rate blend of the dumps (<tag>_play24, maya\play_cache.py blend)
rem   tag names the run: boulders (the finale grid, the default), boulders512, ...
rem   The finale's states cache is written here if missing; another run's caches are
rem   the run's to write (its tag carries its grid and chunking).

setlocal
set "HERE=%~dp0"
set "BIFROST_LIB_CONFIG_FILES=%HERE%pfflip_graphs_config.json"
set "SHOW=%HERE%..\out\collider_showcase"
set "MODE=%~1"
if "%MODE%"=="" set "MODE=states"
set "TAG=%~2"
if "%TAG%"=="" set "TAG=boulders"
set "CACHE=%SHOW%\%TAG%_states"
if /i "%TAG%"=="boulders" set "CACHE=%SHOW%\play_states"
if /i "%MODE%"=="full" set "CACHE=%SHOW%\%TAG%_play"
if /i "%MODE%"=="blend" set "CACHE=%SHOW%\%TAG%_play24"
set "SCENE=%HERE%..\out\play_%TAG%_%MODE%.ma"
set "BOX="
if /i "%TAG%"=="boulders" set "BOX=--box 256 64 128"

if not exist "%CACHE%\cache.json" (
    if /i not "%MODE%_%TAG%"=="states_boulders" (
        echo No cache at %CACHE%: the run writes it.
        exit /b 1
    )
    echo Writing the play cache from the kept states...
    "%HERE%..\.venv\Scripts\python.exe" "%HERE%..\experiments\collider_showcase.py" play || exit /b 1
)
if not exist "%SCENE%" (
    echo Building the Maya scene...
    "C:\Program Files\Autodesk\Maya2027\bin\mayapy.exe" "%HERE%build_particle_play_scene.py" "%CACHE%" "%SCENE%" --obj "%SHOW%\%TAG%.obj" %BOX% || exit /b 1
)
echo Opening Maya...
start "" "C:\Program Files\Autodesk\Maya2027\bin\maya.exe" "%SCENE%"
