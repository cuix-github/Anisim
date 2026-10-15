@echo off
rem Renders chosen frames of a play cache as an Arnold water surface (glass mixed into
rem foam by curvature, sky and sun, the colliders), from a minimized GUI Maya session.
rem The work is in render_surface_preview.ps1 (this is the one-click wrapper).
rem
rem   render_surface_preview.bat <cache_dir> <out_dir> <name> "<frames>" [obj] ["nx ny nz"]
rem     e.g. render_surface_preview.bat out\collider_showcase\boulders512_play24 ^
rem            out\collider_showcase\boulders512_surface boulders512 "52 194" ^
rem            out\collider_showcase\boulders512.obj
rem     obj is "-" for none; the grid comes from the cache's index unless given;
rem     frame size and the surface knobs come from SURF_* environment variables
rem     (see the .ps1).

setlocal
set "HERE=%~dp0"
if "%~4"=="" (
    echo usage: render_surface_preview.bat ^<cache_dir^> ^<out_dir^> ^<name^> "^<frames^>" [obj ^| -] ["nx ny nz"]
    exit /b 2
)
set "OBJ=%~5"
if "%OBJ%"=="-" set "OBJ="
powershell -NoProfile -ExecutionPolicy Bypass -File "%HERE%render_surface_preview.ps1" -Cache "%~1" -Out "%~2" -Name "%~3" -Frames "%~4" -Obj "%OBJ%" -Grid "%~6"
exit /b %ERRORLEVEL%
