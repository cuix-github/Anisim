@echo off
rem Renders a particle play scene to 1920x1080 PNG frames and an mp4: an off-screen
rem viewport playblast from a GUI Maya session, because Maya's Hardware 2.0 batch
rem renderer does not draw the Bifrost points and Arnold batch is not licensed here.
rem The work is in render_play_scene.ps1 (this is the one-click wrapper).
rem
rem   render_play_scene.bat <scene.ma> <out_dir> <name> [first] [last] [fps]
rem     frames land in <out_dir>\<name>.####.png, the video in <out_dir>\<name>.mp4;
rem     first/last default to the scene's timeline, fps to 24.

setlocal
set "HERE=%~dp0"
set "NAME=%~3"
if "%NAME%"=="" set "NAME=play"
set "FPS=%~6"
if "%FPS%"=="" set "FPS=24"
powershell -NoProfile -ExecutionPolicy Bypass -File "%HERE%render_play_scene.ps1" -Scene "%~1" -Out "%~2" -Name "%NAME%" -First "%~4" -Last "%~5" -Fps %FPS%
exit /b %ERRORLEVEL%
