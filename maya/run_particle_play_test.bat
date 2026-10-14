@echo off
rem Runs particle_play_test.json headlessly: the playback compound at one frame of a
rem play cache, its points' positions, colours and sizes written back as .npy.
rem Stock nodes only, so no pack is needed. Driven by gates/validate_particle_play.py.
rem
rem Usage: run_particle_play_test.bat --set-port liquid_pattern <path> ... (any bifcmd args)

setlocal
set "BIFROST=C:\Program Files\Autodesk\Bifrost\Maya2027\3.1.0.8\bifrost"
set "PATH=%BIFROST%\bin;%BIFROST%\thirdparty\bin;%PATH%"
bifcmd "%~dp0particle_play_test.json" %*
