@echo off
rem Runs particle_mist_test.json headlessly: the mist compound at one frame of a play
rem cache, the mist positions and every particle's depth written back as .npy.
rem Stock nodes only, so no pack is needed. Driven by gates/validate_particle_surface.py.
rem
rem Usage: run_particle_mist_test.bat --set-port liquid_pattern <path> ... (any bifcmd args)

setlocal
set "BIFROST=C:\Program Files\Autodesk\Bifrost\Maya2027\3.1.0.8\bifrost"
set "PATH=%BIFROST%\bin;%BIFROST%\thirdparty\bin;%PATH%"
bifcmd "%~dp0particle_mist_test.json" %*
