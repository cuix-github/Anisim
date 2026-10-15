@echo off
rem Runs particle_surface_test.json headlessly: the surface compound at one frame of a
rem play cache, the mesh's vertex positions written back as .npy. Stock nodes only,
rem so no pack is needed. Driven by gates/validate_particle_surface.py.
rem
rem Usage: run_particle_surface_test.bat --set-port liquid_pattern <path> ... (any bifcmd args)

setlocal
set "BIFROST=C:\Program Files\Autodesk\Bifrost\Maya2027\3.1.0.8\bifrost"
set "PATH=%BIFROST%\bin;%BIFROST%\thirdparty\bin;%PATH%"
bifcmd "%~dp0particle_surface_test.json" %*
