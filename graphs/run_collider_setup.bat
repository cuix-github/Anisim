@echo off
rem Runs the mesh-collider setup graph (collider_setup.json) headlessly: meshes from
rem .npy through the shipped collider node and our collider_field adapter, out as
rem solid_phi.npy for the loop graphs. Stock nodes only, so no pack is needed.
rem Driven by graphs/colliders.py.
rem
rem Usage: run_collider_setup.bat --set-port path_points <path> ... (any bifcmd args)

setlocal
set "BIFROST=C:\Program Files\Autodesk\Bifrost\Maya2027\3.1.0.8\bifrost"
set "PATH=%BIFROST%\bin;%BIFROST%\thirdparty\bin;%PATH%"
bifcmd "%~dp0collider_setup.json" %*
