@echo off
rem One click from a run to Maya: play_showcase.py finds the run's cache by name
rem (graphs/runs.py), writes the states cache if missing, builds the scene if missing,
rem then opens Maya on it with the PF-FLIP compound library on BIFROST_LIB_CONFIG_FILES.
rem
rem   play_showcase_in_maya.bat [states / full / blend] [tag] [--no-open] [--rebuild]
rem     states  the run's kept states, one frame per chunk's end (the default)
rem     full    every dumped frame, as the run wrote them
rem     blend   the video-rate blend of the dumps (collider_showcase.py blend --tag <tag>)
rem   tag names the run: boulders (the finale grid, the default), boulders512, ...

setlocal
set "HERE=%~dp0"
set "BIFROST_LIB_CONFIG_FILES=%HERE%pfflip_graphs_config.json"
"%HERE%..\.venv\Scripts\python.exe" "%HERE%play_showcase.py" %*
exit /b %ERRORLEVEL%
