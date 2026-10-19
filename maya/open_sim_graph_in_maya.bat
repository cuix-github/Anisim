@echo off
rem Opens Maya on the simulation graph itself (graphs/sim_3d_fused.json), to look at:
rem builds the node pack if it is missing, builds out\sim_graph\sim_graph.ma if it is
rem missing (or with --rebuild), and starts Maya with the node pack and the graph's
rem library on BIFROST_LIB_CONFIG_FILES. In Maya: Windows > Bifrost Graph Editor, select
rem sim_3d_fused_graph, double-click the node to go inside. Nothing is connected to its
rem outputs, so it never runs.

setlocal
set "HERE=%~dp0"
set "PACK=%HERE%..\nodes\build\PFFlipSpike-1.0.0\PFFlipSpikePackConfig.json"
if not exist "%PACK%" (
    call "%HERE%..\nodes\build.bat" || exit /b 1
)
set "BIFROST_LIB_CONFIG_FILES=%PACK%;%HERE%pfflip_sim_config.json"
set "SCENE=%HERE%..\out\sim_graph\sim_graph.ma"
if "%~1"=="--rebuild" if exist "%SCENE%" del "%SCENE%"
if not exist "%SCENE%" (
    echo Building the Maya scene...
    "C:\Program Files\Autodesk\Maya2027\bin\mayapy.exe" "%HERE%build_sim_graph_scene.py" "%SCENE%" || exit /b 1
)
echo Opening Maya...
start "" "C:\Program Files\Autodesk\Maya2027\bin\maya.exe" "%SCENE%"
