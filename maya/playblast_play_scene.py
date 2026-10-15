"""Playblasts a particle play scene from inside Maya's GUI session (maya.exe, not
mayabatch): Maya's Hardware 2.0 batch renderer draws meshes but not the Bifrost
points of the particle_play compound, while the viewport draws them, so the frames
come from the viewport, off screen, at the size asked for, with multisampling on.
Driven by render_play_scene.bat, which sets the environment this reads:

  PLAY_SCENE  the .ma to open        PLAY_OUT    the folder for the PNG frames
  PLAY_NAME   the frames' name        PLAY_FIRST / PLAY_LAST  the frame range
                                      (the scene's timeline when empty)
  PLAY_WIDTH / PLAY_HEIGHT  the frame size (1920 x 1080 when empty)

Maya quits when the frames are written. Run it as
  maya.exe -command "python(\\"import runpy; runpy.run_path(r'<this file>')\\")"
"""
import os

import maya.cmds as cmds


def main():
    scene, out, name = (os.environ[k] for k in ("PLAY_SCENE", "PLAY_OUT", "PLAY_NAME"))
    width = int(os.environ.get("PLAY_WIDTH") or 1920)
    height = int(os.environ.get("PLAY_HEIGHT") or 1080)
    cmds.file(scene, open=True, force=True)
    first = int(os.environ.get("PLAY_FIRST") or cmds.playbackOptions(q=True, min=True))
    last = int(os.environ.get("PLAY_LAST") or cmds.playbackOptions(q=True, max=True))
    cmds.setAttr("hardwareRenderingGlobals.multiSampleEnable", 1)
    cmds.setAttr("hardwareRenderingGlobals.multiSampleCount", 8)
    # a dark, flat background: the fastest particles are drawn white and vanish on
    # a white one (these are session settings; Maya is closed without saving prefs)
    cmds.displayPref(displayGradient=False)
    cmds.displayRGBColor("background", 0.11, 0.11, 0.13)
    panel = cmds.getPanel(withFocus=True)
    if not panel or cmds.getPanel(typeOf=panel) != "modelPanel":
        panel = cmds.getPanel(type="modelPanel")[0]
    cmds.setFocus(panel)
    cmds.modelEditor(panel, e=True, camera="persp", grid=False, headsUpDisplay=False,
                     displayAppearance="smoothShaded", allObjects=True,
                     rendererName="vp2Renderer")
    os.makedirs(out, exist_ok=True)
    cmds.playblast(startTime=first, endTime=last, format="image", compression="png",
                   filename=os.path.join(out, name), widthHeight=(width, height),
                   percent=100, quality=100, viewer=False, showOrnaments=False,
                   offScreen=True, framePadding=4, forceOverwrite=True, clearCache=True)
    with open(os.path.join(out, name + ".playblast.txt"), "w") as f:
        f.write(f"{scene}\nframes {first}..{last} {width}x{height}\n")
    cmds.quit(force=True)


main()
