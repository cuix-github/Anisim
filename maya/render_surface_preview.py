"""Renders chosen frames of a water-surface scene (build_surface_scene.py) with Arnold,
from inside Maya's GUI session (maya.exe, not mayapy: Arnold batch is a separate
entitlement, rendering in the session is what a Maya subscription licenses). Driven
by render_surface_preview.ps1, which builds the scene first and sets the environment
this reads:

  SURF_SCENE   the .ma to open             SURF_FRAMES  the frames, space separated
  SURF_OUT     the folder for the frames    SURF_NAME    their name (<name>.####.png)
  SURF_WIDTH / SURF_HEIGHT   frame size (1920 x 1080)
  SURF_SAMPLES  Arnold AA samples (4)       SURF_THREADS  render threads (0 = all)
  SURF_GPU      "1" renders on the GPU      SURF_GI_VOLUME_DEPTH  scattering bounces in
  SURF_VOLUME_SAMPLES  volume samples (2)   the mist (3; 0 renders dense mist as smoke)

Each frame goes through Maya's Render Sequence path (arnoldRender -seq), which
writes to <out>/tmp/ and leaves the move to an idle callback, so the file is moved
here; an Arnold log per frame sits beside it. A <name>.done.txt marker is written at
the end, then Maya quits (closed by the launcher if that does not take).
"""
import glob
import os
import shutil
import time
import traceback

import maya.cmds as cmds

E = os.environ
out, name = E["SURF_OUT"], E.get("SURF_NAME", "surface")
os.makedirs(out, exist_ok=True)


def note(msg):
    with open(os.path.join(out, name + ".notes.txt"), "a") as f:
        f.write(msg + "\n")


def env(key, default, cast=float):
    v = E.get(key, "")
    return cast(v) if v else default


def render(frames):
    width, height = env("SURF_WIDTH", 1920, int), env("SURF_HEIGHT", 1080, int)
    cmds.setAttr("defaultRenderGlobals.currentRenderer", "arnold", type="string")
    cmds.setAttr("defaultRenderGlobals.animation", 1)
    cmds.setAttr("defaultRenderGlobals.extensionPadding", 4)
    cmds.setAttr("defaultRenderGlobals.outFormatControl", 0)
    cmds.setAttr("defaultRenderGlobals.putFrameBeforeExt", 1)
    cmds.setAttr("defaultRenderGlobals.imageFilePrefix", name, type="string")
    cmds.setAttr("defaultArnoldDriver.ai_translator", "png", type="string")
    # PNG is a display image. Retaining the driver's EXR default (Use Output
    # Transform) with outputTransformEnabled=False writes linear ACEScg values
    # into an ordinary PNG. Keep the scene's selected view transform explicit.
    color_modes = cmds.attributeQuery("colorManagement", node="defaultArnoldDriver", listEnum=True)[0].split(":")
    cmds.setAttr("defaultArnoldDriver.colorManagement", color_modes.index("Use View Transform"))
    cmds.setAttr("defaultResolution.width", width)
    cmds.setAttr("defaultResolution.height", height)
    cmds.setAttr("defaultResolution.deviceAspectRatio", width / float(height))
    cmds.setAttr("defaultArnoldRenderOptions.AASamples", env("SURF_SAMPLES", 4, int))
    # light must scatter several times inside the mist or a dense cloud renders as
    # dark smoke (Maya's default volume GI depth is 0: single scattering only)
    cmds.setAttr("defaultArnoldRenderOptions.GIVolumeDepth", env("SURF_GI_VOLUME_DEPTH", 3, int))
    cmds.setAttr("defaultArnoldRenderOptions.GIVolumeSamples", env("SURF_VOLUME_SAMPLES", 2, int))
    threads = env("SURF_THREADS", 0, int)
    cmds.setAttr("defaultArnoldRenderOptions.threads_autodetect", 0 if threads else 1)
    if threads:
        cmds.setAttr("defaultArnoldRenderOptions.threads", threads)
    if E.get("SURF_GPU") == "1":
        cmds.setAttr("defaultArnoldRenderOptions.renderDevice", 1)
    cmds.setAttr("defaultArnoldRenderOptions.log_verbosity", 2)
    cmds.setAttr("defaultArnoldRenderOptions.log_to_file", 1)
    cmds.setAttr("defaultArnoldRenderOptions.log_to_console", 0)
    cmds.workspace(fileRule=["images", out.replace("\\", "/")])
    tmp = os.path.join(out, "tmp")
    # Arnold names the log <name>.<frame>.log by itself
    cmds.setAttr("defaultArnoldRenderOptions.log_filename",
                 os.path.join(out, name + ".log").replace("\\", "/"), type="string")
    for f in frames:
        cmds.setAttr("defaultRenderGlobals.startFrame", f)
        cmds.setAttr("defaultRenderGlobals.endFrame", f)
        cmds.currentTime(f)
        t0 = time.time()
        cmds.arnoldRender(seq=str(f), w=width, h=height, cam="persp", srv=False)
        wall = time.time() - t0
        written = sorted(glob.glob(os.path.join(tmp, "*.png")), key=os.path.getmtime)
        final = os.path.join(out, f"{name}.{f:04d}.png")
        if written:
            shutil.move(written[-1], final)
            note(f"frame {f}: {wall:.1f} s -> {os.path.basename(final)}")
        else:
            note(f"frame {f}: {wall:.1f} s, NO IMAGE in {tmp}")


def main():
    frames = [int(v) for v in E["SURF_FRAMES"].split()]
    t0 = time.time()
    try:
        cmds.scriptEditorInfo(writeHistory=True,
                              historyFilename=os.path.join(out, name + ".maya.log"))
        cmds.file(E["SURF_SCENE"], open=True, force=True)
        note(f"scene opened in {time.time() - t0:.1f} s")
        render(frames)
    except Exception:
        note("ERROR\n" + traceback.format_exc())
    finally:
        with open(os.path.join(out, name + ".done.txt"), "w") as f:
            f.write(f"frames {frames}\n{time.time() - t0:.0f} s\n")
        cmds.quit(force=True)


main()
