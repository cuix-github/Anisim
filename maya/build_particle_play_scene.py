"""Builds a ready-to-open Maya scene around the particle_play compound: a Bifrost graph
whose frame port follows Maya's timeline, so scrubbing or playing shows the liquid
particles of one cache frame per Maya frame, coloured by speed; optional meshes (the
colliders, as OBJ) and the domain as a template box, a camera framing it, the
timeline set to the cache's frames.

  "C:\\Program Files\\Autodesk\\Maya2027\\bin\\mayapy.exe" build_particle_play_scene.py
      <cache_dir> <scene.ma> [--obj <file.obj>] [--box nx ny nz] [--fps 24]

The scene references the compound by name, so Maya must start with
BIFROST_LIB_CONFIG_FILES pointing at pfflip_graphs_config.json (the play_*.bat
launchers do this). In the graph, point_size and speed_scale (the speed drawn white)
are on the particle_play node.
"""
import json
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["BIFROST_LIB_CONFIG_FILES"] = os.path.join(HERE, "pfflip_graphs_config.json")
os.environ.setdefault("MAYA_DISABLE_CER", "1")
os.environ.setdefault("MAYA_DISABLE_CIP", "1")
os.environ.setdefault("MAYA_DISABLE_CLIC_IPM", "1")

FPS_UNITS = {24: "film", 25: "pal", 30: "ntsc", 48: "show", 50: "palf", 60: "ntscf"}


def parse(argv):
    cache, scene, opts = argv[0], argv[1], {"obj": None, "box": None, "fps": 24}
    i = 2
    while i < len(argv):
        if argv[i] == "--obj":
            opts["obj"], i = argv[i + 1], i + 2
        elif argv[i] == "--box":
            opts["box"], i = [float(v) for v in argv[i + 1:i + 4]], i + 4
        elif argv[i] == "--fps":
            opts["fps"], i = int(argv[i + 1]), i + 2
        else:
            raise SystemExit("unknown option " + argv[i])
    return os.path.abspath(cache), os.path.abspath(scene), opts


def aim(cmds, cam, eye, target):
    """Place a Maya camera (which looks down its -Z) at eye, looking at target."""
    dx, dy, dz = (t - e for t, e in zip(target, eye))
    cmds.setAttr(cam + ".translate", *eye)
    cmds.setAttr(cam + ".rotate", math.degrees(math.atan2(dy, math.hypot(dx, dz))),
                 math.degrees(math.atan2(-dx, -dz)), 0.0)
    shape = cmds.listRelatives(cam, shapes=True)[0]
    cmds.setAttr(shape + ".centerOfInterest", math.sqrt(dx * dx + dy * dy + dz * dz))
    cmds.setAttr(shape + ".farClipPlane", 1.0e5)


def main():
    cache, scene, opts = parse(sys.argv[1:])
    with open(os.path.join(cache, "cache.json")) as f:
        frames = json.load(f)["frames"]
    import maya.standalone
    try:
        maya.standalone.initialize(name="python")
        standalone = True
    except Exception:
        standalone = False  # already inside Maya
    import maya.cmds as cmds

    for plug in ("bifrostGraph", "objExport"):
        if not cmds.pluginInfo(plug, q=True, loaded=True):
            cmds.loadPlugin(plug)
    cmds.file(new=True, force=True)
    cmds.currentUnit(time=FPS_UNITS.get(opts["fps"], "film"))

    shape = cmds.createNode("bifrostGraphShape", name="particlePlayShape")
    node = cmds.vnnCompound(shape, "/", addNode="PFFlipGraphs,User::PFFlip,particle_play")[0]
    cmds.vnnNode(shape, "/input", createOutputPort=("frame", "float"))
    cmds.vnnNode(shape, "/output", createInputPort=("points", "Object"))
    cmds.vnnConnect(shape, "/input.frame", "/" + node + ".frame")
    cmds.vnnConnect(shape, "/" + node + ".points", "/output.points")
    for port, name in (("liquid_pattern", "liquid.####"), ("speed_pattern", "speed.####")):
        cmds.vnnNode(shape, "/" + node, setPortDefaultValues=(
            port, os.path.join(cache, name).replace("\\", "/")))
    cmds.connectAttr("time1.outTime", shape + ".frame")

    if opts["obj"]:
        before = set(cmds.ls(assemblies=True))
        cmds.file(os.path.abspath(opts["obj"]), i=True, type="OBJ", ignoreVersion=True,
                  options="mo=1", namespace="obj")
        meshes = sorted(set(cmds.ls(assemblies=True)) - before)
        mat = cmds.shadingNode("lambert", asShader=True, name="colliderGrey")
        cmds.setAttr(mat + ".color", 0.42, 0.38, 0.33, type="double3")
        sg = cmds.sets(renderable=True, noSurfaceShader=True, empty=True, name=mat + "SG")
        cmds.connectAttr(mat + ".outColor", sg + ".surfaceShader")
        cmds.sets(meshes, edit=True, forceElement=sg)
        cmds.group(meshes, name="colliders")
    if opts["box"]:
        nx, ny, nz = opts["box"]
        box = cmds.polyCube(width=nx, height=ny, depth=nz, name="domain")[0]
        cmds.setAttr(box + ".translate", nx / 2, ny / 2, nz / 2)
        cmds.setAttr(box + ".template", 1)
        size = max(nx, ny, nz)
        aim(cmds, "persp", (nx * 1.05, ny * 1.9, nz + 0.75 * size), (nx * 0.55, 0.0, nz * 0.45))

    cmds.playbackOptions(min=1, max=frames, ast=1, aet=frames)
    cmds.currentTime(1)
    os.makedirs(os.path.dirname(scene), exist_ok=True)
    cmds.file(rename=scene)
    cmds.file(save=True, type="mayaAscii")
    print("saved", scene, "-", frames, "frames, compound node", node, flush=True)

    if standalone:
        sys.stdout.flush()
        import ctypes
        ctypes.windll.kernel32.TerminateProcess(ctypes.windll.kernel32.GetCurrentProcess(), 0)


if __name__ == "__main__":
    main()
