"""Builds the Arnold water scene around the particle_surface and particle_mist
compounds, headless (mayapy): two Bifrost graphs whose frame ports follow the
timeline, so each Maya frame is one play-cache frame, the bulk meshed as a water
surface and the particles that surface leaves out (spray, mist, one-particle
sheets) splatted into a fog volume; a water shader (glass at the index of water,
tinted with depth) mixed by convex curvature into a white foam shader, the
render-time rule the PF-FLIP and ST-FLIP papers use; a white scattering volume
shader on the mist; a sky dome (the background) and a warm sun; the colliders in
matte grey; a ground plane; the camera framing the domain; Arnold set as the
renderer. Reads the environment render_surface_preview.ps1 sets (SURF_CACHE,
SURF_SCENE, SURF_GRID, SURF_OBJ, SURF_CAMERA, SURF_POINT_SIZE, SURF_DETAIL,
SURF_SMOOTH_ITERS, SURF_MESH_SMOOTH, SURF_MIST_RADIUS, SURF_MIST_DETAIL,
SURF_MIST_DENSITY, SURF_MIST_THRESHOLD, SURF_METRES_PER_CELL, SURF_DEFAULT_RADIUS,
SURF_MCF_ITERS, SURF_MASK_SIGMA, SURF_MASK_THETA, SURF_MASK_ZETA, SURF_MASK_FLOOR,
SURF_MASK_DILATE,
SURF_LOOK, SURF_WATER_DEPTH,
SURF_SKY_INTENSITY, SURF_SKY_VISIBLE, SURF_SUN_INTENSITY, SURF_SUN_ROTATE, SURF_WATER_COLOR,
SURF_CAMERA_PATH (a camera.json: {"frames": [[ex, ey, ez, tx, ty, tz], ...]} in cells, frame
f the list's entry f - 1, keyed on the camera so it tracks; the fixed SURF_CAMERA otherwise);
see the .ps1) and saves SURF_SCENE. Built
here rather than in the GUI session
because a script run at GUI startup sees the compound library before Bifrost's own
libraries are loaded (its ports come back as unknown types); a saved scene opens
fine there.

  "C:\\Program Files\\Autodesk\\Maya2027\\bin\\mayapy.exe" build_surface_scene.py
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
os.environ["BIFROST_LIB_CONFIG_FILES"] = os.path.join(HERE, "pfflip_graphs_config.json")
os.environ.setdefault("MAYA_DISABLE_CER", "1")
os.environ.setdefault("MAYA_DISABLE_CIP", "1")
os.environ.setdefault("MAYA_DISABLE_CLIC_IPM", "1")
sys.path.insert(0, HERE)
from build_particle_play_scene import aim  # noqa: E402

E = os.environ


def env(key, default, cast=float):
    v = E.get(key, "")
    return cast(v) if v else default


def shader(cmds, kind, name, **attrs):
    s = cmds.shadingNode(kind, asShader=True, name=name)
    for k, v in attrs.items():
        if isinstance(v, tuple):
            cmds.setAttr(s + "." + k, *v, type="double3")
        else:
            cmds.setAttr(s + "." + k, v)
    return s


def assign(cmds, shader_node, objects, name):
    sg = cmds.sets(renderable=True, noSurfaceShader=True, empty=True, name=name)
    cmds.connectAttr(shader_node + ".outColor", sg + ".surfaceShader")
    cmds.sets(objects, e=True, forceElement=sg)
    return sg


def graph(cmds, name, compound, out_port, ports):
    """A Bifrost graph shape around one compound: its frame from the timeline, its
    out_port the shape's output, the given ports set."""
    shape = cmds.createNode("bifrostGraphShape", name=name)
    node = cmds.vnnCompound(shape, "/", addNode="PFFlipGraphs,User::PFFlip," + compound)[0]
    cmds.vnnNode(shape, "/input", createOutputPort=("frame", "float"))
    cmds.vnnNode(shape, "/output", createInputPort=(out_port, "Object"))
    cmds.vnnConnect(shape, "/input.frame", "/" + node + ".frame")
    cmds.vnnConnect(shape, "/" + node + "." + out_port, "/output." + out_port)
    for port, value in ports:
        cmds.vnnNode(shape, "/" + node, setPortDefaultValues=(port, str(value)))
    cmds.connectAttr("time1.outTime", shape + ".frame")
    return shape


# the spray channel's regularisation (see particle_mist in generate_particle_surface.py:
# the stored density is the kernel sum over this, the shader multiplies it back)
SPRAY_EPS = 10000.0


def radius_pattern(cache, out):
    """The cache's radius.#### when the run dumped droplet radii; otherwise zeros
    written beside the scene (out/radius/), one file per liquid frame, so the mist
    draws every particle at default_radius."""
    if os.path.exists(os.path.join(cache, "radius.0001.npy")):
        return os.path.join(cache, "radius.####").replace("\\", "/")
    d = os.path.join(out, "radius")
    os.makedirs(d, exist_ok=True)
    import numpy as np
    for name in sorted(os.listdir(cache)):
        if name.startswith("liquid.") and name.endswith(".npy"):
            n = np.load(os.path.join(cache, name), mmap_mode="r").shape[0]
            np.save(os.path.join(d, "radius." + name[7:]), np.zeros(n, np.float32))
    return os.path.join(d, "radius.####").replace("\\", "/")


def build(cmds):
    cache = E["SURF_CACHE"]
    nx, ny, nz = [float(v) for v in E["SURF_GRID"].split()]
    for plug in ("bifrostGraph", "mtoa", "objExport"):
        if not cmds.pluginInfo(plug, q=True, loaded=True):
            cmds.loadPlugin(plug)
    import mtoa.core as core
    cmds.file(new=True, force=True)
    core.createOptions()
    cmds.currentUnit(time="film")

    # the surface and the mist: each compound's frame follows the timeline
    pattern = os.path.join(cache, "liquid.####").replace("\\", "/")
    radius, detail = env("SURF_POINT_SIZE", 0.5), env("SURF_DETAIL", 0.5)
    smooth_iters = env("SURF_SMOOTH_ITERS", 1, int)
    # the escaped droplets' radii: out of the surface, always mist (the papers' split)
    radii = radius_pattern(cache, E["SURF_OUT"])
    # after the Gaussian passes: ST-FLIP's feature-preserving mean-curvature flow
    # (0 = off; the paper's 30), its mask's blur in cells, thresholds and floor
    mcf = (("mcf_iterations", env("SURF_MCF_ITERS", 0, int)),
           ("mask_sigma", env("SURF_MASK_SIGMA", 2.0)), ("mask_theta", env("SURF_MASK_THETA", 2.0)),
           ("mask_zeta", env("SURF_MASK_ZETA", 5.0)), ("mask_floor", env("SURF_MASK_FLOOR", 0.0)),
           ("mask_dilate", env("SURF_MASK_DILATE", 1.0)))
    shape = graph(cmds, "particleSurfaceShape", "particle_surface", "mesh", (
        ("liquid_pattern", pattern), ("point_size", radius), ("detail_size", detail),
        ("smooth_iterations", smooth_iters), ("smoothing", env("SURF_MESH_SMOOTH", 0.1)),
        ("radius_pattern", radii)) + mcf)
    # the mist: each particle the surface leaves out drawn as its cloud of droplets
    # (its radius from the cache, default_radius without one), additive density on a
    # spray grid of SURF_MIST_DETAIL cells, clouds of SURF_MIST_RADIUS cells,
    # SURF_MIST_DENSITY the fraction of a parcel's water drawn (1 = physical)
    mist = graph(cmds, "particleMistShape", "particle_mist", "mist", (
        ("liquid_pattern", pattern), ("radius_pattern", radii),
        ("surface_radius", radius), ("surface_detail", detail),
        ("smooth_iterations", smooth_iters),
        ("mist_threshold", env("SURF_MIST_THRESHOLD", 0.0)),
        ("spray_detail", env("SURF_MIST_DETAIL", 0.5)),
        ("cloud_radius", env("SURF_MIST_RADIUS", 1.0)),
        ("metres_per_cell", env("SURF_METRES_PER_CELL", 1.0)),
        ("default_radius", env("SURF_DEFAULT_RADIUS", 0.0001)),
        ("mist_density", env("SURF_MIST_DENSITY", 1.0)),
        ("spray_epsilon", SPRAY_EPS)) + mcf)

    # water: glass at the index of water, tinted with depth; foam: white, matte;
    # mixed by convex curvature (the papers render high-curvature areas foam-like)
    # SURF_LOOK=ocean: the breaking-wave stills of the PF-FLIP paper, dark water that
    # absorbs within metres, a near-black sky, the sun low behind the spray (backlit,
    # toward a camera on the -z side), a dark floor; the default look is unchanged
    ocean = E.get("SURF_LOOK", "") == "ocean"
    water = shader(cmds, "aiStandardSurface", "water", base=0.0, specular=1.0,
                   specularRoughness=0.03, specularIOR=1.33, transmission=1.0,
                   transmissionColor=tuple(float(v) for v in E.get("SURF_WATER_COLOR", "").split())
                   or ((0.22, 0.42, 0.46) if ocean else (0.7, 0.88, 0.94)),
                   transmissionDepth=env("SURF_WATER_DEPTH", 6.0 if ocean else 60.0))
    foam = shader(cmds, "aiStandardSurface", "foam", base=1.0, baseColor=(0.94, 0.96, 0.98),
                  specular=0.25, specularRoughness=0.55)
    curv = cmds.shadingNode("aiCurvature", asTexture=True, name="foamCurvature")
    for k, v in (("radius", 1.2), ("spread", 1.0), ("threshold", 0.15), ("bias", 0.5),
                 ("multiply", 1.6)):
        cmds.setAttr(curv + "." + k, v)
    mix = cmds.shadingNode("aiMixShader", asShader=True, name="waterFoamMix")
    cmds.connectAttr(water + ".outColor", mix + ".shader1")
    cmds.connectAttr(foam + ".outColor", mix + ".shader2")
    cmds.connectAttr(curv + ".outColorR", mix + ".mix")
    assign(cmds, mix, shape, "waterSG")
    # mist: a white scattering fog read from the spray density channel (its stored
    # value is the extinction per cell over SPRAY_EPS)
    fog = cmds.shadingNode("aiStandardVolume", asShader=True, name="mistFog")
    cmds.setAttr(fog + ".densityChannel", "voxel_fog_density", type="string")
    cmds.setAttr(fog + ".density", SPRAY_EPS)
    cmds.setAttr(fog + ".scatter", 1.0)
    cmds.setAttr(fog + ".scatterColor", 1.0, 1.0, 1.0, type="double3")
    cmds.setAttr(fog + ".scatterAnisotropy", 0.4)
    msg = cmds.sets(renderable=True, noSurfaceShader=True, empty=True, name="mistSG")
    cmds.connectAttr(fog + ".outColor", msg + ".volumeShader")
    cmds.sets(mist, e=True, forceElement=msg)

    # the colliders and the ground
    if E.get("SURF_OBJ"):
        before = set(cmds.ls(assemblies=True))
        cmds.file(os.path.abspath(E["SURF_OBJ"]), i=True, type="OBJ", ignoreVersion=True,
                  options="mo=1", namespace="obj")
        meshes = sorted(set(cmds.ls(assemblies=True)) - before)
        rock = shader(cmds, "aiStandardSurface", "rock", base=1.0,
                      baseColor=(0.08, 0.075, 0.07) if ocean else (0.36, 0.32, 0.28),
                      specular=0.2, specularRoughness=0.7)
        assign(cmds, rock, meshes, "rockSG")
        cmds.group(meshes, name="colliders")
    ground = cmds.polyPlane(w=nx * 6, h=nz * 6, sx=1, sy=1, name="ground")[0]
    cmds.setAttr(ground + ".translate", nx / 2, -1.0, nz / 2)
    sand = shader(cmds, "aiStandardSurface", "sand", base=1.0,
                  baseColor=(0.04, 0.04, 0.045) if ocean else (0.52, 0.47, 0.4),
                  specular=0.1, specularRoughness=0.8)
    assign(cmds, sand, ground, "sandSG")

    # light: a sky dome (seen by the camera as the background) and a warm sun
    sky = cmds.shadingNode("aiSkyDomeLight", asLight=True, name="sky")
    if ocean:
        cmds.setAttr(sky + ".color", 0.03, 0.035, 0.045, type="double3")
        cmds.setAttr(sky + ".intensity", env("SURF_SKY_INTENSITY", 1.0))
        # the sky lights and reflects but the camera sees black behind the water
        cmds.setAttr(sky + ".camera", env("SURF_SKY_VISIBLE", 0.0))
    else:
        cmds.setAttr(sky + ".color", 0.5, 0.66, 0.9, type="double3")
        cmds.setAttr(sky + ".intensity", env("SURF_SKY_INTENSITY", 1.0))
    sun = cmds.directionalLight(name="sun", intensity=env("SURF_SUN_INTENSITY", 4.0 if ocean else 3.0))
    rot = [float(v) for v in E.get("SURF_SUN_ROTATE", "").split()] or ([-18.0, 15.0, 0.0] if ocean
                                                                      else [-48.0, 28.0, 0.0])
    cmds.setAttr(cmds.listRelatives(sun, parent=True)[0] + ".rotate", *rot)
    cmds.setAttr(sun + ".color", 1.0, 0.95, 0.86, type="double3")
    cmds.setAttr(sun + ".aiAngle", 1.5)

    # the camera: the play scene's framing unless given; a camera path keys it per
    # frame (the yaw unwrapped so the keys never spin the long way round)
    size = max(nx, ny, nz)
    cam = E.get("SURF_CAMERA", "").split()
    if len(cam) == 6:
        eye, target = tuple(float(v) for v in cam[:3]), tuple(float(v) for v in cam[3:])
    else:
        eye, target = (nx * 1.05, ny * 1.9, nz + 0.75 * size), (nx * 0.55, 0.0, nz * 0.45)
    aim(cmds, "persp", eye, target)
    if E.get("SURF_CAMERA_PATH"):
        import json
        import math
        with open(E["SURF_CAMERA_PATH"]) as fh:
            frames = json.load(fh)["frames"]
        yaw_prev = None
        for f, row in enumerate(frames, start=1):
            e, t = row[:3], row[3:6]
            aim(cmds, "persp", e, t)
            yaw = cmds.getAttr("persp.rotateY")
            if yaw_prev is not None:
                while yaw - yaw_prev > 180.0:
                    yaw -= 360.0
                while yaw - yaw_prev < -180.0:
                    yaw += 360.0
                cmds.setAttr("persp.rotateY", yaw)
            yaw_prev = yaw
            cmds.setKeyframe("persp", attribute=["translateX", "translateY", "translateZ",
                                                 "rotateX", "rotateY", "rotateZ"], time=f,
                             inTangentType="linear", outTangentType="linear")
        print(f"camera keyed over {len(frames)} frames from {E['SURF_CAMERA_PATH']}", flush=True)
    for c in cmds.ls(cameras=True):
        cmds.setAttr(c + ".renderable", c == "perspShape")
    cmds.setAttr("defaultRenderGlobals.currentRenderer", "arnold", type="string")
    return shape, mist


def main():
    scene = os.path.abspath(E["SURF_SCENE"])
    import maya.standalone
    maya.standalone.initialize(name="python")
    import maya.cmds as cmds
    node = build(cmds)
    os.makedirs(os.path.dirname(scene), exist_ok=True)
    cmds.file(rename=scene)
    cmds.file(save=True, type="mayaAscii")
    print("saved", scene, "- graphs", node, flush=True)
    sys.stdout.flush()
    import ctypes
    ctypes.windll.kernel32.TerminateProcess(ctypes.windll.kernel32.GetCurrentProcess(), 0)


if __name__ == "__main__":
    main()
