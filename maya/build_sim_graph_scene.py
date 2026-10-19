"""Builds a Maya scene that holds the simulation graph itself, to look at, not to run:
one Bifrost graph with the compound User::PFFlip::sim_3d_fused (graphs/sim_3d_fused.json,
the fused 3D loop every 3D run uses through bifcmd) as its only node, nothing connected
to its outputs, so Maya never evaluates it. Double-click the node (or Tab into it) to
see the reads, the loop and the writes; inside the loop, the transfer node (p2g_3d) and
the step node (step_3d), the C++ of nodes/src. Run by open_sim_graph_in_maya.bat, which
puts the node pack and maya/pfflip_sim_config.json on BIFROST_LIB_CONFIG_FILES first.

  "C:\\Program Files\\Autodesk\\Maya2027\\bin\\mayapy.exe" build_sim_graph_scene.py <scene.ma>
"""
import os
import sys

os.environ.setdefault("MAYA_DISABLE_CER", "1")
os.environ.setdefault("MAYA_DISABLE_CIP", "1")
os.environ.setdefault("MAYA_DISABLE_CLIC_IPM", "1")


def main():
    scene = os.path.abspath(sys.argv[1])
    import maya.standalone
    maya.standalone.initialize(name="python")
    import maya.cmds as cmds
    if not cmds.pluginInfo("bifrostGraph", q=True, loaded=True):
        cmds.loadPlugin("bifrostGraph")
    cmds.file(new=True, force=True)
    shape = cmds.createNode("bifrostGraphShape", name="simGraphShape")
    node = cmds.vnnCompound(shape, "/", addNode="PFFlipSim,User::PFFlip,sim_3d_fused")[0]
    parent = cmds.listRelatives(shape, parent=True)[0]
    cmds.rename(parent, "sim_3d_fused_graph")
    os.makedirs(os.path.dirname(scene), exist_ok=True)
    cmds.file(rename=scene)
    cmds.file(save=True, type="mayaAscii")
    print("saved", scene, "- node", node, flush=True)
    sys.stdout.flush()
    import ctypes
    ctypes.windll.kernel32.TerminateProcess(ctypes.windll.kernel32.GetCurrentProcess(), 0)


if __name__ == "__main__":
    main()
