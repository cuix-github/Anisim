"""One click from a run to Maya. Finds the run's play cache by its name (graphs/runs.py,
either layout), writes the states cache if it is missing, builds the Maya scene if it
is missing or points at another cache, then opens Maya on it with the PF-FLIP compound
library on BIFROST_LIB_CONFIG_FILES. Scrub or play: each Maya frame is one cache frame,
the liquid particles coloured by speed, with the colliders and the domain box (the grid,
from the cache's index). play_showcase_in_maya.bat wraps it.

  play_showcase.py [states | full | blend] [tag] [--no-open] [--rebuild]
    states   the run's kept states, one frame per chunk's end (the default)
    full     every dumped frame, as the run wrote them (cache/full)
    blend    the video-rate blend of the dumps (cache/blend24;
             experiments/collider_showcase.py blend --tag <tag> writes it)
    tag      the run: boulders (the finale grid, the default), boulders512, ...
    --no-open  build everything and print the scene's path, without starting Maya
    --rebuild  build the scene again even if it is there
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
from runs import RunPaths, prepare  # noqa: E402

MAYA_BIN = r"C:\Program Files\Autodesk\Maya2027\bin"
MODES = {"states": "states", "full": "full", "blend": "blend24"}


def main(argv):
    opts = {a for a in argv if a.startswith("--")}
    words = [a for a in argv if not a.startswith("--")]
    mode = words[0] if words else "states"
    tag = words[1] if len(words) > 1 else "boulders"
    if mode not in MODES or opts - {"--no-open", "--rebuild"}:
        raise SystemExit(__doc__)
    os.environ["BIFROST_LIB_CONFIG_FILES"] = os.path.join(HERE, "pfflip_graphs_config.json")
    paths = RunPaths(tag)
    cache = os.path.abspath(paths.cache(MODES[mode]))
    index = os.path.join(cache, "cache.json")
    if not os.path.exists(index):
        if mode != "states":
            raise SystemExit(f"No cache at {cache}.\n  full: the run writes it as it goes "
                             f"(collider_showcase.py stage N --tag {tag})\n  blend: "
                             f".venv/Scripts/python.exe experiments/collider_showcase.py blend --tag {tag}")
        print("Writing the play cache from the kept states...", flush=True)
        subprocess.run([sys.executable, os.path.join(HERE, "..", "experiments", "collider_showcase.py"),
                        "play", "--tag", tag], check=True)
    with open(index) as f:
        info = json.load(f)
    scene = os.path.abspath(paths.scene(mode))
    stale = False
    if os.path.exists(scene):          # the scene names its cache: another cache means rebuild
        with open(scene, encoding="utf-8", errors="ignore") as f:
            stale = cache.replace("\\", "/").lower() not in f.read().lower()
    if "--rebuild" in opts or stale or not os.path.exists(scene):
        args = [os.path.join(MAYA_BIN, "mayapy.exe"), os.path.join(HERE, "build_particle_play_scene.py"),
                cache, prepare(scene)]
        obj = paths.input("colliders.obj")
        if os.path.exists(obj):
            args += ["--obj", os.path.abspath(obj)]
        if not info.get("grid") and tag == "boulders":
            args += ["--box", "256", "64", "128"]
        print("Building the Maya scene..." + (" (it pointed at another cache)" if stale else ""), flush=True)
        subprocess.run(args, check=True)
    if "--no-open" in opts:
        print(scene)
        return
    print("Opening Maya...", flush=True)
    subprocess.Popen([os.path.join(MAYA_BIN, "maya.exe"), scene], close_fds=True)


if __name__ == "__main__":
    main(sys.argv[1:])
