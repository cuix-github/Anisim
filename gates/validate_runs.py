"""The run layout (graphs/runs.py): where a run's files go, how the runs made before it
are found where they lie, what the report says is on disk, and that migrating an older
run is a set of renames that loses nothing. Pure Python on a scratch out/ tree it builds
under out/runs_gate/fake_out (its own, rebuilt every time); about a second.

  new      an empty root: every path of a run is under runs/<tag>/
  older    a run laid out the old way (loop3_<tag>_c<N>, collider_showcase/<tag>_*,
           play_<tag>_<mode>.ma): every accessor returns the old path, a chunk that
           does not exist yet goes to the new layout, the finale's play_states is found
  report   discover finds both, describe sizes each part to the byte, the scratch
           folders are listed as everything else, the report prints
  migrate  the plan is exactly the expected renames, folders before the files inside
           them; applied, the bytes and the files are all there under runs/<tag>/, the
           old places are empty, migration.json records every move, a second plan is
           empty; a destination that already exists makes apply refuse before moving

Run from the repo root:  .venv\\Scripts\\python.exe gates\\validate_runs.py
"""
import contextlib
import io
import json
import os
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
import runs  # noqa: E402
from runs import RunPaths, SHOWCASE  # noqa: E402

OUT = os.path.join(HERE, "..", "out", "runs_gate")
FAKE = os.path.join(OUT, "fake_out")
checks = []


def check(name, cond, detail=""):
    checks.append(bool(cond))
    print(f"  {'ok ' if cond else 'FAIL'} {name}" + (f"  ({detail})" if detail and not cond else ""))


def write(path, nbytes=0, text=None):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w" if text is not None else "wb") as f:
        f.write(text if text is not None else b"x" * nbytes)
    return os.path.getsize(path)


def listing(root):
    """Every file under root as (relative path, size)."""
    out = []
    for d, _, fs in os.walk(root):
        for f in fs:
            p = os.path.join(d, f)
            out.append((os.path.relpath(p, root).replace("\\", "/"), os.path.getsize(p)))
    return sorted(out)


def same(a, b):
    return os.path.normpath(a) == os.path.normpath(b)


def older_run(out, tag):
    """A run as the pre-2026-10-16 code left it; returns the bytes written per part."""
    sc = os.path.join(out, SHOWCASE)
    b = {}
    b["dumps"] = (write(os.path.join(out, f"loop3_{tag}_c0", "frame.0000.npy"), 4000)
                  + write(os.path.join(out, f"loop3_{tag}_c0", "positions.npy"), 800)
                  + write(os.path.join(out, f"loop3_{tag}_c1", "final.npy"), 1200))
    b["blend24"] = (write(os.path.join(sc, f"{tag}_play24", "cache.json"), text=json.dumps(
        {"frames": 2, "fps": 24, "blended_from": {"dirs": [f"loop3_{tag}_c0", f"loop3_{tag}_c1"]}}))
        + write(os.path.join(sc, f"{tag}_play24", "liquid.0001.npy"), 400)
        + write(os.path.join(sc, f"{tag}_play24", "speed.0001.npy"), 200))
    b["states"] = (write(os.path.join(sc, f"{tag}_states", "cache.json"), text=json.dumps(
        {"frames": 1, "source": [f"loop3_{tag}_c1"]}))
        + write(os.path.join(sc, f"{tag}_states", "liquid.0001.npy"), 400))
    for n in (f"{tag}_c0.json", f"{tag}_c1.json", f"{tag}_summary.json", f"field_{tag}.json"):
        write(os.path.join(sc, n), text="{}")
    b["inputs"] = (write(os.path.join(sc, f"solid_phi_{tag}.npy"), 2000)
                   + write(os.path.join(sc, f"{tag}.obj"), 300)
                   + write(os.path.join(out, f"collider_{tag}_ab", "solid_phi.npy"), 1000)
                   + write(os.path.join(out, f"collider_{tag}_c", "solid_phi.npy"), 1000))
    b["renders"] = (write(os.path.join(sc, f"{tag}_png", "fr0000.png"), 500)
                    + write(os.path.join(sc, f"{tag}.mp4"), 900)
                    + write(os.path.join(sc, f"{tag}_mist", "mist.0001.png"), 700)
                    + write(os.path.join(sc, f"{tag}_render", "play.0001.png"), 600))
    write(os.path.join(sc, f"{tag}_mist_launch.txt"), text="launch")
    write(os.path.join(out, f"play_{tag}_blend.ma"), text="scene")
    return b


def main():
    shutil.rmtree(FAKE, ignore_errors=True)
    os.makedirs(FAKE)
    out = FAKE

    print("=== new: an empty root")
    rp = RunPaths("fresh", out)
    want = {rp.dumps(0): "runs/fresh/dumps/c0", rp.cache("blend24"): "runs/fresh/cache/blend24",
            rp.cache("states"): "runs/fresh/cache/states", rp.receipt("c0.json"): "runs/fresh/receipts/c0.json",
            rp.receipt("summary.json"): "runs/fresh/receipts/summary.json",
            rp.input("solid_phi.npy"): "runs/fresh/inputs/solid_phi.npy",
            rp.input("collider_ab"): "runs/fresh/inputs/collider_ab",
            rp.renders("audit"): "runs/fresh/renders/audit",
            rp.video("audit"): "runs/fresh/renders/audit/fresh.mp4",
            rp.scene("blend"): "runs/fresh/scenes/play_blend.ma", rp.logs: "runs/fresh/logs"}
    for got, rel in want.items():
        check(rel, same(got, os.path.join(out, rel)), got)
    check("chunks() of nothing is []", rp.chunks() == [])

    print("=== older: a run in the pre-2026-10-16 layout")
    b = older_run(out, "old")
    sc = os.path.join(out, SHOWCASE)
    rp = RunPaths("old", out)
    old = {rp.dumps(0): "loop3_old_c0", rp.dumps(1): "loop3_old_c1",
           rp.cache("blend24"): f"{SHOWCASE}/old_play24", rp.cache("states"): f"{SHOWCASE}/old_states",
           rp.receipt("c1.json"): f"{SHOWCASE}/old_c1.json",
           rp.receipt("summary.json"): f"{SHOWCASE}/old_summary.json",
           rp.receipt("field.json"): f"{SHOWCASE}/field_old.json",
           rp.input("solid_phi.npy"): f"{SHOWCASE}/solid_phi_old.npy",
           rp.input("colliders.obj"): f"{SHOWCASE}/old.obj",
           rp.input("collider_ab"): "collider_old_ab", rp.renders("audit"): f"{SHOWCASE}/old_png",
           rp.renders("playblast"): f"{SHOWCASE}/old_render", rp.renders("mist"): f"{SHOWCASE}/old_mist",
           rp.video("audit"): f"{SHOWCASE}/old.mp4", rp.scene("blend"): "play_old_blend.ma"}
    for got, rel in old.items():
        check(f"found in the older place: {rel}", same(got, os.path.join(out, rel)), got)
    check("a chunk not yet written goes to runs/", same(rp.dumps(2), os.path.join(out, "runs/old/dumps/c2")),
          rp.dumps(2))
    check("the full cache, absent in both, goes to runs/",
          same(rp.cache("full"), os.path.join(out, "runs/old/cache/full")))
    check("chunks() reads the older layout", rp.chunks() == [0, 1], rp.chunks())
    write(os.path.join(sc, "play_states", "cache.json"), text='{"frames": 4}')
    check("the finale's play_states is boulders' states cache",
          same(RunPaths("boulders", out).cache("states"), os.path.join(sc, "play_states")))
    rp2 = RunPaths("mixed", out)
    write(os.path.join(out, "loop3_mixed_c0", "frame.0000.npy"), 400)
    write(os.path.join(out, "runs", "mixed", "dumps", "c1", "frame.0000.npy"), 400)
    check("a run with chunks in both layouts lists both", rp2.chunks() == [0, 1], rp2.chunks())

    print("=== report")
    for n in ("loop3_resume_x_full/a.npy", "step_in/x.npy", "ver_base/a.npy", "play_gate/x.npy"):
        write(os.path.join(out, n), 100)
    write(os.path.join(out, "loose.log"), text="log")
    write(os.path.join(sc, "log_stage0.txt"), text="log")
    tags = runs.discover(out)
    check("discover finds the runs and nothing else", tags == ["boulders", "mixed", "old"], tags)
    d = runs.describe("old", out)
    check("dumps sized to the byte", d["dumps_bytes"] == b["dumps"], (d["dumps_bytes"], b["dumps"]))
    check("caches sized to the byte", {k: v["bytes"] for k, v in d["caches"].items()}
          == {"blend24": b["blend24"], "states": b["states"]}, d["caches"])
    check("the blend's frames and sources read", d["caches"]["blend24"]["frames"] == 2
          and d["caches"]["blend24"]["source_names"] == ["loop3_old_c0", "loop3_old_c1"])
    check("inputs sized to the byte", sum(v["bytes"] for v in d["inputs"].values()) == b["inputs"])
    check("renders: audit (frames + video), mist, playblast", set(d["renders"]) == {"audit", "mist", "playblast"}
          and sum(v["bytes"] for v in d["renders"].values()) == b["renders"], d["renders"])
    check("receipts: two chunks, summary, field",
          set(d["receipts"]) == {"c0.json", "c1.json", "summary.json", "field.json"})
    check("the layout is reported as older", d["layout"] == "older", d["layout"])
    check("the mixed run is reported as mixed", runs.describe("mixed", out)["layout"] == "mixed")
    s = runs.scan(out)
    names = {o["name"] for o in s["others"]}
    check("everything else: the scratch folders and loose files, not the runs' parts",
          names == {"loop3_resume_x_full", "step_in", "ver_base", "play_gate", "loose.log",
                    f"{SHOWCASE}/log_stage0.txt"}, names)
    labels = {o["name"]: o["label"] for o in s["others"]}
    check("scratch labelled", labels["loop3_resume_x_full"].startswith("a gate")
          and labels["play_gate"].startswith("a gate") and labels["loose.log"].startswith("a loose"), labels)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        runs.report(out)
    text = buf.getvalue()
    check("the report prints both tables", "old" in text and "boulders" in text and "Everything else" in text
          and "loose.log" in text)

    print("=== migrate")
    before = listing(out)
    total_before = sum(sz for _, sz in before)
    moves = runs.plan("old", out)
    rel = [(os.path.relpath(s, out).replace("\\", "/"), os.path.relpath(t, out).replace("\\", "/"))
           for s, t in moves]
    expected = [("loop3_old_c0", "runs/old/dumps/c0"), ("loop3_old_c1", "runs/old/dumps/c1"),
                (f"{SHOWCASE}/old_states", "runs/old/cache/states"),
                (f"{SHOWCASE}/old_play24", "runs/old/cache/blend24"),
                ("collider_old_ab", "runs/old/inputs/collider_ab"), ("collider_old_c", "runs/old/inputs/collider_c"),
                (f"{SHOWCASE}/solid_phi_old.npy", "runs/old/inputs/solid_phi.npy"),
                (f"{SHOWCASE}/old.obj", "runs/old/inputs/colliders.obj"),
                (f"{SHOWCASE}/old_c0.json", "runs/old/receipts/c0.json"),
                (f"{SHOWCASE}/old_c1.json", "runs/old/receipts/c1.json"),
                (f"{SHOWCASE}/old_summary.json", "runs/old/receipts/summary.json"),
                (f"{SHOWCASE}/field_old.json", "runs/old/receipts/field.json"),
                (f"{SHOWCASE}/old_mist", "runs/old/renders/mist"), (f"{SHOWCASE}/old_png", "runs/old/renders/audit"),
                (f"{SHOWCASE}/old_render", "runs/old/renders/playblast"),
                (f"{SHOWCASE}/old.mp4", "runs/old/renders/audit/old.mp4"),
                (f"{SHOWCASE}/old_mist_launch.txt", "runs/old/logs/mist_launch.txt"),
                ("play_old_blend.ma", "runs/old/scenes/play_blend.ma")]
    check("the plan is exactly the expected renames, in order", rel == expected,
          "\n    got      " + "\n    got      ".join(map(str, rel)))
    check("migratable picks the showcase runs only (old, and boulders for its play_states)",
          runs.migratable(out) == ["boulders", "old"], runs.migratable(out))
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        done = runs.apply(moves, "old", out)
    after = listing(out)
    receipt = os.path.join(out, "runs", "old", "migration.json")
    check("the bytes are all there", sum(sz for _, sz in after) == total_before + os.path.getsize(receipt))
    check("the file count is unchanged (plus the migration receipt)", len(after) == len(before) + 1)
    gone = [s for s, _ in moves if os.path.exists(s)]
    check("the old places are empty", not gone, gone)
    rp = RunPaths("old", out)
    check("every accessor now finds runs/old", all(os.path.normpath(p).startswith(os.path.normpath(rp.root)) for p in (
        rp.dumps(0), rp.cache("blend24"), rp.cache("states"), rp.receipt("c0.json"), rp.input("solid_phi.npy"),
        rp.renders("audit"), rp.video("audit"), rp.scene("blend"))))
    check("the same files, by size, under their new names (receipts and inputs are renamed)",
          sorted(sz for _, sz in before) == sorted(sz for p, sz in after if not p.endswith("migration.json")))
    with open(receipt) as f:
        rec = json.load(f)
    check("migration.json records every move", len(rec["moves"]) == len(expected) == len(done))
    check("a second plan is empty", runs.plan("old", out) == [])
    d2 = runs.describe("old", out)
    check("described the same afterwards, layout runs/", d2["layout"] == "runs/" and d2["dumps_bytes"] == b["dumps"]
          and {k: v["bytes"] for k, v in d2["caches"].items()} == {k: v["bytes"] for k, v in d["caches"].items()}
          and sum(v["bytes"] for v in d2["renders"].values()) == b["renders"], d2["layout"])
    write(os.path.join(out, "loop3_clash_c0", "a.npy"), 400)
    write(os.path.join(out, "runs", "clash", "dumps", "c0", "b.npy"), 400)
    write(os.path.join(sc, "clash_c0.json"), text="{}")
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            runs.apply(runs.plan("clash", out), "clash", out)
        refused = False
    except SystemExit:
        refused = True
    check("apply refuses an existing destination before moving anything",
          refused and os.path.exists(os.path.join(out, "loop3_clash_c0"))
          and os.path.exists(os.path.join(sc, "clash_c0.json")))

    ok = all(checks)
    print(f"RUN LAYOUT: {'PASS' if ok else 'FAIL'} ({sum(checks)} of {len(checks)} checks)")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
