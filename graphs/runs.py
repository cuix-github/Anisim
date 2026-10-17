"""Where a run's files live under out/, and what is on disk there.

One folder per run, named by the run's tag, holding everything the run produced:

  out/runs/<tag>/
      dumps/c<N>/        chunk N as the loop runner writes it: the arrays it started
                         from (positions.npy ...), the per-substep dumps
                         (frame/vel/esc/scl/ph/diag.####.npy), the resume state (final*.npy)
      cache/states/      the Maya play caches (liquid/speed.####.npy and cache.json):
      cache/full/          one frame per chunk's end; every dumped frame; the 24 fps
      cache/blend24/       blend of the dumps (maya/play_cache.py)
      inputs/            the colliders: solid_phi.npy, colliders.obj, the setup
                         graph's scratch (collider_ab, collider_c)
      receipts/          c<N>.json, summary.json, field.json, cost.json
      renders/<set>/     one render set's frames, logs and mp4 (audit: the painter's
                         frames the run draws as it goes, and its preview video)
      scenes/            the Maya scenes the launcher builds (play_<mode>.ma)
      logs/

The dumps are the only copy of a run; every cache and render is derived from them
and can be rebuilt. Runs made before 2026-10-16 lie in the older layout
(out/loop3_<tag>_c<N>, out/collider_showcase/<tag>_*, out/play_<tag>_<mode>.ma);
RunPaths finds a run's files in either, the new layout first, so nothing had to
move for the code to change. As a script:

  .venv\\Scripts\\python.exe graphs\\runs.py report [--json]
      every run with its parts and sizes and what is derived from what, then
      everything else under out/ (the gates' and studies' scratch, loose files)
  .venv\\Scripts\\python.exe graphs\\runs.py migrate [<tag> ...] [--apply]
      the renames that would bring the older layout's runs into out/runs/: a plan
      only, until --apply, which renames (same volume, nothing copied, nothing
      deleted) and leaves runs/<tag>/migration.json saying what moved from where

ANISIM_OUT overrides the out/ root (the gate points it at a scratch folder).
"""
import json
import os
import re
import shutil
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
SHOWCASE = "collider_showcase"            # the older layout's folder for everything but dumps
CACHES = ("states", "full", "blend24", "crop")   # crop: the liquid around a tracking camera only
LEGACY_CACHE = {"states": "{tag}_states", "full": "{tag}_play", "blend24": "{tag}_play24"}
RUN_DIR = re.compile(r"^loop3_(.+)_c(\d+)$")
# the canonical profiling state (CLAUDE.md): a chunk folder by name, never a run to move
KEEP_IN_PLACE = {"D_finale": "the canonical profiling state: keep, ask before regenerating"}


def out_root(out=None):
    return os.path.normpath(out or os.environ.get("ANISIM_OUT") or os.path.join(ROOT, "out"))


def prepare(path):
    """Makes a file's folder; returns the path."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    return path


def _pick(new, old):
    """The new-layout path, unless only the older one exists."""
    if old is None or os.path.exists(new) or not os.path.exists(old):
        return new
    return old


class RunPaths:
    """Every path of one run, found in the new layout first and the older one second."""

    def __init__(self, tag, out=None):
        self.tag = tag
        self.out = out_root(out)
        self.root = os.path.join(self.out, "runs", tag)
        self.showcase = os.path.join(self.out, SHOWCASE)

    def __repr__(self):
        return f"RunPaths({self.tag!r}, {self.root!r})"

    def dumps(self, chunk):
        return _pick(os.path.join(self.root, "dumps", f"c{chunk}"),
                     os.path.join(self.out, f"loop3_{self.tag}_c{chunk}"))

    def chunks(self):
        """The chunk indices on disk, in order, from both layouts."""
        found = set()
        d = os.path.join(self.root, "dumps")
        if os.path.isdir(d):
            found |= {int(n[1:]) for n in os.listdir(d) if re.fullmatch(r"c\d+", n)}
        if os.path.isdir(self.out):
            for n in os.listdir(self.out):
                m = RUN_DIR.match(n)
                if m and m.group(1) == self.tag:
                    found.add(int(m.group(2)))
        return sorted(found)

    def cache(self, kind):
        """states, full or blend<fps>; the older layout knew blend24 only."""
        assert kind in CACHES or re.fullmatch(r"blend\d+", kind), kind
        old = LEGACY_CACHE.get(kind, "").format(tag=self.tag) or None
        if kind == "states" and self.tag == "boulders":
            old = "play_states"             # the finale run's cache kept its first name
        return _pick(os.path.join(self.root, "cache", kind),
                     os.path.join(self.showcase, old) if old else None)

    def receipt(self, name):
        """c<N>.json, summary.json, field.json or cost.json."""
        legacy = {"summary.json": f"{self.tag}_summary.json", "field.json": f"field_{self.tag}.json",
                  "cost.json": "cost.json"}
        old = legacy.get(name) or (f"{self.tag}_{name}" if re.fullmatch(r"c\d+\.json", name) else None)
        return _pick(os.path.join(self.root, "receipts", name),
                     os.path.join(self.showcase, old) if old else None)

    def input(self, name):
        """solid_phi.npy, colliders.obj, or collider_<x> (the setup graph's scratch)."""
        legacy = {"solid_phi.npy": os.path.join(self.showcase, f"solid_phi_{self.tag}.npy"),
                  "colliders.obj": os.path.join(self.showcase, f"{self.tag}.obj")}
        old = legacy.get(name)
        if old is None and name.startswith("collider_"):
            old = os.path.join(self.out, f"collider_{self.tag}_{name[len('collider_'):]}")
        return _pick(os.path.join(self.root, "inputs", name), old)

    def renders(self, set_name):
        """A render set's folder; audit is the painter's frames, playblast the viewport's."""
        legacy = {"audit": f"{self.tag}_png", "playblast": f"{self.tag}_render"}
        old = legacy.get(set_name, f"{self.tag}_{set_name}")
        return _pick(os.path.join(self.root, "renders", set_name), os.path.join(self.showcase, old))

    def video(self, set_name):
        """The mp4 of a render set (audit's is the run's preview video)."""
        old = os.path.join(self.showcase, f"{self.tag}.mp4") if set_name == "audit" else None
        return _pick(os.path.join(self.root, "renders", set_name, f"{self.tag}.mp4"), old)

    def scene(self, mode):
        return _pick(os.path.join(self.root, "scenes", f"play_{mode}.ma"),
                     os.path.join(self.out, f"play_{self.tag}_{mode}.ma"))

    @property
    def logs(self):
        return os.path.join(self.root, "logs")


# ------------------------------------------------------------------ what is on disk


def du(path):
    """(bytes, files) under path; a file counts itself; a missing path is (0, 0)."""
    if os.path.isfile(path):
        return os.path.getsize(path), 1
    total = files = 0
    stack = [path] if os.path.isdir(path) else []
    while stack:
        try:
            with os.scandir(stack.pop()) as it:
                for e in it:
                    if e.is_dir(follow_symlinks=False):
                        stack.append(e.path)
                    elif e.is_file(follow_symlinks=False):
                        total += e.stat(follow_symlinks=False).st_size
                        files += 1
        except OSError:
            pass
    return total, files


def human(n):
    """Binary units, the ones the file browser shows (1 GB = 1024^3 bytes)."""
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024.0


def _index(cache_dir):
    try:
        with open(os.path.join(cache_dir, "cache.json")) as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def _sources(idx):
    """The dump folders a cache's index says it came from, as names."""
    src = idx.get("blended_from", {}).get("dirs") if isinstance(idx.get("blended_from"), dict) else None
    src = src or idx.get("source")
    return src if isinstance(src, list) else ([src] if src else [])


def discover(out=None):
    """Every run tag with anything on disk, in either layout."""
    out = out_root(out)
    tags = set()
    runs_dir = os.path.join(out, "runs")
    if os.path.isdir(runs_dir):
        tags |= {n for n in os.listdir(runs_dir) if os.path.isdir(os.path.join(runs_dir, n))}
    if os.path.isdir(out):
        for n in os.listdir(out):
            m = RUN_DIR.match(n)
            if m:
                tags.add(m.group(1))
    sc = os.path.join(out, SHOWCASE)
    if os.path.isdir(sc):
        for n in os.listdir(sc):
            if n == "play_states":          # the finale run's states cache, by its first name
                continue
            for pat in (r"^(.+)_(states|play|play24|png|render)$", r"^(.+)_c\d+\.json$",
                        r"^(.+)_summary\.json$", r"^solid_phi_(.+)\.npy$", r"^field_(.+)\.json$"):
                m = re.match(pat, n)
                if m:
                    tags.add(m.group(1))
                    break
        if os.path.exists(os.path.join(sc, "play_states")):
            tags.add("boulders")
    return sorted(tags)


def describe(tag, out=None):
    """One run's parts: where each is, how big, and what it is derived from."""
    out = out_root(out)
    rp = RunPaths(tag, out)
    claimed = []                                   # every path this run accounts for

    def take(p):
        claimed.append(os.path.normpath(p))
        return p
    chunks = rp.chunks()
    dumps = [take(rp.dumps(c)) for c in chunks]
    d_bytes = sum(du(d)[0] for d in dumps)
    caches = {}
    for kind in CACHES:
        p = rp.cache(kind)
        if os.path.isdir(p):
            idx = _index(take(p)) or {}
            caches[kind] = {"path": p, "bytes": du(p)[0], "frames": idx.get("frames"),
                            "from": "the dumps" if kind != "states" else "the chunks' final states",
                            "source_names": _sources(idx), "dumps_present": bool(chunks)}
    receipts = {}
    for d, pat in ((os.path.join(rp.root, "receipts"), r"c(\d+)\.json"),
                   (rp.showcase, re.escape(tag) + r"_c(\d+)\.json")):
        for n in (os.listdir(d) if os.path.isdir(d) else []):
            m = re.fullmatch(pat, n)
            if m and f"c{m.group(1)}.json" not in receipts:
                receipts[f"c{m.group(1)}.json"] = take(rp.receipt(f"c{m.group(1)}.json"))
    for name in ("summary.json", "field.json"):
        p = rp.receipt(name)
        if os.path.isfile(p):
            receipts[name] = take(p)
    inputs = {}
    for name in ("solid_phi.npy", "colliders.obj", "collider_ab", "collider_c"):
        p = rp.input(name)
        if os.path.exists(p):
            inputs[name] = {"path": take(p), "bytes": du(p)[0]}
    renders = {}
    new_r = os.path.join(rp.root, "renders")
    if os.path.isdir(new_r):
        for n in os.listdir(new_r):
            renders[n] = {"path": take(os.path.join(new_r, n)), "bytes": du(os.path.join(new_r, n))[0]}
    if os.path.isdir(rp.showcase):
        names = sorted(os.listdir(rp.showcase), key=lambda n: not os.path.isdir(os.path.join(rp.showcase, n)))
        for n in names:                             # folders first, then the files beside them
            p = os.path.join(rp.showcase, n)
            m = re.fullmatch(re.escape(tag) + r"_(.+)", n)
            if m and os.path.isdir(p) and m.group(1) not in ("states", "play", "play24"):
                name = {"png": "audit", "render": "playblast"}.get(m.group(1), m.group(1))
                if name not in renders:
                    renders[name] = {"path": take(p), "bytes": du(p)[0]}
            elif n == f"{tag}.mp4":
                renders.setdefault("audit", {"path": None, "bytes": 0})
                renders["audit"]["video"] = take(p)
                renders["audit"]["bytes"] += du(p)[0]
            elif m and m.group(1).endswith("_launch.txt"):
                take(p)
    scenes = []
    new_s = os.path.join(rp.root, "scenes")
    if os.path.isdir(new_s):
        scenes += [take(os.path.join(new_s, n)) for n in os.listdir(new_s)]
    if os.path.isdir(out):
        scenes += [take(os.path.join(out, n)) for n in os.listdir(out)
                   if re.fullmatch(re.escape(f"play_{tag}_") + r".+\.ma", n)]
    legacy = [p for p in claimed if not p.startswith(os.path.normpath(rp.root))]
    total = (d_bytes + sum(c["bytes"] for c in caches.values()) + sum(r["bytes"] for r in renders.values())
             + sum(i["bytes"] for i in inputs.values()))
    layout = "mixed" if legacy and len(legacy) < len(claimed) else ("older" if legacy else "runs/")
    return {"tag": tag, "root": rp.root, "chunks": chunks, "dumps_bytes": d_bytes, "caches": caches,
            "receipts": receipts, "inputs": inputs, "renders": renders, "scenes": scenes,
            "bytes": total, "layout": layout, "claimed": claimed, "note": KEEP_IN_PLACE.get(tag, "")}


OTHER_KINDS = (
    (r"^loop3_", "a gate's or study's loop run"),
    (r"^loop_", "a 2D loop gate's run"),
    (r"^(step|step3|ad|esc|sol|stn)_(in|out|out_mg|hydro|hydro_out)$", "a node gate's inputs or outputs"),
    (r"^ver_", "verify_engineering scratch"),
    (r"^collider_", "collider setup scratch"),
    (r"_gate$", "a gate's working folder"),
    (r"^(threading|adaptive_long|adaptive_step|adaptive_3d|adaptive_cost|bore_suite|level_set_node"
     r"|surface_test|seq|footage)$", "a gate's or study's working folder"),
    (r"\.(log|txt|py)$", "a loose log or script"),
    (r"\.npy$", "a stray array"),
    (r"\.ma$", "a Maya scene"),
)


def _other(path, name, label):
    idx = _index(path) if os.path.isdir(path) else None
    if idx is not None:
        label = "a play cache from " + (", ".join(_sources(idx)) or "?")
    if label is None:
        label = next((lab for pat, lab in OTHER_KINDS if re.search(pat, os.path.basename(name))), "")
    b, files = du(path)
    return {"name": name, "bytes": b, "files": files, "label": label}


def scan(out=None):
    """Every run described, then everything else under out/ with a label."""
    out = out_root(out)
    runs = [describe(t, out) for t in discover(out)]
    claimed = {p for r in runs for p in r["claimed"]}
    others = []
    if os.path.isdir(out):
        for n in sorted(os.listdir(out)):
            p = os.path.normpath(os.path.join(out, n))
            if n == "runs" or p in claimed:
                continue
            if n == SHOWCASE:
                for m in sorted(os.listdir(p)):
                    q = os.path.normpath(os.path.join(p, m))
                    if q not in claimed:
                        o = _other(q, f"{SHOWCASE}/{m}", None)
                        o["label"] = o["label"] or "in the older showcase folder, not tied to a run by name"
                        others.append(o)
                continue
            others.append(_other(p, n, None))
    others.sort(key=lambda o: -o["bytes"])
    return {"out": out, "runs": sorted(runs, key=lambda r: -r["bytes"]), "others": others}


def report(out=None, as_json=False):
    s = scan(out)
    if as_json:
        print(json.dumps(s, indent=1, default=str))
        return s
    print(f"Runs under {s['out']}  (layout: runs/ = out/runs/<tag>; older = the pre-2026-10-16 folders)")
    print("The dumps are the only copy of a run; caches and renders are derived from them and can be rebuilt.\n")
    hdr = (f"{'run':<14}{'chunks':>7}{'dumps':>10}  {'caches (derived)':<66}{'renders':>9}  "
           f"{'receipts':<20}{'layout':<8}note")
    print(hdr)
    print("-" * len(hdr))
    for r in s["runs"]:
        caches = ", ".join(f"{k} {human(v['bytes'])}" + (f" ({v['frames']} fr)" if v.get("frames") else "")
                           for k, v in r["caches"].items()) or "-"
        n_chunks = sum(1 for k in r["receipts"] if k.startswith("c"))
        rec = (f"{n_chunks} chunks" + (" + summary" if "summary.json" in r["receipts"] else "")) if n_chunks \
            else ("field only" if "field.json" in r["receipts"] else "-")
        rend = human(sum(v["bytes"] for v in r["renders"].values())) if r["renders"] else "-"
        print(f"{r['tag']:<14}{len(r['chunks']):>7}{human(r['dumps_bytes']):>10}  {caches:<66}{rend:>9}  "
              f"{rec:<20}{r['layout']:<8}{r['note']}")
    total_runs = sum(r["bytes"] for r in s["runs"])
    print(f"\n{'runs total':<33}{human(total_runs):>10}")
    print("\nEverything else under out/ (largest first)")
    hdr2 = f"{'folder or file':<44}{'size':>10}{'files':>8}  what"
    print(hdr2)
    print("-" * len(hdr2))
    for o in s["others"]:
        print(f"{o['name']:<44}{human(o['bytes']):>10}{o['files']:>8}  {o['label']}")
    total_other = sum(o["bytes"] for o in s["others"])
    print(f"\n{'everything else total':<33}{human(total_other):>10}")
    print(f"{'out/ total':<33}{human(total_runs + total_other):>10}")
    try:
        usage = shutil.disk_usage(s["out"])
        print(f"{'disk free':<33}{human(usage.free):>10} of {human(usage.total)}")
    except OSError:
        pass
    return s


# ------------------------------------------------------------------ migration


def plan(tag, out=None):
    """The renames that bring one older-layout run into out/runs/<tag>/, folders
    before the files that go inside them; empty when nothing is in the older layout."""
    out = out_root(out)
    rp = RunPaths(tag, out)
    sc = rp.showcase
    moves = []

    def mv(src, dst):
        if os.path.exists(src) and os.path.normpath(src) != os.path.normpath(dst):
            moves.append((os.path.normpath(src), os.path.normpath(dst)))
    for c in rp.chunks():
        mv(os.path.join(out, f"loop3_{tag}_c{c}"), os.path.join(rp.root, "dumps", f"c{c}"))
    for kind in CACHES:
        old = "play_states" if (kind == "states" and tag == "boulders") else LEGACY_CACHE[kind].format(tag=tag)
        mv(os.path.join(sc, old), os.path.join(rp.root, "cache", kind))
    for suffix in ("ab", "c"):
        mv(os.path.join(out, f"collider_{tag}_{suffix}"), os.path.join(rp.root, "inputs", f"collider_{suffix}"))
    mv(os.path.join(sc, f"solid_phi_{tag}.npy"), os.path.join(rp.root, "inputs", "solid_phi.npy"))
    mv(os.path.join(sc, f"{tag}.obj"), os.path.join(rp.root, "inputs", "colliders.obj"))
    names = sorted(os.listdir(sc)) if os.path.isdir(sc) else []
    for n in names:
        if re.fullmatch(re.escape(tag) + r"_c\d+\.json", n):
            mv(os.path.join(sc, n), os.path.join(rp.root, "receipts", n[len(tag) + 1:]))
    mv(os.path.join(sc, f"{tag}_summary.json"), os.path.join(rp.root, "receipts", "summary.json"))
    mv(os.path.join(sc, f"field_{tag}.json"), os.path.join(rp.root, "receipts", "field.json"))
    for n in names:                                 # render folders first, then the files inside
        m = re.fullmatch(re.escape(tag) + r"_(.+)", n)
        if m and os.path.isdir(os.path.join(sc, n)) and m.group(1) not in ("states", "play", "play24"):
            set_name = {"png": "audit", "render": "playblast"}.get(m.group(1), m.group(1))
            mv(os.path.join(sc, n), os.path.join(rp.root, "renders", set_name))
    mv(os.path.join(sc, f"{tag}.mp4"), os.path.join(rp.root, "renders", "audit", f"{tag}.mp4"))
    for n in names:
        m = re.fullmatch(re.escape(tag) + r"_(.+)_launch\.txt", n)
        if m:
            mv(os.path.join(sc, n), os.path.join(rp.root, "logs", f"{m.group(1)}_launch.txt"))
    for n in (sorted(os.listdir(out)) if os.path.isdir(out) else []):
        m = re.fullmatch(re.escape(f"play_{tag}_") + r"(.+)\.ma", n)
        if m:
            mv(os.path.join(out, n), os.path.join(rp.root, "scenes", f"play_{m.group(1)}.ma"))
    return moves


def migratable(out=None):
    """The older-layout runs worth moving: those with receipts, a cache or a solid
    field (a showcase run), never the canonical profiling state."""
    tags = []
    for r in scan(out)["runs"]:
        if r["tag"] in KEEP_IN_PLACE or r["layout"] == "runs/":
            continue
        if r["receipts"] or r["caches"] or r["inputs"]:
            tags.append(r["tag"])
    return sorted(tags)


def apply(moves, tag, out=None):
    """Renames, after checking that no destination exists; writes migration.json."""
    out = out_root(out)
    for src, dst in moves:
        if os.path.exists(dst):
            raise SystemExit(f"migrate: refusing, destination exists: {dst}")
        if not os.path.exists(src):
            raise SystemExit(f"migrate: refusing, source vanished: {src}")
    done = []
    for src, dst in moves:
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        os.rename(src, dst)
        done.append({"from": os.path.relpath(src, out), "to": os.path.relpath(dst, out)})
        print(f"  {done[-1]['from']}  ->  {done[-1]['to']}")
    rp = RunPaths(tag, out)
    receipt = os.path.join(rp.root, "migration.json")
    previous = []
    if os.path.exists(receipt):
        with open(receipt) as f:
            previous = json.load(f).get("moves", [])
    os.makedirs(rp.root, exist_ok=True)
    with open(receipt, "w") as f:
        json.dump({"tag": tag, "when": time.strftime("%Y-%m-%d %H:%M:%S"), "moves": previous + done},
                  f, indent=1)
    return done


def main(argv):
    as_json, do_apply = "--json" in argv, "--apply" in argv
    words = [a for a in argv if not a.startswith("--")]
    cmd = words[0] if words else "report"
    if cmd == "report":
        report(as_json=as_json)
    elif cmd == "migrate":
        tags = words[1:] or migratable()
        if not tags:
            print("nothing in the older layout to move")
            return
        total = 0
        for tag in tags:
            moves = plan(tag)
            size = sum(du(s)[0] for s, _ in moves)
            total += size
            print(f"{tag}: {len(moves)} renames, {human(size)}" + ("" if do_apply else " (a plan; --apply does it)"))
            if do_apply:
                apply(moves, tag)
            else:
                for s, d in moves:
                    print(f"  {os.path.relpath(s, out_root())}  ->  {os.path.relpath(d, out_root())}")
        print(f"{'moved' if do_apply else 'would move'} {human(total)} in {len(tags)} run(s); "
              "nothing copied, nothing deleted")
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
