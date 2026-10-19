"""Analyze saved reconstruction meshes and closest-point queries; no Maya calls."""
import hashlib
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "out/runs/breaker_surface_diagnosis_20261019_v3"
STAGES = ("raw", "no_mcf", "control")


def quantiles(values):
    return dict(zip(("p50", "p90", "p99", "max"), map(float, np.quantile(values, [.5, .9, .99, 1]))))


def mesh_metrics(mesh):
    points, tri = mesh["points"] * .25, mesh["triangles"]
    corners = points[tri]
    cross = np.cross(corners[:, 1] - corners[:, 0], corners[:, 2] - corners[:, 0])
    lengths = np.linalg.norm(cross, axis=1)
    valid = lengths > 1e-10
    ng = cross / np.maximum(lengths[:, None], 1e-30)
    # Current mesher produces triangles; verify ordering before reading corner normals.
    assert (mesh["face_counts"] == 3).all()
    assert np.array_equal(mesh["face_vertices"].reshape(-1, 3), tri)
    ns = mesh["normals"][mesh["normal_ids"]].reshape(-1, 3, 3)
    nl = np.linalg.norm(ns, axis=-1)
    ns = ns / np.maximum(nl[..., None], 1e-30)
    cosine = (ns * ng[:, None, :]).sum(axis=-1)
    edges = np.sort(np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]]), axis=1)
    unique, inverse, counts = np.unique(edges, axis=0, return_inverse=True, return_counts=True)
    parent = np.arange(len(points))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in unique:
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[rb] = ra
    roots = np.array([find(i) for i in range(len(points))])
    _, components = np.unique(roots, return_inverse=True)
    face_component = components[tri[:, 0]]
    component_faces = np.bincount(face_component)
    face_indices = np.tile(np.arange(len(tri)), 3)
    order = np.argsort(inverse, kind="stable")
    starts = np.r_[0, np.cumsum(counts[:-1])]
    manifold = counts == 2
    fa = face_indices[order[starts[manifold]]]
    fb = face_indices[order[starts[manifold] + 1]]
    good = valid[fa] & valid[fb]
    dihedral = np.degrees(np.arccos(np.clip((ng[fa[good]] * ng[fb[good]]).sum(axis=1), -1, 1)))
    return dict(vertices=len(points), triangles=len(tri), area_m2=float(lengths.sum() / 2),
                components=len(component_faces), components_under_100_faces=int((component_faces < 100).sum()),
                boundary_edges=int((counts == 1).sum()), nonmanifold_edges=int((counts > 2).sum()),
                degenerate_triangles=int((~valid).sum()),
                nonfinite_normals=int((~np.isfinite(mesh["normals"])).any(axis=1).sum()),
                normals_opposed_to_face_fraction=float((cosine[valid] < 0).mean()),
                shading_vs_geometric_angle_deg=quantiles(np.degrees(np.arccos(np.clip(cosine[valid], -1, 1))).ravel()),
                adjacent_face_angle_deg=quantiles(dihedral)), components


def section(mesh, z=68):
    corners = mesh["points"][mesh["triangles"]] * .25
    hit = (corners[:, :, 2].min(axis=1) < z) & (corners[:, :, 2].max(axis=1) >= z)
    t = corners[hit]
    segments = []
    for triangle in t:
        cross = []
        for a, b in ((0, 1), (1, 2), (2, 0)):
            p, q = triangle[a], triangle[b]
            if (p[2] < z) != (q[2] < z):
                cross.append((p + (q - p) * ((z - p[2]) / (q[2] - p[2])))[:2])
        if len(cross) == 2:
            segments.append(cross)
    return np.asarray(segments)


def main():
    worker = json.loads((OUT / "worker_result.json").read_text())
    receipt = json.loads((OUT / "receipt.json").read_text())
    assert worker["status"] == receipt["status"] == "passed"
    probes = np.load(OUT / "probes.npz")
    p, radii = probes["positions"], probes["radii"]
    meshes, closest, metrics = {}, {}, {}
    for stage in STAGES:
        meshes[stage] = np.load(OUT / (stage + "_mesh.npz"))
        closest[stage] = np.load(OUT / (stage + "_closest.npz"))
        metrics[stage], _ = mesh_metrics(meshes[stage])
    local, distances = {}, {}
    for stage in STAGES:
        delta = (p - closest[stage]["points"]) * .25
        distances[stage] = np.linalg.norm(delta, axis=1)
        local[stage] = (delta * closest[stage]["normals"]).sum(axis=1)
    # All probes are measured; restrict this comparison to centres near the raw
    # mesh to avoid treating correctly interior liquid particles as missing detail.
    near_raw = (radii < 1e-9) & (distances["raw"] <= .125)
    comparisons = {}
    for stage in STAGES:
        comparisons[stage] = {
            "near_raw_particle_center_distance_m": quantiles(distances[stage][near_raw]),
            "near_raw_local_outside_over_one_voxel": int((local[stage][near_raw] > .125).sum()),
            "near_raw_closest_point_shift_from_raw_m": quantiles(np.linalg.norm(
                (closest[stage]["points"] - closest["raw"]["points"])[near_raw] * .25, axis=1)),
            "flagged_particle_distance_m": quantiles(distances[stage][radii > 1e-9]),
        }
    # Assert the actual graph changes, including defaults that Maya may omit.
    graphs = {s: json.loads((OUT / (s + "_graph.json")).read_text())["compounds"][0] for s in STAGES}
    changes = {}
    for a, b in (("control", "no_mcf"), ("no_mcf", "raw")):
        av = {v["valueName"]: v["value"] for v in graphs[a]["values"]}
        bv = {v["valueName"]: v["value"] for v in graphs[b]["values"]}
        diff = {k: [av.get(k), bv.get(k)] for k in av.keys() | bv.keys() if av.get(k) != bv.get(k)}
        allowed = "particle_surface.mcf_iterations" if a == "control" else "particle_surface.smoothing"
        assert set(diff) == {allowed}, diff
        assert graphs[a]["connections"] == graphs[b]["connections"]
        changes[a + "_to_" + b] = diff
    result = dict(status="passed", run=receipt, execution=worker, geometry=metrics,
                  near_raw_unflagged_probes=int(near_raw.sum()), flagged_probes=int((radii > 1e-9).sum()),
                  correspondence=comparisons, adjacent_graph_changes=changes,
                  artifacts_sha256={f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                                    for f in OUT.iterdir() if f.suffix in (".npz", ".mb")},
                  limits=["Closest-normal side is a local test, not a global inside/winding test.",
                          "Normals were read through bifrostGeoToMaya; Arnold shading AOVs were not rendered.",
                          "Particle-centre distance is not a mass-conservation measurement.",
                          "No baked-scene beauty equivalence or improved ocean appearance is claimed."])
    (OUT / "analysis.json").write_text(json.dumps(result, indent=2) + "\n")
    # Pure geometric sections, not beauty renders; no lighting or exposure to hide features.
    allp = np.load(ROOT / "out/runs/wavesea4/cache/breaking_section/liquid.0001.npy") * .25
    allr = np.load(ROOT / "out/runs/wavesea4/cache/breaking_section/radius.0001.npy")
    slab = np.abs(allp[:, 2] - 68) <= .125
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    ax = axes.flat[0]
    ax.scatter(allp[slab & (allr < 1e-9), 0], allp[slab & (allr < 1e-9), 1], s=2, color="#117b8b")
    ax.scatter(allp[slab & (allr > 1e-9), 0], allp[slab & (allr > 1e-9), 1], s=12, color="#de6f27", label="Flagged: excluded from water mesh")
    ax.set_title("Particle centres | 0.25 m slab around z=68 m")
    ax.legend(fontsize=8)
    titles = ["Raw reconstruction | no MCF, no mesh smoothing", "Mesher smoothing 0.1 | no MCF", "Rendered control | MCF 10, mesher smoothing 0.1"]
    for ax, stage, title in zip(list(axes.flat)[1:], STAGES, titles):
        ax.add_collection(LineCollection(section(meshes[stage]), colors="#154b66", linewidths=.55))
        ax.set_title(title)
    for ax in axes.flat:
        ax.set_xlim(134, 152); ax.set_ylim(13, 22)
        ax.set_aspect("equal"); ax.set_xlabel("x (m)"); ax.set_ylabel("height (m)")
        ax.grid(alpha=.15)
    fig.suptitle("24.9 s curl: particle slab vs exact mesh sections at z=68 m\nMarker sizes are not particle radii; mesh sections contain no lighting or material", fontsize=13)
    fig.savefig(OUT / "curl_sections.png", dpi=150)
    plt.close(fig)
    print(json.dumps({"geometry": metrics, "correspondence": comparisons,
                      "near_raw_unflagged_probes": int(near_raw.sum()), "flagged_probes": int((radii > 1e-9).sum())}, indent=2))


if __name__ == "__main__":
    main()
