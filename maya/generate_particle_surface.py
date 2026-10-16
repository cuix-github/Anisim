"""Generates the particle surface and mist graphs, stock nodes only: one frame of a
play cache (maya/play_cache.py: liquid.####.npy) becomes a water surface mesh, the
way the PF-FLIP and ST-FLIP papers render theirs (a union of spheres splatted onto a
grid finer than the simulation's, smoothed, contoured at the zero level), and the
particles that surface leaves out (spray and mist: at one cell per metre a lone
particle is a cloud of droplets, not a ball) become a fog-density volume.

particle_surface.json  the compound User::PFFlip::particle_surface: frame picks the
    cache frame; point_size is the sphere radius per particle in cells (0.5 = the
    papers' radius of half a cell at 8 particles a cell), detail_size the voxel size
    of the surface grid in cells (0.5 = twice the simulation resolution),
    smooth_iterations the Gaussian passes on the level set, smoothing the mesher's
    own smoothing, offset pushes the surface out (+) or in (-). Out comes the mesh,
    drawn by the viewport and rendered by Arnold (render_surface_preview.py). A
    sphere survives voxelization only when its radius spans two voxels, so isolated
    particles drop out of this mesh at the defaults: they are the mist's.

particle_mist.json  the compound User::PFFlip::particle_mist: the surface's own
    smoothed level set (surface_radius, surface_detail, smooth_iterations as the
    surface's), sampled at every particle. Negative is inside the water the mesh
    shows; a lone droplet or a one-particle-thick sheet is eroded away by the
    smoothing and samples positive (or the default, far outside, where no level set
    exists). Those above mist_threshold are mist: by construction the water the
    mesh dropped. Each is drawn as the cloud of droplets its water is (PF-FLIP
    Section 7 rasterises spray into a density grid finer than the pressure grid):
    its droplet radius r from radius_pattern (radius.####.npy of the play cache;
    default_radius where it has none), its cloud a sphere of cloud_radius cells, and
    its extinction the cloud's: a parcel of V_p = dx^3 / 8 of water as drops of
    radius r has a total cross-section of 1.5 V_p / r (geometric optics, Q = 2), so
    the extinction per cell at the cloud's centre is k_cloud * dx / (r R^3) with
    k_cloud = 1.5 / (8 c_k 4/3 pi), c_k the splat kernel's mean over its sphere
    (measured by the gate; 0.25 for the linear kernel), dx = metres_per_cell and R =
    cloud_radius. mist_density scales it: the fraction of the parcel's water drawn as
    droplets of its radius, 1 all of it. The clouds are splatted additively on a
    grid of spray_detail cells (splat_points_into_volume divides by the kernel sum
    plus spray_epsilon, so with spray_epsilon far above any kernel sum the stored
    channel is the sum over spray_epsilon: the shader's density multiplier is
    spray_epsilon). Out come the volume (channel voxel_fog_density), the mist positions,
    every particle's depth (that sample) and the mist particles' weights.

particle_surface_test.json / particle_mist_test.json  the bifcmd wrappers the gate
    drives: the same compounds at a given frame, the mesh's vertex positions, or the
    mist positions and every particle's depth, written back as .npy.

Regenerate after editing:  python generate_particle_surface.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
from generate_sim_graph import P  # noqa: E402
from generate_particle_play import HEADER, f3  # noqa: E402

KNOBS = (("point_size", "float", "0.5f"), ("detail_size", "float", "0.5f"),
         ("offset", "float", "0f"), ("smooth_iterations", "int", "1L"),
         ("smoothing", "float", "0.1f"))
MIST_KNOBS = (("surface_radius", "float", "0.5f"), ("surface_detail", "float", "0.5f"),
              ("smooth_iterations", "int", "1L"),
              ("mist_threshold", "float", "0f"), ("spray_detail", "float", "0.5f"),
              ("cloud_radius", "float", "1f"), ("metres_per_cell", "float", "1f"),
              ("default_radius", "float", "0.0001f"), ("mist_density", "float", "1f"),
              ("spray_epsilon", "float", "10000f"))
# the cloud constant k_cloud = 1.5 / (8 c_k 4/3 pi) with c_k = 0.25, the linear
# kernel's mean over its sphere (the gate measures it against this)
K_CLOUD = 1.5 / (8.0 * 0.25 * 4.0 / 3.0 * 3.141592653589793)
F3 = "array<Math::float3>"


def reader(nodes, values, conns):
    """The cache frame's liquid positions: .frame and .liquid_pattern -> rpos.data."""
    nodes += [{"nodeName": "t_pos", "valueType": F3},
              {"nodeName": "f2i", "nodeType": "Core::Type_Conversion::to_int"},
              {"nodeName": "rpos", "nodeType": "File::NumPy::read_NumPy"}]
    values += [{"valueName": "t_pos.value", "valueType": F3, "value": [f3(0, 0, 0)]}]
    conns += [{"source": ".frame", "target": "f2i.from"},
              {"source": "t_pos.output", "target": "rpos.type"},
              {"source": ".liquid_pattern", "target": "rpos.file_path"},
              {"source": "f2i.int", "target": "rpos.frame"}]


def spheres(nodes, values, conns, tag, positions, size):
    """Points at `positions` with point_size from `size`: <tag>_size.out_geometry."""
    nodes += [{"nodeName": tag + "_pts", "nodeType": "Geometry::Points::construct_points"},
              {"nodeName": tag + "_size", "nodeType": "Geometry::Properties::set_geo_property"}]
    values += [{"valueName": tag + "_size.property", "valueType": "string", "value": "point_size"}]
    conns += [{"source": positions, "target": tag + "_pts.point_position"},
              {"source": tag + "_pts.points", "target": tag + "_size.geometry"},
              {"source": size, "target": tag + "_size.data"},
              {"source": size, "target": tag + "_size.default"}]
    return tag + "_size.out_geometry"


def splat(nodes, values, conns, tag, points, detail, level_set):
    """points_to_volume of `points`: a level set or a fog density; <tag>.volume."""
    nodes += [{"nodeName": tag, "nodeType": "Geometry::Converters::points_to_volume"}]
    values += [{"valueName": tag + ".store_level_set", "valueType": "bool",
                "value": "true" if level_set else "false"},
               {"valueName": tag + ".store_fog_density", "valueType": "bool",
                "value": "false" if level_set else "true"}]
    conns += [{"source": points, "target": tag + ".points"},
              {"source": detail, "target": tag + ".detail_size"}]
    return tag + ".volume"


def surface(name):
    ports = [P("liquid_pattern", "input", "string", ""), P("frame", "input", "float", "1f")]
    ports += [P(n, "input", t, d) for n, t, d in KNOBS]
    ports.append(P("mesh", "output", "Object"))
    nodes, values, conns = [], [], []
    reader(nodes, values, conns)
    pts = spheres(nodes, values, conns, "s", "rpos.data", ".point_size")
    vol = splat(nodes, values, conns, "p2v", pts, ".detail_size", level_set=True)
    nodes += [{"nodeName": "sm", "nodeType": "Geometry::Volume::smooth_level_set_property"},
              {"nodeName": "v2m", "nodeType": "Geometry::Converters::volume_to_mesh",
               "fanInPortNames": {"volumes": ["vol"]}},
              {"nodeName": "zero", "valueType": "long"},
              {"nodeName": "first", "nodeType": "Core::Array::get_from_array"}]
    values += [{"valueName": "zero.value", "valueType": "long", "value": "0L"}]
    conns += [{"source": ".offset", "target": "p2v.offset"},
              {"source": vol, "target": "sm.volume"},
              {"source": ".smooth_iterations", "target": "sm.iterations"},
              {"source": "sm.out_volume", "target": "v2m.volumes.vol"},
              {"source": ".smoothing", "target": "v2m.smoothing"},
              {"source": "v2m.meshes", "target": "first.array"},
              {"source": "zero.output", "target": "first.index"},
              {"source": "first.value", "target": ".mesh"}]
    return {"name": name, "ports": ports, "compoundNodes": nodes, "connections": conns,
            "values": values}


def mist(name):
    ports = [P("liquid_pattern", "input", "string", ""), P("radius_pattern", "input", "string", ""),
             P("frame", "input", "float", "1f")]
    ports += [P(n, "input", t, d) for n, t, d in MIST_KNOBS]
    ports += [P("mist", "output", "Object"), P("mist_positions", "output", F3),
              P("depth", "output", "array<float>"), P("spray_weights", "output", "array<float>")]
    nodes, values, conns = [], [], []
    reader(nodes, values, conns)
    # the droplet radii, read like the positions
    nodes += [{"nodeName": "t_flt", "valueType": "array<float>"},
              {"nodeName": "rrad", "nodeType": "File::NumPy::read_NumPy"}]
    values += [{"valueName": "t_flt.value", "valueType": "array<float>", "value": ["0f"]}]
    conns += [{"source": "t_flt.output", "target": "rrad.type"},
              {"source": ".radius_pattern", "target": "rrad.file_path"},
              {"source": "f2i.int", "target": "rrad.frame"}]
    # every particle's depth: the surface's smoothed level set sampled at the
    # particle (far from any surface a sample is the default, well outside)
    pts = spheres(nodes, values, conns, "c", "rpos.data", ".surface_radius")
    vol = splat(nodes, values, conns, "p2v_c", pts, ".surface_detail", level_set=True)
    nodes += [{"nodeName": "sm", "nodeType": "Geometry::Volume::smooth_level_set_property"},
              {"nodeName": "far", "valueType": "float"},
              {"nodeName": "smp", "nodeType": "Geometry::Query::sample_volume"}]
    values += [{"valueName": "far.value", "valueType": "float", "value": "10f"},
               {"valueName": "smp.property", "valueType": "string", "value": "voxel_signed_distance"}]
    conns += [{"source": vol, "target": "sm.volume"},
              {"source": ".smooth_iterations", "target": "sm.iterations"},
              {"source": "sm.out_volume", "target": "smp.volume"},
              {"source": "rpos.data", "target": "smp.positions"},
              {"source": "far.output", "target": "smp.default"}]
    depth = "smp.sampled_data"
    nodes += [{"nodeName": "thin", "nodeType": "Core::Logic::greater"},
              {"nodeName": "pick", "nodeType": "Core::Array::filter_array"},
              {"nodeName": "pickr", "nodeType": "Core::Array::filter_array"}]
    conns += [{"source": depth, "target": ".depth"},
              # mist = depth > threshold; its positions and radii
              {"source": depth, "target": "thin.first"},
              {"source": ".mist_threshold", "target": "thin.second"},
              {"source": "rpos.data", "target": "pick.array"},
              {"source": "thin.output", "target": "pick.boolean_array"},
              {"source": "pick.filtered_array", "target": ".mist_positions"},
              {"source": "rrad.data", "target": "pickr.array"},
              {"source": "thin.output", "target": "pickr.boolean_array"}]
    # r_eff = r + default_radius where r is 0 (a mist particle without a droplet radius)
    nodes += [{"nodeName": "tiny", "valueType": "float"},
              {"nodeName": "none", "nodeType": "Core::Logic::less"},
              {"nodeName": "none_f", "nodeType": "Core::Type_Conversion::to_float"},
              {"nodeName": "fill", "nodeType": "Core::Math::multiply", "multiInPortNames": ["a", "b"]},
              {"nodeName": "r_eff", "nodeType": "Core::Math::add", "multiInPortNames": ["a", "b"]}]
    values += [{"valueName": "tiny.value", "valueType": "float", "value": "0.000000001f"}]
    conns += [{"source": "pickr.filtered_array", "target": "none.first"},
              {"source": "tiny.output", "target": "none.second"},
              {"source": "none.output", "target": "none_f.from"},
              {"source": "none_f.float", "target": "fill.first.a"},
              {"source": ".default_radius", "target": "fill.first.b"},
              {"source": "pickr.filtered_array", "target": "r_eff.first.a"},
              {"source": "fill.output", "target": "r_eff.first.b"}]
    # the weight: mist_density * metres_per_cell * k_cloud / (r_eff * cloud_radius^3)
    nodes += [{"nodeName": "three", "valueType": "float"},
              {"nodeName": "r3", "nodeType": "Core::Math::power"},
              {"nodeName": "k_cloud", "valueType": "float"},
              {"nodeName": "num", "nodeType": "Core::Math::multiply", "multiInPortNames": ["d", "m", "k"]},
              {"nodeName": "den", "nodeType": "Core::Math::multiply", "multiInPortNames": ["r", "r3"]},
              {"nodeName": "w", "nodeType": "Core::Math::divide", "multiInPortNames": ["n", "d"]}]
    values += [{"valueName": "three.value", "valueType": "float", "value": "3f"},
               {"valueName": "k_cloud.value", "valueType": "float", "value": f"{K_CLOUD:.7f}f"}]
    conns += [{"source": ".cloud_radius", "target": "r3.base"},
              {"source": "three.output", "target": "r3.exponent"},
              {"source": ".mist_density", "target": "num.first.d"},
              {"source": ".metres_per_cell", "target": "num.first.m"},
              {"source": "k_cloud.output", "target": "num.first.k"},
              {"source": "r_eff.output", "target": "den.first.r"},
              {"source": "r3.power", "target": "den.first.r3"},
              {"source": "num.output", "target": "w.first.n"},
              {"source": "den.output", "target": "w.first.d"},
              {"source": "w.output", "target": ".spray_weights"}]
    # the clouds: points of cloud_radius carrying the weight, a volume of spray_detail
    # voxels over them, the weight splatted additively (the kernel sum over
    # spray_epsilon) into the channel voxel_fog_density
    mpts = spheres(nodes, values, conns, "m", "pick.filtered_array", ".cloud_radius")
    nodes += [{"nodeName": "setw", "nodeType": "Geometry::Properties::set_geo_property"},
              {"nodeName": "zero_f", "valueType": "float"},
              {"nodeName": "rvox", "nodeType": "Core::Math::divide", "multiInPortNames": ["r", "h"]},
              {"nodeName": "sp_m", "nodeType": "Geometry::Volume::splat_points_into_volume"}]
    values += [{"valueName": "setw.property", "valueType": "string", "value": "voxel_fog_density"},
               {"valueName": "zero_f.value", "valueType": "float", "value": "0f"},
               {"valueName": "sp_m.create_properties", "valueType": "bool", "value": "true"},
               {"valueName": "sp_m.properties", "valueType": "string", "value": "voxel_fog_density"},
               {"valueName": "sp_m.kernel", "valueType": "Geometry::Volume::SplatKernelType",
                "value": "kLinearKernel"},
               {"valueName": "sp_m.add_to_weights", "valueType": "float", "value": "0f"},
               {"valueName": "sp_m.smoothing", "valueType": "float", "value": "0f"},
               {"valueName": "sp_m.coarsest_depth", "valueType": "int", "value": "0"}]
    conns += [{"source": mpts, "target": "setw.geometry"},
              {"source": "w.output", "target": "setw.data"},
              {"source": "zero_f.output", "target": "setw.default"}]
    # the volume the clouds are splatted into is a fog volume with no other property:
    # a level set beside the density is what Arnold reads as the volume's primary
    # channel, and nothing renders (found 2026-10-16); its own ramp is zeroed
    mvol = splat(nodes, values, conns, "p2v_m", "setw.out_geometry", ".spray_detail", level_set=False)
    conns += [{"source": "zero_f.output", "target": "p2v_m.fog_density"},
              {"source": ".cloud_radius", "target": "rvox.first.r"},
              {"source": ".spray_detail", "target": "rvox.first.h"},
              {"source": "rvox.output", "target": "sp_m.radius"},
              {"source": mvol, "target": "sp_m.volume"},
              {"source": "setw.out_geometry", "target": "sp_m.points"},
              {"source": ".spray_epsilon", "target": "sp_m.add_to_denominator"},
              {"source": "sp_m.out_volume", "target": ".mist"}]
    return {"name": name, "ports": ports, "compoundNodes": nodes, "connections": conns,
            "values": values}


def writer(nodes, values, conns, ports, short, source, typ):
    """write_NumPy of `source` to the path port path_<short>; ok_<short> out."""
    w = "w_" + short
    nodes += [{"nodeName": w, "nodeType": "File::NumPy::write_NumPy"}]
    values += [{"valueName": w + ".overwrite", "valueType": "bool", "value": "true"},
               {"valueName": w + ".create_directories", "valueType": "bool", "value": "true"},
               {"valueName": w + ".frame", "valueType": "long", "value": "0L"}]
    conns += [{"source": source, "target": w + ".data"},
              {"source": ".path_" + short, "target": w + ".file_path"},
              {"source": w + ".success", "target": ".ok_" + short}]
    ports += [P("path_" + short, "input", "string", ""), P("ok_" + short, "output", "bool")]


def surface_test():
    """bifcmd: the surface at one frame, the mesh's vertex positions written back."""
    ports = [P("liquid_pattern", "input", "string", ""), P("frame", "input", "float", "1f")]
    ports += [P(n, "input", t, d) for n, t, d in KNOBS]
    nodes = [{"nodeName": "surf", "nodeType": "particle_surface"},
             {"nodeName": "t_f3", "valueType": F3},
             {"nodeName": "g_pos", "nodeType": "Geometry::Properties::get_geo_property"}]
    values = [{"valueName": "t_f3.value", "valueType": F3, "value": [f3(0, 0, 0)]},
              {"valueName": "g_pos.property", "valueType": "string", "value": "point_position"}]
    conns = [{"source": "." + p, "target": "surf." + p}
             for p in ("liquid_pattern", "frame") + tuple(k[0] for k in KNOBS)]
    conns += [{"source": "surf.mesh", "target": "g_pos.geometry"},
              {"source": "t_f3.output", "target": "g_pos.type"}]
    writer(nodes, values, conns, ports, "positions", "g_pos.data", F3)
    return {"name": "User::PFFlip::particle_surface_test", "ports": ports,
            "compounds": [surface("particle_surface")], "compoundNodes": nodes,
            "connections": conns, "values": values}


def mist_test():
    """bifcmd: the mist at one frame, its positions, every particle's depth, the mist
    particles' weights, and the spray density sampled at probe positions (path_probes,
    an (N, 3) .npy) for the gate's integral."""
    ports = [P("liquid_pattern", "input", "string", ""), P("radius_pattern", "input", "string", ""),
             P("frame", "input", "float", "1f"), P("path_probes", "input", "string", "")]
    ports += [P(n, "input", t, d) for n, t, d in MIST_KNOBS]
    nodes = [{"nodeName": "mist", "nodeType": "particle_mist"},
             {"nodeName": "t_pos", "valueType": F3},
             {"nodeName": "rprb", "nodeType": "File::NumPy::read_NumPy"},
             {"nodeName": "zero_l", "valueType": "long"},
             {"nodeName": "zero_f", "valueType": "float"},
             {"nodeName": "smp", "nodeType": "Geometry::Query::sample_volume"}]
    values = [{"valueName": "t_pos.value", "valueType": F3, "value": [f3(0, 0, 0)]},
              {"valueName": "zero_l.value", "valueType": "long", "value": "0L"},
              {"valueName": "zero_f.value", "valueType": "float", "value": "0f"},
              {"valueName": "smp.property", "valueType": "string", "value": "voxel_fog_density"},
              {"valueName": "smp.sampler", "valueType": "Geometry::Query::SamplerType", "value": "kLinear"}]
    conns = [{"source": "." + p, "target": "mist." + p}
             for p in ("liquid_pattern", "radius_pattern", "frame") + tuple(k[0] for k in MIST_KNOBS)]
    conns += [{"source": "t_pos.output", "target": "rprb.type"},
              {"source": ".path_probes", "target": "rprb.file_path"},
              {"source": "zero_l.output", "target": "rprb.frame"},
              {"source": "mist.mist", "target": "smp.volume"},
              {"source": "rprb.data", "target": "smp.positions"},
              {"source": "zero_f.output", "target": "smp.default"}]
    writer(nodes, values, conns, ports, "positions", "mist.mist_positions", F3)
    writer(nodes, values, conns, ports, "depth", "mist.depth", "array<float>")
    writer(nodes, values, conns, ports, "weights", "mist.spray_weights", "array<float>")
    writer(nodes, values, conns, ports, "samples", "smp.sampled_data", "array<float>")
    return {"name": "User::PFFlip::particle_mist_test", "ports": ports,
            "compounds": [mist("particle_mist")], "compoundNodes": nodes,
            "connections": conns, "values": values}


if __name__ == "__main__":
    for fname, c in (("particle_surface.json", surface("User::PFFlip::particle_surface")),
                     ("particle_mist.json", mist("User::PFFlip::particle_mist")),
                     ("particle_surface_test.json", surface_test()),
                     ("particle_mist_test.json", mist_test())):
        with open(os.path.join(HERE, fname), "w", newline="\n") as f:
            json.dump(dict(HEADER, compounds=[c]), f, indent=1)
