"""Generates the particle surface and mist graphs, stock nodes only: one frame of a
play cache (maya/play_cache.py: liquid.####.npy) becomes a water surface mesh, the
way the PF-FLIP and ST-FLIP papers render theirs (a union of spheres splatted onto a
grid finer than the simulation's, smoothed, contoured at the zero level), and the
particles that surface leaves out (spray and mist: at one cell per metre a lone
particle is a cloud of droplets, not a ball) become a fog-density volume.

particle_surface.json  the compound User::PFFlip::particle_surface: frame picks the
    cache frame; the escaped droplets (a radius in radius_pattern, radius.####.npy of
    the play cache) are left out of the surface, the papers' split (their spray
    particles are rendered as droplets, never as liquid; a cache without radii has
    none); point_size is the sphere radius per particle in cells (0.5 = the
    papers' radius of half a cell at 8 particles a cell), detail_size the voxel size
    of the surface grid in cells (0.5 = twice the simulation resolution),
    smooth_iterations the Gaussian passes on the level set, smoothing the mesher's
    own smoothing, offset pushes the surface out (+) or in (-). After the Gaussian
    passes, mcf_iterations of mean-curvature flow under ST-FLIP's feature-preserving
    mask (Appendix B; 0 = off, the paper's k_psi = 30): g = 1 - smoothstep(mask_theta,
    mask_zeta; psi0 / (G * psi0 + 1e-10)), psi0 the union-of-spheres density (1
    inside the raw level set, 0 outside, a half-voxel ramp) and G a Gaussian of
    mask_sigma cells (the paper's 2 dx, theta 2, zeta 5), so g is 1 over the bulk,
    whose density is about twice its blur at the surface, and 0 at a lone drop or
    thin sheet, whose density is many times its blur; the flow smooths where g is 1
    and leaves the rest. mask_dilate (voxels, 1) pushes the ramp of the quotient's numerator out
    (the blur stays the paper's), our one departure from the paper's letter: the paper
    evolves the density itself, whose
    held voxels keep the 0.5 contour between them; we evolve the signed distance,
    where the contour sits between a held inside voxel and a free outside one, and
    the gate showed a one-sphere sheet thinning until the mesher drops it at 0 and
    keeping its thickness at 1 (the bulk's mask unchanged). mask_floor raises g
    everywhere (1: the flow unmasked). Stock nodes: psi0 is a field from the raw
    level set (voxel_field, Core::Math on fields), voxelized by field_to_volume on
    the surface grid for the Gaussian (smooth_voxel_property), and g a field into
    smooth_level_set_property's mask_field in MeanCurvature mode. Out comes the mesh,
    drawn by the viewport and rendered by Arnold (render_surface_preview.py). A
    sphere survives voxelization only when its radius spans two voxels, so isolated
    particles drop out of this mesh at the defaults: they are the mist's.

particle_mist.json  the compound User::PFFlip::particle_mist: the surface's own
    smoothed level set (surface_radius, surface_detail, smooth_iterations and the
    mcf and mask knobs as the surface's), sampled at every particle. Negative is inside the water the mesh
    shows; a lone droplet or a one-particle-thick sheet is eroded away by the
    smoothing and samples positive (or the default, far outside, where no level set
    exists). Those above mist_threshold, and every escaped droplet (a radius), are
    mist: the water the mesh dropped and the spray the papers draw as droplets
    (most droplets sit a voxel deep in their own sphere, which the mesher cannot
    contour: on the breaker's frame 7, 14,653 of 16,438 droplets were neither
    surface nor mist before this rule). Each is drawn as the cloud of droplets its water is (PF-FLIP
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

# the feature-preserving flow after the Gaussian passes (ST-FLIP Appendix B)
MCF_KNOBS = (("mcf_iterations", "int", "0L"), ("mask_sigma", "float", "2f"),
             ("mask_theta", "float", "2f"), ("mask_zeta", "float", "5f"),
             ("mask_floor", "float", "0f"), ("mask_dilate", "float", "1f"))
KNOBS = (("point_size", "float", "0.5f"), ("detail_size", "float", "0.5f"),
         ("offset", "float", "0f"), ("smooth_iterations", "int", "1L"),
         ("smoothing", "float", "0.1f")) + MCF_KNOBS
MIST_KNOBS = (("surface_radius", "float", "0.5f"), ("surface_detail", "float", "0.5f"),
              ("smooth_iterations", "int", "1L"),
              ("mist_threshold", "float", "0f"), ("spray_detail", "float", "0.5f"),
              ("cloud_radius", "float", "1f"), ("metres_per_cell", "float", "1f"),
              ("default_radius", "float", "0.0001f"), ("mist_density", "float", "1f"),
              ("spray_epsilon", "float", "10000f")) + MCF_KNOBS
# the cloud constant k_cloud = 1.5 / (8 c_k 4/3 pi) with c_k = 0.25, the linear
# kernel's mean over its sphere (the gate measures it against this)
K_CLOUD = 1.5 / (8.0 * 0.25 * 4.0 / 3.0 * 3.141592653589793)
# the mask's Gaussian filter width: smooth_voxel_property's voxel_filter_width is
# set to this many sigmas. Measured on the gate's ball (2026-10-17): at 3 sigmas the
# blur is that of a Gaussian of 0.76 sigma, at 6 of 0.92 sigma, at 12 of 1.03 sigma
# (the width is the whole filter, not its half); the gate asserts the measurement
MASK_WIDTH_SIGMAS = 12
F3 = "array<Math::float3>"


def reader(nodes, values, conns):
    """The cache frame's liquid positions and droplet radii: .frame, .liquid_pattern
    and .radius_pattern -> rpos.data and rrad.data; the positions without a radius
    (the liquid that is not spray) -> liq.filtered_array."""
    nodes += [{"nodeName": "t_pos", "valueType": F3},
              {"nodeName": "f2i", "nodeType": "Core::Type_Conversion::to_int"},
              {"nodeName": "rpos", "nodeType": "File::NumPy::read_NumPy"},
              {"nodeName": "t_flt", "valueType": "array<float>"},
              {"nodeName": "rrad", "nodeType": "File::NumPy::read_NumPy"},
              {"nodeName": "tiny", "valueType": "float"},
              {"nodeName": "none", "nodeType": "Core::Logic::less"},
              {"nodeName": "liq", "nodeType": "Core::Array::filter_array"}]
    values += [{"valueName": "t_pos.value", "valueType": F3, "value": [f3(0, 0, 0)]},
               {"valueName": "t_flt.value", "valueType": "array<float>", "value": ["0f"]},
               {"valueName": "tiny.value", "valueType": "float", "value": "0.000000001f"}]
    conns += [{"source": ".frame", "target": "f2i.from"},
              {"source": "t_pos.output", "target": "rpos.type"},
              {"source": ".liquid_pattern", "target": "rpos.file_path"},
              {"source": "f2i.int", "target": "rpos.frame"},
              {"source": "t_flt.output", "target": "rrad.type"},
              {"source": ".radius_pattern", "target": "rrad.file_path"},
              {"source": "f2i.int", "target": "rrad.frame"},
              {"source": "rrad.data", "target": "none.first"},
              {"source": "tiny.output", "target": "none.second"},
              {"source": "rpos.data", "target": "liq.array"},
              {"source": "none.output", "target": "liq.boolean_array"}]


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


def fmath(nodes, conns, name, op, first, second, out="output"):
    """A Core::Math node on two scalar fields (the fan-in form, its field overload)."""
    nodes.append({"nodeName": name, "nodeType": "Core::Math::" + op,
                  "multiInPortNames": ["a", "b"]})
    conns += [{"source": first, "target": name + ".first.a"},
              {"source": second, "target": name + ".first.b"}]
    return name + "." + out


def ffloat(nodes, values, conns, name, source=None, value=None):
    """A scalar field from a float port (`source`) or a constant; returns it."""
    if source is None:
        nodes.append({"nodeName": name + "_v", "valueType": "float"})
        values.append({"valueName": name + "_v.value", "valueType": "float", "value": value})
        source = name + "_v.output"
    nodes.append({"nodeName": name, "nodeType": "Core::Fields::to_field"})
    conns.append({"source": source, "target": name + ".value"})
    return name + ".out_field"


def mcf(nodes, values, conns, vol, raw, pts, detail):
    """ST-FLIP Appendix B after the Gaussian passes: .mcf_iterations of mean-curvature
    flow on the level set `vol` under the mask g = 1 - smoothstep(theta, zeta;
    psi0 / (G * psi0 + 1e-10)), psi0 the union-of-spheres density (1 inside the raw
    level set `raw`, 0 outside, a half-voxel ramp, pushed out by mask_dilate voxels)
    and G a Gaussian of mask_sigma cells; g floored at mask_floor. Returns the
    smoothed volume's port (with mcf_iterations 0 the volume as it was), the density
    volume, its blur and the mask field."""
    # psi0 as a field from the raw level set, clamp(0.5 - sdf / voxel, 0, 1), for the
    # blur; the numerator pushed out by mask_dilate voxels: clamp(0.5 + dilate - sdf / voxel)
    nodes += [{"nodeName": "k_sf", "nodeType": "Core::Fields::scalar_field"},
              {"nodeName": "k_sdf", "nodeType": "Core::Fields::voxel_field"},
              {"nodeName": "k_hv", "valueType": "float"},
              {"nodeName": "k_half", "nodeType": "Core::Math::add", "multiInPortNames": ["a", "b"]},
              {"nodeName": "k_psi0", "nodeType": "Core::Math::clamp"},
              {"nodeName": "k_psi", "nodeType": "Core::Math::clamp"}]
    values += [{"valueName": "k_sdf.property", "valueType": "string", "value": "voxel_signed_distance"},
               {"valueName": "k_hv.value", "valueType": "float", "value": "0.5f"},
               {"valueName": "k_psi0.min", "valueType": "float", "value": "0f"},
               {"valueName": "k_psi0.max", "valueType": "float", "value": "1f"},
               {"valueName": "k_psi.min", "valueType": "float", "value": "0f"},
               {"valueName": "k_psi.max", "valueType": "float", "value": "1f"}]
    conns += [{"source": "k_sf.scalar_field", "target": "k_sdf.field_type"},
              {"source": raw, "target": "k_sdf.volume"},
              {"source": ".mask_dilate", "target": "k_half.first.a"},
              {"source": "k_hv.output", "target": "k_half.first.b"}]
    sdf_vox = fmath(nodes, conns, "k_sv", "divide", "k_sdf.voxel_field",
                    ffloat(nodes, values, conns, "k_detf", source=detail))
    conns.append({"source": fmath(nodes, conns, "k_p00", "subtract",
                                  ffloat(nodes, values, conns, "k_hvf", source="k_hv.output"),
                                  sdf_vox), "target": "k_psi0.value"})
    conns.append({"source": fmath(nodes, conns, "k_p0", "subtract",
                                  ffloat(nodes, values, conns, "k_halff", source="k_half.output"),
                                  sdf_vox), "target": "k_psi.value"})
    psi0, psi = "k_psi0.clamped", "k_psi.clamped"
    # its blur: psi0 voxelized on the surface grid within the raw level set's bounds
    # (volume_bounds; its band reaches past the outermost spheres), then a Gaussian of mask_sigma
    # / detail voxels, the filter MASK_WIDTH_SIGMAS sigmas wide
    nodes += [{"nodeName": "k_bnd", "nodeType": "Geometry::Volume::volume_bounds"},
              {"nodeName": "k_vol", "nodeType": "Geometry::Converters::field_to_volume"},
              {"nodeName": "k_sig", "nodeType": "Core::Math::divide", "multiInPortNames": ["s", "h"]},
              {"nodeName": "k_three", "valueType": "float"},
              {"nodeName": "k_w", "nodeType": "Core::Math::multiply", "multiInPortNames": ["a", "b"]},
              {"nodeName": "k_wi", "nodeType": "Core::Type_Conversion::to_int"},
              {"nodeName": "k_blur", "nodeType": "Geometry::Volume::smooth_voxel_property"},
              {"nodeName": "k_fb", "nodeType": "Core::Fields::voxel_field"}]
    values += [{"valueName": "k_three.value", "valueType": "float", "value": f"{MASK_WIDTH_SIGMAS}f"},
               {"valueName": "k_vol.field_interpretation",
                "valueType": "Geometry::Volume::SamplerInterpretation", "value": "FogDensity"},
               {"valueName": "k_vol.resolution_mode", "valueType": "Geometry::Volume::ResolutionType",
                "value": "Absolute"},
               {"valueName": "k_vol.store_level_set", "valueType": "bool", "value": "false"},
               {"valueName": "k_vol.store_fog_density", "valueType": "bool", "value": "true"},
               {"valueName": "k_blur.property", "valueType": "string", "value": "voxel_fog_density"},
               {"valueName": "k_fb.property", "valueType": "string", "value": "voxel_fog_density"}]
    conns += [{"source": psi0, "target": "k_vol.field"},
              {"source": raw, "target": "k_bnd.volume"},
              {"source": "k_bnd.min_bound", "target": "k_vol.min_bounds"},
              {"source": "k_bnd.max_bound", "target": "k_vol.max_bounds"},
              {"source": detail, "target": "k_vol.detail_size"},
              {"source": ".mask_sigma", "target": "k_sig.first.s"},
              {"source": detail, "target": "k_sig.first.h"},
              {"source": "k_sig.output", "target": "k_w.first.a"},
              {"source": "k_three.output", "target": "k_w.first.b"},
              {"source": "k_w.output", "target": "k_wi.from"},
              {"source": "k_vol.volume", "target": "k_blur.volume"},
              {"source": "k_sig.output", "target": "k_blur.voxel_standard_deviation"},
              {"source": "k_wi.int", "target": "k_blur.voxel_filter_width"},
              {"source": "k_sf.scalar_field", "target": "k_fb.field_type"},
              {"source": "k_blur.out_volume", "target": "k_fb.volume"}]
    dens = "k_vol.volume"
    # q = psi0 / (blur + 1e-10); s = clamp((q - theta) / (zeta - theta), 0, 1);
    # g = 1 - s^2 (3 - 2 s), the cubic smoothstep; floored at mask_floor
    eps = ffloat(nodes, values, conns, "k_eps", value="0.0000000001f")
    den = fmath(nodes, conns, "k_den", "add", "k_fb.voxel_field", eps)
    q = fmath(nodes, conns, "k_q", "divide", psi, den)
    theta = ffloat(nodes, values, conns, "k_th", source=".mask_theta")
    nodes += [{"nodeName": "k_span", "nodeType": "Core::Math::subtract", "multiInPortNames": ["z", "t"]},
              {"nodeName": "k_s", "nodeType": "Core::Math::clamp"}]
    values += [{"valueName": "k_s.min", "valueType": "float", "value": "0f"},
               {"valueName": "k_s.max", "valueType": "float", "value": "1f"}]
    conns += [{"source": ".mask_zeta", "target": "k_span.first.z"},
              {"source": ".mask_theta", "target": "k_span.first.t"}]
    span = ffloat(nodes, values, conns, "k_spanf", source="k_span.output")
    qs = fmath(nodes, conns, "k_qs", "subtract", q, theta)
    conns.append({"source": fmath(nodes, conns, "k_s0", "divide", qs, span), "target": "k_s.value"})
    s2 = fmath(nodes, conns, "k_s2", "multiply", "k_s.clamped", "k_s.clamped")
    s3 = fmath(nodes, conns, "k_s3", "multiply", s2, "k_s.clamped")
    a = fmath(nodes, conns, "k_a", "multiply", ffloat(nodes, values, conns, "k_c3", value="3f"), s2)
    b = fmath(nodes, conns, "k_b", "multiply", ffloat(nodes, values, conns, "k_c2", value="2f"), s3)
    poly = fmath(nodes, conns, "k_poly", "subtract", a, b)
    g = fmath(nodes, conns, "k_g", "subtract", ffloat(nodes, values, conns, "k_c1", value="1f"), poly)
    gf = fmath(nodes, conns, "k_gf", "max", g,
               ffloat(nodes, values, conns, "k_fl", source=".mask_floor"), out="maximum")
    nodes.append({"nodeName": "mcf", "nodeType": "Geometry::Volume::smooth_level_set_property"})
    values += [{"valueName": "mcf.mode", "valueType": "Geometry::Volume::SmoothLevelSetPropertyMode",
                "value": "MeanCurvature"},
               {"valueName": "mcf.mask_interpretation",
                "valueType": "Geometry::Volume::SamplerInterpretation", "value": "FogDensity"},
               {"valueName": "mcf.grid_retopology", "valueType": "bool", "value": "false"}]
    conns += [{"source": vol, "target": "mcf.volume"},
              {"source": ".mcf_iterations", "target": "mcf.iterations"},
              {"source": gf, "target": "mcf.mask_field"}]
    # behind an if on mcf_iterations > 0: with the flow off the volume passes
    # through and none of the mask is built (an if evaluates the branch it takes)
    nodes += [{"nodeName": "k_zero", "valueType": "int"},
              {"nodeName": "k_on", "nodeType": "Core::Logic::greater"},
              {"nodeName": "k_empty", "valueType": "Object"},
              {"nodeName": "k_out", "nodeType": "Core::Logic::if"},
              {"nodeName": "k_if_d", "nodeType": "Core::Logic::if"},
              {"nodeName": "k_if_b", "nodeType": "Core::Logic::if"},
              {"nodeName": "k_if_g", "nodeType": "Core::Logic::if"}]
    values += [{"valueName": "k_zero.value", "valueType": "int", "value": "0"},
               {"valueName": "k_empty.value", "valueType": "Object", "value": {}}]
    conns += [{"source": ".mcf_iterations", "target": "k_on.first"},
              {"source": "k_zero.output", "target": "k_on.second"}]
    for name, yes, no in (("k_out", "mcf.out_volume", vol), ("k_if_d", dens, "k_empty.output"),
                          ("k_if_b", "k_blur.out_volume", "k_empty.output"),
                          ("k_if_g", gf, "k_sf.scalar_field")):
        conns += [{"source": "k_on.output", "target": name + ".condition"},
                  {"source": yes, "target": name + ".true_case"},
                  {"source": no, "target": name + ".false_case"}]
    return "k_out.output", "k_if_d.output", "k_if_b.output", "k_if_g.output"


def surface(name):
    ports = [P("liquid_pattern", "input", "string", ""), P("radius_pattern", "input", "string", ""),
             P("frame", "input", "float", "1f")]
    ports += [P(n, "input", t, d) for n, t, d in KNOBS]
    ports += [P("mesh", "output", "Object"), P("mask", "output", "Core::Fields::ScalarField"),
              P("mask_density", "output", "Object"), P("mask_blur", "output", "Object")]
    nodes, values, conns = [], [], []
    reader(nodes, values, conns)
    pts = spheres(nodes, values, conns, "s", "liq.filtered_array", ".point_size")
    vol = splat(nodes, values, conns, "p2v", pts, ".detail_size", level_set=True)
    nodes += [{"nodeName": "sm", "nodeType": "Geometry::Volume::smooth_level_set_property"},
              {"nodeName": "v2m", "nodeType": "Geometry::Converters::volume_to_mesh",
               "fanInPortNames": {"volumes": ["vol"]}},
              {"nodeName": "zero", "valueType": "long"},
              {"nodeName": "first", "nodeType": "Core::Array::get_from_array"}]
    values += [{"valueName": "zero.value", "valueType": "long", "value": "0L"}]
    smoothed, dens, blur, g = mcf(nodes, values, conns, "sm.out_volume", vol, pts, ".detail_size")
    conns += [{"source": dens, "target": ".mask_density"}, {"source": blur, "target": ".mask_blur"},
              {"source": g, "target": ".mask"},
              {"source": ".offset", "target": "p2v.offset"},
              {"source": vol, "target": "sm.volume"},
              {"source": ".smooth_iterations", "target": "sm.iterations"},
              {"source": smoothed, "target": "v2m.volumes.vol"},
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
    # every particle's depth: the surface's smoothed level set (of the particles
    # without a radius, as the surface's) sampled at the particle (far from any
    # surface a sample is the default, well outside)
    pts = spheres(nodes, values, conns, "c", "liq.filtered_array", ".surface_radius")
    vol = splat(nodes, values, conns, "p2v_c", pts, ".surface_detail", level_set=True)
    nodes += [{"nodeName": "sm", "nodeType": "Geometry::Volume::smooth_level_set_property"},
              {"nodeName": "far", "valueType": "float"},
              {"nodeName": "smp", "nodeType": "Geometry::Query::sample_volume"}]
    values += [{"valueName": "far.value", "valueType": "float", "value": "10f"},
               {"valueName": "smp.property", "valueType": "string", "value": "voxel_signed_distance"}]
    smoothed = mcf(nodes, values, conns, "sm.out_volume", vol, pts, ".surface_detail")[0]
    conns += [{"source": vol, "target": "sm.volume"},
              {"source": ".smooth_iterations", "target": "sm.iterations"},
              {"source": smoothed, "target": "smp.volume"},
              {"source": "rpos.data", "target": "smp.positions"},
              {"source": "far.output", "target": "smp.default"}]
    depth = "smp.sampled_data"
    nodes += [{"nodeName": "thin", "nodeType": "Core::Logic::greater"},
              {"nodeName": "hasr", "nodeType": "Core::Logic::greater"},
              {"nodeName": "ismist", "nodeType": "Core::Logic::or"},
              {"nodeName": "pick", "nodeType": "Core::Array::filter_array"},
              {"nodeName": "pickr", "nodeType": "Core::Array::filter_array"}]
    conns += [{"source": depth, "target": ".depth"},
              # mist = depth > threshold, or a droplet radius; its positions and radii
              {"source": depth, "target": "thin.first"},
              {"source": ".mist_threshold", "target": "thin.second"},
              {"source": "rrad.data", "target": "hasr.first"},
              {"source": "tiny.output", "target": "hasr.second"},
              {"source": "thin.output", "target": "ismist.first"},
              {"source": "hasr.output", "target": "ismist.second"},
              {"source": "rpos.data", "target": "pick.array"},
              {"source": "ismist.output", "target": "pick.boolean_array"},
              {"source": "pick.filtered_array", "target": ".mist_positions"},
              {"source": "rrad.data", "target": "pickr.array"},
              {"source": "ismist.output", "target": "pickr.boolean_array"}]
    # r_eff = r + default_radius where r is 0 (a mist particle without a droplet radius)
    nodes += [{"nodeName": "none_m", "nodeType": "Core::Logic::less"},
              {"nodeName": "none_f", "nodeType": "Core::Type_Conversion::to_float"},
              {"nodeName": "fill", "nodeType": "Core::Math::multiply", "multiInPortNames": ["a", "b"]},
              {"nodeName": "r_eff", "nodeType": "Core::Math::add", "multiInPortNames": ["a", "b"]}]
    conns += [{"source": "pickr.filtered_array", "target": "none_m.first"},
              {"source": "tiny.output", "target": "none_m.second"},
              {"source": "none_m.output", "target": "none_f.from"},
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
    """bifcmd: the surface at one frame, the mesh's vertex positions written back, and
    the mask's density, its blur and the mask sampled at probe positions (path_probes,
    an (N, 3) .npy)."""
    ports = [P("liquid_pattern", "input", "string", ""), P("radius_pattern", "input", "string", ""),
             P("frame", "input", "float", "1f"), P("path_probes", "input", "string", "")]
    ports += [P(n, "input", t, d) for n, t, d in KNOBS]
    nodes = [{"nodeName": "surf", "nodeType": "particle_surface"},
             {"nodeName": "t_f3", "valueType": F3},
             {"nodeName": "g_pos", "nodeType": "Geometry::Properties::get_geo_property"},
             {"nodeName": "rprb", "nodeType": "File::NumPy::read_NumPy"},
             {"nodeName": "zero_l", "valueType": "long"},
             {"nodeName": "zero_f", "valueType": "float"},
             {"nodeName": "smp_d", "nodeType": "Geometry::Query::sample_volume"},
             {"nodeName": "smp_b", "nodeType": "Geometry::Query::sample_volume"},
             {"nodeName": "smp_m", "nodeType": "Core::Fields::sample_field"}]
    values = [{"valueName": "t_f3.value", "valueType": F3, "value": [f3(0, 0, 0)]},
              {"valueName": "g_pos.property", "valueType": "string", "value": "point_position"},
              {"valueName": "zero_l.value", "valueType": "long", "value": "0L"},
              {"valueName": "zero_f.value", "valueType": "float", "value": "0f"}]
    for s in ("smp_d", "smp_b"):
        values += [{"valueName": s + ".property", "valueType": "string", "value": "voxel_fog_density"},
                   {"valueName": s + ".sampler", "valueType": "Geometry::Query::SamplerType",
                    "value": "kLinear"}]
    conns = [{"source": "." + p, "target": "surf." + p}
             for p in ("liquid_pattern", "radius_pattern", "frame") + tuple(k[0] for k in KNOBS)]
    conns += [{"source": "surf.mesh", "target": "g_pos.geometry"},
              {"source": "t_f3.output", "target": "g_pos.type"},
              {"source": "t_f3.output", "target": "rprb.type"},
              {"source": ".path_probes", "target": "rprb.file_path"},
              {"source": "zero_l.output", "target": "rprb.frame"},
              {"source": "surf.mask_density", "target": "smp_d.volume"},
              {"source": "surf.mask_blur", "target": "smp_b.volume"},
              {"source": "rprb.data", "target": "smp_d.positions"},
              {"source": "rprb.data", "target": "smp_b.positions"},
              {"source": "zero_f.output", "target": "smp_d.default"},
              {"source": "zero_f.output", "target": "smp_b.default"},
              {"source": "surf.mask", "target": "smp_m.field"},
              {"source": "rprb.data", "target": "smp_m.position"}]
    writer(nodes, values, conns, ports, "positions", "g_pos.data", F3)
    writer(nodes, values, conns, ports, "density", "smp_d.sampled_data", "array<float>")
    writer(nodes, values, conns, ports, "blur", "smp_b.sampled_data", "array<float>")
    writer(nodes, values, conns, ports, "mask", "smp_m.result", "array<float>")
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
