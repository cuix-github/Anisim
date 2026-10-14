"""Generates the mesh-collider graphs (stage 3 of mesh colliders), stock nodes only.

collider_field.json  the compound User::PFFlip::collider_field: from the objects the
    shipped Simulation::Common::collider node makes to our solver's solid_phi, the
    solid's signed distance at the (nx+1)(ny+1)(nz+1) grid nodes in cells, negative
    inside, in the node's order ((i*(ny+1)+j)*(nz+1)+k). Per collider it reads the
    settings the shipped node stored (its geometry under "geo", then geo_detail_size,
    resolution_mode, volume_offset, volume_thickening, volume_mode, min_hole_radius,
    enable), voxelizes with the same public chain the shipped liquid solver uses
    (compute_voxel_scale, then mesh_to_volume as a level set), samples the result at
    the grid nodes, divides by the cell size and keeps the minimum over colliders.
    The lattice is built in exact integer arithmetic (long indices), so any grid
    size works. For Maya, list this file in maya/pfflip_graphs_config.json.

collider_setup.json  the bifcmd wrapper around it: one or two meshes read from .npy
    (point positions, face vertices, face offsets), the shipped collider node with
    its detail size, offset, thickening and hole radius as ports, the adapter, and
    solid_phi written as .npy, the loop graphs' path_solid_phi. A static solid is
    computed once this way and read by every run. graphs/colliders.py drives it.

Regenerate after editing:  python generate_collider_graph.py
"""
import json
import os

from generate_sim_graph import P

HERE = os.path.dirname(os.path.abspath(__file__))
HEADER = {"header": {"metadata": [{"metaName": "adskFileFormatVersion", "metaValue": "100L"}]},
          "namespaces": [], "types": []}
FAR = "1000000000f"  # the distance before any collider is seen: far outside


def multi(name, op, a="a", b="b"):
    return {"nodeName": name, "nodeType": "Core::Math::" + op, "multiInPortNames": [a, b]}


def per_collider():
    """The iterate body: one collider's level set sampled at the nodes, min'ed into
    the running field. Disabled colliders pass the field through."""
    nodes = [
        {"nodeName": "col", "nodeType": "Core::Array::get_from_array"},
        {"nodeName": "obj_t", "valueType": "Object"},
        {"nodeName": "rmode_t", "valueType": "Geometry::Volume::ResolutionType"},
        {"nodeName": "vmode_t", "valueType": "Geometry::Converters::VolumeMode"},
        {"nodeName": "geo_arr", "nodeType": "Core::Array::build_array", "multiInPortNames": ["g"]},
        {"nodeName": "vscale", "nodeType": "Geometry::Volume::compute_voxel_scale"},
        multi("detail", "multiply"),
        {"nodeName": "m2v", "nodeType": "Geometry::Converters::mesh_to_volume"},
        {"nodeName": "samp", "nodeType": "Geometry::Query::sample_volume"},
        {"nodeName": "far", "valueType": "float"},
        multi("cells", "divide"),
        multi("mn", "min"),
        {"nodeName": "sel", "nodeType": "Core::Logic::if"},
    ]
    values = [
        {"valueName": "obj_t.value", "valueType": "Object", "value": {}},
        {"valueName": "rmode_t.value", "valueType": "Geometry::Volume::ResolutionType",
         "value": "Absolute"},
        {"valueName": "vmode_t.value", "valueType": "Geometry::Converters::VolumeMode",
         "value": "Solid"},
        {"valueName": "far.value", "valueType": "float", "value": FAR},
        {"valueName": "m2v.store_level_set", "valueType": "bool", "value": "true"},
        {"valueName": "m2v.store_fog_density", "valueType": "bool", "value": "false"},
        {"valueName": "m2v.resolution_mode", "valueType": "Geometry::Volume::ResolutionType",
         "value": "Absolute"},
        {"valueName": "m2v.properties", "valueType": "string", "value": ""},
        {"valueName": "samp.property", "valueType": "string", "value": "voxel_signed_distance"},
        {"valueName": "samp.sampler", "valueType": "Geometry::Query::SamplerType",
         "value": "kLinear"},
    ]
    conns = [
        {"source": ".colliders", "target": "col.array"},
        {"source": ".current_index", "target": "col.index"},
    ]
    # the properties the shipped collider node stores on each collider object
    props = [("geo", "geo", "obj_t.output", None),
             ("det", "geo_detail_size", None, ("float", "0.0500000007f")),
             ("rmode", "resolution_mode", "rmode_t.output", None),
             ("off", "volume_offset", None, ("float", "0f")),
             ("thick", "volume_thickening", None, ("float", "1f")),
             ("vmode", "volume_mode", "vmode_t.output", None),
             ("hole", "min_hole_radius", None, ("float", "0f")),
             ("en", "enable", None, ("bool", "true"))]
    for name, key, typed_src, dflt in props:
        nodes.append({"nodeName": "p_" + name, "nodeType": "Core::Object::get_property"})
        values.append({"valueName": f"p_{name}.key", "valueType": "string", "value": key})
        conns.append({"source": "col.value", "target": f"p_{name}.object"})
        if typed_src:
            conns.append({"source": typed_src, "target": f"p_{name}.default_and_type"})
        else:
            values.append({"valueName": f"p_{name}.default_and_type", "valueType": dflt[0],
                           "value": dflt[1]})
    conns += [
        {"source": "p_geo.value", "target": "geo_arr.first.g"},
        {"source": "geo_arr.array", "target": "vscale.geo"},
        {"source": "p_rmode.value", "target": "vscale.method"},
        {"source": "p_det.value", "target": "detail.first.a"},
        {"source": "vscale.voxel_scale", "target": "detail.first.b"},
        {"source": "p_geo.value", "target": "m2v.mesh"},
        {"source": "detail.output", "target": "m2v.detail_size"},
        {"source": "p_off.value", "target": "m2v.offset"},
        {"source": "p_thick.value", "target": "m2v.thickening"},
        {"source": "p_vmode.value", "target": "m2v.volume_mode"},
        {"source": "p_hole.value", "target": "m2v.min_hole_radius"},
        {"source": "m2v.volume", "target": "samp.volume"},
        {"source": ".positions", "target": "samp.positions"},
        {"source": "far.output", "target": "samp.default"},
        {"source": "samp.sampled_data", "target": "cells.first.a"},
        {"source": ".cell_size", "target": "cells.first.b"},
        {"source": ".phi_in", "target": "mn.first.a"},
        {"source": "cells.output", "target": "mn.first.b"},
        {"source": "p_en.value", "target": "sel.condition"},
        {"source": "mn.minimum", "target": "sel.true_case"},
        {"source": ".phi_in", "target": "sel.false_case"},
        {"source": "sel.output", "target": ".phi_out"},
    ]
    return {"name": "per_collider",
            "ports": [P("colliders", "input", "array<Object>"),
                      P("positions", "input", "array<Math::float3>"),
                      P("cell_size", "input", "float"),
                      P("phi_in", "input", "array<float>"), P("phi_out", "output", "array<float>"),
                      P("iteration_limit", "input", "long", portIterationLimit="true"),
                      P("current_index", "input", "long", portIterationCounter="true")],
            "compounds": [], "compoundNodes": nodes, "connections": conns, "values": values,
            "iterateCompound": {"ports": [
                {"portKind": "state", "inputPortName": "phi_in", "outputPortName": "phi_out"}]}}


def adapter(name):
    """User::PFFlip::collider_field (or an inline copy under a local name)."""
    nodes = [
        {"nodeName": "one", "valueType": "int"},
        {"nodeName": "far", "valueType": "float"},
        {"nodeName": "far_arr", "valueType": "array<float>"},
        multi("nx1", "add"), multi("ny1", "add"), multi("nz1", "add"),
        multi("nyz", "multiply"), multi("n", "multiply"),
        {"nodeName": "n_long", "nodeType": "Core::Type_Conversion::to_long"},
        {"nodeName": "filled", "nodeType": "Core::Array::resize_array"},
        {"nodeName": "idx", "nodeType": "Core::Array::get_array_indices"},
        {"nodeName": "nz1_l", "nodeType": "Core::Type_Conversion::to_long"},
        {"nodeName": "ny1_l", "nodeType": "Core::Type_Conversion::to_long"},
        {"nodeName": "nyz_l", "nodeType": "Core::Type_Conversion::to_long"},
        {"nodeName": "k", "nodeType": "Core::Math::modulo"},
        multi("ij", "divide"),
        {"nodeName": "j", "nodeType": "Core::Math::modulo"},
        multi("i", "divide"),
        {"nodeName": "i_f", "nodeType": "Core::Type_Conversion::to_float"},
        {"nodeName": "j_f", "nodeType": "Core::Type_Conversion::to_float"},
        {"nodeName": "k_f", "nodeType": "Core::Type_Conversion::to_float"},
        {"nodeName": "ijk", "nodeType": "Core::Conversion::scalar_to_vector3"},
        multi("scaled", "multiply"), multi("pos", "add"),
        {"nodeName": "count", "nodeType": "Core::Array::array_size"},
        {"nodeName": "each", "nodeType": "per_collider"},
    ]
    values = [
        {"valueName": "one.value", "valueType": "int", "value": "1"},
        {"valueName": "far.value", "valueType": "float", "value": FAR},
        {"valueName": "far_arr.value", "valueType": "array<float>", "value": [FAR]},
    ]
    conns = [
        # the node lattice, in the solver's index order, from exact long indices
        {"source": ".nx", "target": "nx1.first.a"}, {"source": "one.output", "target": "nx1.first.b"},
        {"source": ".ny", "target": "ny1.first.a"}, {"source": "one.output", "target": "ny1.first.b"},
        {"source": ".nz", "target": "nz1.first.a"}, {"source": "one.output", "target": "nz1.first.b"},
        {"source": "ny1.output", "target": "nyz.first.a"}, {"source": "nz1.output", "target": "nyz.first.b"},
        {"source": "nx1.output", "target": "n.first.a"}, {"source": "nyz.output", "target": "n.first.b"},
        {"source": "n.output", "target": "n_long.from"},
        {"source": "far_arr.output", "target": "filled.array"},
        {"source": "n_long.long", "target": "filled.new_size"},
        {"source": "far.output", "target": "filled.default"},
        {"source": "filled.resized", "target": "idx.array"},
        {"source": "nz1.output", "target": "nz1_l.from"},
        {"source": "ny1.output", "target": "ny1_l.from"},
        {"source": "nyz.output", "target": "nyz_l.from"},
        {"source": "idx.indices", "target": "k.value"}, {"source": "nz1_l.long", "target": "k.divisor"},
        {"source": "idx.indices", "target": "ij.first.a"}, {"source": "nz1_l.long", "target": "ij.first.b"},
        {"source": "ij.output", "target": "j.value"}, {"source": "ny1_l.long", "target": "j.divisor"},
        {"source": "idx.indices", "target": "i.first.a"}, {"source": "nyz_l.long", "target": "i.first.b"},
        {"source": "i.output", "target": "i_f.from"},
        {"source": "j.remainder", "target": "j_f.from"},
        {"source": "k.remainder", "target": "k_f.from"},
        {"source": "i_f.float", "target": "ijk.x"}, {"source": "j_f.float", "target": "ijk.y"},
        {"source": "k_f.float", "target": "ijk.z"},
        {"source": "ijk.vector3", "target": "scaled.first.a"},
        {"source": ".cell_size", "target": "scaled.first.b"},
        {"source": "scaled.output", "target": "pos.first.a"}, {"source": ".origin", "target": "pos.first.b"},
        {"source": "pos.output", "target": ".node_positions"},
        {"source": "n_long.long", "target": ".node_count"},
        # the colliders, one after another, into the running minimum
        {"source": ".colliders", "target": "count.array"},
        {"source": "count.size", "target": "each.iteration_limit"},
        {"source": ".colliders", "target": "each.colliders"},
        {"source": "pos.output", "target": "each.positions"},
        {"source": ".cell_size", "target": "each.cell_size"},
        {"source": "filled.resized", "target": "each.phi_in"},
        {"source": "each.phi_out", "target": ".solid_phi"},
    ]
    return {"name": name,
            "ports": [P("colliders", "input", "array<Object>"),
                      P("nx", "input", "int", "32"), P("ny", "input", "int", "32"),
                      P("nz", "input", "int", "32"),
                      P("cell_size", "input", "float", "1f"),
                      P("origin", "input", "Math::float3", {"x": "0f", "y": "0f", "z": "0f"}),
                      P("solid_phi", "output", "array<float>"),
                      P("node_positions", "output", "array<Math::float3>"),
                      P("node_count", "output", "long")],
            "compounds": [per_collider()], "compoundNodes": nodes, "connections": conns,
            "values": values}


def mesh_reader(suffix):
    """Nodes reading one mesh from three .npy files into construct_mesh."""
    s = suffix
    nodes = [
        {"nodeName": "r_pts" + s, "nodeType": "File::NumPy::read_NumPy"},
        {"nodeName": "r_fv" + s, "nodeType": "File::NumPy::read_NumPy"},
        {"nodeName": "r_fo" + s, "nodeType": "File::NumPy::read_NumPy"},
        {"nodeName": "fv_u" + s, "nodeType": "Core::Type_Conversion::to_unsigned_int"},
        {"nodeName": "fo_u" + s, "nodeType": "Core::Type_Conversion::to_unsigned_int"},
        {"nodeName": "mesh" + s, "nodeType": "Geometry::Mesh::construct_mesh"},
    ]
    conns = [
        {"source": "t_pos.output", "target": "r_pts" + s + ".type"},
        {"source": "t_long.output", "target": "r_fv" + s + ".type"},
        {"source": "t_long.output", "target": "r_fo" + s + ".type"},
        {"source": ".path_points" + s, "target": "r_pts" + s + ".file_path"},
        {"source": ".path_face_vertex" + s, "target": "r_fv" + s + ".file_path"},
        {"source": ".path_face_offset" + s, "target": "r_fo" + s + ".file_path"},
        {"source": "r_fv" + s + ".data", "target": "fv_u" + s + ".from"},
        {"source": "r_fo" + s + ".data", "target": "fo_u" + s + ".from"},
        {"source": "r_pts" + s + ".data", "target": "mesh" + s + ".point_position"},
        {"source": "fv_u" + s + ".unsigned_int", "target": "mesh" + s + ".face_vertex"},
        {"source": "fo_u" + s + ".unsigned_int", "target": "mesh" + s + ".face_offset"},
    ]
    ports = [P("path_points" + s, "input", "string", ""),
             P("path_face_vertex" + s, "input", "string", ""),
             P("path_face_offset" + s, "input", "string", "")]
    return nodes, conns, ports


def setup():
    nodes = [
        {"nodeName": "t_pos", "valueType": "array<Math::float3>"},
        {"nodeName": "t_long", "valueType": "array<long>"},
        {"nodeName": "obj_t", "valueType": "Object"},
        {"nodeName": "one_l", "valueType": "long"},
        {"nodeName": "two_l", "valueType": "long"},
        {"nodeName": "geos", "nodeType": "Core::Array::build_array", "multiInPortNames": ["a", "b"]},
        {"nodeName": "n_geo", "nodeType": "Core::Logic::if"},
        {"nodeName": "geos_n", "nodeType": "Core::Array::resize_array"},
        {"nodeName": "collider", "nodeType": "Simulation::Common::collider"},
        {"nodeName": "field", "nodeType": "collider_field"},
        {"nodeName": "origin", "nodeType": "Core::Conversion::scalar_to_vector3"},
        {"nodeName": "w_phi", "nodeType": "File::NumPy::write_NumPy"},
        {"nodeName": "w_pos", "nodeType": "File::NumPy::write_NumPy"},
    ]
    values = [
        {"valueName": "t_pos.value", "valueType": "array<Math::float3>",
         "value": [{"x": "0f", "y": "0f", "z": "0f"}]},
        {"valueName": "t_long.value", "valueType": "array<long>", "value": ["0L"]},
        {"valueName": "obj_t.value", "valueType": "Object", "value": {}},
        {"valueName": "one_l.value", "valueType": "long", "value": "1L"},
        {"valueName": "two_l.value", "valueType": "long", "value": "2L"},
        # the shipped node's own settings we keep fixed here: distances in scene units,
        # a solid (not a shell), enabled, static
        {"valueName": "collider.resolution_mode", "valueType": "Geometry::Volume::ResolutionType",
         "value": "Absolute"},
        {"valueName": "collider.geo_mode", "valueType": "Geometry::Converters::VolumeMode",
         "value": "Solid"},
        {"valueName": "collider.enable", "valueType": "bool", "value": "true"},
    ]
    for w in ("w_phi", "w_pos"):
        values += [{"valueName": w + ".enable", "valueType": "bool", "value": "true"},
                   {"valueName": w + ".overwrite", "valueType": "bool", "value": "true"},
                   {"valueName": w + ".create_directories", "valueType": "bool", "value": "true"},
                   {"valueName": w + ".frame", "valueType": "long", "value": "0L"}]
    ports = [P("use_b", "input", "bool", "false"),
             P("detail_size", "input", "float", "1f"), P("offset", "input", "float", "0f"),
             P("thickening", "input", "float", "0f"), P("min_hole_radius", "input", "float", "0f"),
             P("nx", "input", "int", "32"), P("ny", "input", "int", "32"), P("nz", "input", "int", "32"),
             P("cell_size", "input", "float", "1f"),
             # the origin as three floats: bifcmd's --set-port takes scalars
             P("origin_x", "input", "float", "0f"), P("origin_y", "input", "float", "0f"),
             P("origin_z", "input", "float", "0f"),
             P("path_solid_phi", "input", "string", ""), P("path_positions", "input", "string", ""),
             P("ok", "output", "bool"), P("ok_positions", "output", "bool"),
             P("node_count", "output", "long")]
    conns = [
        {"source": "mesh.mesh", "target": "geos.first.a"},
        {"source": "mesh_b.mesh", "target": "geos.first.b"},
        {"source": ".use_b", "target": "n_geo.condition"},
        {"source": "two_l.output", "target": "n_geo.true_case"},
        {"source": "one_l.output", "target": "n_geo.false_case"},
        {"source": "geos.array", "target": "geos_n.array"},
        {"source": "n_geo.output", "target": "geos_n.new_size"},
        {"source": "obj_t.output", "target": "geos_n.default"},
        {"source": "geos_n.resized", "target": "collider.geometry"},
        {"source": ".detail_size", "target": "collider.detail_size"},
        {"source": ".offset", "target": "collider.offset"},
        {"source": ".thickening", "target": "collider.thickening"},
        {"source": ".min_hole_radius", "target": "collider.min_hole_radius"},
        {"source": "collider.collider", "target": "field.colliders"},
        {"source": ".nx", "target": "field.nx"}, {"source": ".ny", "target": "field.ny"},
        {"source": ".nz", "target": "field.nz"},
        {"source": ".cell_size", "target": "field.cell_size"},
        {"source": ".origin_x", "target": "origin.x"}, {"source": ".origin_y", "target": "origin.y"},
        {"source": ".origin_z", "target": "origin.z"},
        {"source": "origin.vector3", "target": "field.origin"},
        {"source": "field.solid_phi", "target": "w_phi.data"},
        {"source": ".path_solid_phi", "target": "w_phi.file_path"},
        {"source": "w_phi.success", "target": ".ok"},
        {"source": "field.node_positions", "target": "w_pos.data"},
        {"source": ".path_positions", "target": "w_pos.file_path"},
        {"source": "w_pos.success", "target": ".ok_positions"},
        {"source": "field.node_count", "target": ".node_count"},
    ]
    for s in ("", "_b"):
        n2, c2, p2 = mesh_reader(s)
        nodes += n2
        conns += c2
        ports += p2
    return dict(HEADER, compounds=[{"name": "User::PFFlip::collider_setup", "ports": ports,
                                    "compounds": [adapter("collider_field")],
                                    "compoundNodes": nodes, "connections": conns,
                                    "values": values}])


def library():
    return dict(HEADER, compounds=[adapter("User::PFFlip::collider_field")])


if __name__ == "__main__":
    for fname, g in (("collider_field.json", library()), ("collider_setup.json", setup())):
        with open(os.path.join(HERE, fname), "w", newline="\n") as f:
            json.dump(g, f, indent=1)
            f.write("\n")
        print("wrote", fname)
