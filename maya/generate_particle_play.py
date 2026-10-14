"""Generates the 3D particle playback graphs, stock nodes only.

particle_play.json  the compound User::PFFlip::particle_play: Maya's frame picks one
    frame of a play cache (maya/play_cache.py: liquid.####.npy positions and
    speed.####.npy speeds), and out come points the Bifrost viewport draws, sized by
    point_size and coloured from slow_color to fast_color by speed / speed_scale
    (clamped at 1). Listed in pfflip_graphs_config.json.

particle_play_test.json  the bifcmd wrapper the gate drives: the same compound at a
    given frame, its points' positions, colours and sizes written back as .npy.

Regenerate after editing:  python generate_particle_play.py
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "graphs"))
from generate_sim_graph import P  # noqa: E402

HEADER = {"header": {"metadata": [{"metaName": "adskFileFormatVersion", "metaValue": "100L"}]},
          "namespaces": [], "types": []}
F3 = "Math::float3"


def f3(x, y, z):
    return {"x": f"{x}f", "y": f"{y}f", "z": f"{z}f"}


def multi(name, op):
    return {"nodeName": name, "nodeType": "Core::Math::" + op, "multiInPortNames": ["a", "b"]}


def play(name):
    ports = [P("liquid_pattern", "input", "string", ""),
             P("speed_pattern", "input", "string", ""),
             P("frame", "input", "float", "1f"),
             P("point_size", "input", "float", "0.3f"),
             P("speed_scale", "input", "float", "25f"),
             P("slow_color", "input", F3, f3(0.04, 0.22, 0.6)),
             P("fast_color", "input", F3, f3(0.95, 0.97, 1.0)),
             P("points", "output", "Object")]
    nodes = [
        {"nodeName": "t_pos", "valueType": "array<Math::float3>"},
        {"nodeName": "t_flt", "valueType": "array<float>"},
        {"nodeName": "f2i", "nodeType": "Core::Type_Conversion::to_int"},
        {"nodeName": "rpos", "nodeType": "File::NumPy::read_NumPy"},
        {"nodeName": "rspd", "nodeType": "File::NumPy::read_NumPy"},
        {"nodeName": "cpts", "nodeType": "Geometry::Points::construct_points"},
        # colour = slow + min(speed / scale, 1) * (fast - slow), per point
        multi("frac", "divide"),
        {"nodeName": "one", "valueType": "float"},
        multi("clamp", "min"),
        {"nodeName": "frac3", "nodeType": "Core::Conversion::scalar_to_vector3"},
        {"nodeName": "minus", "valueType": "float"},
        multi("neg", "multiply"),
        multi("delta", "add"),
        multi("ramp", "multiply"),
        multi("color", "add"),
        {"nodeName": "s_size", "nodeType": "Geometry::Properties::set_geo_property"},
        {"nodeName": "s_col", "nodeType": "Geometry::Properties::set_geo_property"},
    ]
    values = [
        {"valueName": "t_pos.value", "valueType": "array<Math::float3>", "value": [f3(0, 0, 0)]},
        {"valueName": "t_flt.value", "valueType": "array<float>", "value": ["0f"]},
        {"valueName": "one.value", "valueType": "float", "value": "1f"},
        {"valueName": "minus.value", "valueType": "float", "value": "-1f"},
        {"valueName": "s_size.property", "valueType": "string", "value": "point_size"},
        {"valueName": "s_col.property", "valueType": "string", "value": "point_color"},
    ]
    conns = [
        {"source": ".frame", "target": "f2i.from"},
        {"source": "t_pos.output", "target": "rpos.type"},
        {"source": ".liquid_pattern", "target": "rpos.file_path"},
        {"source": "f2i.int", "target": "rpos.frame"},
        {"source": "t_flt.output", "target": "rspd.type"},
        {"source": ".speed_pattern", "target": "rspd.file_path"},
        {"source": "f2i.int", "target": "rspd.frame"},
        {"source": "rpos.data", "target": "cpts.point_position"},
        {"source": "rspd.data", "target": "frac.first.a"},
        {"source": ".speed_scale", "target": "frac.first.b"},
        {"source": "frac.output", "target": "clamp.first.a"},
        {"source": "one.output", "target": "clamp.first.b"},
        {"source": "clamp.minimum", "target": "frac3.x"},
        {"source": "clamp.minimum", "target": "frac3.y"},
        {"source": "clamp.minimum", "target": "frac3.z"},
        {"source": ".slow_color", "target": "neg.first.a"},
        {"source": "minus.output", "target": "neg.first.b"},
        {"source": ".fast_color", "target": "delta.first.a"},
        {"source": "neg.output", "target": "delta.first.b"},
        {"source": "frac3.vector3", "target": "ramp.first.a"},
        {"source": "delta.output", "target": "ramp.first.b"},
        {"source": ".slow_color", "target": "color.first.a"},
        {"source": "ramp.output", "target": "color.first.b"},
        {"source": "cpts.points", "target": "s_size.geometry"},
        {"source": ".point_size", "target": "s_size.data"},
        {"source": ".point_size", "target": "s_size.default"},
        {"source": "s_size.out_geometry", "target": "s_col.geometry"},
        {"source": "color.output", "target": "s_col.data"},
        {"source": ".slow_color", "target": "s_col.default"},
        {"source": "s_col.out_geometry", "target": ".points"},
    ]
    return {"name": name, "ports": ports, "compoundNodes": nodes, "connections": conns,
            "values": values}


def test():
    """bifcmd: the compound at one frame, its points' properties written back."""
    ports = [P("liquid_pattern", "input", "string", ""),
             P("speed_pattern", "input", "string", ""),
             P("frame", "input", "float", "1f"),
             P("point_size", "input", "float", "0.3f"),
             P("speed_scale", "input", "float", "25f"),
             P("path_positions", "input", "string", ""),
             P("path_colors", "input", "string", ""),
             P("path_sizes", "input", "string", "")]
    nodes = [{"nodeName": "play", "nodeType": "particle_play"},
             {"nodeName": "t_f3", "valueType": "array<Math::float3>"},
             {"nodeName": "t_flt", "valueType": "array<float>"}]
    values = [{"valueName": "t_f3.value", "valueType": "array<Math::float3>", "value": [f3(0, 0, 0)]},
              {"valueName": "t_flt.value", "valueType": "array<float>", "value": ["0f"]}]
    conns = [{"source": "." + p, "target": "play." + p}
             for p in ("liquid_pattern", "speed_pattern", "frame", "point_size", "speed_scale")]
    for prop, short, typ in (("point_position", "positions", "t_f3"),
                             ("point_color", "colors", "t_f3"),
                             ("point_size", "sizes", "t_flt")):
        g, w = "g_" + short, "w_" + short
        nodes += [{"nodeName": g, "nodeType": "Geometry::Properties::get_geo_property"},
                  {"nodeName": w, "nodeType": "File::NumPy::write_NumPy"}]
        values += [{"valueName": g + ".property", "valueType": "string", "value": prop},
                   {"valueName": w + ".overwrite", "valueType": "bool", "value": "true"},
                   {"valueName": w + ".create_directories", "valueType": "bool", "value": "true"},
                   {"valueName": w + ".frame", "valueType": "long", "value": "0L"}]
        conns += [{"source": "play.points", "target": g + ".geometry"},
                  {"source": typ + ".output", "target": g + ".type"},
                  {"source": g + ".data", "target": w + ".data"},
                  {"source": ".path_" + short, "target": w + ".file_path"},
                  {"source": w + ".success", "target": ".ok_" + short}]
        ports.append(P("ok_" + short, "output", "bool"))
    return {"name": "User::PFFlip::particle_play_test", "ports": ports,
            "compounds": [play("particle_play")], "compoundNodes": nodes,
            "connections": conns, "values": values}


if __name__ == "__main__":
    for fname, c in (("particle_play.json", play("User::PFFlip::particle_play")),
                     ("particle_play_test.json", test())):
        with open(os.path.join(HERE, fname), "w", newline="\n") as f:
            json.dump(dict(HEADER, compounds=[c]), f, indent=1)
            f.write("\n")
        print("wrote", fname)
