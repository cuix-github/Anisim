"""Audit the saved Maya scenes and receipts, without evaluating or rendering them."""
import argparse
import hashlib
import json
from pathlib import Path
import re


def digest(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def scene(path):
    text = path.read_text()
    nodes, graphs = {}, {}
    for match in re.finditer(r'^createNode.*?(?=^createNode|^select|^connectAttr|\Z)', text, re.S | re.M):
        block = match.group()
        name = re.search(r'-n "([^"]+)"', block).group(1)
        block = re.sub(r'\s*rename -uid "[^"]+";', '', block)
        serial = re.search(r'setAttr "\.sc" -type "string" \((.*?)\);', block, re.S)
        if serial:
            graph = json.loads(''.join(json.loads(s) for s in re.findall(r'"(?:[^"\\]|\\.)*"', serial.group(1))))
            graphs[name] = {v['valueName']: v['value'] for v in graph['compounds'][0]['values']}
        nodes[name] = block.strip()
    connections = sorted(re.sub(r'\s+', ' ', m.group()) for m in re.finditer(r'^connectAttr.*?;', text, re.S | re.M))
    return nodes, graphs, connections


def compare(left, right):
    a, ga, ca = scene(left)
    b, gb, cb = scene(right)
    return {
        'changed_nodes': sorted(n for n in a.keys() | b.keys() if a.get(n) != b.get(n)),
        'removed_connections': sorted(set(ca) - set(cb)),
        'added_connections': sorted(set(cb) - set(ca)),
        'changed_graph_values': {n: {p: [ga[n].get(p), gb[n].get(p)] for p in ga[n].keys() | gb[n].keys()
                                     if ga[n].get(p) != gb[n].get(p)} for n in ga}
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=Path(__file__).resolve().parents[1] / 'out/runs/breaker_render_study_20261019')
    args = parser.parse_args()
    root = args.out.resolve()
    expected = {'A_to_B': ['foamCurvature'], 'B_to_C': ['particleMistShape'],
                'C_to_D': ['particleMistShape', 'particleSurfaceShape'], 'C_to_E': ['sky'],
                'D_to_candidate_D': ['persp', 'perspShape', 'sky']}
    report = {}
    for pair, nodes in expected.items():
        left, right = pair.split('_to_')
        diff = compare(root / left / 'breaker.ma', root / right / 'breaker.ma')
        assert diff['changed_nodes'] == nodes, (pair, diff)
        assert diff['removed_connections'] == (['connectAttr "waterFoamMix.out" "waterSG.ss";'] if pair == 'A_to_B' else []), (pair, diff)
        assert diff['added_connections'] == (['connectAttr "water.out" "waterSG.ss";'] if pair == 'A_to_B' else []), (pair, diff)
        values = diff['changed_graph_values']
        if pair == 'C_to_D':
            assert values == {
                'particleSurfaceShape': {'particle_surface.mcf_iterations': ['30L', '10L']},
                'particleMistShape': {'particle_mist.mcf_iterations': ['30L', '10L']}}, values
        else:
            assert all(not v for v in values.values()), (pair, values)
        report[pair] = diff
    original = root.parent / 'wavesea4/renders/breaking/tube.ma'
    diff = compare(original, root / 'A/breaker.ma')
    assert diff['changed_nodes'] == ['persp', 'perspShape'], diff
    assert not diff['added_connections'] and not diff['removed_connections'], diff
    assert all(not v for v in diff['changed_graph_values'].values()), diff
    report['original_to_A'] = diff
    receipts = [json.loads((root / arm / 'receipt.json').read_text()) for arm in ['A', 'B', 'C', 'D', 'E', 'candidate_D']]
    cache = root.parent / 'wavesea4/cache/breaking_section'
    assert {p.name: digest(p) for p in sorted(cache.iterdir()) if p.is_file()} == receipts[0]['input_sha256'], 'Source cache changed'
    for receipt in receipts:
        assert receipt['status'] == 'passed' and receipt['exit_code'] == 0, receipt['arm']
        assert receipt['input_sha256'] == receipts[0]['input_sha256'], receipt['arm']
        assert digest(root / receipt['arm'] / 'breaker.0001.png') == receipt['image_sha256'], receipt['arm']
        assert digest(root / receipt['arm'] / 'breaker.ma') == receipt['scene_sha256'], receipt['arm']
        log = (root / receipt['arm'] / 'breaker.1.log').read_text()
        knobs = receipt['overrides']
        assert f"rendering image at {knobs['SURF_WIDTH']} x {knobs['SURF_HEIGHT']}, 4 AA samples" in log, receipt['arm']
        triangles = int(re.search(r'(\d+) tris -- input0_0_0_graph_outputs_0_mesh', log).group(1))
        assert triangles == (292516 if receipt['arm'] in ['D', 'candidate_D'] else 274090), (receipt['arm'], triangles)
        receipt['water_triangles'] = triangles
    result = {'status': 'passed', 'scene_differences': report, 'renders': receipts}
    (root / 'verified_study.json').write_text(json.dumps(result, indent=2) + '\n')
    print('PASS: six render receipts, source cache hashes, baked scene changes, resolutions, AA samples and water triangle counts')


if __name__ == '__main__':
    main()
