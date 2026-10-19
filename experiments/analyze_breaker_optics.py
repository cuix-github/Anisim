"""Verify and inspect saved optical A/B outputs without launching Maya or rendering."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import subprocess

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
OIIO = Path('C:/Program Files/Autodesk/Arnold/Maya2027/bin/oiiotool.exe')
OCIO = Path('C:/Program Files/Autodesk/Maya2027/resources/OCIO-configs/Maya2022-default/config.ocio')
CURL = (200, 175, 800, 530)
LAYERS = ['beauty', 'beauty_raw', 'specular', 'specular_direct', 'specular_indirect',
          'transmission', 'transmission_direct', 'transmission_indirect', 'diffuse',
          'sss', 'coat', 'sheen', 'emission', 'volume', 'background', 'N', 'P']
COMPONENTS = ['specular', 'transmission', 'diffuse', 'sss', 'coat', 'sheen',
              'emission', 'volume', 'background']


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def run(*args):
    return subprocess.check_output([str(OIIO), *map(str, args)], text=True)


def channels(layer):
    if layer == 'beauty':
        return 'R,G,B'
    return ','.join(layer + '.' + suffix for suffix in ('XYZ' if layer in ('N', 'P') else 'RGB'))


def pfm(path):
    with Path(path).open('rb') as stream:
        assert stream.readline().strip() == b'PF'
        width, height = map(int, stream.readline().split())
        scale = float(stream.readline())
        data = np.frombuffer(stream.read(), dtype='<f4' if scale < 0 else '>f4')
    assert data.size == width * height * 3
    data = np.flipud(data.reshape(height, width, 3)).copy()
    assert np.isfinite(data).all()
    return data


def read_layers(folder):
    exr = folder / 'optics.0001.exr'
    arrays = {}
    for layer in LAYERS:
        target = folder / (layer + '.pfm')
        run(exr, '--ch', channels(layer), '-o', target)
        arrays[layer] = pfm(target)
    # Display every radiometric panel through the same OCIO view as the Maya PNG.
    for layer in ['beauty', 'beauty_raw', 'specular', 'transmission', 'specular_indirect']:
        run('--colorconfig', OCIO, exr, '--ch', channels(layer), '--chnames', 'R,G,B',
            '--ociodisplay:from=ACEScg', 'sRGB', 'ACES 1.0 SDR-video', '-d', 'uint8',
            '-o', folder / (layer + '_view.png'))
    n = arrays['N']
    visible = np.linalg.norm(n, axis=2) > .5
    mapped = np.where(visible[..., None], (n + 1) * .5, 0)
    Image.fromarray(np.uint8(np.clip(mapped, 0, 1) * 255)).save(folder / 'normals_view.png')
    return arrays


def state_difference(a, b, path=''):
    changes = {}
    for key in a.keys() | b.keys():
        label = path + key
        if isinstance(a.get(key), dict) and isinstance(b.get(key), dict):
            changes.update(state_difference(a[key], b[key], label + '.'))
        elif a.get(key) != b.get(key):
            changes[label] = [a.get(key), b.get(key)]
    return changes


def region(array):
    x0, y0, x1, y1 = CURL
    return array[y0:y1, x0:x1]


def panel_image(items, destination, columns=2, width=720):
    height = width * 540 // 960
    rows = (len(items) + columns - 1) // columns
    sheet = Image.new('RGB', (columns * width, rows * (height + 34)), (20, 23, 28))
    draw = ImageDraw.Draw(sheet)
    font = ImageFont.truetype('C:/Windows/Fonts/arial.ttf', 19)
    for i, (file, label) in enumerate(items):
        x, y = (i % columns) * width, (i // columns) * (height + 34)
        draw.text((x + 12, y + 7), label, fill='white', font=font)
        picture = Image.open(file).convert('RGB').resize((width, height), Image.Resampling.LANCZOS)
        sheet.paste(picture, (x, y + 34))
    sheet.save(destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=ROOT / 'out/runs/breaker_optical_aovs_20261019_v3')
    args = parser.parse_args()
    root = args.out.resolve()
    assert not subprocess.check_output(['git', 'status', '--porcelain'], cwd=ROOT).strip()
    receipts = {arm: json.loads((root / arm / 'receipt.json').read_text()) for arm in 'AB'}
    states = {arm: json.loads((root / arm / 'state.json').read_text()) for arm in 'AB'}
    assert receipts['A']['git_sha'] == receipts['B']['git_sha']
    assert receipts['A']['baked_sha256'] == receipts['B']['baked_sha256']
    cache = ROOT / 'out/runs/wavesea4/cache/breaking_section'
    assert {p.name: digest(p) for p in cache.iterdir() if p.is_file()} == receipts['A']['cache_sha256'] == receipts['B']['cache_sha256']
    assert receipts['A']['control_check']['status'] == 'passed'
    changes = state_difference(states['A'], states['B'])
    # Public Maya node docs define lightBlindData as an addr, not a light knob:
    # https://help.autodesk.com/cloudhelp/2024/ENU/Maya-Tech-Docs/Nodes/defaultLightList.html
    runtime_addresses = {key: changes.pop(key) for key in list(changes) if key.endswith('.lightBlindData')}
    assert changes == {'defaultArnoldRenderOptions.GISpecularDepth': [1, 3]}, changes
    arrays, checks = {}, {}
    for arm in 'AB':
        folder, receipt = root / arm, receipts[arm]
        assert receipt['status'] == 'passed'
        assert digest(folder / 'beauty.png') == receipt['png_sha256']
        assert digest(folder / 'optics.0001.exr') == receipt['exr_sha256']
        if receipt['exit_code'] != 0:
            # This run's collector failed after successful image writes. Preserve
            # that failure and verify the specific recovery, never ignore exits.
            original = json.loads((folder / 'failed_collection_receipt.json').read_text())
            worker = json.loads((folder / 'worker_result.json').read_text())
            assert original['status'] == 'failed' and original['exit_code'] == 1
            assert receipt['artifact_validation'] == 'passed_after_offline_collection_recovery'
            assert 'assert len(images) == 1' in worker['error'] and 'AssertionError: []' in worker['error']
        info = run('--info', '-v', folder / 'optics.0001.exr')
        assert '960 x  540, 54 channel, float openexr' in info
        assert f'arnold/specular_depth: {1 if arm == "A" else 3}' in info
        assert 'arnold/transmission_depth: 8' in info and 'arnold/license_state: "valid"' in info
        log = (folder / 'arnold.1.log').read_text()
        assert 'ERROR' not in log and '| render done' in log
        data = arrays[arm] = read_layers(folder)
        total = sum(data[layer].astype(np.float64) for layer in COMPONENTS)
        error = np.abs(total - data['beauty_raw'])
        # Relative normalization uses the total absolute raw RGB signal, not black pixels.
        relative = float(error.sum() / np.abs(data['beauty_raw']).sum())
        assert relative < 1e-5, relative
        assert np.max(np.abs(data['specular'] - data['specular_direct'] - data['specular_indirect'])) < 1e-5
        assert np.max(np.abs(data['transmission'] - data['transmission_direct'] - data['transmission_indirect'])) < 1e-5
        png = np.asarray(Image.open(folder / 'beauty.png').convert('RGB'), dtype=float)
        converted = np.asarray(Image.open(folder / 'beauty_view.png').convert('RGB'), dtype=float)
        display_error = float(np.abs(png - converted).mean())
        assert display_error < 1, 'EXR display transform disagrees with Maya PNG'
        signal = float(region(data['beauty_raw']).sum())
        shares = {name: float(region(data[name]).sum() / signal) for name in COMPONENTS}
        checks[arm] = {'exr_channels': 54, 'aov_sum_relative_l1_error': relative,
                       'aov_sum_max_absolute_error': float(error.max()),
                       'exr_display_vs_png_mae_levels': display_error,
                       'curl_component_linear_rgb_sum_fractions': shares,
                       'render_call_s': receipt['worker']['render_call_s'],
                       'process_wall_s': receipt['process_wall_s'],
                       'actual_process_exit_code': receipt['exit_code'],
                       'arnold_render_s': float(re.search(r'arnold/stats/time/render: ([\d.]+)', info).group(1))}
    assert np.array_equal(arrays['A']['N'], arrays['B']['N'])
    assert np.array_equal(arrays['A']['P'], arrays['B']['P'])
    delta = arrays['B']['beauty_raw'] - arrays['A']['beauty_raw']
    reference = region(arrays['A']['beauty_raw'])
    a_png = np.asarray(Image.open(root / 'A/beauty.png').convert('RGB'), dtype=float)
    b_png = np.asarray(Image.open(root / 'B/beauty.png').convert('RGB'), dtype=float)
    display_delta = np.abs(b_png - a_png)
    metrics = {'curl_raw_rgb_sum_change_fraction': float(region(delta).sum() / reference.sum()),
               'curl_raw_relative_l1_change': float(np.abs(region(delta)).sum() / reference.sum()),
               'curl_display_mae_levels': float(region(display_delta).mean()),
               'curl_display_fraction_pixels_changed_over_8_levels': float((region(display_delta).max(axis=2) > 8).mean()),
               'full_display_mae_levels': float(display_delta.mean()),
               'normal_aov_exact': True, 'position_aov_exact': True}
    for layer in ['specular', 'specular_direct', 'specular_indirect', 'transmission']:
        values = region(arrays['B'][layer] - arrays['A'][layer])
        metrics[layer + '_delta_over_control_raw_rgb'] = float(values.sum() / reference.sum())
    result = {'status': 'passed', 'render_git_sha': receipts['A']['git_sha'],
              'analysis_git_sha': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
              'state_differences': changes, 'excluded_process_addresses': runtime_addresses,
              'control': receipts['A']['control_check'],
              'curl_xyxy': CURL, 'checks': checks, 'comparison': metrics,
              'input_sha256': receipts['A']['cache_sha256'], 'baked_sha256': receipts['A']['baked_sha256'],
              'limits': ['The rectangular curl ROI is not an object-isolated water mask.',
                         'RGB sums are scene-linear channel sums, not perceptual luminance or material reflectance.',
                         'Raw components compose before denoising; OIDN-denoised components need not compose.',
                         'The two original Maya process exits were 1 due solely to the PNG collection path; retained failures and offline validation document recovery.',
                         'AOVs identify the first scattering category; transmission includes any later reflection events.',
                         'The comparison tests specular depth 1 versus 3 only, not convergence at unlimited depth.']}
    (root / 'analysis.json').write_text(json.dumps(result, indent=2) + '\n')
    panel_image([(root / 'A/beauty.png', 'A | specular depth 1'),
                 (root / 'B/beauty.png', 'B | specular depth 3')], root / 'comparison.png')
    panel_image([(root / 'B/beauty.png', 'Beauty | original display transform'),
                 (root / 'B/beauty_raw_view.png', 'Raw beauty | same display transform'),
                 (root / 'B/specular_view.png', 'Specular | same exposure'),
                 (root / 'B/transmission_view.png', 'Transmission | same exposure'),
                 (root / 'B/normals_view.png', 'World normals | N mapped from [-1,1] to [0,1]'),
                 (root / 'B/specular_indirect_view.png', 'Indirect specular | same exposure')], root / 'aov_panels.png')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
