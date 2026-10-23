"""Read-only analysis of original full-spray causal/precision controls."""
import argparse
import json
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from inspect_spray_events_2d import panel


def load(path, cpp, step=None):
    if not cpp:
        with np.load(path / ('final.npz' if step is None else f'frame.{step:04d}.npz')) as z:
            return {k: z[k] for k in z.files}
    prefixes = dict(pos='frame', vel='vel', typ='ph', spr='spr', esc='esc', rad='rad')
    final = dict(pos='final', vel='final_velocities', typ='final_phase',
                 spr='final_spray', rad='final_radius')
    a = {k: np.load(path / (final[k]+'.npy' if step is None else f'{v}.{step:04d}.npy'))
         for k, v in prefixes.items() if step is not None or k != 'esc'}
    a['pos'], a['vel'] = a['pos'][:, :2], a['vel'][:, :2]
    for k in ('spr', 'esc'):
        if k in a:
            a[k] = a[k] > .5
    return a


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--root', type=Path, required=True)
    p.add_argument('--original', type=Path, required=True)
    a = p.parse_args(); root = a.root.resolve(); dest = root / 'analysis'
    dest.mkdir(exist_ok=False)
    arms = {
        'Original NumPy': (a.original/'python', False),
        'Original C++': (root/'trace/bifrost', True),
        'NumPy, shared pressure': (root/'shared_pressure', False),
        'NumPy, matched independent': (root/'independent_matched', False),
        'C++, tighter pressure': (root/'cpp_strict/bifrost', True),
    }
    all_results = dict(scope='Original 20s, escape/emission/flight ON; phase 0 tag law in matched arms',
        timing='Some causal runs overlap. Recorded wall times are not performance benchmarks.', arms={})
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    colours = ['#888888', '#c43c39', '#2b73bc', '#2d995d', '#a57423']
    for (label, (path, cpp)), color in zip(arms.items(), colours):
        ledger = json.loads((path/'ledger.json').read_text())
        chunks = []
        for lo in range(0, 8000, 800):
            rr = ledger[lo:lo+800]
            chunks.append(dict(start_s=lo/400, end_s=(lo+800)/400,
                **{k: sum(r[k] for r in rr) for k in ('accepted', 'evaporated', 'reentered')},
                mean_candidates=float(np.mean([r['candidates'] for r in rr])), alive=rr[-1]['alive']))
        fin = load(path, cpp)
        water = (fin['typ'] == 1) & ~fin['spr']; air = (fin['typ'] == 0) & ~fin['spr']
        result = dict(chunks=chunks, final_count=len(fin['pos']), final_spray=int(fin['spr'].sum()),
            total_births=sum(r['accepted'] for r in ledger),
            water_rms=float(np.sqrt(np.mean(np.sum(fin['vel'][water].astype(float)**2, axis=1)))),
            air_rms=float(np.sqrt(np.mean(np.sum(fin['vel'][air].astype(float)**2, axis=1)))))
        all_results['arms'][label] = result
        stride = 40
        axs[0].plot((np.arange(8000)[::stride]+1)/400, [r['alive'] for r in ledger[::stride]],
                    label=label, color=color, linestyle='--' if 'shared' in label else '-')
        axs[1].plot([r['end_s'] for r in chunks], [r['accepted']/2 for r in chunks],
                    label=label, color=color, linestyle='--' if 'shared' in label else '-')
    shared = load(root/'shared_pressure', False)
    cpp_final = load(root/'trace/bifrost', True)
    all_results['shared_final_bit_identical'] = {
        k: bool(np.array_equal(shared[k], cpp_final[k])) for k in cpp_final}
    for name in ('shared_pressure', 'independent_matched'):
        r = json.loads((root/name/'receipt.json').read_text())
        all_results[name] = {k: v for k, v in r.items() if k != 'rows'}
        events = json.loads((root/name/'source_events.json').read_text())
        groups = []
        for lo in range(0, 8000, 800):
            ev = [e for e in events if lo <= e['step'] < lo+800]
            weights = np.array([e['events'] for e in ev]); n = int(weights.sum())
            side = np.array([min(e['x'], 64-e['x']) < 1.5 for e in ev])
            groups.append(dict(start_s=lo/400, end_s=(lo+800)/400, events=n,
                side_wall_events=int(weights[side].sum()) if n else 0,
                mean_source_gradient=float(np.average([e['gradient'] for e in ev], weights=weights)) if n else 0))
        all_results[name]['source_event_groups'] = groups
    axs[0].set(ylabel='Alive auxiliary spray particles', title='Original 20-second scene: causal controls')
    axs[0].legend(fontsize=8); axs[0].grid(alpha=.2)
    axs[1].set(xlabel='Simulation time (s)', ylabel='Births / second (2 s windows)')
    axs[1].grid(alpha=.2); fig.tight_layout(); fig.savefig(dest/'spray_counts.png', dpi=150); plt.close(fig)
    # Actual cached 18.0025-second states, with the same marker rules as the movie.
    from PIL import Image, ImageDraw
    titles = ['Original C++', 'NumPy, shared pressure', 'NumPy, matched independent', 'C++, tighter pressure']
    tiles = [panel(load(*arms[t], step=7200), t, 18.0025, scale=7) for t in titles]
    sheet = Image.new('RGB', (2*tiles[0].width, 2*tiles[0].height), (15,20,27))
    for i, tile in enumerate(tiles):
        sheet.paste(tile, ((i%2)*tile.width, (i//2)*tile.height))
    sheet.save(dest/'at_18s.png')
    (dest/'summary.json').write_text(json.dumps(all_results, indent=2))
    print(json.dumps(all_results, indent=2))


if __name__ == '__main__':
    main()
