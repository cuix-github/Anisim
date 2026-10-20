"""Index existing local media. Never invoke Maya, Arnold, a solver or a renderer.

Only out/gallery is written. Images are resized for browsing; source artifacts
and experimental receipts stay at their original paths. Works offline/file://.
"""
import argparse
from collections import defaultdict
import hashlib
import json
import os
from pathlib import Path
import re
from urllib.parse import quote

from PIL import Image, ImageOps

ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / 'out/gallery'
REGISTRY = ROOT / 'maya/presets/breaker_validated_20261019.json'
IMAGES = {'.png', '.jpg', '.jpeg', '.webp', '.gif'}
VIDEOS = {'.mp4', '.mov', '.webm', '.avi'}


def read(path):
    try:
        return json.loads(path.read_text(encoding='utf-8'))
    except (OSError, ValueError):
        return {}


def url(path):
    return quote(os.path.relpath(path, OUTPUT).replace('\\', '/'), safe='/')


def discover():
    media = []
    for top in [ROOT / 'out', ROOT / 'results']:
        for parent, dirs, files in os.walk(top):
            dirs[:] = [d for d in dirs if d not in {'gallery', 'repo_organization_20261019', '__pycache__', 'fake_out'}
                       and not (Path(parent) / d).is_symlink()]
            for name in files:
                path = Path(parent) / name
                if path.suffix.lower() in IMAGES | VIDEOS:
                    media.append(path)
    return sorted(media)


def study(path):
    parts = path.relative_to(ROOT).parts
    if parts[:2] == ('out', 'runs'):
        return parts[2]
    if parts[0] == 'results':
        return '/'.join(parts[:3])
    return '/'.join(parts[:2])


def study_doc(name):
    if name.startswith('wave2d_'):
        return ROOT / 'docs/wave-feature-diagnostics.md'
    mapping = [('mix60_perspective', 'breaker-perspective-study.md'),
               ('six_views', 'breaker-six-view-study.md'),
               ('mix60', 'breaker-six-view-study.md'),
               ('foam_mix', 'breaker-foam-view-study.md'),
               ('foam_views', 'breaker-foam-view-study.md'),
               ('foam_layer', 'breaker-foam-layer-study.md'),
               ('interface', 'breaker-interface-study.md'),
               ('deep_water', 'breaker-deep-water-study.md'),
               ('neutral_mist', 'breaker-neutral-mist-study.md'),
               ('correspondence', 'breaker-surface-correspondence.md'),
               ('render_study', 'breaker-render-study.md')]
    for token, document in mapping:
        if token in name:
            return ROOT / 'docs/rendering/studies' / document
    return ROOT / 'docs/evidence.md'


def build():
    OUTPUT.mkdir(parents=True, exist_ok=True)
    thumbs = OUTPUT / 'thumbnails'
    thumbs.mkdir(exist_ok=True)
    registry = read(REGISTRY)
    chosen = {v['display_image']: v for v in registry['variants']}
    for path in sorted((ROOT / 'out/runs/breaker_six_views_20261019/shots').glob('*/C/beauty_view.png')):
        rel = path.relative_to(ROOT).as_posix()
        chosen.setdefault(rel, dict(id=path.parents[1].name, uniform_weight=.2, clustered_weight=.8,
                                   assessment='Six-view observation; no final look selected.'))
    media = discover()
    directories = defaultdict(list)
    for path in media:
        if path.suffix.lower() in IMAGES:
            directories[path.parent].append(path)
    # Collapse numbered sequences into one browsable entry, retaining every
    # member in the manifest and detail view rather than listing thousands of cards.
    sequence_members = {}
    numbered = re.compile(r'(?:fr|frame|img|image|render)?[_.-]?\d{3,}$', re.I)
    for directory, paths in directories.items():
        members = [p for p in paths if numbered.fullmatch(p.stem)]
        if len(members) >= 8:
            sequence_members[members[0]] = members
    consumed = {p for members in sequence_members.values() for p in members[1:]}
    entries, errors = [], []
    for path in media:
        rel = path.relative_to(ROOT).as_posix()
        if path in consumed:
            continue
        # The verified display conversion is the browse image. The original PNG
        # is retained as a link; it is never deleted or silently recolored.
        if path.name == 'beauty.png' and path.with_name('beauty_view.png').exists():
            continue
        selected = chosen.get(rel)
        members = sequence_members.get(path, [])
        kind = ('video' if path.suffix.lower() in VIDEOS else 'sequence' if members
                else 'comparison' if 'comparison' in path.stem else
                'render' if path.stem in {'beauty', 'beauty_view'} else 'diagnostic')
        state = read(path.parent / 'state.json')
        check = read(path.parent / 'image_check.json') if kind == 'render' else {}
        validation = 'passed' if check.get('status') == 'passed' else 'failed' if check.get('status') == 'failed' else 'unreviewed'
        camera = state.get('view', {})
        if not isinstance(camera, dict):
            camera = {}
        ratio = None
        if selected:
            ratio = [selected['uniform_weight'], selected['clustered_weight']]
        elif state.get('alpha_B') is not None:
            ratio = [1 - state['alpha_B'], state['alpha_B']]
        elif state.get('foam_arm') in {'A', 'B'}:
            ratio = [1, 0] if state['foam_arm'] == 'A' else [0, 1]
        identity = hashlib.sha256(rel.encode()).hexdigest()[:20]
        thumb = None
        if kind != 'video':
            signature = hashlib.sha256((rel + str(path.stat().st_mtime_ns)).encode()).hexdigest()[:20]
            target = thumbs / (signature + '.jpg')
            try:
                if not target.exists():
                    with Image.open(path) as source:
                        preview = ImageOps.exif_transpose(source).convert('RGB')
                        preview.thumbnail((640, 360), Image.Resampling.LANCZOS)
                        preview.save(target, quality=88)
                thumb = url(target)
            except (OSError, ValueError) as exc:
                errors.append(dict(path=rel, error=str(exc)))
                validation = 'unreadable'
        links = [{'label': '研究记录', 'url': url(study_doc(study(path)))}]
        for name, label in [('state.json', '场景参数'), ('image_check.json', '图像检查'),
                            ('receipt.json', '运行记录'), ('candidate.mb', 'Maya 场景'),
                            ('optical.mb', 'Maya 场景'), ('foam.ass', '泡沫 ASS'),
                            ('beauty.png', '原始 PNG'), ('optics.0001.exr', 'EXR'),
                            ('scenes/wave_2d.ma', 'Maya 粒子场景'),
                            ('receipts/particle_audit.json', '粒子检查')]:
            if (path.parent / name).is_file():
                links.append(dict(label=label, url=url(path.parent / name)))
        entries.append(dict(id=identity, path=rel, url=url(path), thumbnail=thumb,
                            title=selected['id'] if selected else (path.parent.name + ' · ' + path.stem),
                            study=study(path), kind=kind, validation=validation, selected=bool(selected),
                            time_s=camera.get('time_s'),
                            width_m=camera.get('width_m', camera.get('target_plane_width_m')), ratio=ratio,
                            note=selected.get('assessment', '') if selected else '', links=links,
                            members=[dict(name=p.name, url=url(p)) for p in members]))
    entries.sort(key=lambda e: (not e['selected'], e['study'], e['path']))
    manifest = dict(schema=1, source_media_count=len(media), entry_count=len(entries),
                    selected_count=sum(e['selected'] for e in entries), entries=entries, thumbnail_errors=errors,
                    registry=url(REGISTRY), setup=url(ROOT / 'maya/presets/README.md'),
                    guide=url(ROOT / 'docs/rendering/README.md'),
                    policy='Index only. Existing media, scenes, simulation data and receipts remain in place.')
    (OUTPUT / 'manifest.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2), encoding='utf-8')
    template = (Path(__file__).parent / 'gallery.html').read_text(encoding='utf-8')
    payload = json.dumps(manifest, ensure_ascii=False).replace('<', '\\u003c')
    (OUTPUT / 'index.html').write_text(template.replace('/*__MANIFEST__*/', payload), encoding='utf-8')
    print(json.dumps({k: v for k, v in manifest.items() if k != 'entries'}, ensure_ascii=False, indent=2))
    return manifest


if __name__ == '__main__':
    argparse.ArgumentParser(description=__doc__).parse_args()
    build()
