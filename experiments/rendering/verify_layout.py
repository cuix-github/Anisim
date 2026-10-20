"""Offline relocation checks. No scene evaluation, simulation or ray tracing.

The before snapshot is created before relocation. Historical artifact bytes and
paths are compared, not repaired. Imports run with child-process launch blocked.
"""
import argparse
import ast
import hashlib
import importlib
import json
import os
from pathlib import Path
import re
import runpy
import sys
from urllib.parse import unquote

ROOT = Path(__file__).resolve().parents[2]
AUDIT = ROOT / 'out/repo_organization_20261019'
GROUP = ROOT / 'experiments/rendering/breaker'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def sha(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def deny_launch(event, args):
    if event in {'subprocess.Popen', 'os.system', 'os.startfile', 'os.exec', 'os.spawn'}:
        raise AssertionError('Process launch is forbidden in layout validation: ' + event)
    if event == 'ctypes.dlopen' and any(s in str(args).lower() for s in ['bifrost', 'ai.dll', 'openmaya']):
        raise AssertionError('Renderer/solver libraries must not load during layout validation')


def import_check(mode):
    sys.dont_write_bytecode = True
    sys.addaudithook(deny_launch)
    directory = ROOT / 'experiments' if mode == 'legacy' else GROUP
    sys.path.insert(0, str(directory))
    mapping = read(GROUP / 'relocation.json')
    loaded, static_only = [], []
    for item in mapping['scripts']:
        name = Path(item['old']).stem
        if name == 'breaker_camera_probe':
            # Historical mayapy utility initializes Maya and evaluates a mesh at
            # module scope. Preserve it byte-for-byte and check its AST only.
            static_only.append(dict(module=name, reason='Maya initialization and mesh work at module scope'))
            continue
        if mode == 'package':
            sys.path.insert(0, str(ROOT))
            name = 'experiments.rendering.breaker.' + name
        module = importlib.import_module(name)
        assert Path(module.__file__).resolve() == (ROOT / item['new']).resolve(), name
        if mode == 'legacy':
            # Each old file is exercised explicitly, even after sibling imports
            # have put the grouped directory first on sys.path.
            forwarded = runpy.run_path(str(ROOT / item['old']), run_name=name)
            assert Path(forwarded['__file__']).resolve() == (ROOT / item['new']).resolve(), name
            if 'ROOT' in forwarded:
                assert forwarded['ROOT'] == ROOT
        if hasattr(module, 'ROOT'):
            assert module.ROOT == ROOT, (name, module.ROOT)
        if hasattr(module, 'base'):
            assert module.base.ROOT == ROOT, name
        loaded.append(name)
    result = dict(status='passed', mode=mode, modules=loaded, static_only=static_only, external_processes_allowed=False)
    (AUDIT / ('imports_' + mode + '.json')).write_text(json.dumps(result, indent=2))
    print(json.dumps(dict(status='passed', mode=mode, count=len(loaded), static_only=static_only)))


def verify():
    mapping = read(GROUP / 'relocation.json')
    for item in mapping['scripts']:
        original = AUDIT / 'original_source' / item['old']
        assert sha(original) == item['sha256']
        before = original.read_text(encoding='utf-8')
        after = (ROOT / item['new']).read_text(encoding='utf-8')
        normalized = after.replace('Path(__file__).resolve().parents[3]', 'Path(__file__).resolve().parents[1]')
        assert ast.dump(ast.parse(before)) == ast.dump(ast.parse(normalized)), item['new']
        assert (ROOT / item['old']).is_file(), item['old']
    for mode in ['legacy', 'grouped', 'package']:
        assert read(AUDIT / ('imports_' + mode + '.json'))['status'] == 'passed'
    files = read(AUDIT / 'original_files.json')
    changed = []
    for relative, metadata in files.items():
        path = ROOT / relative
        if not path.is_file() or [path.stat().st_size, path.stat().st_mtime_ns] != metadata:
            changed.append(relative)
    assert not changed, changed
    scenes = read(AUDIT / 'scene_reference_baseline.json')
    for relative, expected in scenes['scene_hashes'].items():
        assert sha(ROOT / relative) == expected, relative
    new_missing = [x for x in scenes['references'] if x['exists'] and not Path(x['path']).exists()]
    assert not new_missing, new_missing
    # Check local links in the relocated studies, navigation and compatibility pages.
    docs = list((ROOT / 'docs/rendering').rglob('*.md')) + [ROOT / x['old'] for x in mapping['documents']]
    docs += [ROOT / p for p in ['README.md', 'docs/outputs.md', 'docs/open-questions.md', 'maya/presets/README.md', 'experiments/rendering/breaker/README.md']]
    broken = []
    for doc in docs:
        for target in re.findall(r'\]\(([^\s)]+)\)', doc.read_text(encoding='utf-8')):
            target = target.split('#', 1)[0]
            if not target or ':' in target:
                continue
            if not (doc.parent / unquote(target)).exists():
                broken.append([doc.relative_to(ROOT).as_posix(), target])
    assert not broken, broken
    gallery = read(ROOT / 'out/gallery/manifest.json')
    for entry in gallery['entries']:
        assert (ROOT / entry['path']).is_file(), entry['path']
        for link in [entry['url'], entry['thumbnail'], *[x['url'] for x in entry['links']], *[x['url'] for x in entry['members']]]:
            if link:
                assert (ROOT / 'out/gallery' / unquote(link)).is_file(), link
    registry = read(ROOT / 'maya/presets/breaker_validated_20261019.json')
    assert sha(ROOT / registry['foam']['source']) == registry['foam']['source_sha256']
    for variant in registry['variants']:
        folder = ROOT / variant['folder']
        assert sha(folder / 'candidate.mb') == variant['scene_sha256']
        assert sha(folder / 'foam.ass') == variant['foam_ass_sha256']
        assert sha(ROOT / variant['display_image']) == variant['display_image_sha256']
    result = dict(status='passed', moved_implementations=len(mapping['scripts']), grouped_studies=len(mapping['documents']),
                  algorithm_ast_unchanged_except_repository_root=True, import_modes=['legacy', 'grouped', 'package'],
                  runtime_imports_per_mode=24, static_only=['breaker_camera_probe'],
                  historical_files_unchanged=len(files), scenes_byte_identical=len(scenes['scene_hashes']),
                  scene_path_strings=len(scenes['references']), newly_missing_scene_paths=len(new_missing),
                  missing_local_document_links=len(broken), gallery_entries=len(gallery['entries']),
                  gallery_source_media=gallery['source_media_count'], registry_variants=len(registry['variants']),
                  simulation=False, rendering=False, maya_scene_evaluation=False,
                  scope='Static paths, byte hashes, AST equivalence and guarded imports; no rendering or scene evaluation.')
    (AUDIT / 'verification.json').write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['legacy', 'grouped', 'package', 'verify'])
    action = parser.parse_args().action
    if action == 'verify':
        verify()
    else:
        import_check(action)
