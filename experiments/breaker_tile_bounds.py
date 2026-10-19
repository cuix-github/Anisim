"""Inspect public stock-volume bounds without density blur, MCF or rendering."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import time

import breaker_surface_correspondence as surface
import breaker_six_views as study
import breaker_deep_water as base


def inspect(tile):
    src = study.OUT/'frames/pre/tiles'/tile
    out = src/'bounds_audit'
    assert not out.exists()
    out.mkdir()
    g = surface.graph(1)
    surf = g['compounds'][0]['compounds'][0]
    needed = {'k_bnd'}
    while True:
        old = set(needed)
        for connection in surf['connections']:
            if connection['target'].split('.')[0] in needed and not connection['source'].startswith('.'):
                needed.add(connection['source'].split('.')[0])
        if needed == old:
            break
    surf['compoundNodes'] = [n for n in surf['compoundNodes'] if n['nodeName'] in needed]
    surf['values'] = [v for v in surf['values'] if v['valueName'].split('.')[0] in needed]
    surf['connections'] = [c for c in surf['connections'] if c['target'].split('.')[0] in needed]
    surf['ports'] = [p for p in surf['ports'] if p['portDirection'] == 'input']
    for key in ['min_bound', 'max_bound']:
        surf['ports'].append(surface.gen.P(key, 'output', 'Math::float3'))
        surf['connections'].append(dict(source='k_bnd.'+key, target='.'+key))
    inputs = [p for p in surf['ports'] if p['portDirection'] == 'input']
    wrapper = dict(name='User::PFFlip::raw_bounds', ports=inputs+[surface.gen.P(k,'output','Math::float3') for k in ['min_bound','max_bound']],
                   compounds=[surf], compoundNodes=[dict(nodeName='surf',nodeType='particle_surface')],
                   connections=[dict(source='.'+p['portName'],target='surf.'+p['portName']) for p in inputs]
                   +[dict(source='surf.'+k,target='.'+k) for k in ['min_bound','max_bound']])
    base.save(out/'graph.json',dict(surface.gen.HEADER, compounds=[wrapper]))
    env = os.environ.copy()
    env['PATH'] = str(surface.BIF/'bin')+';'+str(surface.BIF/'thirdparty/bin')+';'+env['PATH']
    command = [str(surface.BIF/'bin/bifcmd.exe'),str(out/'graph.json')]
    values = dict(surface.KNOBS,liquid_pattern=(src/'cache/liquid.####').as_posix(),
                  radius_pattern=(src/'cache/radius.####').as_posix(),frame=1)
    for key,value in values.items():command += ['--set-port',key,str(value)]
    start = time.perf_counter()
    with (out/'bifcmd.log').open('w') as log:
        process=subprocess.Popen(command,env=env,stdout=log,stderr=subprocess.STDOUT)
        while process.poll() is None:
            if study.free_gb()<5 or time.perf_counter()-start>90:
                process.kill();process.wait();raise RuntimeError('Bounds-only guard')
            time.sleep(1)
    assert process.returncode==0
    text=(out/'bifcmd.log').read_text()
    lines=[line for line in text.splitlines() if 'Output port' in line]
    base.save(out/'receipt.json',dict(status='passed',wall_s=time.perf_counter()-start,output=lines,
              retained_nodes=sorted(needed),simulation=False,render=False))
    print(json.dumps(study.read(out/'receipt.json'),indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('tile')
    inspect(parser.parse_args().tile)
