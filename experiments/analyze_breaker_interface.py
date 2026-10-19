"""Validate the two SSS images without any rendering or simulation."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
import breaker_deep_water as base
from analyze_breaker_deep_water import pfm, convert, OIIO, MAYAPY, COMPONENTS, LAYERS
OUT=base.ROOT/'out/runs/breaker_interface_sss_20261019_v3'


def check(folder):
    from PIL import Image
    receipt=json.loads((folder/'receipt.json').read_text());assert receipt['status']=='passed'
    state=json.loads((folder/'state.json').read_text())
    assert state['mist'] is False
    assert state['interfaceEnable']['input2']==int(folder.name=='B')
    assert base.sha(folder/'beauty.png')==receipt['png_sha256']
    exr=folder/receipt['worker']['exr'];assert base.sha(exr)==receipt['exr_sha256']
    info=subprocess.check_output([OIIO,'--info','-v',str(exr)],text=True)
    (folder/'exr_info.txt').write_text(info)
    assert '960 x  540, 54 channel, float openexr' in info
    assert 'arnold/license_state: "valid"' in info
    log=(folder/'arnold.1.log').read_text();assert 'ERROR' not in log and '| render done' in log
    data={}
    for layer in LAYERS:
        ch='R,G,B' if layer=='beauty' else ','.join(layer+'.'+s for s in ('XYZ' if layer in ['N','P'] else 'RGB'))
        target=folder/(layer+'.pfm')
        subprocess.run([OIIO,str(exr),'--ch',ch,'-o',str(target)],check=True)
        data[layer]=pfm(target)
    subprocess.run([MAYAPY,str(Path(__file__).resolve()),'--convert',str(folder)],check=True)
    total=data['beauty_raw'].sum(dtype=float)
    closure=float(np.abs(sum(data[k].astype(float) for k in COMPONENTS)-data['beauty_raw']).sum()/total)
    assert closure<1e-5,closure
    png=np.asarray(Image.open(folder/'beauty.png').convert('RGB'),dtype=float)
    view=np.asarray(Image.open(folder/'beauty_view.png').convert('RGB'),dtype=float)
    delta=np.abs(png-view);assert delta.mean()<=1 and (delta.max(2)>8).mean()<.01
    shares={k:float(data[k].sum(dtype=float)/total) for k in COMPONENTS}
    assert shares['volume']==0,'Mist/volume must be absent'
    if folder.name=='A':assert np.max(np.abs(data['sss']))<1e-8
    else:assert shares['sss']>1e-5,'SSS mask did not reach rendered shading'
    result=dict(status='passed',aov_closure_relative_l1=closure,component_raw_rgb_share=shares,
                display_mean_abs_levels=float(delta.mean()),render_call_s=receipt['worker']['render_call_s'],
                process_s=receipt['process_s'],renderer_license='valid')
    base.save(folder/('control_check.json' if folder.name=='A' else 'image_check.json'),result)
    print(json.dumps(result,indent=2));return result


def report(out):
    from PIL import Image,ImageDraw,ImageFont
    states=[json.loads((out/arm/'state.json').read_text()) for arm in ['A','B']]
    differences=[k for k in states[0] if states[0][k]!=states[1][k]]
    assert differences==['interfaceEnable'],differences
    assert states[0]['geometry']==json.loads((base.PRIOR/'B/state.json').read_text())['geometry']
    normals=[pfm(out/arm/'N.pfm') for arm in ['A','B']]
    positions=[pfm(out/arm/'P.pfm') for arm in ['A','B']]
    checks={arm:json.loads((out/arm/('control_check.json' if arm=='A' else 'image_check.json')).read_text()) for arm in ['A','B']}
    result=dict(status='passed',simulation_rerun=False,mesh_rebuilt=False,mist=False,render_count=2,
                changed_nodes=differences,N_exact=bool(np.array_equal(*normals)),P_exact=bool(np.array_equal(*positions)),
                N_max_abs_difference=float(np.abs(normals[0]-normals[1]).max()),arms=checks,
                approximation=states[0]['approximation'],mask_sha256=states[0]['mask_sha256'])
    font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',20)
    sheet=Image.new('RGB',(960,1148),(20,23,28));draw=ImageDraw.Draw(sheet)
    for i,(arm,label) in enumerate([('A','A | Corrected baseline / close side view / mist off'),
                                    ('B','B | Same scene + curvature-gated randomwalk SSS')]):
        draw.text((12,i*574+6),label,font=font,fill='white')
        sheet.paste(Image.open(out/arm/'beauty_view.png').convert('RGB'),(0,i*574+34))
    sheet.save(out/'comparison.png');base.save(out/'analysis.json',result);print(json.dumps(result,indent=2))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--convert':convert(Path(sys.argv[2]))
    else:
        p=argparse.ArgumentParser();p.add_argument('mode',choices=['A','B','report']);p.add_argument('--out',type=Path,default=OUT)
        a=p.parse_args()
        if a.mode=='report':report(a.out)
        else:check(a.out/a.mode)
