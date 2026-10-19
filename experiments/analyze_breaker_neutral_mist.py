"""Offline evidence for the two frozen-cache neutral-water / finer-mist previews."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'out/runs/breaker_neutral_mist_20261019_v2'
PRIOR = ROOT / 'out/runs/breaker_optical_aovs_20261019_v3'
OIIO = 'C:/Program Files/Autodesk/Arnold/Maya2027/bin/oiiotool.exe'
MAYAPY = 'C:/Program Files/Autodesk/Maya2027/bin/mayapy.exe'
COMPONENTS = ['specular','transmission','diffuse','sss','coat','sheen','emission','volume','background']
LAYERS = ['beauty','beauty_raw',*COMPONENTS,'N','P']


def pfm(path):
    with path.open('rb') as f:
        assert f.readline().strip() == b'PF'
        w,h = map(int,f.readline().split())
        scale = float(f.readline())
        a = np.frombuffer(f.read(),dtype='<f4' if scale<0 else '>f4').reshape(h,w,3)
    assert np.isfinite(a).all()
    return np.flipud(a).copy()


def convert(folder):
    import PyOpenColorIO as ocio
    config = ocio.Config.CreateFromFile('C:/Program Files/Autodesk/Maya2027/resources/OCIO-configs/Maya2022-default/config.ocio')
    processor = config.getProcessor(ocio.DisplayViewTransform(src='ACEScg',display='sRGB',view='ACES 1.0 SDR-video')).getDefaultCPUProcessor()
    for layer in ['beauty','beauty_raw','specular','transmission']:
        a = np.ascontiguousarray(np.flipud(pfm(folder/(layer+'.pfm'))))
        h,w,_ = a.shape
        processor.apply(ocio.PackedImageDesc(a,w,h,3))
        target = folder/(layer+'_display.pfm')
        with target.open('wb') as f:
            f.write(f'PF\n{w} {h}\n-1.0\n'.encode())
            f.write(a.astype('<f4').tobytes())
        subprocess.run([OIIO,str(target),'--attrib','oiio:ColorSpace','sRGB','-d','uint8','-o',str(folder/(layer+'_view.png'))],check=True)


def sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f,'sha256').hexdigest()


def main():
    from PIL import Image, ImageDraw, ImageFont
    result = dict(status='running', render_count=2, simulation_rerun=False,
                  water_surface_rebuilt=False, mist_classification_and_volume_built=True, arms={})
    states = {}
    for arm in ['original','side']:
        folder=OUT/arm
        receipt=json.loads((folder/'receipt.json').read_text())
        if receipt['status']!='passed':
            recovery=json.loads((folder/'recovery.json').read_text())
            assert recovery['status']=='passed'
            assert recovery['receipt_sha256']==sha(folder/'receipt.json')
            receipt={**receipt,**recovery}
        assert receipt['status']=='passed'
        assert receipt['exit_code']==0 or receipt.get('kind')=='output_collection_recovery_no_render'
        assert sha(folder/'beauty.png')==receipt['png_sha256']
        assert sha(folder/'optics.0001.exr')==receipt['exr_sha256']
        state=states[arm]=json.loads((folder/'state.json').read_text())
        assert state['defaultArnoldDriver']['colorManagement']==1
        assert state['defaultColorMgtGlobals']['wsn']=='ACEScg'
        assert state['defaultColorMgtGlobals']['vtn']=='ACES 1.0 SDR-video (sRGB)'
        info=subprocess.check_output([OIIO,'--info','-v',str(folder/'optics.0001.exr')],text=True)
        (folder/'exr_info.txt').write_text(info)
        assert '960 x  540, 54 channel, float openexr' in info
        assert 'arnold/license_state: "valid"' in info
        log=(folder/'arnold.1.log').read_text()
        assert 'ERROR' not in log and '| render done' in log
        data={}
        for layer in LAYERS:
            ch='R,G,B' if layer=='beauty' else ','.join(layer+'.'+s for s in ('XYZ' if layer in ['N','P'] else 'RGB'))
            target=folder/(layer+'.pfm')
            subprocess.run([OIIO,str(folder/'optics.0001.exr'),'--ch',ch,'-o',str(target)],check=True)
            data[layer]=pfm(target)
        subprocess.run([MAYAPY,str(Path(__file__).resolve()),'--convert',str(folder)],check=True)
        rawsum=sum(data[k].astype(np.float64) for k in COMPONENTS)
        closure=float(np.abs(rawsum-data['beauty_raw']).sum()/np.abs(data['beauty_raw']).sum())
        assert closure<1e-5, closure
        png=np.asarray(Image.open(folder/'beauty.png').convert('RGB'),dtype=float)
        view=np.asarray(Image.open(folder/'beauty_view.png').convert('RGB'),dtype=float)
        diff=np.abs(png-view)
        display=dict(mean_abs_levels=float(diff.mean()),fraction_over_8=float((diff.max(axis=2)>8).mean()),max_abs_levels=float(diff.max()))
        assert display['mean_abs_levels']<=1 and display['fraction_over_8']<=.01, display
        n=data['N']; mask=np.linalg.norm(n,axis=2)>.5
        Image.fromarray(np.uint8(np.clip(np.where(mask[...,None],(n+1)*.5,0),0,1)*255)).save(folder/'normals_view.png')
        record=dict(process_exit_code=receipt['exit_code'],collection_recovered=receipt.get('kind')=='output_collection_recovery_no_render',render_source_sha=receipt['git_sha'],render_call_s=receipt['worker']['render_call_s'],process_s=receipt['process_s'],
                    aov_closure_relative_l1=closure,png_matches_ocio_view=display,
                    component_raw_rgb_share={k:float(data[k].sum()/data['beauty_raw'].sum()) for k in COMPONENTS})
        if arm=='original':
            record['prior_N_exact']=bool(np.array_equal(data['N'],pfm(PRIOR/'B/N.pfm')))
            record['prior_P_exact']=bool(np.array_equal(data['P'],pfm(PRIOR/'B/P.pfm')))
            normal_difference=np.abs(data['N']-pfm(PRIOR/'B/N.pfm'))
            record['prior_N_max_abs_difference']=float(normal_difference.max())
            record['prior_N_mean_abs_difference']=float(normal_difference.mean())
            record['mean_linear_specular_rgb']=data['specular'].mean((0,1)).tolist()
            record['mean_linear_transmission_rgb']=data['transmission'].mean((0,1)).tolist()
            record['baseline_mean_linear_specular_rgb']=pfm(PRIOR/'B/specular.pfm').mean((0,1)).tolist()
            record['baseline_mean_linear_transmission_rgb']=pfm(PRIOR/'B/transmission.pfm').mean((0,1)).tolist()
        result['arms'][arm]=record
    assert states['original']['geometry']==states['side']['geometry']
    assert states['original']['mist_visible']
    assert states['original']['mist_ports']['spray_detail']==.25
    assert states['original']['mist_ports']['cloud_radius']==.5
    assert result['arms']['original']['component_raw_rgb_share']['volume'] > 0
    differing=[key for key in states['original'] if states['original'][key]!=states['side'][key]]
    assert sorted(differing)==['persp','perspShape'], differing
    camera_rotation=states['side']['persp']['rotate'][0]
    assert camera_rotation[0]==0 and camera_rotation[2]==0
    assert states['side']['perspShape']['orthographic']
    theta=np.radians(12); yaw=np.radians(camera_rotation[1])
    perpendicular_dot=float(-np.sin(yaw)*np.cos(theta)-np.cos(yaw)*np.sin(theta))
    assert abs(perpendicular_dot)<1e-12
    result.update(status='passed',between_images_changed_nodes=differing,
                  geometry=states['original']['geometry'],material=states['original']['water'],
                  sky=states['original']['candidateSky'],mist_ports=states['original']['mist_ports'],mist_shader=states['original']['candidateMistFog'],
                  side_view_dot_dominant_wave=perpendicular_dot,
                  render_source_sha=json.loads((OUT/'original/receipt.json').read_text())['git_sha'])
    items=[(PRIOR/'B/beauty_view.png','Corrected baseline / existing camera'),
           (OUT/'original/beauty_view.png','Neutral absorption + finer cached mist / same camera'),
           (OUT/'side/beauty_view.png','Same candidate / close side view (15 m wide)')]
    font=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',20)
    sheet=Image.new('RGB',(960,3*574),(20,23,28));draw=ImageDraw.Draw(sheet)
    for i,(path,label) in enumerate(items):
        draw.text((12,i*574+5),label,font=font,fill='white')
        sheet.paste(Image.open(path).convert('RGB'),(0,i*574+34))
    sheet.save(OUT/'comparison.png')
    panels=Image.new('RGB',(1440,2*439),(20,23,28));draw=ImageDraw.Draw(panels)
    for i,name in enumerate(['beauty','beauty_raw','specular','transmission']):
        x=(i%2)*720;y=(i//2)*439
        draw.text((x+12,y+5),name,font=font,fill='white')
        panels.paste(Image.open(OUT/'original'/(name+'_view.png')).convert('RGB').resize((720,405)),(x,y+34))
    panels.save(OUT/'aov_panels.png')
    (OUT/'analysis.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    if len(sys.argv)>1 and sys.argv[1]=='--convert':
        convert(Path(sys.argv[2]))
    else:
        main()
