"""Read-only physical/visual comparison. Identical display markers in every arm."""
import argparse
import json
from pathlib import Path
import subprocess

import numpy as np
from PIL import Image, ImageDraw
from inspect_spray_events_2d import panel, COLOURS, SMALL
from spray_events_2d import load_frame, NX, NY


def final(root, arm):
    d=root/arm
    if arm=='python':
        with np.load(d/'final.npz') as z: return {k:z[k] for k in z.files}
    return dict(pos=np.load(d/'final.npy')[:,:2],vel=np.load(d/'final_velocities.npy')[:,:2],
        typ=np.load(d/'final_phase.npy')>.5,spr=np.load(d/'final_spray.npy')>.5,
        esc=np.zeros(len(np.load(d/'final.npy')),bool),rad=np.load(d/'final_radius.npy'))


def measure(a):
    primary=~a['spr']; water=primary & (a['typ']==1)
    p=a['pos'][primary].astype(float);v=a['vel'][primary].astype(float)
    m=np.where(a['typ'][primary]==1,1000.,1.)
    return dict(water_count=int(water.sum()),primary_count=int(primary.sum()),spray=int(a['spr'].sum()),
        energy=float(np.sum(m*(.5*np.sum(v*v,axis=1)+float(np.float32(9.8))*p[:,1]))),
        water_com=a['pos'][water].astype(float).mean(0).tolist(),
        height_p99=float(np.quantile(a['pos'][water,1],.99)),
        water_rms=float(np.sqrt(np.mean(np.sum(a['vel'][water].astype(float)**2,axis=1)))))


def compare(pair):
    waters=[a['pos'][(a['typ']==1)&~a['spr']].astype(float) for a in pair]
    h=[np.histogram2d(p[:,0],p[:,1],bins=(16,20),range=((0,NX),(0,NY)))[0] for p in waters]
    return dict(tv=float(.5*abs(h[0]/h[0].sum()-h[1]/h[1].sum()).sum()),
        com=float(np.linalg.norm(waters[0].mean(0)-waters[1].mean(0))),
        height_p99=float(abs(np.quantile(waters[0][:,1],.99)-np.quantile(waters[1][:,1],.99))))


def analyse(root):
    job=json.loads((root/'job.json').read_text());rows=[]
    for j in range(0,job['steps'],job['dump_every']):
        a=[load_frame(root,arm,j) for arm in ('python','bifrost')]
        rows.append(dict(time=(j+1)*job['dt'],difference=compare(a),
                         **{arm:measure(z) for arm,z in zip(('python','bifrost'),a)}))
    a=[final(root,arm) for arm in ('python','bifrost')]
    rows.append(dict(time=job['steps']*job['dt'],difference=compare(a),
                     **{arm:measure(z) for arm,z in zip(('python','bifrost'),a)}))
    worst={k:max(r['difference'][k] for r in rows) for k in ('tv','com','height_p99')}
    initial=measure(dict(np.load(root/'initial.npz')))
    physics={arm:dict(max_energy_ratio=max(r[arm]['energy'] for r in rows)/initial['energy'],
        final_energy_ratio=rows[-1][arm]['energy']/initial['energy'],final=rows[-1][arm]) for arm in ('python','bifrost')}
    result=dict(job=job,max_difference=worst,
        similarity_pass=worst['tv']<.15 and worst['com']<1 and worst['height_p99']<4,
        physics=physics,rows=rows)
    (root/'comparison.json').write_text(json.dumps(result,indent=2))
    return result


def movie(old,new,dest,kind='phase'):
    job=json.loads((new/'job.json').read_text());oldjob=json.loads((old/'job.json').read_text())
    assert all(job[k]==oldjob[k] for k in ('nx','ny','steps','dt','dump_every'))
    dest.mkdir(parents=True,exist_ok=False)
    frames=dest/'frames';frames.mkdir();cppframes=dest/'cpp_frames';cppframes.mkdir()
    images={}
    for i,j in enumerate(range(0,job['steps'],job['dump_every'])):
        t=(j+1)*job['dt']
        labels=['OLD NumPy: mass-only','OLD C++: tag mean','NEW NumPy: Eq.7 + pure faces','NEW C++: Eq.7 + pure faces']
        if kind=='escape':
            labels=['Count criterion | NumPy','Count criterion | C++',
                    'Phase criterion | NumPy','Phase criterion | C++']
        arrays=[load_frame(r,a,j) for r in (old,new) for a in ('python','bifrost')]
        tiles=[panel(a,l,t) for a,l in zip(arrays,labels)]
        w,h=tiles[0].size
        im=Image.new('RGB',(2*w+8,2*h+80),(15,20,27))
        for k,tile in enumerate(tiles):im.paste(tile,((k%2)*(w+8),(k//2)*(h+4)))
        c=Image.new('RGB',(2*w+8,h+76),(15,20,27));c.paste(tiles[1],(0,0));c.paste(tiles[3],(w+8,0))
        for target,y in ((im,2*h+12),(c,h+10)):
            d=ImageDraw.Draw(target);x=12
            for label,rgb in COLOURS:
                d.rectangle((x,y,x+12,y+12),fill=rgb);d.text((x+18,y-4),label,font=SMALL,fill='white');x+=220
            d.text((12,y+24),'Identical seed, view, colors, marker sizes and playback rate | 1 simulated second / second',font=SMALL,fill='white')
            d.text((12,y+44),'Display markers are not physical droplet sizes. Independent pressure solves.',font=SMALL,fill=(180,190,205))
        im.save(frames/f'frame.{i:04d}.png');c.save(cppframes/f'frame.{i:04d}.png')
        if any(abs(t-x)<job['dt']*job['dump_every']*.51 for x in (3.2,6,12,18,19.9)):
            im.save(dest/f'compare_{t:.2f}s.png');c.save(dest/f'cpp_{t:.2f}s.png');images[str(t)]=str(dest/f'compare_{t:.2f}s.png')
    for f,name in ((frames,'old_new_four_way.mp4'),(cppframes,'cpp_before_after.mp4')):
        subprocess.run(['ffmpeg','-hide_banner','-loglevel','error','-framerate',str(1/(job['dt']*job['dump_every'])),
            '-i',str(f/'frame.%04d.png'),'-c:v','libx264','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart',str(dest/name)],check=True)
    (dest/'job.json').write_text(json.dumps(dict(kind='Read-only video',old=str(old),new=str(new),
        comparison_kind=kind,marker_scale=8,markers='water/air 1px; escaped/spray 2px',images=images),indent=2))


def main():
    p=argparse.ArgumentParser();p.add_argument('--run',required=True,type=Path)
    p.add_argument('--old',type=Path);p.add_argument('--visuals',type=Path)
    p.add_argument('--comparison-kind',choices=('phase','escape'),default='phase');a=p.parse_args()
    r=analyse(a.run.resolve());print(json.dumps({k:v for k,v in r.items() if k not in ('rows','job')},indent=2))
    if a.old:movie(a.old.resolve(),a.run.resolve(),a.visuals.resolve(),a.comparison_kind)
    if not r['similarity_pass']:raise SystemExit('Independent similarity gate failed; outputs retained')


if __name__=='__main__':main()
