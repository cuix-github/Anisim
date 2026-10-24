"""Matched dark-background movies from retained caches; never simulates."""
import argparse
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image,ImageDraw
from inspect_spray_events_2d import panel,COLOURS,SMALL
from inspect_settling_phase_2d import composite
from current_plunge_2d import final_frame
from settling_phase_2d import frame,measure,write,DT,EVERY,CHUNK_STEPS


def load(p):return json.loads(p.read_text())


def samples(root,limit=None):
    progress=load(root/'progress.json')
    yield 0.,dict(np.load(root/'initial.npz'))
    for chunk in progress['chunks']:
        for j in range(0,CHUNK_STEPS,EVERY):
            t=(chunk['offset']+j+1)*DT
            if chunk['offset']==0 and j==0:continue
            if limit is not None and t>=limit:break
            yield t,frame(root/chunk['directory'],j)


def pair(left,right,t,labels):
    imgs=[panel(a,label,t,scale=10) for a,label in zip((left,right),labels)]
    im=Image.new('RGB',(1288,966),(15,20,27))
    im.paste(imgs[0],(0,0));im.paste(imgs[1],(648,0))
    d=ImageDraw.Draw(im);x=12
    for label,rgb in COLOURS:
        d.rectangle((x,877,x+12,889),fill=rgb)
        d.text((x+18,874),label,font=SMALL,fill='white');x+=250
    d.text((12,900),'Same initial water, view, markers and playback. Cumulative changes; NOT a causal ablation.',font=SMALL,fill='white')
    d.text((12,925),'Current: phase 2 / escape 2 / transport 3 / detached primary / finite density / evaporation OFF',font=SMALL,fill=(180,195,210))
    return im


def movie(path,images,fps):
    assert not path.exists(),path
    proc=None;n=0
    try:
        for im in images:
            if proc is None:
                proc=subprocess.Popen(['ffmpeg','-hide_banner','-loglevel','error',
                    '-f','rawvideo','-pixel_format','rgb24','-video_size',f'{im.width}x{im.height}',
                    '-framerate',str(fps),'-i','-','-an','-c:v','libx264','-crf','18',
                    '-pix_fmt','yuv420p','-movflags','+faststart',str(path)],stdin=subprocess.PIPE)
            proc.stdin.write(im.tobytes());n+=1
    finally:
        if proc:
            proc.stdin.close()
            if proc.wait()!=0:raise RuntimeError('ffmpeg failed')
    return n


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);a=ap.parse_args()
    root=a.run.resolve();job=load(root/'job.json');progress=load(root/'progress.json')
    assert progress['status'] in ('near_rest','limit_reached_not_settled')
    dest=root/'visuals';dest.mkdir(exist_ok=False)
    previous=Path(job['old_long_run']);original=Path(job['original_20s'])
    rows=load(root/'measurements.json');oldrows=load(previous/'measurements.json')
    fps=1/(DT*EVERY)
    def full(speed=1):
        for t,a in samples(root):
            im=composite(a,t,measure(a,t));d=ImageDraw.Draw(im)
            d.rectangle((0,0,im.width,31),fill=(0,0,0))
            d.text((12,5),f'CURRENT C++ | {speed}x playback | simulation t = {t:.2f} s',font=SMALL,fill='white')
            d.rectangle((0,938,im.width,966),fill=(15,20,27))
            d.text((12,944),'Current cumulative candidate | evaporation OFF | fixed dt/damping | diagnostic, not certification',font=SMALL,fill='white')
            if speed==1 and any(abs(t-v)<.01 for v in (0,3.2,6,20,40,80,120,180)):
                im.save(dest/f'current_{t:07.2f}s.png')
            yield im
    n=movie(dest/'current_full_realtime.mp4',full(),fps)
    movie(dest/'current_full_4x.mp4',full(4),fps*4)
    def earliest():
        for i,(t,a) in enumerate(samples(root,20)):
            b=dict(np.load(root/'initial.npz')) if i==0 else frame(original/'bifrost',i*EVERY)
            im=pair(b,a,t,('Original C++ (Oct 23)','Current C++ candidate'))
            if any(abs(t-v)<.01 for v in (3.2,6,14,19.8)):
                im.save(dest/f'original_comparison_{t:05.2f}s.png')
            yield im
    movie(dest/'original_vs_current_20s.mp4',earliest(),fps)
    def old_long():
        for (t,a),(u,b) in zip(samples(root),samples(previous)):
            assert abs(t-u)<1e-9
            yield pair(b,a,t,('Previous phase/escape repair','Current cumulative candidate'))
    movie(dest/'previous_vs_current_long.mp4',old_long(),fps)
    last=root/progress['chunks'][-1]['directory'];end=progress['simulated_seconds']
    af=final_frame(last);lastrow=measure(af,end)
    composite(af,end,lastrow).save(dest/'current_exact_final.png')
    op=load(previous/'progress.json');index=len(progress['chunks'])-1
    bf=final_frame(previous/op['chunks'][index]['directory'])
    pair(bf,af,end,('Previous phase/escape repair','Current cumulative candidate')).save(dest/'final_comparison.png')
    fig,axs=plt.subplots(2,2,figsize=(13,8),layout='constrained')
    for series,label in [(oldrows,'Previous phase/escape repair'),(rows,'Current candidate')]:
        ts=[r['time'] for r in series];peak=max(r['water_kinetic'] for r in series)
        for ax,key in [(axs[0,0],'water_rms'),(axs[0,1],'water_speed_p99')]:
            ax.semilogy(ts,[max(r[key],1e-8) for r in series],label=label)
        axs[1,0].plot(ts,[r['spray'] for r in series],label=label)
        axs[1,1].semilogy(ts,[max(r['water_kinetic']/peak,1e-9) for r in series],label=label)
    axs[0,0].axhline(.05,color='grey',ls='--');axs[0,1].axhline(.2,color='grey',ls='--')
    axs[1,1].axhline(1e-4,color='grey',ls='--')
    for ax,title in zip(axs.flat,['Primary water RMS speed (cells/s)','Primary water speed P99 (cells/s)',
        'Living auxiliary spray (evaporation differs)','Kinetic energy / each run sampled peak']):
        ax.set_title(title);ax.set_xlabel('Simulation time (s)');ax.grid(alpha=.2);ax.legend()
    fig.suptitle('Same plunge, cumulative implementation changes; no retuned damping')
    fig.savefig(dest/'decay_comparison.png',dpi=140);plt.close(fig)
    write(dest/'job.json',dict(kind='Read-only matched cache visualization',source=str(root),
        frames=n,fps=fps,markers='Same original colors and pixel markers; NOT droplet radii',
        original=str(original),previous=str(previous),limits=job['caveats'],
        final_checkpoint_displayed=True,no_holds_or_interpolation=True))
    print(str(dest),flush=True)


if __name__=='__main__':main()
