"""Read-only particle movies, decay curves and spray lifecycle accounting."""
import argparse
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw

from inspect_spray_events_2d import panel, COLOURS, SMALL
from settling_phase_2d import frame, measure, rest_check, DT, EVERY, CHUNK_STEPS, write


def records(root):
    job = json.loads((root/'job.json').read_text())
    progress = json.loads((root/'progress.json').read_text())
    rows = json.loads((root/'measurements.json').read_text())
    ledger = []
    previous = 0
    for chunk in progress['chunks']:
        ev = json.loads((root/chunk['directory']/'ledger.json').read_text())
        for e in ev:
            assert e['alive'] == previous+e['accepted']-e['reentered']-e['evaporated'], e
            previous = e['alive']
        ledger.extend(ev)
    assert len(ledger) == len(progress['chunks'])*CHUNK_STEPS
    last = root/progress['chunks'][-1]['directory']
    phases = np.load(last/'final_phase.npy') > .5
    spray = np.load(last/'final_spray.npy') > .5
    p = np.load(last/'final.npy')[:, :2]
    v = np.load(last/'final_velocities.npy')[:, :2].astype(float)
    water = phases & ~spray
    radius = np.load(last/'final_radius.npy')
    # In this enabled droplet model a surviving escaped water particle has a
    # positive radius; rejoining resets it to zero. Air escape flags are unused
    # by the water metrics. Include the exact checkpoint, not only the last dump.
    exact_final = measure(dict(pos=p, vel=v, typ=phases, spr=spray,
        esc=phases & (radius > 0), rad=radius), progress['simulated_seconds'])
    final_rest = rest_check(rows+[exact_final])
    peaks = {k:max(rows, key=lambda r:r[k]) for k in ('water_rms', 'water_kinetic', 'primary_energy')}
    summary = dict(job=job, status=progress['status'], simulation_seconds=progress['simulated_seconds'],
        elapsed_seconds=progress['elapsed_seconds'], near_rest=progress['near_rest'],
        near_rest_including_final_checkpoint=final_rest,
        max_pressure_residual=max(c['max_pressure_residual'] for c in progress['chunks']),
        max_particle_cfl=max(c['max_particle_cfl'] for c in progress['chunks']),
        peak_water_rms=dict(time=peaks['water_rms']['time'],value=peaks['water_rms']['water_rms']),
        peak_kinetic=peaks['water_kinetic']['water_kinetic'],
        max_primary_energy_ratio=peaks['primary_energy']['primary_energy']/rows[0]['primary_energy'],
        final=dict(water_count=int(water.sum()), primary_count=int((~spray).sum()),
             spray=int(spray.sum()), water_rms=float(np.sqrt(np.mean(np.sum(v[water]**2,axis=1)))),
             water_com=p[water].astype(float).mean(0).tolist(),
             height_p99=float(np.quantile(p[water,1], .99))),
        lifecycle={k:sum(e[k] for e in ledger) for k in
            ('events', 'requested', 'accepted', 'rejected_cap', 'rejected_geometry', 'evaporated', 'reentered')},
        last_auxiliary_birth_time=max(((i+1)*DT for i,e in enumerate(ledger) if e['accepted']), default=None),
        last_auxiliary_alive_time=max(((i+1)*DT for i,e in enumerate(ledger) if e['alive']), default=None),
        lifecycle_closes_every_step=True, last_sample=rows[-1],
        caveats=['Near-rest criteria are observational, not a correctness theorem.',
                 'Surface proxy excludes escaped water, but speed and kinetic criteria include it.',
                 'Particle mechanical energy excludes lost primary mass and zero-mass auxiliary spray.',
                 'Air, water and spray use visibility markers, not physical particle radii.'])
    bins = []
    for start in range(0,len(ledger),400):
        ev = ledger[start:start+400]
        bins.append(dict(time=(start+len(ev))*DT, births=sum(e['accepted'] for e in ev),
            reentered=sum(e['reentered'] for e in ev), evaporated=sum(e['evaporated'] for e in ev),
            alive=ev[-1]['alive'], candidates_mean=float(np.mean([e['candidates'] for e in ev]))))
    summary['spray_per_second'] = bins
    write(root/'summary.json', summary)
    return rows, summary, progress


def charts(dest, rows, summary):
    t=np.array([r['time'] for r in rows]);peak=summary['peak_kinetic']
    fig,axs=plt.subplots(2,2,figsize=(13.5,8.4),layout='constrained')
    ax=axs[0,0]
    for key,label,style in [('water_rms','All primary water RMS','-'),('bulk_water_rms','Non-escaped water RMS','--')]:
        ax.semilogy(t,[r[key] if r[key] and r[key]>0 else np.nan for r in rows],style,label=label)
    ax.axhline(.05,color='grey',ls='--',label='Near-rest RMS criterion')
    ax.set_ylim(.005,20)
    ax.set_ylabel('Speed (cells/s)');ax.legend();ax.set_title('Water motion; includes the difficult particles')
    ax=axs[0,1]
    ax.semilogy(t,[r['water_kinetic']/peak if r['water_kinetic']>0 else np.nan for r in rows],label='Water kinetic / sampled peak')
    ax.set_ylim(1e-6,2)
    ax.axhline(1e-4,color='grey',ls='--');ax.set_ylabel('Fraction of peak');ax.legend()
    ax.set_title('Motion energy (not total-energy conservation)')
    ax=axs[1,0];b=summary['spray_per_second'];bt=[r['time'] for r in b]
    for k,label in [('births','Births / 1 s'),('reentered','Returned / 1 s'),('evaporated','Evaporated / 1 s')]:
        ax.plot(bt,[r[k] for r in b],label=label,alpha=.85)
    ax.legend();ax.set_ylabel('Auxiliary particle count');ax.set_title('White-point population is births minus removals')
    ax=axs[1,1]
    ax.plot(t,[r['height_p99'] for r in rows],label='All-water height P99')
    ax.plot(t,[r['water_com'][1] for r in rows],label='All-water centre of mass')
    ax.set_ylabel('Height (cells)');ax.legend();ax.set_title('Height and centre of mass (not a volume measurement)')
    for ax in axs.flat: ax.set_xlabel('Simulation time (s)');ax.grid(alpha=.22)
    fig.suptitle('Current C++ plunge: fixed equations, dt and damping throughout',fontsize=16)
    fig.savefig(dest/'decay_curves.png',dpi=150);plt.close(fig)


def composite(a, t, row):
    left=panel(a,'C++ | all particles',t,scale=10)
    keep=a['typ'] & ~a['spr']
    water={k:v[keep] for k,v in a.items()}
    right=panel(water,'Same cache | primary water only',t,scale=10)
    im=Image.new('RGB',(left.width+right.width+8,left.height+100),(15,20,27))
    im.paste(left,(0,0));im.paste(right,(left.width+8,0))
    draw=ImageDraw.Draw(im);x=12;y=left.height+10
    for label,rgb in COLOURS:
        draw.rectangle((x,y,x+12,y+12),fill=rgb);draw.text((x+18,y-3),label,font=SMALL,fill='white');x+=250
    draw.text((12,y+24),f"Water RMS {row['water_rms']:.4f} cells/s | bulk {row['bulk_water_rms']:.4f} | "
              f"water {row['water_count']} | same view throughout",font=SMALL,fill='white')
    draw.text((12,y+46),'Markers are NOT droplet radii. Right pane hides air/spray only; no simulation data removed.',font=SMALL,fill=(180,195,210))
    draw.text((12,y+68),'Phase 2 / escape 2 | fixed dt 0.0025 s | no additional damping | time shown is simulation time',font=SMALL,fill=(180,195,210))
    return im


def movies(root,dest,rows,progress):
    dest.mkdir(exist_ok=True)
    samples=[]
    for chunk in progress['chunks']:
        samples.extend((root/chunk['directory'],j,(chunk['offset']+j+1)*DT)
                       for j in range(0,CHUNK_STEPS,EVERY))
    # Replace the first post-step sample with the exact initial condition. The
    # remaining samples are dt after nominal 0.2s intervals; duration stays T.
    # No repeated hold that could imply settling.
    samples[0]=(None,None,0.)
    process=None;selected=[];targets=[0.,3.2,20.,40.,80.,120.,180.,239.8]
    output=dest/'cpp_settling_realtime.mp4'
    assert not output.exists(), output
    try:
        for idx,(d,j,t) in enumerate(samples):
            a=dict(np.load(root/'initial.npz')) if d is None else frame(d,j)
            im=composite(a,t,rows[0 if idx==0 else idx+1])
            if process is None:
                process=subprocess.Popen(['ffmpeg','-hide_banner','-loglevel','error',
                    '-f','rawvideo','-pixel_format','rgb24','-video_size',f'{im.width}x{im.height}',
                    '-framerate',str(1/(DT*EVERY)),'-i','-','-an','-c:v','libx264','-crf','18',
                    '-pix_fmt','yuv420p','-movflags','+faststart',str(output)],stdin=subprocess.PIPE)
            process.stdin.write(im.tobytes())
            if any(abs(t-v)<.01 for v in targets) or idx==len(samples)-1:
                im.save(dest/f'state_{t:07.2f}s.png');selected.append((t,im.copy()))
    finally:
        if process:
            process.stdin.close()
            if process.wait()!=0: raise RuntimeError('ffmpeg failed')
    # Explicitly labeled shorter overview; the original is retained above.
    subprocess.run(['ffmpeg','-hide_banner','-loglevel','error','-i',str(output),
        '-vf',"setpts=PTS/4,drawbox=x=0:y=0:w=iw:h=32:color=black:t=fill,drawtext=fontfile='C\\:/Windows/Fonts/arial.ttf':text='4x SPEED overview - see time labels':fontcolor=white:fontsize=22:x=10:y=4",
        '-r','20','-frames:v',str(len(samples)),'-c:v','libx264','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart',str(dest/'cpp_settling_4x.mp4')],check=True)
    chosen=selected[-4:]
    sheet=Image.new('RGB',(1288,round(chosen[0][1].height/2)*len(chosen)),(15,20,27))
    for i,(_,im) in enumerate(chosen):
        thumb=im.resize((644,round(im.height/2)))
        # one row per time; retain a readable native-resolution last image too
        sheet.paste(thumb,(0,i*thumb.height))
    # Do not add a redundant wide empty half: crop to the actual panel width.
    sheet.crop((0,0,644,sheet.height)).save(dest/'late_contact_sheet.png')
    write(dest/'job.json',dict(kind='Read-only particle visualization',source=str(root),
        realtime=str(output),frame_count=len(samples),fps=1/(DT*EVERY),
        playback='1x and explicitly labeled 4x',markers='Same markers in every frame',
        comparison='Same cache; right hides air and auxiliary spray for inspection'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--run',required=True,type=Path)
    parser.add_argument('--video',action='store_true');a=parser.parse_args()
    root=a.run.resolve();dest=root/'visuals';dest.mkdir(exist_ok=True)
    rows,summary,progress=records(root);charts(dest,rows,summary)
    if a.video: movies(root,dest,rows,progress)
    print(json.dumps({k:v for k,v in summary.items() if k not in ('job','last_sample','spray_per_second')},indent=2))


if __name__=='__main__': main()
