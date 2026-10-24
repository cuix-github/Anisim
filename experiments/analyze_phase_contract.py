"""Read-only phase experiment assessment and particle visualization; no solver."""
import argparse
import json
from pathlib import Path
import subprocess

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image, ImageDraw, ImageFont

FONT=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',20)
SMALL=ImageFont.truetype('C:/Windows/Fonts/arial.ttf',16)
ARMS=('python','bifrost')
CASES=('rest_raw','rest_pure','plunge_raw','plunge_pure','plunge_pure_momentum',
       'plunge_pure_strict','plunge_pure_halfstep')


def load(d,arm,j):
    p=d/arm
    if arm=='python':
        with np.load(p/f'frame.{j:04d}.npz') as a:
            return {k:a[k] for k in a.files}
    return dict(pos=np.load(p/f'frame.{j:04d}.npy')[:,:2],
                vel=np.load(p/f'vel.{j:04d}.npy')[:,:2],
                typ=(np.load(p/f'ph.{j:04d}.npy')>.5).astype(int))


def panel(a,name,t,scale=8):
    w,h=64*scale,80*scale
    pix=np.full((h,w,3),(255,255,255),np.uint8)
    xy=np.rint(a['pos']*scale).astype(int)
    inside=(xy[:,0]>=0)&(xy[:,0]<w)&(xy[:,1]>=0)&(xy[:,1]<h)
    for typ,rgb in [(0,(208,215,225)),(1,(20,67,200))]:
        q=xy[inside&(a['typ']==typ)]
        for dx,dy in ((0,0),(1,0),(0,1),(1,1)):
            pix[np.maximum(h-1-q[:,1]-dy,0),np.minimum(q[:,0]+dx,w-1)]=rgb
    im=Image.new('RGB',(w,h+72),(22,27,34));im.paste(Image.fromarray(pix),(0,72))
    d=ImageDraw.Draw(im);d.text((10,6),f'{name} | t={t:.3f}s',font=FONT,fill='white')
    d.text((10,36),'Blue: water | Grey: air | Escape / spray OFF',font=SMALL,fill=(190,200,211))
    return im


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',required=True,type=Path)
    p.add_argument('--video-case',default='plunge_pure_halfstep')
    a=p.parse_args();root=a.root.resolve();dest=root/'analysis';dest.mkdir(exist_ok=True)
    result=dict(cases={},scope='Uniform 2D, no escape/emission/ST/adaptivity; not production acceptance')
    loaded={}
    for case in CASES:
        d=root/case
        if not (d/'comparison.json').exists(): continue
        job=json.loads((d/'job.json').read_text());comp=json.loads((d/'comparison.json').read_text())
        runtime=json.loads((d/'runtime.json').read_text())
        metrics={arm:json.loads((d/arm/'metrics.json').read_text()) for arm in ARMS}
        loaded[case]=(job,comp,metrics,runtime)
        rec=dict(commit=job['commit'],dt=job['dt'],tolerance=job['tolerance'],
            max_difference=comp['max_difference'],similarity_pass=comp['similarity_pass'],
            momentum_port_matched=job.get('independent_momentum_port',False),
            physics={},alpha_liquid=job.get('alpha_liquid',float(np.float32(.97))),
            alpha_air=job.get('alpha_air',float(np.float32(.9))))
        for arm,rows in metrics.items():
            final=runtime[arm]['final'];allrows=rows+[final]
            rec['physics'][arm]=dict(final=final,
                peak_energy_ratio=max([job['initial']['energy']]+[r['energy'] for r in allrows])/job['initial']['energy'],
                peak_air_rms=max(r['air_rms'] for r in allrows),
                peak_particle_component=max(r['max_component'] for r in allrows),
                wall_seconds_including_startup_and_io=runtime[arm]['wall_seconds'])
        diag=np.stack([np.load(f).ravel() for f in sorted((d/'bifrost').glob('diag.*.npy'))])
        rec['cpp_all_step_diagnostics']=dict(iterations_max=float(diag[:,0].max()),
            reported_residual_max=float(diag[:,1].max()),divergence_max=float(diag[:,2].max()),
            particle_component_max=float(diag[:,3].max()),finite=bool(np.isfinite(diag).all()))
        result['cases'][case]=rec

    chosen=('plunge_raw','plunge_pure_momentum','plunge_pure_strict','plunge_pure_halfstep')
    labels=('Raw phase, initial harness','Pure phase, corrected ports','Pure phase, tighter solve','Pure phase, half dt / scaled blend')
    colours=('#777777','#da7727','#277cb9','#3c9561')
    fig,axes=plt.subplots(2,3,figsize=(16,8),layout='constrained')
    for case,label,col in zip(chosen,labels,colours):
        if case not in loaded:continue
        job,comp,metrics,runtime=loaded[case]
        for arm,style in zip(ARMS,('--','-')):
            rows=metrics[arm]+[dict(time=job['steps']*job['dt'],**runtime[arm]['final'])]
            t=[r['time'] for r in rows]
            for ax,key,den in zip(axes[0],('energy','air_rms','water_rms'),(job['initial']['energy'],1,1)):
                ax.plot(t,[r[key]/den for r in rows],style,c=col,lw=1.2,
                        label=label if arm=='bifrost' else None)
        t=[(r['step']+1)*job['dt'] for r in comp['states']]
        for ax,key in zip(axes[1],('tv','com','height_p99')):
            ax.plot(t,[r[key] for r in comp['states']],c=col,lw=1.2)
    titles=('Particle mechanical energy / initial','Air RMS speed (cells/s)','Water RMS speed (cells/s)',
            'Water histogram difference (TV)','Water COM difference (cells)','Water 99th-height difference (cells)')
    for ax,title in zip(axes.ravel(),titles):
        ax.set_title(title);ax.set_xlabel('Simulation time (s)');ax.grid(alpha=.25)
    for ax,threshold in zip(axes[1],(.05,.2,1)):
        ax.axhline(threshold,c='red',ls=':',label='Predeclared limit')
    axes[0,0].legend(fontsize=8);axes[1,2].legend(fontsize=8)
    fig.suptitle('Uniform 2D diagnostics | solid: Bifrost, dashed: independent NumPy | escape / spray OFF')
    fig.savefig(dest/'diagnostics.png',dpi=150);plt.close(fig)

    # Time refinement changes both dt and the derived blend coefficients to
    # preserve (1-alpha)/dt approximately. It is not a grid-convergence proof.
    coarse,fine='plunge_pure_strict','plunge_pure_halfstep'
    if fine in loaded:
        result['time_refinement']={}
        for arm in ARMS:
            rows=[]
            for k in range(60):
                # Match exact 6s final separately; dump offsets differ by .00125s.
                x=load(root/coarse,arm,k*40);y=load(root/fine,arm,k*80)
                x=x['pos'][x['typ']==1];y=y['pos'][y['typ']==1]
                hs=[np.histogram2d(q[:,0],q[:,1],bins=(16,20),range=((0,64),(0,80)))[0] for q in (x,y)]
                rows.append(dict(time=k*.1,tv=float(.5*abs(hs[0]/hs[0].sum()-hs[1]/hs[1].sum()).sum()),
                    com=float(np.linalg.norm(x.mean(0)-y.mean(0))),
                    p99=float(abs(np.quantile(x[:,1],.99)-np.quantile(y[:,1],.99)))))
            result['time_refinement'][arm]={key:max(r[key] for r in rows) for key in ('tv','com','p99')}
        result['time_refinement']['sampling_note']='Coarse dumps t=.0025+.1k; fine dumps t=.00125+.1k. Compare qualitative sensitivity, not a convergence order.'
    (dest/'summary.json').write_text(json.dumps(result,indent=2))

    case=a.video_case
    if case in loaded:
        job=loaded[case][0];frames=dest/(case+'_frames');frames.mkdir(exist_ok=True)
        selected=[]
        for i,j in enumerate(range(0,job['steps'],job['dump_every'])):
            t=(j+1)*job['dt'];panels=[panel(load(root/case,arm,j),label,t) for arm,label in zip(ARMS,('Independent NumPy','Actual Bifrost graph'))]
            im=Image.new('RGB',(panels[0].width*2+8,panels[0].height+42),(14,19,26))
            for n,pn in enumerate(panels):im.paste(pn,(n*(pn.width+8),0))
            ImageDraw.Draw(im).text((10,im.height-32),'Diagnostic: gates remain open. Equal 2-pixel markers for visibility; not physical radii.',font=SMALL,fill='white')
            im.save(frames/f'{i:04d}.png')
            if i in (20,32,59):
                im.save(dest/f'{case}_{t:.3f}s.png');selected.append(im)
        if selected:
            sheet=Image.new('RGB',(len(selected)*688,502),(14,19,26))
            for n,im in enumerate(selected):sheet.paste(im.resize((688,502)),(n*688,0))
            sheet.save(dest/'particle_sheet.png')
        subprocess.run(['ffmpeg','-hide_banner','-loglevel','error','-y','-framerate',str(1/(job['dt']*job['dump_every'])),
            '-i',str(frames/'%04d.png'),'-c:v','libx264','-crf','18','-pix_fmt','yuv420p','-movflags','+faststart',
            str(dest/'phase_diagnostic.mp4')],check=True)
    print(json.dumps({n:dict(pass_similarity=r['similarity_pass'],difference=r['max_difference']) for n,r in result['cases'].items()},indent=2))


if __name__=='__main__':main()
