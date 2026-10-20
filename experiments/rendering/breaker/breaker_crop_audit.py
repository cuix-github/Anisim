"""Inspect existing position/normal AOVs for visible crop-side surfaces; no render."""
import json
from pathlib import Path
import sys
import numpy as np

sys.path.insert(0,str(Path(__file__).resolve().parent))
import breaker_deep_water as base
from analyze_breaker_deep_water import pfm

OUT=base.ROOT/'out/runs/breaker_foam_views_20261019'


def main():
    assert not base.git('status','--porcelain')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    views=['wide','oblique','later26']
    fig,axes=plt.subplots(3,1,figsize=(12,18),facecolor='#171b21')
    records={}
    for ax,view in zip(axes,views):
        folder=OUT/view/'A'
        assert json.loads((folder/'image_check.json').read_text())['status']=='passed'
        p=pfm(folder/'P.pfm')*.25;n=pfm(folder/'N.pfm');tr=pfm(folder/'transmission.pfm')
        # Restrict to transmitted surface samples in the water crop. This is
        # deliberately not an object-ID matte or a count of all water pixels.
        water=(tr.sum(2)>1e-6)&(p[:,:,0]>119)&(p[:,:,0]<173)&(p[:,:,1]>5)&(p[:,:,1]<25)&(p[:,:,2]>59)&(p[:,:,2]<79)
        near=(p[:,:,2]<60.25)&(np.abs(n[:,:,2])>.7)
        category=np.zeros(p.shape[:2]+(3,),np.float32)
        category[water]=[.1,.45,.65];category[water&near]=[1,.65,.1]
        ax.imshow(category);ax.axis('off')
        fraction=float((water&near).sum()/water.sum())
        ax.set_title(f'{view} | {fraction:.1%} of selected transmitted samples near front crop side',color='white')
        records[view]=dict(selected_transmitted_samples=int(water.sum()),front_crop_samples=int((water&near).sum()),
                           front_fraction=fraction,z_quantiles_m=np.quantile(p[:,:,2][water],[0,.25,.5,.75,1]).tolist(),
                           source_sha256={k:base.sha(folder/(k+'.pfm')) for k in ['P','N','transmission']})
    fig.suptitle('Existing AOV data only: orange = near z=60 m cut, |Nz| > 0.7; blue = other selected water\n'
                 'No geometry change or new render; excludes pixels without transmission contribution',color='white')
    fig.tight_layout(rect=(0,0,1,.96));fig.savefig(OUT/'crop_audit.png',dpi=100);plt.close(fig)
    base.save(OUT/'crop_audit.json',dict(status='passed',git_sha=base.git('rev-parse','HEAD'),
              method='World-space Arnold P/N AOVs, simulation metres = scene units / 4',
              selection='transmission RGB sum > 1e-6; 119<x<173, 5<y<25, 59<z<79 metres',
              crop_side_selection='z < 60.25 m and abs(Nz) > 0.7; crop starts at z=60 m',
              limitation='Thresholded AOV samples, not exact object-ID attribution; no claim to explain all scale perception.',
              views=records))


if __name__=='__main__':main()
