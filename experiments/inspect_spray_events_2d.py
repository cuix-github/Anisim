"""Plot existing 2D particle caches and comparison metrics; never simulates."""
import argparse
import json
from pathlib import Path
import subprocess
import sys

import numpy as np
from PIL import Image, ImageDraw, ImageFont
from spray_events_2d import NX, NY, load_frame

FONT = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", 20)
SMALL = ImageFont.truetype("C:/Windows/Fonts/arial.ttf", 16)
COLOURS = [("air", (63,72,84)), ("water",(48,125,230)),
           ("escaped water",(255,153,40)), ("emitted spray",(238,248,250))]


def panel(a, title, t, scale=8, bounds=(0,NX,0,NY)):
    x0,x1,y0,y1 = bounds
    W,H = round((x1-x0)*scale), round((y1-y0)*scale)
    pix = np.full((H,W,3), (15,20,27), np.uint8)
    xy = np.rint((a["pos"]-[x0,y0])*scale).astype(int)
    inside = (xy[:,0]>=0)&(xy[:,0]<W)&(xy[:,1]>=0)&(xy[:,1]<H)
    masks = [(a["typ"]==0)&~a["spr"], (a["typ"]==1)&~a["spr"]&~a["esc"],
             (a["typ"]==1)&~a["spr"]&a["esc"], a["spr"]]
    for (_,rgb), mask in zip(COLOURS,masks):
        p = xy[inside&mask]
        pix[H-1-p[:,1],p[:,0]] = rgb
        # Water/air use one pixel, spray two: display markers, not droplet radii.
        if rgb in (COLOURS[2][1],COLOURS[3][1]):
            pix[H-1-p[:,1],np.minimum(p[:,0]+1,W-1)] = rgb
    im=Image.new("RGB",(W,H+66),(22,27,34)); im.paste(Image.fromarray(pix),(0,66))
    d=ImageDraw.Draw(im)
    d.text((10,6),f"{title} | t = {t:.2f} s",font=FONT,fill="white")
    d.text((10,34),f"escaped: {int(masks[2].sum())}   spray: {int(a['spr'].sum())}",font=SMALL,fill=(195,204,217))
    return im


def main():
    p=argparse.ArgumentParser(); p.add_argument("--out",required=True,type=Path)
    a=p.parse_args(); root=a.out.resolve(); job=json.loads((root/"job.json").read_text())
    dt,every,steps=job["dt"],job["dump_every"],job["steps"]
    dest=root/"visuals"; dest.mkdir(exist_ok=True)
    frames=dest/"frames"; frames.mkdir(exist_ok=True)
    selected=[]; targets=[1.8,3.2,6.,14.]
    for i,j in enumerate(range(0,steps,every)):
        left=panel(load_frame(root,"python",j),"Independent NumPy",(j+1)*dt)
        right=panel(load_frame(root,"bifrost",j),"Bifrost C++ graph",(j+1)*dt)
        im=Image.new("RGB",(left.width+right.width+8,left.height+64),(15,20,27))
        im.paste(left,(0,0)); im.paste(right,(left.width+8,0))
        d=ImageDraw.Draw(im)
        x=12
        for label,rgb in COLOURS:
            d.rectangle((x,left.height+10,x+12,left.height+22),fill=rgb)
            d.text((x+18,left.height+6),label,font=SMALL,fill="white")
            x+=220
        d.text((12,left.height+34),"Same seed and parameters | marker sizes are for visibility | no surface reconstruction",font=SMALL,fill=(175,188,205))
        im.save(frames/f"frame.{i:04d}.png")
        if any(abs((j+1)*dt-t)<every*dt*.51 for t in targets):
            im.save(dest/f"comparison_{(j+1)*dt:.2f}s.png")
            selected.append(im)
    fps=1/(dt*every)
    subprocess.run(["ffmpeg","-hide_banner","-loglevel","error","-y","-framerate",str(fps),
        "-i",str(frames/"frame.%04d.png"),"-c:v","libx264","-crf","18","-pix_fmt","yuv420p",
        "-movflags","+faststart",str(dest/"python_bifrost_20s.mp4")],check=True)
    if selected:
        thumbs=[im.resize((774,578)) for im in selected]
        sheet=Image.new("RGB",(1548,578*((len(thumbs)+1)//2)),(15,20,27))
        for i,im in enumerate(thumbs): sheet.paste(im,((i%2)*774,(i//2)*578))
        sheet.save(dest/"comparison_sheet.png")
    print(str(dest/"python_bifrost_20s.mp4"))


if __name__=="__main__": main()
