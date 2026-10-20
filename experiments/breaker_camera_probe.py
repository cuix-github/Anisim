"""Offline primary-ray mesh intersection check; no Arnold render or shading."""
import sys
import numpy as np
import maya.standalone
maya.standalone.initialize(name='python')
import maya.api.OpenMaya as om
from PIL import Image, ImageDraw
import breaker_six_views as study
import breaker_deep_water as base

name=sys.argv[1]
root=study.OUT/'frames'/name
data=np.load(root/'mesh.npz');p=data['points'];t=data['triangles']
fn=om.MFnMesh();fn.create(om.MPointArray(p.tolist()),[3]*len(t),t.ravel().tolist())
accel=fn.autoUniformGridParams()
W,H=240,135
sheet=Image.new('RGB',(W*2,H+24),(30,30,35));draw=ImageDraw.Draw(sheet)
result={}
for col,view in enumerate(['wide','detail']):
    spec=study.view_spec(name,view);pitch,yaw,_=np.radians(spec['rotation'])
    back=np.array([np.sin(yaw)*np.cos(pitch),-np.sin(pitch),np.cos(yaw)*np.cos(pitch)])
    right=np.array([np.cos(yaw),0,-np.sin(yaw)]);up=np.cross(back,right)
    eye=np.array(spec['target_m'])*4+500*back
    direction=om.MFloatVector(*(-back))
    colors=np.zeros((H,W,3),np.uint8);hits=0;cut=0
    for y in range(H):
        for x in range(W):
            origin=eye+4*spec['width_m']*((x+.5)/W-.5)*right+4*spec['width_m']*(.5-(y+.5)/H)*H/W*up
            hit=fn.closestIntersection(om.MFloatPoint(*origin),direction,om.MSpace.kObject,2000,False,accelParams=accel)
            if hit is None:continue
            xyz=np.array(list(hit[0])[:3])*.25;hits+=1
            boundary=(xyz[0]<study.FRAMES[name]['xmin']+.4 or xyz[2]<10.4 or xyz[2]>77.6) and xyz[1]>(xyz[0]-36)*.1+.3
            if boundary:colors[y,x]=[255,30,30];cut+=1
            else:colors[y,x]=[40,100+int(np.clip((xyz[1]-10)*20,0,150)),180]
    sheet.paste(Image.fromarray(colors),(col*W,24));draw.text((col*W+5,5),name+' '+view,fill='white')
    result[view]=dict(rays=W*H,water_hits=hits,visible_cut_hits=cut,fraction=cut/max(hits,1))
sheet.resize((960,318)).save(root/'camera_ray_check.png')
base.save(root/'camera_ray_check.json',result)
print(result,flush=True)
base.exit_worker(0)
