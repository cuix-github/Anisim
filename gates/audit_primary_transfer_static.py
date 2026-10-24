"""Static analytical Eq.6/7/P2G probe, not a simulation or repository import."""
import argparse
import json
import math
from pathlib import Path

K = sum(max(1-x*x-y*y, 0)**3
        for x in (-.75, -.25, .25, .75)
        for y in (-.75, -.25, .25, .75))
RHO_L = 1000.0
RHO_G = 1.0

def kernel(x, y, q):
    return max(1-(x-q[0])**2-(y-q[1])**2, 0)**3

def phi(raw):
    return min(1.0, math.sqrt(max(raw-math.log(1000)*K, 0)/(1000*K)))

def cell(i, j, q):
    return phi(K+1000*kernel(i+.5, j+.5, q)) >= .5

def at_face(i, j, axis, q):
    x, y = (i,j+.5) if axis == 0 else (i+.5,j)
    w = kernel(x,y,q)
    own = 1000*w/(K+1000*w)
    a = cell(i-(axis==0),j-(axis==1),q)
    b = cell(i,j,q)
    raw = phi(K+1000*w)
    pure = float(a) if a == b else raw
    return own,raw,pure

def sample(axis,q):
    ox,oy = (0,.5) if axis == 0 else (.5,0)
    x,y=q[0]-ox,q[1]-oy
    i,j=math.floor(x),math.floor(y)
    f,g=x-i,y-j
    ans=[0.,0.,0.]
    for di,wx in ((0,1-f),(1,f)):
        for dj,wy in ((0,1-g),(1,g)):
            vals=at_face(i+di,j+dj,axis,q)
            for k in range(3):ans[k]+=wx*wy*vals[k]
    return ans

rows=[]
for name,q in [('quarter',(4.25,4.25)),('cell_center',(4.5,4.5)),
               ('u_face_center',(4,4.5)),('grid_vertex',(4,4))]:
    u,v=sample(0,q),sample(1,q)
    rows.append(dict(position=name,own_u_fraction=u[0],own_v_fraction=v[0],
        raw_sample_phi=.5*(u[1]+v[1]),pure_sample_phi=.5*(u[2]+v[2]),
        escape_pure=.5*(u[2]+v[2])<.3))
scan=[]
for i in range(41):
    for j in range(41):
        q=(4+(i+.5)/41,4+(j+.5)/41)
        u,v=sample(0,q),sample(1,q)
        scan.append((.5*(u[2]+v[2]),u[0],v[0]))
escaped=[s for s in scan if s[0]<.3]
report=dict(scope='One added liquid primary in stationary, uniformly seeded air; no projection or evolution',
    regular_kernel_sum=K,examples=rows,scan_samples=len(scan),escape_samples=len(escaped),
    escaped_own_component_fraction_min=min(min(s[1:]) for s in escaped),
    escaped_own_component_fraction_max=max(max(s[1:]) for s in escaped),
    density_global_normalization='log(n/(c*n_ref))-mean(log(n/(c*n_ref))) equals log(n/n_ref)-mean(log(n/n_ref)) exactly algebraically')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--out',type=Path,required=True)
output=parser.parse_args().out
if output.exists():raise FileExistsError(output)
output.parent.mkdir(parents=True,exist_ok=True)
output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
print(json.dumps(report,indent=2))
