"""Area of a declared piecewise-linear phi=.5 reconstruction, in cell^2.

This measures an interface geometry, not particle mass. At domain edges the
nearest cell-centre value is extended constantly through the last half cell.
Average both triangle diagonals, and report their disagreement as an ambiguity.
The quadrature is exact on each linear triangle; the field itself is still a
grid/particle-resolution approximation. Never call this exact material volume.
"""
import numpy as np
from phase_contract import compress


def triangle_fraction(a,b,c,level):
    v=np.sort(np.stack([a,b,c],axis=-1)-level,axis=-1)
    lo,mid,hi=v[...,0],v[...,1],v[...,2]
    f=np.zeros_like(lo);f[lo>=0]=1
    one=(hi>0)&(mid<0)
    f[one]=hi[one]**2/((hi[one]-lo[one])*(hi[one]-mid[one]))
    two=(lo<0)&(mid>=0)&(hi>0)
    f[two]=1-lo[two]**2/((mid[two]-lo[two])*(hi[two]-lo[two]))
    return f


def area(field,level=.5):
    field=np.asarray(field,dtype=float)
    if field.ndim!=2 or not np.isfinite(field).all():raise ValueError('invalid phase field')
    nx,ny=field.shape
    f=np.pad(field,1,mode='edge')
    x=np.r_[0,np.arange(nx)+.5,nx];y=np.r_[0,np.arange(ny)+.5,ny]
    w=np.diff(x)[:,None]*np.diff(y)[None,:]/2
    a,b,c,d=f[:-1,:-1],f[1:,:-1],f[1:,1:],f[:-1,1:]
    first=float(np.sum(w*(triangle_fraction(a,b,c,level)+triangle_fraction(a,c,d,level))))
    second=float(np.sum(w*(triangle_fraction(a,b,d,level)+triangle_fraction(b,c,d,level))))
    return dict(area=.5*(first+second),diagonal_difference=abs(first-second))


def particle_field(s):
    # Independent diagnostic: rebuild cell-centred raw mass, never overwrite
    # the solver's current face fields or loop-carried momentum.
    q=s.pos-.5;b=np.floor(q+.5).astype(np.int64)
    raw,=s._splat((s.nx,s.ny),(b,q-b),[s.masses()])
    return compress(raw,s.rho0_face,s.rho_l,s.rho_g,s.phase_alpha)
