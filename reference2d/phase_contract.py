"""Opt-in uniform-grid PF-FLIP phase contract, independent of the C++ code.

Model 1: Eq.7 on raw face mass. Model 2: also enforce pure-cell faces (§3.3).
Cell labels use a separate cell-centred Eq.6 mass splat; this interpolation
choice is explicit because the paper does not specify how it samples cells.
ST, particle adaptivity and solid-support correction require separate gates.
"""
import numpy as np
from pfflip2d import Sim


def compress(raw, rho0, rho_l=1000., rho_g=1., alpha=1.):
    if not (np.isfinite([rho0,rho_l,rho_g,alpha]).all()
            and rho0>0 and rho_l>rho_g>0 and alpha>0):
        raise ValueError('invalid phase calibration / density / alpha')
    floor = np.log(rho_l/rho_g)*rho_g*rho0
    return np.minimum(np.sqrt(np.maximum(np.asarray(raw)-floor,0)/(alpha*rho0*rho_l)),1)


def pure_faces(faces, cell_phi):
    liquid = cell_phi >= .5
    result = [p.copy() for p in faces]
    for axis,p in enumerate(result):
        lo,hi,inside=[slice(None)]*liquid.ndim,[slice(None)]*liquid.ndim,[slice(None)]*liquid.ndim
        lo[axis],hi[axis],inside[axis]=slice(None,-1),slice(1,None),slice(1,-1)
        a,b=liquid[tuple(lo)],liquid[tuple(hi)]
        p[tuple(inside)]=np.where(a==b,a.astype(float),p[tuple(inside)])
    return tuple(result)


class PhaseContractSim(Sim):
    def __init__(self,*args,phase_model=2,phase_alpha=1.,port_rounding=True,pressure_tolerance=1e-6,**kw):
        super().__init__(*args,**kw)
        if phase_model not in (1,2) or not np.isfinite(phase_alpha) or phase_alpha<=0:
            raise ValueError('phase model must be 1/2 with positive alpha')
        if self.st or self.adapt or self.has_solids or self.zone is not None:
            raise ValueError('phase contract currently requires uniform, non-ST, unobstructed domain')
        self.phase_model,self.phase_alpha=phase_model,phase_alpha
        self.port_rounding,self.pressure_tolerance=port_rounding,pressure_tolerance
        self._contract_faces=None

    def _round(self,a):
        return a.astype(np.float32).astype(float) if self.port_rounding else a

    def p2g(self):
        # Public C++ particle momentum ports are float32 before deposition;
        # mean mass/momentum and reconstructed phi cross float32 ports too.
        m=self._round(self.masses())
        valid=(~self.spray).astype(float)
        faces=[]
        for axis,name,shape,frame in zip(range(2),('u','v'),(self.u.shape,self.v.shape),self._face_frames()):
            mom=self._round(m*self.vel[:,axis])
            mass,momentum,support=self._splat(shape,frame,[m,mom,valid])
            setattr(self,'m'+name,mass)
            mm=self._round(mass/(support+float(np.float32(1e-6))))
            mo=self._round(momentum/(support+float(np.float32(1e-6))))
            setattr(self,name,np.divide(mo,mm,out=np.zeros_like(mm),where=mm>0))
            faces.append(compress(mass,self.rho0_face,self.rho_l,self.rho_g,self.phase_alpha))
        if self.phase_model==2:
            q=self.pos-.5;b=np.floor(q+.5).astype(np.int64)
            cm,=self._splat((self.nx,self.ny),(b,q-b),[m])
            cell=compress(cm,self.rho0_face,self.rho_l,self.rho_g,self.phase_alpha)
            faces=pure_faces(faces,cell)
            self.cell_phase=cell
        self._contract_faces=tuple(self._round(p) for p in faces)
        self._enforce_walls()

    def _face_phase(self):
        if self._contract_faces is None:
            raise RuntimeError('p2g must prepare the explicit phase field')
        return self._contract_faces

    def project(self,dt,tol=None,max_iter=4000):
        return super().project(dt,self.pressure_tolerance if tol is None else tol,max_iter)
