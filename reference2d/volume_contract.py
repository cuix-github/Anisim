"""Opt-in log number-density relaxation; see docs/verification-route.md.

Four primary particles/cell, closed uniform 2D, no escaped/auxiliary particles.
Reference boundary support uses the same .51-cell particle clamp. The spatial
mean is removed to satisfy the closed Neumann pressure compatibility condition.
These are explicit implementation choices, not claimed author tuning.
"""
import numpy as np
from phase_contract import PhaseContractSim


def target_divergence(current, reference, tau):
    current,reference=np.asarray(current),np.asarray(reference)
    if (current.shape!=reference.shape or current.size==0 or not np.isfinite(tau) or tau<=0
        or not np.isfinite(current).all() or not np.isfinite(reference).all()
        or (current<=0).any() or (reference<=0).any()):
        raise ValueError('invalid or empty volume number-density support')
    d=np.log(current/reference)/tau
    return d-d.mean()


class VolumeContractSim(PhaseContractSim):
    def __init__(self,*args,volume_tau=0.,**kwargs):
        super().__init__(*args,**kwargs)
        if not np.isfinite(volume_tau) or volume_tau<0:
            raise ValueError('invalid volume relaxation time')
        if volume_tau and (self.phase_model!=2 or self.escape or self.droplet or self.emit):
            raise ValueError('volume contract requires phase 2 without escape/spray')
        self.volume_tau=volume_tau
        self._density_reference=None

    def number_density(self):
        q=self.pos-.5; b=np.floor(q+.5).astype(np.int64)
        return self._splat((self.nx,self.ny),(b,q-b),[np.ones(len(self.pos))])[0]

    def density_reference(self):
        if self._density_reference is None:
            probe=PhaseContractSim(self.nx,self.ny)
            probe.seed(lambda x,y:y<0,jitter=0)
            probe.pos=np.maximum(.51,np.minimum([self.nx-.51,self.ny-.51],probe.pos))
            probe.pos=probe.pos.astype(np.float32).astype(float)
            q=probe.pos-.5; b=np.floor(q+.5).astype(np.int64)
            self._density_reference=probe._splat((self.nx,self.ny),(b,q-b),[np.ones(len(q))])[0]
        return self._density_reference

    def project(self,dt,tol=None,max_iter=4000):
        if self.volume_tau and (len(self.pos)!=4*self.nx*self.ny or self.spray.any()
            or (self.radius>0).any() or not np.isfinite(self.pos).all()):
            raise ValueError('volume contract requires finite primary particles, four per initial cell, with no droplets')
        self._target_divergence=(target_divergence(self.number_density(),self.density_reference(),self.volume_tau)
                                 if self.volume_tau else None)
        return super().project(dt,tol,max_iter)
