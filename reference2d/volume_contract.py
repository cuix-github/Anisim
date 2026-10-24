"""Opt-in log number-density relaxation; see docs/verification-route.md.

Four primary particles/cell, closed uniform 2D, no escaped/auxiliary particles.
Reference boundary support uses the same .51-cell particle clamp. The spatial
mean is removed to satisfy the closed Neumann pressure compatibility condition.
These are explicit implementation choices, not claimed author tuning.
"""
import numpy as np
from phase_contract import PhaseContractSim
from transport_contract import TransportContractSim


def target_divergence(current, reference, tau):
    current,reference=np.asarray(current),np.asarray(reference)
    if (current.shape!=reference.shape or current.size==0 or not np.isfinite(tau) or tau<=0
        or not np.isfinite(current).all() or not np.isfinite(reference).all()
        or (current<=0).any() or (reference<=0).any()):
        raise ValueError('invalid or empty volume number-density support')
    d=np.log(current/reference)/tau
    return d-d.mean()


def finite_target_divergence(current, reference, tau):
    """Finite Eq. (5) candidate, independently evaluated using NumPy arrays."""
    current,reference=np.asarray(current,dtype=float),np.asarray(reference,dtype=float)
    if (current.shape!=reference.shape or current.size==0 or not np.isfinite(tau) or tau<=0
        or not np.isfinite(current).all() or not np.isfinite(reference).all()
        or (current<0).any() or (reference<=0).any() or not (current>0).any()):
        raise ValueError('invalid finite-density support')
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        raw=(current/reference-1.)/tau
        result=raw-raw.mean()
    if not np.isfinite(result).all():
        raise ValueError('nonfinite finite-density target')
    return result


class VolumeContractSim(TransportContractSim):
    def __init__(self,*args,volume_tau=0.,volume_model=0,**kwargs):
        super().__init__(*args,**kwargs)
        if volume_model not in (0,1) or not np.isfinite(volume_tau) or volume_tau<0:
            raise ValueError('invalid volume relaxation time')
        if (volume_tau or volume_model) and (self.phase_model!=2 or
            ((self.escape or self.droplet or self.emit) and not self.droplet_transfer_model)):
            raise ValueError('volume contract requires phase 2 without escape/spray')
        self.volume_tau=volume_tau
        self.volume_model=volume_model
        self._density_reference=None

    def number_density(self):
        q=self.pos-.5; b=np.floor(q+.5).astype(np.int64)
        active=(self.masses()>0).astype(float) if self.droplet_transfer_model else np.ones(len(self.pos))
        return self._splat((self.nx,self.ny),(b,q-b),[active])[0]

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
        if self.volume_tau and (not np.isfinite(self.pos).all() or
            (not self.droplet_transfer_model and (len(self.pos)!=4*self.nx*self.ny
             or self.spray.any() or (self.radius>0).any()))):
            raise ValueError('volume contract requires finite primary particles, four per initial cell, with no droplets')
        formula=finite_target_divergence if self.volume_model else target_divergence
        self._target_divergence=(formula(self.number_density(),self.density_reference(),self.volume_tau)
                                 if self.volume_tau else None)
        return super().project(dt,tol,max_iter)
