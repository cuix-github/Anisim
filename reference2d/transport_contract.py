"""Experimental PF-FLIP §3.7 transport, independent of the production node.

0: retained midpoint transport; 1: local RK3; 2: water extension + midpoint;
3: water extension + local RK3. Pressure, G2P and phase tagging do not change.
The paper does not specify the extension stencil: simultaneous nearest-layer
averaging of face velocities adjacent to liquid is an explicit local choice.
"""
import numpy as np
from phase_contract import PhaseContractSim


def extend_faces(values, valid, layers=3):
    out=values.copy(); known=valid.copy()
    for _ in range(layers):
        total=np.zeros_like(out); count=np.zeros(out.shape,np.int32)
        for axis in range(out.ndim):
            a=[slice(None)]*out.ndim; b=a.copy(); a[axis]=slice(1,None); b[axis]=slice(None,-1)
            a,b=tuple(a),tuple(b)
            total[a]+=np.where(known[b],out[b],0); count[a]+=known[b]
            total[b]+=np.where(known[a],out[a],0); count[b]+=known[a]
        fill=(~known)&(count>0)
        out[fill]=total[fill]/count[fill]; known|=fill
    return out,known


def rk3_local(pos,dt,sample,cfl=1.0):
    """Third-order Ralston RK, stage-velocity CFL checks, no discarded time."""
    p=pos.copy(); remaining=np.full(len(p),dt); advances=0
    while np.any(remaining>1e-12):
        ids=np.flatnonzero(remaining>1e-12); q=p[ids]
        k1=sample(q,ids)
        h=np.minimum(remaining[ids],cfl/np.maximum(np.max(abs(k1),axis=1),1e-30))
        for retry in range(32):
            k2=sample(q+.5*h[:,None]*k1,ids)
            k3=sample(q+.75*h[:,None]*k2,ids)
            peak=np.maximum.reduce([np.max(abs(k),axis=1) for k in (k1,k2,k3)])
            bad=h*peak>cfl*(1+1e-12)
            if not np.any(bad):break
            h[bad]=.95*cfl/peak[bad]
        else: raise RuntimeError('RK3 stage CFL refinement did not converge')
        p[ids]=q+h[:,None]*(2/9*k1+1/3*k2+4/9*k3)
        remaining[ids]-=h; advances+=1
        if advances>4096:raise RuntimeError('RK3 exceeded local step budget; no time discarded')
    return p


class TransportContractSim(PhaseContractSim):
    def __init__(self,*a,transport_model=0,**kw):
        super().__init__(*a,**kw)
        if transport_model not in (0,1,2,3):raise ValueError('unknown transport model')
        self.transport_model=transport_model

    def transport_fields(self):
        values=[self.u,self.v]
        if self.transport_model<2:return values
        # An interior face touching a liquid cell has nonzero phi in mode 2.
        # Wall-normal samples stay exactly zero rather than being extrapolated.
        fields=[]
        for axis,(v,phi) in enumerate(zip(values,self._face_phase())):
            valid=phi>0
            wall=[slice(None)]*2;wall[axis]=[0,-1];valid[tuple(wall)]=True
            f,known=extend_faces(v,valid,3); fields.append(f)
        return fields

    def advect(self,dt,ballistic=None,fly=None):
        if self.transport_model==0:return super().advect(dt,ballistic,fly)
        if ballistic is not None or fly is not None or self.escape or self.st or self.adapt:
            raise ValueError('transport experiment requires no escape/ST/adapt')
        fields=self.transport_fields()
        def sample(q,ids):
            raw=self.sample_velocity(q)
            if self.transport_model>=2:
                mask=self.typ[ids]==1
                for axis,off in enumerate(((0,.5),(.5,0))):
                    raw[mask,axis]=self._sample_faces(fields[axis],np.array(off),q[mask])
            return raw
        if self.transport_model in (1,3):
            self.pos=rk3_local(self.pos,dt,sample)
        else:
            ids=np.arange(len(self.pos))
            mid=self.pos+.5*dt*sample(self.pos,ids)
            self.pos=self.pos+dt*sample(np.clip(mid,.51,None),ids)
        for axis,hi in enumerate((self.nx-.51,self.ny-.51)):
            hit=(self.pos[:,axis]<.51)|(self.pos[:,axis]>hi)
            self.pos[:,axis]=np.clip(self.pos[:,axis],.51,hi);self.vel[hit,axis]=0
