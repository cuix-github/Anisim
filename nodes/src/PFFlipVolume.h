#pragma once
#include <algorithm>
#include <cmath>
#include <vector>

// Experimental closed, uniform, four-particles/cell 2D density relaxation.
// The log law is from Losasso 2008; wall calibration and removal of its
// constant mode are explicit local choices, NOT the unreleased PF solver.
namespace PFFlip { namespace Volume {
using Vec = std::vector<double>;
inline void deposit(Vec& n, int nx, int ny, double x, double y) {
    const double qx=x-.5, qy=y-.5;
    const int bx=int(std::floor(qx+.5)), by=int(std::floor(qy+.5));
    for (int di=-1;di<=1;++di) for (int dj=-1;dj<=1;++dj) {
        const int i=bx+di, j=by+dj;
        if (i<0 || i>=nx || j<0 || j>=ny) continue;
        const double dx=qx-i, dy=qy-j, w=std::max(1-dx*dx-dy*dy,0.);
        n[size_t(i)*ny+j]+=w*w*w;
    }
}
inline Vec reference(int nx, int ny) {
    Vec n(size_t(nx)*ny,0.);
    for (double oy : {.25,.75}) for (double ox : {.25,.75})
        for (int i=0;i<nx;++i) for (int j=0;j<ny;++j) {
            // Match float32 initial/checkpoint ports at the admissible boundary.
            const double x=float(std::max(.51,std::min(nx-.51,i+ox)));
            const double y=float(std::max(.51,std::min(ny-.51,j+oy)));
            deposit(n,nx,ny,x,y);
        }
    return n;
}
inline bool target(const Vec& n, const Vec& ref, double tau, Vec& divergence) {
    if (n.empty() || n.size()!=ref.size() || !std::isfinite(tau) || tau<=0) return false;
    divergence.resize(n.size()); double mean=0;
    for (size_t i=0;i<n.size();++i) {
        if (!std::isfinite(n[i]) || !std::isfinite(ref[i]) || n[i]<=0 || ref[i]<=0) return false;
        divergence[i]=std::log(n[i]/ref[i])/tau;
        mean+=divergence[i];
    }
    mean/=double(n.size());
    for (double& x : divergence) x-=mean;
    return true;
}

// Eq. (5) finite-difference alternative with a seconds-based relaxation time.
// Local zero sampling is permitted; no epsilon or change to the legacy log law.
inline bool targetFinite(const Vec& n, const Vec& ref, double tau, Vec& divergence) {
    if (n.empty() || n.size()!=ref.size() || !std::isfinite(tau) || tau<=0) return false;
    divergence.resize(n.size()); double mean=0; bool supported=false;
    for (size_t i=0;i<n.size();++i) {
        if (!std::isfinite(n[i]) || !std::isfinite(ref[i]) || n[i]<0 || ref[i]<=0) return false;
        supported |= n[i]>0;
        divergence[i]=(n[i]/ref[i]-1.)/tau;
        if (!std::isfinite(divergence[i])) return false;
        mean+=divergence[i];
    }
    if (!supported || !std::isfinite(mean)) return false;
    mean/=double(n.size());
    for (double& x : divergence) {
        x-=mean;
        if (!std::isfinite(x)) return false;
    }
    return true;
}
}}
