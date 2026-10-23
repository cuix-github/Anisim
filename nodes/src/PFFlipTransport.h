#pragma once
#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

// Experimental uniform 2D PF-FLIP §3.7 transport. Pressure and G2P are untouched.
// The simultaneous face-layer averaging is our stated implementation choice;
// the paper specifies a few cells of liquid extension, not its precise stencil.
namespace PFFlip { namespace Transport {
using Vec = std::vector<double>;
inline Vec extend2(const Vec& values, const Vec& phase, int nx, int ny, int axis) {
    Vec out=values;
    std::vector<char> known(values.size());
    for (int i=0;i<nx;++i) for (int j=0;j<ny;++j) {
        const size_t k=static_cast<size_t>(i)*ny+j;
        known[k]=phase[k]>0 || (axis==0 ? i==0 || i==nx-1 : j==0 || j==ny-1);
    }
    for (int layer=0;layer<3;++layer) {
        Vec next=out; auto mask=known;
        for (int i=0;i<nx;++i) for (int j=0;j<ny;++j) {
            const size_t k=static_cast<size_t>(i)*ny+j;
            if (known[k]) continue;
            double sum=0; int count=0;
            auto add=[&](int x,int y) {
                if (x<0 || x>=nx || y<0 || y>=ny) return;
                const size_t q=static_cast<size_t>(x)*ny+y;
                if (known[q]) {sum+=out[q];++count;}
            };
            add(i-1,j);add(i+1,j);add(i,j-1);add(i,j+1);
            if (count) {next[k]=sum/count;mask[k]=1;}
        }
        out.swap(next); known.swap(mask);
    }
    return out;
}

template<class Sample>
bool rk3(double& x,double& y,double dt,const Sample& sample) {
    double rem=dt; int advances=0;
    while (rem>1e-12) {
        const auto k1=sample(x,y);
        double h=std::min(rem,1.0/std::max({std::abs(k1[0]),std::abs(k1[1]),1e-30}));
        std::array<double,2> k2{},k3{}; bool accepted=false;
        for (int retry=0;retry<32;++retry) {
            k2=sample(x+.5*h*k1[0],y+.5*h*k1[1]);
            k3=sample(x+.75*h*k2[0],y+.75*h*k2[1]);
            double peak=0;
            for (const auto& k : {k1,k2,k3}) for (double v:k) {
                if (!std::isfinite(v)) return false;
                peak=std::max(peak,std::abs(v));
            }
            if (h*peak<=1.0*(1+1e-12)) {accepted=true;break;}
            h=.95/peak;
        }
        if (!accepted) return false;
        x+=h*((2.0/9)*k1[0]+(1.0/3)*k2[0]+(4.0/9)*k3[0]);
        y+=h*((2.0/9)*k1[1]+(1.0/3)*k2[1]+(4.0/9)*k3[1]);
        rem-=h;
        if (++advances>4096) return false; // refuse, never silently discard time
    }
    return true;
}
}} // namespace PFFlip::Transport
