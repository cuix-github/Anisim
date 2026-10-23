#include "../nodes/src/PFFlipTransport.h"
#include <iostream>
#include <limits>
#include <stdexcept>
using namespace PFFlip::Transport;
int main() {
    const int nx=14,ny=13;Vec v(nx*ny,99),phi(nx*ny,0);
    for (int i=0;i<nx;++i) for (int j=0;j<ny;++j) {
        const int k=i*ny+j;
        if (i>=4 && i<8 && j>=4 && j<8) {phi[k]=1;v[k]=2.75;}
        if (i==0 || i==nx-1) v[k]=2.75;
    }
    const auto ext=extend2(v,phi,nx,ny,0);
    for (int i=0;i<nx;++i) for (int j=0;j<ny;++j) {
        int dist=std::max({4-i,0,i-7})+std::max({4-j,0,j-7});
        dist=std::min({dist,i,nx-1-i});
        if (ext[i*ny+j]!=(dist<=3?2.75:99)) throw std::runtime_error("extension");
    }
    double x=2,y=3;
    if (!rk3(x,y,2.,[](double,double){return std::array<double,2>{2,-.3};}) ||
        std::abs(x-6)>1e-13 || std::abs(y-2.4)>1e-13) throw std::runtime_error("constant transport");
    double e[3];int n=0;
    for (double h : {.4,.2,.1}) {
        x=1;y=0;
        for (int j=0;j<int(std::round(2/h));++j)
            if (!rk3(x,y,h,[](double a,double b){return std::array<double,2>{-b,a};})) return 1;
        e[n++]=std::hypot(x-std::cos(2.),y-std::sin(2.));
    }
    if (e[0]/e[1]<6 || e[1]/e[2]<6) throw std::runtime_error("third order");
    if (rk3(x,y,.1,[](double,double){return std::array<double,2>{std::numeric_limits<double>::quiet_NaN(),0};}))
        throw std::runtime_error("nonfinite accepted");
    std::cout.precision(17);
    std::cout<<"{\"rotation_errors\":["<<e[0]<<","<<e[1]<<","<<e[2]<<"],\"passed\":true}";
}
