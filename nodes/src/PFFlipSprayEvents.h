#pragma once
// Public, independently specified sampler; see docs/spray-events-2d.md.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

namespace SprayEvents {
constexpr uint64_t gold = UINT64_C(0x9E3779B97F4A7C15);
inline uint64_t mix(uint64_t s) {
    uint64_t z = s + gold;
    z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
    return z ^ (z >> 31);
}
inline uint64_t key(uint64_t seed, uint64_t step, uint64_t source,
                    uint64_t event=0, uint64_t child=0) {
    uint64_t k=mix(seed);
    for (uint64_t v : {step,source,event,child}) k=mix(k ^ mix(v));
    return k;
}
inline double uniform(uint64_t k, int channel) {
    return static_cast<double>(mix(k + channel*gold) >> 11) * 0x1.0p-53;
}
inline std::array<double, 3> ball(uint64_t k, int ch, bool flat) {
    const double a=uniform(k,ch), theta=6.2831853071795864769*uniform(k,ch+1);
    if (flat) { const double r=std::sqrt(a); return {r*std::cos(theta),r*std::sin(theta),0}; }
    const double z=2*a-1, xy=std::sqrt(std::max(1-z*z,0.0));
    const double r=std::pow(uniform(k,ch+2),1.0/3.0);
    return {r*xy*std::cos(theta),r*xy*std::sin(theta),r*z};
}
inline double rawRadius(uint64_t k, int n, double sigma, double radius) {
    const double z=std::sqrt(-2*std::log(1-uniform(k,1))) *
                   std::cos(6.2831853071795864769*uniform(k,2));
    return radius*std::pow(static_cast<double>(n),-1.0/3.0)*std::exp(sigma*z-1.5*sigma*sigma);
}
} // namespace SprayEvents
