#pragma once
#include <wcns/core/types.hpp>
#include <array>
#include <algorithm>
#include <cmath>

namespace wcns {
// Fixed physical inlet for the supplied chapter-5 ramp. The table is a
// compressible laminar similarity solution, NOT a turbulent wall-law profile.
inline std::array<Real,5> ramp_inlet(Real y) {
    static const std::array<Real,3> table[] = {
#include <wcns/physics/ramp_inlet_table.inc>
    };
    const Real a=std::clamp(y*256,Real(0),Real(1024));
    const int i=std::min(1023,static_cast<int>(a)); const Real w=a-i;
    std::array<Real,3> q{};
    for(int m=0;m<3;++m) q[m]=(1-w)*table[i][m]+w*table[i+1][m];
    return {1/q[2],q[0],q[1],0,q[2]};
}

// Smooth, counter-flow trip, with reproducible spanwise/time modulation.
// This is an explicitly documented implementation choice, not an unpublished
// reconstruction of Porter & Poggie's exact forcing waveform.
inline Real ramp_trip(Real x,Real y,Real z,Real time,Real amplitude,Real span) {
    if(x<=2.33 || x>=2.67 || y<=0 || y>=.01) return 0;
    constexpr Real pi=3.14159265358979323846;
    const Real a=std::sin(pi*(x-2.33)/.34), b=std::sin(pi*y/.01);
    const Real phase=2*pi*z/span;
    return -amplitude*a*a*b*b*(1+.2*std::sin(3*phase+2*pi*.7*time)
        +.2*std::sin(7*phase-2*pi*1.1*time));
}
} // namespace wcns
