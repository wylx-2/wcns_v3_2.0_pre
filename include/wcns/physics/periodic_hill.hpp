#pragma once

#include <wcns/core/types.hpp>
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

namespace wcns {

// Almeida/TMR geometry, h=1, crest-to-crest period 9. Original polynomial
// coordinates are in units in which h=28. Return height and dy/dx.
inline std::array<Real, 2> periodic_hill_geometry(Real x)
{
    if (!std::isfinite(x)) throw std::invalid_argument("hill coordinate must be finite");
    x -= 9.0 * std::floor(x / 9.0);
    const Real sign = x <= 4.5 ? 1.0 : -1.0;
    const Real s = 28.0 * std::min(x, 9.0 - x);
    if (s >= 54.0) return {{0.0, 0.0}};
    constexpr Real end[] = {9, 14, 20, 30, 40, 54};
    constexpr Real c[][4] = {
        {28, 0, 6.775070969851e-3, -2.124527775800e-3},
        {25.07355893131, .9754803562315, -.1016116352781, .001889794677828},
        {25.79601052357, .8206693007457, -.09055370274339, .001626510569859},
        {40.46435022819, -1.379581654948, .01945884504128, -.000207031893219},
        {17.92461334664, .8743920332081, -.05567361123058, .0006277731764683},
        {56.39011190988, -2.010520359035, .01644919857549, .00002674976141766}};
    int n = 0;
    while (s > end[n]) ++n;
    const auto& a = c[n];
    const Real y = a[0] + s * (a[1] + s * (a[2] + s * a[3]));
    if (n == 0 && y >= 28.0) return {{1.0, 0.0}};
    if (n == 5 && y <= 0.0) return {{0.0, 0.0}};
    return {{y / 28.0, sign * (a[1] + s * (2 * a[2] + 3 * s * a[3]))}};
}

// Quadratic one-sided derivative at a no-slip wall; d is inward distance.
inline Real hill_wall_derivative(Real q1, Real q2, Real d1, Real d2)
{
    if (!(d2 > d1 && d1 > 0)) throw std::invalid_argument("invalid hill wall distances");
    return (d2 * d2 * q1 - d1 * d1 * q2) / (d1 * d2 * (d2 - d1));
}

// Incremental PI control of crest bulk mass velocity, frozen within each RK step.
// Clipping limits the accumulated force as well, providing anti-windup.
inline Real hill_feedback_force(Real force, Real error, Real previous_error,
                                Real dt, Real tau, Real density, Real limit)
{
    if (!(dt > 0 && tau > 0 && density > 0 && limit > 0)
        || !std::isfinite(force + error + previous_error + dt + tau + density + limit))
        throw std::invalid_argument("invalid periodic-hill feedback state");
    return std::clamp(force + density * ((error - previous_error) / tau
                                        + dt * error / (tau * tau)), -limit, limit);
}
} // namespace wcns
