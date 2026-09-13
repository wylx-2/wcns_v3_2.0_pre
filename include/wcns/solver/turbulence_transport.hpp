#pragma once

#include <wcns/core/types.hpp>
#include <wcns/solver/turbulence_fields.hpp>

#include <vector>

namespace wcns {

struct TurbulenceFaceFlux {
    std::vector<Real> advective;
    std::vector<Real> diffusive;
    std::vector<Real> net;
};

// Transported variables are specific model scalars. Their conservative
// advective flux is mass_flux * phi_upwind and never enters Euler eigenvectors.
[[nodiscard]] TurbulenceFaceFlux turbulence_face_flux(
    Real mass_flux,
    const std::vector<Real>& left_specific_values,
    const std::vector<Real>& right_specific_values,
    const std::vector<Real>& diffusion_coefficients,
    const std::vector<Real>& normal_gradients);

enum class TurbulenceBoundaryRule {
    Extrapolate,
    Dirichlet,
};

[[nodiscard]] std::vector<Real>
turbulence_boundary_ghost(const std::vector<Real>& interior,
                          const std::vector<Real>& boundary_values,
                          TurbulenceBoundaryRule rule);

[[nodiscard]] std::vector<Real> admissible_turbulence_update(
    const std::vector<Real>& current,
    const std::vector<Real>& increment,
    const std::vector<TurbulenceFieldDescriptor>& descriptors,
    Real relaxation = 1.0);

} // namespace wcns
