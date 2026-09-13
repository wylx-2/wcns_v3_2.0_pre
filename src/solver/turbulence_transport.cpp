#include <wcns/solver/turbulence_transport.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

void require_finite(const std::vector<Real>& values, const char* label)
{
    if (!std::all_of(values.begin(), values.end(), [](Real value) { return std::isfinite(value); })) {
        throw std::invalid_argument(std::string(label) + " contains non-finite values");
    }
}

} // namespace

TurbulenceFaceFlux turbulence_face_flux(
    Real mass_flux,
    const std::vector<Real>& left_specific_values,
    const std::vector<Real>& right_specific_values,
    const std::vector<Real>& diffusion_coefficients,
    const std::vector<Real>& normal_gradients)
{
    const auto count = left_specific_values.size();
    if (!std::isfinite(mass_flux) || right_specific_values.size() != count
        || diffusion_coefficients.size() != count || normal_gradients.size() != count) {
        throw std::invalid_argument("turbulence face-flux metadata is invalid");
    }
    require_finite(left_specific_values, "left turbulence trace");
    require_finite(right_specific_values, "right turbulence trace");
    require_finite(diffusion_coefficients, "turbulence diffusion coefficient");
    require_finite(normal_gradients, "turbulence normal gradient");
    if (!std::all_of(diffusion_coefficients.begin(),
                     diffusion_coefficients.end(),
                     [](Real value) { return value >= 0.0; })) {
        throw std::invalid_argument("turbulence diffusion coefficient must be non-negative");
    }

    TurbulenceFaceFlux result;
    result.advective.resize(count);
    result.diffusive.resize(count);
    result.net.resize(count);
    const auto& upwind = mass_flux >= 0.0 ? left_specific_values : right_specific_values;
    for (std::size_t variable = 0; variable < count; ++variable) {
        result.advective[variable] = mass_flux * upwind[variable];
        result.diffusive[variable]
            = diffusion_coefficients[variable] * normal_gradients[variable];
        result.net[variable] = result.advective[variable] - result.diffusive[variable];
    }
    return result;
}

std::vector<Real> turbulence_boundary_ghost(const std::vector<Real>& interior,
                                            const std::vector<Real>& boundary_values,
                                            TurbulenceBoundaryRule rule)
{
    if (interior.size() != boundary_values.size()) {
        throw std::invalid_argument("turbulence boundary vectors have different sizes");
    }
    require_finite(interior, "turbulence boundary interior");
    require_finite(boundary_values, "turbulence boundary value");
    if (rule == TurbulenceBoundaryRule::Extrapolate) return interior;
    std::vector<Real> result(interior.size());
    for (std::size_t variable = 0; variable < result.size(); ++variable) {
        result[variable] = 2.0 * boundary_values[variable] - interior[variable];
    }
    require_finite(result, "turbulence boundary ghost");
    return result;
}

std::vector<Real> admissible_turbulence_update(
    const std::vector<Real>& current,
    const std::vector<Real>& increment,
    const std::vector<TurbulenceFieldDescriptor>& descriptors,
    Real relaxation)
{
    if (current.size() != increment.size() || current.size() != descriptors.size()
        || !std::isfinite(relaxation) || relaxation <= 0.0 || relaxation > 1.0) {
        throw std::invalid_argument("turbulence update metadata is invalid");
    }
    require_finite(current, "current turbulence state");
    require_finite(increment, "turbulence increment");
    std::vector<Real> result(current.size());
    for (std::size_t variable = 0; variable < result.size(); ++variable) {
        descriptors[variable].validate();
        result[variable] = current[variable] + relaxation * increment[variable];
        const bool admissible = std::isfinite(result[variable])
            && (descriptors[variable].strictly_positive
                    ? result[variable] > descriptors[variable].lower_bound
                    : result[variable] >= descriptors[variable].lower_bound);
        if (!admissible) {
            throw std::runtime_error("turbulence update violates bound for "
                                     + descriptors[variable].name);
        }
    }
    return result;
}

} // namespace wcns
