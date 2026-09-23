#pragma once

#include <wcns/solver/turbulence_model.hpp>

#include <array>

namespace wcns {

using VelocityGradient = std::array<std::array<Real, 3>, 3>;

struct LesPointEvaluation {
    SymmetricStress deviatoric_strain;
    SymmetricStress physical_sgs_stress;
    Real strain_magnitude = 0.0;
    Real eddy_viscosity = 0.0;
    Real energy_transfer = 0.0;
};

[[nodiscard]] SymmetricStress les_deviatoric_strain(const VelocityGradient& gradient);
[[nodiscard]] Real les_symmetric_inner_product(const SymmetricStress& lhs,
                                               const SymmetricStress& rhs);
[[nodiscard]] SymmetricStress les_deviatoric_part(const SymmetricStress& tensor);
[[nodiscard]] Real les_strain_magnitude(const SymmetricStress& deviatoric_strain);
[[nodiscard]] Real les_van_driest_factor(Real wall_y_plus, Real a_plus = 26.0);

[[nodiscard]] LesPointEvaluation evaluate_smagorinsky(Real density,
                                                      Real filter_width,
                                                      Real coefficient,
                                                      const VelocityGradient& gradient,
                                                      Real damping_factor = 1.0);

[[nodiscard]] LesPointEvaluation evaluate_scale_similarity(
    const VelocityGradient& gradient,
    const SymmetricStress& leonard_stress,
    Real similarity_coefficient);

[[nodiscard]] LesPointEvaluation evaluate_mixed_smagorinsky_similarity(
    Real density,
    Real filter_width,
    Real smagorinsky_coefficient,
    const VelocityGradient& gradient,
    const SymmetricStress& leonard_stress,
    Real similarity_coefficient,
    Real damping_factor = 1.0);

[[nodiscard]] Real evaluate_dynamic_smagorinsky_coefficient(Real averaged_l_dot_m,
                                                            Real averaged_m_dot_m,
                                                            Real denominator_floor,
                                                            Real minimum,
                                                            Real maximum);

[[nodiscard]] LesPointEvaluation evaluate_dynamic_smagorinsky(
    Real density,
    Real filter_width,
    Real dynamic_coefficient,
    const VelocityGradient& gradient);

[[nodiscard]] LesPointEvaluation evaluate_wale(Real density,
                                               Real filter_width,
                                               Real coefficient,
                                               const VelocityGradient& gradient,
                                               Real denominator_floor = 1.0e-30);

} // namespace wcns
