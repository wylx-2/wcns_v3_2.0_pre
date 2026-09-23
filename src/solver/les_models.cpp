#include <wcns/solver/les_models.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

void validate_input(Real density, Real filter_width, const VelocityGradient& gradient)
{
    if (!std::isfinite(density) || density <= 0.0 || !std::isfinite(filter_width)
        || filter_width <= 0.0) {
        throw std::invalid_argument("LES density and filter width must be positive and finite");
    }
    for (const auto& row : gradient) {
        for (const Real value : row) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("LES velocity gradient contains a non-finite value");
            }
        }
    }
}

SymmetricStress scaled(const SymmetricStress& value, Real factor)
{
    return {factor * value.xx,
            factor * value.yy,
            factor * value.zz,
            factor * value.xy,
            factor * value.xz,
            factor * value.yz};
}

SymmetricStress added(const SymmetricStress& lhs, const SymmetricStress& rhs)
{
    return {lhs.xx + rhs.xx,
            lhs.yy + rhs.yy,
            lhs.zz + rhs.zz,
            lhs.xy + rhs.xy,
            lhs.xz + rhs.xz,
            lhs.yz + rhs.yz};
}

LesPointEvaluation finish(const VelocityGradient& gradient,
                          SymmetricStress stress,
                          Real eddy_viscosity)
{
    LesPointEvaluation result;
    result.deviatoric_strain = les_deviatoric_strain(gradient);
    result.strain_magnitude = les_strain_magnitude(result.deviatoric_strain);
    result.physical_sgs_stress = stress;
    result.eddy_viscosity = eddy_viscosity;
    result.energy_transfer
        = -les_symmetric_inner_product(stress, result.deviatoric_strain);
    if (!std::isfinite(result.energy_transfer) || !std::isfinite(eddy_viscosity)) {
        throw std::invalid_argument("LES evaluation produced a non-finite value");
    }
    return result;
}

} // namespace

SymmetricStress les_deviatoric_strain(const VelocityGradient& gradient)
{
    for (const auto& row : gradient) {
        for (const Real value : row) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("LES velocity gradient contains a non-finite value");
            }
        }
    }
    const Real divergence = gradient[0][0] + gradient[1][1] + gradient[2][2];
    return {gradient[0][0] - divergence / 3.0,
            gradient[1][1] - divergence / 3.0,
            gradient[2][2] - divergence / 3.0,
            0.5 * (gradient[0][1] + gradient[1][0]),
            0.5 * (gradient[0][2] + gradient[2][0]),
            0.5 * (gradient[1][2] + gradient[2][1])};
}

Real les_symmetric_inner_product(const SymmetricStress& lhs, const SymmetricStress& rhs)
{
    return lhs.xx * rhs.xx + lhs.yy * rhs.yy + lhs.zz * rhs.zz
        + 2.0 * (lhs.xy * rhs.xy + lhs.xz * rhs.xz + lhs.yz * rhs.yz);
}

SymmetricStress les_deviatoric_part(const SymmetricStress& tensor)
{
    const Real trace_third = (tensor.xx + tensor.yy + tensor.zz) / 3.0;
    return {tensor.xx - trace_third,
            tensor.yy - trace_third,
            tensor.zz - trace_third,
            tensor.xy,
            tensor.xz,
            tensor.yz};
}

Real les_strain_magnitude(const SymmetricStress& deviatoric_strain)
{
    const Real square = 2.0
        * les_symmetric_inner_product(deviatoric_strain, deviatoric_strain);
    return std::sqrt(std::max(Real {0.0}, square));
}

Real les_van_driest_factor(Real wall_y_plus, Real a_plus)
{
    if (!std::isfinite(wall_y_plus) || wall_y_plus < 0.0 || !std::isfinite(a_plus)
        || a_plus <= 0.0) {
        throw std::invalid_argument("van-Driest damping input is invalid");
    }
    const Real factor = 1.0 - std::exp(-wall_y_plus / a_plus);
    return factor * factor;
}

LesPointEvaluation evaluate_smagorinsky(Real density,
                                        Real filter_width,
                                        Real coefficient,
                                        const VelocityGradient& gradient,
                                        Real damping_factor)
{
    validate_input(density, filter_width, gradient);
    if (!std::isfinite(coefficient) || coefficient < 0.0 || coefficient > 0.3
        || !std::isfinite(damping_factor) || damping_factor < 0.0
        || damping_factor > 1.0) {
        throw std::invalid_argument("Smagorinsky coefficient or damping is invalid");
    }
    const auto strain = les_deviatoric_strain(gradient);
    const Real magnitude = les_strain_magnitude(strain);
    const Real length = coefficient * filter_width;
    const Real viscosity = density * length * length * magnitude * damping_factor;
    return finish(gradient, scaled(strain, -2.0 * viscosity), viscosity);
}

LesPointEvaluation evaluate_scale_similarity(const VelocityGradient& gradient,
                                             const SymmetricStress& leonard_stress,
                                             Real similarity_coefficient)
{
    if (!std::isfinite(similarity_coefficient)) {
        throw std::invalid_argument("scale-similarity coefficient is non-finite");
    }
    const auto deviatoric = les_deviatoric_part(leonard_stress);
    return finish(gradient, scaled(deviatoric, similarity_coefficient), 0.0);
}

LesPointEvaluation evaluate_mixed_smagorinsky_similarity(
    Real density,
    Real filter_width,
    Real smagorinsky_coefficient,
    const VelocityGradient& gradient,
    const SymmetricStress& leonard_stress,
    Real similarity_coefficient,
    Real damping_factor)
{
    const auto eddy = evaluate_smagorinsky(
        density, filter_width, smagorinsky_coefficient, gradient, damping_factor);
    const auto similarity
        = evaluate_scale_similarity(gradient, leonard_stress, similarity_coefficient);
    return finish(gradient,
                  added(eddy.physical_sgs_stress, similarity.physical_sgs_stress),
                  eddy.eddy_viscosity);
}

Real evaluate_dynamic_smagorinsky_coefficient(Real averaged_l_dot_m,
                                              Real averaged_m_dot_m,
                                              Real denominator_floor,
                                              Real minimum,
                                              Real maximum)
{
    if (!std::isfinite(averaged_l_dot_m) || !std::isfinite(averaged_m_dot_m)
        || averaged_m_dot_m < 0.0 || !std::isfinite(denominator_floor)
        || denominator_floor <= 0.0 || !std::isfinite(minimum)
        || !std::isfinite(maximum) || minimum > maximum) {
        throw std::invalid_argument("dynamic Smagorinsky least-squares input is invalid");
    }
    return std::clamp(averaged_l_dot_m / (averaged_m_dot_m + denominator_floor),
                      minimum,
                      maximum);
}

LesPointEvaluation evaluate_dynamic_smagorinsky(Real density,
                                                Real filter_width,
                                                Real dynamic_coefficient,
                                                const VelocityGradient& gradient)
{
    validate_input(density, filter_width, gradient);
    if (!std::isfinite(dynamic_coefficient)) {
        throw std::invalid_argument("dynamic Smagorinsky coefficient is non-finite");
    }
    const auto strain = les_deviatoric_strain(gradient);
    const Real magnitude = les_strain_magnitude(strain);
    const Real viscosity
        = density * dynamic_coefficient * filter_width * filter_width * magnitude;
    return finish(gradient, scaled(strain, -2.0 * viscosity), viscosity);
}

LesPointEvaluation evaluate_wale(Real density,
                                 Real filter_width,
                                 Real coefficient,
                                 const VelocityGradient& gradient,
                                 Real denominator_floor)
{
    validate_input(density, filter_width, gradient);
    if (!std::isfinite(coefficient) || coefficient <= 0.0 || coefficient > 1.0
        || !std::isfinite(denominator_floor) || denominator_floor <= 0.0) {
        throw std::invalid_argument("WALE coefficient or denominator floor is invalid");
    }
    VelocityGradient square {};
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            for (int inner = 0; inner < 3; ++inner) {
                square[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
                    += gradient[static_cast<std::size_t>(row)][static_cast<std::size_t>(inner)]
                    * gradient[static_cast<std::size_t>(inner)]
                              [static_cast<std::size_t>(column)];
            }
        }
    }
    const Real trace_square = square[0][0] + square[1][1] + square[2][2];
    SymmetricStress sd {
        square[0][0] - trace_square / 3.0,
        square[1][1] - trace_square / 3.0,
        square[2][2] - trace_square / 3.0,
        0.5 * (square[0][1] + square[1][0]),
        0.5 * (square[0][2] + square[2][0]),
        0.5 * (square[1][2] + square[2][1]),
    };
    const SymmetricStress full_strain {
        gradient[0][0],
        gradient[1][1],
        gradient[2][2],
        0.5 * (gradient[0][1] + gradient[1][0]),
        0.5 * (gradient[0][2] + gradient[2][0]),
        0.5 * (gradient[1][2] + gradient[2][1]),
    };
    const Real strain_square = les_symmetric_inner_product(full_strain, full_strain);
    const Real sd_square = std::max(
        Real {0.0}, les_symmetric_inner_product(sd, sd));
    const Real numerator = std::pow(sd_square, 1.5);
    const Real denominator = std::pow(std::max(Real {0.0}, strain_square), 2.5)
        + std::pow(sd_square, 1.25) + denominator_floor;
    const Real length = coefficient * filter_width;
    const Real viscosity = density * length * length * numerator / denominator;
    const auto deviatoric_strain = les_deviatoric_strain(gradient);
    return finish(
        gradient, scaled(deviatoric_strain, -2.0 * viscosity), viscosity);
}

} // namespace wcns
