#include <wcns/solver/two_equation_models.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

bool positive_finite(Real value)
{
    return std::isfinite(value) && value > 0.0;
}

void validate_context(const TurbulenceCellContext& context)
{
    if ((context.dimension != 2 && context.dimension != 3)
        || context.model_values.size() != 2 || context.model_gradients.size() != 2
        || !positive_finite(context.mean_state[temperature_density])
        || !positive_finite(context.molecular_kinematic_viscosity)
        || !positive_finite(context.reference_reynolds)
        || !positive_finite(context.reference_mach)
        || !positive_finite(context.wall_distance)) {
        throw std::invalid_argument("two-equation turbulence context is invalid");
    }
    for (const Real value : context.model_values) {
        if (!positive_finite(value)) {
            throw PhysicsError("two-equation transported value is not strictly positive");
        }
    }
    for (const auto& gradient : context.model_gradients) {
        for (const Real value : gradient) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("two-equation model gradient is non-finite");
            }
        }
    }
}

Real strain_contraction(const TurbulenceCellContext& context)
{
    const auto& du = context.primitive_gradients[0];
    const auto& dv = context.primitive_gradients[1];
    const auto& dw = context.primitive_gradients[2];
    const Real sxy = 0.5 * (du[1] + dv[0]);
    const Real sxz = 0.5 * (du[2] + dw[0]);
    const Real syz = 0.5 * (dv[2] + dw[1]);
    return 2.0 * (du[0] * du[0] + dv[1] * dv[1] + dw[2] * dw[2])
        + 4.0 * (sxy * sxy + sxz * sxz + syz * syz);
}

Real dot(const std::array<Real, 3>& left, const std::array<Real, 3>& right)
{
    return left[0] * right[0] + left[1] * right[1] + left[2] * right[2];
}

} // namespace

void SstConstants::validate() const
{
    const std::array<Real, 13> values {{sigma_k1,
                                        sigma_k2,
                                        sigma_omega1,
                                        sigma_omega2,
                                        beta1,
                                        beta2,
                                        beta_star,
                                        kappa,
                                        a1,
                                        gamma1,
                                        gamma2,
                                        production_limiter,
                                        cross_diffusion_floor}};
    for (const Real value : values) {
        if (!positive_finite(value)) throw std::invalid_argument("invalid SST constant");
    }
}

void StandardKEpsilonConstants::validate() const
{
    const std::array<Real, 6> values {
        {c_mu, c_epsilon1, c_epsilon2, sigma_k, sigma_epsilon, kappa}};
    for (const Real value : values) {
        if (!positive_finite(value)) {
            throw std::invalid_argument("invalid standard k-epsilon constant");
        }
    }
}

TwoEquationEvaluation evaluate_k_omega_sst(const TurbulenceCellContext& context,
                                            const SstConstants& constants)
{
    constants.validate();
    validate_context(context);
    const Real rho = context.mean_state[temperature_density];
    const Real k = context.model_values[0];
    const Real omega = context.model_values[1];
    const Real distance = context.wall_distance;
    const Real nu = context.molecular_kinematic_viscosity;
    const Real cross_raw = 2.0 * rho * constants.sigma_omega2 / omega
        * dot(context.model_gradients[0], context.model_gradients[1]);
    const Real cross_denominator = std::max(cross_raw, constants.cross_diffusion_floor);
    const Real wall_argument
        = 500.0 * nu / (distance * distance * omega);
    const Real k_omega_argument
        = std::sqrt(k) / (constants.beta_star * omega * distance);
    const Real arg2 = std::max(2.0 * k_omega_argument, wall_argument);
    const Real arg1 = std::min(
        std::max(k_omega_argument, wall_argument),
        4.0 * rho * constants.sigma_omega2 * k
            / (cross_denominator * distance * distance));
    TwoEquationEvaluation result;
    result.blending_f1 = std::tanh(std::pow(arg1, 4.0));
    result.blending_f2 = std::tanh(arg2 * arg2);
    const Real blend = result.blending_f1;
    const Real sigma_k = blend * constants.sigma_k1 + (1.0 - blend) * constants.sigma_k2;
    const Real sigma_omega
        = blend * constants.sigma_omega1 + (1.0 - blend) * constants.sigma_omega2;
    const Real beta = blend * constants.beta1 + (1.0 - blend) * constants.beta2;
    const Real gamma = blend * constants.gamma1 + (1.0 - blend) * constants.gamma2;
    const Real omega_scale
        = std::max(constants.a1 * omega,
                   std::sqrt(strain_contraction(context)) * result.blending_f2);
    result.eddy_kinematic_viscosity = constants.a1 * k / omega_scale;
    const Real production_unlimited
        = result.eddy_kinematic_viscosity * strain_contraction(context);
    result.production = std::min(production_unlimited,
                                 constants.production_limiter * constants.beta_star * omega * k);
    result.destruction_k = constants.beta_star * omega * k;
    result.destruction_second = beta * omega * omega;
    result.cross_diffusion = 2.0 * (1.0 - blend) * constants.sigma_omega2 / omega
        * dot(context.model_gradients[0], context.model_gradients[1]);
    result.source = {{rho * (result.production - result.destruction_k),
                      rho * (gamma * result.production
                                  / std::max(result.eddy_kinematic_viscosity, 1.0e-30)
                              - result.destruction_second + result.cross_diffusion)}};
    const Real mu_t = rho * context.reference_reynolds
        * result.eddy_kinematic_viscosity;
    result.diffusion_coefficients = {{
        rho * nu + sigma_k * mu_t / context.reference_reynolds,
        rho * nu + sigma_omega * mu_t / context.reference_reynolds,
    }};
    // The transported unknowns are rho*k and rho*omega.  At fixed density,
    // d[rho*f(q)]/d(rho*q) = df/dq, so no extra density factor belongs here.
    result.source_jacobian = {{-constants.beta_star * omega,
                               -constants.beta_star * k,
                               0.0,
                               -2.0 * beta * omega}};
    return result;
}

TwoEquationEvaluation evaluate_standard_k_epsilon(
    const TurbulenceCellContext& context,
    const StandardKEpsilonConstants& constants)
{
    constants.validate();
    validate_context(context);
    const Real rho = context.mean_state[temperature_density];
    const Real k = context.model_values[0];
    const Real epsilon = context.model_values[1];
    TwoEquationEvaluation result;
    result.blending_f1 = 1.0;
    result.blending_f2 = 1.0;
    result.eddy_kinematic_viscosity = constants.c_mu * k * k / epsilon;
    result.production = result.eddy_kinematic_viscosity * strain_contraction(context);
    result.destruction_k = epsilon;
    result.destruction_second = constants.c_epsilon2 * epsilon * epsilon / k;
    result.source = {{rho * (result.production - epsilon),
                      rho * (constants.c_epsilon1 * epsilon / k * result.production
                              - result.destruction_second)}};
    const Real mu_t = rho * context.reference_reynolds
        * result.eddy_kinematic_viscosity;
    result.diffusion_coefficients = {{
        rho * context.molecular_kinematic_viscosity
            + mu_t / (constants.sigma_k * context.reference_reynolds),
        rho * context.molecular_kinematic_viscosity
            + mu_t / (constants.sigma_epsilon * context.reference_reynolds),
    }};
    result.source_jacobian = {{-epsilon / k,
                               -1.0,
                               0.0,
                               -2.0 * constants.c_epsilon2 * epsilon / k}};
    return result;
}

std::array<Real, 2> two_equation_farfield_values(const TurbulenceModelConfig& config)
{
    config.validate();
    if (config.kind != TurbulenceModelKind::KOmegaSst
        && config.kind != TurbulenceModelKind::KEpsilon) {
        throw std::invalid_argument("two-equation farfield requires SST or k-epsilon");
    }
    const Real k = 1.5 * config.freestream_turbulence_intensity
        * config.freestream_turbulence_intensity;
    if (config.kind == TurbulenceModelKind::KOmegaSst) {
        const SstConstants constants;
        const Real omega = std::sqrt(k)
            / (std::pow(constants.beta_star, 0.25) * config.freestream_length_scale);
        return {{k, omega}};
    }
    const StandardKEpsilonConstants constants;
    const Real epsilon = std::pow(constants.c_mu, 0.75) * std::pow(k, 1.5)
        / config.freestream_length_scale;
    return {{k, epsilon}};
}

} // namespace wcns
