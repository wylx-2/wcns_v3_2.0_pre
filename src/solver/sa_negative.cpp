#include <wcns/solver/turbulence_model.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

struct Dual {
    Real value = 0.0;
    Real derivative = 0.0;
};

Dual operator+(Dual lhs, Dual rhs)
{
    return {lhs.value + rhs.value, lhs.derivative + rhs.derivative};
}

Dual operator-(Dual lhs, Dual rhs)
{
    return {lhs.value - rhs.value, lhs.derivative - rhs.derivative};
}

Dual operator*(Dual lhs, Dual rhs)
{
    return {lhs.value * rhs.value,
            lhs.derivative * rhs.value + lhs.value * rhs.derivative};
}

Dual operator/(Dual lhs, Dual rhs)
{
    const Real denominator = rhs.value * rhs.value;
    return {lhs.value / rhs.value,
            (lhs.derivative * rhs.value - lhs.value * rhs.derivative) / denominator};
}

Dual operator-(Dual value) { return {-value.value, -value.derivative}; }

Dual operator+(Dual lhs, Real rhs) { return lhs + Dual {rhs, 0.0}; }
Dual operator+(Real lhs, Dual rhs) { return Dual {lhs, 0.0} + rhs; }
Dual operator-(Real lhs, Dual rhs) { return Dual {lhs, 0.0} - rhs; }
Dual operator*(Dual lhs, Real rhs) { return lhs * Dual {rhs, 0.0}; }
Dual operator*(Real lhs, Dual rhs) { return Dual {lhs, 0.0} * rhs; }
Dual operator/(Dual lhs, Real rhs) { return lhs / Dual {rhs, 0.0}; }
Dual operator/(Real lhs, Dual rhs) { return Dual {lhs, 0.0} / rhs; }

Dual square(Dual value) { return value * value; }

Dual cube(Dual value) { return value * value * value; }

Dual sixth(Dual value)
{
    const Dual value3 = cube(value);
    return value3 * value3;
}

Dual exponential(Dual value)
{
    const Real result = std::exp(value.value);
    return {result, result * value.derivative};
}

Dual sixth_root(Dual value)
{
    const Real result = std::pow(value.value, 1.0 / 6.0);
    return {result, value.derivative / (6.0 * std::pow(value.value, 5.0 / 6.0))};
}

bool positive_finite(Real value)
{
    return std::isfinite(value) && value > 0.0;
}

void validate_context(const TurbulenceCellContext& context)
{
    if ((context.dimension != 2 && context.dimension != 3)
        || context.model_values.size() != 1 || context.model_gradients.size() != 1
        || !positive_finite(context.mean_state[temperature_density])
        || !positive_finite(context.molecular_kinematic_viscosity)
        || !positive_finite(context.wall_distance)) {
        throw std::invalid_argument("SA-neg cell context is incomplete or invalid");
    }
    for (const auto& gradient : context.primitive_gradients) {
        for (const Real value : gradient) {
            if (!std::isfinite(value)) {
                throw std::invalid_argument("SA-neg velocity gradient is non-finite");
            }
        }
    }
    for (const Real value : context.model_gradients.front()) {
        if (!std::isfinite(value)) {
            throw std::invalid_argument("SA-neg working-variable gradient is non-finite");
        }
    }
    if (!std::isfinite(context.model_values.front())) {
        throw std::invalid_argument("SA-neg working variable is non-finite");
    }
}

Real vorticity_magnitude(const TurbulenceCellContext& context)
{
    const auto& du = context.primitive_gradients[0];
    const auto& dv = context.primitive_gradients[1];
    const auto& dw = context.primitive_gradients[2];
    const Real wz = dv[0] - du[1];
    if (context.dimension == 2) return std::abs(wz);
    const Real wy = du[2] - dw[0];
    const Real wx = dw[1] - dv[2];
    return std::sqrt(wx * wx + wy * wy + wz * wz);
}

} // namespace

Real SaNegativeConstants::cw1() const noexcept
{
    return cb1 / (kappa * kappa) + (1.0 + cb2) / sigma;
}

void SaNegativeConstants::validate() const
{
    const Real values[] = {cb1, sigma, cb2, kappa, cw2, cw3, cv1, ct3, ct4, cn1, c2, c3};
    if (!std::all_of(std::begin(values), std::end(values), positive_finite)
        || !(c2 < c3) || !positive_finite(cw1())) {
        throw std::invalid_argument("SA-neg constants are invalid");
    }
}

SaNegativeEvaluation evaluate_sa_negative(const TurbulenceCellContext& context,
                                          const SaNegativeConstants& constants)
{
    constants.validate();
    validate_context(context);
    const Dual nu_tilde {context.model_values.front(), 1.0};
    const Real nu = context.molecular_kinematic_viscosity;
    const Real distance = context.wall_distance;
    const Real omega = vorticity_magnitude(context);
    const Dual chi = nu_tilde / nu;
    const Dual chi3 = cube(chi);

    SaNegativeEvaluation result;
    result.chi = chi.value;
    result.vorticity_magnitude = omega;
    result.negative_branch = nu_tilde.value < 0.0;

    Dual diffusion;
    Dual production;
    Dual destruction;
    if (!result.negative_branch) {
        const Dual fv1
            = chi3 / (chi3 + constants.cv1 * constants.cv1 * constants.cv1);
        const Dual ft2 = constants.ct3 * exponential(-constants.ct4 * square(chi));
        const Dual fv2 = 1.0 - chi / (1.0 + chi * fv1);
        const Dual sbar = nu_tilde * fv2
            / (constants.kappa * constants.kappa * distance * distance);
        Dual modified;
        if (sbar.value >= -constants.c2 * omega) {
            modified = omega + sbar;
        } else if (omega == 0.0) {
            modified = {0.0, 0.0};
        } else {
            modified = omega
                + omega * (constants.c2 * constants.c2 * omega + constants.c3 * sbar)
                    / ((constants.c3 - 2.0 * constants.c2) * omega - sbar);
        }
        Dual r {10.0, 0.0};
        if (modified.value != 0.0) {
            const Dual raw = nu_tilde
                / (modified * (constants.kappa * constants.kappa * distance * distance));
            if (raw.value < 10.0) r = raw;
        }
        const Dual g = r + constants.cw2 * (sixth(r) - r);
        const Dual fw = g
            * sixth_root((1.0 + std::pow(constants.cw3, 6.0))
                         / (sixth(g) + std::pow(constants.cw3, 6.0)));
        production = constants.cb1 * (1.0 - ft2) * modified * nu_tilde;
        destruction = -(constants.cw1() * fw
                         - constants.cb1 / (constants.kappa * constants.kappa) * ft2)
            * square(nu_tilde / distance);
        diffusion = (nu + nu_tilde) / constants.sigma;
        result.fv1 = fv1.value;
        result.ft2 = ft2.value;
        result.fv2 = fv2.value;
        result.modified_vorticity = modified.value;
        result.r = r.value;
        result.fw = fw.value;
        result.fn = 1.0;
        result.eddy_kinematic_viscosity = nu_tilde.value * fv1.value;
    } else {
        const Dual fn = (constants.cn1 + chi3) / (constants.cn1 - chi3);
        production = constants.cb1 * (1.0 - constants.ct3) * omega * nu_tilde;
        destruction = constants.cw1() * square(nu_tilde / distance);
        diffusion = (nu + nu_tilde * fn) / constants.sigma;
        result.fn = fn.value;
        result.fv1 = 0.0;
        result.ft2 = 0.0;
        result.modified_vorticity = omega;
        result.eddy_kinematic_viscosity = 0.0;
    }

    Real gradient_squared = 0.0;
    for (int direction = 0; direction < context.dimension; ++direction) {
        const Real value = context.model_gradients.front()[static_cast<std::size_t>(direction)];
        gradient_squared += value * value;
    }
    const Real cross = constants.cb2 * gradient_squared / constants.sigma;
    const Dual source = production + destruction + cross;
    const Real finite_values[] = {result.chi,
                                  result.fv1,
                                  result.fv2,
                                  result.ft2,
                                  result.fn,
                                  result.vorticity_magnitude,
                                  result.modified_vorticity,
                                  result.r,
                                  result.fw,
                                  diffusion.value,
                                  result.eddy_kinematic_viscosity,
                                  production.value,
                                  destruction.value,
                                  cross,
                                  source.value,
                                  source.derivative};
    if (!std::all_of(std::begin(finite_values), std::end(finite_values), [](Real value) {
            return std::isfinite(value);
        })
        || !positive_finite(diffusion.value) || result.eddy_kinematic_viscosity < 0.0) {
        throw PhysicsError("SA-neg evaluation produced an inadmissible result");
    }
    result.diffusion_coefficient = diffusion.value;
    result.production_source = production.value;
    result.destruction_source = destruction.value;
    result.cross_diffusion_source = cross;
    result.source = source.value;
    result.source_derivative = source.derivative;
    return result;
}

Real sa_negative_farfield_value(const TurbulenceModelConfig& config, Real reference_reynolds)
{
    config.validate();
    if (config.kind != TurbulenceModelKind::SaNegative || !positive_finite(reference_reynolds)) {
        throw std::invalid_argument("SA-neg farfield value requires SA-neg and positive Reynolds number");
    }
    return config.sa_farfield_nu_tilde_ratio / reference_reynolds;
}

Real local_implicit_turbulence_increment(Real explicit_residual,
                                         Real source_jacobian,
                                         Real time_step)
{
    if (!std::isfinite(explicit_residual) || !std::isfinite(source_jacobian)
        || !positive_finite(time_step)) {
        throw std::invalid_argument("local implicit turbulence update inputs are invalid");
    }
    const Real denominator = 1.0 - time_step * source_jacobian;
    if (!positive_finite(denominator)) {
        throw PhysicsError("local implicit turbulence source diagonal is not positive");
    }
    return time_step * explicit_residual / denominator;
}

} // namespace wcns
