#include <wcns/solver/riemann_solver.hpp>
#include <wcns/solver/wcns_reconstruction.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace wcns {
namespace {

bool valid_algorithm_name(std::string_view name)
{
    if (name.empty()) return false;
    for (const char character : name) {
        const bool lower = character >= 'a' && character <= 'z';
        const bool digit = character >= '0' && character <= '9';
        if (!lower && !digit && character != '_') return false;
    }
    return true;
}

Normal3 checked_unit_normal(Normal3 normal)
{
    const Real norm = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1.0e-12) {
        throw PhysicsError("Riemann solver requires a unit normal");
    }
    return normal;
}

Real normal_velocity(const PressurePrimitiveState& state, Normal3 normal)
{
    return state[1] * normal.x + state[2] * normal.y + state[3] * normal.z;
}

bool finite_state(const ConservativeState& state)
{
    return std::all_of(state.begin(), state.end(), [](Real value) { return std::isfinite(value); });
}

Matrix5 shifted_matrix(const Matrix5& matrix, Real shift)
{
    Matrix5 result = matrix;
    for (int index = 0; index < euler_components; ++index) {
        result[static_cast<std::size_t>(index)][static_cast<std::size_t>(index)] -= shift;
    }
    return result;
}

void add_scaled(Matrix5& destination, const Matrix5& source, Real scale)
{
    for (int row = 0; row < euler_components; ++row) {
        for (int column = 0; column < euler_components; ++column) {
            destination[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
                += scale
                * source[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)];
        }
    }
}

ConservativeState matrix_vector_product(const Matrix5& matrix,
                                        const ConservativeState& vector)
{
    ConservativeState result {};
    for (int row = 0; row < euler_components; ++row) {
        for (int column = 0; column < euler_components; ++column) {
            result[static_cast<std::size_t>(row)]
                += matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
                * vector[static_cast<std::size_t>(column)];
        }
    }
    return result;
}

Real entropy_fixed_magnitude(Real eigenvalue, Real delta)
{
    Real magnitude = std::abs(eigenvalue);
    if (magnitude < delta) {
        magnitude = 0.5 * (magnitude * magnitude / delta + delta);
    }
    return magnitude;
}

Matrix5 roe_flux_jacobian(const EulerCharacteristicBasis& basis,
                          const std::array<Real, euler_components>& eigenvalues)
{
    Matrix5 result {};
    for (int column = 0; column < euler_components; ++column) {
        ConservativeState unit {};
        unit[static_cast<std::size_t>(column)] = 1.0;
        auto characteristic = project_characteristic(unit, basis);
        for (int wave = 0; wave < euler_components; ++wave) {
            characteristic[static_cast<std::size_t>(wave)]
                *= eigenvalues[static_cast<std::size_t>(wave)];
        }
        const auto physical = restore_characteristic(characteristic, basis);
        for (int row = 0; row < euler_components; ++row) {
            result[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
                = physical[static_cast<std::size_t>(row)];
        }
    }
    return result;
}

// A_p has the three distinct eigenvalues (u_n, lambda_-, lambda_+).  Its
// entropy-fixed absolute value is therefore the quadratic interpolation of
// |lambda|_delta evaluated at A_p.  This avoids importing ordinary Roe
// eigenvectors into the preconditioned system.
Matrix5 absolute_preconditioned_jacobian(
    const Matrix5& matrix,
    const std::array<Real, 3>& eigenvalues,
    const std::array<Real, 3>& magnitudes)
{
    Matrix5 result {};
    for (int index = 0; index < 3; ++index) {
        const int first = (index + 1) % 3;
        const int second = (index + 2) % 3;
        const Real denominator = (eigenvalues[static_cast<std::size_t>(index)]
                                  - eigenvalues[static_cast<std::size_t>(first)])
            * (eigenvalues[static_cast<std::size_t>(index)]
               - eigenvalues[static_cast<std::size_t>(second)]);
        const Real scale = std::max(
            {1.0,
             std::abs(eigenvalues[static_cast<std::size_t>(index)]),
             std::abs(eigenvalues[static_cast<std::size_t>(first)]),
             std::abs(eigenvalues[static_cast<std::size_t>(second)])});
        if (!std::isfinite(denominator)
            || std::abs(denominator) <= std::numeric_limits<Real>::epsilon() * scale * scale) {
            throw PhysicsError("Weiss-Smith eigenvalues are numerically degenerate");
        }
        const auto polynomial
            = multiply(shifted_matrix(matrix, eigenvalues[static_cast<std::size_t>(first)]),
                       shifted_matrix(matrix, eigenvalues[static_cast<std::size_t>(second)]));
        add_scaled(result,
                   polynomial,
                   magnitudes[static_cast<std::size_t>(index)] / denominator);
    }
    return result;
}

struct RoeAverage {
    Normal3 velocity {};
    Real normal_velocity = 0.0;
    Real enthalpy = 0.0;
    Real sound_speed = 0.0;
};

RoeAverage roe_average(const PressurePrimitiveState& left,
                       const PressurePrimitiveState& right,
                       Normal3 normal,
                       const IdealGas& gas)
{
    const auto left_conservative = to_conservative(left, gas);
    const auto right_conservative = to_conservative(right, gas);
    const Real root_left = std::sqrt(left[0]);
    const Real root_right = std::sqrt(right[0]);
    const Real denominator = root_left + root_right;
    if (!std::isfinite(denominator) || denominator <= 0.0) {
        throw PhysicsError("Roe average has an invalid density denominator");
    }
    const auto average = [&](Real left_value, Real right_value) {
        return (root_left * left_value + root_right * right_value) / denominator;
    };
    RoeAverage result;
    result.velocity = {
        average(left[1], right[1]),
        average(left[2], right[2]),
        average(left[3], right[3]),
    };
    const Real left_enthalpy = (left_conservative[4] + left[4]) / left[0];
    const Real right_enthalpy = (right_conservative[4] + right[4]) / right[0];
    result.enthalpy = average(left_enthalpy, right_enthalpy);
    result.normal_velocity = result.velocity.x * normal.x + result.velocity.y * normal.y
        + result.velocity.z * normal.z;
    const Real speed_squared = result.velocity.x * result.velocity.x
        + result.velocity.y * result.velocity.y + result.velocity.z * result.velocity.z;
    const Real sound_squared = (gas.gamma - 1.0) * (result.enthalpy - 0.5 * speed_squared);
    if (!std::isfinite(sound_squared) || sound_squared <= 0.0) {
        throw PhysicsError("Roe average has an invalid sound speed");
    }
    result.sound_speed = std::sqrt(sound_squared);
    return result;
}

RiemannResult rusanov_result(const PressurePrimitiveState& left,
                             const PressurePrimitiveState& right,
                             Normal3 normal,
                             const GasModel& gas,
                             const NumericalFloors& floors,
                             std::string requested,
                             RiemannFallbackReason reason)
{
    const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
    const Real spectral_radius
        = std::max(std::abs(normal_velocity(left, normal)) + sound_speed(left, ideal),
                   std::abs(normal_velocity(right, normal)) + sound_speed(right, ideal));
    RiemannResult result {
        rusanov_flux(left, right, normal, ideal),
        spectral_radius,
        requested,
        "rusanov",
        reason,
        {},
    };
    if (requested != "rusanov") {
        result.fallback_path.push_back({std::move(requested), "rusanov", reason});
    }
    return result;
}

Real state_speed(const PressurePrimitiveState& state)
{
    return std::sqrt(state[1] * state[1] + state[2] * state[2] + state[3] * state[3]);
}

RiemannResult preconditioned_rusanov_result(const PressurePrimitiveState& left,
                                            const PressurePrimitiveState& right,
                                            Normal3 normal,
                                            const GasModel& gas,
                                            const NumericalFloors& floors,
                                            Real viscous_speed,
                                            const RiemannSolverParameters& parameters,
                                            RiemannFallbackReason reason)
{
    const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
    const auto left_conservative = to_conservative(left, ideal);
    const auto right_conservative = to_conservative(right, ideal);
    const auto left_flux = euler_flux(left, normal, ideal);
    const auto right_flux = euler_flux(right, normal, ideal);
    PressurePrimitiveState average {};
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        average[index] = 0.5 * (left[index] + right[index]);
    }
    const auto preconditioned = weiss_smith_pressure_state(
        average, normal, gas, floors, parameters.preconditioner, viscous_speed, 3);
    auto full_speed_parameters = parameters.preconditioner;
    full_speed_parameters.mach_cutoff = 1.0;
    full_speed_parameters.viscous_cutoff = 0.0;
    const auto physical = weiss_smith_pressure_state(
        average, normal, gas, floors, full_speed_parameters, 0.0, 3);
    Real spectral_radius = 0.0;
    for (const Real eigenvalue : preconditioned.eigenvalues) {
        spectral_radius = std::max(spectral_radius, std::abs(eigenvalue));
    }
    ConservativeState jump {};
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        jump[index] = right_conservative[index] - left_conservative[index];
    }
    auto primitive_jump = matrix_vector_product(physical.inverse, jump);
    for (Real& value : primitive_jump) value *= spectral_radius;
    const auto dissipation
        = matrix_vector_product(preconditioned.gamma, primitive_jump);
    ConservativeState flux {};
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        flux[index] = 0.5 * (left_flux[index] + right_flux[index])
            - 0.5 * dissipation[index];
    }
    return {flux,
            spectral_radius,
            "roe",
            "preconditioned_rusanov",
            reason,
            {{"roe", "preconditioned_rusanov", reason}}};
}

RiemannResult preconditioned_roe_result(const PressurePrimitiveState& left,
                                        const PressurePrimitiveState& right,
                                        Normal3 normal,
                                        const GasModel& gas,
                                        const NumericalFloors& floors,
                                        Real viscous_speed,
                                        const RiemannSolverParameters& parameters)
{
    const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
    const auto fallback = [&](RiemannFallbackReason reason) {
        return preconditioned_rusanov_result(
            left, right, normal, gas, floors, viscous_speed, parameters, reason);
    };
    try {
        const auto left_conservative = to_conservative(left, ideal);
        const auto right_conservative = to_conservative(right, ideal);
        const auto left_flux = euler_flux(left, normal, ideal);
        const auto right_flux = euler_flux(right, normal, ideal);
        const auto average = roe_average(left, right, normal, ideal);
        const Real left_sound = sound_speed(left, ideal);
        const Real right_sound = sound_speed(right, ideal);
        const Real left_reference = weiss_smith_reference_speed(
            state_speed(left), left_sound, viscous_speed, parameters.preconditioner);
        const Real right_reference = weiss_smith_reference_speed(
            state_speed(right), right_sound, viscous_speed, parameters.preconditioner);
        const auto basis = make_roe_characteristic_basis(left, right, normal, gas, floors, 3);
        ConservativeState jump {};
        for (int component = 0; component < euler_components; ++component) {
            const auto index = static_cast<std::size_t>(component);
            jump[index] = right_conservative[index] - left_conservative[index];
        }
        const Real roe_density = std::sqrt(left[0] * right[0]);
        const Real roe_pressure
            = roe_density * average.sound_speed * average.sound_speed / gas.gamma();
        const PressurePrimitiveState roe_state {{roe_density,
                                                 average.velocity.x,
                                                 average.velocity.y,
                                                 average.velocity.z,
                                                 roe_pressure}};
        const auto preconditioned_state = weiss_smith_pressure_state(
            roe_state, normal, gas, floors, parameters.preconditioner, viscous_speed, 3);
        const Real beta = preconditioned_state.beta;
        const Real delta = parameters.entropy_fix_coefficient
            * std::max({left_reference, right_reference, preconditioned_state.reference_speed});
        ConservativeState dissipation {};
        Real spectral_radius = 0.0;
        if (beta >= 1.0 - 16.0 * std::numeric_limits<Real>::epsilon()) {
            const std::array<Real, euler_components> eigenvalues {{
                average.normal_velocity - average.sound_speed,
                average.normal_velocity,
                average.normal_velocity,
                average.normal_velocity,
                average.normal_velocity + average.sound_speed,
            }};
            const auto strengths = project_characteristic(jump, basis);
            ConservativeState scaled_strengths {};
            for (int wave = 0; wave < euler_components; ++wave) {
                const Real magnitude = entropy_fixed_magnitude(
                    eigenvalues[static_cast<std::size_t>(wave)], delta);
                spectral_radius = std::max(spectral_radius, magnitude);
                scaled_strengths[static_cast<std::size_t>(wave)]
                    = magnitude * strengths[static_cast<std::size_t>(wave)];
            }
            dissipation = restore_characteristic(scaled_strengths, basis);
        } else {
            const auto preconditioned
                = weiss_smith_eigenvalues(average.normal_velocity, average.sound_speed, beta);
            const std::array<Real, 3> distinct {{
                preconditioned[0], preconditioned[3], preconditioned[4]}};
            std::array<Real, 3> magnitudes {};
            for (int wave = 0; wave < 3; ++wave) {
                magnitudes[static_cast<std::size_t>(wave)] = entropy_fixed_magnitude(
                    distinct[static_cast<std::size_t>(wave)], delta);
                spectral_radius
                    = std::max(spectral_radius, magnitudes[static_cast<std::size_t>(wave)]);
            }

            auto full_speed_parameters = parameters.preconditioner;
            full_speed_parameters.mach_cutoff = 1.0;
            full_speed_parameters.viscous_cutoff = 0.0;
            const auto physical_state = weiss_smith_pressure_state(
                roe_state, normal, gas, floors, full_speed_parameters, 0.0, 3);

            const std::array<Real, euler_components> physical_eigenvalues {{
                average.normal_velocity - average.sound_speed,
                average.normal_velocity,
                average.normal_velocity,
                average.normal_velocity,
                average.normal_velocity + average.sound_speed,
            }};
            const auto physical_jacobian = roe_flux_jacobian(basis, physical_eigenvalues);
            const auto flux_primitive = multiply(physical_jacobian, physical_state.gamma);
            const auto preconditioned_jacobian
                = multiply(preconditioned_state.inverse, flux_primitive);
            const auto absolute_jacobian = absolute_preconditioned_jacobian(
                preconditioned_jacobian, distinct, magnitudes);
            const auto primitive_jump
                = matrix_vector_product(physical_state.inverse, jump);
            dissipation = matrix_vector_product(
                preconditioned_state.gamma,
                matrix_vector_product(absolute_jacobian, primitive_jump));
        }
        ConservativeState flux {};
        for (int component = 0; component < euler_components; ++component) {
            const auto index = static_cast<std::size_t>(component);
            flux[index]
                = 0.5 * (left_flux[index] + right_flux[index]) - 0.5 * dissipation[index];
        }
        if (!finite_state(flux)) return fallback(RiemannFallbackReason::NonFiniteFlux);
        return {flux,
                spectral_radius,
                "roe",
                "roe",
                RiemannFallbackReason::None,
                {}};
    } catch (const PhysicsError&) {
        return fallback(RiemannFallbackReason::InvalidRoeAverage);
    }
}

class RusanovStrategy final : public IRiemannSolver {
public:
    [[nodiscard]] std::string_view name() const noexcept override { return "rusanov"; }

    [[nodiscard]] RiemannResult solve(const PressurePrimitiveState& left,
                                      const PressurePrimitiveState& right,
                                      Normal3 unit_normal,
                                      const GasModel& gas,
                                      const NumericalFloors& floors) const override
    {
        const auto normal = checked_unit_normal(unit_normal);
        floors.validate();
        return rusanov_result(
            left, right, normal, gas, floors, "rusanov", RiemannFallbackReason::None);
    }
};

class HllcStrategy final : public IRiemannSolver {
public:
    explicit HllcStrategy(RiemannSolverParameters parameters)
        : parameters_(parameters)
    {
        parameters_.validate();
    }

    [[nodiscard]] std::string_view name() const noexcept override { return "hllc"; }

    [[nodiscard]] RiemannResult solve(const PressurePrimitiveState& left,
                                      const PressurePrimitiveState& right,
                                      Normal3 unit_normal,
                                      const GasModel& gas,
                                      const NumericalFloors& floors) const override
    {
        const auto normal = checked_unit_normal(unit_normal);
        floors.validate();
        const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
        const auto left_conservative = to_conservative(left, ideal);
        const auto right_conservative = to_conservative(right, ideal);
        const auto left_flux = euler_flux(left, normal, ideal);
        const auto right_flux = euler_flux(right, normal, ideal);
        const Real un_left = normal_velocity(left, normal);
        const Real un_right = normal_velocity(right, normal);
        const Real sound_left = sound_speed(left, ideal);
        const Real sound_right = sound_speed(right, ideal);
        RoeAverage average;
        try {
            average = roe_average(left, right, normal, ideal);
        } catch (const PhysicsError&) {
            return rusanov_result(
                left, right, normal, gas, floors, "hllc", RiemannFallbackReason::InvalidWaveSpeed);
        }
        const Real speed_left
            = std::min(un_left - sound_left, average.normal_velocity - average.sound_speed);
        const Real speed_right
            = std::max(un_right + sound_right, average.normal_velocity + average.sound_speed);
        const Real spectral_radius = std::max(std::abs(speed_left), std::abs(speed_right));
        if (!std::isfinite(speed_left) || !std::isfinite(speed_right)
            || speed_left >= speed_right) {
            return rusanov_result(
                left, right, normal, gas, floors, "hllc", RiemannFallbackReason::InvalidWaveSpeed);
        }
        if (speed_left >= 0.0) {
            return {left_flux, spectral_radius, "hllc", "hllc", RiemannFallbackReason::None, {}};
        }
        if (speed_right <= 0.0) {
            return {right_flux, spectral_radius, "hllc", "hllc", RiemannFallbackReason::None, {}};
        }
        const Real left_term = left[0] * (speed_left - un_left);
        const Real right_term = right[0] * (speed_right - un_right);
        const Real denominator = left_term - right_term;
        const Real denominator_scale = std::max({1.0, std::abs(left_term), std::abs(right_term)});
        if (!std::isfinite(denominator)
            || std::abs(denominator) <= parameters_.denominator_tolerance * denominator_scale) {
            return rusanov_result(
                left, right, normal, gas, floors, "hllc", RiemannFallbackReason::InvalidWaveSpeed);
        }
        const Real speed_middle
            = (right[4] - left[4] + left_term * un_left - right_term * un_right) / denominator;
        if (!std::isfinite(speed_middle) || speed_middle <= speed_left
            || speed_middle >= speed_right) {
            return rusanov_result(
                left, right, normal, gas, floors, "hllc", RiemannFallbackReason::InvalidWaveSpeed);
        }
        const auto star_flux = [&](const PressurePrimitiveState& state,
                                   const ConservativeState& conservative,
                                   const ConservativeState& physical_flux,
                                   Real wave_speed,
                                   Real un) {
            const Real star_denominator = wave_speed - speed_middle;
            const Real scale = std::max({1.0, std::abs(wave_speed), std::abs(speed_middle)});
            if (!std::isfinite(star_denominator)
                || std::abs(star_denominator) <= parameters_.denominator_tolerance * scale) {
                throw PhysicsError("HLLC star-state denominator is invalid");
            }
            const Real star_density = state[0] * (wave_speed - un) / star_denominator;
            const Real star_pressure
                = state[4] + state[0] * (wave_speed - un) * (speed_middle - un);
            const Normal3 star_velocity {
                state[1] + (speed_middle - un) * normal.x,
                state[2] + (speed_middle - un) * normal.y,
                state[3] + (speed_middle - un) * normal.z,
            };
            ConservativeState star {{
                star_density,
                star_density * star_velocity.x,
                star_density * star_velocity.y,
                star_density * star_velocity.z,
                ((wave_speed - un) * conservative[4] - state[4] * un + star_pressure * speed_middle)
                    / star_denominator,
            }};
            if (!finite_state(star)) {
                throw PhysicsError("HLLC star state is non-finite");
            }
            static_cast<void>(to_primitive(star, ideal));
            ConservativeState flux = physical_flux;
            for (int component = 0; component < euler_components; ++component) {
                const auto index = static_cast<std::size_t>(component);
                flux[index] += wave_speed * (star[index] - conservative[index]);
            }
            if (!finite_state(flux)) throw PhysicsError("HLLC star flux is non-finite");
            return flux;
        };
        try {
            const auto flux = speed_middle >= 0.0
                ? star_flux(left, left_conservative, left_flux, speed_left, un_left)
                : star_flux(right, right_conservative, right_flux, speed_right, un_right);
            return {flux, spectral_radius, "hllc", "hllc", RiemannFallbackReason::None, {}};
        } catch (const PhysicsError&) {
            return rusanov_result(left,
                                  right,
                                  normal,
                                  gas,
                                  floors,
                                  "hllc",
                                  RiemannFallbackReason::InvalidIntermediateState);
        }
    }

private:
    RiemannSolverParameters parameters_ {};
};

class RoeStrategy final : public IRiemannSolver {
public:
    explicit RoeStrategy(RiemannSolverParameters parameters)
        : parameters_(parameters)
    {
        parameters_.validate();
    }

    [[nodiscard]] std::string_view name() const noexcept override { return "roe"; }

    [[nodiscard]] RiemannResult solve(const PressurePrimitiveState& left,
                                      const PressurePrimitiveState& right,
                                      Normal3 unit_normal,
                                      const GasModel& gas,
                                      const NumericalFloors& floors) const override
    {
        const auto normal = checked_unit_normal(unit_normal);
        floors.validate();
        const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
        const auto fallback = [&](RiemannFallbackReason reason) {
            auto result = HllcStrategy(parameters_).solve(left, right, normal, gas, floors);
            result.fallback_path.insert(result.fallback_path.begin(), {"roe", "hllc", reason});
            result.requested_solver = "roe";
            result.fallback_reason = reason;
            return result;
        };
        try {
            const auto left_conservative = to_conservative(left, ideal);
            const auto right_conservative = to_conservative(right, ideal);
            const auto left_flux = euler_flux(left, normal, ideal);
            const auto right_flux = euler_flux(right, normal, ideal);
            const auto average = roe_average(left, right, normal, ideal);
            const auto basis = make_roe_characteristic_basis(left, right, normal, gas, floors, 3);
            ConservativeState jump {};
            for (int component = 0; component < euler_components; ++component) {
                const auto index = static_cast<std::size_t>(component);
                jump[index] = right_conservative[index] - left_conservative[index];
            }
            const auto strengths = project_characteristic(jump, basis);
            std::array<Real, euler_components> eigenvalues {{
                average.normal_velocity - average.sound_speed,
                average.normal_velocity,
                average.normal_velocity,
                average.normal_velocity,
                average.normal_velocity + average.sound_speed,
            }};
            const Real delta
                = parameters_.entropy_fix_coefficient
                * std::max(
                      {sound_speed(left, ideal), sound_speed(right, ideal), average.sound_speed});
            ConservativeState scaled_strengths {};
            for (int wave = 0; wave < euler_components; ++wave) {
                Real magnitude = std::abs(eigenvalues[static_cast<std::size_t>(wave)]);
                if (magnitude < delta) {
                    magnitude = 0.5 * (magnitude * magnitude / delta + delta);
                }
                scaled_strengths[static_cast<std::size_t>(wave)]
                    = magnitude * strengths[static_cast<std::size_t>(wave)];
            }
            const auto dissipation = restore_characteristic(scaled_strengths, basis);
            ConservativeState flux {};
            for (int component = 0; component < euler_components; ++component) {
                const auto index = static_cast<std::size_t>(component);
                flux[index]
                    = 0.5 * (left_flux[index] + right_flux[index]) - 0.5 * dissipation[index];
            }
            if (!finite_state(flux)) return fallback(RiemannFallbackReason::NonFiniteFlux);
            return {
                flux,
                std::abs(average.normal_velocity) + average.sound_speed,
                "roe",
                "roe",
                RiemannFallbackReason::None,
                {},
            };
        } catch (const PhysicsError&) {
            return fallback(RiemannFallbackReason::InvalidRoeAverage);
        }
    }

private:
    RiemannSolverParameters parameters_ {};
};

} // namespace

std::string_view riemann_fallback_reason_name(RiemannFallbackReason reason)
{
    switch (reason) {
    case RiemannFallbackReason::None: return "none";
    case RiemannFallbackReason::InvalidWaveSpeed: return "invalid_wave_speed";
    case RiemannFallbackReason::InvalidIntermediateState: return "invalid_intermediate_state";
    case RiemannFallbackReason::InvalidRoeAverage: return "invalid_roe_average";
    case RiemannFallbackReason::NonFiniteFlux: return "non_finite_flux";
    }
    throw std::invalid_argument("unknown Riemann fallback reason");
}

void RiemannDiagnostics::record(const RiemannResult& result, FaceDiagnosticLocation location)
{
    if (result.requested_solver.empty() || result.used_solver.empty()) {
        throw std::invalid_argument("Riemann diagnostic requires solver names");
    }
    ++total_faces;
    ++requested_faces[result.requested_solver];
    ++used_faces[result.used_solver];
    if (result.fallback_reason == RiemannFallbackReason::None) {
        if (result.requested_solver != result.used_solver) {
            throw std::logic_error("Riemann solver changed without a fallback reason");
        }
        if (!result.fallback_path.empty()) {
            throw std::logic_error("Riemann result has a path without a fallback");
        }
        return;
    }
    if (result.requested_solver == result.used_solver || result.fallback_path.empty()) {
        throw std::logic_error("Riemann fallback diagnostic is inconsistent");
    }
    std::string current = result.requested_solver;
    for (const auto& step : result.fallback_path) {
        const auto index = static_cast<std::size_t>(step.reason);
        if (step.from_solver != current || step.to_solver.empty()
            || step.from_solver == step.to_solver || step.reason == RiemannFallbackReason::None
            || index >= fallback_reasons.size()) {
            throw std::logic_error("Riemann fallback path is inconsistent");
        }
        ++fallback_reasons[index];
        fallback_events.push_back({location, step.from_solver, step.to_solver, step.reason});
        current = step.to_solver;
    }
    if (current != result.used_solver) {
        throw std::logic_error("Riemann fallback path does not reach the used solver");
    }
}

std::size_t RiemannDiagnostics::fallback_count() const noexcept
{
    return fallback_events.size();
}

std::size_t RiemannDiagnostics::fallback_count(RiemannFallbackReason reason) const
{
    const auto index = static_cast<std::size_t>(reason);
    if (index >= fallback_reasons.size()) {
        throw std::invalid_argument("unknown Riemann fallback reason");
    }
    return fallback_reasons[index];
}

void RiemannSolverRegistry::register_solver(std::string name, Factory factory)
{
    if (!valid_algorithm_name(name) || !factory) {
        throw std::invalid_argument("Riemann registration has an invalid name or factory");
    }
    auto probe = factory();
    if (!probe || probe->name() != name) {
        throw std::invalid_argument("Riemann factory name does not match its registry key");
    }
    if (!factories_.emplace(std::move(name), std::move(factory)).second) {
        throw std::invalid_argument("duplicate Riemann solver registration");
    }
}

bool RiemannSolverRegistry::contains(std::string_view name) const noexcept
{
    return factories_.find(std::string(name)) != factories_.end();
}

std::unique_ptr<IRiemannSolver> RiemannSolverRegistry::create(std::string_view name) const
{
    const auto iterator = factories_.find(std::string(name));
    if (iterator == factories_.end()) {
        throw std::invalid_argument("unknown Riemann solver: " + std::string(name));
    }
    auto result = iterator->second();
    if (!result || result->name() != iterator->first) {
        throw std::logic_error("registered Riemann factory returned an invalid strategy");
    }
    return result;
}

std::vector<std::string> RiemannSolverRegistry::names() const
{
    std::vector<std::string> result;
    result.reserve(factories_.size());
    for (const auto& [name, factory] : factories_) {
        static_cast<void>(factory);
        result.push_back(name);
    }
    std::sort(result.begin(), result.end());
    return result;
}

void RiemannSolverParameters::validate() const
{
    if (!std::isfinite(entropy_fix_coefficient) || entropy_fix_coefficient <= 0.0) {
        throw std::invalid_argument("Roe entropy-fix coefficient must be positive and finite");
    }
    if (!std::isfinite(denominator_tolerance) || denominator_tolerance <= 0.0
        || denominator_tolerance >= 1.0) {
        throw std::invalid_argument(
            "Riemann denominator tolerance must be finite and lie in (0,1)");
    }
    preconditioner.validate();
}

void RiemannConfig::validate() const
{
    if (!valid_algorithm_name(scheme)) {
        throw std::invalid_argument("Riemann scheme name is invalid");
    }
    parameters.validate();
    if (parameters.weiss_smith && scheme != "roe") {
        throw std::invalid_argument("Weiss-Smith preconditioning requires the Roe solver");
    }
}

void RiemannConfig::validate(const RiemannSolverRegistry& registry) const
{
    validate();
    if (!registry.contains(scheme)) {
        throw std::invalid_argument("unknown Riemann solver: " + scheme);
    }
}

std::string RiemannConfig::summary() const
{
    validate();
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<Real>::max_digits10)
           << "riemann_solver=" << scheme
           << ";roe_entropy_fix=" << parameters.entropy_fix_coefficient
           << ";riemann_denominator_tolerance=" << parameters.denominator_tolerance;
    if (parameters.weiss_smith) {
        stream << ";preconditioner=weiss_smith"
               << ";preconditioner_mach_cutoff=" << parameters.preconditioner.mach_cutoff
               << ";preconditioner_viscous_cutoff="
               << parameters.preconditioner.viscous_cutoff;
    }
    return stream.str();
}

std::string RiemannConfig::restart_signature() const
{
    return "riemann_config_v1;" + summary();
}

RiemannSolverRegistry
RiemannSolverRegistry::with_builtins(const RiemannSolverParameters& parameters)
{
    parameters.validate();
    RiemannSolverRegistry result;
    result.register_solver("rusanov", [] { return std::make_unique<RusanovStrategy>(); });
    result.register_solver("hllc",
                           [parameters] { return std::make_unique<HllcStrategy>(parameters); });
    result.register_solver("roe",
                           [parameters] { return std::make_unique<RoeStrategy>(parameters); });
    return result;
}

std::string_view riemann_solver_name(RiemannSolverKind kind)
{
    switch (kind) {
    case RiemannSolverKind::Rusanov: return "rusanov";
    case RiemannSolverKind::Hllc: return "hllc";
    case RiemannSolverKind::Roe: return "roe";
    }
    throw std::invalid_argument("unknown Riemann solver kind");
}

RiemannSolver::RiemannSolver(RiemannSolverKind kind, RiemannSolverParameters parameters)
    : parameters_(parameters)
{
    parameters_.validate();
    const auto registry = RiemannSolverRegistry::with_builtins(parameters_);
    implementation_ = registry.create(riemann_solver_name(kind));
}

RiemannSolver::RiemannSolver(std::string_view name,
                             const RiemannSolverRegistry& registry,
                             RiemannSolverParameters parameters)
    : implementation_(registry.create(name))
    , parameters_(parameters)
{
    parameters_.validate();
}

std::string_view RiemannSolver::name() const noexcept
{
    return implementation_->name();
}

std::string RiemannSolver::summary() const
{
    RiemannConfig config;
    config.scheme = std::string(name());
    config.parameters = parameters_;
    return config.summary();
}

std::string RiemannSolver::restart_signature() const
{
    return "riemann_v3;" + summary();
}

RiemannResult RiemannSolver::solve(const PressurePrimitiveState& left,
                                   const PressurePrimitiveState& right,
                                   Normal3 unit_normal,
                                   const GasModel& gas,
                                   const NumericalFloors& floors,
                                   RiemannFaceContext context) const
{
    const auto normal = checked_unit_normal(unit_normal);
    if (!std::isfinite(context.viscous_speed) || context.viscous_speed < 0.0) {
        throw std::invalid_argument("Riemann viscous reference speed is invalid");
    }
    auto result = parameters_.weiss_smith
        ? preconditioned_roe_result(left,
                                    right,
                                    normal,
                                    gas,
                                    floors,
                                    context.viscous_speed,
                                    parameters_)
        : implementation_->solve(left, right, normal, gas, floors);
    if (result.requested_solver.empty()) result.requested_solver = std::string(name());
    if (result.used_solver.empty()) result.used_solver = result.requested_solver;
    if (result.requested_solver != name()) {
        throw std::logic_error("Riemann result requested-solver diagnostic is inconsistent");
    }
    if (result.fallback_reason == RiemannFallbackReason::None) {
        if (result.used_solver != result.requested_solver || !result.fallback_path.empty()) {
            throw std::logic_error("Riemann result changed solver without a fallback path");
        }
    } else {
        if (result.fallback_path.empty()) {
            throw std::logic_error("Riemann fallback result has no path");
        }
        std::string current = result.requested_solver;
        for (const auto& step : result.fallback_path) {
            if (step.from_solver != current || step.to_solver.empty()
                || step.from_solver == step.to_solver
                || step.reason == RiemannFallbackReason::None) {
                throw std::logic_error("Riemann fallback path is inconsistent");
            }
            current = step.to_solver;
        }
        if (current != result.used_solver) {
            throw std::logic_error("Riemann fallback path does not reach the used solver");
        }
    }
    if (!std::isfinite(result.spectral_radius) || result.spectral_radius < 0.0) {
        throw PhysicsError("Riemann solver returned an invalid spectral radius");
    }
    for (const auto value : result.flux_per_unit_area) {
        if (!std::isfinite(value)) {
            throw PhysicsError("Riemann solver returned a non-finite flux");
        }
    }
    return result;
}

ConservativeState RiemannSolver::flux(const PressurePrimitiveState& left,
                                      const PressurePrimitiveState& right,
                                      Normal3 unit_normal,
                                      const GasModel& gas,
                                      const NumericalFloors& floors,
                                      RiemannFaceContext context) const
{
    return solve(left, right, unit_normal, gas, floors, context).flux_per_unit_area;
}

} // namespace wcns
