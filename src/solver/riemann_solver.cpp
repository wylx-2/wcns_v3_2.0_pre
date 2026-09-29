#include <wcns/solver/riemann_solver.hpp>
#include <wcns/solver/wcns_reconstruction.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iterator>
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

struct EulerFaceData {
    ConservativeState conservative {};
    ConservativeState flux {};
    Real normal_velocity = 0.0;
    Real sound_speed = 0.0;
    Real enthalpy = 0.0;
};

EulerFaceData make_face_data(const PressurePrimitiveState& state,
                             Normal3 normal,
                             const IdealGas& gas)
{
    // to_conservative() performs the public-state and gas validation once.
    // The old Riemann paths called it again through euler_flux(), and then
    // repeated the same validation in sound_speed().
    EulerFaceData result;
    result.conservative = to_conservative(state, gas);
    const Real rho = state[0];
    const Real pressure = state[4];
    result.normal_velocity = normal_velocity(state, normal);
    result.sound_speed = std::sqrt(gas.gamma * pressure / rho);
    result.enthalpy = (result.conservative[4] + pressure) / rho;
    result.flux = {{
        rho * result.normal_velocity,
        result.conservative[1] * result.normal_velocity + pressure * normal.x,
        result.conservative[2] * result.normal_velocity + pressure * normal.y,
        result.conservative[3] * result.normal_velocity + pressure * normal.z,
        (result.conservative[4] + pressure) * result.normal_velocity,
    }};
    if (!std::isfinite(result.sound_speed) || !std::isfinite(result.enthalpy)
        || !std::all_of(result.flux.begin(), result.flux.end(), [](Real value) {
               return std::isfinite(value);
           })) {
        throw PhysicsError("Euler face data is non-finite");
    }
    return result;
}

bool finite_state(const ConservativeState& state)
{
    return std::all_of(state.begin(), state.end(), [](Real value) { return std::isfinite(value); });
}

Real dot(Normal3 lhs, Normal3 rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Normal3 cross(Normal3 lhs, Normal3 rhs)
{
    return {
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x,
    };
}

Normal3 normalized(Normal3 vector)
{
    const Real magnitude = std::sqrt(dot(vector, vector));
    if (!std::isfinite(magnitude) || magnitude <= std::numeric_limits<Real>::epsilon()) {
        throw PhysicsError("cannot normalize a degenerate direction");
    }
    return {vector.x / magnitude, vector.y / magnitude, vector.z / magnitude};
}

ConservativeState conservative_jump(const ConservativeState& left,
                                     const ConservativeState& right)
{
    ConservativeState result {};
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        result[index] = right[index] - left[index];
    }
    return result;
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
    Real density = 0.0;
    Normal3 velocity {};
    Real normal_velocity = 0.0;
    Real enthalpy = 0.0;
    Real sound_speed = 0.0;
};

RoeAverage roe_average(const PressurePrimitiveState& left,
                       const PressurePrimitiveState& right,
                       const EulerFaceData& left_data,
                       const EulerFaceData& right_data,
                       Normal3 normal,
                       const IdealGas& gas)
{
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
    result.density = root_left * root_right;
    result.velocity = {
        average(left[1], right[1]),
        average(left[2], right[2]),
        average(left[3], right[3]),
    };
    result.enthalpy = average(left_data.enthalpy, right_data.enthalpy);
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

struct RoePikeStrengths {
    Real acoustic_minus = 0.0;
    Real entropy = 0.0;
    Normal3 shear_momentum {};
    Real acoustic_plus = 0.0;
};

RoePikeStrengths roe_pike_strengths(const ConservativeState& jump,
                                    const RoeAverage& average,
                                    Normal3 normal,
                                    Real gamma)
{
    const Real density_jump = jump[0];
    const Normal3 momentum_jump {jump[1], jump[2], jump[3]};
    const Real normal_momentum_jump = dot(momentum_jump, normal);
    const Real velocity_dot_momentum = dot(average.velocity, momentum_jump);
    const Real kinetic = 0.5 * dot(average.velocity, average.velocity);
    const Real inverse_sound_squared
        = 1.0 / (average.sound_speed * average.sound_speed);
    const Real pressure_over_sound_squared = (gamma - 1.0) * inverse_sound_squared
        * (jump[4] - velocity_dot_momentum + kinetic * density_jump);
    const Real normal_wave = (normal_momentum_jump
                              - average.normal_velocity * density_jump)
        / average.sound_speed;
    const Normal3 tangential_velocity {
        average.velocity.x - average.normal_velocity * normal.x,
        average.velocity.y - average.normal_velocity * normal.y,
        average.velocity.z - average.normal_velocity * normal.z,
    };
    const Normal3 tangential_momentum_jump {
        momentum_jump.x - normal_momentum_jump * normal.x,
        momentum_jump.y - normal_momentum_jump * normal.y,
        momentum_jump.z - normal_momentum_jump * normal.z,
    };
    return {
        0.5 * (pressure_over_sound_squared - normal_wave),
        density_jump - pressure_over_sound_squared,
        {
            tangential_momentum_jump.x - tangential_velocity.x * density_jump,
            tangential_momentum_jump.y - tangential_velocity.y * density_jump,
            tangential_momentum_jump.z - tangential_velocity.z * density_jump,
        },
        0.5 * (pressure_over_sound_squared + normal_wave),
    };
}

ConservativeState roe_pike_dissipation(const RoePikeStrengths& strengths,
                                       const RoeAverage& average,
                                       Normal3 normal,
                                       Real acoustic_minus_magnitude,
                                       Real contact_magnitude,
                                       Real acoustic_plus_magnitude)
{
    const Real minus = acoustic_minus_magnitude * strengths.acoustic_minus;
    const Real entropy = contact_magnitude * strengths.entropy;
    const Real plus = acoustic_plus_magnitude * strengths.acoustic_plus;
    const Real mass = minus + entropy + plus;
    const Real acoustic_momentum = average.sound_speed * (plus - minus);
    const Normal3 shear {
        contact_magnitude * strengths.shear_momentum.x,
        contact_magnitude * strengths.shear_momentum.y,
        contact_magnitude * strengths.shear_momentum.z,
    };
    const Normal3 tangential_velocity {
        average.velocity.x - average.normal_velocity * normal.x,
        average.velocity.y - average.normal_velocity * normal.y,
        average.velocity.z - average.normal_velocity * normal.z,
    };
    const Real kinetic = 0.5 * dot(average.velocity, average.velocity);
    return {{
        mass,
        average.velocity.x * mass + normal.x * acoustic_momentum + shear.x,
        average.velocity.y * mass + normal.y * acoustic_momentum + shear.y,
        average.velocity.z * mass + normal.z * acoustic_momentum + shear.z,
        (average.enthalpy - average.normal_velocity * average.sound_speed) * minus
            + kinetic * entropy
            + (average.enthalpy + average.normal_velocity * average.sound_speed) * plus
            + dot(tangential_velocity, shear),
    }};
}

RiemannResult rusanov_result(const EulerFaceData& left,
                             const EulerFaceData& right,
                             std::string requested,
                             RiemannFallbackReason reason)
{
    const Real spectral_radius
        = std::max(std::abs(left.normal_velocity) + left.sound_speed,
                   std::abs(right.normal_velocity) + right.sound_speed);
    ConservativeState flux {};
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        flux[index] = 0.5 * (left.flux[index] + right.flux[index])
            - 0.5 * spectral_radius * (right.conservative[index] - left.conservative[index]);
    }
    RiemannResult result {
        flux,
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

RiemannResult rusanov_result(const PressurePrimitiveState& left,
                             const PressurePrimitiveState& right,
                             Normal3 normal,
                             const GasModel& gas,
                             const NumericalFloors& floors,
                             std::string requested,
                             RiemannFallbackReason reason)
{
    const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
    return rusanov_result(make_face_data(left, normal, ideal),
                          make_face_data(right, normal, ideal),
                          std::move(requested),
                          reason);
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
    const auto left_data = make_face_data(left, normal, ideal);
    const auto right_data = make_face_data(right, normal, ideal);
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
        jump[index] = right_data.conservative[index] - left_data.conservative[index];
    }
    auto primitive_jump = matrix_vector_product(physical.inverse, jump);
    for (Real& value : primitive_jump) value *= spectral_radius;
    const auto dissipation
        = matrix_vector_product(preconditioned.gamma, primitive_jump);
    ConservativeState flux {};
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        flux[index] = 0.5 * (left_data.flux[index] + right_data.flux[index])
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
        const auto left_data = make_face_data(left, normal, ideal);
        const auto right_data = make_face_data(right, normal, ideal);
        const auto average = roe_average(left, right, left_data, right_data, normal, ideal);
        const Real left_sound = left_data.sound_speed;
        const Real right_sound = right_data.sound_speed;
        const Real left_reference = weiss_smith_reference_speed(
            state_speed(left), left_sound, viscous_speed, parameters.preconditioner);
        const Real right_reference = weiss_smith_reference_speed(
            state_speed(right), right_sound, viscous_speed, parameters.preconditioner);
        const auto basis = make_roe_characteristic_basis(left, right, normal, gas, floors, 3);
        ConservativeState jump {};
        for (int component = 0; component < euler_components; ++component) {
            const auto index = static_cast<std::size_t>(component);
            jump[index] = right_data.conservative[index] - left_data.conservative[index];
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
                = 0.5 * (left_data.flux[index] + right_data.flux[index])
                - 0.5 * dissipation[index];
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

class HllStrategy final : public IRiemannSolver {
public:
    explicit HllStrategy(RiemannSolverParameters parameters)
        : parameters_(parameters)
    {
        parameters_.validate();
    }

    [[nodiscard]] std::string_view name() const noexcept override { return "hll"; }

    [[nodiscard]] RiemannResult solve(const PressurePrimitiveState& left,
                                      const PressurePrimitiveState& right,
                                      Normal3 unit_normal,
                                      const GasModel& gas,
                                      const NumericalFloors& floors) const override
    {
        const auto normal = checked_unit_normal(unit_normal);
        floors.validate();
        const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
        const auto left_data = make_face_data(left, normal, ideal);
        const auto right_data = make_face_data(right, normal, ideal);
        const Real un_left = left_data.normal_velocity;
        const Real un_right = right_data.normal_velocity;

        try {
            const auto average = roe_average(left, right, left_data, right_data, normal, ideal);
            const Real speed_left = std::min(
                un_left - left_data.sound_speed,
                average.normal_velocity - average.sound_speed);
            const Real speed_right = std::max(
                un_right + right_data.sound_speed,
                average.normal_velocity + average.sound_speed);
            const Real spectral_radius = std::max(std::abs(speed_left), std::abs(speed_right));
            const Real denominator = speed_right - speed_left;
            const Real scale = std::max({1.0, std::abs(speed_left), std::abs(speed_right)});
            if (!std::isfinite(speed_left) || !std::isfinite(speed_right)
                || !std::isfinite(denominator) || speed_left >= speed_right
                || denominator <= parameters_.denominator_tolerance * scale) {
                return rusanov_result(left,
                                      right,
                                      normal,
                                      gas,
                                      floors,
                                      "hll",
                                      RiemannFallbackReason::InvalidWaveSpeed);
            }
            if (speed_left >= 0.0) {
                return {left_data.flux,
                        spectral_radius,
                        "hll",
                        "hll",
                        RiemannFallbackReason::None,
                        {}};
            }
            if (speed_right <= 0.0) {
                return {
                    right_data.flux,
                    spectral_radius,
                    "hll",
                    "hll",
                    RiemannFallbackReason::None,
                    {}};
            }

            ConservativeState flux {};
            for (int component = 0; component < euler_components; ++component) {
                const auto index = static_cast<std::size_t>(component);
                flux[index] = (speed_right * left_data.flux[index]
                               - speed_left * right_data.flux[index]
                               + speed_left * speed_right
                                    * (right_data.conservative[index]
                                       - left_data.conservative[index]))
                    / denominator;
            }
            if (!finite_state(flux)) {
                return rusanov_result(left,
                                      right,
                                      normal,
                                      gas,
                                      floors,
                                      "hll",
                                      RiemannFallbackReason::NonFiniteFlux);
            }
            return {flux, spectral_radius, "hll", "hll", RiemannFallbackReason::None, {}};
        } catch (const PhysicsError&) {
            return rusanov_result(left,
                                  right,
                                  normal,
                                  gas,
                                  floors,
                                  "hll",
                                  RiemannFallbackReason::InvalidWaveSpeed);
        }
    }

private:
    RiemannSolverParameters parameters_ {};
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
        const auto left_data = make_face_data(left, normal, ideal);
        const auto right_data = make_face_data(right, normal, ideal);
        const Real un_left = left_data.normal_velocity;
        const Real un_right = right_data.normal_velocity;
        const Real sound_left = left_data.sound_speed;
        const Real sound_right = right_data.sound_speed;
        RoeAverage average;
        try {
            average = roe_average(left, right, left_data, right_data, normal, ideal);
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
            return {left_data.flux,
                    spectral_radius,
                    "hllc",
                    "hllc",
                    RiemannFallbackReason::None,
                    {}};
        }
        if (speed_right <= 0.0) {
            return {right_data.flux,
                    spectral_radius,
                    "hllc",
                    "hllc",
                    RiemannFallbackReason::None,
                    {}};
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
                ? star_flux(left,
                            left_data.conservative,
                            left_data.flux,
                            speed_left,
                            un_left)
                : star_flux(right,
                            right_data.conservative,
                            right_data.flux,
                            speed_right,
                            un_right);
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

enum class RoeVariant {
    Standard,
    LiGuAllSpeed,
    RieperAllSpeed,
};

class RoeStrategy final : public IRiemannSolver {
public:
    RoeStrategy(RiemannSolverParameters parameters, RoeVariant variant)
        : parameters_(parameters)
        , variant_(variant)
    {
        parameters_.validate();
    }

    [[nodiscard]] std::string_view name() const noexcept override
    {
        switch (variant_) {
        case RoeVariant::Standard: return "roe";
        case RoeVariant::LiGuAllSpeed: return "roe_all_speed";
        case RoeVariant::RieperAllSpeed: return "roe_all_speed_rieper";
        }
        return "roe";
    }

    [[nodiscard]] RiemannResult solve(const PressurePrimitiveState& left,
                                      const PressurePrimitiveState& right,
                                      Normal3 unit_normal,
                                      const GasModel& gas,
                                      const NumericalFloors& floors) const override
    {
        const auto normal = checked_unit_normal(unit_normal);
        floors.validate();
        const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
        const std::string requested(name());
        const auto fallback = [&](RiemannFallbackReason reason) {
            auto result = HllcStrategy(parameters_).solve(left, right, normal, gas, floors);
            result.fallback_path.insert(
                result.fallback_path.begin(), {requested, "hllc", reason});
            result.requested_solver = requested;
            result.fallback_reason = reason;
            return result;
        };
        try {
            const auto left_data = make_face_data(left, normal, ideal);
            const auto right_data = make_face_data(right, normal, ideal);
            const auto average
                = roe_average(left, right, left_data, right_data, normal, ideal);
            const auto jump
                = conservative_jump(left_data.conservative, right_data.conservative);
            auto strengths = roe_pike_strengths(jump, average, normal, ideal.gamma);
            if (variant_ == RoeVariant::RieperAllSpeed) {
                // Rieper-type low-Mach Roe fix: scale only the velocity-jump
                // contribution to the two acoustic wave strengths.  The
                // physical Euler flux and physical wave speeds are unchanged,
                // so this remains an explicit conservative spatial flux.
                const Real speed = std::sqrt(average.velocity.x * average.velocity.x
                                             + average.velocity.y * average.velocity.y
                                             + average.velocity.z * average.velocity.z);
                const Real velocity_factor = std::min(1.0, speed / average.sound_speed);
                const Real delta_normal_velocity
                    = normal_velocity(right, normal) - normal_velocity(left, normal);
                const Real acoustic_correction = (1.0 - velocity_factor) * average.density
                    * delta_normal_velocity / (2.0 * average.sound_speed);
                strengths.acoustic_minus += acoustic_correction;
                strengths.acoustic_plus -= acoustic_correction;
            }
            Real acoustic_speed = average.sound_speed;
            Real flux_dissipation_scale = 0.5;
            Real all_speed_factor = 1.0;
            if (variant_ == RoeVariant::LiGuAllSpeed) {
                // Li--Gu all-speed Roe, test_pdf_3.pdf equations (5.83) and
                // (5.92).  The eigenvectors remain the Roe--Pike vectors;
                // only the two acoustic eigenvalues are changed.
                const Real local_normal_mach
                    = std::abs(average.normal_velocity) / average.sound_speed;
                const Real mach
                    = std::max(local_normal_mach, parameters_.all_speed.reference_mach);
                if (mach < 1.0) {
                    const Real mach_squared = mach * mach;
                    all_speed_factor
                        = mach
                        * std::sqrt((4.0 + (1.0 - mach_squared) * (1.0 - mach_squared))
                                    / (1.0 + mach_squared));
                }
                all_speed_factor = std::min(1.0, all_speed_factor);
                acoustic_speed = all_speed_factor * average.sound_speed;
                flux_dissipation_scale = parameters_.all_speed.dissipation_scale;
            }
            const Real delta
                = parameters_.entropy_fix_coefficient
                * std::max(
                      {left_data.sound_speed, right_data.sound_speed, average.sound_speed});
            const auto wave_magnitude = [&](Real eigenvalue) {
                Real magnitude = std::abs(eigenvalue);
                if (variant_ != RoeVariant::LiGuAllSpeed && magnitude < delta) {
                    magnitude = 0.5 * (magnitude * magnitude / delta + delta);
                }
                return magnitude;
            };
            const auto dissipation = roe_pike_dissipation(
                strengths,
                average,
                normal,
                wave_magnitude(average.normal_velocity - acoustic_speed),
                wave_magnitude(average.normal_velocity),
                wave_magnitude(average.normal_velocity + acoustic_speed));
            ConservativeState flux {};
            for (int component = 0; component < euler_components; ++component) {
                const auto index = static_cast<std::size_t>(component);
                flux[index] = 0.5 * (left_data.flux[index] + right_data.flux[index])
                    - flux_dissipation_scale * dissipation[index];
            }
            if (variant_ == RoeVariant::LiGuAllSpeed) {
                ConservativeState averaged_conservative {};
                for (int component = 0; component < euler_components; ++component) {
                    const auto index = static_cast<std::size_t>(component);
                    averaged_conservative[index]
                        = 0.5
                        * (left_data.conservative[index] + right_data.conservative[index]);
                }
                const Real averaged_pressure = to_primitive(averaged_conservative, ideal)[4];
                ConservativeState pressure_vector = averaged_conservative;
                pressure_vector[4] += averaged_pressure;
                const Real pressure_correction
                    = (1.0 - all_speed_factor) * parameters_.all_speed.pressure_coefficient
                    * (right[4] - left[4])
                    / (parameters_.all_speed.reference_mach * average.density
                       * average.sound_speed);
                for (int component = 0; component < euler_components; ++component) {
                    const auto index = static_cast<std::size_t>(component);
                    flux[index] -= pressure_correction * pressure_vector[index];
                }
            }
            if (!finite_state(flux)) return fallback(RiemannFallbackReason::NonFiniteFlux);
            return {
                flux,
                std::abs(average.normal_velocity) + acoustic_speed,
                requested,
                requested,
                RiemannFallbackReason::None,
                {},
            };
        } catch (const PhysicsError&) {
            return fallback(RiemannFallbackReason::InvalidRoeAverage);
        }
    }

private:
    RiemannSolverParameters parameters_ {};
    RoeVariant variant_ = RoeVariant::Standard;
};

class RotatedRoeStrategy final : public IRiemannSolver {
public:
    explicit RotatedRoeStrategy(RiemannSolverParameters parameters)
        : parameters_(parameters)
    {
        parameters_.validate();
    }

    [[nodiscard]] std::string_view name() const noexcept override { return "roe_rotated"; }

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
            result.fallback_path.insert(
                result.fallback_path.begin(), {"roe_rotated", "hllc", reason});
            result.requested_solver = "roe_rotated";
            result.fallback_reason = reason;
            return result;
        };

        try {
            const auto left_data = make_face_data(left, normal, ideal);
            const auto right_data = make_face_data(right, normal, ideal);
            const auto jump
                = conservative_jump(left_data.conservative, right_data.conservative);
            const auto average
                = roe_average(left, right, left_data, right_data, normal, ideal);
            const Normal3 velocity_jump {
                right[1] - left[1],
                right[2] - left[2],
                right[3] - left[3],
            };
            const Real jump_norm = std::sqrt(dot(velocity_jump, velocity_jump));
            const Real velocity_scale = std::max(
                {1.0,
                 std::abs(left[1]),
                 std::abs(left[2]),
                 std::abs(left[3]),
                 std::abs(right[1]),
                 std::abs(right[2]),
                 std::abs(right[3])});
            const Normal3 first = jump_norm > parameters_.denominator_tolerance * velocity_scale
                ? normalized(velocity_jump)
                : normal;

            const std::array<Real, 3> absolute_components {
                std::abs(first.x), std::abs(first.y), std::abs(first.z)};
            const auto least_aligned = static_cast<std::size_t>(std::distance(
                absolute_components.begin(),
                std::min_element(absolute_components.begin(), absolute_components.end())));
            const std::array<Normal3, 3> coordinate_axes {{
                {1.0, 0.0, 0.0},
                {0.0, 1.0, 0.0},
                {0.0, 0.0, 1.0},
            }};
            const Normal3 second = normalized(cross(coordinate_axes[least_aligned], first));
            const Normal3 third = normalized(cross(first, second));
            const std::array<Normal3, 3> directions {{first, second, third}};

            const Real pressure_jump = right[4] - left[4];
            const ConservativeState roe_enthalpy_state {{
                average.density,
                average.density * average.velocity.x,
                average.density * average.velocity.y,
                average.density * average.velocity.z,
                average.density * average.enthalpy,
            }};
            ConservativeState dissipation {};
            Real spectral_radius = 0.0;
            for (const auto direction : directions) {
                const Real weight = std::abs(dot(normal, direction));
                if (weight <= std::numeric_limits<Real>::epsilon()) continue;
                const Real directional_velocity = dot(average.velocity, direction);
                const Real directional_jump = dot(velocity_jump, direction);
                const Real velocity_magnitude = std::abs(directional_velocity);
                const Real sign = directional_velocity > 0.0
                    ? 1.0
                    : (directional_velocity < 0.0 ? -1.0 : 0.0);
                const Real t = 2.0 * sign
                    * std::min(velocity_magnitude, average.sound_speed);
                const Real s
                    = 2.0 * std::max(0.0, average.sound_speed - velocity_magnitude);
                const Real delta_velocity = t * directional_jump / (2.0 * average.sound_speed)
                    + s * pressure_jump
                        / (2.0 * average.density * average.sound_speed * average.sound_speed);
                const Real delta_pressure = 0.5 * s * average.density * directional_jump
                    + t * pressure_jump / (2.0 * average.sound_speed);
                const ConservativeState direction_vector {{
                    0.0,
                    direction.x,
                    direction.y,
                    direction.z,
                    directional_velocity,
                }};
                for (int component = 0; component < euler_components; ++component) {
                    const auto index = static_cast<std::size_t>(component);
                    dissipation[index] += weight
                        * (velocity_magnitude * jump[index]
                           + delta_velocity * roe_enthalpy_state[index]
                           + delta_pressure * direction_vector[index]);
                }
                spectral_radius
                    += weight * (velocity_magnitude + average.sound_speed);
            }

            ConservativeState flux {};
            for (int component = 0; component < euler_components; ++component) {
                const auto index = static_cast<std::size_t>(component);
                flux[index]
                    = 0.5 * (left_data.flux[index] + right_data.flux[index])
                    - 0.5 * dissipation[index];
            }
            if (!finite_state(flux) || !std::isfinite(spectral_radius)) {
                return fallback(RiemannFallbackReason::NonFiniteFlux);
            }
            return {flux,
                    spectral_radius,
                    "roe_rotated",
                    "roe_rotated",
                    RiemannFallbackReason::None,
                    {}};
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

void LiGuAllSpeedRoeParameters::validate() const
{
    if (!std::isfinite(reference_mach) || reference_mach <= 0.0 || reference_mach > 1.0) {
        throw std::invalid_argument(
            "Li-Gu all-speed Roe reference Mach number must lie in (0,1]");
    }
    if (!std::isfinite(dissipation_scale) || dissipation_scale <= 0.0
        || dissipation_scale > 1.0) {
        throw std::invalid_argument(
            "Li-Gu all-speed Roe dissipation scale must lie in (0,1]");
    }
    if (!std::isfinite(pressure_coefficient) || pressure_coefficient < 0.0
        || pressure_coefficient > 1.0) {
        throw std::invalid_argument(
            "Li-Gu all-speed Roe pressure coefficient must lie in [0,1]");
    }
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
    all_speed.validate();
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
    if (scheme == "roe_all_speed") {
        stream << ";all_speed_reference_mach=" << parameters.all_speed.reference_mach
               << ";all_speed_dissipation_scale="
               << parameters.all_speed.dissipation_scale
               << ";all_speed_pressure_coefficient="
               << parameters.all_speed.pressure_coefficient;
    }
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
    return "riemann_config_v2;" + summary();
}

RiemannSolverRegistry
RiemannSolverRegistry::with_builtins(const RiemannSolverParameters& parameters)
{
    parameters.validate();
    RiemannSolverRegistry result;
    result.register_solver("rusanov", [] { return std::make_unique<RusanovStrategy>(); });
    result.register_solver("hll",
                           [parameters] { return std::make_unique<HllStrategy>(parameters); });
    result.register_solver("hllc",
                           [parameters] { return std::make_unique<HllcStrategy>(parameters); });
    result.register_solver("roe",
                           [parameters] {
                               return std::make_unique<RoeStrategy>(parameters,
                                                                    RoeVariant::Standard);
                           });
    result.register_solver("roe_all_speed",
                           [parameters] {
                               return std::make_unique<RoeStrategy>(parameters,
                                                                    RoeVariant::LiGuAllSpeed);
                           });
    result.register_solver("roe_all_speed_rieper",
                           [parameters] {
                               return std::make_unique<RoeStrategy>(parameters,
                                                                    RoeVariant::RieperAllSpeed);
                           });
    result.register_solver("roe_rotated",
                           [parameters] {
                               return std::make_unique<RotatedRoeStrategy>(parameters);
                           });
    return result;
}

std::string_view riemann_solver_name(RiemannSolverKind kind)
{
    switch (kind) {
    case RiemannSolverKind::Rusanov: return "rusanov";
    case RiemannSolverKind::Hll: return "hll";
    case RiemannSolverKind::Hllc: return "hllc";
    case RiemannSolverKind::Roe: return "roe";
    case RiemannSolverKind::AllSpeedRoe: return "roe_all_speed";
    case RiemannSolverKind::RieperAllSpeedRoe: return "roe_all_speed_rieper";
    case RiemannSolverKind::RotatedRoe: return "roe_rotated";
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
    return "riemann_v4;" + summary();
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
