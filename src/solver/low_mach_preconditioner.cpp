#include <wcns/solver/low_mach_preconditioner.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

Matrix5 invert(Matrix5 matrix)
{
    Matrix5 inverse {};
    for (int row = 0; row < 5; ++row)
        inverse[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)] = 1.0;
    for (int pivot = 0; pivot < 5; ++pivot) {
        int selected = pivot;
        for (int row = pivot + 1; row < 5; ++row) {
            if (std::abs(matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(pivot)])
                > std::abs(matrix[static_cast<std::size_t>(selected)]
                                  [static_cast<std::size_t>(pivot)])) {
                selected = row;
            }
        }
        const Real value
            = matrix[static_cast<std::size_t>(selected)][static_cast<std::size_t>(pivot)];
        if (!std::isfinite(value) || std::abs(value) <= 1.0e-14) {
            throw PhysicsError("Weiss-Smith matrix is singular");
        }
        if (selected != pivot) {
            std::swap(matrix[static_cast<std::size_t>(selected)],
                      matrix[static_cast<std::size_t>(pivot)]);
            std::swap(inverse[static_cast<std::size_t>(selected)],
                      inverse[static_cast<std::size_t>(pivot)]);
        }
        const Real scale
            = 1.0 / matrix[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(pivot)];
        for (int column = 0; column < 5; ++column) {
            matrix[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(column)] *= scale;
            inverse[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(column)] *= scale;
        }
        for (int row = 0; row < 5; ++row) {
            if (row == pivot) continue;
            const Real factor
                = matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(pivot)];
            for (int column = 0; column < 5; ++column) {
                matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] -= factor
                    * matrix[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(column)];
                inverse[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)] -= factor
                    * inverse[static_cast<std::size_t>(pivot)][static_cast<std::size_t>(column)];
            }
        }
    }
    return inverse;
}

} // namespace

void WeissSmithParameters::validate() const
{
    if (!std::isfinite(mach_cutoff) || mach_cutoff <= 0.0 || mach_cutoff > 1.0) {
        throw std::invalid_argument("Weiss-Smith Mach cutoff must lie in (0,1]");
    }
    if (!std::isfinite(viscous_cutoff) || viscous_cutoff < 0.0
        || viscous_cutoff > 10.0) {
        throw std::invalid_argument("Weiss-Smith viscous cutoff must lie in [0,10]");
    }
}

Matrix5 multiply(const Matrix5& left, const Matrix5& right)
{
    Matrix5 result {};
    for (int row = 0; row < 5; ++row) {
        for (int column = 0; column < 5; ++column) {
            for (int inner = 0; inner < 5; ++inner) {
                result[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
                    += left[static_cast<std::size_t>(row)][static_cast<std::size_t>(inner)]
                    * right[static_cast<std::size_t>(inner)][static_cast<std::size_t>(column)];
            }
        }
    }
    return result;
}

Real weiss_smith_reference_speed(Real velocity_magnitude,
                                 Real sound_speed,
                                 Real viscous_speed,
                                 const WeissSmithParameters& parameters)
{
    parameters.validate();
    if (!std::isfinite(velocity_magnitude) || velocity_magnitude < 0.0
        || !std::isfinite(sound_speed) || sound_speed <= 0.0
        || !std::isfinite(viscous_speed) || viscous_speed < 0.0) {
        throw std::invalid_argument("invalid Weiss-Smith reference-speed input");
    }
    return std::min(sound_speed,
                    std::max({velocity_magnitude,
                              parameters.mach_cutoff * sound_speed,
                              parameters.viscous_cutoff * viscous_speed}));
}

std::array<Real, 5>
weiss_smith_eigenvalues(Real normal_velocity, Real sound_speed, Real beta)
{
    if (!std::isfinite(normal_velocity) || !std::isfinite(sound_speed)
        || sound_speed <= 0.0 || !std::isfinite(beta) || beta <= 0.0 || beta > 1.0) {
        throw std::invalid_argument("invalid Weiss-Smith eigenvalue input");
    }
    const Real discriminant = std::sqrt(
        (1.0 - beta) * (1.0 - beta) * normal_velocity * normal_velocity
        + 4.0 * beta * sound_speed * sound_speed);
    return {{normal_velocity,
             normal_velocity,
             normal_velocity,
             0.5 * ((1.0 + beta) * normal_velocity - discriminant),
             0.5 * ((1.0 + beta) * normal_velocity + discriminant)}};
}

WeissSmithState weiss_smith_state(const TemperaturePrimitiveState& state,
                                  Normal3 unit_normal,
                                  const GasModel& gas,
                                  const ReferenceScales& reference,
                                  const NumericalFloors& floors,
                                  const WeissSmithParameters& parameters,
                                  Real viscous_speed,
                                  int dimension)
{
    parameters.validate();
    if (!std::isfinite(viscous_speed) || viscous_speed < 0.0) {
        throw std::invalid_argument("Weiss-Smith viscous reference speed is invalid");
    }
    const Real normal_norm = std::sqrt(unit_normal.x * unit_normal.x
                                       + unit_normal.y * unit_normal.y
                                       + unit_normal.z * unit_normal.z);
    if (!std::isfinite(normal_norm) || std::abs(normal_norm - 1.0) > 1.0e-12) {
        throw std::invalid_argument("Weiss-Smith preconditioner requires a unit normal");
    }
    const auto pressure = pressure_primitive(state, gas, reference, floors, dimension);
    static_cast<void>(pressure);
    const Real rho = state[temperature_density];
    const Real u = state[temperature_velocity_x];
    const Real v = state[temperature_velocity_y];
    const Real w = state[temperature_velocity_z];
    const Real temperature = state[temperature_value];
    const Real sound = thermodynamic_sound_speed(state, gas, reference, floors, dimension);
    const Real speed = std::sqrt(u * u + v * v + w * w);
    const Real reference_speed
        = weiss_smith_reference_speed(speed, sound, viscous_speed, parameters);
    if (!std::isfinite(reference_speed) || reference_speed <= 0.0) {
        throw PhysicsError("Weiss-Smith reference speed is invalid");
    }
    const Real beta = (reference_speed / sound) * (reference_speed / sound);
    const Real rho_temperature = -rho / temperature;
    const Real cp = 1.0
        / ((gas.gamma() - 1.0) * reference.mach() * reference.mach());
    const Real enthalpy
        = thermodynamic_total_enthalpy(state, gas, reference, floors, dimension);
    const Real theta = 1.0 / (reference_speed * reference_speed)
        - rho_temperature / (rho * cp);

    WeissSmithState result;
    result.reference_speed = reference_speed;
    result.beta = beta;
    result.gamma = {{{{theta, 0.0, 0.0, 0.0, rho_temperature}},
                     {{theta * u, rho, 0.0, 0.0, rho_temperature * u}},
                     {{theta * v, 0.0, rho, 0.0, rho_temperature * v}},
                     {{theta * w, 0.0, 0.0, rho, rho_temperature * w}},
                     {{theta * enthalpy - 1.0,
                       rho * u,
                       rho * v,
                       rho * w,
                       rho_temperature * enthalpy + rho * cp}}}};
    result.inverse = invert(result.gamma);
    const Real normal_velocity
        = u * unit_normal.x + v * unit_normal.y + w * unit_normal.z;
    result.eigenvalues = weiss_smith_eigenvalues(normal_velocity, sound, beta);
    for (const auto value : result.eigenvalues) {
        if (!std::isfinite(value)) {
            throw PhysicsError("Weiss-Smith eigenvalue is non-finite");
        }
    }
    return result;
}

} // namespace wcns
