#include <wcns/physics/thermodynamics.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace wcns {
namespace {

bool positive_finite(Real value)
{
    return std::isfinite(value) && value > 0.0;
}

bool nonnegative_finite(Real value)
{
    return std::isfinite(value) && value >= 0.0;
}

void require_dimension(int dimension)
{
    if (dimension != 2 && dimension != 3) {
        throw PhysicsConfigurationError("fluid state dimension must be two or three");
    }
}

template <class State> void require_finite_state(const State& state, const char* label)
{
    for (const auto value : state) {
        if (!std::isfinite(value)) {
            throw PhysicsConfigurationError(std::string(label) + " contains a non-finite value");
        }
    }
}

void require_two_dimensional_velocity(Real velocity_z, int dimension)
{
    require_dimension(dimension);
    if (dimension == 2 && velocity_z != 0.0) {
        throw PhysicsConfigurationError("two-dimensional state must have zero z velocity");
    }
}

std::string number(Real value)
{
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<Real>::max_digits10) << value;
    return stream.str();
}

void validate_temperature_state(const TemperaturePrimitiveState& state,
                                const NumericalFloors& floors,
                                int dimension)
{
    floors.validate();
    require_finite_state(state, "temperature primitive state");
    require_two_dimensional_velocity(state[temperature_velocity_z], dimension);
    if (state[temperature_density] <= floors.density) {
        throw PhysicsConfigurationError("temperature primitive density is below its floor");
    }
    if (state[temperature_value] <= floors.temperature) {
        throw PhysicsConfigurationError("temperature primitive temperature is below its floor");
    }
}

void validate_pressure_state(const PressurePrimitiveState& state,
                             const NumericalFloors& floors,
                             int dimension)
{
    floors.validate();
    require_finite_state(state, "pressure primitive state");
    require_two_dimensional_velocity(state[3], dimension);
    if (state[0] <= floors.density) {
        throw PhysicsConfigurationError("pressure primitive density is below its floor");
    }
    if (state[4] <= floors.pressure) {
        throw PhysicsConfigurationError("pressure primitive pressure is below its floor");
    }
}

} // namespace

GasModel GasModel::from_input(const GasModelInput& input)
{
    if (!positive_finite(input.gamma) || input.gamma <= 1.0) {
        throw PhysicsConfigurationError("gas gamma must be finite and greater than one");
    }
    if (input.molar_mass.has_value() == input.specific_gas_constant.has_value()) {
        throw PhysicsConfigurationError(
            "exactly one of molar_mass and specific_gas_constant is required");
    }
    if (input.molar_mass) {
        if (!positive_finite(*input.molar_mass)) {
            throw PhysicsConfigurationError("gas molar mass must be positive and finite");
        }
        return {input.gamma, *input.molar_mass, universal_gas_constant / *input.molar_mass};
    }
    if (!positive_finite(*input.specific_gas_constant)) {
        throw PhysicsConfigurationError("specific gas constant must be positive and finite");
    }
    return {
        input.gamma,
        universal_gas_constant / *input.specific_gas_constant,
        *input.specific_gas_constant,
    };
}

std::string GasModel::summary() const
{
    return "gas_model=calorically_perfect gamma=" + number(gamma_) + " molar_mass="
        + number(molar_mass_) + " specific_gas_constant=" + number(specific_gas_constant_);
}

std::string GasModel::restart_signature() const
{
    return summary();
}

ReferenceScales ReferenceScales::derive(const ReferenceInput& input, const GasModel& gas)
{
    if (input.reynolds || input.mach) {
        throw PhysicsConfigurationError("Re and Ma are derived values and must not be provided");
    }
    const std::array<Real, 5> values {{
        input.velocity,
        input.density,
        input.temperature,
        input.length,
        input.viscosity,
    }};
    if (!std::all_of(values.begin(), values.end(), positive_finite)) {
        throw PhysicsConfigurationError("all five reference scales must be positive and finite");
    }
    const Real sound = std::sqrt(gas.gamma() * gas.specific_gas_constant() * input.temperature);
    const Real reynolds = input.density * input.velocity * input.length / input.viscosity;
    const Real mach = input.velocity / sound;
    const Real dynamic_pressure = input.density * input.velocity * input.velocity;
    const Real time = input.length / input.velocity;
    if (!positive_finite(sound) || !positive_finite(reynolds) || !positive_finite(mach)
        || !positive_finite(dynamic_pressure) || !positive_finite(time)) {
        throw PhysicsConfigurationError("derived reference values are not positive and finite");
    }
    return {
        input.velocity,
        input.density,
        input.temperature,
        input.length,
        input.viscosity,
        reynolds,
        mach,
        dynamic_pressure,
        time,
    };
}

std::string ReferenceScales::summary() const
{
    return "U_ref=" + number(velocity_) + " rho_ref=" + number(density_)
        + " T_ref=" + number(temperature_) + " L_ref=" + number(length_)
        + " mu_ref=" + number(viscosity_) + " Re=" + number(reynolds_) + " Ma=" + number(mach_)
        + " dynamic_pressure_ref=" + number(dynamic_pressure_) + " time_ref=" + number(time_);
}

std::string ReferenceScales::restart_signature() const
{
    return summary();
}

void NumericalFloors::validate() const
{
    if (!positive_finite(density) || !positive_finite(pressure) || !positive_finite(temperature)
        || !nonnegative_finite(jacobian_absolute) || !nonnegative_finite(jacobian_relative)
        || !nonnegative_finite(face_area_absolute) || !nonnegative_finite(face_area_relative)
        || !positive_finite(reconstruction_scale) || !positive_finite(reconstruction_epsilon)) {
        throw PhysicsConfigurationError("numerical floors are invalid");
    }
}

Real NumericalFloors::jacobian_floor(Real reference_volume) const
{
    validate();
    if (!positive_finite(reference_volume)) {
        throw PhysicsConfigurationError("reference volume must be positive and finite");
    }
    return std::max(jacobian_absolute, jacobian_relative * reference_volume);
}

Real NumericalFloors::face_area_floor(Real reference_area) const
{
    validate();
    if (!positive_finite(reference_area)) {
        throw PhysicsConfigurationError("reference area must be positive and finite");
    }
    return std::max(face_area_absolute, face_area_relative * reference_area);
}

void ReconstructionScaling::validate() const
{
    if (!positive_finite(scale_floor) || !positive_finite(epsilon)
        || !std::all_of(component.begin(), component.end(), positive_finite)) {
        throw PhysicsConfigurationError("reconstruction scaling is invalid");
    }
}

Real ReconstructionScaling::normalized_smoothness(Real beta, int component_index) const
{
    validate();
    if (!nonnegative_finite(beta) || component_index < 0 || component_index >= fluid_components) {
        throw PhysicsConfigurationError("smoothness input is invalid");
    }
    const Real scale = std::max(component[static_cast<std::size_t>(component_index)], scale_floor);
    return beta / (scale * scale);
}

PressurePrimitiveState pressure_primitive(const TemperaturePrimitiveState& state,
                                          const GasModel& gas,
                                          const ReferenceScales& reference,
                                          const NumericalFloors& floors,
                                          int dimension)
{
    validate_temperature_state(state, floors, dimension);
    const Real pressure = state[temperature_density] * state[temperature_value]
        / (gas.gamma() * reference.mach() * reference.mach());
    if (!std::isfinite(pressure) || pressure <= floors.pressure) {
        throw PhysicsConfigurationError("derived pressure is below its floor");
    }
    return {
        state[temperature_density],
        state[temperature_velocity_x],
        state[temperature_velocity_y],
        state[temperature_velocity_z],
        pressure,
    };
}

TemperaturePrimitiveState temperature_primitive(const PressurePrimitiveState& state,
                                                const GasModel& gas,
                                                const ReferenceScales& reference,
                                                const NumericalFloors& floors,
                                                int dimension)
{
    validate_pressure_state(state, floors, dimension);
    const Real temperature
        = gas.gamma() * reference.mach() * reference.mach() * state[4] / state[0];
    if (!std::isfinite(temperature) || temperature <= floors.temperature) {
        throw PhysicsConfigurationError("derived temperature is below its floor");
    }
    return {state[0], state[1], state[2], state[3], temperature};
}

ThermodynamicConservativeState thermodynamic_conservative(const TemperaturePrimitiveState& state,
                                                          const GasModel& gas,
                                                          const ReferenceScales& reference,
                                                          const NumericalFloors& floors,
                                                          int dimension)
{
    const auto primitive = pressure_primitive(state, gas, reference, floors, dimension);
    const Real kinetic = 0.5 * primitive[0]
        * (primitive[1] * primitive[1] + primitive[2] * primitive[2] + primitive[3] * primitive[3]);
    return {
        primitive[0],
        primitive[0] * primitive[1],
        primitive[0] * primitive[2],
        primitive[0] * primitive[3],
        primitive[4] / (gas.gamma() - 1.0) + kinetic,
    };
}

TemperaturePrimitiveState
temperature_primitive_from_conservative(const ThermodynamicConservativeState& state,
                                        const GasModel& gas,
                                        const ReferenceScales& reference,
                                        const NumericalFloors& floors,
                                        int dimension)
{
    floors.validate();
    require_finite_state(state, "thermodynamic conservative state");
    const Real rho = state[0];
    if (rho <= floors.density) {
        throw PhysicsConfigurationError("conservative density is below its floor");
    }
    const Real u = state[1] / rho;
    const Real v = state[2] / rho;
    const Real w = state[3] / rho;
    require_two_dimensional_velocity(w, dimension);
    const Real kinetic = 0.5 * rho * (u * u + v * v + w * w);
    const Real pressure = (gas.gamma() - 1.0) * (state[4] - kinetic);
    return temperature_primitive(
        PressurePrimitiveState {rho, u, v, w, pressure}, gas, reference, floors, dimension);
}

Real thermodynamic_sound_speed(const TemperaturePrimitiveState& state,
                               const GasModel& gas,
                               const ReferenceScales& reference,
                               const NumericalFloors& floors,
                               int dimension)
{
    const auto primitive = pressure_primitive(state, gas, reference, floors, dimension);
    return std::sqrt(gas.gamma() * primitive[4] / primitive[0]);
}

Real thermodynamic_total_enthalpy(const TemperaturePrimitiveState& state,
                                  const GasModel& gas,
                                  const ReferenceScales& reference,
                                  const NumericalFloors& floors,
                                  int dimension)
{
    const auto primitive = pressure_primitive(state, gas, reference, floors, dimension);
    const Real velocity_squared = primitive[1] * primitive[1] + primitive[2] * primitive[2]
        + primitive[3] * primitive[3];
    // The state has already been validated by pressure_primitive.  Evaluating
    // H directly avoids converting and validating the same state a second time.
    return gas.gamma() * primitive[4] / ((gas.gamma() - 1.0) * primitive[0])
        + 0.5 * velocity_squared;
}

} // namespace wcns
