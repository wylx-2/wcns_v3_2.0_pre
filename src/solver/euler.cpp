#include <wcns/solver/euler.hpp>

#include <algorithm>
#include <cmath>
#include <string>

namespace wcns {
namespace {

void validate_finite(const std::array<Real, euler_components>& state, const char* label)
{
    for (const auto value : state) {
        if (!std::isfinite(value)) {
            throw PhysicsError(std::string(label) + " contains a non-finite value");
        }
    }
}

void validate_primitive(const PrimitiveState& primitive, const IdealGas& gas)
{
    validate_finite(primitive, "primitive state");
    if (primitive[primitive_density] <= gas.density_floor) {
        throw PhysicsError("primitive density is below its floor");
    }
    if (primitive[pressure] <= gas.pressure_floor) {
        throw PhysicsError("primitive pressure is below its floor");
    }
}

Normal3 checked_normal(Normal3 normal)
{
    const Real magnitude
        = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!std::isfinite(magnitude) || magnitude <= 0.0) {
        throw PhysicsError("Euler flux normal must be finite and non-zero");
    }
    normal.x /= magnitude;
    normal.y /= magnitude;
    normal.z /= magnitude;
    return normal;
}

Real normal_velocity(const PrimitiveState& state, Normal3 normal)
{
    return state[velocity_x] * normal.x + state[velocity_y] * normal.y
        + state[velocity_z] * normal.z;
}

ConservativeState to_conservative_unchecked(const PrimitiveState& primitive,
                                            const IdealGas& gas)
{
    const Real rho = primitive[primitive_density];
    const Real u = primitive[velocity_x];
    const Real v = primitive[velocity_y];
    const Real w = primitive[velocity_z];
    const Real kinetic = 0.5 * rho * (u * u + v * v + w * w);
    return {
        rho,
        rho * u,
        rho * v,
        rho * w,
        primitive[pressure] / (gas.gamma - 1.0) + kinetic,
    };
}

ConservativeState euler_flux_unchecked(const PrimitiveState& primitive,
                                       const ConservativeState& conservative,
                                       Normal3 normal)
{
    const Real rho = primitive[primitive_density];
    const Real pressure_value = primitive[pressure];
    const Real un = normal_velocity(primitive, normal);
    return {
        rho * un,
        conservative[momentum_x] * un + pressure_value * normal.x,
        conservative[momentum_y] * un + pressure_value * normal.y,
        conservative[momentum_z] * un + pressure_value * normal.z,
        (conservative[total_energy] + pressure_value) * un,
    };
}

Real sound_speed_unchecked(const PrimitiveState& primitive, const IdealGas& gas)
{
    return std::sqrt(gas.gamma * primitive[pressure] / primitive[primitive_density]);
}

} // namespace

void IdealGas::validate() const
{
    if (!std::isfinite(gamma) || gamma <= 1.0 || !std::isfinite(density_floor)
        || density_floor < 0.0 || !std::isfinite(pressure_floor) || pressure_floor < 0.0) {
        throw PhysicsError("ideal-gas parameters are invalid");
    }
}

ConservativeState to_conservative(const PrimitiveState& primitive, const IdealGas& gas)
{
    gas.validate();
    validate_primitive(primitive, gas);
    return to_conservative_unchecked(primitive, gas);
}

PrimitiveState to_primitive(const ConservativeState& conservative, const IdealGas& gas)
{
    gas.validate();
    validate_finite(conservative, "conservative state");
    const Real rho = conservative[density];
    if (rho <= gas.density_floor) {
        throw PhysicsError("conservative density is below its floor");
    }
    const Real inverse_density = 1.0 / rho;
    const Real u = conservative[momentum_x] * inverse_density;
    const Real v = conservative[momentum_y] * inverse_density;
    const Real w = conservative[momentum_z] * inverse_density;
    const Real kinetic = 0.5 * rho * (u * u + v * v + w * w);
    const Real p = (gas.gamma - 1.0) * (conservative[total_energy] - kinetic);
    PrimitiveState result {rho, u, v, w, p};
    validate_primitive(result, gas);
    return result;
}

Real sound_speed(const PrimitiveState& primitive, const IdealGas& gas)
{
    gas.validate();
    validate_primitive(primitive, gas);
    return sound_speed_unchecked(primitive, gas);
}

ConservativeState
euler_flux(const PrimitiveState& primitive, Normal3 unit_normal, const IdealGas& gas)
{
    gas.validate();
    validate_primitive(primitive, gas);
    const auto normal = checked_normal(unit_normal);
    return euler_flux_unchecked(primitive, to_conservative_unchecked(primitive, gas), normal);
}

ConservativeState rusanov_flux(const PrimitiveState& left,
                               const PrimitiveState& right,
                               Normal3 unit_normal,
                               const IdealGas& gas)
{
    gas.validate();
    const auto normal = checked_normal(unit_normal);
    validate_primitive(left, gas);
    validate_primitive(right, gas);
    const auto left_conservative = to_conservative_unchecked(left, gas);
    const auto right_conservative = to_conservative_unchecked(right, gas);
    const auto left_flux = euler_flux_unchecked(left, left_conservative, normal);
    const auto right_flux = euler_flux_unchecked(right, right_conservative, normal);
    const Real speed
        = std::max(std::abs(normal_velocity(left, normal)) + sound_speed_unchecked(left, gas),
                   std::abs(normal_velocity(right, normal)) + sound_speed_unchecked(right, gas));
    ConservativeState result {};
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        result[index] = 0.5 * (left_flux[index] + right_flux[index])
            - 0.5 * speed * (right_conservative[index] - left_conservative[index]);
    }
    return result;
}

ConservativeState load_conservative(const Field<Real>& field, Index3 index)
{
    if (field.components() != euler_components) {
        throw std::invalid_argument("Euler field must have five components");
    }
    ConservativeState result {};
    for (int component = 0; component < euler_components; ++component) {
        result[static_cast<std::size_t>(component)] = field(index.i, index.j, index.k, component);
    }
    return result;
}

PrimitiveState load_primitive(const Field<Real>& field, Index3 index)
{
    return load_conservative(field, index);
}

void store_state(Field<Real>& field, Index3 index, const ConservativeState& state)
{
    if (field.components() != euler_components) {
        throw std::invalid_argument("Euler field must have five components");
    }
    for (int component = 0; component < euler_components; ++component) {
        field(index.i, index.j, index.k, component) = state[static_cast<std::size_t>(component)];
    }
}

} // namespace wcns
