#include <wcns/solver/physical_boundary.hpp>
#include <wcns/physics/chapter5.hpp>

#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace wcns {
namespace {

bool finite(Real value)
{
    return std::isfinite(value);
}

void require_unit_normal(Normal3 normal)
{
    const Real norm = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!finite(norm) || std::abs(norm - 1.0) > 1.0e-12) {
        throw PhysicsConfigurationError("boundary face normal must be a unit vector");
    }
}

TemperaturePrimitiveState load_temperature(const Field<Real>& field, Index3 index)
{
    TemperaturePrimitiveState result {};
    for (int component = 0; component < fluid_components; ++component) {
        result[static_cast<std::size_t>(component)] = field(index.i, index.j, index.k, component);
    }
    return result;
}

void store_temperature(Field<Real>& field, Index3 index, const TemperaturePrimitiveState& state)
{
    for (int component = 0; component < fluid_components; ++component) {
        field(index.i, index.j, index.k, component) = state[static_cast<std::size_t>(component)];
    }
}

Index3 mirror_index(Index3 face, FaceLocation location, int layer, Extent3 extent)
{
    const auto axis = static_cast<std::size_t>(location.axis);
    face[axis] = location.side == Side::Lower ? layer - 1 : extent[axis] - layer;
    return face;
}

Index3 ghost_index(Index3 face, FaceLocation location, int layer, Extent3 extent)
{
    const auto axis = static_cast<std::size_t>(location.axis);
    face[axis] = location.side == Side::Lower ? -layer : extent[axis] - 1 + layer;
    return face;
}

TemperaturePrimitiveState reflected_velocity(TemperaturePrimitiveState state,
                                             Normal3 normal,
                                             const std::array<Real, 3>& wall_velocity,
                                             bool no_slip)
{
    if (no_slip) {
        state[temperature_velocity_x] = 2.0 * wall_velocity[0] - state[temperature_velocity_x];
        state[temperature_velocity_y] = 2.0 * wall_velocity[1] - state[temperature_velocity_y];
        state[temperature_velocity_z] = 2.0 * wall_velocity[2] - state[temperature_velocity_z];
        return state;
    }
    const Real relative_x = state[temperature_velocity_x] - wall_velocity[0];
    const Real relative_y = state[temperature_velocity_y] - wall_velocity[1];
    const Real relative_z = state[temperature_velocity_z] - wall_velocity[2];
    const Real normal_velocity
        = relative_x * normal.x + relative_y * normal.y + relative_z * normal.z;
    state[temperature_velocity_x] -= 2.0 * normal_velocity * normal.x;
    state[temperature_velocity_y] -= 2.0 * normal_velocity * normal.y;
    state[temperature_velocity_z] -= 2.0 * normal_velocity * normal.z;
    return state;
}

PressurePrimitiveState characteristic_boundary_state(const PressurePrimitiveState& interior,
                                                     const PressurePrimitiveState& target,
                                                     Normal3 outward_normal,
                                                     const GasModel& gas,
                                                     const NumericalFloors& floors,
                                                     int dimension)
{
    const Real rho = interior[0];
    const Real pressure = interior[4];
    const Real sound = std::sqrt(gas.gamma() * pressure / rho);
    const Real normal_velocity = interior[1] * outward_normal.x + interior[2] * outward_normal.y
        + interior[3] * outward_normal.z;
    if (normal_velocity <= -sound) return target;
    if (normal_velocity >= sound) return interior;

    const Real density_delta = target[0] - interior[0];
    const Real pressure_delta = target[4] - interior[4];
    const std::array<Real, 3> velocity_delta {{
        target[1] - interior[1],
        target[2] - interior[2],
        target[3] - interior[3],
    }};
    const Real normal_delta = velocity_delta[0] * outward_normal.x
        + velocity_delta[1] * outward_normal.y + velocity_delta[2] * outward_normal.z;
    Real acoustic_minus = 0.5 * (pressure_delta / (sound * sound) - rho * normal_delta / sound);
    Real acoustic_plus = 0.5 * (pressure_delta / (sound * sound) + rho * normal_delta / sound);
    Real entropy = density_delta - pressure_delta / (sound * sound);
    std::array<Real, 3> tangential {{
        velocity_delta[0] - normal_delta * outward_normal.x,
        velocity_delta[1] - normal_delta * outward_normal.y,
        velocity_delta[2] - normal_delta * outward_normal.z,
    }};
    if (normal_velocity - sound >= 0.0) acoustic_minus = 0.0;
    if (normal_velocity + sound >= 0.0) acoustic_plus = 0.0;
    if (normal_velocity >= 0.0) {
        entropy = 0.0;
        tangential = {{0.0, 0.0, 0.0}};
    }
    const Real selected_normal = sound * (acoustic_plus - acoustic_minus) / rho;
    PressurePrimitiveState result {
        interior[0] + acoustic_minus + acoustic_plus + entropy,
        interior[1] + selected_normal * outward_normal.x + tangential[0],
        interior[2] + selected_normal * outward_normal.y + tangential[1],
        interior[3] + selected_normal * outward_normal.z + tangential[2],
        interior[4] + sound * sound * (acoustic_minus + acoustic_plus),
    };
    if (dimension == 2) result[3] = 0.0;
    if (!finite(result[0]) || !finite(result[4]) || result[0] <= floors.density
        || result[4] <= floors.pressure) {
        throw PhysicsConfigurationError("characteristic boundary state is non-physical");
    }
    return result;
}

Normal3 patch_outward_normal(const StructuredBlock& block, const BoundaryPatch& patch, Index3 face)
{
    const FaceMetric* metric = nullptr;
    switch (patch.face.axis) {
    case Axis::I: metric = &block.face_metrics.i_faces; break;
    case Axis::J: metric = &block.face_metrics.j_faces; break;
    case Axis::K: metric = &block.face_metrics.k_faces; break;
    }
    const Real sign = patch.face.side == Side::Lower ? -1.0 : 1.0;
    Normal3 normal {
        sign * metric->normal_x(face.i, face.j, face.k),
        sign * metric->normal_y(face.i, face.j, face.k),
        sign * metric->normal_z(face.i, face.j, face.k),
    };
    const Real norm = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!finite(norm) || norm <= 0.0) {
        throw PhysicsConfigurationError("physical boundary has an invalid face normal");
    }
    normal.x /= norm;
    normal.y /= norm;
    normal.z /= norm;
    return normal;
}

TemperaturePrimitiveState wall_ghost(const TemperaturePrimitiveState& interior,
                                     Normal3 normal,
                                     const BoundaryData& data,
                                     BoundaryType type,
                                     const GasModel& gas,
                                     const ReferenceScales& reference,
                                     const NumericalFloors& floors,
                                     int dimension)
{
    const auto interior_pressure = pressure_primitive(interior, gas, reference, floors, dimension);
    const bool no_slip
        = type == BoundaryType::NoSlipAdiabaticWall || type == BoundaryType::NoSlipIsothermalWall;
    auto ghost = reflected_velocity(interior, normal, data.wall_velocity, no_slip);
    if (type == BoundaryType::NoSlipIsothermalWall) {
        ghost[temperature_value] = 2.0 * *data.wall_temperature - interior[temperature_value];
    } else {
        ghost[temperature_value] = interior[temperature_value];
    }
    ghost[temperature_density] = gas.gamma() * reference.mach() * reference.mach()
        * interior_pressure[4] / ghost[temperature_value];
    if (dimension == 2) {
        ghost[temperature_velocity_z] = 0.0;
    }
    static_cast<void>(pressure_primitive(ghost, gas, reference, floors, dimension));
    return ghost;
}

TemperaturePrimitiveState make_ghost(const TemperaturePrimitiveState& interior,
                                     Normal3 normal,
                                     const BoundaryData& data,
                                     const BoundaryPatch& patch,
                                     std::array<Real, 3> face_coordinates,
                                     Real time,
                                     const GasModel& gas,
                                     const ReferenceScales& reference,
                                     const NumericalFloors& floors,
                                     int dimension)
{
    switch (patch.type) {
    case BoundaryType::Farfield:
    case BoundaryType::Inflow:
    case BoundaryType::Outflow: {
        if (!data.target_state) return interior;
        const auto interior_pressure
            = pressure_primitive(interior, gas, reference, floors, dimension);
        const auto target_pressure = pressure_primitive(
            farfield_target_at(data, face_coordinates, dimension),
            gas,
            reference,
            floors,
            dimension);
        return temperature_primitive(
            characteristic_boundary_state(
                interior_pressure, target_pressure, normal, gas, floors, dimension),
            gas,
            reference,
            floors,
            dimension);
    }
    case BoundaryType::SlipWall:
    case BoundaryType::Symmetry:
    case BoundaryType::NoSlipAdiabaticWall:
    case BoundaryType::NoSlipIsothermalWall:
        return wall_ghost(interior, normal, data, patch.type, gas, reference, floors, dimension);
    case BoundaryType::DoubleMachReflection: {
        if (!data.double_mach_reflection) {
            throw PhysicsConfigurationError("double-Mach-reflection boundary data are missing");
        }
        const auto& model = *data.double_mach_reflection;
        model.validate(gas.gamma(), dimension);
        if (patch.face.axis == Axis::I && patch.face.side == Side::Lower) {
            return temperature_primitive(
                model.post_shock_state(), gas, reference, floors, dimension);
        }
        if (patch.face.axis == Axis::J && patch.face.side == Side::Upper) {
            return temperature_primitive(
                model.exact_state(face_coordinates[0], face_coordinates[1], time),
                gas,
                reference,
                floors,
                dimension);
        }
        if (patch.face.axis == Axis::J && patch.face.side == Side::Lower) {
            if (face_coordinates[0] < model.shock_foot()) {
                return temperature_primitive(
                    model.post_shock_state(), gas, reference, floors, dimension);
            }
            return wall_ghost(
                interior, normal, data, BoundaryType::SlipWall, gas, reference, floors, dimension);
        }
        throw PhysicsConfigurationError("double-Mach-reflection boundary is valid only on i-lower, "
                                        "j-lower or j-upper faces");
    }
    case BoundaryType::Periodic:
        throw PhysicsConfigurationError("periodic boundary must use a connectivity");
    case BoundaryType::Undefined:
        throw PhysicsConfigurationError("undefined physical boundary type");
    }
    throw PhysicsConfigurationError("unsupported physical boundary type");
}

PressurePrimitiveState reflected_face_trace(PressurePrimitiveState state,
                                            Normal3 normal,
                                            const std::array<Real, 3>& wall_velocity)
{
    const Real relative_x = state[1] - wall_velocity[0];
    const Real relative_y = state[2] - wall_velocity[1];
    const Real relative_z = state[3] - wall_velocity[2];
    const Real normal_velocity
        = relative_x * normal.x + relative_y * normal.y + relative_z * normal.z;
    state[1] -= 2.0 * normal_velocity * normal.x;
    state[2] -= 2.0 * normal_velocity * normal.y;
    state[3] -= 2.0 * normal_velocity * normal.z;
    return state;
}

} // namespace

void FarfieldPointVortex::validate() const
{
    if (!finite(lift_coefficient) || !finite(center[0]) || !finite(center[1])
        || !finite(chord) || chord <= 0.0) {
        throw PhysicsConfigurationError("farfield point vortex parameters are invalid");
    }
}

TemperaturePrimitiveState farfield_target_at(const BoundaryData& data,
                                              std::array<Real, 3> face_coordinates,
                                              int dimension)
{
    if (!data.target_state) {
        throw PhysicsConfigurationError("farfield target state is missing");
    }
    if (data.compression_ramp_inlet) return ramp_inlet(face_coordinates[1]);
    auto target = *data.target_state;
    if (!data.farfield_point_vortex) return target;
    if (dimension != 2) {
        throw PhysicsConfigurationError("farfield point vortex is available only in 2-D");
    }
    const auto& vortex = *data.farfield_point_vortex;
    vortex.validate();
    const Real dx = face_coordinates[0] - vortex.center[0];
    const Real dy = face_coordinates[1] - vortex.center[1];
    const Real radius_squared = dx * dx + dy * dy;
    if (!finite(radius_squared) || radius_squared <= std::numeric_limits<Real>::epsilon()) {
        throw PhysicsConfigurationError("farfield point vortex is singular at a boundary face");
    }
    const Real speed = std::hypot(target[temperature_velocity_x],
                                  target[temperature_velocity_y]);
    const Real circulation = 0.5 * vortex.lift_coefficient * speed * vortex.chord;
    const Real factor = circulation / (2.0 * std::acos(-1.0) * radius_squared);
    // Positive lift uses the aerodynamic (clockwise) circulation convention.
    target[temperature_velocity_x] += factor * dy;
    target[temperature_velocity_y] -= factor * dx;
    return target;
}

void BoundaryData::validate(BoundaryType type, int dimension) const
{
    if (dimension != 2 && dimension != 3) {
        throw PhysicsConfigurationError("boundary data dimension must be two or three");
    }
    for (const Real value : wall_velocity) {
        if (!finite(value)) {
            throw PhysicsConfigurationError("wall velocity must be finite");
        }
    }
    if (dimension == 2 && wall_velocity[2] != 0.0) {
        throw PhysicsConfigurationError("two-dimensional wall z velocity must be zero");
    }
    if(compression_ramp_inlet && (type!=BoundaryType::Inflow || dimension!=3))
        throw std::invalid_argument("ramp profile requires 3D inflow");
    const bool needs_target = type == BoundaryType::Farfield || type == BoundaryType::Inflow;
    const bool permits_target = needs_target || type == BoundaryType::Outflow;
    if ((needs_target && !target_state) || (!permits_target && target_state)) {
        throw PhysicsConfigurationError(
            needs_target ? "boundary target state is required"
                         : "boundary target state is not valid for this boundary type");
    }
    const bool needs_temperature = type == BoundaryType::NoSlipIsothermalWall;
    if (needs_temperature != wall_temperature.has_value()) {
        throw PhysicsConfigurationError(
            needs_temperature ? "isothermal wall temperature is required"
                              : "wall temperature is only valid for an isothermal wall");
    }
    if (wall_temperature && (!finite(*wall_temperature) || *wall_temperature <= 0.0)) {
        throw PhysicsConfigurationError("wall temperature must be positive and finite");
    }
    const bool needs_double_mach = type == BoundaryType::DoubleMachReflection;
    if (needs_double_mach != double_mach_reflection.has_value()) {
        throw PhysicsConfigurationError(
            needs_double_mach
                ? "double-Mach-reflection model data are required"
                : "double-Mach-reflection model data are invalid for this boundary type");
    }
    if (needs_double_mach) {
        if (wall_velocity != std::array<Real, 3> {{0.0, 0.0, 0.0}}) {
            throw PhysicsConfigurationError(
                "double-Mach-reflection boundary requires a stationary wall");
        }
    }
    if (farfield_point_vortex) {
        if (type != BoundaryType::Farfield) {
            throw PhysicsConfigurationError(
                "farfield point vortex is only valid on a farfield boundary");
        }
        if (dimension != 2) {
            throw PhysicsConfigurationError("farfield point vortex is available only in 2-D");
        }
        farfield_point_vortex->validate();
    }
}

std::array<Real, 3>
boundary_face_coordinates(const StructuredBlock& block, const BoundaryPatch& patch, Index3 face)
{
    const int dimension = block.cell_dimension();
    const int tangential_count = 1 << (dimension - 1);
    std::array<Real, 3> result {{0.0, 0.0, 0.0}};
    for (int mask = 0; mask < tangential_count; ++mask) {
        Index3 vertex = face;
        int bit = 0;
        for (int axis = 0; axis < dimension; ++axis) {
            if (axis == static_cast<int>(patch.face.axis)) continue;
            vertex[static_cast<std::size_t>(axis)] += (mask >> bit) & 1;
            ++bit;
        }
        result[0] += block.coordinates.x(vertex.i, vertex.j, vertex.k);
        result[1] += block.coordinates.y(vertex.i, vertex.j, vertex.k);
        result[2] += block.coordinates.z(vertex.i, vertex.j, vertex.k);
    }
    for (auto& value : result)
        value /= static_cast<Real>(tangential_count);
    return result;
}

PhysicalGhostFillResult PhysicalGhostStateOperator::fill(StructuredBlock& block,
                                                         const BoundaryDataMap& boundary_data,
                                                         const GasModel& gas,
                                                         const ReferenceScales& reference,
                                                         const NumericalFloors& floors,
                                                         std::uint64_t version,
                                                         Real time)
{
    if (version == 0) {
        throw PhysicsConfigurationError("physical ghost version must be non-zero");
    }
    if (block.ghost_width() < 3) {
        throw PhysicsConfigurationError("WCNS physical ghost fill requires three layers");
    }
    if (!finite(time) || time < 0.0) {
        throw PhysicsConfigurationError("physical boundary time must be finite and non-negative");
    }
    PhysicalGhostFillResult result {version, 0};
    for (const auto& patch : block.boundaries) {
        const auto data_iterator = boundary_data.find(patch.name);
        if (data_iterator == boundary_data.end()) {
            throw PhysicsConfigurationError("physical boundary data is missing for patch "
                                            + patch.name);
        }
        const auto& data = data_iterator->second;
        data.validate(patch.type, block.cell_dimension());
        const auto counts = patch.boundary_face_range.counts();
        for (int ok = 0; ok < counts.nk; ++ok) {
            for (int oj = 0; oj < counts.nj; ++oj) {
                for (int oi = 0; oi < counts.ni; ++oi) {
                    const auto face = patch.boundary_face_range.at({oi, oj, ok});
                    const auto normal = patch_outward_normal(block, patch, face);
                    const auto coordinates = boundary_face_coordinates(block, patch, face);
                    for (int layer = 1; layer <= 3; ++layer) {
                        const auto interior_index
                            = mirror_index(face, patch.face, layer, block.cell_extent());
                        const auto ghost = make_ghost(
                            load_temperature(block.flow.temperature_primitive, interior_index),
                            normal,
                            data,
                            patch,
                            coordinates,
                            time,
                            gas,
                            reference,
                            floors,
                            block.cell_dimension());
                        const auto destination
                            = ghost_index(face, patch.face, layer, block.cell_extent());
                        const auto pressure = pressure_primitive(
                            ghost, gas, reference, floors, block.cell_dimension());
                        const auto conservative = thermodynamic_conservative(
                            ghost, gas, reference, floors, block.cell_dimension());
                        store_temperature(block.flow.temperature_primitive, destination, ghost);
                        store_state(block.flow.primitive, destination, pressure);
                        store_state(block.flow.conservative, destination, conservative);
                        ++result.state_count;
                    }
                }
            }
        }
    }
    block.flow.physical_ghost_version = version;
    return result;
}

std::string InviscidBoundaryOptions::summary() const
{
    return std::string("strong_boundary_face_state=")
        + (strong_boundary_face_state ? "true" : "false");
}

std::string InviscidBoundaryOptions::restart_signature() const
{
    return "inviscid_boundary_v1;" + summary();
}

PressurePrimitiveState
apply_inviscid_boundary_face_state(const BoundaryPatch& patch,
                                   const PressurePrimitiveState& interior_trace,
                                   const PressurePrimitiveState& reconstructed_exterior_trace,
                                   Normal3 outward_unit_normal,
                                   const BoundaryData& data,
                                   const InviscidBoundaryOptions& options,
                                   const GasModel& gas,
                                   const ReferenceScales& reference,
                                   const NumericalFloors& floors,
                                   int dimension,
                                   std::array<Real, 3> face_coordinates,
                                   Real time)
{
    require_unit_normal(outward_unit_normal);
    data.validate(patch.type, dimension);
    static_cast<void>(temperature_primitive(interior_trace, gas, reference, floors, dimension));
    static_cast<void>(
        temperature_primitive(reconstructed_exterior_trace, gas, reference, floors, dimension));
    if (!options.strong_boundary_face_state) {
        return reconstructed_exterior_trace;
    }
    switch (patch.type) {
    case BoundaryType::SlipWall:
    case BoundaryType::Symmetry:
    case BoundaryType::NoSlipAdiabaticWall:
    case BoundaryType::NoSlipIsothermalWall:
        return reflected_face_trace(interior_trace, outward_unit_normal, data.wall_velocity);
    case BoundaryType::Farfield:
    case BoundaryType::Inflow:
    case BoundaryType::Outflow:
        if (!data.target_state) return interior_trace;
        return characteristic_boundary_state(
            interior_trace,
            pressure_primitive(farfield_target_at(data, face_coordinates, dimension),
                               gas,
                               reference,
                               floors,
                               dimension),
            outward_unit_normal,
            gas,
            floors,
            dimension);
    case BoundaryType::DoubleMachReflection: {
        const auto& model = *data.double_mach_reflection;
        model.validate(gas.gamma(), dimension);
        if (patch.face.axis == Axis::I && patch.face.side == Side::Lower) {
            return model.post_shock_state();
        }
        if (patch.face.axis == Axis::J && patch.face.side == Side::Upper) {
            return model.exact_state(face_coordinates[0], face_coordinates[1], time);
        }
        if (patch.face.axis == Axis::J && patch.face.side == Side::Lower) {
            if (face_coordinates[0] < model.shock_foot()) {
                return model.post_shock_state();
            }
            return reflected_face_trace(interior_trace, outward_unit_normal, data.wall_velocity);
        }
        throw PhysicsConfigurationError("double-Mach-reflection boundary is valid only on i-lower, "
                                        "j-lower or j-upper faces");
    }
    case BoundaryType::Periodic:
        throw PhysicsConfigurationError("periodic boundary must use a connectivity");
    case BoundaryType::Undefined:
        throw PhysicsConfigurationError("undefined physical boundary type");
    }
    throw PhysicsConfigurationError("unsupported inviscid boundary type");
}

void update_temperature_primitive_interior(StructuredBlock& block,
                                           const GasModel& gas,
                                           const ReferenceScales& reference,
                                           const NumericalFloors& floors)
{
    const auto extent = block.cell_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                update_temperature_primitive_cell(block, {i, j, k}, gas, reference, floors);
            }
        }
    }
}

void update_temperature_primitive_cell(StructuredBlock& block,
                                       Index3 index,
                                       const GasModel& gas,
                                       const ReferenceScales& reference,
                                       const NumericalFloors& floors)
{
    const auto conservative = load_conservative(block.flow.conservative, index);
    try {
        const auto temperature = temperature_primitive_from_conservative(
            conservative, gas, reference, floors, block.cell_dimension());
        const auto pressure
            = pressure_primitive(temperature, gas, reference, floors, block.cell_dimension());
        store_temperature(block.flow.temperature_primitive, index, temperature);
        store_state(block.flow.primitive, index, pressure);
    } catch (const PhysicsConfigurationError& error) {
        std::ostringstream message;
        message << std::setprecision(17) << "block=" << block.id() << " cell=(" << index.i << ','
                << index.j << ',' << index.k << ") conservative=(";
        for (std::size_t component = 0; component < conservative.size(); ++component) {
            if (component != 0) message << ',';
            message << conservative[component];
        }
        message << "): " << error.what();
        throw PhysicsError(message.str());
    }
}

} // namespace wcns
