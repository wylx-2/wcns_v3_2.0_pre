#include <wcns/solver/viscous_boundary.hpp>

#include <wcns/mesh/linear_operators.hpp>

#include <array>
#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

std::array<Real, 3> as_array(Normal3 normal)
{
    return {{normal.x, normal.y, normal.z}};
}

Real dot(const std::array<Real, 3>& lhs, const std::array<Real, 3>& rhs)
{
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

void require_unit_normal(Normal3 normal)
{
    const Real norm = std::sqrt(normal.x * normal.x + normal.y * normal.y + normal.z * normal.z);
    if (!std::isfinite(norm) || std::abs(norm - 1.0) > 1.0e-12) {
        throw PhysicsConfigurationError("viscous boundary requires a unit normal");
    }
}

std::array<Real, 3>
face_centroid(const StructuredBlock& block, const BoundaryPatch& patch, Index3 face)
{
    std::array<Real, 3> result {{0.0, 0.0, 0.0}};
    int count = 0;
    const auto add = [&](Index3 vertex) {
        result[0] += block.coordinates.x(vertex.i, vertex.j, vertex.k);
        result[1] += block.coordinates.y(vertex.i, vertex.j, vertex.k);
        result[2] += block.coordinates.z(vertex.i, vertex.j, vertex.k);
        ++count;
    };
    if (block.cell_dimension() == 2) {
        if (patch.face.axis == Axis::I) {
            add({face.i, face.j, 0});
            add({face.i, face.j + 1, 0});
        } else {
            add({face.i, face.j, 0});
            add({face.i + 1, face.j, 0});
        }
    } else {
        const int axis = static_cast<int>(patch.face.axis);
        const int tangent_a = (axis + 1) % 3;
        const int tangent_b = (axis + 2) % 3;
        for (int a = 0; a <= 1; ++a) {
            for (int b = 0; b <= 1; ++b) {
                auto vertex = face;
                vertex[static_cast<std::size_t>(tangent_a)] += a;
                vertex[static_cast<std::size_t>(tangent_b)] += b;
                add(vertex);
            }
        }
    }
    for (auto& value : result)
        value /= static_cast<Real>(count);
    return result;
}

Index3
inward_cell(const StructuredBlock& block, const BoundaryPatch& patch, Index3 face, int ordinal)
{
    auto result = face;
    const auto axis = static_cast<std::size_t>(patch.face.axis);
    result[axis]
        = patch.face.side == Side::Lower ? ordinal : block.cell_extent()[axis] - 1 - ordinal;
    return result;
}

Real normal_physical_scale(const StructuredBlock& block,
                           const MetricField& metric,
                           const BoundaryPatch& patch,
                           Index3 face,
                           Normal3 outward,const AlgorithmProfile& profile)
{
    const auto center = inward_cell(block, patch, face, 0);
    const auto& coordinates = metric.cell_coordinates();
    const std::array<Real, 3> cell {{
        coordinates.x(center.i, center.j, center.k),
        coordinates.y(center.i, center.j, center.k),
        coordinates.z(center.i, center.j, center.k),
    }};
    const auto boundary = viscous_wall_coordinates(block, metric, patch, face, profile);
    const auto outward_array = as_array(outward);
    std::array<Real, 3> inward {{-outward_array[0], -outward_array[1], -outward_array[2]}};
    const std::array<Real, 3> delta {
        {cell[0] - boundary[0], cell[1] - boundary[1], cell[2] - boundary[2]}};
    const Real distance = dot(delta, inward);
    if (!std::isfinite(distance) || distance <= 0.0) {
        throw PhysicsError("wall first cell center is not inside the boundary");
    }
    return 0.5 / distance;
}

std::vector<Real> inward_values(
    const StructuredBlock& block, const BoundaryPatch& patch, Index3 face, int component, int count)
{
    std::vector<Real> result;
    result.reserve(static_cast<std::size_t>(count));
    for (int ordinal = 0; ordinal < count; ++ordinal) {
        const auto cell = inward_cell(block, patch, face, ordinal);
        result.push_back(block.flow.temperature_primitive(cell.i, cell.j, cell.k, component));
    }
    return result;
}

CartesianGradient normal_gradient(Real inward_derivative, Normal3 outward)
{
    return {{
        -inward_derivative * outward.x,
        -inward_derivative * outward.y,
        -inward_derivative * outward.z,
    }};
}

CartesianGradient remove_normal(CartesianGradient gradient, Normal3 outward)
{
    const Real component
        = gradient[0] * outward.x + gradient[1] * outward.y + gradient[2] * outward.z;
    gradient[0] -= component * outward.x;
    gradient[1] -= component * outward.y;
    gradient[2] -= component * outward.z;
    return gradient;
}

} // namespace

std::array<Real,3> viscous_wall_coordinates(const StructuredBlock& block,const MetricField& metric,
    const BoundaryPatch& patch,Index3 face,const AlgorithmProfile& profile)
{
    if(profile.kind()!=AlgorithmProfileKind::Scmm6Wcns) return face_centroid(block,patch,face);
    const auto axis=static_cast<std::size_t>(patch.face.axis);
    const auto& operators=cached_line_operators(profile,block.cell_extent()[axis]);
    const auto& row=operators.interpolation_rows()[static_cast<std::size_t>(face[axis])];
    const auto& c=metric.cell_coordinates();std::array<Real,3> result{};
    for(const auto [j,w]:row) {
        auto p=face;p[axis]=j;
        result[0]+=w*c.x(p.i,p.j,p.k);result[1]+=w*c.y(p.i,p.j,p.k);result[2]+=w*c.z(p.i,p.j,p.k);
    }
    return result;
}

Real wall_dirichlet_computational_derivative(Real wall_value,
                                             const std::vector<Real>& inward_center_values,
                                             const AlgorithmProfile& profile)
{
    const int required = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 4 : 6;
    if (static_cast<int>(inward_center_values.size()) != required || !std::isfinite(wall_value)) {
        throw std::invalid_argument("wall Dirichlet stencil has the wrong size or value");
    }
    std::vector<Real> nodes(static_cast<std::size_t>(required + 1));
    std::vector<Real> values(static_cast<std::size_t>(required + 1));
    nodes[0] = 0.0;
    values[0] = wall_value;
    for (int index = 0; index < required; ++index) {
        nodes[static_cast<std::size_t>(index + 1)] = static_cast<Real>(index) + 0.5;
        values[static_cast<std::size_t>(index + 1)]
            = inward_center_values[static_cast<std::size_t>(index)];
        if (!std::isfinite(values[static_cast<std::size_t>(index + 1)])) {
            throw PhysicsError("wall Dirichlet stencil contains a non-finite value");
        }
    }
    Real derivative = 0.0;
    for (int basis = 0; basis <= required; ++basis) {
        Real basis_derivative = 0.0;
        for (int differentiated = 0; differentiated <= required; ++differentiated) {
            if (differentiated == basis) continue;
            Real term = 1.0
                / (nodes[static_cast<std::size_t>(basis)]
                   - nodes[static_cast<std::size_t>(differentiated)]);
            for (int factor = 0; factor <= required; ++factor) {
                if (factor == basis || factor == differentiated) continue;
                term *= (0.0 - nodes[static_cast<std::size_t>(factor)])
                    / (nodes[static_cast<std::size_t>(basis)]
                       - nodes[static_cast<std::size_t>(factor)]);
            }
            basis_derivative += term;
        }
        derivative += basis_derivative * values[static_cast<std::size_t>(basis)];
    }
    return derivative;
}

Real interpolate_internal_pressure_trace(const StructuredBlock& block,
                                         const AlgorithmProfile& profile,
                                         Axis axis,
                                         Index3 face)
{
    const int count = block.cell_extent()[static_cast<std::size_t>(axis)];
    const int normal = face[static_cast<std::size_t>(axis)];
    const auto& operators = cached_line_operators(profile, count);
    const auto& row = operators.interpolation_rows()[static_cast<std::size_t>(normal)];
    Real result = 0.0;
    for (const auto [center_index, coefficient] : row) {
        auto center = face;
        center[static_cast<std::size_t>(axis)] = center_index;
        result += coefficient * block.flow.primitive(center.i, center.j, center.k, pressure);
    }
    if (!std::isfinite(result)) {
        throw PhysicsError("internal wall pressure trace is non-finite");
    }
    return result;
}

ViscousFaceTrace apply_viscous_boundary_trace(const StructuredBlock& block,
                                              const MetricField& metric,
                                              const BoundaryPatch& patch,
                                              Index3 face,
                                              const ViscousFaceTrace& raw_trace,
                                              Real internal_pressure_trace,
                                              Normal3 outward_unit_normal,
                                              const BoundaryData& data,
                                              const AlgorithmProfile& profile,
                                              const GasModel& gas,
                                              const ReferenceScales& reference,
                                              const NumericalFloors& floors)
{
    require_unit_normal(outward_unit_normal);
    data.validate(patch.type, block.cell_dimension());
    (void)pressure_primitive(raw_trace.state, gas, reference, floors, block.cell_dimension());
    if (patch.type != BoundaryType::NoSlipAdiabaticWall
        && patch.type != BoundaryType::NoSlipIsothermalWall) {
        return raw_trace;
    }
    const auto normal = as_array(outward_unit_normal);
    if (std::abs(dot(data.wall_velocity, normal)) > 1.0e-12) {
        throw PhysicsConfigurationError(
            "static no-slip wall velocity must have zero normal component");
    }
    const int stencil_count = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 4 : 6;
    if (block.cell_extent()[static_cast<std::size_t>(patch.face.axis)] < stencil_count) {
        throw ProfileError("wall Dirichlet stencil does not fit the block");
    }
    const Real scale = normal_physical_scale(block, metric, patch, face, outward_unit_normal, profile);
    ViscousFaceTrace result = raw_trace;
    for (int velocity = 0; velocity < 3; ++velocity) {
        result.state[static_cast<std::size_t>(temperature_velocity_x + velocity)]
            = data.wall_velocity[static_cast<std::size_t>(velocity)];
        const auto values
            = inward_values(block, patch, face, temperature_velocity_x + velocity, stencil_count);
        const Real derivative
            = scale
            * wall_dirichlet_computational_derivative(
                  data.wall_velocity[static_cast<std::size_t>(velocity)], values, profile);
        result.gradients[static_cast<std::size_t>(velocity)]
            = normal_gradient(derivative, outward_unit_normal);
    }
    if (patch.type == BoundaryType::NoSlipAdiabaticWall) {
        result.gradients[static_cast<int>(ViscousPrimitive::Temperature)]
            = remove_normal(raw_trace.gradients[static_cast<int>(ViscousPrimitive::Temperature)],
                            outward_unit_normal);
    } else {
        result.state[temperature_value] = *data.wall_temperature;
        const auto values = inward_values(block, patch, face, temperature_value, stencil_count);
        const Real derivative = scale
            * wall_dirichlet_computational_derivative(*data.wall_temperature, values, profile);
        result.gradients[static_cast<int>(ViscousPrimitive::Temperature)]
            = normal_gradient(derivative, outward_unit_normal);
    }
    const Real temperature = result.state[temperature_value];
    result.state[temperature_density]
        = gas.gamma() * reference.mach() * reference.mach() * internal_pressure_trace / temperature;
    (void)pressure_primitive(result.state, gas, reference, floors, block.cell_dimension());
    if (block.cell_dimension() == 2) {
        result.state[temperature_velocity_z] = 0.0;
        result.gradients[static_cast<int>(ViscousPrimitive::VelocityZ)] = {{0.0, 0.0, 0.0}};
        for (auto& gradient : result.gradients)
            gradient[2] = 0.0;
    }
    return result;
}

} // namespace wcns
