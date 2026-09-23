#include <wcns/solver/sa_negative_transport.hpp>

#include <wcns/mesh/linear_operators.hpp>
#include <wcns/solver/wall_distance.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace wcns {
namespace {

bool contains(const IndexRange3& range, Index3 index)
{
    for (int axis = 0; axis < 3; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        const int lower = std::min(range.begin[a], range.end[a]);
        const int upper = std::max(range.begin[a], range.end[a]);
        if (index[a] < lower || index[a] > upper) return false;
    }
    return true;
}

bool connection_covers(const StructuredBlock& block, Axis axis, Side side, Index3 face)
{
    for (const auto& connection : block.connectivities) {
        if (connection.receiver_face.axis == axis && connection.receiver_face.side == side
            && contains(connection.shared_face_range.untyped(), face)) {
            return true;
        }
    }
    return false;
}

const BoundaryPatch* physical_patch(const StructuredBlock& block, Axis axis, Index3 face)
{
    for (const auto& patch : block.boundaries) {
        if (patch.face.axis == axis && contains(patch.boundary_face_range.untyped(), face)) {
            return &patch;
        }
    }
    return nullptr;
}

const FaceAreaVectors& faces(const MetricField& metric, Axis axis)
{
    if (axis == Axis::I) return metric.i_faces();
    if (axis == Axis::J) return metric.j_faces();
    return metric.k_faces();
}

std::array<Real, 3> area_vector(const MetricField& metric, Axis axis, Index3 face)
{
    const auto& values = faces(metric, axis);
    return {{values.x(face.i, face.j, face.k),
             values.y(face.i, face.j, face.k),
             values.z(face.i, face.j, face.k)}};
}

Real area_magnitude(const std::array<Real, 3>& area)
{
    return std::sqrt(area[0] * area[0] + area[1] * area[1] + area[2] * area[2]);
}

Index3 adjacent_cell(Index3 face, FaceLocation location, Extent3 extent)
{
    face[static_cast<std::size_t>(location.axis)]
        = location.side == Side::Lower ? 0
                                       : extent[static_cast<std::size_t>(location.axis)] - 1;
    return face;
}

Index3 mirror_cell(Index3 face, FaceLocation location, int layer, Extent3 extent)
{
    face[static_cast<std::size_t>(location.axis)]
        = location.side == Side::Lower ? layer - 1
                                       : extent[static_cast<std::size_t>(location.axis)] - layer;
    return face;
}

Index3 ghost_cell(Index3 face, FaceLocation location, int layer, Extent3 extent)
{
    face[static_cast<std::size_t>(location.axis)]
        = location.side == Side::Lower ? -layer
                                       : extent[static_cast<std::size_t>(location.axis)] - 1 + layer;
    return face;
}

bool is_no_slip(BoundaryType type)
{
    return type == BoundaryType::NoSlipAdiabaticWall
        || type == BoundaryType::NoSlipIsothermalWall;
}

bool uses_farfield_value(BoundaryType type)
{
    return type == BoundaryType::Farfield || type == BoundaryType::Inflow;
}

Real interpolate_model_face(const StructuredBlock& block,
                            const AlgorithmProfile& profile,
                            Axis axis,
                            Index3 face)
{
    const auto& field = block.turbulence.storage();
    const int component = block.turbulence.component("nu_tilde");
    Real result = 0.0;
    if (profile.kind() == AlgorithmProfileKind::PhengleiWcns) {
        constexpr std::array<int, 4> offsets {{-2, -1, 0, 1}};
        constexpr std::array<Real, 4> coefficients {
            {-1.0 / 16.0, 9.0 / 16.0, 9.0 / 16.0, -1.0 / 16.0}};
        for (std::size_t entry = 0; entry < offsets.size(); ++entry) {
            auto cell = face;
            cell[static_cast<std::size_t>(axis)] += offsets[entry];
            result += coefficients[entry] * field(cell.i, cell.j, cell.k, component);
        }
    } else {
        constexpr std::array<int, 6> offsets {{-3, -2, -1, 0, 1, 2}};
        constexpr std::array<Real, 6> coefficients {
            {3.0 / 256.0,
             -25.0 / 256.0,
             150.0 / 256.0,
             150.0 / 256.0,
             -25.0 / 256.0,
             3.0 / 256.0}};
        for (std::size_t entry = 0; entry < offsets.size(); ++entry) {
            auto cell = face;
            cell[static_cast<std::size_t>(axis)] += offsets[entry];
            result += coefficients[entry] * field(cell.i, cell.j, cell.k, component);
        }
    }
    if (!std::isfinite(result)) throw PhysicsError("SA-neg face interpolation is non-finite");
    return result;
}

Real second_order_upwind(const Field<Real>& field,
                         int component,
                         Axis axis,
                         Index3 face,
                         bool left)
{
    auto first = face;
    auto second = face;
    if (left) {
        first[static_cast<std::size_t>(axis)] -= 1;
        second[static_cast<std::size_t>(axis)] -= 2;
    } else {
        second[static_cast<std::size_t>(axis)] += 1;
    }
    const Real result = 1.5 * field(first.i, first.j, first.k, component)
        - 0.5 * field(second.i, second.j, second.k, component);
    if (!std::isfinite(result)) throw PhysicsError("SA-neg upwind trace is non-finite");
    return result;
}

TemperaturePrimitiveState load_mean_state(const StructuredBlock& block, Index3 cell)
{
    TemperaturePrimitiveState result {};
    for (int component = 0; component < fluid_components; ++component) {
        result[static_cast<std::size_t>(component)]
            = block.flow.temperature_primitive(cell.i, cell.j, cell.k, component);
    }
    return result;
}

TurbulenceCellContext make_cell_context(const StructuredBlock& block,
                                        Index3 cell,
                                        const PrimitiveGradientField& mean_gradients,
                                        const PrimitiveGradientField& model_gradients,
                                        const TransportModel& transport,
                                        const GasModel& gas,
                                        const ReferenceScales& reference)
{
    TurbulenceCellContext context;
    context.mean_state = load_mean_state(block, cell);
    for (int variable = 0; variable < viscous_primitive_components; ++variable) {
        for (int direction = 0; direction < 3; ++direction) {
            context.primitive_gradients[static_cast<std::size_t>(variable)]
                                       [static_cast<std::size_t>(direction)]
                = mean_gradients(cell,
                                 static_cast<ViscousPrimitive>(variable),
                                 direction);
        }
    }
    context.model_values = {block.turbulence.at(cell, "nu_tilde")};
    context.model_gradients = {{{model_gradients(cell, ViscousPrimitive::Temperature, 0),
                                 model_gradients(cell, ViscousPrimitive::Temperature, 1),
                                 model_gradients(cell, ViscousPrimitive::Temperature, 2)}}};
    const Real rho = context.mean_state[temperature_density];
    context.molecular_kinematic_viscosity
        = transport.viscosity(context.mean_state[temperature_value])
        / (rho * reference.reynolds());
    context.wall_distance = block.turbulence.at(cell, "wall_distance");
    context.reference_reynolds = reference.reynolds();
    context.reference_mach = reference.mach();
    context.heat_capacity_ratio = gas.gamma();
    context.dimension = block.cell_dimension();
    return context;
}

Real centered_face_derivative(const Field<Real>& field,
                              Axis axis,
                              Index3 cell,
                              AlgorithmProfileKind profile)
{
    const auto value = [&](int offset) {
        auto face = cell;
        face[static_cast<std::size_t>(axis)] += offset;
        const Real result = field(face.i, face.j, face.k, density);
        if (!std::isfinite(result)) {
            throw PhysicsError("SA-neg face-flux halo is incomplete");
        }
        return result;
    };
    if (profile == AlgorithmProfileKind::PhengleiWcns) {
        return (value(-1) - 27.0 * value(0) + 27.0 * value(1) - value(2)) / 24.0;
    }
    return (-9.0 * value(-2) + 125.0 * value(-1) - 2250.0 * value(0)
            + 2250.0 * value(1) - 125.0 * value(2) + 9.0 * value(3))
        / 1920.0;
}

} // namespace

void initialize_sa_negative_fields(const MpiRuntime& mpi,
                                   LocalBlockSet& local_blocks,
                                   const TurbulenceModelConfig& config,
                                   Real reference_reynolds)
{
    if (config.kind != TurbulenceModelKind::SaNegative) {
        throw std::invalid_argument("SA-neg initialization requires model=sa_neg");
    }
    const auto model = TurbulenceModelRegistry::create_builtin().create(config);
    const auto local_walls = extract_wall_primitives(local_blocks.blocks());
    const auto global_walls = collect_global_wall_primitives(mpi, local_walls);
    if (global_walls.empty()) {
        throw PhysicsConfigurationError("SA-neg resolved-wall mode requires a no-slip wall");
    }
    const WallDistanceIndex wall_index(global_walls);
    const Real farfield = sa_negative_farfield_value(config, reference_reynolds);
    for (auto& block : local_blocks.blocks()) {
        block.turbulence.reset(block.cell_extent(), block.ghost_width(), model->fields());
        block.turbulence.fill(0.0);
        wall_index.fill_cell_field(
            block, block.turbulence.storage(), block.turbulence.component("wall_distance"));
        const auto cells = block.cell_extent();
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    block.turbulence.at({i, j, k}, "nu_tilde") = farfield;
                }
            }
        }
        block.turbulence.validate_interior();
    }
}

void synchronize_sa_negative_fields(const HaloExchanger& exchanger,
                                    LocalBlockSet& local_blocks,
                                    Real farfield_nu_tilde)
{
    if (!std::isfinite(farfield_nu_tilde) || farfield_nu_tilde <= 0.0) {
        throw std::invalid_argument("SA-neg farfield value must be positive and finite");
    }
    if (local_blocks.blocks().empty()) return;
    const int components = local_blocks.blocks().front().turbulence.components();
    BlockFieldRegistry fields(components);
    for (auto& block : local_blocks.blocks()) {
        if (block.turbulence.components() != components
            || !block.turbulence.contains("nu_tilde")) {
            throw std::invalid_argument("SA-neg block fields are inconsistent");
        }
        fields.add(block.id(), block.turbulence.storage());
    }
    exchanger.exchange(fields);

    for (auto& block : local_blocks.blocks()) {
        const auto extent = block.cell_extent();
        auto& storage = block.turbulence.storage();
        const int transported = block.turbulence.component("nu_tilde");
        for (const auto& patch : block.boundaries) {
            const bool dirichlet = is_no_slip(patch.type) || uses_farfield_value(patch.type);
            const Real boundary_value = is_no_slip(patch.type) ? 0.0 : farfield_nu_tilde;
            const auto counts = patch.boundary_face_range.counts();
            for (int ok = 0; ok < counts.nk; ++ok) {
                for (int oj = 0; oj < counts.nj; ++oj) {
                    for (int oi = 0; oi < counts.ni; ++oi) {
                        const auto face = patch.boundary_face_range.at({oi, oj, ok});
                        for (int layer = 1; layer <= block.ghost_width(); ++layer) {
                            const auto mirror = mirror_cell(face, patch.face, layer, extent);
                            const auto ghost = ghost_cell(face, patch.face, layer, extent);
                            for (int component = 0; component < components; ++component) {
                                const Real interior
                                    = storage(mirror.i, mirror.j, mirror.k, component);
                                storage(ghost.i, ghost.j, ghost.k, component) = component == transported
                                        && dirichlet
                                    ? 2.0 * boundary_value - interior
                                    : interior;
                            }
                        }
                    }
                }
            }
        }
    }
}

void compute_sa_negative_gradient(PrimitiveGradientField& result,
                                  const StructuredBlock& block,
                                  const MetricField& metric,
                                  const AlgorithmProfile& profile,
                                  std::uint64_t version)
{
    if (metric.profile() != profile.kind() || metric.dimension() != block.cell_dimension()
        || result.profile() != profile.kind() || result.dimension() != block.cell_dimension()) {
        throw ProfileError("SA-neg gradient metadata are incompatible");
    }
    result.reset(version);
    const auto cells = block.cell_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const Index3 cell {i, j, k};
                const Real jacobian = metric.jacobian()(i, j, k);
                if (!std::isfinite(jacobian) || jacobian <= 0.0) {
                    throw PhysicsError("SA-neg gradient has an invalid Jacobian");
                }
                std::array<Real, 3> gradient {{0.0, 0.0, 0.0}};
                for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                    const auto axis = static_cast<Axis>(logical);
                    auto lower = cell;
                    auto upper = cell;
                    ++upper[static_cast<std::size_t>(axis)];
                    const Real lower_value = interpolate_model_face(block, profile, axis, lower);
                    const Real upper_value = interpolate_model_face(block, profile, axis, upper);
                    const auto lower_area = area_vector(metric, axis, lower);
                    const auto upper_area = area_vector(metric, axis, upper);
                    for (int direction = 0; direction < 3; ++direction) {
                        gradient[static_cast<std::size_t>(direction)]
                            += upper_value * upper_area[static_cast<std::size_t>(direction)]
                            - lower_value * lower_area[static_cast<std::size_t>(direction)];
                    }
                }
                for (int variable = 0; variable < viscous_primitive_components; ++variable) {
                    for (int direction = 0; direction < 3; ++direction) {
                        result(cell, static_cast<ViscousPrimitive>(variable), direction)
                            = variable == static_cast<int>(ViscousPrimitive::Temperature)
                            ? gradient[static_cast<std::size_t>(direction)] / jacobian
                            : 0.0;
                    }
                }
            }
        }
    }
}

void compute_sa_negative_flux_and_source(
    ViscousFaceFluxField& flux,
    Field<Real>& residual,
    Field<Real>& source_jacobian,
    StructuredBlock& block,
    const MetricField& metric,
    const InviscidFaceFluxField& mean_flux,
    const PrimitiveGradientField& mean_gradients,
    const PrimitiveGradientField& model_gradients,
    const AlgorithmProfile& profile,
    const ITurbulenceModel& model,
    const TransportModel& transport,
    const GasModel& gas,
    const ReferenceScales& reference,
    std::uint64_t version)
{
    if (model.config().kind != TurbulenceModelKind::SaNegative
        || residual.components() != 1 || source_jacobian.components() != 1
        || residual.interior_extent() != block.cell_extent()
        || source_jacobian.interior_extent() != block.cell_extent()
        || mean_flux.version() != version || mean_gradients.version() != version
        || model_gradients.version() != version) {
        throw std::invalid_argument("SA-neg transport inputs are incompatible");
    }
    if (metric.profile() != profile.kind()) {
        throw ProfileError("SA-neg transport profile differs from metric profile");
    }
    residual.fill(0.0);
    source_jacobian.fill(0.0);
    flux.reset(version);
    const auto cells = block.cell_extent();
    const int nu_component = block.turbulence.component("nu_tilde");

    const auto compute_axis = [&](Axis axis) {
        const auto extent = faces(metric, axis).x.interior_extent();
        auto& output = flux.field(axis);
        const auto& mass = mean_flux.field(axis);
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 face {i, j, k};
                    const Real mass_flux = mass(i, j, k, density);
                    const Real upwind = second_order_upwind(block.turbulence.storage(),
                                                            nu_component,
                                                            axis,
                                                            face,
                                                            mass_flux >= 0.0);
                    auto left = face;
                    --left[static_cast<std::size_t>(axis)];
                    auto right = face;
                    const auto clamp_cell = [&](Index3 cell) {
                        for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                            const auto a = static_cast<std::size_t>(logical);
                            cell[a] = std::clamp(cell[a], 0, cells[a] - 1);
                        }
                        return cell;
                    };
                    const auto left_cell = clamp_cell(left);
                    const auto right_cell = clamp_cell(right);
                    const auto left_context = make_cell_context(block,
                                                                left_cell,
                                                                mean_gradients,
                                                                model_gradients,
                                                                transport,
                                                                gas,
                                                                reference);
                    const auto right_context = make_cell_context(block,
                                                                 right_cell,
                                                                 mean_gradients,
                                                                 model_gradients,
                                                                 transport,
                                                                 gas,
                                                                 reference);
                    const Real diffusion = 0.5
                        * (model.diffusion_coefficients(left_context).front()
                           + model.diffusion_coefficients(right_context).front());
                    const auto area = area_vector(metric, axis, face);
                    Real normal_gradient = 0.0;
                    if (const auto* patch = physical_patch(block, axis, face);
                        patch != nullptr && is_no_slip(patch->type)) {
                        const auto inside = adjacent_cell(face, patch->face, cells);
                        const Real value = block.turbulence.at(inside, "nu_tilde");
                        const Real distance = block.turbulence.at(inside, "wall_distance");
                        const Real sign = patch->face.side == Side::Lower ? -1.0 : 1.0;
                        normal_gradient = -sign * value * area_magnitude(area) / distance;
                    } else {
                        for (int direction = 0; direction < 3; ++direction) {
                            const Real gradient = 0.5
                                * (model_gradients(left_cell, ViscousPrimitive::Temperature, direction)
                                   + model_gradients(right_cell,
                                                     ViscousPrimitive::Temperature,
                                                     direction));
                            normal_gradient
                                += gradient * area[static_cast<std::size_t>(direction)];
                        }
                    }
                    const Real value = mass_flux * upwind - diffusion * normal_gradient;
                    if (!std::isfinite(value)) {
                        throw PhysicsError("SA-neg computational face flux is non-finite");
                    }
                    for (int component = 0; component < euler_components; ++component) {
                        output(i, j, k, component) = component == density ? value : 0.0;
                    }
                }
            }
        }
    };
    compute_axis(Axis::I);
    compute_axis(Axis::J);
    if (block.cell_dimension() == 3) compute_axis(Axis::K);

    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const Index3 cell {i, j, k};
                const auto context = make_cell_context(block,
                                                       cell,
                                                       mean_gradients,
                                                       model_gradients,
                                                       transport,
                                                       gas,
                                                       reference);
                const auto evaluation = evaluate_sa_negative(context);
                const auto source = model.source_linearization(context);
                residual(i, j, k, 0) = source.source.front();
                source_jacobian(i, j, k, 0) = source.jacobian.front();
                const Real mu = transport.viscosity(context.mean_state[temperature_value]);
                block.turbulence.at(cell, "mu_t_over_mu")
                    = context.mean_state[temperature_density] * reference.reynolds()
                    * evaluation.eddy_kinematic_viscosity / mu;
                block.turbulence.at(cell, "sa_production")
                    = context.mean_state[temperature_density] * evaluation.production_source;
                block.turbulence.at(cell, "sa_destruction")
                    = context.mean_state[temperature_density] * evaluation.destruction_source;
                block.turbulence.at(cell, "sa_negative_branch")
                    = evaluation.negative_branch ? 1.0 : 0.0;
            }
        }
    }
}

void assemble_sa_negative_residual(Field<Real>& residual,
                                   const StructuredBlock& block,
                                   const MetricField& metric,
                                   const ViscousFaceFluxField& flux,
                                   const AlgorithmProfile& profile)
{
    if (residual.components() != 1 || residual.interior_extent() != block.cell_extent()
        || metric.profile() != profile.kind() || flux.profile() != profile.kind()) {
        throw std::invalid_argument("SA-neg residual inputs are incompatible");
    }
    const auto cells = block.cell_extent();
    const auto accumulate_axis = [&](Axis axis) {
        const int count = cells[static_cast<std::size_t>(axis)];
        const auto& rows = cached_line_operators(profile, count).derivative_rows();
        const auto& values = flux.field(axis);
        const int boundary_width = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 1 : 2;
        const bool lower_complete = connection_side_is_fully_covered(block, axis, Side::Lower);
        const bool upper_complete = connection_side_is_fully_covered(block, axis, Side::Upper);
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    const Index3 cell {i, j, k};
                    const int normal = cell[static_cast<std::size_t>(axis)];
                    Index3 lower = cell;
                    lower[static_cast<std::size_t>(axis)] = 0;
                    Index3 upper = cell;
                    upper[static_cast<std::size_t>(axis)] = count;
                    Real derivative = 0.0;
                    if ((lower_complete && normal < boundary_width
                         && connection_covers(block, axis, Side::Lower, lower))
                        || (upper_complete && normal >= count - boundary_width
                            && connection_covers(block, axis, Side::Upper, upper))) {
                        derivative = centered_face_derivative(values, axis, cell, profile.kind());
                    } else {
                        for (const auto [face_index, coefficient] : rows[static_cast<std::size_t>(normal)]) {
                            auto face = cell;
                            face[static_cast<std::size_t>(axis)] = face_index;
                            derivative += coefficient
                                * values(face.i, face.j, face.k, density);
                        }
                    }
                    const Real jacobian = metric.jacobian()(i, j, k);
                    if (!std::isfinite(derivative) || !std::isfinite(jacobian)
                        || jacobian <= 0.0) {
                        throw PhysicsError("SA-neg flux divergence is invalid");
                    }
                    residual(i, j, k, 0) -= derivative / jacobian;
                }
            }
        }
    };
    accumulate_axis(Axis::I);
    accumulate_axis(Axis::J);
    if (block.cell_dimension() == 3) accumulate_axis(Axis::K);
}

} // namespace wcns
