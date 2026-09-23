#include <wcns/solver/rans_two_equation_transport.hpp>

#include <wcns/solver/two_equation_models.hpp>
#include <wcns/solver/wall_distance.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace wcns {
namespace {

bool is_wall(BoundaryType type)
{
    return type == BoundaryType::NoSlipAdiabaticWall
        || type == BoundaryType::NoSlipIsothermalWall;
}

bool uses_farfield_value(BoundaryType type)
{
    return type == BoundaryType::Farfield || type == BoundaryType::Inflow;
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

const FaceAreaVectors& face_vectors(const MetricField& metric, Axis axis)
{
    if (axis == Axis::I) return metric.i_faces();
    if (axis == Axis::J) return metric.j_faces();
    return metric.k_faces();
}

std::array<Real, 3> area_vector(const MetricField& metric, Axis axis, Index3 face)
{
    const auto& values = face_vectors(metric, axis);
    return {{values.x(face.i, face.j, face.k),
             values.y(face.i, face.j, face.k),
             values.z(face.i, face.j, face.k)}};
}

std::vector<TurbulenceFieldDescriptor> transported_fields(const ITurbulenceModel& model)
{
    std::vector<TurbulenceFieldDescriptor> result;
    for (const auto& descriptor : model.fields()) {
        if (descriptor.role == TurbulenceFieldRole::Transported) result.push_back(descriptor);
    }
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
        --first[static_cast<std::size_t>(axis)];
        second = first;
        --second[static_cast<std::size_t>(axis)];
    } else {
        second[static_cast<std::size_t>(axis)] += 1;
    }
    return 1.5 * field(first.i, first.j, first.k, component)
        - 0.5 * field(second.i, second.j, second.k, component);
}

Index3 clamp_cell(Index3 cell, Extent3 extent)
{
    for (int axis = 0; axis < 3; ++axis) {
        const auto entry = static_cast<std::size_t>(axis);
        cell[entry] = std::clamp(cell[entry], 0, extent[entry] - 1);
    }
    return cell;
}

TurbulenceCellContext make_context(const StructuredBlock& block,
                                   Index3 cell,
                                   const PrimitiveGradientField& mean_gradients,
                                   const Field<Real>& model_gradients,
                                   const std::vector<TurbulenceFieldDescriptor>& descriptors,
                                   const TransportModel& transport,
                                   const GasModel& gas,
                                   const ReferenceScales& reference)
{
    TurbulenceCellContext context;
    for (int component = 0; component < fluid_components; ++component) {
        context.mean_state[static_cast<std::size_t>(component)]
            = block.flow.temperature_primitive(cell.i, cell.j, cell.k, component);
    }
    for (int variable = 0; variable < viscous_primitive_components; ++variable) {
        for (int direction = 0; direction < 3; ++direction) {
            context.primitive_gradients[static_cast<std::size_t>(variable)]
                                       [static_cast<std::size_t>(direction)]
                = mean_gradients(cell, static_cast<ViscousPrimitive>(variable), direction);
        }
    }
    for (std::size_t variable = 0; variable < descriptors.size(); ++variable) {
        context.model_values.push_back(block.turbulence.at(cell, descriptors[variable].name));
        context.model_gradients.push_back({{
            model_gradients(cell.i, cell.j, cell.k, static_cast<int>(3 * variable)),
            model_gradients(cell.i, cell.j, cell.k, static_cast<int>(3 * variable + 1)),
            model_gradients(cell.i, cell.j, cell.k, static_cast<int>(3 * variable + 2)),
        }});
    }
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

Real equilibrium_friction_velocity(Real velocity, Real distance, Real nu)
{
    constexpr Real kappa = 0.41;
    constexpr Real e = 9.8;
    Real result = std::max(std::sqrt(std::max(nu * velocity / distance, 1.0e-30)), 1.0e-12);
    for (int iteration = 0; iteration < 30; ++iteration) {
        const Real y_plus = std::max(result * distance / nu, 1.0 + 1.0e-12);
        const Real logarithm = std::log(e * y_plus);
        const Real function = velocity / result - logarithm / kappa;
        const Real derivative = -velocity / (result * result) - 1.0 / (kappa * result);
        const Real candidate = result - function / derivative;
        if (!std::isfinite(candidate) || candidate <= 0.0) break;
        if (std::abs(candidate - result) <= 1.0e-12 * std::max(Real {1.0}, result)) {
            result = candidate;
            break;
        }
        result = candidate;
    }
    return result;
}

} // namespace

void initialize_two_equation_fields(const MpiRuntime& mpi,
                                    LocalBlockSet& local_blocks,
                                    const TurbulenceModelConfig& config)
{
    if (config.kind != TurbulenceModelKind::KOmegaSst
        && config.kind != TurbulenceModelKind::KEpsilon) {
        throw std::invalid_argument("two-equation initialization requires SST or k-epsilon");
    }
    const auto model = TurbulenceModelRegistry::create_builtin().create(config);
    const auto walls = collect_global_wall_primitives(mpi, extract_wall_primitives(local_blocks.blocks()));
    if (walls.empty()) {
        throw PhysicsConfigurationError("two-equation RANS requires a no-slip wall");
    }
    const WallDistanceIndex wall_index(walls);
    const auto farfield = two_equation_farfield_values(config);
    const auto transported = transported_fields(*model);
    for (auto& block : local_blocks.blocks()) {
        block.turbulence.reset(block.cell_extent(), block.ghost_width(), model->fields());
        block.turbulence.fill(0.0);
        wall_index.fill_cell_field(
            block, block.turbulence.storage(), block.turbulence.component("wall_distance"));
        const auto extent = block.cell_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 cell {i, j, k};
                    block.turbulence.at(cell, transported[0].name) = farfield[0];
                    block.turbulence.at(cell, transported[1].name) = farfield[1];
                }
            }
        }
        block.turbulence.validate_interior();
    }
}

void synchronize_two_equation_fields(const HaloExchanger& exchanger,
                                     LocalBlockSet& local_blocks,
                                     const TurbulenceModelConfig& config,
                                     const TransportModel& transport,
                                     const ReferenceScales& reference)
{
    const auto model = TurbulenceModelRegistry::create_builtin().create(config);
    const auto descriptors = transported_fields(*model);
    const auto farfield = two_equation_farfield_values(config);
    if (local_blocks.blocks().empty()) return;
    const int components = local_blocks.blocks().front().turbulence.components();
    BlockFieldRegistry fields(components);
    for (auto& block : local_blocks.blocks()) fields.add(block.id(), block.turbulence.storage());
    exchanger.exchange(fields);

    for (auto& block : local_blocks.blocks()) {
        auto& storage = block.turbulence.storage();
        const auto extent = block.cell_extent();
        for (const auto& patch : block.boundaries) {
            const auto counts = patch.boundary_face_range.counts();
            for (int ok = 0; ok < counts.nk; ++ok) {
                for (int oj = 0; oj < counts.nj; ++oj) {
                    for (int oi = 0; oi < counts.ni; ++oi) {
                        const auto face = patch.boundary_face_range.at({oi, oj, ok});
                        const auto first = mirror_cell(face, patch.face, 1, extent);
                        std::array<Real, 2> boundary {{
                            block.turbulence.at(first, descriptors[0].name),
                            block.turbulence.at(first, descriptors[1].name),
                        }};
                        bool dirichlet = false;
                        if (uses_farfield_value(patch.type)) {
                            boundary = farfield;
                            dirichlet = true;
                        } else if (is_wall(patch.type)) {
                            dirichlet = true;
                            const Real distance = block.turbulence.at(first, "wall_distance");
                            const Real rho = block.flow.temperature_primitive(
                                first.i, first.j, first.k, temperature_density);
                            const Real temperature = block.flow.temperature_primitive(
                                first.i, first.j, first.k, temperature_value);
                            const Real nu = transport.viscosity(temperature)
                                / (rho * reference.reynolds());
                            if (config.wall_treatment == WallTreatment::Resolved) {
                                const SstConstants constants;
                                boundary = {{0.0,
                                             60.0 * nu
                                                 / (constants.beta1 * distance * distance)}};
                            } else {
                                const Real u = block.flow.temperature_primitive(
                                    first.i, first.j, first.k, temperature_velocity_x);
                                const Real v = block.flow.temperature_primitive(
                                    first.i, first.j, first.k, temperature_velocity_y);
                                const Real w = block.flow.temperature_primitive(
                                    first.i, first.j, first.k, temperature_velocity_z);
                                const Real friction = equilibrium_friction_velocity(
                                    std::sqrt(u * u + v * v + w * w), distance, nu);
                                const Real y_plus = friction * distance / nu;
                                block.turbulence.at(first, "wall_y_plus") = y_plus;
                                const StandardKEpsilonConstants constants;
                                const Real k_wall
                                    = friction * friction / std::sqrt(constants.c_mu);
                                if (config.kind == TurbulenceModelKind::KEpsilon) {
                                    boundary = {{k_wall,
                                                 friction * friction * friction
                                                     / (constants.kappa * distance)}};
                                } else {
                                    boundary = {{k_wall,
                                                 friction
                                                     / (constants.kappa * distance
                                                        * std::pow(constants.c_mu, 0.25))}};
                                }
                            }
                        }
                        for (int layer = 1; layer <= block.ghost_width(); ++layer) {
                            const auto mirror = mirror_cell(face, patch.face, layer, extent);
                            const auto ghost = ghost_cell(face, patch.face, layer, extent);
                            for (int component = 0; component < components; ++component) {
                                storage(ghost.i, ghost.j, ghost.k, component)
                                    = storage(mirror.i, mirror.j, mirror.k, component);
                            }
                            if (dirichlet) {
                                for (int variable = 0; variable < 2; ++variable) {
                                    const int component
                                        = block.turbulence.component(descriptors[variable].name);
                                    storage(ghost.i, ghost.j, ghost.k, component)
                                        = 2.0 * boundary[static_cast<std::size_t>(variable)]
                                        - storage(mirror.i, mirror.j, mirror.k, component);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

WallFunctionDiagnostics global_wall_function_diagnostics(
    const MpiRuntime& mpi,
    const LocalBlockSet& local_blocks,
    const TurbulenceModelConfig& config)
{
    WallFunctionDiagnostics local;
    if (config.wall_treatment == WallTreatment::WallFunction) {
        for (const auto& block : local_blocks.blocks()) {
            const auto extent = block.cell_extent();
            for (const auto& patch : block.boundaries) {
                if (!is_wall(patch.type)) continue;
                const auto counts = patch.boundary_face_range.counts();
                for (int ok = 0; ok < counts.nk; ++ok) {
                    for (int oj = 0; oj < counts.nj; ++oj) {
                        for (int oi = 0; oi < counts.ni; ++oi) {
                            const auto face = patch.boundary_face_range.at({oi, oj, ok});
                            const auto cell = mirror_cell(face, patch.face, 1, extent);
                            const Real y_plus = block.turbulence.at(cell, "wall_y_plus");
                            if (!std::isfinite(y_plus) || y_plus <= 0.0) {
                                throw PhysicsError("wall-function y+ diagnostic is invalid");
                            }
                            ++local.face_count;
                            if (y_plus < config.wall_function_y_plus_min
                                || y_plus > config.wall_function_y_plus_max) {
                                ++local.out_of_range_count;
                            }
                            local.minimum_y_plus = std::min(local.minimum_y_plus, y_plus);
                            local.maximum_y_plus = std::max(local.maximum_y_plus, y_plus);
                        }
                    }
                }
            }
        }
    }
    WallFunctionDiagnostics result;
    result.face_count = static_cast<std::size_t>(mpi.sum(static_cast<Real>(local.face_count)));
    result.out_of_range_count
        = static_cast<std::size_t>(mpi.sum(static_cast<Real>(local.out_of_range_count)));
    result.minimum_y_plus = mpi.min(local.minimum_y_plus);
    result.maximum_y_plus = mpi.max(local.maximum_y_plus);
    return result;
}

void compute_two_equation_gradients_and_source(
    Field<Real>& residual,
    Field<Real>& source_jacobian,
    Field<Real>& face_workspace,
    StructuredBlock& block,
    const MetricField& metric,
    const PrimitiveGradientField& mean_gradients,
    const ITurbulenceModel& model,
    const TransportModel& transport,
    const GasModel& gas,
    const ReferenceScales& reference)
{
    const auto descriptors = transported_fields(model);
    if (descriptors.size() != 2 || residual.components() != 2
        || source_jacobian.components() != 4
        || face_workspace.components() != two_equation_face_workspace_components
        || face_workspace.ghost_width() < 1
        || residual.interior_extent() != block.cell_extent()
        || source_jacobian.interior_extent() != block.cell_extent()
        || face_workspace.interior_extent() != block.cell_extent()) {
        throw std::invalid_argument("two-equation transport workspace is incompatible");
    }
    const auto extent = block.cell_extent();
    residual.fill(0.0);
    source_jacobian.fill(0.0);
    face_workspace.fill(0.0);

    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const Index3 cell {i, j, k};
                const Real volume = metric.jacobian()(i, j, k);
                for (std::size_t variable = 0; variable < descriptors.size(); ++variable) {
                    std::array<Real, 3> gradient {{0.0, 0.0, 0.0}};
                    for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                        const auto axis = static_cast<Axis>(logical);
                        auto lower = cell;
                        auto upper = cell;
                        ++upper[static_cast<std::size_t>(axis)];
                        auto left = lower;
                        --left[static_cast<std::size_t>(axis)];
                        auto right = upper;
                        const int component
                            = block.turbulence.component(descriptors[variable].name);
                        const Real lower_value = 0.5
                            * (block.turbulence.storage()(left.i, left.j, left.k, component)
                               + block.turbulence.storage()(cell.i, cell.j, cell.k, component));
                        const Real upper_value = 0.5
                            * (block.turbulence.storage()(cell.i, cell.j, cell.k, component)
                               + block.turbulence.storage()(right.i, right.j, right.k, component));
                        const auto lower_area = area_vector(metric, axis, lower);
                        const auto upper_area = area_vector(metric, axis, upper);
                        for (int direction = 0; direction < 3; ++direction) {
                            gradient[static_cast<std::size_t>(direction)]
                                += upper_value * upper_area[static_cast<std::size_t>(direction)]
                                - lower_value * lower_area[static_cast<std::size_t>(direction)];
                        }
                    }
                    for (int direction = 0; direction < 3; ++direction) {
                        face_workspace(i, j, k, static_cast<int>(3 * variable + direction))
                            = gradient[static_cast<std::size_t>(direction)] / volume;
                    }
                }
            }
        }
    }

    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const Index3 cell {i, j, k};
                const auto context = make_context(block,
                                                  cell,
                                                  mean_gradients,
                                                  face_workspace,
                                                  descriptors,
                                                  transport,
                                                  gas,
                                                  reference);
                const auto source = model.source_linearization(context);
                for (int variable = 0; variable < 2; ++variable) {
                    residual(i, j, k, variable) = source.source[static_cast<std::size_t>(variable)];
                }
                for (int entry = 0; entry < 4; ++entry) {
                    source_jacobian(i, j, k, entry)
                        = source.jacobian[static_cast<std::size_t>(entry)];
                }
                const Real mu = transport.viscosity(context.mean_state[temperature_value]);
                TwoEquationEvaluation evaluation;
                if (model.config().kind == TurbulenceModelKind::KOmegaSst) {
                    evaluation = evaluate_k_omega_sst(context);
                    block.turbulence.at(cell, "sst_f1") = evaluation.blending_f1;
                    block.turbulence.at(cell, "sst_f2") = evaluation.blending_f2;
                    block.turbulence.at(cell, "cross_diffusion")
                        = evaluation.cross_diffusion;
                } else {
                    evaluation = evaluate_standard_k_epsilon(context);
                }
                block.turbulence.at(cell, "mu_t_over_mu")
                    = context.mean_state[temperature_density] * reference.reynolds()
                    * evaluation.eddy_kinematic_viscosity / mu;
                block.turbulence.at(cell, "turbulence_production")
                    = context.mean_state[temperature_density] * evaluation.production;
                block.turbulence.at(cell, "turbulence_destruction")
                    = context.mean_state[temperature_density] * evaluation.destruction_k;
                face_workspace(i, j, k, 6)
                    = evaluation.diffusion_coefficients[0];
                face_workspace(i, j, k, 7)
                    = evaluation.diffusion_coefficients[1];
            }
        }
    }
}

void fill_two_equation_workspace_physical_ghosts(const StructuredBlock& block,
                                                 Field<Real>& face_workspace)
{
    if (face_workspace.components() != two_equation_face_workspace_components
        || face_workspace.ghost_width() < 1
        || face_workspace.interior_extent() != block.cell_extent()) {
        throw std::invalid_argument("two-equation face workspace is incompatible");
    }
    const auto extent = block.cell_extent();
    for (const auto& patch : block.boundaries) {
        const auto counts = patch.boundary_face_range.counts();
        for (int ok = 0; ok < counts.nk; ++ok) {
            for (int oj = 0; oj < counts.nj; ++oj) {
                for (int oi = 0; oi < counts.ni; ++oi) {
                    const auto face = patch.boundary_face_range.at({oi, oj, ok});
                    for (int layer = 1; layer <= face_workspace.ghost_width(); ++layer) {
                        const auto mirror = mirror_cell(face, patch.face, layer, extent);
                        const auto ghost = ghost_cell(face, patch.face, layer, extent);
                        for (int component = 0;
                             component < two_equation_face_workspace_components;
                             ++component) {
                            face_workspace(ghost.i, ghost.j, ghost.k, component)
                                = face_workspace(mirror.i, mirror.j, mirror.k, component);
                        }
                    }
                }
            }
        }
    }
}

void assemble_two_equation_flux_residual(Field<Real>& residual,
                                         const Field<Real>& face_workspace,
                                         const StructuredBlock& block,
                                         const MetricField& metric,
                                         const InviscidFaceFluxField& mean_flux,
                                         const ITurbulenceModel& model)
{
    const auto descriptors = transported_fields(model);
    if (descriptors.size() != 2 || residual.components() != 2
        || face_workspace.components() != two_equation_face_workspace_components
        || face_workspace.ghost_width() < 1
        || residual.interior_extent() != block.cell_extent()
        || face_workspace.interior_extent() != block.cell_extent()) {
        throw std::invalid_argument("two-equation flux workspace is incompatible");
    }
    const auto extent = block.cell_extent();

    const auto face_flux = [&](Axis axis, Index3 face, int variable) {
        const int component = block.turbulence.component(
            descriptors[static_cast<std::size_t>(variable)].name);
        const Real mass_flux = mean_flux.field(axis)(face.i, face.j, face.k, density);
        Real advected = second_order_upwind(
            block.turbulence.storage(), component, axis, face, mass_flux >= 0.0);
        auto left = face;
        --left[static_cast<std::size_t>(axis)];
        const auto left_cell = clamp_cell(left, extent);
        const auto right_cell = clamp_cell(face, extent);
        const auto donor = mass_flux >= 0.0 ? left_cell : right_cell;
        const auto& descriptor = descriptors[static_cast<std::size_t>(variable)];
        if (!std::isfinite(advected)
            || (descriptor.strictly_positive && !(advected > descriptor.lower_bound))) {
            // A high-order trace may overshoot even when both cell averages are admissible.
            // Fall back locally to the upwind cell; the transported state itself is not clipped.
            advected = block.turbulence.storage()(donor.i, donor.j, donor.k, component);
        }
        const Real diffusion = 0.5
            * (face_workspace(left.i, left.j, left.k, 6 + variable)
               + face_workspace(face.i, face.j, face.k, 6 + variable));
        const auto area = area_vector(metric, axis, face);
        Real gradient_dot_area = 0.0;
        for (int direction = 0; direction < 3; ++direction) {
            gradient_dot_area += 0.5
                * (face_workspace(left.i,
                                   left.j,
                                   left.k,
                                   3 * variable + direction)
                   + face_workspace(face.i,
                                     face.j,
                                     face.k,
                                     3 * variable + direction))
                * area[static_cast<std::size_t>(direction)];
        }
        return mass_flux * advected - diffusion * gradient_dot_area;
    };

    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const Index3 cell {i, j, k};
                const Real volume = metric.jacobian()(i, j, k);
                for (int variable = 0; variable < 2; ++variable) {
                    Real divergence = 0.0;
                    for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                        const auto axis = static_cast<Axis>(logical);
                        auto lower = cell;
                        auto upper = cell;
                        ++upper[static_cast<std::size_t>(axis)];
                        divergence += face_flux(axis, upper, variable)
                            - face_flux(axis, lower, variable);
                    }
                    residual(i, j, k, variable) -= divergence / volume;
                    if (!std::isfinite(residual(i, j, k, variable))) {
                        throw PhysicsError("two-equation turbulence residual is non-finite");
                    }
                }
            }
        }
    }
}

} // namespace wcns
