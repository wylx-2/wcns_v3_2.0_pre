#include <wcns/solver/implicit_time_integrator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace wcns {
namespace {

int lower_component(int axis)
{
    return 2 * axis;
}

int upper_component(int axis)
{
    return 2 * axis + 1;
}

void validate_system(const Field<Real>& right_hand_side,
                     const Field<Real>& diagonal,
                     const Field<Real>& coupling,
                     int dimension)
{
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("LU-SGS dimension must be two or three");
    }
    const auto extent = right_hand_side.interior_extent();
    if (extent != diagonal.interior_extent() || extent != coupling.interior_extent()
        || right_hand_side.components() <= 0
        || (diagonal.components() != 1
            && diagonal.components() != right_hand_side.components())
        || coupling.components() != 2 * dimension || right_hand_side.ghost_width() != 0
        || diagonal.ghost_width() != 0 || coupling.ghost_width() != 0) {
        throw std::invalid_argument("LU-SGS field metadata differ");
    }
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                for (int variable = 0; variable < right_hand_side.components(); ++variable) {
                    const Real rhs = right_hand_side(i, j, k, variable);
                    const Real d = diagonal(i,
                                            j,
                                            k,
                                            diagonal.components() == 1 ? 0 : variable);
                    if (!std::isfinite(rhs) || !std::isfinite(d) || d <= 0.0) {
                        throw std::invalid_argument("LU-SGS right-hand side or diagonal is invalid");
                    }
                }
                for (int component = 0; component < coupling.components(); ++component) {
                    const Real value = coupling(i, j, k, component);
                    if (!std::isfinite(value) || value < 0.0) {
                        throw std::invalid_argument("LU-SGS coupling must be finite and non-negative");
                    }
                }
            }
        }
    }
}

void solve_dense_block(const Field<Real>& blocks,
                       Index3 cell,
                       const std::vector<Real>& right_hand_side,
                       std::vector<Real>& matrix,
                       std::vector<Real>& result)
{
    const int variables = static_cast<int>(right_hand_side.size());
    if (blocks.components() != variables * variables || variables == 0) {
        throw std::invalid_argument("LU-SGS dense diagonal block has an invalid size");
    }
    matrix.resize(static_cast<std::size_t>(variables * variables));
    result = right_hand_side;
    for (int row = 0; row < variables; ++row) {
        for (int column = 0; column < variables; ++column) {
            matrix[static_cast<std::size_t>(row * variables + column)]
                = blocks(cell.i, cell.j, cell.k, row * variables + column);
        }
    }
    for (int pivot = 0; pivot < variables; ++pivot) {
        int best = pivot;
        for (int row = pivot + 1; row < variables; ++row) {
            if (std::abs(matrix[static_cast<std::size_t>(row * variables + pivot)])
                > std::abs(matrix[static_cast<std::size_t>(best * variables + pivot)])) {
                best = row;
            }
        }
        const Real pivot_value
            = matrix[static_cast<std::size_t>(best * variables + pivot)];
        if (!std::isfinite(pivot_value)
            || std::abs(pivot_value) <= 64.0 * std::numeric_limits<Real>::epsilon()) {
            throw std::runtime_error("LU-SGS dense diagonal block is singular");
        }
        if (best != pivot) {
            for (int column = pivot; column < variables; ++column) {
                std::swap(matrix[static_cast<std::size_t>(pivot * variables + column)],
                          matrix[static_cast<std::size_t>(best * variables + column)]);
            }
            std::swap(result[static_cast<std::size_t>(pivot)],
                      result[static_cast<std::size_t>(best)]);
        }
        for (int row = pivot + 1; row < variables; ++row) {
            const Real factor
                = matrix[static_cast<std::size_t>(row * variables + pivot)]
                / matrix[static_cast<std::size_t>(pivot * variables + pivot)];
            matrix[static_cast<std::size_t>(row * variables + pivot)] = 0.0;
            for (int column = pivot + 1; column < variables; ++column) {
                matrix[static_cast<std::size_t>(row * variables + column)] -= factor
                    * matrix[static_cast<std::size_t>(pivot * variables + column)];
            }
            result[static_cast<std::size_t>(row)]
                -= factor * result[static_cast<std::size_t>(pivot)];
        }
    }
    for (int row = variables - 1; row >= 0; --row) {
        for (int column = row + 1; column < variables; ++column) {
            result[static_cast<std::size_t>(row)]
                -= matrix[static_cast<std::size_t>(row * variables + column)]
                * result[static_cast<std::size_t>(column)];
        }
        result[static_cast<std::size_t>(row)]
            /= matrix[static_cast<std::size_t>(row * variables + row)];
        if (!std::isfinite(result[static_cast<std::size_t>(row)])) {
            throw std::runtime_error("LU-SGS dense block solve is non-finite");
        }
    }
}

void multiply_dense_block(const Field<Real>& blocks,
                          Index3 cell,
                          const Field<Real>& values,
                          std::vector<Real>& result)
{
    const int variables = values.components();
    result.assign(static_cast<std::size_t>(variables), 0.0);
    for (int row = 0; row < variables; ++row) {
        for (int column = 0; column < variables; ++column) {
            result[static_cast<std::size_t>(row)]
                += blocks(cell.i, cell.j, cell.k, row * variables + column)
                * values(cell.i, cell.j, cell.k, column);
        }
    }
}

Matrix5 euler_normal_flux_jacobian(const PressurePrimitiveState& state,
                                   Normal3 normal,
                                   Real gamma)
{
    const Real rho = state[0];
    const std::array<Real, 3> velocity {{state[1], state[2], state[3]}};
    const std::array<Real, 3> direction {{normal.x, normal.y, normal.z}};
    const Real normal_velocity = velocity[0] * direction[0]
        + velocity[1] * direction[1] + velocity[2] * direction[2];
    const Real speed_squared = velocity[0] * velocity[0] + velocity[1] * velocity[1]
        + velocity[2] * velocity[2];
    const Real energy = state[4] / (gamma - 1.0) + 0.5 * rho * speed_squared;
    const Real enthalpy = (energy + state[4]) / rho;
    Matrix5 result {};
    result[0] = {{0.0, direction[0], direction[1], direction[2], 0.0}};
    for (int component = 0; component < 3; ++component) {
        const int row = component + 1;
        result[static_cast<std::size_t>(row)][0]
            = -velocity[static_cast<std::size_t>(component)] * normal_velocity
            + 0.5 * (gamma - 1.0) * speed_squared
                * direction[static_cast<std::size_t>(component)];
        for (int momentum = 0; momentum < 3; ++momentum) {
            result[static_cast<std::size_t>(row)]
                  [static_cast<std::size_t>(momentum + 1)]
                = (component == momentum ? normal_velocity : 0.0)
                + velocity[static_cast<std::size_t>(component)]
                    * direction[static_cast<std::size_t>(momentum)]
                - (gamma - 1.0) * velocity[static_cast<std::size_t>(momentum)]
                    * direction[static_cast<std::size_t>(component)];
        }
        result[static_cast<std::size_t>(row)][4]
            = (gamma - 1.0) * direction[static_cast<std::size_t>(component)];
    }
    result[4][0]
        = normal_velocity * (0.5 * (gamma - 1.0) * speed_squared - enthalpy);
    for (int momentum = 0; momentum < 3; ++momentum) {
        result[4][static_cast<std::size_t>(momentum + 1)]
            = enthalpy * direction[static_cast<std::size_t>(momentum)]
            - (gamma - 1.0) * velocity[static_cast<std::size_t>(momentum)]
                * normal_velocity;
    }
    result[4][4] = gamma * normal_velocity;
    return result;
}

void multiply_face_block(const Field<Real>& face_blocks,
                         Index3 cell,
                         int face_component,
                         const Field<Real>& values,
                         Index3 neighbor,
                         std::vector<Real>& result)
{
    const int variables = values.components();
    result.assign(static_cast<std::size_t>(variables), 0.0);
    const int offset = face_component * variables * variables;
    for (int row = 0; row < variables; ++row) {
        for (int column = 0; column < variables; ++column) {
            result[static_cast<std::size_t>(row)]
                += face_blocks(cell.i,
                               cell.j,
                               cell.k,
                               offset + row * variables + column)
                * values(neighbor.i, neighbor.j, neighbor.k, column);
        }
    }
}

Real neighbor_value(const Field<Real>& field,
                    int i,
                    int j,
                    int k,
                    int variable,
                    int axis,
                    int offset)
{
    if (axis == 0) i += offset;
    if (axis == 1) j += offset;
    if (axis == 2) k += offset;
    return field(i, j, k, variable);
}

bool has_neighbor(Extent3 extent, int i, int j, int k, int axis, int offset)
{
    if (axis == 0) i += offset;
    if (axis == 1) j += offset;
    if (axis == 2) k += offset;
    return i >= 0 && i < extent.ni && j >= 0 && j < extent.nj && k >= 0
        && k < extent.nk;
}

const FaceAreaVectors& faces(const MetricField& metric, Axis axis)
{
    if (axis == Axis::I) return metric.i_faces();
    if (axis == Axis::J) return metric.j_faces();
    return metric.k_faces();
}

std::array<Real, 3> area_vector(const MetricField& metric, Axis axis, Index3 face)
{
    const auto& field = faces(metric, axis);
    return {{field.x(face.i, face.j, face.k),
             field.y(face.i, face.j, face.k),
             field.z(face.i, face.j, face.k)}};
}

} // namespace

void LuSgsIterationConfig::validate() const
{
    if (sweeps < 1 || sweeps > 4) {
        throw std::invalid_argument("LU-SGS forward/backward sweeps must be in [1,4]");
    }
    if (!std::isfinite(relaxation) || relaxation <= 0.0 || relaxation > 1.0) {
        throw std::invalid_argument("LU-SGS relaxation must lie in (0,1]");
    }
}

void DualTimeIterationConfig::validate() const
{
    if (max_iterations == 0 || !std::isfinite(absolute_tolerance)
        || absolute_tolerance <= 0.0 || !std::isfinite(relative_tolerance)
        || relative_tolerance <= 0.0 || !std::isfinite(cfl) || cfl <= 0.0) {
        throw std::invalid_argument("invalid dual-time iteration configuration");
    }
    lu_sgs.validate();
}

ScalarSpectralSystem build_scalar_spectral_system(
    const StructuredBlock& block,
    const MetricField& metric,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    Real pseudo_cfl,
    Real physical_time_diagonal,
    const WeissSmithParameters* preconditioner,
    Real viscous_preconditioner_scale,
    const Field<Real>* additional_face_coupling)
{
    if (metric.dimension() != block.cell_dimension()
        || metric.jacobian().interior_extent() != block.cell_extent()
        || !std::isfinite(pseudo_cfl) || pseudo_cfl <= 0.0
        || !std::isfinite(physical_time_diagonal) || physical_time_diagonal < 0.0
        || !std::isfinite(viscous_preconditioner_scale)
        || viscous_preconditioner_scale < 0.0) {
        throw std::invalid_argument("invalid scalar-spectral system input");
    }
    if (preconditioner != nullptr) preconditioner->validate();
    if (additional_face_coupling != nullptr
        && (additional_face_coupling->interior_extent() != block.cell_extent()
            || additional_face_coupling->components() != 2 * block.cell_dimension()
            || additional_face_coupling->ghost_width() != 0)) {
        throw std::invalid_argument("additional implicit face coupling is incompatible");
    }
    const auto extent = block.cell_extent();
    ScalarSpectralSystem result {
        Field<Real>(extent, 1, 0, 0.0),
        Field<Real>(extent, 2 * block.cell_dimension(), 0, 0.0),
    };
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const Index3 cell {i, j, k};
                TemperaturePrimitiveState state {};
                for (int component = 0; component < fluid_components; ++component) {
                    state[static_cast<std::size_t>(component)]
                        = block.flow.temperature_primitive(i, j, k, component);
                }
                const Real u = state[temperature_velocity_x];
                const Real v = state[temperature_velocity_y];
                const Real w = state[temperature_velocity_z];
                const Real speed = std::sqrt(u * u + v * v + w * w);
                const Real sound = thermodynamic_sound_speed(
                    state, gas, reference, floors, block.cell_dimension());
                const Real volume = metric.jacobian()(i, j, k);
                Real area_sum = 0.0;
                for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                    const auto axis = static_cast<Axis>(logical);
                    for (int side = 0; side <= 1; ++side) {
                        auto face = cell;
                        face[static_cast<std::size_t>(axis)] += side;
                        const auto area = area_vector(metric, axis, face);
                        area_sum += std::sqrt(
                            area[0] * area[0] + area[1] * area[1] + area[2] * area[2]);
                    }
                }
                const Real length = 2.0 * volume / area_sum;
                const Real viscous_speed = viscous_preconditioner_scale / length;
                const Real reference_speed = preconditioner == nullptr
                    ? sound
                    : weiss_smith_reference_speed(
                          speed, sound, viscous_speed, *preconditioner);
                const Real beta
                    = (reference_speed / sound) * (reference_speed / sound);
                Real spectral_sum = 0.0;
                for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                    const auto axis = static_cast<Axis>(logical);
                    for (int side = 0; side <= 1; ++side) {
                        auto face = cell;
                        face[static_cast<std::size_t>(axis)] += side;
                        const auto area = area_vector(metric, axis, face);
                        const Real magnitude = std::sqrt(
                            area[0] * area[0] + area[1] * area[1] + area[2] * area[2]);
                        const Real normal_velocity
                            = (u * area[0] + v * area[1] + w * area[2]) / magnitude;
                        Real face_radius = std::abs(normal_velocity) + sound;
                        if (preconditioner != nullptr) {
                            face_radius = 0.0;
                            for (const Real eigenvalue :
                                 weiss_smith_eigenvalues(normal_velocity, sound, beta)) {
                                face_radius = std::max(face_radius, std::abs(eigenvalue));
                            }
                        }
                        // Scalar spectral splitting retains the convective
                        // direction: A+ ~ (u_n+a)I/2 on the lower side and
                        // -A- ~ (a-u_n)I/2 on the upper side.
                        const Real signed_velocity = side == 0 ? normal_velocity
                                                               : -normal_velocity;
                        Real coupling
                            = 0.5 * (face_radius + signed_velocity) * magnitude / volume;
                        if (additional_face_coupling != nullptr) {
                            const Real extra = (*additional_face_coupling)(
                                i, j, k, 2 * logical + side);
                            if (!std::isfinite(extra) || extra < 0.0) {
                                throw std::invalid_argument(
                                    "additional implicit face coupling is invalid");
                            }
                            coupling += extra;
                        }
                        result.coupling(i, j, k, 2 * logical + side) = coupling;
                        spectral_sum += coupling;
                    }
                }
                result.diagonal(i, j, k, 0)
                    = physical_time_diagonal + spectral_sum * (1.0 + 1.0 / pseudo_cfl);
            }
        }
    }
    return result;
}

Field<Real> build_preconditioned_time_diagonal_blocks(
    const StructuredBlock& block,
    const MetricField& metric,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    const Field<Real>& scalar_diagonal,
    Real physical_time_diagonal,
    const WeissSmithParameters& preconditioner,
    Real viscous_preconditioner_scale)
{
    const auto extent = block.cell_extent();
    if (metric.dimension() != block.cell_dimension()
        || metric.jacobian().interior_extent() != extent
        || scalar_diagonal.interior_extent() != extent || scalar_diagonal.components() != 1
        || scalar_diagonal.ghost_width() != 0 || !std::isfinite(physical_time_diagonal)
        || physical_time_diagonal <= 0.0 || !std::isfinite(viscous_preconditioner_scale)
        || viscous_preconditioner_scale < 0.0) {
        throw std::invalid_argument("invalid preconditioned physical-time diagonal input");
    }
    preconditioner.validate();
    Field<Real> result(extent, euler_components * euler_components, 0, 0.0);
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const Index3 cell {i, j, k};
                TemperaturePrimitiveState temperature_state {};
                for (int component = 0; component < euler_components; ++component) {
                    temperature_state[static_cast<std::size_t>(component)]
                        = block.flow.temperature_primitive(i, j, k, component);
                }
                Real area_sum = 0.0;
                for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                    const auto axis = static_cast<Axis>(logical);
                    for (int side = 0; side <= 1; ++side) {
                        auto face = cell;
                        face[static_cast<std::size_t>(axis)] += side;
                        const auto area = area_vector(metric, axis, face);
                        area_sum += std::sqrt(
                            area[0] * area[0] + area[1] * area[1] + area[2] * area[2]);
                    }
                }
                const Real length = 2.0 * metric.jacobian()(i, j, k) / area_sum;
                if (!std::isfinite(length) || length <= 0.0) {
                    throw PhysicsError("preconditioner characteristic length is invalid");
                }
                const auto pressure_state = pressure_primitive(temperature_state,
                                                               gas,
                                                               reference,
                                                               floors,
                                                               block.cell_dimension());
                const auto matrix = weiss_smith_conservative_preconditioner(
                    pressure_state,
                    gas,
                    floors,
                    preconditioner,
                    viscous_preconditioner_scale / length,
                    block.cell_dimension());
                for (int row = 0; row < euler_components; ++row) {
                    for (int column = 0; column < euler_components; ++column) {
                        const int entry = row * euler_components + column;
                        result(i, j, k, entry)
                            = physical_time_diagonal
                            * matrix[static_cast<std::size_t>(row)]
                                    [static_cast<std::size_t>(column)];
                        if (row == column) {
                            result(i, j, k, entry) += scalar_diagonal(i, j, k, 0);
                        }
                    }
                }
            }
        }
    }
    return result;
}

Field<Real> solve_scalar_lu_sgs(const Field<Real>& right_hand_side,
                                const Field<Real>& diagonal,
                                const Field<Real>& coupling,
                                int dimension,
                                const LuSgsIterationConfig& config)
{
    config.validate();
    validate_system(right_hand_side, diagonal, coupling, dimension);
    const auto extent = right_hand_side.interior_extent();
    Field<Real> solution(extent, right_hand_side.components(), 0, 0.0);
    Field<Real> defect(extent, right_hand_side.components(), 0, 0.0);
    Field<Real> intermediate(extent, right_hand_side.components(), 0, 0.0);
    Field<Real> correction(extent, right_hand_side.components(), 0, 0.0);

    for (int sweep = 0; sweep < config.sweeps; ++sweep) {
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    for (int variable = 0; variable < right_hand_side.components(); ++variable) {
                        const int diagonal_component
                            = diagonal.components() == 1 ? 0 : variable;
                        Real value = right_hand_side(i, j, k, variable)
                            - diagonal(i, j, k, diagonal_component)
                                * solution(i, j, k, variable);
                        for (int axis = 0; axis < dimension; ++axis) {
                            for (const int side : {-1, 1}) {
                                if (has_neighbor(extent, i, j, k, axis, side)) {
                                    value += coupling(i,
                                                      j,
                                                      k,
                                                      side < 0 ? lower_component(axis)
                                                               : upper_component(axis))
                                        * neighbor_value(
                                            solution, i, j, k, variable, axis, side);
                                }
                            }
                        }
                        defect(i, j, k, variable) = value;
                    }
                }
            }
        }
        intermediate.fill(0.0);
        correction.fill(0.0);
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    for (int variable = 0; variable < right_hand_side.components(); ++variable) {
                        Real value = defect(i, j, k, variable);
                        for (int axis = 0; axis < dimension; ++axis) {
                            if (has_neighbor(extent, i, j, k, axis, -1)) {
                                value += coupling(i, j, k, lower_component(axis))
                                    * neighbor_value(
                                        intermediate, i, j, k, variable, axis, -1);
                            }
                        }
                        const int diagonal_component
                            = diagonal.components() == 1 ? 0 : variable;
                        intermediate(i, j, k, variable)
                            = value / diagonal(i, j, k, diagonal_component);
                    }
                }
            }
        }
        for (int k = extent.nk - 1; k >= 0; --k) {
            for (int j = extent.nj - 1; j >= 0; --j) {
                for (int i = extent.ni - 1; i >= 0; --i) {
                    for (int variable = 0; variable < right_hand_side.components(); ++variable) {
                        Real value = intermediate(i, j, k, variable);
                        const int diagonal_component
                            = diagonal.components() == 1 ? 0 : variable;
                        for (int axis = 0; axis < dimension; ++axis) {
                            if (has_neighbor(extent, i, j, k, axis, 1)) {
                                value += coupling(i, j, k, upper_component(axis))
                                    * neighbor_value(
                                          correction, i, j, k, variable, axis, 1)
                                    / diagonal(i, j, k, diagonal_component);
                            }
                        }
                        correction(i, j, k, variable) = value;
                        solution(i, j, k, variable) += config.relaxation * value;
                        if (!std::isfinite(solution(i, j, k, variable))) {
                            throw std::runtime_error("LU-SGS produced a non-finite increment");
                        }
                    }
                }
            }
        }
    }
    return solution;
}

Field<Real> solve_block_lu_sgs(const Field<Real>& right_hand_side,
                               const Field<Real>& diagonal_blocks,
                               const Field<Real>& coupling,
                               int dimension,
                               const LuSgsIterationConfig& config)
{
    config.validate();
    if (right_hand_side.components() <= 0
        || diagonal_blocks.components()
            != right_hand_side.components() * right_hand_side.components()) {
        throw std::invalid_argument("LU-SGS block diagonal size is incompatible");
    }
    validate_system(right_hand_side, Field<Real>(right_hand_side.interior_extent(), 1, 0, 1.0),
                    coupling, dimension);
    if (diagonal_blocks.interior_extent() != right_hand_side.interior_extent()
        || diagonal_blocks.ghost_width() != 0) {
        throw std::invalid_argument("LU-SGS block diagonal metadata are incompatible");
    }
    const auto extent = right_hand_side.interior_extent();
    const int variables = right_hand_side.components();
    Field<Real> solution(extent, variables, 0, 0.0);
    Field<Real> defect(extent, variables, 0, 0.0);
    Field<Real> intermediate(extent, variables, 0, 0.0);
    Field<Real> correction(extent, variables, 0, 0.0);
    // These buffers are reused for every cell.  Allocating them in the sweeps made the
    // implicit path scale its heap traffic with cell count without changing the algebra.
    std::vector<Real> product(static_cast<std::size_t>(variables));
    std::vector<Real> value(static_cast<std::size_t>(variables));
    std::vector<Real> solved(static_cast<std::size_t>(variables));
    std::vector<Real> matrix(static_cast<std::size_t>(variables * variables));

    for (int sweep = 0; sweep < config.sweeps; ++sweep) {
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 cell {i, j, k};
                    multiply_dense_block(diagonal_blocks, cell, solution, product);
                    for (int variable = 0; variable < variables; ++variable) {
                        Real value = right_hand_side(i, j, k, variable)
                            - product[static_cast<std::size_t>(variable)];
                        for (int axis = 0; axis < dimension; ++axis) {
                            for (const int side : {-1, 1}) {
                                if (has_neighbor(extent, i, j, k, axis, side)) {
                                    value += coupling(i,
                                                      j,
                                                      k,
                                                      side < 0 ? lower_component(axis)
                                                               : upper_component(axis))
                                        * neighbor_value(
                                            solution, i, j, k, variable, axis, side);
                                }
                            }
                        }
                        defect(i, j, k, variable) = value;
                    }
                }
            }
        }
        intermediate.fill(0.0);
        correction.fill(0.0);
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 cell {i, j, k};
                    for (int variable = 0; variable < variables; ++variable) {
                        value[static_cast<std::size_t>(variable)]
                            = defect(i, j, k, variable);
                        for (int axis = 0; axis < dimension; ++axis) {
                            if (has_neighbor(extent, i, j, k, axis, -1)) {
                                value[static_cast<std::size_t>(variable)]
                                    += coupling(i, j, k, lower_component(axis))
                                    * neighbor_value(
                                        intermediate, i, j, k, variable, axis, -1);
                            }
                        }
                    }
                    solve_dense_block(diagonal_blocks, cell, value, matrix, solved);
                    for (int variable = 0; variable < variables; ++variable) {
                        intermediate(i, j, k, variable)
                            = solved[static_cast<std::size_t>(variable)];
                    }
                }
            }
        }
        for (int k = extent.nk - 1; k >= 0; --k) {
            for (int j = extent.nj - 1; j >= 0; --j) {
                for (int i = extent.ni - 1; i >= 0; --i) {
                    const Index3 cell {i, j, k};
                    multiply_dense_block(diagonal_blocks, cell, intermediate, value);
                    for (int variable = 0; variable < variables; ++variable) {
                        for (int axis = 0; axis < dimension; ++axis) {
                            if (has_neighbor(extent, i, j, k, axis, 1)) {
                                value[static_cast<std::size_t>(variable)]
                                    += coupling(i, j, k, upper_component(axis))
                                    * neighbor_value(
                                        correction, i, j, k, variable, axis, 1);
                            }
                        }
                    }
                    solve_dense_block(diagonal_blocks, cell, value, matrix, solved);
                    for (int variable = 0; variable < variables; ++variable) {
                        const Real delta = solved[static_cast<std::size_t>(variable)];
                        correction(i, j, k, variable) = delta;
                        solution(i, j, k, variable) += config.relaxation * delta;
                        if (!std::isfinite(solution(i, j, k, variable))) {
                            throw std::runtime_error("LU-SGS produced a non-finite increment");
                        }
                    }
                }
            }
        }
    }
    return solution;
}

Field<Real> build_euler_face_coupling_blocks(
    const StructuredBlock& block,
    const MetricField& metric,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    const Field<Real>* additional_face_coupling)
{
    const auto extent = block.cell_extent();
    const int dimension = block.cell_dimension();
    if (metric.dimension() != dimension || metric.jacobian().interior_extent() != extent) {
        throw std::invalid_argument("Euler face-block metric is incompatible");
    }
    if (additional_face_coupling != nullptr
        && (additional_face_coupling->interior_extent() != extent
            || additional_face_coupling->components() != 2 * dimension
            || additional_face_coupling->ghost_width() != 0)) {
        throw std::invalid_argument("Euler face-block viscous coupling is incompatible");
    }
    Field<Real> result(extent,
                       2 * dimension * euler_components * euler_components,
                       0,
                       0.0);
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                TemperaturePrimitiveState temperature_state {};
                for (int component = 0; component < euler_components; ++component) {
                    temperature_state[static_cast<std::size_t>(component)]
                        = block.flow.temperature_primitive(i, j, k, component);
                }
                const auto pressure_state = pressure_primitive(
                    temperature_state, gas, reference, floors, dimension);
                const Real sound = thermodynamic_sound_speed(
                    temperature_state, gas, reference, floors, dimension);
                const Real volume = metric.jacobian()(i, j, k);
                const Index3 cell {i, j, k};
                for (int logical = 0; logical < dimension; ++logical) {
                    const auto axis = static_cast<Axis>(logical);
                    for (int side = 0; side <= 1; ++side) {
                        auto face = cell;
                        face[static_cast<std::size_t>(axis)] += side;
                        const auto area = area_vector(metric, axis, face);
                        const Real magnitude = std::sqrt(
                            area[0] * area[0] + area[1] * area[1] + area[2] * area[2]);
                        const Normal3 normal {
                            area[0] / magnitude,
                            area[1] / magnitude,
                            area[2] / magnitude,
                        };
                        const Real normal_velocity = pressure_state[1] * normal.x
                            + pressure_state[2] * normal.y + pressure_state[3] * normal.z;
                        const Real radius = std::abs(normal_velocity) + sound;
                        const auto jacobian
                            = euler_normal_flux_jacobian(pressure_state, normal, gas.gamma());
                        const Real scale = 0.5 * magnitude / volume;
                        const Real diffusion = additional_face_coupling == nullptr
                            ? 0.0
                            : (*additional_face_coupling)(i, j, k, 2 * logical + side);
                        if (!std::isfinite(diffusion) || diffusion < 0.0) {
                            throw std::invalid_argument(
                                "Euler face-block viscous coupling is invalid");
                        }
                        const int face_component = 2 * logical + side;
                        const int offset
                            = face_component * euler_components * euler_components;
                        for (int row = 0; row < euler_components; ++row) {
                            for (int column = 0; column < euler_components; ++column) {
                                const Real identity = row == column ? radius : 0.0;
                                const Real split = side == 0
                                    ? identity
                                        + jacobian[static_cast<std::size_t>(row)]
                                                  [static_cast<std::size_t>(column)]
                                    : identity
                                        - jacobian[static_cast<std::size_t>(row)]
                                                  [static_cast<std::size_t>(column)];
                                result(i, j, k, offset + row * euler_components + column)
                                    = scale * split
                                    + (row == column ? diffusion : 0.0);
                            }
                        }
                    }
                }
            }
        }
    }
    return result;
}

Field<Real> solve_face_block_lu_sgs(const Field<Real>& right_hand_side,
                                    const Field<Real>& diagonal,
                                    const Field<Real>& face_blocks,
                                    int dimension,
                                    const LuSgsIterationConfig& config)
{
    config.validate();
    const auto extent = right_hand_side.interior_extent();
    const int variables = right_hand_side.components();
    if (dimension < 1 || dimension > 3 || variables <= 0
        || diagonal.interior_extent() != extent || diagonal.ghost_width() != 0
        || (diagonal.components() != 1 && diagonal.components() != variables)
        || face_blocks.interior_extent() != extent || face_blocks.ghost_width() != 0
        || face_blocks.components() != 2 * dimension * variables * variables) {
        throw std::invalid_argument("LU-SGS face-block system metadata are incompatible");
    }
    Field<Real> solution(extent, variables, 0, 0.0);
    Field<Real> defect(extent, variables, 0, 0.0);
    Field<Real> intermediate(extent, variables, 0, 0.0);
    Field<Real> correction(extent, variables, 0, 0.0);
    std::vector<Real> product(static_cast<std::size_t>(variables));
    std::vector<Real> value(static_cast<std::size_t>(variables));
    for (int sweep = 0; sweep < config.sweeps; ++sweep) {
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 cell {i, j, k};
                    for (int variable = 0; variable < variables; ++variable) {
                        const int diagonal_component
                            = diagonal.components() == 1 ? 0 : variable;
                        defect(i, j, k, variable) = right_hand_side(i, j, k, variable)
                            - diagonal(i, j, k, diagonal_component)
                                * solution(i, j, k, variable);
                    }
                    for (int axis = 0; axis < dimension; ++axis) {
                        for (const int side : {-1, 1}) {
                            if (!has_neighbor(extent, i, j, k, axis, side)) continue;
                            auto neighbor = cell;
                            neighbor[static_cast<std::size_t>(axis)] += side;
                            multiply_face_block(
                                face_blocks,
                                cell,
                                side < 0 ? lower_component(axis) : upper_component(axis),
                                solution,
                                neighbor,
                                product);
                            for (int variable = 0; variable < variables; ++variable) {
                                defect(i, j, k, variable)
                                    += product[static_cast<std::size_t>(variable)];
                            }
                        }
                    }
                }
            }
        }
        intermediate.fill(0.0);
        correction.fill(0.0);
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 cell {i, j, k};
                    for (int variable = 0; variable < variables; ++variable) {
                        value[static_cast<std::size_t>(variable)]
                            = defect(i, j, k, variable);
                    }
                    for (int axis = 0; axis < dimension; ++axis) {
                        if (!has_neighbor(extent, i, j, k, axis, -1)) continue;
                        auto neighbor = cell;
                        --neighbor[static_cast<std::size_t>(axis)];
                        multiply_face_block(
                            face_blocks,
                            cell,
                            lower_component(axis),
                            intermediate,
                            neighbor,
                            product);
                        for (int variable = 0; variable < variables; ++variable) {
                            value[static_cast<std::size_t>(variable)]
                                += product[static_cast<std::size_t>(variable)];
                        }
                    }
                    for (int variable = 0; variable < variables; ++variable) {
                        const int diagonal_component
                            = diagonal.components() == 1 ? 0 : variable;
                        intermediate(i, j, k, variable)
                            = value[static_cast<std::size_t>(variable)]
                            / diagonal(i, j, k, diagonal_component);
                    }
                }
            }
        }
        for (int k = extent.nk - 1; k >= 0; --k) {
            for (int j = extent.nj - 1; j >= 0; --j) {
                for (int i = extent.ni - 1; i >= 0; --i) {
                    const Index3 cell {i, j, k};
                    for (int variable = 0; variable < variables; ++variable) {
                        const int diagonal_component
                            = diagonal.components() == 1 ? 0 : variable;
                        value[static_cast<std::size_t>(variable)]
                            = diagonal(i, j, k, diagonal_component)
                            * intermediate(i, j, k, variable);
                    }
                    for (int axis = 0; axis < dimension; ++axis) {
                        if (!has_neighbor(extent, i, j, k, axis, 1)) continue;
                        auto neighbor = cell;
                        ++neighbor[static_cast<std::size_t>(axis)];
                        multiply_face_block(
                            face_blocks,
                            cell,
                            upper_component(axis),
                            correction,
                            neighbor,
                            product);
                        for (int variable = 0; variable < variables; ++variable) {
                            value[static_cast<std::size_t>(variable)]
                                += product[static_cast<std::size_t>(variable)];
                        }
                    }
                    for (int variable = 0; variable < variables; ++variable) {
                        const int diagonal_component
                            = diagonal.components() == 1 ? 0 : variable;
                        const Real delta = value[static_cast<std::size_t>(variable)]
                            / diagonal(i, j, k, diagonal_component);
                        correction(i, j, k, variable) = delta;
                        solution(i, j, k, variable) += config.relaxation * delta;
                        if (!std::isfinite(solution(i, j, k, variable))) {
                            throw std::runtime_error(
                                "LU-SGS face-block solve produced a non-finite increment");
                        }
                    }
                }
            }
        }
    }
    return solution;
}

Real bdf_time_diagonal(BdfOrder order, Real physical_time_step)
{
    if (!std::isfinite(physical_time_step) || physical_time_step <= 0.0) {
        throw std::invalid_argument("BDF physical time step must be positive and finite");
    }
    if (order == BdfOrder::First) return 1.0 / physical_time_step;
    if (order == BdfOrder::Second) return 1.5 / physical_time_step;
    throw std::invalid_argument("unknown BDF order");
}

Real bdf_time_residual(BdfOrder order,
                       Real candidate,
                       Real current,
                       Real previous,
                       Real physical_time_step)
{
    if (!std::isfinite(candidate) || !std::isfinite(current) || !std::isfinite(previous)) {
        throw std::invalid_argument("BDF state contains a non-finite value");
    }
    static_cast<void>(bdf_time_diagonal(order, physical_time_step));
    if (order == BdfOrder::First) return (candidate - current) / physical_time_step;
    if (order == BdfOrder::Second) {
        return (3.0 * candidate - 4.0 * current + previous)
            / (2.0 * physical_time_step);
    }
    throw std::invalid_argument("unknown BDF order");
}

} // namespace wcns
