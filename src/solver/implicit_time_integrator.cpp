#include <wcns/solver/implicit_time_integrator.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>

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
    if (sweeps != 1) {
        throw std::invalid_argument("stage AA supports exactly one LU-SGS forward/backward sweep");
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
    const Field<Real>* additional_cell_spectral_radius)
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
    if (additional_cell_spectral_radius != nullptr
        && (additional_cell_spectral_radius->interior_extent() != block.cell_extent()
            || additional_cell_spectral_radius->components() != 1
            || additional_cell_spectral_radius->ghost_width() != 0)) {
        throw std::invalid_argument("additional implicit spectral field is incompatible");
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
                        const Real coupling = 0.5 * face_radius * magnitude / volume;
                        result.coupling(i, j, k, 2 * logical + side) = coupling;
                        spectral_sum += coupling;
                    }
                }
                if (additional_cell_spectral_radius != nullptr) {
                    const Real extra = (*additional_cell_spectral_radius)(i, j, k, 0);
                    if (!std::isfinite(extra) || extra < 0.0) {
                        throw std::invalid_argument(
                            "additional implicit spectral radius is invalid");
                    }
                    spectral_sum += extra;
                }
                result.diagonal(i, j, k, 0)
                    = physical_time_diagonal + spectral_sum * (1.0 + 1.0 / pseudo_cfl);
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
    Field<Real> intermediate(extent, right_hand_side.components(), 0, 0.0);
    Field<Real> increment(extent, right_hand_side.components(), 0, 0.0);

    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                for (int variable = 0; variable < right_hand_side.components(); ++variable) {
                    Real value = right_hand_side(i, j, k, variable);
                    for (int axis = 0; axis < dimension; ++axis) {
                        if (has_neighbor(extent, i, j, k, axis, -1)) {
                            value += coupling(i, j, k, lower_component(axis))
                                * neighbor_value(intermediate, i, j, k, variable, axis, -1);
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
                    Real correction = 0.0;
                    for (int axis = 0; axis < dimension; ++axis) {
                        if (has_neighbor(extent, i, j, k, axis, 1)) {
                            correction += coupling(i, j, k, upper_component(axis))
                                * neighbor_value(increment, i, j, k, variable, axis, 1);
                        }
                    }
                    const int diagonal_component
                        = diagonal.components() == 1 ? 0 : variable;
                    increment(i, j, k, variable)
                        = config.relaxation
                        * (intermediate(i, j, k, variable)
                           + correction / diagonal(i, j, k, diagonal_component));
                    if (!std::isfinite(increment(i, j, k, variable))) {
                        throw std::runtime_error("LU-SGS produced a non-finite increment");
                    }
                }
            }
        }
    }
    return increment;
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
