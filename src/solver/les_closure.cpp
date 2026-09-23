#include <wcns/solver/les_closure.hpp>

#include <wcns/solver/les_filter.hpp>
#include <wcns/solver/les_models.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace wcns {
namespace {

constexpr int resolved_density = 0;
constexpr int resolved_velocity = 1;
constexpr int resolved_gradient = 4;
constexpr int resolved_delta = 13;

Real stencil_weight(int di, int dj, int dk)
{
    return les_box3_weights[static_cast<std::size_t>(di + 1)]
        * les_box3_weights[static_cast<std::size_t>(dj + 1)]
        * les_box3_weights[static_cast<std::size_t>(dk + 1)];
}

SymmetricStress matrix_to_symmetric(const std::array<std::array<Real, 3>, 3>& matrix)
{
    return {matrix[0][0],
            matrix[1][1],
            matrix[2][2],
            0.5 * (matrix[0][1] + matrix[1][0]),
            0.5 * (matrix[0][2] + matrix[2][0]),
            0.5 * (matrix[1][2] + matrix[2][1])};
}

VelocityGradient resolved_gradient_at(const Field<Real>& resolved, Index3 cell)
{
    VelocityGradient result {};
    for (int velocity = 0; velocity < 3; ++velocity) {
        for (int direction = 0; direction < 3; ++direction) {
            result[static_cast<std::size_t>(velocity)][static_cast<std::size_t>(direction)]
                = les_fixed_stencil_value(
                    resolved, cell, resolved_gradient + 3 * velocity + direction);
        }
    }
    return result;
}

SymmetricStress leonard_at(const Field<Real>& resolved, Index3 cell)
{
    const auto sample = favre_box3_tensor_sample(
        resolved, cell, resolved_density, {{1, 2, 3}});
    return les_deviatoric_part(matrix_to_symmetric(sample.leonard_stress));
}

SymmetricStress germano_m_at(const Field<Real>& resolved,
                             Index3 center,
                             const TurbulenceModelConfig& config)
{
    const auto favre = favre_box3_tensor_sample(
        resolved, center, resolved_density, {{1, 2, 3}});
    VelocityGradient test_gradient {};
    SymmetricStress filtered_term;
    for (int dk = -1; dk <= 1; ++dk) {
        for (int dj = -1; dj <= 1; ++dj) {
            for (int di = -1; di <= 1; ++di) {
                const Index3 cell {center.i + di, center.j + dj, center.k + dk};
                const Real weight = stencil_weight(di, dj, dk);
                const Real density
                    = les_fixed_stencil_value(resolved, cell, resolved_density);
                const auto gradient = resolved_gradient_at(resolved, cell);
                for (int velocity = 0; velocity < 3; ++velocity) {
                    for (int direction = 0; direction < 3; ++direction) {
                        test_gradient[static_cast<std::size_t>(velocity)]
                                     [static_cast<std::size_t>(direction)]
                            += weight * density
                            * gradient[static_cast<std::size_t>(velocity)]
                                      [static_cast<std::size_t>(direction)]
                            / favre.density;
                    }
                }
                const auto strain = les_deviatoric_strain(gradient);
                const Real magnitude = les_strain_magnitude(strain);
                const Real delta = les_fixed_stencil_value(resolved, cell, resolved_delta);
                const Real factor = weight * density * delta * delta * magnitude;
                filtered_term.xx += factor * strain.xx;
                filtered_term.yy += factor * strain.yy;
                filtered_term.zz += factor * strain.zz;
                filtered_term.xy += factor * strain.xy;
                filtered_term.xz += factor * strain.xz;
                filtered_term.yz += factor * strain.yz;
            }
        }
    }
    const auto test_strain = les_deviatoric_strain(test_gradient);
    const Real test_magnitude = les_strain_magnitude(test_strain);
    const Real delta = les_fixed_stencil_value(resolved, center, resolved_delta);
    const Real test_delta = config.les_test_filter_ratio * delta;
    const Real test_factor = favre.density * test_delta * test_delta * test_magnitude;
    return {-2.0 * (test_factor * test_strain.xx - filtered_term.xx),
            -2.0 * (test_factor * test_strain.yy - filtered_term.yy),
            -2.0 * (test_factor * test_strain.zz - filtered_term.zz),
            -2.0 * (test_factor * test_strain.xy - filtered_term.xy),
            -2.0 * (test_factor * test_strain.xz - filtered_term.xz),
            -2.0 * (test_factor * test_strain.yz - filtered_term.yz)};
}

void require_compatible(const Field<Real>& field, int components, const char* label)
{
    if (field.components() != components || field.ghost_width() < 1
        || field.interior_extent().ni <= 0 || field.interior_extent().nj <= 0
        || field.interior_extent().nk <= 0) {
        throw std::invalid_argument(std::string(label) + " field metadata is invalid");
    }
}

} // namespace

void populate_les_resolved_field(Field<Real>& resolved,
                                 const StructuredBlock& block,
                                 const MetricField& metric,
                                 const PrimitiveGradientField& gradients,
                                 const TurbulenceModelConfig& config)
{
    config.validate();
    require_compatible(resolved, les_resolved_components, "LES resolved");
    if (block.cell_dimension() != 3 || metric.dimension() != 3 || gradients.dimension() != 3
        || resolved.interior_extent() != block.cell_extent()
        || metric.jacobian().interior_extent() != block.cell_extent()
        || gradients.values().interior_extent() != block.cell_extent()) {
        throw std::invalid_argument("LES resolved inputs require matching three-dimensional fields");
    }
    resolved.fill(std::numeric_limits<Real>::quiet_NaN());
    const auto extent = block.cell_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const Real density
                    = block.flow.temperature_primitive(i, j, k, temperature_density);
                const Real volume = metric.jacobian()(i, j, k);
                if (!std::isfinite(density) || density <= 0.0 || !std::isfinite(volume)
                    || volume <= 0.0) {
                    throw PhysicsError("LES resolved field has invalid density or volume");
                }
                resolved(i, j, k, resolved_density) = density;
                for (int velocity = 0; velocity < 3; ++velocity) {
                    resolved(i, j, k, resolved_velocity + velocity)
                        = block.flow.temperature_primitive(
                            i, j, k, temperature_velocity_x + velocity);
                    for (int direction = 0; direction < 3; ++direction) {
                        resolved(i,
                                 j,
                                 k,
                                 resolved_gradient + 3 * velocity + direction)
                            = gradients({i, j, k},
                                        static_cast<ViscousPrimitive>(velocity),
                                        direction);
                    }
                }
                resolved(i, j, k, resolved_delta)
                    = config.les_filter_width_ratio * std::cbrt(volume);
            }
        }
    }
}

void complete_les_fixed_stencil_ghosts(Field<Real>& field)
{
    if (field.ghost_width() < 1) {
        throw std::invalid_argument("LES fixed stencil requires at least one ghost layer");
    }
    const auto extent = field.interior_extent();
    const int ghost = field.ghost_width();
    const Field<Real> exchanged = field;
    // Halo exchange supplies connected face ghosts but not their edge/corner
    // extensions.  Extend one physical direction per pass so a valid exchanged
    // face value is preserved at a partition/physical-boundary intersection.
    for (int pass = 0; pass < 3; ++pass) {
        for (int k = -ghost; k < extent.nk + ghost; ++k) {
            for (int j = -ghost; j < extent.nj + ghost; ++j) {
                for (int i = -ghost; i < extent.ni + ghost; ++i) {
                    const Index3 cell {i, j, k};
                    for (int component = 0; component < field.components(); ++component) {
                        Real& value = field(i, j, k, component);
                        if (std::isfinite(value)) continue;
                        for (const bool connected_axis : {false, true}) {
                            for (int axis = 0; axis < 3; ++axis) {
                                const auto direction = static_cast<std::size_t>(axis);
                                const int upper = extent[direction] - 1;
                                if (cell[direction] >= 0 && cell[direction] <= upper) continue;
                                auto anchor = cell;
                                for (int other = 0; other < 3; ++other) {
                                    if (other == axis) continue;
                                    const auto other_direction
                                        = static_cast<std::size_t>(other);
                                    anchor[other_direction] = std::clamp(
                                        anchor[other_direction],
                                        0,
                                        extent[other_direction] - 1);
                                }
                                const bool axis_has_exchange = std::isfinite(exchanged(
                                    anchor.i, anchor.j, anchor.k, component));
                                if (axis_has_exchange != connected_axis) continue;
                                auto source = cell;
                                source[direction]
                                    = std::clamp(cell[direction], 0, upper);
                                const Real candidate
                                    = field(source.i, source.j, source.k, component);
                                if (std::isfinite(candidate)) {
                                    value = candidate;
                                    break;
                                }
                            }
                            if (std::isfinite(value)) break;
                        }
                    }
                }
            }
        }
    }
    for (int k = -ghost; k < extent.nk + ghost; ++k) {
        for (int j = -ghost; j < extent.nj + ghost; ++j) {
            for (int i = -ghost; i < extent.ni + ghost; ++i) {
                for (int component = 0; component < field.components(); ++component) {
                    if (!std::isfinite(field(i, j, k, component))) {
                        throw PhysicsError("LES fixed-stencil ghost extension is non-finite");
                    }
                }
            }
        }
    }
}

void compute_les_dynamic_moments(Field<Real>& moments,
                                 const Field<Real>& resolved,
                                 const TurbulenceModelConfig& config)
{
    config.validate();
    require_compatible(moments, les_dynamic_moment_components, "LES dynamic moment");
    require_compatible(resolved, les_resolved_components, "LES resolved");
    if (moments.interior_extent() != resolved.interior_extent()) {
        throw std::invalid_argument("LES dynamic moment extent differs from resolved field");
    }
    moments.fill(std::numeric_limits<Real>::quiet_NaN());
    const auto extent = resolved.interior_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const Index3 cell {i, j, k};
                const auto leonard = leonard_at(resolved, cell);
                const auto germano = germano_m_at(resolved, cell, config);
                moments(i, j, k, 0) = les_symmetric_inner_product(leonard, germano);
                moments(i, j, k, 1) = les_symmetric_inner_product(germano, germano);
            }
        }
    }
}

void compute_les_closure_field(Field<Real>& closure,
                               const Field<Real>& resolved,
                               const Field<Real>& dynamic_moments,
                               const TurbulenceModelConfig& config)
{
    config.validate();
    require_compatible(closure, les_closure_components, "LES closure");
    require_compatible(resolved, les_resolved_components, "LES resolved");
    require_compatible(
        dynamic_moments, les_dynamic_moment_components, "LES dynamic moment");
    if (closure.interior_extent() != resolved.interior_extent()
        || dynamic_moments.interior_extent() != resolved.interior_extent()) {
        throw std::invalid_argument("LES closure field extents differ");
    }
    closure.fill(std::numeric_limits<Real>::quiet_NaN());
    const auto extent = resolved.interior_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const Index3 cell {i, j, k};
                const auto leonard = leonard_at(resolved, cell);
                closure(i, j, k, les_leonard_xx) = leonard.xx;
                closure(i, j, k, les_leonard_yy) = leonard.yy;
                closure(i, j, k, les_leonard_zz) = leonard.zz;
                closure(i, j, k, les_leonard_xy) = leonard.xy;
                closure(i, j, k, les_leonard_xz) = leonard.xz;
                closure(i, j, k, les_leonard_yz) = leonard.yz;
                closure(i, j, k, les_dynamic_coefficient)
                    = evaluate_dynamic_smagorinsky_coefficient(
                        les_box3_tensor_sample(dynamic_moments, cell, 0),
                        les_box3_tensor_sample(dynamic_moments, cell, 1),
                        config.les_dynamic_denominator_floor,
                        config.les_dynamic_coefficient_minimum,
                        config.les_dynamic_coefficient_maximum);
                closure(i, j, k, les_filter_width)
                    = les_fixed_stencil_value(resolved, cell, resolved_delta);
            }
        }
    }
}

Real interpolate_les_closure_face(const Field<Real>& closure,
                                  Axis axis,
                                  Index3 face,
                                  int component)
{
    if (component < 0 || component >= closure.components()) {
        throw std::invalid_argument("LES face interpolation component is invalid");
    }
    auto left = face;
    --left[static_cast<std::size_t>(axis)];
    const Real value = 0.5
        * (les_fixed_stencil_value(closure, left, component)
           + les_fixed_stencil_value(closure, face, component));
    if (!std::isfinite(value)) {
        throw PhysicsError("LES face interpolation produced a non-finite value");
    }
    return value;
}

} // namespace wcns
