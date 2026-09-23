#pragma once

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/mesh/structured_block.hpp>
#include <wcns/solver/turbulence_model.hpp>
#include <wcns/solver/viscous_gradient.hpp>

namespace wcns {

inline constexpr int les_resolved_components = 14;
inline constexpr int les_dynamic_moment_components = 2;
inline constexpr int les_closure_components = 8;

enum LesClosureComponent {
    les_leonard_xx = 0,
    les_leonard_yy = 1,
    les_leonard_zz = 2,
    les_leonard_xy = 3,
    les_leonard_xz = 4,
    les_leonard_yz = 5,
    les_dynamic_coefficient = 6,
    les_filter_width = 7,
};

void populate_les_resolved_field(Field<Real>& resolved,
                                 const StructuredBlock& block,
                                 const MetricField& metric,
                                 const PrimitiveGradientField& gradients,
                                 const TurbulenceModelConfig& config);

void complete_les_fixed_stencil_ghosts(Field<Real>& field);

void compute_les_dynamic_moments(Field<Real>& moments,
                                 const Field<Real>& resolved,
                                 const TurbulenceModelConfig& config);

void compute_les_closure_field(Field<Real>& closure,
                               const Field<Real>& resolved,
                               const Field<Real>& dynamic_moments,
                               const TurbulenceModelConfig& config);

[[nodiscard]] Real interpolate_les_closure_face(const Field<Real>& closure,
                                                Axis axis,
                                                Index3 face,
                                                int component);

} // namespace wcns
