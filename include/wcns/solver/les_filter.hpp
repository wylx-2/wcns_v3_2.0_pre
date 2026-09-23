#pragma once

#include <wcns/core/field.hpp>
#include <wcns/core/types.hpp>

#include <array>

namespace wcns {

inline constexpr std::array<Real, 3> les_box3_weights {{0.25, 0.5, 0.25}};

struct FavreFilterSample {
    Real density = 0.0;
    std::array<Real, 3> velocity {{0.0, 0.0, 0.0}};
    std::array<std::array<Real, 3>, 3> second_moment {};
    std::array<std::array<Real, 3>, 3> leonard_stress {};
};

// A non-finite or unavailable physical ghost is replaced by the nearest
// interior cell. Connected/MPI ghosts remain part of the fixed stencil.
[[nodiscard]] Real les_box3_tensor_sample(const Field<Real>& field,
                                          Index3 center,
                                          int component);

[[nodiscard]] Real les_fixed_stencil_value(const Field<Real>& field,
                                           Index3 index,
                                           int component);

void les_box3_tensor_filter(const Field<Real>& input, Field<Real>& output);

[[nodiscard]] FavreFilterSample
favre_box3_tensor_sample(const Field<Real>& primitive,
                         Index3 center,
                         int density_component,
                         std::array<int, 3> velocity_components);

[[nodiscard]] Real les_box3_fourier_response(std::array<Real, 3> wave_angles);

} // namespace wcns
