#include <wcns/solver/les_filter.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

Index3 nearest_interior(Index3 index, Extent3 extent)
{
    index.i = std::clamp(index.i, 0, extent.ni - 1);
    index.j = std::clamp(index.j, 0, extent.nj - 1);
    index.k = std::clamp(index.k, 0, extent.nk - 1);
    return index;
}

bool stored(const Field<Real>& field, Index3 index)
{
    const auto extent = field.interior_extent();
    const int ghost = field.ghost_width();
    return index.i >= -ghost && index.i < extent.ni + ghost && index.j >= -ghost
        && index.j < extent.nj + ghost && index.k >= -ghost
        && index.k < extent.nk + ghost;
}

Real extended_value(const Field<Real>& field, Index3 index, int component)
{
    if (stored(field, index)) {
        const Real candidate = field(index.i, index.j, index.k, component);
        if (std::isfinite(candidate)) return candidate;
    }
    const auto cell = nearest_interior(index, field.interior_extent());
    const Real candidate = field(cell.i, cell.j, cell.k, component);
    if (!std::isfinite(candidate)) {
        throw std::invalid_argument("LES filter input contains a non-finite interior value");
    }
    return candidate;
}

template <class Sample> Real filtered_sample(Sample&& sample)
{
    Real result = 0.0;
    for (int dk = -1; dk <= 1; ++dk) {
        for (int dj = -1; dj <= 1; ++dj) {
            for (int di = -1; di <= 1; ++di) {
                const Real weight = les_box3_weights[static_cast<std::size_t>(di + 1)]
                    * les_box3_weights[static_cast<std::size_t>(dj + 1)]
                    * les_box3_weights[static_cast<std::size_t>(dk + 1)];
                result += weight * sample(di, dj, dk);
            }
        }
    }
    return result;
}

} // namespace

Real les_box3_tensor_sample(const Field<Real>& field, Index3 center, int component)
{
    if (field.interior_extent().ni <= 0 || field.interior_extent().nj <= 0
        || field.interior_extent().nk <= 0 || component < 0
        || component >= field.components() || center.i < 0
        || center.i >= field.interior_extent().ni || center.j < 0
        || center.j >= field.interior_extent().nj || center.k < 0
        || center.k >= field.interior_extent().nk) {
        throw std::invalid_argument("LES box filter sample metadata is invalid");
    }
    return filtered_sample([&](int di, int dj, int dk) {
        return extended_value(
            field, {center.i + di, center.j + dj, center.k + dk}, component);
    });
}

Real les_fixed_stencil_value(const Field<Real>& field, Index3 index, int component)
{
    if (component < 0 || component >= field.components()) {
        throw std::invalid_argument("LES fixed-stencil component is invalid");
    }
    return extended_value(field, index, component);
}

void les_box3_tensor_filter(const Field<Real>& input, Field<Real>& output)
{
    if (input.interior_extent() != output.interior_extent()
        || input.components() != output.components()) {
        throw std::invalid_argument("LES box filter fields have incompatible metadata");
    }
    const auto extent = input.interior_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                for (int component = 0; component < input.components(); ++component) {
                    output(i, j, k, component)
                        = les_box3_tensor_sample(input, {i, j, k}, component);
                }
            }
        }
    }
}

FavreFilterSample favre_box3_tensor_sample(const Field<Real>& primitive,
                                           Index3 center,
                                           int density_component,
                                           std::array<int, 3> velocity_components)
{
    if (density_component < 0 || density_component >= primitive.components()
        || std::any_of(velocity_components.begin(), velocity_components.end(), [&](int entry) {
               return entry < 0 || entry >= primitive.components();
           })) {
        throw std::invalid_argument("Favre filter component map is invalid");
    }
    FavreFilterSample result;
    result.density = filtered_sample([&](int di, int dj, int dk) {
        return extended_value(primitive,
                              {center.i + di, center.j + dj, center.k + dk},
                              density_component);
    });
    if (!std::isfinite(result.density) || result.density <= 0.0) {
        throw std::invalid_argument("Favre filter density must be positive and finite");
    }
    for (int row = 0; row < 3; ++row) {
        const Real momentum = filtered_sample([&](int di, int dj, int dk) {
            const Index3 index {center.i + di, center.j + dj, center.k + dk};
            return extended_value(primitive, index, density_component)
                * extended_value(
                    primitive, index, velocity_components[static_cast<std::size_t>(row)]);
        });
        result.velocity[static_cast<std::size_t>(row)] = momentum / result.density;
    }
    for (int row = 0; row < 3; ++row) {
        for (int column = 0; column < 3; ++column) {
            const Real moment = filtered_sample([&](int di, int dj, int dk) {
                const Index3 index {center.i + di, center.j + dj, center.k + dk};
                return extended_value(primitive, index, density_component)
                    * extended_value(
                        primitive, index, velocity_components[static_cast<std::size_t>(row)])
                    * extended_value(
                        primitive, index, velocity_components[static_cast<std::size_t>(column)]);
            });
            result.second_moment[static_cast<std::size_t>(row)]
                                [static_cast<std::size_t>(column)]
                = moment;
            result.leonard_stress[static_cast<std::size_t>(row)]
                                  [static_cast<std::size_t>(column)]
                = moment
                - result.density * result.velocity[static_cast<std::size_t>(row)]
                    * result.velocity[static_cast<std::size_t>(column)];
        }
    }
    return result;
}

Real les_box3_fourier_response(std::array<Real, 3> wave_angles)
{
    Real response = 1.0;
    for (const Real angle : wave_angles) {
        if (!std::isfinite(angle)) {
            throw std::invalid_argument("LES Fourier response angle is non-finite");
        }
        response *= 0.5 * (1.0 + std::cos(angle));
    }
    return response;
}

} // namespace wcns
