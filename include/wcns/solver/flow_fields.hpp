#pragma once

#include <wcns/core/field.hpp>
#include <wcns/core/types.hpp>

#include <cstdint>
#include <limits>

namespace wcns {

inline constexpr int euler_components = 5;

struct ImplicitTimeHistory {
    bool valid = false;
    Real physical_time_step = 0.0;
    Field<Real> previous_mean;
    Field<Real> previous_model;

    void prepare(Extent3 extent, int model_components)
    {
        previous_mean = Field<Real>(extent, euler_components, 0, 0.0);
        previous_model = model_components > 0
            ? Field<Real>(extent, model_components, 0, 0.0)
            : Field<Real> {};
    }
};

struct FlowFields {
    FlowFields(Extent3 cell_extent, int ghost_width)
        : conservative(cell_extent, euler_components, ghost_width)
        , primitive(cell_extent, euler_components, ghost_width)
        , temperature_primitive(cell_extent, euler_components, ghost_width)
        , residual(cell_extent, euler_components, 0)
    {
        const Real nan = std::numeric_limits<Real>::quiet_NaN();
        conservative.fill(nan);
        primitive.fill(nan);
        temperature_primitive.fill(nan);
        residual.fill(0.0);
    }

    Field<Real> conservative;
    Field<Real> primitive;
    Field<Real> temperature_primitive;
    Field<Real> residual;
    ImplicitTimeHistory implicit_history;
    std::uint64_t physical_ghost_version = 0;
};

} // namespace wcns
