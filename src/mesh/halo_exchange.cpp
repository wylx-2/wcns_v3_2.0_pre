#include <wcns/mesh/halo_exchange.hpp>

#include <cstdlib>

namespace wcns {
namespace {

int axis_index(Axis axis)
{
    return static_cast<int>(axis);
}

int outward_step(Side side)
{
    return side == Side::Lower ? -1 : 1;
}

int inward_step(Side side)
{
    return -outward_step(side);
}

template <class Range>
void validate_cell_range(const Range& range, Extent3 extent, int dimension, const char* label)
{
    for (int axis = 0; axis < dimension; ++axis) {
        const auto axis_index = static_cast<std::size_t>(axis);
        if (range.begin[axis_index] < 0 || range.begin[axis_index] >= extent[axis_index]
            || range.end[axis_index] < 0 || range.end[axis_index] >= extent[axis_index]) {
            throw TopologyError(std::string(label) + " is outside its block cell extent");
        }
    }
}

} // namespace

HaloExchangePlan make_halo_exchange_plan(const ConnectivityPatch& connection,
                                         Extent3 receiver_cell_extent,
                                         Extent3 donor_cell_extent,
                                         int dimension)
{
    if (!connection.transform.valid(dimension)) {
        throw TopologyError("cannot build a halo plan from an invalid index transform");
    }
    if (connection.ghost_width < 0) {
        throw TopologyError("connectivity ghost width must be non-negative");
    }

    const int receiver_normal = axis_index(connection.receiver_face.axis);
    const int donor_normal = axis_index(connection.donor_face.axis);
    const int mapped_normal
        = std::abs(
              connection.transform.receiver_to_donor[static_cast<std::size_t>(receiver_normal)])
        - 1;
    if (mapped_normal != donor_normal) {
        throw TopologyError("connectivity transform does not map the receiver normal axis to the "
                            "donor normal axis");
    }
    validate_cell_range(connection.receiver_adjacent_cell_range,
                        receiver_cell_extent,
                        dimension,
                        "receiver cell-face range");
    validate_cell_range(connection.donor_adjacent_cell_range,
                        donor_cell_extent,
                        dimension,
                        "donor cell-face range");
    const int expected_receiver_face = connection.receiver_face.side == Side::Lower
        ? 0
        : receiver_cell_extent[static_cast<std::size_t>(receiver_normal)] - 1;
    const int expected_donor_face = connection.donor_face.side == Side::Lower
        ? 0
        : donor_cell_extent[static_cast<std::size_t>(donor_normal)] - 1;
    if (connection.receiver_adjacent_cell_range.begin[static_cast<std::size_t>(receiver_normal)]
            != expected_receiver_face
        || connection.receiver_adjacent_cell_range.end[static_cast<std::size_t>(receiver_normal)]
            != expected_receiver_face
        || connection.donor_adjacent_cell_range.begin[static_cast<std::size_t>(donor_normal)]
            != expected_donor_face
        || connection.donor_adjacent_cell_range.end[static_cast<std::size_t>(donor_normal)]
            != expected_donor_face) {
        throw TopologyError("connectivity cell ranges are not located on their declared faces");
    }
    if (connection.ghost_width > donor_cell_extent[static_cast<std::size_t>(donor_normal)]) {
        throw TopologyError("donor block is too thin for the requested ghost width");
    }

    const auto receiver_counts = connection.receiver_adjacent_cell_range.counts();
    const auto donor_counts = connection.donor_adjacent_cell_range.counts();
    for (int receiver_axis = 0; receiver_axis < dimension; ++receiver_axis) {
        const int donor_axis
            = std::abs(
                  connection.transform.receiver_to_donor[static_cast<std::size_t>(receiver_axis)])
            - 1;
        if (receiver_counts[static_cast<std::size_t>(receiver_axis)]
            != donor_counts[static_cast<std::size_t>(donor_axis)]) {
            throw TopologyError(
                "receiver and donor cell ranges have incompatible transformed extents");
        }
    }

    HaloExchangePlan plan {
        connection.name,
        connection.receiver_block,
        connection.donor_block,
        connection.donor_rank,
        receiver_cell_extent,
        donor_cell_extent,
        dimension,
        connection.transform,
        {},
    };
    plan.cell_pairs.reserve(connection.receiver_adjacent_cell_range.size()
                            * static_cast<std::size_t>(connection.ghost_width));

    for (int k = 0; k < receiver_counts.nk; ++k) {
        for (int j = 0; j < receiver_counts.nj; ++j) {
            for (int i = 0; i < receiver_counts.ni; ++i) {
                const Index3 receiver_ordinal {i, j, k};
                Index3 donor_ordinal;
                for (int receiver_axis = 0; receiver_axis < dimension; ++receiver_axis) {
                    const int donor_axis
                        = std::abs(connection.transform
                                       .receiver_to_donor[static_cast<std::size_t>(receiver_axis)])
                        - 1;
                    donor_ordinal[static_cast<std::size_t>(donor_axis)]
                        = receiver_ordinal[static_cast<std::size_t>(receiver_axis)];
                }

                const auto receiver_face_cell
                    = connection.receiver_adjacent_cell_range.at(receiver_ordinal);
                const auto donor_face_cell = connection.donor_adjacent_cell_range.at(donor_ordinal);
                for (int layer = 0; layer < connection.ghost_width; ++layer) {
                    auto receiver_ghost = receiver_face_cell;
                    receiver_ghost[static_cast<std::size_t>(receiver_normal)]
                        += outward_step(connection.receiver_face.side) * (layer + 1);
                    auto donor_interior = donor_face_cell;
                    donor_interior[static_cast<std::size_t>(donor_normal)]
                        += inward_step(connection.donor_face.side) * layer;
                    plan.cell_pairs.push_back({receiver_ghost, donor_interior});
                }
            }
        }
    }
    return plan;
}

} // namespace wcns
