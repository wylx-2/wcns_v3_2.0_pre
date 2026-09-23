#include <wcns/mesh/structured_mesh.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace wcns {
namespace {

std::array<Real, 3> coordinate(const StructuredBlock& block, Index3 index)
{
    return {
        block.coordinates.x(index.i, index.j, index.k),
        block.coordinates.y(index.i, index.j, index.k),
        block.coordinates.z(index.i, index.j, index.k),
    };
}

template <class LeftRange, class RightRange>
bool same_undirected_range(const LeftRange& lhs, const RightRange& rhs)
{
    return (lhs.begin == rhs.begin && lhs.end == rhs.end)
        || (lhs.begin == rhs.end && lhs.end == rhs.begin);
}

bool is_reciprocal(const ConnectivityPatch& connection,
                   const ConnectivityPatch& candidate,
                   int dimension)
{
    return candidate.receiver_block == connection.donor_block
        && candidate.donor_block == connection.receiver_block
        && candidate.receiver_face == connection.donor_face
        && candidate.donor_face == connection.receiver_face
        && same_undirected_range(candidate.receiver_vertex_range, connection.donor_vertex_range)
        && same_undirected_range(candidate.donor_vertex_range, connection.receiver_vertex_range)
        && same_undirected_range(candidate.receiver_adjacent_cell_range,
                                 connection.donor_adjacent_cell_range)
        && same_undirected_range(candidate.donor_adjacent_cell_range,
                                 connection.receiver_adjacent_cell_range)
        && candidate.transform == connection.transform.inverse(dimension)
        && candidate.periodic == connection.periodic.inverse();
}

bool coordinates_match(Real receiver, Real donor, bool cgns_periodic_transform)
{
    const Real scale = std::max({Real {1}, std::abs(receiver), std::abs(donor)});
    const Real transform_tolerance = cgns_periodic_transform
        ? Real {8} * std::numeric_limits<float>::epsilon() * scale
        : Real {0};
    return std::abs(receiver - donor)
        <= std::max(Real {256} * std::numeric_limits<Real>::epsilon() * scale, transform_tolerance);
}

bool ranges_overlap(const IndexRange3& lhs, const IndexRange3& rhs)
{
    for (int axis = 0; axis < 3; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        const int lhs_lower = std::min(lhs.begin[a], lhs.end[a]);
        const int lhs_upper = std::max(lhs.begin[a], lhs.end[a]);
        const int rhs_lower = std::min(rhs.begin[a], rhs.end[a]);
        const int rhs_upper = std::max(rhs.begin[a], rhs.end[a]);
        if (lhs_upper < rhs_lower || rhs_upper < lhs_lower) return false;
    }
    return true;
}

void validate_interface_coordinates(const StructuredBlock& receiver,
                                    const StructuredBlock& donor,
                                    const ConnectivityPatch& connection)
{
    const auto counts = connection.receiver_vertex_range.counts();
    const bool cgns_periodic_transform = !(connection.periodic == PeriodicTransform {});
    for (int k = 0; k < counts.nk; ++k) {
        for (int j = 0; j < counts.nj; ++j) {
            for (int i = 0; i < counts.ni; ++i) {
                const auto receiver_index = connection.receiver_vertex_range.at({i, j, k});
                const auto donor_index
                    = connection.transform.map(receiver_index,
                                               connection.receiver_vertex_range.begin,
                                               connection.donor_vertex_range.begin,
                                               receiver.cell_dimension());
                const auto expected_donor
                    = connection.periodic.apply_point(coordinate(receiver, receiver_index));
                const auto actual_donor = coordinate(donor, donor_index);
                if (!coordinates_match(expected_donor[0], actual_donor[0], cgns_periodic_transform)
                    || !coordinates_match(
                        expected_donor[1], actual_donor[1], cgns_periodic_transform)
                    || !coordinates_match(
                        expected_donor[2], actual_donor[2], cgns_periodic_transform)) {
                    throw TopologyError("connectivity " + connection.name
                                        + " maps vertices with different physical coordinates");
                }
            }
        }
    }
}

} // namespace

StructuredMesh::StructuredMesh(std::vector<StructuredBlock> blocks)
    : blocks_(std::move(blocks))
{
    block_index_.reserve(blocks_.size());
    for (std::size_t index = 0; index < blocks_.size(); ++index) {
        const auto [iterator, inserted] = block_index_.emplace(blocks_[index].id(), index);
        static_cast<void>(iterator);
        if (!inserted) {
            throw TopologyError("duplicate structured block id "
                                + std::to_string(blocks_[index].id()));
        }
    }
}

bool StructuredMesh::contains(BlockId id) const noexcept
{
    return block_index_.find(id) != block_index_.end();
}

StructuredBlock& StructuredMesh::block(BlockId id)
{
    const auto iterator = block_index_.find(id);
    if (iterator == block_index_.end()) {
        throw std::out_of_range("structured block id is not present in the mesh");
    }
    return blocks_[iterator->second];
}

const StructuredBlock& StructuredMesh::block(BlockId id) const
{
    const auto iterator = block_index_.find(id);
    if (iterator == block_index_.end()) {
        throw std::out_of_range("structured block id is not present in the mesh");
    }
    return blocks_[iterator->second];
}

void StructuredMesh::validate_connectivities(bool validate_coordinates) const
{
    for (const auto& receiver : blocks_) {
        for (const auto& connection : receiver.connectivities) {
            if (connection.receiver_block != receiver.id()) {
                throw TopologyError("connectivity " + connection.name
                                    + " has a receiver id inconsistent with its owning block");
            }
            if (!contains(connection.donor_block)) {
                throw TopologyError("connectivity " + connection.name
                                    + " references an unknown donor block");
            }
            const auto& donor = block(connection.donor_block);
            if (receiver.cell_dimension() != donor.cell_dimension()) {
                throw TopologyError("connectivity " + connection.name
                                    + " joins blocks with different cell dimensions");
            }
            const int dimension = receiver.cell_dimension();
            if (connection.donor_rank != donor.owner_rank()) {
                throw TopologyError("connectivity " + connection.name
                                    + " has a donor rank inconsistent with its donor block");
            }
            if (connection.ghost_width != receiver.ghost_width()) {
                throw TopologyError("connectivity " + connection.name
                                    + " has a ghost width inconsistent with its receiver block");
            }
            if (!connection.transform.valid(dimension)) {
                throw TopologyError("connectivity " + connection.name
                                    + " has an invalid index transform");
            }
            if (!connection.periodic.valid(dimension)) {
                throw TopologyError("connectivity " + connection.name
                                    + " has a non-orthogonal or improper periodic transform");
            }
            if (connection.transform.map(connection.receiver_vertex_range.end,
                                         connection.receiver_vertex_range.begin,
                                         connection.donor_vertex_range.begin,
                                         dimension)
                != connection.donor_vertex_range.end) {
                throw TopologyError(
                    "connectivity " + connection.name
                    + " transform does not map the receiver range onto the donor range");
            }
            const auto reciprocal
                = std::find_if(donor.connectivities.begin(),
                               donor.connectivities.end(),
                               [&](const ConnectivityPatch& candidate) {
                                   return is_reciprocal(connection, candidate, dimension);
                               });
            if (reciprocal == donor.connectivities.end()) {
                throw TopologyError("connectivity " + connection.name
                                    + " has no reciprocal donor record");
            }
            if (validate_coordinates) {
                validate_interface_coordinates(receiver, donor, connection);
            }
        }
    }
}

bool connection_side_is_fully_covered(const StructuredBlock& block, Axis axis, Side side)
{
    const auto cells = block.cell_extent();
    const auto normal = static_cast<std::size_t>(axis);
    std::size_t expected = 1;
    for (int logical = 0; logical < block.cell_dimension(); ++logical) {
        if (static_cast<std::size_t>(logical) != normal) {
            expected *= static_cast<std::size_t>(cells[static_cast<std::size_t>(logical)]);
        }
    }

    std::size_t covered = 0;
    std::vector<IndexRange3> ranges;
    for (const auto& connection : block.connectivities) {
        if (connection.receiver_face.axis != axis || connection.receiver_face.side != side) {
            continue;
        }
        const auto range = connection.shared_face_range.untyped();
        for (const auto& existing : ranges) {
            if (ranges_overlap(existing, range)) return false;
        }
        covered += range.size();
        ranges.push_back(range);
    }
    return !ranges.empty() && covered == expected;
}

} // namespace wcns
