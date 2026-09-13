#include <wcns/mesh/conservation_weights.hpp>

#include <wcns/mesh/linear_operators.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace wcns {
namespace {

std::vector<Real> solve_dense(std::vector<std::vector<Real>> matrix, std::vector<Real> rhs)
{
    const int count = static_cast<int>(rhs.size());
    for (int column = 0; column < count; ++column) {
        int pivot_row = column;
        for (int row = column + 1; row < count; ++row) {
            if (std::abs(matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)])
                > std::abs(matrix[static_cast<std::size_t>(pivot_row)]
                                 [static_cast<std::size_t>(column)])) {
                pivot_row = row;
            }
        }
        const Real pivot
            = matrix[static_cast<std::size_t>(pivot_row)][static_cast<std::size_t>(column)];
        if (!std::isfinite(pivot)
            || std::abs(pivot) < 1024.0 * std::numeric_limits<Real>::epsilon()) {
            throw ConservationWeightError(
                "conservation-weight system does not have a unique solution");
        }
        std::swap(matrix[static_cast<std::size_t>(column)],
                  matrix[static_cast<std::size_t>(pivot_row)]);
        std::swap(rhs[static_cast<std::size_t>(column)], rhs[static_cast<std::size_t>(pivot_row)]);
        for (int row = column + 1; row < count; ++row) {
            const Real factor
                = matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
                / matrix[static_cast<std::size_t>(column)][static_cast<std::size_t>(column)];
            for (int inner = column; inner < count; ++inner) {
                matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(inner)] -= factor
                    * matrix[static_cast<std::size_t>(column)][static_cast<std::size_t>(inner)];
            }
            rhs[static_cast<std::size_t>(row)] -= factor * rhs[static_cast<std::size_t>(column)];
        }
    }
    std::vector<Real> result(static_cast<std::size_t>(count), 0.0);
    for (int row = count - 1; row >= 0; --row) {
        Real value = rhs[static_cast<std::size_t>(row)];
        for (int column = row + 1; column < count; ++column) {
            value -= matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
                * result[static_cast<std::size_t>(column)];
        }
        result[static_cast<std::size_t>(row)]
            = value / matrix[static_cast<std::size_t>(row)][static_cast<std::size_t>(row)];
    }
    return result;
}

LineConservationWeights
build_uncached_line(const AlgorithmProfile& profile, int cell_count, bool periodic)
{
    LineConservationWeights result;
    result.profile = profile.kind();
    result.cell_count = cell_count;
    result.periodic = periodic;
    if (cell_count <= 0) {
        throw ConservationWeightError("conservation-weight line must contain cells");
    }
    if (periodic) {
        result.cell_weights.assign(static_cast<std::size_t>(cell_count), 1.0);
        return result;
    }

    const auto operators = LineOperators::build(profile, cell_count);
    const auto& rows = operators.derivative_rows();
    std::vector<std::vector<Real>> transpose(
        static_cast<std::size_t>(cell_count + 1),
        std::vector<Real>(static_cast<std::size_t>(cell_count), 0.0));
    for (int cell = 0; cell < cell_count; ++cell) {
        for (const auto [face, coefficient] : rows[static_cast<std::size_t>(cell)]) {
            transpose[static_cast<std::size_t>(face)][static_cast<std::size_t>(cell)] = coefficient;
        }
    }
    std::vector<Real> boundary(static_cast<std::size_t>(cell_count + 1), 0.0);
    boundary.front() = -1.0;
    boundary.back() = 1.0;

    std::vector<std::vector<Real>> normal(
        static_cast<std::size_t>(cell_count),
        std::vector<Real>(static_cast<std::size_t>(cell_count), 0.0));
    std::vector<Real> normal_rhs(static_cast<std::size_t>(cell_count), 0.0);
    for (int row = 0; row < cell_count; ++row) {
        for (int column = 0; column < cell_count; ++column) {
            for (int face = 0; face <= cell_count; ++face) {
                normal[static_cast<std::size_t>(row)][static_cast<std::size_t>(column)]
                    += transpose[static_cast<std::size_t>(face)][static_cast<std::size_t>(row)]
                    * transpose[static_cast<std::size_t>(face)][static_cast<std::size_t>(column)];
            }
        }
        for (int face = 0; face <= cell_count; ++face) {
            normal_rhs[static_cast<std::size_t>(row)]
                += transpose[static_cast<std::size_t>(face)][static_cast<std::size_t>(row)]
                * boundary[static_cast<std::size_t>(face)];
        }
    }
    result.cell_weights = solve_dense(std::move(normal), std::move(normal_rhs));
    for (const Real weight : result.cell_weights) {
        if (!std::isfinite(weight) || weight <= 0.0) {
            throw ConservationWeightError(
                "conservation-weight system produced a non-positive cell weight");
        }
    }
    for (int face = 0; face <= cell_count; ++face) {
        Real actual = 0.0;
        for (int cell = 0; cell < cell_count; ++cell) {
            actual += transpose[static_cast<std::size_t>(face)][static_cast<std::size_t>(cell)]
                * result.cell_weights[static_cast<std::size_t>(cell)];
        }
        result.maximum_residual = std::max(
            result.maximum_residual, std::abs(actual - boundary[static_cast<std::size_t>(face)]));
    }
    if (result.maximum_residual > 1.0e-10) {
        throw ConservationWeightError(
            "conservation-weight residual exceeds the initialization tolerance: count="
            + std::to_string(cell_count) + ", residual=" + std::to_string(result.maximum_residual));
    }
    return result;
}

Real tangential_weight(const std::array<LineConservationWeights, 3>& lines,
                       FaceLocation face,
                       Index3 index,
                       int dimension)
{
    Real result = 1.0;
    for (int axis = 0; axis < dimension; ++axis) {
        if (axis != static_cast<int>(face.axis)) {
            result *= lines[static_cast<std::size_t>(axis)].cell_weights[static_cast<std::size_t>(
                index[static_cast<std::size_t>(axis)])];
        }
    }
    return result;
}

Array3D<Real>& boundary_field(BlockConservationWeights& block, Axis axis)
{
    switch (axis) {
    case Axis::I: return block.physical_i_face;
    case Axis::J: return block.physical_j_face;
    case Axis::K: return block.physical_k_face;
    }
    throw ConservationWeightError("invalid boundary-weight face axis");
}

bool range_contains(const IndexRange3& range, Index3 index)
{
    for (int axis = 0; axis < 3; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        const int lower = std::min(range.begin[a], range.end[a]);
        const int upper = std::max(range.begin[a], range.end[a]);
        if (index[a] < lower || index[a] > upper) return false;
    }
    return true;
}

std::vector<const ConnectivityPatch*>
side_connections(const StructuredBlock& block, Axis axis, Side side)
{
    const auto cells = block.cell_extent();
    std::vector<const ConnectivityPatch*> result;
    for (const auto& connection : block.connectivities) {
        if (connection.receiver_face.axis == axis && connection.receiver_face.side == side) {
            result.push_back(&connection);
        }
    }
    if (result.empty()) return result;
    // A side that mixes physical boundary segments and connectivity segments
    // cannot share one tensor-product normal-line continuation. Treat its
    // normal direction as locally bounded; paired connection fluxes still
    // cancel through the independently checked tangential weights below.
    const int normal = side == Side::Lower ? 0 : cells[static_cast<std::size_t>(axis)];
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                Index3 face {i, j, k};
                face[static_cast<std::size_t>(axis)] = normal;
                bool connected = false;
                for (const auto* connection : result) {
                    if (range_contains(connection->shared_face_range.untyped(), face)) {
                        connected = true;
                        break;
                    }
                }
                if (!connected) {
                    return {};
                }
            }
        }
    }
    return result;
}

struct BoundaryTrace {
    bool cycle = false;
    int distance = 0;
};

BoundaryTrace trace_boundary(const StructuredMesh& mesh,
                             BlockId block_id,
                             Axis axis,
                             Side outward_side,
                             std::set<std::tuple<BlockId, int, int>> path)
{
    const auto key = std::tuple {block_id, static_cast<int>(axis), static_cast<int>(outward_side)};
    if (!path.insert(key).second) return {true, 0};
    const auto& block = mesh.block(block_id);
    const auto connections = side_connections(block, axis, outward_side);
    if (connections.empty()) return {false, 0};

    bool first = true;
    BoundaryTrace merged;
    for (const auto* connection : connections) {
        const int donor_axis_index
            = std::abs(connection->transform.receiver_to_donor[static_cast<std::size_t>(axis)]) - 1;
        if (donor_axis_index < 0 || donor_axis_index >= block.cell_dimension()) {
            throw ConservationWeightError(
                "connection transform maps a normal axis outside the mesh dimension");
        }
        const auto donor_axis = static_cast<Axis>(donor_axis_index);
        const Side next_side
            = connection->donor_face.side == Side::Lower ? Side::Upper : Side::Lower;
        auto branch = trace_boundary(mesh, connection->donor_block, donor_axis, next_side, path);
        if (!branch.cycle) {
            const int cells = mesh.block(connection->donor_block)
                                  .cell_extent()[static_cast<std::size_t>(donor_axis)];
            if (branch.distance > std::numeric_limits<int>::max() - cells) {
                throw ConservationWeightError("connected line length exceeds int range");
            }
            branch.distance += cells;
        }
        if (first) {
            merged = branch;
            first = false;
        } else if (merged.cycle != branch.cycle
                   || (!merged.cycle && merged.distance != branch.distance)) {
            throw ConservationWeightError(
                "split connection branches have incompatible normal line lengths");
        }
    }
    return merged;
}

} // namespace

LineConservationWeights
build_line_conservation_weights(const AlgorithmProfile& profile, int cell_count, bool periodic)
{
    using CacheKey = std::tuple<int, int, bool>;
    static std::map<CacheKey, LineConservationWeights> cache;
    const CacheKey key {static_cast<int>(profile.kind()), cell_count, periodic};
    const auto iterator = cache.find(key);
    if (iterator != cache.end()) {
        return iterator->second;
    }
    auto result = build_uncached_line(profile, cell_count, periodic);
    cache.emplace(key, result);
    return result;
}

GlobalConservationWeights GlobalConservationWeights::build(const StructuredMesh& mesh,
                                                           const AlgorithmProfile& profile)
{
    mesh.validate_connectivities(false);
    GlobalConservationWeights result;
    result.profile_ = profile.kind();
    std::unordered_map<BlockId, std::array<LineConservationWeights, 3>> line_weights;
    for (const auto& block : mesh.blocks()) {
        std::array<LineConservationWeights, 3> lines;
        for (int axis = 0; axis < block.cell_dimension(); ++axis) {
            const auto direction = static_cast<Axis>(axis);
            const auto lower = trace_boundary(mesh, block.id(), direction, Side::Lower, {});
            const auto upper = trace_boundary(mesh, block.id(), direction, Side::Upper, {});
            if (lower.cycle != upper.cycle) {
                throw ConservationWeightError("connected line reaches a cycle on only one side");
            }
            const int local_count = block.cell_extent()[static_cast<std::size_t>(axis)];
            if (lower.cycle) {
                lines[static_cast<std::size_t>(axis)]
                    = build_line_conservation_weights(profile, local_count, true);
            } else {
                if (lower.distance
                    > std::numeric_limits<int>::max() - local_count - upper.distance) {
                    throw ConservationWeightError("connected line length exceeds int range");
                }
                const int global_count = lower.distance + local_count + upper.distance;
                const auto global = build_line_conservation_weights(profile, global_count, false);
                auto& local = lines[static_cast<std::size_t>(axis)];
                local.profile = profile.kind();
                local.cell_count = local_count;
                local.cell_weights.assign(global.cell_weights.begin() + lower.distance,
                                          global.cell_weights.begin() + lower.distance
                                              + local_count);
                local.maximum_residual = global.maximum_residual;
            }
            result.maximum_line_residual_
                = std::max(result.maximum_line_residual_,
                           lines[static_cast<std::size_t>(axis)].maximum_residual);
        }
        if (block.cell_dimension() == 2) {
            lines[2] = build_line_conservation_weights(profile, 1, true);
        }
        line_weights.emplace(block.id(), lines);
        auto [iterator, inserted]
            = result.blocks_.emplace(block.id(), BlockConservationWeights(block.cell_extent()));
        static_cast<void>(inserted);
        auto& block_weights = iterator->second;
        const auto extent = block.cell_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    block_weights.cell(i, j, k) = lines[0].cell_weights[static_cast<std::size_t>(i)]
                        * lines[1].cell_weights[static_cast<std::size_t>(j)]
                        * lines[2].cell_weights[static_cast<std::size_t>(k)];
                }
            }
        }
        for (const auto& patch : block.boundaries) {
            auto& faces = boundary_field(block_weights, patch.face.axis);
            const auto counts = patch.boundary_face_range.counts();
            for (int k = 0; k < counts.nk; ++k) {
                for (int j = 0; j < counts.nj; ++j) {
                    for (int i = 0; i < counts.ni; ++i) {
                        const auto face = patch.boundary_face_range.at({i, j, k});
                        faces(face.i, face.j, face.k)
                            = tangential_weight(lines, patch.face, face, block.cell_dimension());
                    }
                }
            }
        }
    }

    for (const auto& block : mesh.blocks()) {
        for (const auto& connection : block.connectivities) {
            const auto& receiver_lines = line_weights.at(connection.receiver_block);
            const auto& donor_lines = line_weights.at(connection.donor_block);
            const auto counts = connection.shared_face_range.counts();
            for (int k = 0; k < counts.nk; ++k) {
                for (int j = 0; j < counts.nj; ++j) {
                    for (int i = 0; i < counts.ni; ++i) {
                        const Index3 receiver_ordinal {i, j, k};
                        Index3 donor_ordinal;
                        for (int receiver_axis = 0; receiver_axis < block.cell_dimension();
                             ++receiver_axis) {
                            const int donor_axis
                                = std::abs(connection.transform
                                               .receiver_to_donor[static_cast<std::size_t>(
                                                   receiver_axis)])
                                - 1;
                            donor_ordinal[static_cast<std::size_t>(donor_axis)]
                                = receiver_ordinal[static_cast<std::size_t>(receiver_axis)];
                        }
                        const auto receiver_face
                            = connection.shared_face_range.at(receiver_ordinal);
                        auto donor_face
                            = connection.donor_adjacent_cell_range.at(donor_ordinal);
                        const auto donor_axis
                            = static_cast<std::size_t>(connection.donor_face.axis);
                        donor_face[donor_axis]
                            = connection.donor_face.side == Side::Lower
                            ? 0
                            : mesh.block(connection.donor_block).cell_extent()[donor_axis];
                        const Real receiver_weight = tangential_weight(receiver_lines,
                                                                       connection.receiver_face,
                                                                       receiver_face,
                                                                       block.cell_dimension());
                        const Real donor_weight = tangential_weight(
                            donor_lines, connection.donor_face, donor_face, block.cell_dimension());
                        result.maximum_shared_face_mismatch_
                            = std::max(result.maximum_shared_face_mismatch_,
                                       std::abs(receiver_weight - donor_weight));
                    }
                }
            }
        }
    }
    if (result.maximum_shared_face_mismatch_ > 1.0e-10) {
        throw ConservationWeightError(
            "global conservation weights do not cancel on a shared connection face");
    }
    return result;
}

const BlockConservationWeights& GlobalConservationWeights::block(BlockId id) const
{
    const auto iterator = blocks_.find(id);
    if (iterator == blocks_.end()) {
        throw std::out_of_range("conservation weights do not contain the requested block");
    }
    return iterator->second;
}

} // namespace wcns
