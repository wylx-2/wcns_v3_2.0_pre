#include <wcns/parallel/distributed_topology.hpp>

#include <algorithm>
#include <limits>
#include <tuple>

namespace wcns {
namespace {

template <class LeftRange, class RightRange>
bool same_undirected_range(const LeftRange& lhs, const RightRange& rhs)
{
    return (lhs.begin == rhs.begin && lhs.end == rhs.end)
        || (lhs.begin == rhs.end && lhs.end == rhs.begin);
}

bool reciprocal(const ConnectivityPatch& lhs, const ConnectivityPatch& rhs, int dimension)
{
    return rhs.receiver_block == lhs.donor_block && rhs.donor_block == lhs.receiver_block
        && rhs.receiver_face == lhs.donor_face && rhs.donor_face == lhs.receiver_face
        && same_undirected_range(rhs.receiver_vertex_range, lhs.donor_vertex_range)
        && same_undirected_range(rhs.donor_vertex_range, lhs.receiver_vertex_range)
        && rhs.transform == lhs.transform.inverse(dimension)
        && rhs.periodic == lhs.periodic.inverse();
}

auto connection_key(const ConnectivityPatch& connection)
{
    return std::tuple {
        connection.receiver_block,
        connection.donor_block,
        static_cast<int>(connection.receiver_face.axis),
        static_cast<int>(connection.receiver_face.side),
        connection.receiver_vertex_range.begin.i,
        connection.receiver_vertex_range.begin.j,
        connection.receiver_vertex_range.begin.k,
        connection.receiver_vertex_range.end.i,
        connection.receiver_vertex_range.end.j,
        connection.receiver_vertex_range.end.k,
    };
}

} // namespace

int DirectedExchange::message_tag(int tag_base) const
{
    if (tag_base < 0 || connection < 0) {
        throw TopologyError("message tags require non-negative base and connection id");
    }
    const int direction = halo.receiver_block < halo.donor_block ? 0 : 1;
    const auto tag
        = static_cast<long long>(tag_base) + 2LL * static_cast<long long>(connection) + direction;
    if (tag > std::numeric_limits<int>::max()) {
        throw TopologyError("message tag exceeds int range");
    }
    return static_cast<int>(tag);
}

DistributedTopology DistributedTopology::build(const StructuredMesh& mesh,
                                               const BlockDistribution& distribution,
                                               bool validate_coordinates)
{
    if (mesh.block_count() != distribution.assignments().size()) {
        throw TopologyError("mesh and distribution contain different block counts");
    }
    mesh.validate_connectivities(validate_coordinates);

    std::vector<const ConnectivityPatch*> canonical;
    for (const auto& block : mesh.blocks()) {
        for (const auto& connection : block.connectivities) {
            if (connection.receiver_block < connection.donor_block) {
                canonical.push_back(&connection);
            } else if (connection.receiver_block == connection.donor_block) {
                const auto reciprocal_iterator = std::find_if(
                    block.connectivities.begin(),
                    block.connectivities.end(),
                    [&](const ConnectivityPatch& candidate) {
                        return &candidate != &connection
                            && reciprocal(connection, candidate, block.cell_dimension());
                    });
                if (reciprocal_iterator != block.connectivities.end()) {
                    if (connection_key(connection) < connection_key(*reciprocal_iterator)) {
                        canonical.push_back(&connection);
                    }
                } else if (reciprocal(connection, connection, block.cell_dimension())) {
                    canonical.push_back(&connection);
                }
            }
        }
    }
    std::sort(canonical.begin(), canonical.end(), [](const auto* lhs, const auto* rhs) {
        return connection_key(*lhs) < connection_key(*rhs);
    });

    DistributedTopology result;
    result.exchanges_.reserve(canonical.size() * 2);
    for (std::size_t index = 0; index < canonical.size(); ++index) {
        if (index > static_cast<std::size_t>(std::numeric_limits<ConnectionId>::max())) {
            throw TopologyError("connection count exceeds ConnectionId range");
        }
        const auto& forward = *canonical[index];
        const auto& donor_block = mesh.block(forward.donor_block);
        const auto reverse_iterator = std::find_if(
            donor_block.connectivities.begin(),
            donor_block.connectivities.end(),
            [&](const ConnectivityPatch& candidate) {
                return reciprocal(
                    forward, candidate, mesh.block(forward.receiver_block).cell_dimension());
            });
        if (reverse_iterator == donor_block.connectivities.end()) {
            throw TopologyError("canonical connectivity has no reciprocal record");
        }

        const auto connection_id = static_cast<ConnectionId>(index);
        std::vector<const ConnectivityPatch*> directed_connections {&forward};
        if (&forward != &*reverse_iterator) directed_connections.push_back(&*reverse_iterator);
        for (const auto* directed : directed_connections) {
            auto connection = *directed;
            connection.id = connection_id;
            connection.donor_rank = distribution.owner(connection.donor_block);
            const auto receiver_rank = distribution.owner(connection.receiver_block);
            const auto donor_rank = distribution.owner(connection.donor_block);
            result.exchanges_.push_back({
                connection_id,
                receiver_rank,
                donor_rank,
                make_halo_exchange_plan(connection,
                                        mesh.block(connection.receiver_block).cell_extent(),
                                        mesh.block(connection.donor_block).cell_extent(),
                                        mesh.block(connection.receiver_block).cell_dimension()),
            });
        }
    }
    std::sort(result.exchanges_.begin(),
              result.exchanges_.end(),
              [](const DirectedExchange& lhs, const DirectedExchange& rhs) {
                  return std::tuple {lhs.connection,
                                     lhs.halo.receiver_block,
                                     lhs.halo.connectivity_name}
                      < std::tuple {rhs.connection,
                                     rhs.halo.receiver_block,
                                     rhs.halo.connectivity_name};
              });
    return result;
}

std::vector<const DirectedExchange*> DistributedTopology::receives(RankId rank) const
{
    std::vector<const DirectedExchange*> result;
    for (const auto& exchange : exchanges_) {
        if (exchange.receiver_rank == rank && exchange.donor_rank != rank) {
            result.push_back(&exchange);
        }
    }
    return result;
}

std::vector<const DirectedExchange*> DistributedTopology::sends(RankId rank) const
{
    std::vector<const DirectedExchange*> result;
    for (const auto& exchange : exchanges_) {
        if (exchange.donor_rank == rank && exchange.receiver_rank != rank) {
            result.push_back(&exchange);
        }
    }
    return result;
}

std::vector<const DirectedExchange*> DistributedTopology::local_copies(RankId rank) const
{
    std::vector<const DirectedExchange*> result;
    for (const auto& exchange : exchanges_) {
        if (exchange.receiver_rank == rank && exchange.donor_rank == rank) {
            result.push_back(&exchange);
        }
    }
    return result;
}

} // namespace wcns
