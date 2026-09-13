#include "test_support.hpp"

#include <wcns/mesh/halo_exchange.hpp>
#include <wcns/mesh/structured_mesh.hpp>
#include <wcns/parallel/block_distribution.hpp>
#include <wcns/parallel/distributed_topology.hpp>

#include <stdexcept>
#include <utility>
#include <vector>

namespace {

wcns::ConnectivityPatch receiver_connection()
{
    return {
        "left-to-right",
        0,
        1,
        1,
        {wcns::Axis::I, wcns::Side::Upper},
        {wcns::Axis::J, wcns::Side::Lower},
        {{4, 0, 0}, {4, 3, 0}},
        {{3, 0, 0}, {0, 0, 0}},
        {{3, 0, 0}, {3, 2, 0}},
        {{2, 0, 0}, {0, 0, 0}},
        {{4, 0, 0}, {4, 2, 0}},
        {{{2, -1, 3}}},
        3,
    };
}

wcns::ConnectivityPatch donor_connection()
{
    return {
        "right-to-left",
        1,
        0,
        0,
        {wcns::Axis::J, wcns::Side::Lower},
        {wcns::Axis::I, wcns::Side::Upper},
        {{3, 0, 0}, {0, 0, 0}},
        {{4, 0, 0}, {4, 3, 0}},
        {{2, 0, 0}, {0, 0, 0}},
        {{3, 0, 0}, {3, 2, 0}},
        {{2, 0, 0}, {0, 0, 0}},
        {{{-2, 1, 3}}},
        3,
    };
}

} // namespace

// 验收多块连接变换、互反记录、坐标一致性和 halo 映射。
void test_topology()
{
    using namespace wcns;

    const IndexTransform transform {{{2, -1, 3}}};
    WCNS_REQUIRE(transform.valid(2));
    WCNS_REQUIRE(transform.map({4, 3, 0}, {4, 0, 0}, {3, 0, 0}, 2) == (Index3 {0, 0, 0}));
    WCNS_REQUIRE(transform.inverse(2) == (IndexTransform {{{-2, 1, 3}}}));
    WCNS_REQUIRE((!IndexTransform {{{2, 2, 3}}}.valid(2)));
    WCNS_REQUIRE((!IndexTransform {{{1, 2, -3}}}.valid(2)));
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        (IndexTransform {{{2, 2, 3}}}.map({0, 0, 0}, {0, 0, 0}, {0, 0, 0}, 2)));

    PeriodicTransform periodic {
        {{{{0.0, -1.0, 0.0}}, {{1.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}}},
        {{3.0, -2.0, 0.0}},
    };
    WCNS_REQUIRE(periodic.valid(2));
    const auto mapped_point = periodic.apply_point({{2.0, 1.0, 0.0}});
    WCNS_REQUIRE(mapped_point == (std::array<Real, 3> {{2.0, 0.0, 0.0}}));
    const auto recovered_point = periodic.inverse().apply_point(mapped_point);
    WCNS_REQUIRE(recovered_point == (std::array<Real, 3> {{2.0, 1.0, 0.0}}));
    auto invalid_periodic = periodic;
    invalid_periodic.rotation[0][0] = 0.5;
    WCNS_REQUIRE(!invalid_periodic.valid(2));

    StructuredBlock left(0, "left", 0, 2, 2, {5, 4, 1}, 3);
    StructuredBlock right(1, "right", 1, 2, 2, {4, 5, 1}, 3);
    left.connectivities.push_back(receiver_connection());
    right.connectivities.push_back(donor_connection());

    std::vector<StructuredBlock> blocks;
    blocks.push_back(std::move(left));
    blocks.push_back(std::move(right));
    StructuredMesh mesh(std::move(blocks));
    WCNS_REQUIRE(mesh.block_count() == 2);
    WCNS_REQUIRE(mesh.contains(0));
    WCNS_REQUIRE(mesh.block(1).name() == "right");
    mesh.validate_connectivities();

    const auto halo = make_halo_exchange_plan(mesh.block(0).connectivities.front(),
                                              mesh.block(0).cell_extent(),
                                              mesh.block(1).cell_extent(),
                                              2);
    WCNS_REQUIRE(halo.receiver_block == 0);
    WCNS_REQUIRE(halo.donor_block == 1);
    WCNS_REQUIRE(halo.donor_rank == 1);
    WCNS_REQUIRE(halo.cell_pairs.size() == 9);
    WCNS_REQUIRE(halo.cell_pairs[0] == (HaloCellPair {{4, 0, 0}, {2, 0, 0}}));
    WCNS_REQUIRE(halo.cell_pairs[1] == (HaloCellPair {{5, 0, 0}, {2, 1, 0}}));
    WCNS_REQUIRE(halo.cell_pairs[3] == (HaloCellPair {{4, 1, 0}, {1, 0, 0}}));
    WCNS_REQUIRE_THROWS(
        TopologyError,
        make_halo_exchange_plan(
            mesh.block(0).connectivities.front(), mesh.block(0).cell_extent(), {3, 2, 1}, 2));

    mesh.block(0).connectivities.front().donor_rank = 7;
    WCNS_REQUIRE_THROWS(TopologyError, mesh.validate_connectivities());
    mesh.block(0).connectivities.front().donor_rank = 1;
    mesh.block(0).connectivities.front().ghost_width = 2;
    WCNS_REQUIRE_THROWS(TopologyError, mesh.validate_connectivities());
    mesh.block(0).connectivities.front().ghost_width = 3;
    mesh.block(1).coordinates.x(3, 0, 0) = 1.0;
    WCNS_REQUIRE_THROWS(TopologyError, mesh.validate_connectivities());

    StructuredBlock incomplete_left(0, "left", 0, 2, 2, {5, 4, 1}, 3);
    StructuredBlock incomplete_right(1, "right", 1, 2, 2, {4, 5, 1}, 3);
    incomplete_left.connectivities.push_back(receiver_connection());
    std::vector<StructuredBlock> incomplete_blocks;
    incomplete_blocks.push_back(std::move(incomplete_left));
    incomplete_blocks.push_back(std::move(incomplete_right));
    StructuredMesh incomplete(std::move(incomplete_blocks));
    WCNS_REQUIRE_THROWS(TopologyError, incomplete.validate_connectivities());
    WCNS_REQUIRE_THROWS(std::out_of_range, incomplete.block(7));

    auto first_cut = receiver_connection();
    first_cut.donor_block = 0;
    first_cut.donor_rank = 0;
    auto second_cut = donor_connection();
    second_cut.receiver_block = 0;
    second_cut.donor_block = 0;
    second_cut.donor_rank = 0;
    StructuredBlock o_grid(0, "o-grid", 0, 2, 2, {5, 4, 1}, 3);
    o_grid.connectivities = {first_cut, second_cut};
    StructuredMesh self_mesh(std::vector<StructuredBlock> {std::move(o_grid)});
    self_mesh.validate_connectivities(false);
    const auto self_distribution
        = BlockDistribution::balanced({{0, self_mesh.block(0).cell_extent().size()}}, 1);
    const auto self_topology
        = DistributedTopology::build(self_mesh, self_distribution, false);
    WCNS_REQUIRE(self_topology.exchanges().size() == 2);
    WCNS_REQUIRE(self_topology.local_copies(0).size() == 2);
}
