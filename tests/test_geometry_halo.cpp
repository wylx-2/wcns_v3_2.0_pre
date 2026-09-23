#include "test_support.hpp"

#include <wcns/mesh/geometry_halo.hpp>
#include <wcns/mesh/conservation_weights.hpp>
#include <wcns/solver/inviscid_flux.hpp>

#include <set>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

wcns::StructuredMesh make_geometry_mesh()
{
    using namespace wcns;
    StructuredBlock left(0, "geometry-left", 0, 2, 2, {7, 7, 1}, 3);
    StructuredBlock right(1, "geometry-right", 1, 2, 2, {7, 7, 1}, 3);
    for (int j = 0; j < 7; ++j) {
        for (int i = 0; i < 7; ++i) {
            left.coordinates.x(i, j, 0) = static_cast<Real>(i);
            left.coordinates.y(i, j, 0) = static_cast<Real>(j);
            left.coordinates.z(i, j, 0) = 0.0;
            right.coordinates.x(i, j, 0) = static_cast<Real>(j + 6);
            right.coordinates.y(i, j, 0) = static_cast<Real>(6 - i);
            right.coordinates.z(i, j, 0) = 0.0;
        }
    }
    left.connectivities.push_back({
        "left-to-right",
        0,
        1,
        1,
        {Axis::I, Side::Upper},
        {Axis::J, Side::Lower},
        {{6, 0, 0}, {6, 6, 0}},
        {{6, 0, 0}, {0, 0, 0}},
        {{5, 0, 0}, {5, 5, 0}},
        {{5, 0, 0}, {0, 0, 0}},
        {{6, 0, 0}, {6, 5, 0}},
        {{{2, -1, 3}}},
        3,
    });
    right.connectivities.push_back({
        "right-to-left",
        1,
        0,
        0,
        {Axis::J, Side::Lower},
        {Axis::I, Side::Upper},
        {{6, 0, 0}, {0, 0, 0}},
        {{6, 0, 0}, {6, 6, 0}},
        {{5, 0, 0}, {0, 0, 0}},
        {{5, 0, 0}, {5, 5, 0}},
        {{5, 0, 0}, {0, 0, 0}},
        {{{-2, 1, 3}}},
        3,
    });
    std::vector<StructuredBlock> blocks;
    blocks.push_back(std::move(left));
    blocks.push_back(std::move(right));
    return StructuredMesh(std::move(blocks));
}

wcns::StructuredMesh make_periodic_geometry_mesh()
{
    using namespace wcns;
    StructuredBlock left(0, "periodic-left", 0, 2, 2, {7, 7, 1}, 3);
    StructuredBlock right(1, "periodic-right", 1, 2, 2, {7, 7, 1}, 3);
    const PeriodicTransform periodic {
        {{{{0.0, -1.0, 0.0}}, {{1.0, 0.0, 0.0}}, {{0.0, 0.0, 1.0}}}},
        {{10.0, -3.0, 0.0}},
    };
    for (int j = 0; j < 7; ++j) {
        for (int i = 0; i < 7; ++i) {
            left.coordinates.x(i, j, 0) = static_cast<Real>(i);
            left.coordinates.y(i, j, 0) = static_cast<Real>(j);
            left.coordinates.z(i, j, 0) = 0.0;
            const auto mapped
                = periodic.apply_point({{static_cast<Real>(i + 6), static_cast<Real>(j), 0.0}});
            right.coordinates.x(i, j, 0) = mapped[0];
            right.coordinates.y(i, j, 0) = mapped[1];
            right.coordinates.z(i, j, 0) = mapped[2];
        }
    }
    left.connectivities.push_back({
        "periodic-forward",
        0,
        1,
        1,
        {Axis::I, Side::Upper},
        {Axis::I, Side::Lower},
        {{6, 0, 0}, {6, 6, 0}},
        {{0, 0, 0}, {0, 6, 0}},
        {{5, 0, 0}, {5, 5, 0}},
        {{0, 0, 0}, {0, 5, 0}},
        {{6, 0, 0}, {6, 5, 0}},
        {{{1, 2, 3}}},
        3,
        invalid_connection_id,
        periodic,
    });
    right.connectivities.push_back({
        "periodic-reverse",
        1,
        0,
        0,
        {Axis::I, Side::Lower},
        {Axis::I, Side::Upper},
        {{0, 0, 0}, {0, 6, 0}},
        {{6, 0, 0}, {6, 6, 0}},
        {{0, 0, 0}, {0, 5, 0}},
        {{5, 0, 0}, {5, 5, 0}},
        {{0, 0, 0}, {0, 5, 0}},
        {{{1, 2, 3}}},
        3,
        invalid_connection_id,
        periodic.inverse(),
    });
    std::vector<StructuredBlock> blocks;
    blocks.push_back(std::move(left));
    blocks.push_back(std::move(right));
    return StructuredMesh(std::move(blocks));
}

wcns::StructuredMesh make_partial_c_grid_mesh()
{
    using namespace wcns;
    StructuredBlock block(0, "partial-c-grid", 0, 2, 2, {17, 9, 1}, 3);
    for (int j = 0; j < 9; ++j) {
        for (int i = 0; i < 17; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<Real>(i);
            block.coordinates.y(i, j, 0) = static_cast<Real>(j);
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    block.connectivities.push_back({"wake-left-to-right",
                                    0,
                                    0,
                                    0,
                                    {Axis::J, Side::Lower},
                                    {Axis::J, Side::Lower},
                                    {{0, 0, 0}, {4, 0, 0}},
                                    {{16, 0, 0}, {12, 0, 0}},
                                    {{0, 0, 0}, {3, 0, 0}},
                                    {{15, 0, 0}, {12, 0, 0}},
                                    {{0, 0, 0}, {3, 0, 0}},
                                    {{{-1, 2, 3}}},
                                    3});
    block.connectivities.push_back({"wake-right-to-left",
                                    0,
                                    0,
                                    0,
                                    {Axis::J, Side::Lower},
                                    {Axis::J, Side::Lower},
                                    {{16, 0, 0}, {12, 0, 0}},
                                    {{0, 0, 0}, {4, 0, 0}},
                                    {{15, 0, 0}, {12, 0, 0}},
                                    {{0, 0, 0}, {3, 0, 0}},
                                    {{15, 0, 0}, {12, 0, 0}},
                                    {{{-1, 2, 3}}},
                                    3});
    std::vector<StructuredBlock> blocks;
    blocks.push_back(std::move(block));
    return StructuredMesh(std::move(blocks));
}

} // namespace

// 验收分阶段几何消息的种类、层宽、donor 路径、唯一标签和共享面所有者。
void test_geometry_halo_plan()
{
    using namespace wcns;
    const auto mesh = make_geometry_mesh();
    const auto ph = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
    const auto ph_plan = GeometryHaloPlan::build(mesh, ph);
    WCNS_REQUIRE(ph_plan.exchanges().size() == 12);
    std::set<int> tags;
    for (const auto& exchange : ph_plan.exchanges()) {
        WCNS_REQUIRE(exchange.shared_face_owner == 0);
        WCNS_REQUIRE(exchange.donor_path.size() == 1);
        WCNS_REQUIRE(tags.insert(exchange.message_tag()).second);
        if (exchange.kind == GeometryMessageKind::GeometryVertex) {
            WCNS_REQUIRE(exchange.halo_width == 2);
            if (exchange.receiver_block == 0) {
                WCNS_REQUIRE(exchange.index_transform == (IndexTransform {{{2, -1, 3}}}));
            }
        } else {
            WCNS_REQUIRE(exchange.halo_width == 3);
        }
    }

    const auto scmm = ProfileFactory::create(AlgorithmProfileKind::Scmm6Wcns);
    const auto scmm_plan = GeometryHaloPlan::build(mesh, scmm);
    for (const auto& exchange : scmm_plan.exchanges()) {
        if (exchange.kind == GeometryMessageKind::GeometryOperand) {
            WCNS_REQUIRE(exchange.halo_width == 5);
        }
    }
}

// 验收较小块号发布唯一共享面矢量后两侧面积和外法向严格对应。
void test_shared_metric_synchronization()
{
    using namespace wcns;
    auto mesh = make_geometry_mesh();
    const auto profile = ProfileFactory::create(AlgorithmProfileKind::Scmm6Wcns);
    std::unordered_map<BlockId, MetricField> metrics;
    for (const BlockId id : {BlockId {0}, BlockId {1}}) {
        auto& block = mesh.block(id);
        auto result = initialize_metric_field(block, profile);
        metrics.emplace(block.id(), std::move(result.metric));
    }
    SharedMetricSynchronizer::synchronize(mesh, metrics);
    const auto& left = metrics.at(0).i_faces();
    const auto& right = metrics.at(1).j_faces();
    for (int j = 0; j < 6; ++j) {
        const int donor_i = 5 - j;
        WCNS_REQUIRE(left.x(6, j, 0) == right.x(donor_i, 0, 0));
        WCNS_REQUIRE(left.y(6, j, 0) == right.y(donor_i, 0, 0));
        WCNS_REQUIRE(left.area(6, j, 0) == right.area(donor_i, 0, 0));
        WCNS_REQUIRE_NEAR(left.x(6, j, 0), 1.0, 2.0e-12);
    }
}

// 验收旋转周期连接按 Q/d 校验坐标并按 Q 旋转唯一共享面积矢量。
void test_periodic_shared_metric_synchronization()
{
    using namespace wcns;
    auto mesh = make_periodic_geometry_mesh();
    mesh.validate_connectivities();
    const auto profile = ProfileFactory::create(AlgorithmProfileKind::Scmm6Wcns);
    const auto plan = GeometryHaloPlan::build(mesh, profile);
    const auto rotated = plan.exchanges().front().periodic.apply_vector({{1.0, 0.0, 0.0}});
    WCNS_REQUIRE_NEAR(rotated[0], 0.0, 0.0);
    WCNS_REQUIRE_NEAR(rotated[1], 1.0, 0.0);

    std::unordered_map<BlockId, MetricField> metrics;
    for (const BlockId id : {BlockId {0}, BlockId {1}}) {
        auto result = initialize_metric_field(mesh.block(id), profile);
        metrics.emplace(id, std::move(result.metric));
    }
    SharedMetricSynchronizer::synchronize(mesh, metrics);
    const auto& owner = metrics.at(0).i_faces();
    const auto& donor = metrics.at(1).i_faces();
    for (int j = 0; j < 6; ++j) {
        WCNS_REQUIRE_NEAR(owner.x(6, j, 0), 1.0, 2.0e-12);
        WCNS_REQUIRE_NEAR(donor.x(0, j, 0), 0.0, 2.0e-12);
        WCNS_REQUIRE_NEAR(donor.y(0, j, 0), 1.0, 2.0e-12);
    }
}

// 验收多块轴置换连接按整条物理线切分守恒权重，且共享面加权贡献抵消。
void test_global_conservation_weights()
{
    using namespace wcns;
    const auto mesh = make_geometry_mesh();
    for (const auto kind : {AlgorithmProfileKind::PhengleiWcns, AlgorithmProfileKind::Scmm6Wcns}) {
        const auto profile = ProfileFactory::create(kind);
        const auto weights = GlobalConservationWeights::build(mesh, profile);
        const auto composite = build_line_conservation_weights(profile, 12);
        const auto transverse = build_line_conservation_weights(profile, 6);
        WCNS_REQUIRE(weights.maximum_line_residual() < 1.0e-11);
        WCNS_REQUIRE(weights.maximum_shared_face_mismatch() < 1.0e-11);
        for (const BlockId id : {BlockId {0}, BlockId {1}}) {
            const auto& block_weights = weights.block(id).cell;
            const auto extent = block_weights.interior_extent();
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    WCNS_REQUIRE(block_weights(i, j, 0) > 0.0);
                    const Real expected = id == 0
                        ? composite.cell_weights[static_cast<std::size_t>(i)]
                            * transverse.cell_weights[static_cast<std::size_t>(j)]
                        : transverse.cell_weights[static_cast<std::size_t>(i)]
                            * composite.cell_weights[static_cast<std::size_t>(j + 6)];
                    WCNS_REQUIRE_NEAR(block_weights(i, j, 0), expected, 2.0e-13);
                }
            }
        }
        WCNS_REQUIRE_THROWS(std::out_of_range, weights.block(99));
    }

    for (const auto kind : {AlgorithmProfileKind::PhengleiWcns,
                            AlgorithmProfileKind::Scmm6Wcns}) {
        auto partial = make_partial_c_grid_mesh();
        auto& block = partial.block(0);
        const auto profile = ProfileFactory::create(kind);
        WCNS_REQUIRE(!connection_side_is_fully_covered(block, Axis::J, Side::Lower));

        const auto weights = GlobalConservationWeights::build(partial, profile);
        const auto normal = build_line_conservation_weights(profile, 8);
        const auto tangent = build_line_conservation_weights(profile, 16);
        const auto& cells = weights.block(0).cell;
        for (int j = 0; j < 8; ++j) {
            for (int i = 0; i < 16; ++i) {
                WCNS_REQUIRE_NEAR(cells(i, j, 0),
                                  tangent.cell_weights[static_cast<std::size_t>(i)]
                                      * normal.cell_weights[static_cast<std::size_t>(j)],
                                  2.0e-13);
            }
        }
        WCNS_REQUIRE(weights.maximum_shared_face_mismatch() < 1.0e-13);

        // A connected sub-face on a mixed side must use the same bounded row
        // that generated the tensor-product conservation weights.
        const auto bounded = cached_line_operators(profile, 8).derivative_rows().front();
        WCNS_REQUIRE(inviscid_residual_stencil(block,
                                               profile,
                                               FluxDifferenceMode::Profile,
                                               Axis::J,
                                               {0, 0, 0})
                     == bounded);

        const auto metric = initialize_metric_field(block, profile).metric;
        InviscidFaceFluxField flux(block.cell_extent(), 2, kind, 19);
        flux.field(Axis::I).fill(0.0);
        flux.field(Axis::J).fill(0.0);
        for (int i = 0; i < 16; ++i) {
            flux.field(Axis::J)(i, 8, 0, density) = 0.2 + 0.013 * i;
        }
        for (int i = 0; i < 4; ++i) {
            const Real value = 0.07 * (i + 1);
            flux.field(Axis::J)(i, 0, 0, density) = value;
            flux.field(Axis::J)(15 - i, 0, 0, density) = -value;
        }
        for (int i = 4; i < 12; ++i) {
            flux.field(Axis::J)(i, 0, 0, density) = -0.15 + 0.009 * i;
        }

        compute_wcns_inviscid_residual(block, metric, flux, profile);
        Real actual = 0.0;
        for (int j = 0; j < 8; ++j) {
            for (int i = 0; i < 16; ++i) {
                actual += cells(i, j, 0) * metric.jacobian()(i, j, 0)
                    * block.flow.residual(i, j, 0, density);
            }
        }
        Real expected = 0.0;
        for (int i = 0; i < 16; ++i) {
            const Real lower_physical
                = i >= 4 && i < 12 ? flux.field(Axis::J)(i, 0, 0, density) : 0.0;
            expected -= tangent.cell_weights[static_cast<std::size_t>(i)]
                * (flux.field(Axis::J)(i, 8, 0, density) - lower_physical);
        }
        WCNS_REQUIRE_NEAR(actual, expected, 2.0e-12);

        const auto flux_plan = FaceFluxHaloPlan::build(partial, profile, 17);
        WCNS_REQUIRE(flux_plan.exchanges().size() == 2);
        WCNS_REQUIRE(flux_plan.exchanges()[0].receiver_block == 0);
        WCNS_REQUIRE(flux_plan.exchanges()[0].donor_block == 0);
        WCNS_REQUIRE(flux_plan.exchanges()[0].direction == 0);
        WCNS_REQUIRE(flux_plan.exchanges()[1].direction == 1);
        WCNS_REQUIRE(flux_plan.exchanges()[0].pairs.size()
                     == static_cast<std::size_t>(kind == AlgorithmProfileKind::PhengleiWcns ? 4
                                                                                           : 8));
        WCNS_REQUIRE(flux_plan.exchanges()[1].pairs.size()
                     == static_cast<std::size_t>(kind == AlgorithmProfileKind::PhengleiWcns ? 8
                                                                                           : 12));
        WCNS_REQUIRE(flux_plan.exchanges()[0].pairs.front().layer == 1);
        WCNS_REQUIRE(flux_plan.exchanges()[1].pairs.front().layer == 0);
        WCNS_REQUIRE(flux_plan.exchanges()[0].message_tag()
                     != flux_plan.exchanges()[1].message_tag());
    }

    const auto complete = make_geometry_mesh();
    WCNS_REQUIRE(
        connection_side_is_fully_covered(complete.block(0), Axis::I, Side::Upper));
    const auto complete_profile = ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns);
    const auto centered = inviscid_residual_stencil(complete.block(0),
                                                    complete_profile,
                                                    FluxDifferenceMode::Profile,
                                                    Axis::I,
                                                    {5, 2, 0});
    WCNS_REQUIRE(centered.back().first > 6);
}
