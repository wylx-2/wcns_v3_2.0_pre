#include "test_support.hpp"

#include <wcns/mesh/metrics.hpp>
#include <wcns/solver/physical_boundary.hpp>

#include <cmath>
#include <limits>

namespace {

wcns::GasModel test_gas()
{
    wcns::GasModelInput input;
    input.specific_gas_constant = 287.0;
    return wcns::GasModel::from_input(input);
}

wcns::ReferenceScales test_reference(const wcns::GasModel& gas)
{
    return wcns::ReferenceScales::derive({340.0, 1.2, 288.0, 1.0, 1.8e-5, {}, {}}, gas);
}

wcns::StructuredBlock make_wall_block()
{
    using namespace wcns;
    StructuredBlock block(0, "wall", 0, 2, 2, {7, 6, 1}, 3);
    const auto vertices = block.vertex_extent();
    for (int j = 0; j < vertices.nj; ++j) {
        for (int i = 0; i < vertices.ni; ++i) {
            block.coordinates.x(i, j, 0) = static_cast<Real>(i);
            block.coordinates.y(i, j, 0) = static_cast<Real>(j);
            block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    compute_metrics(block);
    block.boundaries.push_back({"i-lower-wall",
                                BoundaryType::SlipWall,
                                {Axis::I, Side::Lower},
                                {{0, 0, 0}, {0, vertices.nj - 1, 0}},
                                {{0, 0, 0}, {0, vertices.nj - 2, 0}},
                                {{0, 0, 0}, {0, vertices.nj - 2, 0}},
                                {}});
    return block;
}

} // namespace

// 验收温度型物理 ghost 的三层面状事务、派生状态版本及边角不写入约束。
void test_physical_ghost_state()
{
    using namespace wcns;
    const auto gas = test_gas();
    const auto reference = test_reference(gas);
    const NumericalFloors floors;
    auto block = make_wall_block();
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    block.flow.temperature_primitive.fill(nan);
    block.flow.primitive.fill(nan);
    block.flow.conservative.fill(nan);

    const TemperaturePrimitiveState interior {1.0, 2.0, 1.0, 0.0, 1.0};
    const auto conservative = thermodynamic_conservative(interior, gas, reference, floors, 2);
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            store_state(block.flow.conservative, {i, j, 0}, conservative);
        }
    }
    update_temperature_primitive_interior(block, gas, reference, floors);
    BoundaryDataMap data {{"i-lower-wall", BoundaryData {}}};
    const auto result = PhysicalGhostStateOperator::fill(block, data, gas, reference, floors, 7);
    WCNS_REQUIRE(result.version == 7);
    WCNS_REQUIRE(result.state_count == 3U * static_cast<std::size_t>(cells.nj));
    WCNS_REQUIRE(block.flow.physical_ghost_version == 7);
    for (int layer = 1; layer <= 3; ++layer) {
        const auto ghost = load_primitive(block.flow.primitive, {-layer, 2, 0});
        WCNS_REQUIRE_NEAR(ghost[0], 1.0, 1.0e-14);
        WCNS_REQUIRE_NEAR(ghost[1], -2.0, 1.0e-14);
        WCNS_REQUIRE_NEAR(ghost[2], 1.0, 1.0e-14);
        WCNS_REQUIRE(ghost[3] == 0.0);
        WCNS_REQUIRE_NEAR(block.flow.temperature_primitive(-layer, 2, 0, 4), 1.0, 1.0e-14);
        const auto recovered = temperature_primitive_from_conservative(
            load_conservative(block.flow.conservative, {-layer, 2, 0}), gas, reference, floors, 2);
        WCNS_REQUIRE_NEAR(recovered[1], -2.0, 1.0e-14);
    }
    for (int component = 0; component < fluid_components; ++component) {
        WCNS_REQUIRE(std::isnan(block.flow.temperature_primitive(-1, -1, 0, component)));
        WCNS_REQUIRE(std::isnan(block.flow.primitive(-1, -1, 0, component)));
        WCNS_REQUIRE(std::isnan(block.flow.conservative(-1, -1, 0, component)));
    }
    WCNS_REQUIRE_THROWS(PhysicsConfigurationError,
                        PhysicalGhostStateOperator::fill(block, data, gas, reference, floors, 0));
}

// 验收无粘物理边界面强开关仅反射固壁法向速度，弱模式保留重构外侧迹。
void test_inviscid_boundary_face_state()
{
    using namespace wcns;
    const auto gas = test_gas();
    const auto reference = test_reference(gas);
    const NumericalFloors floors;
    const BoundaryPatch patch {
        "wall", BoundaryType::NoSlipAdiabaticWall, {Axis::I, Side::Lower}, {}, {}, {}, {}};
    const PressurePrimitiveState interior {1.0, 2.0, 3.0, 0.0, 1.0};
    const PressurePrimitiveState reconstructed {0.9, -1.5, -4.0, 0.0, 0.8};
    const BoundaryData data;
    const auto strong = apply_inviscid_boundary_face_state(
        patch, interior, reconstructed, {-1.0, 0.0, 0.0}, data, {}, gas, reference, floors, 2);
    WCNS_REQUIRE(strong[0] == interior[0]);
    WCNS_REQUIRE(strong[1] == -interior[1]);
    WCNS_REQUIRE(strong[2] == interior[2]);
    WCNS_REQUIRE(strong[4] == interior[4]);

    InviscidBoundaryOptions weak;
    weak.strong_boundary_face_state = false;
    const auto unchanged = apply_inviscid_boundary_face_state(
        patch, interior, reconstructed, {-1.0, 0.0, 0.0}, data, weak, gas, reference, floors, 2);
    WCNS_REQUIRE(unchanged == reconstructed);
    WCNS_REQUIRE(weak.summary() == "strong_boundary_face_state=false");
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        apply_inviscid_boundary_face_state(
            patch, interior, reconstructed, {-2.0, 0.0, 0.0}, data, {}, gas, reference, floors, 2));

    BoundaryPatch farfield = patch;
    farfield.type = BoundaryType::Farfield;
    BoundaryData farfield_data;
    auto subsonic_interior = interior;
    subsonic_interior[1] = 0.2;
    const auto interior_temperature
        = temperature_primitive(subsonic_interior, gas, reference, floors, 2);
    auto target_temperature = interior_temperature;
    target_temperature[temperature_density] *= 1.1;
    target_temperature[temperature_velocity_x] += 0.15;
    target_temperature[temperature_value] *= 0.9;
    farfield_data.target_state = target_temperature;
    const auto target_pressure = pressure_primitive(target_temperature, gas, reference, floors, 2);
    const auto characteristic = apply_inviscid_boundary_face_state(farfield,
                                                                   subsonic_interior,
                                                                   reconstructed,
                                                                   {1.0, 0.0, 0.0},
                                                                   farfield_data,
                                                                   {},
                                                                   gas,
                                                                   reference,
                                                                   floors,
                                                                   2);
    WCNS_REQUIRE(characteristic != subsonic_interior);
    WCNS_REQUIRE(characteristic != target_pressure);

    farfield_data.farfield_point_vortex
        = FarfieldPointVortex {1.0, {{0.25, 0.0}}, 1.0};
    const auto corrected = farfield_target_at(farfield_data, {{0.25, 2.0, 0.0}}, 2);
    const double speed = std::hypot(target_temperature[temperature_velocity_x],
                                    target_temperature[temperature_velocity_y]);
    WCNS_REQUIRE_NEAR(corrected[temperature_velocity_x]
                          - target_temperature[temperature_velocity_x],
                      0.5 * speed / (4.0 * std::acos(-1.0)),
                      1.0e-14);
    WCNS_REQUIRE_NEAR(corrected[temperature_velocity_y],
                      target_temperature[temperature_velocity_y],
                      1.0e-14);
    WCNS_REQUIRE_THROWS(
        PhysicsConfigurationError,
        farfield_target_at(farfield_data, {{0.25, 2.0, 0.0}}, 3));
}
