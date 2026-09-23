#include "test_support.hpp"

#include <wcns/runtime/boundary_output.hpp>

void test_boundary_face_physics()
{
    using namespace wcns;
    BoundaryOutputConfig config;
    config.reference_pressure = 1.0;
    config.reference_density = 2.0;
    config.reference_velocity = {{3.0, 0.0, 0.0}};
    config.reference_area = 1.0;
    config.reference_length = 1.0;
    config.drag_direction = {{1.0, 0.0, 0.0}};
    config.lift_direction = {{0.0, 1.0, 0.0}};
    config.tangent_direction = {{0.0, 1.0, 0.0}};

    const auto face = evaluate_boundary_face_physics(2.0,
                                                     1.5,
                                                     1.2,
                                                     {{1.0, 0.0, 0.0}},
                                                     {{3.0, 4.0, 0.0}},
                                                     2.0,
                                                     {{5.0, 0.0, 0.0}},
                                                     10.0,
                                                     true,
                                                     config);
    WCNS_REQUIRE_NEAR(face.pressure_coefficient, 1.0 / 9.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(face.pressure_traction[0], 2.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(face.pressure_traction[1], 0.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(face.viscous_traction[0], -0.3, 1.0e-15);
    WCNS_REQUIRE_NEAR(face.viscous_traction[1], -0.4, 1.0e-15);
    WCNS_REQUIRE_NEAR(face.total_traction[0], 1.7, 1.0e-15);
    WCNS_REQUIRE_NEAR(face.skin_friction_coefficient, -0.4 / 9.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(face.heat_flux_into_wall, -1.0, 1.0e-15);

    auto resolved_wall = face;
    populate_boundary_wall_units(
        resolved_wall, {{1.0, 0.0, 0.0}}, 1.25, 0.01, 100.0);
    const Real expected_friction_velocity = std::sqrt(0.4 / 1.25);
    WCNS_REQUIRE_NEAR(
        resolved_wall.friction_velocity, expected_friction_velocity, 1.0e-15);
    WCNS_REQUIRE_NEAR(resolved_wall.wall_y_plus,
                      expected_friction_velocity * 0.01 * 1.25 * 100.0 / 1.2,
                      1.0e-15);
    WCNS_REQUIRE_NEAR(resolved_wall.wall_y_plus_class, 0.0, 0.0);
    auto wall_function = face;
    populate_boundary_wall_units(
        wall_function, {{1.0, 0.0, 0.0}}, 1.25, 0.01, 100.0, 50.0);
    WCNS_REQUIRE_NEAR(wall_function.wall_y_plus, 50.0, 0.0);
    WCNS_REQUIRE_NEAR(wall_function.friction_velocity, 48.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(wall_function.wall_y_plus_class, 2.0, 0.0);
    WCNS_REQUIRE_THROWS(
        PhysicsError,
        populate_boundary_wall_units(
            wall_function, {{1.0, 0.0, 0.0}}, 1.25, 0.0, 100.0));

    const auto flipped = evaluate_boundary_face_physics(2.0,
                                                        1.5,
                                                        1.2,
                                                        {{-1.0, 0.0, 0.0}},
                                                        {{-3.0, -4.0, 0.0}},
                                                        2.0,
                                                        {{5.0, 0.0, 0.0}},
                                                        10.0,
                                                        true,
                                                        config);
    for (std::size_t component = 0; component < 3; ++component) {
        WCNS_REQUIRE_NEAR(
            flipped.pressure_traction[component], -face.pressure_traction[component], 1.0e-15);
        WCNS_REQUIRE_NEAR(
            flipped.viscous_traction[component], -face.viscous_traction[component], 1.0e-15);
    }
    WCNS_REQUIRE_NEAR(flipped.skin_friction_coefficient, -face.skin_friction_coefficient, 1.0e-15);
    WCNS_REQUIRE_NEAR(flipped.heat_flux_into_wall, -face.heat_flux_into_wall, 1.0e-15);

    const auto inviscid = evaluate_boundary_face_physics(2.0,
                                                         1.5,
                                                         1.2,
                                                         {{1.0, 0.0, 0.0}},
                                                         {{0.0, 0.0, 0.0}},
                                                         2.0,
                                                         {{0.0, 0.0, 0.0}},
                                                         10.0,
                                                         false,
                                                         config);
    WCNS_REQUIRE_NEAR(inviscid.viscous_traction[0], 0.0, 0.0);
    WCNS_REQUIRE_NEAR(inviscid.heat_flux_into_wall, 0.0, 0.0);

    GasModelInput gas_input;
    gas_input.specific_gas_constant = 4.0;
    const auto gas = GasModel::from_input(gas_input);
    const auto reference = ReferenceScales::derive({2.0, 2.0, 3.0, 5.0, 7.0, {}, {}}, gas);
    QuantityContext dimensional {
        gas, reference, NumericalFloors {}, TransportModel(TransportConfig {}), true};
    const auto scales_2d = boundary_output_scales(dimensional, 2);
    WCNS_REQUIRE_NEAR(scales_2d.coordinate, 5.0, 0.0);
    WCNS_REQUIRE_NEAR(scales_2d.area, 5.0, 0.0);
    WCNS_REQUIRE_NEAR(scales_2d.pressure, 8.0, 0.0);
    WCNS_REQUIRE_NEAR(scales_2d.temperature, 3.0, 0.0);
    WCNS_REQUIRE_NEAR(scales_2d.viscosity, 7.0, 0.0);
    WCNS_REQUIRE_NEAR(scales_2d.velocity, 2.0, 0.0);
    WCNS_REQUIRE_NEAR(scales_2d.force, 40.0, 0.0);
    WCNS_REQUIRE_NEAR(scales_2d.moment, 200.0, 0.0);
    const Real nondimensional_force = 1.23456789;
    WCNS_REQUIRE_NEAR(
        nondimensional_force * scales_2d.force / scales_2d.force, nondimensional_force, 1.0e-15);
    const auto scales_3d = boundary_output_scales(dimensional, 3);
    WCNS_REQUIRE_NEAR(scales_3d.area, 25.0, 0.0);
    WCNS_REQUIRE_NEAR(scales_3d.force, 200.0, 0.0);
    WCNS_REQUIRE_NEAR(scales_3d.moment, 1000.0, 0.0);
    QuantityContext nondimensional = dimensional;
    nondimensional.dimensional = false;
    const auto unit_scales = boundary_output_scales(nondimensional, 2);
    WCNS_REQUIRE_NEAR(unit_scales.force, 1.0, 0.0);
}
