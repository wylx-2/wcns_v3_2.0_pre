#include "test_support.hpp"

#include <wcns/runtime/quantity_registry.hpp>
#include <wcns/solver/turbulence_fields.hpp>
#include <wcns/solver/turbulence_model.hpp>
#include <wcns/solver/turbulence_transport.hpp>
#include <wcns/solver/viscous_flux.hpp>
#include <wcns/solver/wall_distance.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

wcns::GasModel test_gas()
{
    wcns::GasModelInput input;
    input.specific_gas_constant = 287.0;
    return wcns::GasModel::from_input(input);
}

wcns::ReferenceScales test_reference(const wcns::GasModel& gas)
{
    wcns::ReferenceInput input;
    input.velocity = 100.0;
    input.density = 1.0;
    input.temperature = 300.0;
    input.length = 2.0;
    input.viscosity = 1.0e-5;
    return wcns::ReferenceScales::derive(input, gas);
}

wcns::Real manufactured_transport_error(int cell_count, int component_count)
{
    using namespace wcns;
    const Real two_pi = 2.0 * std::acos(-1.0);
    const Real spacing = two_pi / static_cast<Real>(cell_count);
    const Real velocity = 0.7;
    const Real diffusion = 0.08;
    std::vector<std::vector<Real>> state(
        static_cast<std::size_t>(component_count),
        std::vector<Real>(static_cast<std::size_t>(cell_count)));
    for (int component = 0; component < component_count; ++component) {
        const Real wave_number = static_cast<Real>(component + 1);
        for (int cell = 0; cell < cell_count; ++cell) {
            const Real x = (static_cast<Real>(cell) + 0.5) * spacing;
            state[static_cast<std::size_t>(component)][static_cast<std::size_t>(cell)]
                = 2.0 + std::sin(wave_number * x);
        }
    }

    std::vector<std::vector<Real>> face_flux(
        static_cast<std::size_t>(component_count),
        std::vector<Real>(static_cast<std::size_t>(cell_count)));
    for (int face = 0; face < cell_count; ++face) {
        const int left = face;
        const int right = (face + 1) % cell_count;
        std::vector<Real> reconstructed(static_cast<std::size_t>(component_count));
        std::vector<Real> gradients(static_cast<std::size_t>(component_count));
        for (int component = 0; component < component_count; ++component) {
            const auto& values = state[static_cast<std::size_t>(component)];
            reconstructed[static_cast<std::size_t>(component)]
                = 0.5 * (values[static_cast<std::size_t>(left)]
                         + values[static_cast<std::size_t>(right)]);
            gradients[static_cast<std::size_t>(component)]
                = (values[static_cast<std::size_t>(right)]
                   - values[static_cast<std::size_t>(left)])
                / spacing;
        }
        const auto flux = turbulence_face_flux(
            velocity,
            reconstructed,
            reconstructed,
            std::vector<Real>(static_cast<std::size_t>(component_count), diffusion),
            gradients);
        for (int component = 0; component < component_count; ++component) {
            face_flux[static_cast<std::size_t>(component)][static_cast<std::size_t>(face)]
                = flux.net[static_cast<std::size_t>(component)];
        }
    }

    Real error_squared = 0.0;
    for (int component = 0; component < component_count; ++component) {
        const Real wave_number = static_cast<Real>(component + 1);
        for (int cell = 0; cell < cell_count; ++cell) {
            const int lower_face = (cell + cell_count - 1) % cell_count;
            const Real x = (static_cast<Real>(cell) + 0.5) * spacing;
            const Real source = velocity * wave_number * std::cos(wave_number * x)
                + diffusion * wave_number * wave_number * std::sin(wave_number * x);
            const auto& fluxes = face_flux[static_cast<std::size_t>(component)];
            const Real residual
                = -(fluxes[static_cast<std::size_t>(cell)]
                    - fluxes[static_cast<std::size_t>(lower_face)])
                    / spacing
                + source;
            error_squared += residual * residual;
        }
    }
    return std::sqrt(error_squared
                     / static_cast<Real>(cell_count * component_count));
}

wcns::Real ssprk3_decay_error(wcns::Real time_step)
{
    using namespace wcns;
    constexpr Real decay = 0.7;
    constexpr Real final_time = 1.0;
    const int step_count = static_cast<int>(std::llround(final_time / time_step));
    time_step = final_time / static_cast<Real>(step_count);
    Real value = 1.0;
    for (int step = 0; step < step_count; ++step) {
        const Real first = value - time_step * decay * value;
        const Real second
            = 0.75 * value + 0.25 * (first - time_step * decay * first);
        value = value / 3.0 + 2.0 / 3.0 * (second - time_step * decay * second);
    }
    return std::abs(value - std::exp(-decay * final_time));
}

std::vector<wcns::WallPrimitive> circle_segments(int segment_count)
{
    using namespace wcns;
    const Real two_pi = 2.0 * std::acos(-1.0);
    std::vector<WallPrimitive> result;
    for (int segment = 0; segment < segment_count; ++segment) {
        const Real first_angle = two_pi * static_cast<Real>(segment)
            / static_cast<Real>(segment_count);
        const Real second_angle = two_pi * static_cast<Real>(segment + 1)
            / static_cast<Real>(segment_count);
        result.push_back(WallPrimitive::segment(
            {{std::cos(first_angle), std::sin(first_angle), 0.0}},
            {{std::cos(second_angle), std::sin(second_angle), 0.0}}));
    }
    return result;
}

wcns::WallPoint sphere_point(wcns::Real latitude, wcns::Real longitude)
{
    const wcns::Real radius = std::cos(latitude);
    return {{radius * std::cos(longitude),
             radius * std::sin(longitude),
             std::sin(latitude)}};
}

std::vector<wcns::WallPrimitive> sphere_triangles(int latitude_bands)
{
    using namespace wcns;
    const Real pi = std::acos(-1.0);
    const int longitude_count = 2 * latitude_bands;
    std::vector<WallPrimitive> result;
    const WallPoint south {{0.0, 0.0, -1.0}};
    const WallPoint north {{0.0, 0.0, 1.0}};
    const auto ring_point = [&](int latitude_index, int longitude_index) {
        const Real latitude = -0.5 * pi
            + pi * static_cast<Real>(latitude_index) / static_cast<Real>(latitude_bands);
        const Real longitude = 2.0 * pi * static_cast<Real>(longitude_index)
            / static_cast<Real>(longitude_count);
        return sphere_point(latitude, longitude);
    };
    for (int longitude = 0; longitude < longitude_count; ++longitude) {
        const int next = (longitude + 1) % longitude_count;
        result.push_back(WallPrimitive::triangle(
            south, ring_point(1, next), ring_point(1, longitude)));
        for (int latitude = 1; latitude + 1 < latitude_bands; ++latitude) {
            const auto q00 = ring_point(latitude, longitude);
            const auto q10 = ring_point(latitude, next);
            const auto q11 = ring_point(latitude + 1, next);
            const auto q01 = ring_point(latitude + 1, longitude);
            result.push_back(WallPrimitive::triangle(q00, q10, q11));
            result.push_back(WallPrimitive::triangle(q00, q11, q01));
        }
        result.push_back(WallPrimitive::triangle(
            ring_point(latitude_bands - 1, longitude),
            ring_point(latitude_bands - 1, next),
            north));
    }
    return result;
}

} // namespace

// 验收模型场不改变平均流五分量布局，并可确定性打包、恢复和注册输出。
void test_turbulence_fields()
{
    using namespace wcns;
    const std::vector<TurbulenceFieldDescriptor> descriptors {
        {"k", TurbulenceFieldRole::Transported, TurbulenceFieldScale::VelocitySquared, true, 0.0},
        {"omega", TurbulenceFieldRole::Transported, TurbulenceFieldScale::InverseTime, true, 0.0},
        {"mu_model", TurbulenceFieldRole::Diagnostic, TurbulenceFieldScale::DynamicViscosity, false, 0.0},
    };
    StructuredBlock block(0, "model-fields", 0, 3, 3, {3, 2, 2}, 2);
    WCNS_REQUIRE(block.flow.conservative.components() == euler_components);
    WCNS_REQUIRE(block.turbulence.empty());
    block.turbulence.reset(block.cell_extent(), block.ghost_width(), descriptors);
    WCNS_REQUIRE(block.turbulence.components() == 3);
    WCNS_REQUIRE(block.flow.conservative.components() == euler_components);

    for (int k = 0; k < block.cell_extent().nk; ++k) {
        for (int j = 0; j < block.cell_extent().nj; ++j) {
            for (int i = 0; i < block.cell_extent().ni; ++i) {
                block.turbulence.at({i, j, k}, "k") = 1.0 + i;
                block.turbulence.at({i, j, k}, "omega") = 2.0 + j;
                block.turbulence.at({i, j, k}, "mu_model") = 3.0 + k;
            }
        }
    }
    block.turbulence.validate_interior();
    const auto payload = block.turbulence.pack_interior();
    TurbulenceFieldSet restored(block.cell_extent(), block.ghost_width(), descriptors);
    restore_turbulence_restart(capture_turbulence_restart(block.turbulence), restored);
    WCNS_REQUIRE(restored.pack_interior() == payload);
    WCNS_REQUIRE(restored.descriptor_signature() == block.turbulence.descriptor_signature());
    WCNS_REQUIRE_THROWS(std::invalid_argument, restored.unpack_interior({1.0}));
    auto wrong_restart = capture_turbulence_restart(block.turbulence);
    wrong_restart.descriptor_signature += ";different";
    WCNS_REQUIRE_THROWS(std::runtime_error,
                        restore_turbulence_restart(wrong_restart, restored));
    restored.at({0, 0, 0}, "k") = 0.0;
    WCNS_REQUIRE_THROWS(std::runtime_error, restored.validate_interior());

    FieldQuantityRegistry registry;
    for (const auto& descriptor : descriptors)
        registry.register_turbulence_field(descriptor);
    WCNS_REQUIRE(registry.contains("k"));
    WCNS_REQUIRE(registry.contains("omega"));
    WCNS_REQUIRE(registry.descriptor("mu_model").scale == QuantityScale::Viscosity);
    WCNS_REQUIRE_THROWS(std::invalid_argument, registry.register_turbulence_field(descriptors[0]));

    auto duplicate = descriptors;
    duplicate[1].name = "k";
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        (TurbulenceFieldSet(block.cell_extent(), 1, duplicate)));
    auto invalid_name = descriptors[0];
    invalid_name.name = "K bad";
    WCNS_REQUIRE_THROWS(std::invalid_argument, invalid_name.validate());
}

// 验收 none registry、模型族分类和供 LU-SGS 使用的源项分块对角符号。
void test_turbulence_model_interface()
{
    using namespace wcns;
    TurbulenceModelConfig none;
    const auto registry = TurbulenceModelRegistry::create_builtin();
    const auto model = registry.create(none);
    WCNS_REQUIRE(model->family() == TurbulenceModelFamily::None);
    WCNS_REQUIRE(model->fields().empty());
    WCNS_REQUIRE(model->viscous_contribution({}).eddy_viscosity == 0.0);
    WCNS_REQUIRE(model->source_linearization({}).source.empty());
    WCNS_REQUIRE(turbulence_model_family(TurbulenceModelKind::SaNegative)
                 == TurbulenceModelFamily::RansTransport);
    WCNS_REQUIRE(turbulence_model_family(TurbulenceModelKind::Wale)
                 == TurbulenceModelFamily::LesAlgebraic);
    WCNS_REQUIRE_THROWS(std::invalid_argument, turbulence_model_kind("unknown"));

    TurbulenceModelConfig unavailable;
    unavailable.kind = TurbulenceModelKind::SaNegative;
    WCNS_REQUIRE_THROWS(std::invalid_argument, registry.create(unavailable));
    TurbulenceModelConfig invalid_kepsilon;
    invalid_kepsilon.kind = TurbulenceModelKind::KEpsilon;
    WCNS_REQUIRE_THROWS(std::invalid_argument, invalid_kepsilon.validate());

    TurbulenceSourceLinearization source {{1.0, -2.0}, {3.0, 4.0, 5.0, 6.0}};
    source.validate(2);
    WCNS_REQUIRE(source.implicit_diagonal_block(10.0)
                 == std::vector<Real>({7.0, -4.0, -5.0, 4.0}));
    source.jacobian.pop_back();
    WCNS_REQUIRE_THROWS(std::invalid_argument, source.validate(2));
}

void test_turbulence_transport_interface()
{
    using namespace wcns;
    const auto positive
        = turbulence_face_flux(2.0, {1.0, 3.0}, {5.0, 7.0}, {0.5, 2.0}, {4.0, -1.0});
    WCNS_REQUIRE(positive.advective == std::vector<Real>({2.0, 6.0}));
    WCNS_REQUIRE(positive.diffusive == std::vector<Real>({2.0, -2.0}));
    WCNS_REQUIRE(positive.net == std::vector<Real>({0.0, 8.0}));
    const auto negative
        = turbulence_face_flux(-2.0, {1.0}, {5.0}, {0.0}, {9.0});
    WCNS_REQUIRE(negative.advective == std::vector<Real>({-10.0}));
    WCNS_REQUIRE(turbulence_boundary_ghost(
                     {2.0}, {1.0}, TurbulenceBoundaryRule::Dirichlet)
                 == std::vector<Real>({0.0}));
    WCNS_REQUIRE(turbulence_boundary_ghost(
                     {2.0}, {9.0}, TurbulenceBoundaryRule::Extrapolate)
                 == std::vector<Real>({2.0}));

    const std::vector<TurbulenceFieldDescriptor> descriptors {
        {"a", TurbulenceFieldRole::Transported, TurbulenceFieldScale::Dimensionless, true, 0.0},
        {"b", TurbulenceFieldRole::Transported, TurbulenceFieldScale::Dimensionless, false, -1.0},
    };
    WCNS_REQUIRE(admissible_turbulence_update({1.0, 0.0}, {-0.5, -1.0}, descriptors)
                 == std::vector<Real>({0.5, -1.0}));
    WCNS_REQUIRE_THROWS(
        std::runtime_error,
        admissible_turbulence_update({1.0, 0.0}, {-1.0, 0.0}, descriptors));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        turbulence_face_flux(1.0, {1.0}, {1.0}, {-1.0}, {0.0}));

    for (int component_count = 1; component_count <= 3; ++component_count) {
        const Real coarse = manufactured_transport_error(32, component_count);
        const Real fine = manufactured_transport_error(64, component_count);
        WCNS_REQUIRE(std::log(coarse / fine) / std::log(2.0) > 1.95);
    }
    const Real coarse_time = ssprk3_decay_error(0.1);
    const Real fine_time = ssprk3_decay_error(0.05);
    WCNS_REQUIRE(std::log(coarse_time / fine_time) / std::log(2.0) > 2.9);
}

// 验收零闭合逐位恢复层流结果，非零应力/热流只经唯一黏性通量入口相加。
void test_turbulence_viscous_coupling()
{
    using namespace wcns;
    const auto gas = test_gas();
    const auto reference = test_reference(gas);
    TransportModel transport(TransportConfig {});
    NumericalFloors floors;
    ViscousFaceTrace trace;
    trace.state = {{1.0, 0.2, -0.1, 0.0, 1.0}};
    trace.gradients[static_cast<int>(ViscousPrimitive::VelocityX)] = {{0.1, 0.2, 0.0}};
    trace.gradients[static_cast<int>(ViscousPrimitive::VelocityY)] = {{-0.05, 0.3, 0.0}};
    trace.gradients[static_cast<int>(ViscousPrimitive::VelocityZ)] = {{0.0, 0.0, 0.0}};
    trace.gradients[static_cast<int>(ViscousPrimitive::Temperature)] = {{0.4, -0.2, 0.0}};

    const auto laminar
        = compute_viscous_cartesian_flux(trace, transport, gas, reference, floors, 2);
    const auto explicit_zero = compute_viscous_cartesian_flux(
        trace, transport, gas, reference, floors, TurbulenceViscousContribution {}, 2);
    WCNS_REQUIRE(laminar.x == explicit_zero.x);
    WCNS_REQUIRE(laminar.y == explicit_zero.y);
    WCNS_REQUIRE(laminar.z == explicit_zero.z);
    WCNS_REQUIRE(laminar.viscosity == explicit_zero.viscosity);

    TurbulenceViscousContribution closure;
    closure.stress = {1.0, 2.0, 3.0, 0.5, 0.0, 0.0};
    closure.energy_heat_flux = {{4.0, -2.0, 0.0}};
    closure.eddy_viscosity = 7.0;
    const auto coupled
        = compute_viscous_cartesian_flux(trace, transport, gas, reference, floors, closure, 2);
    WCNS_REQUIRE_NEAR(coupled.x[momentum_x] - laminar.x[momentum_x], 1.0, 0.0);
    WCNS_REQUIRE_NEAR(coupled.x[momentum_y] - laminar.x[momentum_y], 0.5, 0.0);
    WCNS_REQUIRE_NEAR(coupled.y[momentum_y] - laminar.y[momentum_y], 2.0, 0.0);
    WCNS_REQUIRE_NEAR(coupled.x[total_energy] - laminar.x[total_energy], 4.15, 1.0e-14);
    WCNS_REQUIRE_NEAR(coupled.y[total_energy] - laminar.y[total_energy], -2.1, 1.0e-14);
    WCNS_REQUIRE_NEAR(coupled.viscosity - laminar.viscosity, 7.0, 0.0);
    closure.stress.xz = 1.0;
    WCNS_REQUIRE_THROWS(std::invalid_argument, closure.validate(2));
}

// 验收二维 segment、三维 triangle、BVH 最近元和等距时的确定性 primitive 选择。
void test_wall_distance()
{
    using namespace wcns;
    WallDistanceIndex segments({WallPrimitive::segment({{-1.0, 0.0, 0.0}},
                                                       {{1.0, 0.0, 0.0}}),
                                WallPrimitive::segment({{-1.0, 2.0, 0.0}},
                                                       {{1.0, 2.0, 0.0}})});
    WCNS_REQUIRE_NEAR(segments.query({{0.25, 0.5, 0.0}}).distance, 0.5, 1.0e-15);
    const auto tie = segments.query({{0.0, 1.0, 0.0}});
    WCNS_REQUIRE_NEAR(tie.distance, 1.0, 1.0e-15);
    WCNS_REQUIRE(tie.primitive == 0);

    WallDistanceIndex triangle({WallPrimitive::triangle({{0.0, 0.0, 0.0}},
                                                        {{1.0, 0.0, 0.0}},
                                                        {{0.0, 1.0, 0.0}})});
    WCNS_REQUIRE_NEAR(triangle.query({{0.25, 0.25, 2.0}}).distance, 2.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(triangle.query({{2.0, 0.0, 0.0}}).distance, 1.0, 1.0e-15);
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        triangle.query({{std::numeric_limits<Real>::quiet_NaN(), 0.0, 0.0}}));

    const Real circle_error_16
        = 1.0 - WallDistanceIndex(circle_segments(16)).query({{0.0, 0.0, 0.0}}).distance;
    const Real circle_error_32
        = 1.0 - WallDistanceIndex(circle_segments(32)).query({{0.0, 0.0, 0.0}}).distance;
    WCNS_REQUIRE(circle_error_16 / circle_error_32 > 3.9);
    const Real sphere_error_8
        = 1.0 - WallDistanceIndex(sphere_triangles(8)).query({{0.0, 0.0, 0.0}}).distance;
    const Real sphere_error_16
        = 1.0 - WallDistanceIndex(sphere_triangles(16)).query({{0.0, 0.0, 0.0}}).distance;
    WCNS_REQUIRE(sphere_error_8 / sphere_error_16 > 3.7);

    StructuredBlock block(0, "distance", 0, 2, 2, {3, 2, 1}, 1);
    block.cell_metrics.center_x(0, 0, 0) = -0.5;
    block.cell_metrics.center_y(0, 0, 0) = 0.25;
    block.cell_metrics.center_z(0, 0, 0) = 0.0;
    block.cell_metrics.center_x(1, 0, 0) = 0.5;
    block.cell_metrics.center_y(1, 0, 0) = 0.75;
    block.cell_metrics.center_z(1, 0, 0) = 0.0;
    Field<Real> distance(block.cell_extent(), 1, 1);
    segments.fill_cell_field(block, distance);
    WCNS_REQUIRE_NEAR(distance(0, 0, 0, 0), 0.25, 1.0e-15);
    WCNS_REQUIRE_NEAR(distance(1, 0, 0, 0), 0.75, 1.0e-15);

    StructuredBlock wall_block(1, "wall-block", 0, 2, 2, {3, 2, 1}, 1);
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < 3; ++i) {
            wall_block.coordinates.x(i, j, 0) = static_cast<Real>(i);
            wall_block.coordinates.y(i, j, 0) = static_cast<Real>(j);
            wall_block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    wall_block.boundaries.push_back({"wall",
                                     BoundaryType::NoSlipAdiabaticWall,
                                     {Axis::J, Side::Lower},
                                     {{0, 0, 0}, {2, 0, 0}},
                                     {{0, 0, 0}, {1, 0, 0}},
                                     {{0, 0, 0}, {1, 0, 0}},
                                     {}});
    const auto extracted
        = extract_wall_primitives(std::vector<StructuredBlock> {wall_block});
    WCNS_REQUIRE(extracted.size() == 2);
    WCNS_REQUIRE_NEAR(WallDistanceIndex(extracted).query({{1.5, 0.5, 0.0}}).distance,
                      0.5,
                      1.0e-15);
}
