#include "test_support.hpp"

#include <wcns/runtime/quantity_registry.hpp>
#include <wcns/solver/turbulence_fields.hpp>
#include <wcns/solver/turbulence_model.hpp>
#include <wcns/solver/turbulence_transport.hpp>
#include <wcns/solver/two_equation_models.hpp>
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

struct ErrorNorms {
    wcns::Real l1 = 0.0;
    wcns::Real l2 = 0.0;
    wcns::Real linf = 0.0;
};

ErrorNorms sa_manufactured_error(int cell_count)
{
    using namespace wcns;
    const SaNegativeConstants constants;
    const Real two_pi = 2.0 * std::acos(-1.0);
    const Real spacing = two_pi / static_cast<Real>(cell_count);
    constexpr Real velocity = 0.4;
    constexpr Real molecular_nu = 1.5e-5;
    constexpr Real base = 5.0e-5;
    constexpr Real amplitude = 1.0e-5;
    std::vector<Real> value(static_cast<std::size_t>(cell_count));
    for (int cell = 0; cell < cell_count; ++cell) {
        const Real x = (static_cast<Real>(cell) + 0.5) * spacing;
        value[static_cast<std::size_t>(cell)] = base + amplitude * std::sin(x);
    }

    std::vector<Real> flux(static_cast<std::size_t>(cell_count));
    for (int face = 0; face < cell_count; ++face) {
        const int left = face;
        const int right = (face + 1) % cell_count;
        const Real left_value = value[static_cast<std::size_t>(left)];
        const Real right_value = value[static_cast<std::size_t>(right)];
        const Real face_value = 0.5 * (left_value + right_value);
        const Real diffusion = (molecular_nu + face_value) / constants.sigma;
        flux[static_cast<std::size_t>(face)] = velocity * face_value
            - diffusion * (right_value - left_value) / spacing;
    }

    ErrorNorms result;
    for (int cell = 0; cell < cell_count; ++cell) {
        const int lower_face = (cell + cell_count - 1) % cell_count;
        const Real x = (static_cast<Real>(cell) + 0.5) * spacing;
        const Real exact_value = base + amplitude * std::sin(x);
        const Real first = amplitude * std::cos(x);
        const Real second = -amplitude * std::sin(x);

        TurbulenceCellContext context;
        context.mean_state = {{1.0, velocity, 0.0, 0.0, 1.0}};
        context.primitive_gradients[0] = {{0.0, -2.0, 0.0}};
        context.model_values = {exact_value};
        context.model_gradients = {{{first, 0.0, 0.0}}};
        context.molecular_kinematic_viscosity = molecular_nu;
        context.wall_distance = 0.7;
        context.dimension = 2;
        const Real model_source = evaluate_sa_negative(context, constants).source;
        const Real exact_flux_derivative = velocity * first
            - (first * first + (molecular_nu + exact_value) * second) / constants.sigma;
        const Real manufactured_source = exact_flux_derivative - model_source;
        const Real residual
            = -(flux[static_cast<std::size_t>(cell)]
                - flux[static_cast<std::size_t>(lower_face)])
                    / spacing
            + model_source + manufactured_source;
        const Real error = std::abs(residual);
        result.l1 += error;
        result.l2 += error * error;
        result.linf = std::max(result.linf, error);
    }
    result.l1 /= static_cast<Real>(cell_count);
    result.l2 = std::sqrt(result.l2 / static_cast<Real>(cell_count));
    return result;
}

ErrorNorms two_equation_manufactured_error(int cell_count,
                                            wcns::TurbulenceModelKind kind)
{
    using namespace wcns;
    const Real two_pi = 2.0 * std::acos(-1.0);
    const Real spacing = two_pi / static_cast<Real>(cell_count);
    constexpr Real velocity = 0.4;
    constexpr Real molecular_nu = 1.5e-5;
    const Real second_base = kind == TurbulenceModelKind::KOmegaSst ? 2.0 : 8.0e-3;
    const Real second_amplitude
        = kind == TurbulenceModelKind::KOmegaSst ? 0.1 : 1.0e-3;

    const auto exact_values = [&](Real x) {
        return std::array<Real, 2> {{0.04 + 0.005 * std::sin(x),
                                    second_base + second_amplitude * std::cos(2.0 * x)}};
    };
    const auto exact_gradients = [&](Real x) {
        return std::array<Real, 2> {{0.005 * std::cos(x),
                                    -2.0 * second_amplitude * std::sin(2.0 * x)}};
    };
    const auto context_at = [&](Real x,
                                const std::array<Real, 2>& values,
                                const std::array<Real, 2>& gradients) {
        TurbulenceCellContext context;
        context.mean_state = {{1.0, velocity, 0.0, 0.0, 1.0}};
        context.primitive_gradients[0] = {{0.08, 0.02, 0.0}};
        context.primitive_gradients[1] = {{-0.01, -0.03, 0.0}};
        context.model_values = {values[0], values[1]};
        context.model_gradients = {{{gradients[0], 0.0, 0.0}},
                                   {{gradients[1], 0.0, 0.0}}};
        context.molecular_kinematic_viscosity = molecular_nu;
        context.wall_distance = 0.4 + 0.02 * std::cos(x);
        context.reference_reynolds = 1.0e6;
        context.reference_mach = 0.2;
        context.heat_capacity_ratio = 1.4;
        context.dimension = 2;
        return context;
    };
    const auto evaluation_at = [&](Real x,
                                   const std::array<Real, 2>& values,
                                   const std::array<Real, 2>& gradients) {
        const auto context = context_at(x, values, gradients);
        return kind == TurbulenceModelKind::KOmegaSst
            ? evaluate_k_omega_sst(context)
            : evaluate_standard_k_epsilon(context);
    };
    const auto exact_flux = [&](Real x, int variable) {
        const auto values = exact_values(x);
        const auto gradients = exact_gradients(x);
        const auto evaluation = evaluation_at(x, values, gradients);
        return velocity * values[static_cast<std::size_t>(variable)]
            - evaluation.diffusion_coefficients[static_cast<std::size_t>(variable)]
            * gradients[static_cast<std::size_t>(variable)];
    };
    const auto exact_flux_derivative = [&](Real x, int variable) {
        constexpr Real delta = 1.0e-4;
        return (-exact_flux(x + 2.0 * delta, variable)
                + 8.0 * exact_flux(x + delta, variable)
                - 8.0 * exact_flux(x - delta, variable)
                + exact_flux(x - 2.0 * delta, variable))
            / (12.0 * delta);
    };

    std::array<std::vector<Real>, 2> values;
    std::array<std::vector<Real>, 2> gradients;
    for (auto& component : values) component.resize(static_cast<std::size_t>(cell_count));
    for (auto& component : gradients) component.resize(static_cast<std::size_t>(cell_count));
    for (int cell = 0; cell < cell_count; ++cell) {
        const Real x = (static_cast<Real>(cell) + 0.5) * spacing;
        const auto exact = exact_values(x);
        for (int variable = 0; variable < 2; ++variable) {
            values[static_cast<std::size_t>(variable)][static_cast<std::size_t>(cell)]
                = exact[static_cast<std::size_t>(variable)];
        }
    }
    for (int cell = 0; cell < cell_count; ++cell) {
        const int lower = (cell + cell_count - 1) % cell_count;
        const int upper = (cell + 1) % cell_count;
        for (int variable = 0; variable < 2; ++variable) {
            const auto& component = values[static_cast<std::size_t>(variable)];
            gradients[static_cast<std::size_t>(variable)][static_cast<std::size_t>(cell)]
                = (component[static_cast<std::size_t>(upper)]
                   - component[static_cast<std::size_t>(lower)])
                / (2.0 * spacing);
        }
    }

    std::array<std::vector<Real>, 2> fluxes;
    for (auto& component : fluxes) component.resize(static_cast<std::size_t>(cell_count));
    for (int face = 0; face < cell_count; ++face) {
        const int left = face;
        const int left_left = (face + cell_count - 1) % cell_count;
        const int right = (face + 1) % cell_count;
        const Real left_x = (static_cast<Real>(left) + 0.5) * spacing;
        const Real right_x = (static_cast<Real>(right) + 0.5) * spacing;
        std::array<Real, 2> left_values;
        std::array<Real, 2> right_values;
        std::array<Real, 2> left_gradients;
        std::array<Real, 2> right_gradients;
        for (int variable = 0; variable < 2; ++variable) {
            left_values[static_cast<std::size_t>(variable)]
                = values[static_cast<std::size_t>(variable)][static_cast<std::size_t>(left)];
            right_values[static_cast<std::size_t>(variable)]
                = values[static_cast<std::size_t>(variable)][static_cast<std::size_t>(right)];
            left_gradients[static_cast<std::size_t>(variable)]
                = gradients[static_cast<std::size_t>(variable)][static_cast<std::size_t>(left)];
            right_gradients[static_cast<std::size_t>(variable)]
                = gradients[static_cast<std::size_t>(variable)][static_cast<std::size_t>(right)];
        }
        const auto left_evaluation
            = evaluation_at(left_x, left_values, left_gradients);
        const auto right_evaluation
            = evaluation_at(right_x, right_values, right_gradients);
        for (int variable = 0; variable < 2; ++variable) {
            const auto entry = static_cast<std::size_t>(variable);
            const Real upwind = 1.5 * values[entry][static_cast<std::size_t>(left)]
                - 0.5 * values[entry][static_cast<std::size_t>(left_left)];
            const Real diffusion = 0.5
                * (left_evaluation.diffusion_coefficients[entry]
                   + right_evaluation.diffusion_coefficients[entry]);
            const Real face_gradient
                = 0.5 * (left_gradients[entry] + right_gradients[entry]);
            fluxes[entry][static_cast<std::size_t>(face)]
                = velocity * upwind - diffusion * face_gradient;
        }
    }

    ErrorNorms result;
    for (int cell = 0; cell < cell_count; ++cell) {
        const int lower_face = (cell + cell_count - 1) % cell_count;
        const Real x = (static_cast<Real>(cell) + 0.5) * spacing;
        std::array<Real, 2> cell_values;
        std::array<Real, 2> cell_gradients;
        for (int variable = 0; variable < 2; ++variable) {
            cell_values[static_cast<std::size_t>(variable)]
                = values[static_cast<std::size_t>(variable)][static_cast<std::size_t>(cell)];
            cell_gradients[static_cast<std::size_t>(variable)]
                = gradients[static_cast<std::size_t>(variable)][static_cast<std::size_t>(cell)];
        }
        const auto discrete_evaluation = evaluation_at(x, cell_values, cell_gradients);
        const auto exact_evaluation
            = evaluation_at(x, exact_values(x), exact_gradients(x));
        for (int variable = 0; variable < 2; ++variable) {
            const auto entry = static_cast<std::size_t>(variable);
            const Real manufactured_source
                = exact_flux_derivative(x, variable) - exact_evaluation.source[entry];
            const Real residual
                = -(fluxes[entry][static_cast<std::size_t>(cell)]
                    - fluxes[entry][static_cast<std::size_t>(lower_face)])
                    / spacing
                + discrete_evaluation.source[entry] + manufactured_source;
            const Real error = std::abs(residual);
            result.l1 += error;
            result.l2 += error * error;
            result.linf = std::max(result.linf, error);
        }
    }
    result.l1 /= static_cast<Real>(2 * cell_count);
    result.l2 = std::sqrt(result.l2 / static_cast<Real>(2 * cell_count));
    return result;
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
    auto positive_descriptor = descriptors[0];
    positive_descriptor.lower_bound = 1.0e-12;
    const auto unchanged
        = project_positive_turbulence_state(4.0e-6, 2.0, positive_descriptor);
    WCNS_REQUIRE(!unchanged.repaired);
    WCNS_REQUIRE_NEAR(unchanged.specific, 2.0e-6, 0.0);
    const auto projected
        = project_positive_turbulence_state(-3.0e-6, 2.0, positive_descriptor);
    WCNS_REQUIRE(projected.repaired);
    WCNS_REQUIRE(projected.specific > positive_descriptor.lower_bound);
    WCNS_REQUIRE(projected.specific
                 == std::nextafter(positive_descriptor.lower_bound,
                                   std::numeric_limits<Real>::infinity()));
    WCNS_REQUIRE_NEAR(projected.conservative, 2.0 * projected.specific, 0.0);
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        project_positive_turbulence_state(1.0, 1.0, descriptors[2]));
    WCNS_REQUIRE_THROWS(
        std::invalid_argument,
        project_positive_turbulence_state(
            std::numeric_limits<Real>::quiet_NaN(), 1.0, positive_descriptor));
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
    WCNS_REQUIRE(model->diffusion_coefficients({}).empty());
    WCNS_REQUIRE(model->source_linearization({}).source.empty());
    WCNS_REQUIRE(turbulence_model_family(TurbulenceModelKind::SaNegative)
                 == TurbulenceModelFamily::RansTransport);
    WCNS_REQUIRE(turbulence_model_family(TurbulenceModelKind::Wale)
                 == TurbulenceModelFamily::LesAlgebraic);
    WCNS_REQUIRE_THROWS(std::invalid_argument, turbulence_model_kind("unknown"));

    TurbulenceModelConfig sst_config;
    sst_config.kind = TurbulenceModelKind::KOmegaSst;
    const auto sst_model = registry.create(sst_config);
    WCNS_REQUIRE(sst_model->fields().size() == 10);
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

// 独立标量式检查 SST-2003m 基线和标准高 Re k-epsilon 的源、扩散与远场换算。
void test_two_equation_turbulence_models()
{
    using namespace wcns;
    TurbulenceCellContext context;
    context.mean_state = {{1.2, 0.3, 0.0, 0.0, 1.0}};
    context.model_values = {0.04, 2.0};
    context.model_gradients = {{{0.0, 0.0, 0.0}}, {{0.0, 0.0, 0.0}}};
    context.molecular_kinematic_viscosity = 1.0e-5;
    context.wall_distance = 0.1;
    context.reference_reynolds = 1.0e6;
    context.reference_mach = 0.2;
    context.heat_capacity_ratio = 1.4;
    context.dimension = 2;

    const auto sst = evaluate_k_omega_sst(context);
    WCNS_REQUIRE_NEAR(sst.blending_f1, 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(sst.blending_f2, 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(sst.eddy_kinematic_viscosity, 0.02, 2.0e-15);
    WCNS_REQUIRE_NEAR(sst.production, 0.0, 0.0);
    WCNS_REQUIRE_NEAR(sst.source[0], -0.00864, 2.0e-15);
    WCNS_REQUIRE_NEAR(sst.source[1], -0.36, 2.0e-14);
    WCNS_REQUIRE_NEAR(sst.diffusion_coefficients[0], 0.020412, 2.0e-15);
    WCNS_REQUIRE_NEAR(sst.diffusion_coefficients[1], 0.012012, 2.0e-15);
    WCNS_REQUIRE_NEAR(sst.source_jacobian[0], -0.18, 2.0e-15);
    WCNS_REQUIRE_NEAR(sst.source_jacobian[1], -0.0036, 2.0e-15);

    const SstConstants sst_constants;
    WCNS_REQUIRE_NEAR(sst_constants.production_limiter, 10.0, 0.0);
    WCNS_REQUIRE_NEAR(sst_constants.cross_diffusion_floor, 1.0e-10, 0.0);
    WCNS_REQUIRE_NEAR(sst_constants.gamma1, 5.0 / 9.0, 0.0);
    WCNS_REQUIRE_NEAR(sst_constants.gamma2, 0.44, 0.0);

    auto blend_context = context;
    blend_context.mean_state[temperature_density] = 1.0;
    blend_context.model_values = {1.0e-4, 1.0};
    blend_context.molecular_kinematic_viscosity = 1.0e-8;
    blend_context.wall_distance = 1.0;
    const auto blended = evaluate_k_omega_sst(blend_context);
    const Real expected_arg1 = std::sqrt(1.0e-4) / 0.09;
    const Real expected_arg2 = 2.0 * expected_arg1;
    WCNS_REQUIRE_NEAR(blended.blending_f1,
                      std::tanh(std::pow(expected_arg1, 4.0)),
                      2.0e-16);
    WCNS_REQUIRE_NEAR(blended.blending_f2,
                      std::tanh(expected_arg2 * expected_arg2),
                      2.0e-16);

    auto rotation_context = context;
    rotation_context.primitive_gradients[0][1] = 4.0;
    rotation_context.primitive_gradients[1][0] = -4.0;
    const auto rotation = evaluate_k_omega_sst(rotation_context);
    WCNS_REQUIRE_NEAR(rotation.eddy_kinematic_viscosity, 0.02, 2.0e-15);
    WCNS_REQUIRE_NEAR(rotation.production, 0.0, 0.0);

    auto limited_context = context;
    limited_context.primitive_gradients[0][1] = 100.0;
    const auto limited = evaluate_k_omega_sst(limited_context);
    WCNS_REQUIRE_NEAR(limited.production, 10.0 * 0.09 * 2.0 * 0.04, 2.0e-15);
    WCNS_REQUIRE_NEAR(sst.source_jacobian[3], -0.3, 2.0e-15);

    TurbulenceModelConfig sst_config;
    sst_config.kind = TurbulenceModelKind::KOmegaSst;
    const auto sst_farfield = two_equation_farfield_values(sst_config);
    WCNS_REQUIRE_NEAR(sst_farfield[0], 1.5e-4, 2.0e-18);
    WCNS_REQUIRE(sst_farfield[1] > 0.0);
    const auto sst_model = TurbulenceModelRegistry::create_builtin().create(sst_config);
    WCNS_REQUIRE(sst_model->source_linearization(context).jacobian.size() == 4);
    WCNS_REQUIRE(sst_model->diffusion_coefficients(context).size() == 2);
    WCNS_REQUIRE(sst_model->viscous_contribution(context).eddy_viscosity > 0.0);

    context.model_values = {0.04, 0.008};
    const auto ke = evaluate_standard_k_epsilon(context);
    WCNS_REQUIRE_NEAR(ke.eddy_kinematic_viscosity, 0.018, 2.0e-15);
    WCNS_REQUIRE_NEAR(ke.source[0], -0.0096, 2.0e-15);
    WCNS_REQUIRE_NEAR(ke.source[1], -0.0036864, 2.0e-15);
    WCNS_REQUIRE_NEAR(ke.diffusion_coefficients[0], 0.021612, 2.0e-15);
    WCNS_REQUIRE_NEAR(ke.diffusion_coefficients[1], 0.016627384615384615, 2.0e-15);
    TurbulenceModelConfig ke_config;
    ke_config.kind = TurbulenceModelKind::KEpsilon;
    ke_config.wall_treatment = WallTreatment::WallFunction;
    ke_config.experimental = true;
    const auto ke_model = TurbulenceModelRegistry::create_builtin().create(ke_config);
    WCNS_REQUIRE(ke_model->fields().size() == 7);
    WCNS_REQUIRE(two_equation_farfield_values(ke_config)[1] > 0.0);

    const auto order = [](Real coarse, Real fine) {
        return std::log(coarse / fine) / std::log(2.0);
    };
    for (const auto kind : {TurbulenceModelKind::KOmegaSst,
                            TurbulenceModelKind::KEpsilon}) {
        const auto coarse = two_equation_manufactured_error(32, kind);
        const auto medium = two_equation_manufactured_error(64, kind);
        const auto fine = two_equation_manufactured_error(128, kind);
        WCNS_REQUIRE(order(coarse.l1, medium.l1) > 1.9);
        WCNS_REQUIRE(order(medium.l1, fine.l1) > 1.9);
        WCNS_REQUIRE(order(coarse.l2, medium.l2) > 1.9);
        WCNS_REQUIRE(order(medium.l2, fine.l2) > 1.9);
        WCNS_REQUIRE(order(coarse.linf, medium.linf) > 1.8);
        WCNS_REQUIRE(order(medium.linf, fine.linf) > 1.8);
    }

    context.model_values[0] = 0.0;
    WCNS_REQUIRE_THROWS(PhysicsError, evaluate_standard_k_epsilon(context));
}

// 对照 NASA TMR SA-neg 的正/负分支参考点，并检查自动微分源 Jacobian。
void test_sa_negative_model()
{
    using namespace wcns;
    SaNegativeConstants constants;
    WCNS_REQUIRE_NEAR(constants.cw1(), 3.2390678167757287, 2.0e-15);

    TurbulenceCellContext context;
    context.mean_state = {{1.2, 0.3, -0.1, 0.0, 1.0}};
    context.primitive_gradients[0] = {{0.0, 2.0, 0.0}};
    context.primitive_gradients[1] = {{-1.0, 0.0, 0.0}};
    context.model_values = {4.5e-5};
    context.model_gradients = {{{1.0e-3, -2.0e-3, 0.0}}};
    context.molecular_kinematic_viscosity = 1.5e-5;
    context.wall_distance = 0.03;
    context.reference_reynolds = 6.0e6;
    context.reference_mach = 0.15;
    context.heat_capacity_ratio = 1.4;
    context.dimension = 2;

    const auto positive = evaluate_sa_negative(context);
    WCNS_REQUIRE(!positive.negative_branch);
    WCNS_REQUIRE_NEAR(positive.chi, 3.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(positive.fv1, 0.07014608571851676, 2.0e-15);
    WCNS_REQUIRE_NEAR(positive.fv2, -1.4784411615093869, 3.0e-15);
    WCNS_REQUIRE_NEAR(positive.ft2, 0.013330795845890767, 2.0e-15);
    WCNS_REQUIRE_NEAR(positive.modified_vorticity, 2.560249505797327, 3.0e-15);
    WCNS_REQUIRE_NEAR(positive.r, 0.11617695780693105, 2.0e-15);
    WCNS_REQUIRE_NEAR(positive.fw, 0.08153502509866944, 2.0e-15);
    WCNS_REQUIRE_NEAR(positive.diffusion_coefficient, 9.0e-5, 1.0e-19);
    WCNS_REQUIRE_NEAR(positive.eddy_kinematic_viscosity,
                      3.1565738573332544e-6,
                      1.0e-20);
    WCNS_REQUIRE_NEAR(positive.production_source, 1.5403012689802298e-5, 2.0e-19);
    WCNS_REQUIRE_NEAR(positive.destruction_source, -5.700418880258732e-7, 2.0e-20);
    WCNS_REQUIRE_NEAR(positive.cross_diffusion_source, 4.665e-6, 1.0e-20);
    WCNS_REQUIRE_NEAR(positive.source, 1.9497970801776424e-5, 2.0e-19);

    const auto derivative_check = [&](Real value) {
        context.model_values[0] = value;
        const auto center = evaluate_sa_negative(context);
        const Real step = 1.0e-7 * std::max(std::abs(value), context.molecular_kinematic_viscosity);
        context.model_values[0] = value + step;
        const Real upper = evaluate_sa_negative(context).source;
        context.model_values[0] = value - step;
        const Real lower = evaluate_sa_negative(context).source;
        context.model_values[0] = value;
        WCNS_REQUIRE_NEAR(center.source_derivative,
                          (upper - lower) / (2.0 * step),
                          2.0e-7 * std::max(1.0, std::abs(center.source_derivative)));
    };
    derivative_check(4.5e-5);

    context.model_values[0] = -2.25e-5;
    const auto negative = evaluate_sa_negative(context);
    WCNS_REQUIRE(negative.negative_branch);
    WCNS_REQUIRE_NEAR(negative.fn, 0.6516129032258065, 2.0e-15);
    WCNS_REQUIRE_NEAR(negative.diffusion_coefficient, 5.080645161290326e-7, 2.0e-21);
    WCNS_REQUIRE(negative.eddy_kinematic_viscosity == 0.0);
    WCNS_REQUIRE_NEAR(negative.production_source, 1.82925e-6, 2.0e-20);
    WCNS_REQUIRE_NEAR(negative.destruction_source, 1.8219756469363482e-6, 2.0e-20);
    WCNS_REQUIRE_NEAR(negative.source, 8.316225646936349e-6, 2.0e-20);
    derivative_check(-2.25e-5);

    context.model_values[0] = -7.1 * context.molecular_kinematic_viscosity;
    const auto formerly_singular_negative = evaluate_sa_negative(context);
    WCNS_REQUIRE(formerly_singular_negative.negative_branch);
    WCNS_REQUIRE(formerly_singular_negative.eddy_kinematic_viscosity == 0.0);
    WCNS_REQUIRE(std::isfinite(formerly_singular_negative.source));
    context.model_values[0] = -2.25e-5;

    TurbulenceModelConfig config;
    config.kind = TurbulenceModelKind::SaNegative;
    const auto model = TurbulenceModelRegistry::create_builtin().create(config);
    WCNS_REQUIRE(model->fields().size() == 6);
    WCNS_REQUIRE(model->diffusion_coefficients(context).size() == 1);
    WCNS_REQUIRE(model->viscous_contribution(context).eddy_viscosity == 0.0);
    const auto source = model->source_linearization(context);
    WCNS_REQUIRE_NEAR(source.source[0], 1.2 * negative.source, 2.0e-20);
    WCNS_REQUIRE_NEAR(source.jacobian[0], negative.source_derivative, 2.0e-12);
    WCNS_REQUIRE_NEAR(sa_negative_farfield_value(config, 6.0e6), 5.0e-7, 1.0e-21);

    WCNS_REQUIRE_NEAR(local_implicit_turbulence_increment(-10.0, -100.0, 0.1),
                      -1.0 / 11.0,
                      2.0e-16);
    WCNS_REQUIRE_THROWS(PhysicsError,
                        local_implicit_turbulence_increment(1.0, 20.0, 0.1));

    const auto coarse = sa_manufactured_error(32);
    const auto medium = sa_manufactured_error(64);
    const auto fine = sa_manufactured_error(128);
    const auto order = [](Real first, Real second) {
        return std::log(first / second) / std::log(2.0);
    };
    WCNS_REQUIRE(order(coarse.l1, medium.l1) > 1.9);
    WCNS_REQUIRE(order(medium.l1, fine.l1) > 1.9);
    WCNS_REQUIRE(order(coarse.l2, medium.l2) > 1.9);
    WCNS_REQUIRE(order(medium.l2, fine.l2) > 1.9);
    WCNS_REQUIRE(order(coarse.linf, medium.linf) > 1.9);
    WCNS_REQUIRE(order(medium.linf, fine.linf) > 1.9);
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
