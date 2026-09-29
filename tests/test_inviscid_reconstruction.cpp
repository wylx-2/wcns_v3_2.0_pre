#include "test_support.hpp"

#include <wcns/solver/riemann_solver.hpp>
#include <wcns/solver/wcns_reconstruction.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <memory>
#include <string_view>
#include <vector>

namespace {

wcns::GasModel reconstruction_gas()
{
    wcns::GasModelInput input;
    input.specific_gas_constant = 287.0;
    return wcns::GasModel::from_input(input);
}

wcns::ReferenceScales reconstruction_reference(const wcns::GasModel& gas)
{
    return wcns::ReferenceScales::derive({340.0, 1.2, 288.0, 1.0, 1.8e-5, {}, {}}, gas);
}

class CustomAverageReconstruction final : public wcns::IReconstructionScheme {
public:
    [[nodiscard]] std::string_view name() const noexcept override { return "custom_average"; }

    [[nodiscard]] wcns::StencilRequirement stencil_requirement() const noexcept override
    {
        return {};
    }

    [[nodiscard]] wcns::Real reconstruct_scalar(wcns::ScalarStencilView stencil,
                                                wcns::TraceSide side,
                                                const wcns::ReconstructionContext&) const override
    {
        if (stencil.size() != 6) throw std::invalid_argument("custom stencil size");
        return side == wcns::TraceSide::Left ? 0.5 * (stencil[2] + stencil[3])
                                             : 0.5 * (stencil[3] + stencil[2]);
    }
};

class CustomCentralRiemann final : public wcns::IRiemannSolver {
public:
    [[nodiscard]] std::string_view name() const noexcept override { return "custom_central"; }

    [[nodiscard]] wcns::RiemannResult solve(const wcns::PressurePrimitiveState& left,
                                            const wcns::PressurePrimitiveState& right,
                                            wcns::Normal3 unit_normal,
                                            const wcns::GasModel& gas,
                                            const wcns::NumericalFloors& floors) const override
    {
        const wcns::IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
        const auto left_flux = wcns::euler_flux(left, unit_normal, ideal);
        const auto right_flux = wcns::euler_flux(right, unit_normal, ideal);
        wcns::RiemannResult result;
        for (int component = 0; component < wcns::euler_components; ++component) {
            result.flux_per_unit_area[static_cast<std::size_t>(component)] = 0.5
                * (left_flux[static_cast<std::size_t>(component)]
                   + right_flux[static_cast<std::size_t>(component)]);
        }
        result.requested_solver = "custom_central";
        result.used_solver = "custom_central";
        return result;
    }
};

} // namespace

// 验收阶段 L 的重构/Riemann 注册表可增加自定义策略且拒绝未知或重复名称。
void test_stage_l_algorithm_registries()
{
    using namespace wcns;
    auto reconstruction_registry = ReconstructionRegistry::with_builtins();
    WCNS_REQUIRE(reconstruction_registry.names()
                 == std::vector<std::string>(
                     {"linear5", "mdcd_hybrid", "mdcd_linear", "weno_js", "weno_z", "zero_order"}));
    reconstruction_registry.register_scheme(
        "custom_average", [] { return std::make_unique<CustomAverageReconstruction>(); });
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        reconstruction_registry.register_scheme("custom_average", [] {
                            return std::make_unique<CustomAverageReconstruction>();
                        }));
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        reconstruction_registry.create("missing_reconstruction"));
    auto custom = reconstruction_registry.create("custom_average");
    const std::array<Real, 6> stencil {{0.0, 1.0, 2.0, 4.0, 8.0, 16.0}};
    WCNS_REQUIRE_NEAR(custom->reconstruct_scalar(stencil, TraceSide::Left, {}), 3.0, 0.0);
    const std::array<Real, 5> wrong_stencil {{0.0, 1.0, 2.0, 3.0, 4.0}};
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        reconstruction_registry.create("linear5")->reconstruct_scalar(
                            wrong_stencil, TraceSide::Left, {}));
    ReconstructionConfig custom_config;
    custom_config.scheme = "custom_average";
    custom_config.validate(reconstruction_registry);
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        custom_config.validate(ReconstructionRegistry::with_builtins()));

    auto riemann_registry = RiemannSolverRegistry::with_builtins();
    WCNS_REQUIRE(riemann_registry.names()
                 == std::vector<std::string>(
                     {"hll",
                      "hllc",
                      "roe",
                      "roe_all_speed",
                      "roe_all_speed_rieper",
                      "roe_rotated",
                      "rusanov"}));
    riemann_registry.register_solver("custom_central",
                                     [] { return std::make_unique<CustomCentralRiemann>(); });
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        riemann_registry.register_solver("custom_central", [] {
                            return std::make_unique<CustomCentralRiemann>();
                        }));
    WCNS_REQUIRE_THROWS(std::invalid_argument, riemann_registry.create("missing_riemann"));
    const auto gas = reconstruction_gas();
    const PressurePrimitiveState state {1.0, 0.2, 0.0, 0.0, 1.0};
    const RiemannSolver custom_riemann("custom_central", riemann_registry);
    const auto result = custom_riemann.solve(state, state, {1.0, 0.0, 0.0}, gas, {});
    WCNS_REQUIRE(result.requested_solver == "custom_central");
    WCNS_REQUIRE(result.used_solver == "custom_central");
    WCNS_REQUIRE(custom_riemann.summary().find("riemann_solver=custom_central;") == 0);
}

// 验收 WENO-JS/WENO-Z/MDCD-LINEAR/MDCD-HYBRID 的常数保持、尺度及镜像约定。
void test_stage_l_scalar_reconstruction_schemes()
{
    using namespace wcns;
    const auto registry = ReconstructionRegistry::with_builtins();
    const std::array<Real, 6> asymmetric {{-9.0, -3.0, 1.25, 7.5, 20.0, 80.0}};
    const auto zero_order = registry.create("zero_order");
    WCNS_REQUIRE(zero_order->stencil_requirement().point_count == 6);
    WCNS_REQUIRE_NEAR(
        zero_order->reconstruct_scalar(asymmetric, TraceSide::Left, {}), asymmetric[2], 0.0);
    WCNS_REQUIRE_NEAR(
        zero_order->reconstruct_scalar(asymmetric, TraceSide::Right, {}), asymmetric[3], 0.0);
    WCNS_REQUIRE(reconstruction_name(ReconstructionKind::ZeroOrder) == "zero_order");
    const std::array<Real, 6> constant {{2.5, 2.5, 2.5, 2.5, 2.5, 2.5}};
    for (const auto* name : {"weno_js", "weno_z", "mdcd_linear", "mdcd_hybrid"}) {
        const auto scheme = registry.create(name);
        WCNS_REQUIRE_NEAR(scheme->reconstruct_scalar(constant, TraceSide::Left, {}), 2.5, 1.0e-14);
        WCNS_REQUIRE_NEAR(scheme->reconstruct_scalar(constant, TraceSide::Right, {}), 2.5, 1.0e-14);
    }
    WCNS_REQUIRE_NEAR(mdcd_six_point_smoothness(constant, 2.5), 0.0, 0.0);

    const std::array<Real, 6> smooth {{0.7, 0.9, 1.2, 1.6, 2.1, 2.7}};
    auto reversed = smooth;
    std::reverse(reversed.begin(), reversed.end());
    for (const auto* name : {"weno_js", "weno_z", "mdcd_hybrid"}) {
        const auto scheme = registry.create(name);
        ReconstructionContext base;
        base.scale = 1.0;
        std::array<Real, 6> scaled {};
        for (std::size_t index = 0; index < smooth.size(); ++index) {
            scaled[index] = 1.0e6 * smooth[index];
        }
        ReconstructionContext large = base;
        large.scale = 1.0e6;
        const Real left = scheme->reconstruct_scalar(smooth, TraceSide::Left, base);
        const Real right = scheme->reconstruct_scalar(smooth, TraceSide::Right, base);
        WCNS_REQUIRE_NEAR(
            scheme->reconstruct_scalar(scaled, TraceSide::Left, large) / 1.0e6, left, 3.0e-15);
        WCNS_REQUIRE_NEAR(
            scheme->reconstruct_scalar(reversed, TraceSide::Left, base), right, 3.0e-15);
    }

    WcnsParameters parameters;
    const auto mdcd = mdcd_linear_reconstruct(smooth, parameters);
    const Real dispersion = parameters.mdcd_dispersion;
    const Real dissipation = parameters.mdcd_dissipation;
    const std::array<Real, 6> coefficients {{
        3.0 * (dispersion + dissipation) / 8.0,
        (-18.0 * dispersion - 30.0 * dissipation - 1.0) / 16.0,
        (12.0 * dispersion + 60.0 * dissipation + 9.0) / 16.0,
        (12.0 * dispersion - 60.0 * dissipation + 9.0) / 16.0,
        (-18.0 * dispersion + 30.0 * dissipation - 1.0) / 16.0,
        3.0 * (dispersion - dissipation) / 8.0,
    }};
    Real expected = 0.0;
    for (std::size_t index = 0; index < smooth.size(); ++index) {
        expected += coefficients[index] * smooth[index];
    }
    WCNS_REQUIRE_NEAR(mdcd.left, expected, 1.0e-15);

    const std::array<Real, 6> linear_data {{0.0, 1.0, 2.0, 3.0, 4.0, 5.0}};
    const auto hybrid_smooth = mdcd_hybrid_reconstruct_scaled(linear_data, 1.0, parameters);
    const auto linear_smooth = mdcd_linear_reconstruct(linear_data, parameters);
    WCNS_REQUIRE_NEAR(hybrid_smooth.left, linear_smooth.left, 0.0);
    WCNS_REQUIRE_NEAR(hybrid_smooth.right, linear_smooth.right, 0.0);

    const std::array<Real, 6> discontinuity {{1.0, 1.0, 1.0, 2.0, 2.0, 2.0}};
    WCNS_REQUIRE_NEAR(mdcd_six_point_smoothness(discontinuity, 1.0), 279739.0 / 5040.0, 1.0e-13);
    const auto hybrid_jump = mdcd_hybrid_reconstruct_scaled(discontinuity, 1.0, parameters);
    WCNS_REQUIRE_NEAR(hybrid_jump.left, 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(hybrid_jump.right, 2.0, 1.0e-14);

    auto invalid = parameters;
    invalid.mdcd_dissipation = invalid.mdcd_dispersion;
    WCNS_REQUIRE_THROWS(std::invalid_argument, invalid.validate());
}

// 验收同一面 Roe 特征基的 L*R=I、守恒量往返变换和特征重构回退次序。
void test_stage_l_characteristic_reconstruction()
{
    using namespace wcns;
    const auto gas = reconstruction_gas();
    const auto reference = reconstruction_reference(gas);
    const NumericalFloors floors;
    const PressurePrimitiveState left_state {1.0, 0.6, -0.2, 0.0, 1.0};
    const PressurePrimitiveState right_state {0.8, -0.1, 0.3, 0.0, 0.7};
    const Normal3 normal {0.6, 0.8, 0.0};
    const auto basis
        = make_roe_characteristic_basis(left_state, right_state, normal, gas, floors, 2);
    for (int row = 0; row < euler_components; ++row) {
        for (int column = 0; column < euler_components; ++column) {
            Real product = 0.0;
            for (int inner = 0; inner < euler_components; ++inner) {
                product
                    += basis.left[static_cast<std::size_t>(row)][static_cast<std::size_t>(inner)]
                    * basis
                          .right[static_cast<std::size_t>(inner)][static_cast<std::size_t>(column)];
            }
            WCNS_REQUIRE_NEAR(product, row == column ? 1.0 : 0.0, 2.0e-14);
        }
    }
    const auto conservative = to_conservative(left_state);
    const auto restored
        = restore_characteristic(project_characteristic(conservative, basis), basis);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(restored[static_cast<std::size_t>(component)],
                          conservative[static_cast<std::size_t>(component)],
                          3.0e-14);
    }

    const PressurePrimitiveState left_3d {1.1, 0.7, -0.4, 0.2, 0.9};
    const PressurePrimitiveState right_3d {0.75, -0.2, 0.5, -0.3, 0.65};
    const auto basis_3d = make_roe_characteristic_basis(
        left_3d, right_3d, {1.0 / 3.0, 2.0 / 3.0, 2.0 / 3.0}, gas, floors, 3);
    for (int row = 0; row < euler_components; ++row) {
        for (int column = 0; column < euler_components; ++column) {
            Real product = 0.0;
            for (int inner = 0; inner < euler_components; ++inner) {
                product
                    += basis_3d.left[static_cast<std::size_t>(row)][static_cast<std::size_t>(inner)]
                    * basis_3d
                          .right[static_cast<std::size_t>(inner)][static_cast<std::size_t>(column)];
            }
            WCNS_REQUIRE_NEAR(product, row == column ? 1.0 : 0.0, 3.0e-14);
        }
    }
    const auto conservative_3d = to_conservative(left_3d);
    const auto restored_3d
        = restore_characteristic(project_characteristic(conservative_3d, basis_3d), basis_3d);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(restored_3d[static_cast<std::size_t>(component)],
                          conservative_3d[static_cast<std::size_t>(component)],
                          4.0e-14);
    }

    Field<Real> conservative_field({4, 1, 1}, euler_components, 3);
    Field<Real> primitive_field({4, 1, 1}, euler_components, 3);
    for (int i = -3; i < 7; ++i) {
        store_state(primitive_field, {i, 0, 0}, left_state);
        store_state(conservative_field, {i, 0, 0}, conservative);
    }
    ReconstructionConfig config;
    config.scheme = "weno_z";
    config.variables = ReconstructionVariables::Characteristic;
    ReconstructionDiagnostics diagnostics;
    const auto uniform = reconstruct_thermodynamic_face(conservative_field,
                                                        primitive_field,
                                                        Axis::I,
                                                        {2, 0, 0},
                                                        config,
                                                        gas,
                                                        reference,
                                                        diagnostics,
                                                        2,
                                                        normal);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(uniform.left[static_cast<std::size_t>(component)],
                          left_state[static_cast<std::size_t>(component)],
                          3.0e-14);
        WCNS_REQUIRE_NEAR(uniform.right[static_cast<std::size_t>(component)],
                          left_state[static_cast<std::size_t>(component)],
                          3.0e-14);
    }
    WCNS_REQUIRE(diagnostics.characteristic_faces == 1);
    WCNS_REQUIRE(diagnostics.characteristic_fallbacks == 0);

    const PressurePrimitiveState valid {1.0, 0.2, 0.0, 0.0, 1.0};
    for (int i = -3; i < 7; ++i) {
        store_state(primitive_field, {i, 0, 0}, valid);
        store_state(conservative_field, {i, 0, 0}, to_conservative(valid));
    }
    for (const int i : {0, 3}) {
        auto remote = valid;
        remote[0] = 100.0;
        remote[4] = 100.0;
        store_state(primitive_field, {i, 0, 0}, remote);
        store_state(conservative_field, {i, 0, 0}, to_conservative(remote));
    }
    config.scheme = "linear5";
    diagnostics = {};
    const FaceDiagnosticLocation location {7, 2, Axis::I, {2, 0, 0}, 41, 3};
    const auto fallback = reconstruct_thermodynamic_face(conservative_field,
                                                         primitive_field,
                                                         Axis::I,
                                                         {2, 0, 0},
                                                         config,
                                                         gas,
                                                         reference,
                                                         diagnostics,
                                                         2,
                                                         {1.0, 0.0, 0.0},
                                                         location);
    WCNS_REQUIRE(fallback.left == valid);
    WCNS_REQUIRE(fallback.right == valid);
    WCNS_REQUIRE(diagnostics.characteristic_fallbacks == 1);
    WCNS_REQUIRE(diagnostics.primitive_fallbacks == 1);
    WCNS_REQUIRE(diagnostics.first_order_fallbacks == 1);
    WCNS_REQUIRE(diagnostics.fallback_events.size() == 2);
    WCNS_REQUIRE(diagnostics.fallback_events[0].location.block == 7);
    WCNS_REQUIRE(diagnostics.fallback_events[0].location.rank == 2);
    WCNS_REQUIRE(diagnostics.fallback_events[0].location.rk_stage == 3);
    WCNS_REQUIRE(diagnostics.fallback_events[0].from_strategy == "linear5:characteristic");
    WCNS_REQUIRE(diagnostics.fallback_events[0].to_strategy == "linear5:primitive");
    WCNS_REQUIRE(diagnostics.fallback_events[0].reason
                 == ReconstructionFallbackReason::InvalidCharacteristicState);
    WCNS_REQUIRE(diagnostics.fallback_events[1].to_strategy == "first_order:primitive");
}

// 验收五阶线性重构的四次多项式精确性及 WCNS-JS 的尺度不变性。
void test_stage_j_scalar_reconstruction()
{
    using namespace wcns;
    std::array<Real, 6> polynomial {};
    for (int point = 0; point < 6; ++point) {
        const Real x = static_cast<Real>(point - 2);
        polynomial[static_cast<std::size_t>(point)] = x * x * x * x;
    }
    const auto linear = linear5_reconstruct(polynomial);
    WCNS_REQUIRE_NEAR(linear.left, 0.0625, 1.0e-14);
    WCNS_REQUIRE_NEAR(linear.right, 0.0625, 1.0e-14);

    const std::array<Real, 6> smooth {{0.7, 0.9, 1.2, 1.6, 2.1, 2.7}};
    std::array<Real, 6> scaled {};
    for (std::size_t index = 0; index < smooth.size(); ++index) {
        scaled[index] = 1.0e6 * smooth[index];
    }
    const auto base = wcns5_reconstruct_scaled(smooth, 1.0);
    const auto large = wcns5_reconstruct_scaled(scaled, 1.0e6);
    WCNS_REQUIRE_NEAR(large.left / 1.0e6, base.left, 2.0e-15);
    WCNS_REQUIRE_NEAR(large.right / 1.0e6, base.right, 2.0e-15);
    WCNS_REQUIRE_THROWS(std::invalid_argument, wcns5_reconstruct_scaled(smooth, 0.0));
    ReconstructionConfig inconsistent;
    inconsistent.nonlinear.epsilon = 1.0e-5;
    WCNS_REQUIRE_THROWS(std::invalid_argument, inconsistent.validate());
}

// 验收重构出现非物理状态时确定性回退到相邻单元一阶状态并记录计数。
void test_reconstruction_positivity_fallback()
{
    using namespace wcns;
    const auto gas = reconstruction_gas();
    const auto reference = reconstruction_reference(gas);
    Field<Real> conservative({4, 1, 1}, euler_components, 3);
    Field<Real> primitive({4, 1, 1}, euler_components, 3);
    const PressurePrimitiveState valid {1.0, 0.2, 0.0, 0.0, 1.0};
    for (int i = -3; i < 7; ++i) {
        store_state(primitive, {i, 0, 0}, valid);
        store_state(conservative, {i, 0, 0}, to_conservative(valid));
    }
    for (const int i : {0, 3}) {
        auto remote = valid;
        remote[0] = 100.0;
        remote[4] = 100.0;
        store_state(primitive, {i, 0, 0}, remote);
    }
    ReconstructionConfig config;
    config.scheme = std::string(reconstruction_name(ReconstructionKind::Linear5));
    ReconstructionDiagnostics diagnostics;
    const auto result = reconstruct_thermodynamic_face(
        conservative, primitive, Axis::I, {2, 0, 0}, config, gas, reference, diagnostics, 2);
    WCNS_REQUIRE(result.left == valid);
    WCNS_REQUIRE(result.right == valid);
    WCNS_REQUIRE(diagnostics.linear_faces == 1);
    WCNS_REQUIRE(diagnostics.first_order_fallbacks == 1);
    WCNS_REQUIRE(diagnostics.fallback_events.size() == 1);
    WCNS_REQUIRE(diagnostics.fallback_events[0].from_strategy == "linear5:primitive");
    WCNS_REQUIRE(diagnostics.fallback_events[0].to_strategy == "first_order:primitive");
}

// 验收内置 Riemann 求解器的相容性、接触保持、迎风极限和法向反转对称性。
void test_stage_l_riemann_solvers()
{
    using namespace wcns;
    const auto gas = reconstruction_gas();
    const NumericalFloors floors;
    const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
    const Normal3 normal {0.6, 0.8, 0.0};
    const PressurePrimitiveState uniform {1.1, 0.7, -0.2, 0.0, 0.9};
    const auto exact = euler_flux(uniform, normal, ideal);
    for (const auto kind :
         {RiemannSolverKind::Rusanov,
          RiemannSolverKind::Hll,
          RiemannSolverKind::Hllc,
          RiemannSolverKind::Roe,
          RiemannSolverKind::AllSpeedRoe,
          RiemannSolverKind::RieperAllSpeedRoe,
          RiemannSolverKind::RotatedRoe}) {
        const RiemannSolver solver(kind);
        const auto result = solver.solve(uniform, uniform, normal, gas, floors);
        WCNS_REQUIRE(result.requested_solver == solver.name());
        WCNS_REQUIRE(result.used_solver == solver.name());
        WCNS_REQUIRE(result.fallback_reason == RiemannFallbackReason::None);
        WCNS_REQUIRE(result.spectral_radius > 0.0);
        for (int component = 0; component < euler_components; ++component) {
            WCNS_REQUIRE_NEAR(result.flux_per_unit_area[static_cast<std::size_t>(component)],
                              exact[static_cast<std::size_t>(component)],
                              2.0e-13);
        }
    }

    // The production Roe path uses the matrix-free Roe-Pike formulas.  Keep a
    // matrix reference here so future optimizations remain numerically
    // equivalent to R*|Lambda|*L*delta(U) for an oblique three-dimensional face.
    const PressurePrimitiveState matrix_left {1.15, 0.8, -0.35, 0.22, 1.1};
    const PressurePrimitiveState matrix_right {0.72, -0.18, 0.41, -0.09, 0.63};
    const Normal3 matrix_normal {2.0 / 3.0, -1.0 / 3.0, 2.0 / 3.0};
    const auto matrix_basis = make_roe_characteristic_basis(
        matrix_left, matrix_right, matrix_normal, gas, floors, 3);
    const auto matrix_left_conservative = to_conservative(matrix_left, ideal);
    const auto matrix_right_conservative = to_conservative(matrix_right, ideal);
    ConservativeState matrix_jump {};
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        matrix_jump[index] = matrix_right_conservative[index]
            - matrix_left_conservative[index];
    }
    auto matrix_strengths = project_characteristic(matrix_jump, matrix_basis);
    const Real matrix_roe_un = matrix_basis.right[1][1];
    const Real matrix_roe_sound
        = 0.5 * (matrix_basis.right[1][4] - matrix_basis.right[1][0]);
    const Real entropy_delta = 0.1
        * std::max(
              {sound_speed(matrix_left, ideal), sound_speed(matrix_right, ideal), matrix_roe_sound});
    const std::array<Real, euler_components> matrix_eigenvalues {{
        matrix_roe_un - matrix_roe_sound,
        matrix_roe_un,
        matrix_roe_un,
        matrix_roe_un,
        matrix_roe_un + matrix_roe_sound}};
    for (int wave = 0; wave < euler_components; ++wave) {
        Real magnitude = std::abs(matrix_eigenvalues[static_cast<std::size_t>(wave)]);
        if (magnitude < entropy_delta) {
            magnitude = 0.5 * (magnitude * magnitude / entropy_delta + entropy_delta);
        }
        matrix_strengths[static_cast<std::size_t>(wave)] *= magnitude;
    }
    const auto matrix_dissipation = restore_characteristic(matrix_strengths, matrix_basis);
    const auto matrix_left_flux = euler_flux(matrix_left, matrix_normal, ideal);
    const auto matrix_right_flux = euler_flux(matrix_right, matrix_normal, ideal);
    const auto matrix_free_result = RiemannSolver(RiemannSolverKind::Roe).solve(
        matrix_left, matrix_right, matrix_normal, gas, floors);
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        const Real expected
            = 0.5 * (matrix_left_flux[index] + matrix_right_flux[index])
            - 0.5 * matrix_dissipation[index];
        WCNS_REQUIRE_NEAR(matrix_free_result.flux_per_unit_area[index], expected, 3.0e-13);
    }

    const PressurePrimitiveState contact_left {1.0, 0.0, 0.0, 0.0, 1.0};
    const PressurePrimitiveState contact_right {2.0, 0.0, 0.0, 0.0, 1.0};
    const ConservativeState contact_flux {{0.0, 1.0, 0.0, 0.0, 0.0}};
    const auto hllc_contact = RiemannSolver(RiemannSolverKind::Hllc)
                                  .flux(contact_left, contact_right, {1.0, 0.0, 0.0}, gas, floors);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(hllc_contact[static_cast<std::size_t>(component)],
                          contact_flux[static_cast<std::size_t>(component)],
                          3.0e-13);
    }
    const auto rusanov_contact
        = RiemannSolver(RiemannSolverKind::Rusanov)
              .flux(contact_left, contact_right, {1.0, 0.0, 0.0}, gas, floors);
    WCNS_REQUIRE(std::abs(hllc_contact[0]) < std::abs(rusanov_contact[0]));

    RiemannSolverParameters strong_entropy_fix;
    strong_entropy_fix.entropy_fix_coefficient = 0.2;
    const auto roe_contact = RiemannSolver(RiemannSolverKind::Roe, strong_entropy_fix)
                                 .flux(contact_left, contact_right, {1.0, 0.0, 0.0}, gas, floors);
    WCNS_REQUIRE(roe_contact[0] < 0.0);

    RiemannSolverParameters preconditioned_parameters;
    preconditioned_parameters.weiss_smith = true;
    preconditioned_parameters.preconditioner.mach_cutoff = 1.0e-3;
    preconditioned_parameters.preconditioner.viscous_cutoff = 1.0;
    const PressurePrimitiveState low_mach_left {1.0, 0.01, 0.0, 0.0, 1.0};
    const PressurePrimitiveState low_mach_right {1.0, 0.011, 0.0, 0.0, 1.000001};
    const RiemannSolver low_mach_roe(RiemannSolverKind::Roe, preconditioned_parameters);
    const auto low_mach_result = low_mach_roe.solve(
        low_mach_left, low_mach_right, {1.0, 0.0, 0.0}, gas, floors);
    const auto ordinary_result = RiemannSolver(RiemannSolverKind::Roe).solve(
        low_mach_left, low_mach_right, {1.0, 0.0, 0.0}, gas, floors);
    WCNS_REQUIRE(low_mach_result.used_solver == "roe");
    WCNS_REQUIRE(low_mach_result.spectral_radius < ordinary_result.spectral_radius);
    WCNS_REQUIRE(low_mach_roe.summary().find("preconditioner=weiss_smith")
                 != std::string::npos);
    const auto low_mach_reverse = low_mach_roe.solve(
        low_mach_right, low_mach_left, {-1.0, 0.0, 0.0}, gas, floors);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            low_mach_result.flux_per_unit_area[static_cast<std::size_t>(component)],
            -low_mach_reverse.flux_per_unit_area[static_cast<std::size_t>(component)],
            2.0e-11);
    }
    const auto low_mach_uniform = low_mach_roe.solve(
        low_mach_left, low_mach_left, {1.0, 0.0, 0.0}, gas, floors);
    const auto low_mach_uniform_exact
        = euler_flux(low_mach_left, {1.0, 0.0, 0.0}, ideal);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            low_mach_uniform.flux_per_unit_area[static_cast<std::size_t>(component)],
            low_mach_uniform_exact[static_cast<std::size_t>(component)],
            5.0e-12);
    }
    Real scaled_pressure_reference = 0.0;
    for (const Real mach : {1.0e-1, 1.0e-2, 1.0e-3}) {
        constexpr Real pressure_jump = 0.01;
        const Real background_pressure = 1.0 / (gas.gamma() * mach * mach);
        const PressurePrimitiveState asymptotic_left {
            {1.0, 1.0, 0.0, 0.0, background_pressure - pressure_jump}};
        const PressurePrimitiveState asymptotic_right {
            {1.0, 1.0, 0.0, 0.0, background_pressure + pressure_jump}};
        const auto result = low_mach_roe.solve(
            asymptotic_left, asymptotic_right, {1.0, 0.0, 0.0}, gas, floors);
        WCNS_REQUIRE(result.used_solver == "roe");
        WCNS_REQUIRE(result.spectral_radius > 1.0);
        WCNS_REQUIRE(result.spectral_radius < 2.0);
        for (const Real value : result.flux_per_unit_area) {
            WCNS_REQUIRE(std::isfinite(value));
        }
        const Real scaled_pressure
            = ((asymptotic_right[4] - asymptotic_left[4]) / background_pressure)
            / (mach * mach);
        if (scaled_pressure_reference == 0.0) scaled_pressure_reference = scaled_pressure;
        WCNS_REQUIRE_NEAR(scaled_pressure, scaled_pressure_reference, 2.0e-8);
    }

    const PressurePrimitiveState supersonic_left {1.0, 4.0, 0.1, 0.0, 1.0};
    const PressurePrimitiveState supersonic_right {0.8, 3.5, -0.2, 0.0, 0.7};
    const auto supersonic_exact = euler_flux(supersonic_left, {1.0, 0.0, 0.0}, ideal);
    const auto supersonic_hllc
        = RiemannSolver(RiemannSolverKind::Hllc)
              .flux(supersonic_left, supersonic_right, {1.0, 0.0, 0.0}, gas, floors);
    const auto supersonic_hll = RiemannSolver(RiemannSolverKind::Hll)
                                    .flux(supersonic_left,
                                          supersonic_right,
                                          {1.0, 0.0, 0.0},
                                          gas,
                                          floors);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(supersonic_hllc[static_cast<std::size_t>(component)],
                          supersonic_exact[static_cast<std::size_t>(component)],
                          0.0);
        WCNS_REQUIRE_NEAR(supersonic_hll[static_cast<std::size_t>(component)],
                          supersonic_exact[static_cast<std::size_t>(component)],
                          0.0);
    }
    const auto ordinary_roe = RiemannSolver(RiemannSolverKind::Roe).solve(
        supersonic_left, supersonic_right, {1.0, 0.0, 0.0}, gas, floors);
    const auto recovered_roe = low_mach_roe.solve(
        supersonic_left, supersonic_right, {1.0, 0.0, 0.0}, gas, floors);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(recovered_roe.flux_per_unit_area[static_cast<std::size_t>(component)],
                          ordinary_roe.flux_per_unit_area[static_cast<std::size_t>(component)],
                          5.0e-13);
    }
    WCNS_REQUIRE_NEAR(recovered_roe.spectral_radius,
                      ordinary_roe.spectral_radius,
                      5.0e-13);

    const auto all_speed_low_mach
        = RiemannSolver(RiemannSolverKind::RieperAllSpeedRoe).solve(
        low_mach_left, low_mach_right, {1.0, 0.0, 0.0}, gas, floors);
    WCNS_REQUIRE(all_speed_low_mach.requested_solver == "roe_all_speed_rieper");
    WCNS_REQUIRE(all_speed_low_mach.used_solver == "roe_all_speed_rieper");
    WCNS_REQUIRE_NEAR(all_speed_low_mach.spectral_radius,
                      ordinary_result.spectral_radius,
                      5.0e-13);
    Real low_mach_flux_difference = 0.0;
    for (int component = 0; component < euler_components; ++component) {
        low_mach_flux_difference
            = std::max(low_mach_flux_difference,
                       std::abs(all_speed_low_mach.flux_per_unit_area[static_cast<std::size_t>(
                                    component)]
                                - ordinary_result.flux_per_unit_area[static_cast<std::size_t>(
                                    component)]));
    }
    WCNS_REQUIRE(low_mach_flux_difference > 1.0e-8);

    const PressurePrimitiveState velocity_jump_right {1.0, 0.011, 0.0, 0.0, 1.0};
    const auto ordinary_velocity_jump = RiemannSolver(RiemannSolverKind::Roe).solve(
        low_mach_left, velocity_jump_right, {1.0, 0.0, 0.0}, gas, floors);
    const auto all_speed_velocity_jump
        = RiemannSolver(RiemannSolverKind::RieperAllSpeedRoe)
              .solve(low_mach_left,
                     velocity_jump_right,
                     {1.0, 0.0, 0.0},
                     gas,
                     floors);
    const auto velocity_jump_basis = make_roe_characteristic_basis(
        low_mach_left, velocity_jump_right, {1.0, 0.0, 0.0}, gas, floors, 3);
    const auto left_velocity_flux
        = euler_flux(low_mach_left, {1.0, 0.0, 0.0}, ideal);
    const auto right_velocity_flux
        = euler_flux(velocity_jump_right, {1.0, 0.0, 0.0}, ideal);
    const auto dissipation_coefficients = [&](const RiemannResult& result) {
        ConservativeState dissipation {};
        for (int component = 0; component < euler_components; ++component) {
            const auto index = static_cast<std::size_t>(component);
            dissipation[index] = left_velocity_flux[index] + right_velocity_flux[index]
                - 2.0 * result.flux_per_unit_area[index];
        }
        return project_characteristic(dissipation, velocity_jump_basis);
    };
    const auto ordinary_dissipation = dissipation_coefficients(ordinary_velocity_jump);
    const auto all_speed_dissipation = dissipation_coefficients(all_speed_velocity_jump);
    for (const auto wave : {0U, 4U}) {
        const Real ratio
            = std::abs(all_speed_dissipation[wave] / ordinary_dissipation[wave]);
        WCNS_REQUIRE(ratio > 0.0);
        WCNS_REQUIRE(ratio < 0.02);
    }

    const auto all_speed_supersonic
        = RiemannSolver(RiemannSolverKind::RieperAllSpeedRoe)
              .solve(supersonic_left,
                     supersonic_right,
                     {1.0, 0.0, 0.0},
                     gas,
                     floors);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(
            all_speed_supersonic.flux_per_unit_area[static_cast<std::size_t>(component)],
            ordinary_roe.flux_per_unit_area[static_cast<std::size_t>(component)],
            5.0e-13);
    }
    WCNS_REQUIRE_NEAR(all_speed_supersonic.spectral_radius,
                      ordinary_roe.spectral_radius,
                      5.0e-13);

    RiemannSolverParameters li_gu_parameters;
    li_gu_parameters.all_speed.reference_mach = 0.1;
    li_gu_parameters.all_speed.dissipation_scale = 0.02;
    li_gu_parameters.all_speed.pressure_coefficient = 0.05;
    const RiemannSolver li_gu(RiemannSolverKind::AllSpeedRoe, li_gu_parameters);
    const auto li_gu_low_mach = li_gu.solve(
        low_mach_left, low_mach_right, {1.0, 0.0, 0.0}, gas, floors);
    WCNS_REQUIRE(li_gu_low_mach.requested_solver == "roe_all_speed");
    WCNS_REQUIRE(li_gu_low_mach.used_solver == "roe_all_speed");
    WCNS_REQUIRE(li_gu_low_mach.spectral_radius < ordinary_result.spectral_radius);
    WCNS_REQUIRE(li_gu.summary().find("all_speed_reference_mach=0.10000000000000001")
                 != std::string::npos);

    const Real root_left = std::sqrt(low_mach_left[0]);
    const Real root_right = std::sqrt(low_mach_right[0]);
    const Real root_sum = root_left + root_right;
    const Real roe_un
        = (root_left * low_mach_left[1] + root_right * low_mach_right[1]) / root_sum;
    const auto low_mach_left_conservative = to_conservative(low_mach_left, ideal);
    const auto low_mach_right_conservative = to_conservative(low_mach_right, ideal);
    const Real left_enthalpy
        = (low_mach_left_conservative[4] + low_mach_left[4]) / low_mach_left[0];
    const Real right_enthalpy
        = (low_mach_right_conservative[4] + low_mach_right[4]) / low_mach_right[0];
    const Real roe_enthalpy
        = (root_left * left_enthalpy + root_right * right_enthalpy) / root_sum;
    const Real roe_sound = std::sqrt(
        (gas.gamma() - 1.0) * (roe_enthalpy - 0.5 * roe_un * roe_un));
    const Real mach = std::max(std::abs(roe_un) / roe_sound, 0.1);
    const Real mach_squared = mach * mach;
    const Real factor = mach
        * std::sqrt((4.0 + (1.0 - mach_squared) * (1.0 - mach_squared))
                    / (1.0 + mach_squared));
    WCNS_REQUIRE_NEAR(li_gu_low_mach.spectral_radius,
                      std::abs(roe_un) + factor * roe_sound,
                      2.0e-13);

    auto no_pressure_parameters = li_gu_parameters;
    no_pressure_parameters.all_speed.pressure_coefficient = 0.0;
    const auto no_pressure = RiemannSolver(RiemannSolverKind::AllSpeedRoe,
                                           no_pressure_parameters)
                                 .solve(low_mach_left,
                                        low_mach_right,
                                        {1.0, 0.0, 0.0},
                                        gas,
                                        floors);
    ConservativeState averaged_conservative {};
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        averaged_conservative[index] = 0.5
            * (low_mach_left_conservative[index] + low_mach_right_conservative[index]);
    }
    auto pressure_vector = averaged_conservative;
    pressure_vector[4] += to_primitive(averaged_conservative, ideal)[4];
    const Real expected_pressure_coefficient = (1.0 - factor) * 0.05
        * (low_mach_right[4] - low_mach_left[4])
        / (0.1 * std::sqrt(low_mach_left[0] * low_mach_right[0]) * roe_sound);
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        WCNS_REQUIRE_NEAR(li_gu_low_mach.flux_per_unit_area[index]
                              - no_pressure.flux_per_unit_area[index],
                          -expected_pressure_coefficient * pressure_vector[index],
                          2.0e-13);
    }

    auto high_mach_recovery_parameters = li_gu_parameters;
    high_mach_recovery_parameters.all_speed.dissipation_scale = 0.5;
    const auto li_gu_supersonic
        = RiemannSolver(RiemannSolverKind::AllSpeedRoe, high_mach_recovery_parameters)
              .solve(supersonic_left,
                     supersonic_right,
                     {1.0, 0.0, 0.0},
                     gas,
                     floors);
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        WCNS_REQUIRE_NEAR(li_gu_supersonic.flux_per_unit_area[index],
                          ordinary_roe.flux_per_unit_area[index],
                          5.0e-13);
    }

    const PressurePrimitiveState aligned_left {1.0, 4.0, 0.2, -0.1, 1.0};
    const PressurePrimitiveState aligned_right {0.9, 3.7, 0.2, -0.1, 0.8};
    const auto aligned_roe = RiemannSolver(RiemannSolverKind::Roe).solve(
        aligned_left, aligned_right, {1.0, 0.0, 0.0}, gas, floors);
    const auto aligned_rotated = RiemannSolver(RiemannSolverKind::RotatedRoe).solve(
        aligned_left, aligned_right, {1.0, 0.0, 0.0}, gas, floors);
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        WCNS_REQUIRE_NEAR(aligned_rotated.flux_per_unit_area[index],
                          aligned_roe.flux_per_unit_area[index],
                          2.0e-12);
    }

    const PressurePrimitiveState left {1.0, 0.9, -0.3, 0.1, 1.2};
    const PressurePrimitiveState right {0.7, -0.2, 0.4, -0.1, 0.6};
    const Normal3 reverse_normal {-normal.x, -normal.y, -normal.z};
    for (const auto kind :
         {RiemannSolverKind::Rusanov,
          RiemannSolverKind::Hll,
          RiemannSolverKind::Hllc,
          RiemannSolverKind::Roe,
          RiemannSolverKind::AllSpeedRoe,
          RiemannSolverKind::RieperAllSpeedRoe,
          RiemannSolverKind::RotatedRoe}) {
        const RiemannSolver solver(kind);
        const auto forward = solver.flux(left, right, normal, gas, floors);
        const auto reverse = solver.flux(right, left, reverse_normal, gas, floors);
        for (int component = 0; component < euler_components; ++component) {
            WCNS_REQUIRE_NEAR(forward[static_cast<std::size_t>(component)],
                              -reverse[static_cast<std::size_t>(component)],
                              5.0e-13);
        }
    }

    RiemannConfig config;
    config.scheme = "hllc";
    config.validate(RiemannSolverRegistry::with_builtins(config.parameters));
    WCNS_REQUIRE(config.restart_signature().find("riemann_config_v2;") == 0);
    config.scheme = "missing";
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        config.validate(RiemannSolverRegistry::with_builtins(config.parameters)));

    RiemannSolverParameters fallback_parameters;
    fallback_parameters.denominator_tolerance = 0.5;
    const PressurePrimitiveState near_floor {1.0e-8, 0.0, 0.0, 0.0, 1.0e-8};
    const auto fallback_result = RiemannSolver(RiemannSolverKind::Hllc, fallback_parameters)
                                     .solve(near_floor, near_floor, {1.0, 0.0, 0.0}, gas, floors);
    WCNS_REQUIRE(fallback_result.requested_solver == "hllc");
    WCNS_REQUIRE(fallback_result.used_solver == "rusanov");
    WCNS_REQUIRE(fallback_result.fallback_reason == RiemannFallbackReason::InvalidWaveSpeed);
    RiemannDiagnostics diagnostics;
    const FaceDiagnosticLocation location {9, 3, Axis::J, {4, 5, 0}, 17, 2};
    diagnostics.record(fallback_result, location);
    WCNS_REQUIRE(diagnostics.total_faces == 1);
    WCNS_REQUIRE(diagnostics.requested_faces.at("hllc") == 1);
    WCNS_REQUIRE(diagnostics.used_faces.at("rusanov") == 1);
    WCNS_REQUIRE(diagnostics.fallback_count() == 1);
    WCNS_REQUIRE(diagnostics.fallback_count(RiemannFallbackReason::InvalidWaveSpeed) == 1);
    WCNS_REQUIRE(diagnostics.fallback_events[0].location.axis == Axis::J);
    WCNS_REQUIRE(diagnostics.fallback_events[0].location.face.i == 4);
    WCNS_REQUIRE(diagnostics.fallback_events[0].location.face.j == 5);
    WCNS_REQUIRE(diagnostics.fallback_events[0].location.face.k == 0);
    WCNS_REQUIRE(riemann_fallback_reason_name(diagnostics.fallback_events[0].reason)
                 == "invalid_wave_speed");

    RiemannResult two_level;
    two_level.requested_solver = "roe";
    two_level.used_solver = "rusanov";
    two_level.fallback_reason = RiemannFallbackReason::InvalidRoeAverage;
    two_level.fallback_path = {
        {"roe", "hllc", RiemannFallbackReason::InvalidRoeAverage},
        {"hllc", "rusanov", RiemannFallbackReason::InvalidIntermediateState},
    };
    RiemannDiagnostics two_level_diagnostics;
    two_level_diagnostics.record(two_level, location);
    WCNS_REQUIRE(two_level_diagnostics.fallback_count() == 2);
    WCNS_REQUIRE(two_level_diagnostics.fallback_count(RiemannFallbackReason::InvalidRoeAverage)
                 == 1);
    WCNS_REQUIRE(
        two_level_diagnostics.fallback_count(RiemannFallbackReason::InvalidIntermediateState) == 1);
    WCNS_REQUIRE(two_level_diagnostics.fallback_events[0].from_solver == "roe");
    WCNS_REQUIRE(two_level_diagnostics.fallback_events[0].to_solver == "hllc");
    WCNS_REQUIRE(two_level_diagnostics.fallback_events[1].from_solver == "hllc");
    WCNS_REQUIRE(two_level_diagnostics.fallback_events[1].to_solver == "rusanov");
}

// 验收阶段 L 的光滑特征收敛、激波/高 Mach 健壮性和非有限模板失败语义。
void test_stage_l_algorithm_benchmarks()
{
    using namespace wcns;
    const auto gas = reconstruction_gas();
    const auto reference = reconstruction_reference(gas);
    const NumericalFloors floors;
    const auto smooth_error = [&](const char* scheme, Real spacing) {
        Field<Real> conservative({4, 1, 1}, euler_components, 3);
        Field<Real> primitive({4, 1, 1}, euler_components, 3);
        const auto analytic = [](Real x) {
            return PressurePrimitiveState {
                1.0 + 0.1 * std::exp(0.4 * x),
                0.25 + 0.04 * std::exp(0.3 * x),
                -0.1 + 0.02 * std::exp(0.2 * x),
                0.0,
                1.0 + 0.08 * std::exp(0.5 * x),
            };
        };
        for (int i = -3; i < 7; ++i) {
            const Real x = (static_cast<Real>(i) - 1.5) * spacing;
            const auto state = analytic(x);
            store_state(primitive, {i, 0, 0}, state);
            store_state(conservative, {i, 0, 0}, to_conservative(state));
        }
        ReconstructionConfig config;
        config.scheme = scheme;
        config.variables = ReconstructionVariables::Characteristic;
        ReconstructionDiagnostics diagnostics;
        const auto reconstructed = reconstruct_thermodynamic_face(conservative,
                                                                  primitive,
                                                                  Axis::I,
                                                                  {2, 0, 0},
                                                                  config,
                                                                  gas,
                                                                  reference,
                                                                  diagnostics,
                                                                  2,
                                                                  {1.0, 0.0, 0.0});
        WCNS_REQUIRE(diagnostics.fallback_events.empty());
        const auto exact = analytic(0.0);
        Real error = 0.0;
        for (int component = 0; component < euler_components; ++component) {
            const auto index = static_cast<std::size_t>(component);
            error = std::max(error, std::abs(reconstructed.left[index] - exact[index]));
            error = std::max(error, std::abs(reconstructed.right[index] - exact[index]));
        }
        return error;
    };
    for (const auto* scheme : {"weno_js", "weno_z", "mdcd_linear", "mdcd_hybrid"}) {
        const Real coarse = smooth_error(scheme, 0.4);
        const Real fine = smooth_error(scheme, 0.2);
        WCNS_REQUIRE(coarse > fine);
        WCNS_REQUIRE(coarse / fine > 12.0);
    }

    auto registry = ReconstructionRegistry::with_builtins();
    std::array<Real, 6> non_finite {{0.0, 1.0, 2.0, 3.0, 4.0, 5.0}};
    non_finite[2] = std::numeric_limits<Real>::infinity();
    for (const auto* scheme : {"weno_js", "weno_z", "mdcd_linear", "mdcd_hybrid"}) {
        WCNS_REQUIRE_THROWS(
            PhysicsError,
            registry.create(scheme)->reconstruct_scalar(non_finite, TraceSide::Left, {}));
    }
    for (int sample = 0; sample < 50; ++sample) {
        std::array<Real, 6> values {};
        for (std::size_t point = 0; point < values.size(); ++point) {
            values[point]
                = std::sin(0.37 * static_cast<Real>((sample + 1) * static_cast<int>(point + 1)));
        }
        WCNS_REQUIRE(mdcd_six_point_smoothness(values, 1.0) >= 0.0);
    }

    const PressurePrimitiveState sod_left {1.0, 0.0, 0.0, 0.0, 1.0};
    const PressurePrimitiveState sod_right {0.125, 0.0, 0.0, 0.0, 0.1};
    const PressurePrimitiveState high_mach_left {1.0, 50.0, 0.0, 0.0, 1.0};
    const PressurePrimitiveState high_mach_right {0.8, 45.0, 0.0, 0.0, 0.7};
    for (const auto kind : {RiemannSolverKind::Hll,
                            RiemannSolverKind::Hllc,
                            RiemannSolverKind::Roe,
                            RiemannSolverKind::AllSpeedRoe,
                            RiemannSolverKind::RieperAllSpeedRoe,
                            RiemannSolverKind::RotatedRoe}) {
        const RiemannSolver solver(kind);
        for (const auto result :
             {solver.solve(sod_left, sod_right, {1.0, 0.0, 0.0}, gas, floors),
              solver.solve(high_mach_left, high_mach_right, {1.0, 0.0, 0.0}, gas, floors)}) {
            WCNS_REQUIRE(result.fallback_reason == RiemannFallbackReason::None);
            WCNS_REQUIRE(result.spectral_radius > 0.0);
            for (const auto component : result.flux_per_unit_area) {
                WCNS_REQUIRE(std::isfinite(component));
            }
        }
    }

    const IdealGas ideal {gas.gamma(), floors.density, floors.pressure};
    const PressurePrimitiveState wave_left {1.0, 0.5, 0.0, 0.0, 1.0};
    const auto wave_basis
        = make_roe_characteristic_basis(wave_left, wave_left, {1.0, 0.0, 0.0}, gas, floors, 3);
    ConservativeState characteristic_jump {};
    characteristic_jump[4] = 1.0e-6;
    const auto conservative_left = to_conservative(wave_left, ideal);
    const auto physical_jump = restore_characteristic(characteristic_jump, wave_basis);
    ConservativeState conservative_right = conservative_left;
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        conservative_right[index] += physical_jump[index];
    }
    const auto wave_right = to_primitive(conservative_right, ideal);
    const auto roe_wave = RiemannSolver(RiemannSolverKind::Roe)
                              .solve(wave_left, wave_right, {1.0, 0.0, 0.0}, gas, floors);
    const auto wave_left_flux = euler_flux(wave_left, {1.0, 0.0, 0.0}, ideal);
    const auto wave_right_flux = euler_flux(wave_right, {1.0, 0.0, 0.0}, ideal);
    const Real acoustic_speed = 0.5 + std::sqrt(gas.gamma());
    for (int component = 0; component < euler_components; ++component) {
        const auto index = static_cast<std::size_t>(component);
        const Real expected = 0.5 * (wave_left_flux[index] + wave_right_flux[index])
            - 0.5 * acoustic_speed * physical_jump[index];
        WCNS_REQUIRE_NEAR(roe_wave.flux_per_unit_area[index], expected, 2.0e-12);
    }
}

// 验收默认 Rusanov 只接收单位法向并保持阶段 J 的解析 Euler 通量基线。
void test_rusanov_riemann_solver()
{
    using namespace wcns;
    const auto gas = reconstruction_gas();
    const NumericalFloors floors;
    const PressurePrimitiveState state {1.1, 0.7, -0.2, 0.0, 0.9};
    const Normal3 normal {0.6, 0.8, 0.0};
    const RiemannSolver solver;
    const auto flux = solver.flux(state, state, normal, gas, floors);
    const auto exact = euler_flux(state, normal, {gas.gamma(), floors.density, floors.pressure});
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(flux[static_cast<std::size_t>(component)],
                          exact[static_cast<std::size_t>(component)],
                          1.0e-14);
    }
    WCNS_REQUIRE(solver.summary().find("riemann_solver=rusanov;") == 0);
    WCNS_REQUIRE_THROWS(PhysicsError, solver.flux(state, state, {2.0, 0.0, 0.0}, gas, floors));
}
