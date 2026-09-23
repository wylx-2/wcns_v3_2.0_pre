#include "test_support.hpp"

#include <wcns/solver/implicit_time_integrator.hpp>
#include <wcns/solver/low_mach_preconditioner.hpp>

#include <cmath>

namespace {

wcns::Real bdf2_decay_error(wcns::Real time_step)
{
    using namespace wcns;
    constexpr Real decay = 0.7;
    constexpr Real final_time = 1.0;
    const int step_count = static_cast<int>(std::llround(final_time / time_step));
    time_step = final_time / static_cast<Real>(step_count);
    Real previous = 1.0;
    // Use an exact first layer so this test isolates the BDF2 recurrence rather
    // than measuring the deliberately first-order BDF1 startup policy.
    Real current = std::exp(-decay * time_step);
    for (int step = 1; step < step_count; ++step) {
        const Real next = (4.0 * current - previous) / (3.0 + 2.0 * decay * time_step);
        previous = current;
        current = next;
    }
    return std::abs(current - std::exp(-decay * final_time));
}

} // namespace

// 验收 LU-SGS 扫掠、BDF 离散和 Weiss-Smith 矩阵的核心代数性质。
void test_implicit_time_integrator()
{
    using namespace wcns;

    const Extent3 extent {3, 1, 1};
    Field<Real> rhs(extent, 1, 0, 0.0);
    Field<Real> diagonal(extent, 1, 0, 4.0);
    Field<Real> coupling(extent, 4, 0, 0.0);
    for (int i = 0; i < extent.ni; ++i) {
        rhs(i, 0, 0, 0) = static_cast<Real>(i + 1);
        coupling(i, 0, 0, 0) = 1.0;
        coupling(i, 0, 0, 1) = 1.0;
    }
    const auto increment = solve_scalar_lu_sgs(rhs, diagonal, coupling, 2);
    WCNS_REQUIRE_NEAR(increment(0, 0, 0, 0), 0.4462890625, 1.0e-15);
    WCNS_REQUIRE_NEAR(increment(1, 0, 0, 0), 0.78515625, 1.0e-15);
    WCNS_REQUIRE_NEAR(increment(2, 0, 0, 0), 0.890625, 1.0e-15);

    auto residual_norm = [&](const Field<Real>& solution) {
        Real sum = 0.0;
        for (int i = 0; i < extent.ni; ++i) {
            Real residual = rhs(i, 0, 0, 0) - diagonal(i, 0, 0, 0)
                * solution(i, 0, 0, 0);
            if (i > 0) residual += coupling(i, 0, 0, 0) * solution(i - 1, 0, 0, 0);
            if (i + 1 < extent.ni) {
                residual += coupling(i, 0, 0, 1) * solution(i + 1, 0, 0, 0);
            }
            sum += residual * residual;
        }
        return std::sqrt(sum);
    };
    LuSgsIterationConfig two_sweeps;
    two_sweeps.sweeps = 2;
    const auto twice = solve_scalar_lu_sgs(rhs, diagonal, coupling, 2, two_sweeps);
    WCNS_REQUIRE(residual_norm(twice) < residual_norm(increment));
    LuSgsIterationConfig invalid_sweeps;
    invalid_sweeps.sweeps = 0;
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        solve_scalar_lu_sgs(rhs, diagonal, coupling, 2, invalid_sweeps));
    invalid_sweeps.sweeps = 5;
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        solve_scalar_lu_sgs(rhs, diagonal, coupling, 2, invalid_sweeps));
    const auto face_scalar_increment
        = solve_face_block_lu_sgs(rhs, diagonal, coupling, 2);
    for (int i = 0; i < extent.ni; ++i) {
        WCNS_REQUIRE_NEAR(face_scalar_increment(i, 0, 0, 0),
                          increment(i, 0, 0, 0),
                          1.0e-15);
    }

    Field<Real> zero_rhs(extent, 2, 0, 0.0);
    const auto zero_increment = solve_scalar_lu_sgs(zero_rhs, diagonal, coupling, 2);
    for (int i = 0; i < extent.ni; ++i) {
        for (int variable = 0; variable < zero_increment.components(); ++variable) {
            WCNS_REQUIRE_NEAR(zero_increment(i, 0, 0, variable), 0.0, 1.0e-14);
        }
    }

    StructuredBlock spectral_block(17, "implicit-spectral", 0, 2, 2, {9, 9, 1}, 3);
    for (int j = 0; j < spectral_block.vertex_extent().nj; ++j) {
        for (int i = 0; i < spectral_block.vertex_extent().ni; ++i) {
            spectral_block.coordinates.x(i, j, 0) = 0.125 * i;
            spectral_block.coordinates.y(i, j, 0) = 0.125 * j;
            spectral_block.coordinates.z(i, j, 0) = 0.0;
        }
    }
    const auto metric = initialize_metric_field(
                            spectral_block,
                            ProfileFactory::create(AlgorithmProfileKind::PhengleiWcns))
                            .metric;
    for (int j = 0; j < spectral_block.cell_extent().nj; ++j) {
        for (int i = 0; i < spectral_block.cell_extent().ni; ++i) {
            spectral_block.flow.temperature_primitive(i, j, 0, temperature_density) = 1.0;
            spectral_block.flow.temperature_primitive(i, j, 0, temperature_velocity_x) = 0.2;
            spectral_block.flow.temperature_primitive(i, j, 0, temperature_velocity_y) = 0.0;
            spectral_block.flow.temperature_primitive(i, j, 0, temperature_velocity_z) = 0.0;
            spectral_block.flow.temperature_primitive(i, j, 0, temperature_value) = 1.0;
        }
    }
    const auto spectral_gas = GasModel::from_input(GasModelInput {1.4, {}, 287.0});
    const auto spectral_reference = ReferenceScales::derive(
        ReferenceInput {34.0, 1.2, 300.0, 1.0, 1.8e-5, {}, {}}, spectral_gas);
    Field<Real> face_viscous_coupling(spectral_block.cell_extent(), 4, 0, 0.0);
    for (int component = 0; component < face_viscous_coupling.components(); ++component) {
        face_viscous_coupling(3, 3, 0, component) = static_cast<Real>(component + 1);
    }
    const auto inviscid_system = build_scalar_spectral_system(spectral_block,
                                                               metric,
                                                               spectral_gas,
                                                               spectral_reference,
                                                               NumericalFloors {},
                                                               1.0);
    WCNS_REQUIRE_NEAR(inviscid_system.coupling(3, 3, 0, 0)
                          - inviscid_system.coupling(3, 3, 0, 1),
                      1.6,
                      1.0e-12);
    const auto viscous_system = build_scalar_spectral_system(spectral_block,
                                                              metric,
                                                              spectral_gas,
                                                              spectral_reference,
                                                              NumericalFloors {},
                                                              1.0,
                                                              0.0,
                                                              nullptr,
                                                              0.0,
                                                              &face_viscous_coupling);
    for (int component = 0; component < face_viscous_coupling.components(); ++component) {
        WCNS_REQUIRE_NEAR(viscous_system.coupling(3, 3, 0, component)
                              - inviscid_system.coupling(3, 3, 0, component),
                          static_cast<Real>(component + 1),
                          1.0e-13);
    }
    WCNS_REQUIRE_NEAR(viscous_system.diagonal(3, 3, 0, 0)
                          - inviscid_system.diagonal(3, 3, 0, 0),
                      20.0,
                      1.0e-12);
    const auto euler_face_blocks = build_euler_face_coupling_blocks(
        spectral_block,
        metric,
        spectral_gas,
        spectral_reference,
        NumericalFloors {},
        &face_viscous_coupling);
    constexpr int entries_per_face = euler_components * euler_components;
    WCNS_REQUIRE_NEAR(euler_face_blocks(3, 3, 0, 0 * entries_per_face + 1)
                          - euler_face_blocks(3, 3, 0, 1 * entries_per_face + 1),
                      8.0,
                      1.0e-12);
    WCNS_REQUIRE_NEAR(euler_face_blocks(3, 3, 0, 0)
                          - euler_face_blocks(3, 3, 0, 1 * entries_per_face),
                      -1.0,
                      1.0e-12);
    Field<Real> invalid_face_coupling(spectral_block.cell_extent(), 1, 0, 0.0);
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        build_scalar_spectral_system(spectral_block,
                                                     metric,
                                                     spectral_gas,
                                                     spectral_reference,
                                                     NumericalFloors {},
                                                     1.0,
                                                     0.0,
                                                     nullptr,
                                                     0.0,
                                                     &invalid_face_coupling));

    const Extent3 block_extent {1, 1, 1};
    Field<Real> block_rhs(block_extent, 2, 0, 0.0);
    Field<Real> block_diagonal(block_extent, 4, 0, 0.0);
    Field<Real> block_coupling(block_extent, 4, 0, 0.0);
    block_rhs(0, 0, 0, 0) = 1.0;
    block_rhs(0, 0, 0, 1) = 2.0;
    block_diagonal(0, 0, 0, 0) = 4.0;
    block_diagonal(0, 0, 0, 1) = 1.0;
    block_diagonal(0, 0, 0, 2) = 2.0;
    block_diagonal(0, 0, 0, 3) = 3.0;
    const auto block_increment
        = solve_block_lu_sgs(block_rhs, block_diagonal, block_coupling, 2);
    WCNS_REQUIRE_NEAR(block_increment(0, 0, 0, 0), 0.1, 2.0e-11);
    WCNS_REQUIRE_NEAR(block_increment(0, 0, 0, 1), 0.6, 2.0e-11);
    block_diagonal(0, 0, 0, 3) = 0.5;
    WCNS_REQUIRE_THROWS(std::runtime_error,
                        solve_block_lu_sgs(
                            block_rhs, block_diagonal, block_coupling, 2));

    WCNS_REQUIRE_NEAR(bdf_time_diagonal(BdfOrder::First, 0.25), 4.0, 0.0);
    WCNS_REQUIRE_NEAR(bdf_time_diagonal(BdfOrder::Second, 0.25), 6.0, 0.0);
    WCNS_REQUIRE_NEAR(
        bdf_time_residual(BdfOrder::First, 2.0, 1.0, 0.5, 0.25), 4.0, 0.0);
    WCNS_REQUIRE_NEAR(
        bdf_time_residual(BdfOrder::Second, 2.0, 1.0, 0.5, 0.25), 5.0, 0.0);
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        bdf_time_diagonal(BdfOrder::First, 0.0));
    const Real coarse_bdf2_error = bdf2_decay_error(0.1);
    const Real medium_bdf2_error = bdf2_decay_error(0.05);
    const Real fine_bdf2_error = bdf2_decay_error(0.025);
    WCNS_REQUIRE(std::log(coarse_bdf2_error / medium_bdf2_error) / std::log(2.0)
                 > 1.9);
    WCNS_REQUIRE(std::log(medium_bdf2_error / fine_bdf2_error) / std::log(2.0)
                 > 1.9);

    const auto gas = GasModel::from_input(GasModelInput {1.4, {}, 287.0});
    const auto reference = ReferenceScales::derive(
        ReferenceInput {34.0, 1.2, 300.0, 1.0, 1.8e-5, {}, {}}, gas);
    const NumericalFloors floors;
    WeissSmithParameters parameters;
    parameters.mach_cutoff = 1.0e-3;
    parameters.viscous_cutoff = 1.0;
    const TemperaturePrimitiveState low_mach {{1.0, 0.1, 0.0, 0.0, 1.0}};
    const auto preconditioned = weiss_smith_state(
        low_mach, {1.0, 0.0, 0.0}, gas, reference, floors, parameters, 0.0, 2);
    WCNS_REQUIRE(preconditioned.beta < 1.0e-2);
    WCNS_REQUIRE(std::abs(preconditioned.eigenvalues[4])
                 < thermodynamic_sound_speed(low_mach, gas, reference, floors, 2));

    const auto identity = multiply(preconditioned.gamma, preconditioned.inverse);
    for (int row = 0; row < 5; ++row) {
        for (int column = 0; column < 5; ++column) {
            WCNS_REQUIRE_NEAR(identity[static_cast<std::size_t>(row)]
                                     [static_cast<std::size_t>(column)],
                              row == column ? 1.0 : 0.0,
                              5.0e-12);
        }
    }

    const auto pressure_low_mach
        = pressure_primitive(low_mach, gas, reference, floors, 2);
    const auto pressure_preconditioned = weiss_smith_pressure_state(
        pressure_low_mach, {1.0, 0.0, 0.0}, gas, floors, parameters, 0.0, 2);
    const auto pressure_identity
        = multiply(pressure_preconditioned.gamma, pressure_preconditioned.inverse);
    for (int row = 0; row < 5; ++row) {
        for (int column = 0; column < 5; ++column) {
            WCNS_REQUIRE_NEAR(pressure_identity[static_cast<std::size_t>(row)]
                                               [static_cast<std::size_t>(column)],
                              row == column ? 1.0 : 0.0,
                              5.0e-12);
        }
    }
    const ConservativeState residual {{0.2, -0.4, 0.1, 0.0, 0.7}};
    const auto transformed = weiss_smith_precondition_residual(
        pressure_low_mach, residual, gas, floors, parameters, 0.0, 2);
    WCNS_REQUIRE(std::abs(transformed[0] - residual[0]) > 1.0e-4);
    for (const Real value : transformed) WCNS_REQUIRE(std::isfinite(value));

    const TemperaturePrimitiveState sonic {{1.0, 20.0, 0.0, 0.0, 1.0}};
    const auto unpreconditioned = weiss_smith_state(
        sonic, {1.0, 0.0, 0.0}, gas, reference, floors, parameters, 0.0, 2);
    WCNS_REQUIRE_NEAR(unpreconditioned.beta, 1.0, 1.0e-15);
    const auto pressure_sonic = pressure_primitive(sonic, gas, reference, floors, 2);
    const auto identity_residual = weiss_smith_precondition_residual(
        pressure_sonic, residual, gas, floors, parameters, 0.0, 2);
    for (int component = 0; component < euler_components; ++component) {
        WCNS_REQUIRE_NEAR(identity_residual[static_cast<std::size_t>(component)],
                          residual[static_cast<std::size_t>(component)],
                          2.0e-12);
    }
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        weiss_smith_state(low_mach,
                                          {2.0, 0.0, 0.0},
                                          gas,
                                          reference,
                                          floors,
                                          parameters,
                                          0.0,
                                          2));
}
