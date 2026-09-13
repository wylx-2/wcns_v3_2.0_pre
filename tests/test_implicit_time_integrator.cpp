#include "test_support.hpp"

#include <wcns/solver/implicit_time_integrator.hpp>
#include <wcns/solver/low_mach_preconditioner.hpp>

#include <cmath>

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

    WCNS_REQUIRE_NEAR(bdf_time_diagonal(BdfOrder::First, 0.25), 4.0, 0.0);
    WCNS_REQUIRE_NEAR(bdf_time_diagonal(BdfOrder::Second, 0.25), 6.0, 0.0);
    WCNS_REQUIRE_NEAR(
        bdf_time_residual(BdfOrder::First, 2.0, 1.0, 0.5, 0.25), 4.0, 0.0);
    WCNS_REQUIRE_NEAR(
        bdf_time_residual(BdfOrder::Second, 2.0, 1.0, 0.5, 0.25), 5.0, 0.0);
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        bdf_time_diagonal(BdfOrder::First, 0.0));

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
                              1.0e-10);
        }
    }

    const TemperaturePrimitiveState sonic {{1.0, 20.0, 0.0, 0.0, 1.0}};
    const auto unpreconditioned = weiss_smith_state(
        sonic, {1.0, 0.0, 0.0}, gas, reference, floors, parameters, 0.0, 2);
    WCNS_REQUIRE_NEAR(unpreconditioned.beta, 1.0, 1.0e-15);
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
