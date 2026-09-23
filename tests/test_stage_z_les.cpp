#include "test_support.hpp"

#include <wcns/runtime/weighted_statistics.hpp>
#include <wcns/solver/les_closure.hpp>
#include <wcns/solver/les_filter.hpp>
#include <wcns/solver/les_models.hpp>

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace {

void test_filter()
{
    using namespace wcns;
    Field<Real> constant({4, 3, 2}, 1, 1, 7.25);
    Field<Real> filtered({4, 3, 2}, 1, 0, 0.0);
    les_box3_tensor_filter(constant, filtered);
    for (std::size_t entry = 0; entry < filtered.size(); ++entry) {
        WCNS_REQUIRE_NEAR(filtered.data()[entry], 7.25, 1.0e-15);
    }

    Field<Real> linear({5, 4, 3}, 1, 1, std::numeric_limits<Real>::quiet_NaN());
    for (int k = -1; k <= 3; ++k) {
        for (int j = -1; j <= 4; ++j) {
            for (int i = -1; i <= 5; ++i) {
                linear(i, j, k, 0) = 2.0 + 3.0 * i - 0.5 * j + 0.25 * k;
            }
        }
    }
    WCNS_REQUIRE_NEAR(les_box3_tensor_sample(linear, {2, 2, 1}, 0),
                      linear(2, 2, 1, 0),
                      1.0e-14);

    Field<Real> boundary({2, 2, 2}, 1, 1, std::numeric_limits<Real>::quiet_NaN());
    for (int k = 0; k < 2; ++k) {
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) boundary(i, j, k, 0) = 1.0 + i;
        }
    }
    // nearest-constant physical extension: 3/4 left value + 1/4 right value.
    WCNS_REQUIRE_NEAR(les_box3_tensor_sample(boundary, {0, 0, 0}, 0), 1.25, 1.0e-15);

    const Real pi = std::acos(-1.0);
    WCNS_REQUIRE_NEAR(les_box3_fourier_response({pi / 2.0, 0.0, 0.0}), 0.5, 1.0e-15);
    WCNS_REQUIRE_NEAR(les_box3_fourier_response({pi, 0.0, 0.0}), 0.0, 1.0e-15);

    Field<Real> primitive({3, 3, 3}, 4, 1, 0.0);
    for (int k = -1; k <= 3; ++k) {
        for (int j = -1; j <= 3; ++j) {
            for (int i = -1; i <= 3; ++i) {
                primitive(i, j, k, 0) = 2.0;
                primitive(i, j, k, 1) = static_cast<Real>(i);
                primitive(i, j, k, 2) = 1.0;
                primitive(i, j, k, 3) = -2.0;
            }
        }
    }
    const auto favre = favre_box3_tensor_sample(primitive, {1, 1, 1}, 0, {{1, 2, 3}});
    WCNS_REQUIRE_NEAR(favre.density, 2.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(favre.velocity[0], 1.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(favre.leonard_stress[0][0], 1.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(favre.leonard_stress[1][1], 0.0, 1.0e-15);
}

void test_models()
{
    using namespace wcns;
    const VelocityGradient strain {{
        {{1.0, 0.0, 0.0}},
        {{0.0, -0.5, 0.0}},
        {{0.0, 0.0, -0.5}},
    }};
    const auto deviatoric = les_deviatoric_strain(strain);
    WCNS_REQUIRE_NEAR(deviatoric.xx + deviatoric.yy + deviatoric.zz, 0.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(les_strain_magnitude(deviatoric), std::sqrt(3.0), 1.0e-15);

    const auto smag = evaluate_smagorinsky(2.0, 0.25, 0.2, strain);
    const Real expected_mu = 2.0 * std::pow(0.2 * 0.25, 2) * std::sqrt(3.0);
    WCNS_REQUIRE_NEAR(smag.eddy_viscosity, expected_mu, 1.0e-15);
    WCNS_REQUIRE_NEAR(smag.physical_sgs_stress.xx, -2.0 * expected_mu, 1.0e-15);
    WCNS_REQUIRE(smag.energy_transfer > 0.0);
    const auto zero = evaluate_smagorinsky(2.0, 0.25, 0.0, strain);
    WCNS_REQUIRE(zero.eddy_viscosity == 0.0);
    WCNS_REQUIRE(zero.physical_sgs_stress.xx == 0.0);

    SymmetricStress leonard {3.0, 0.0, 0.0, 0.5, -0.25, 0.75};
    const auto similarity = evaluate_scale_similarity(strain, leonard, 1.0);
    WCNS_REQUIRE_NEAR(similarity.physical_sgs_stress.xx, 2.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(similarity.physical_sgs_stress.yy, -1.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(similarity.physical_sgs_stress.zz, -1.0, 1.0e-15);
    const auto mixed
        = evaluate_mixed_smagorinsky_similarity(2.0, 0.25, 0.2, strain, leonard, 1.0);
    WCNS_REQUIRE_NEAR(mixed.physical_sgs_stress.xy, 0.5, 1.0e-15);
    WCNS_REQUIRE_NEAR(mixed.eddy_viscosity, expected_mu, 1.0e-15);

    WCNS_REQUIRE_NEAR(evaluate_dynamic_smagorinsky_coefficient(
                          2.0, 4.0, 1.0e-20, -0.05, 0.09),
                      0.09,
                      1.0e-15);
    WCNS_REQUIRE_NEAR(evaluate_dynamic_smagorinsky_coefficient(
                          -1.0, 10.0, 1.0e-20, -0.05, 0.09),
                      -0.05,
                      1.0e-15);
    const auto dynamic = evaluate_dynamic_smagorinsky(1.0, 0.2, -0.01, strain);
    WCNS_REQUIRE(dynamic.eddy_viscosity < 0.0);
    WCNS_REQUIRE(dynamic.energy_transfer < 0.0);

    const VelocityGradient simple_shear {{
        {{0.0, 1.0, 0.0}},
        {{0.0, 0.0, 0.0}},
        {{0.0, 0.0, 0.0}},
    }};
    const auto wale_shear = evaluate_wale(1.0, 0.5, 0.325, simple_shear);
    WCNS_REQUIRE_NEAR(wale_shear.eddy_viscosity, 0.0, 1.0e-15);
    const VelocityGradient general {{
        {{0.2, 0.7, 0.1}},
        {{-0.3, -0.1, 0.4}},
        {{0.2, -0.5, -0.1}},
    }};
    const auto wale_general = evaluate_wale(1.2, 0.3, 0.325, general);
    WCNS_REQUIRE(wale_general.eddy_viscosity > 0.0);
    WCNS_REQUIRE(wale_general.energy_transfer >= 0.0);
    WCNS_REQUIRE_NEAR(les_van_driest_factor(0.0), 0.0, 1.0e-15);
    WCNS_REQUIRE(les_van_driest_factor(1000.0) > 0.999999999999999);

    TurbulenceModelConfig config;
    config.kind = TurbulenceModelKind::Smagorinsky;
    const auto model = TurbulenceModelRegistry::create_builtin().create(config);
    WCNS_REQUIRE(model->family() == TurbulenceModelFamily::LesAlgebraic);
    WCNS_REQUIRE(model->fields().empty());
    TurbulenceCellContext context;
    context.mean_state = {{2.0, 0.0, 0.0, 0.0, 1.0}};
    context.primitive_gradients[0] = strain[0];
    context.primitive_gradients[1] = strain[1];
    context.primitive_gradients[2] = strain[2];
    context.filter_width = 0.25;
    context.reference_mach = 0.2;
    context.heat_capacity_ratio = 1.4;
    context.dimension = 3;
    const auto contribution = model->viscous_contribution(context);
    const auto reference = evaluate_smagorinsky(2.0, 0.25, 0.17, strain);
    WCNS_REQUIRE_NEAR(contribution.stress.xx,
                      -reference.physical_sgs_stress.xx,
                      1.0e-15);
    WCNS_REQUIRE_NEAR(contribution.sgs_energy_transfer,
                      reference.energy_transfer,
                      1.0e-15);
}

void test_closure()
{
    using namespace wcns;
    const Extent3 extent {4, 4, 4};
    Field<Real> resolved(extent, les_resolved_components, 1, 0.0);
    for (int k = -1; k <= 4; ++k) {
        for (int j = -1; j <= 4; ++j) {
            for (int i = -1; i <= 4; ++i) {
                resolved(i, j, k, 0) = 1.0;
                resolved(i, j, k, 1) = 2.0;
                resolved(i, j, k, 2) = -1.0;
                resolved(i, j, k, 3) = 0.5;
                for (int component = 4; component < 13; ++component) {
                    resolved(i, j, k, component) = 0.0;
                }
                resolved(i, j, k, 13) = 0.25;
            }
        }
    }
    Field<Real> moments(extent, les_dynamic_moment_components, 1, 0.0);
    Field<Real> closure(extent, les_closure_components, 1, 0.0);
    TurbulenceModelConfig config;
    config.kind = TurbulenceModelKind::DynamicSmagorinsky;
    compute_les_dynamic_moments(moments, resolved, config);
    complete_les_fixed_stencil_ghosts(moments);
    compute_les_closure_field(closure, resolved, moments, config);
    complete_les_fixed_stencil_ghosts(closure);
    for (int component = les_leonard_xx; component <= les_leonard_yz; ++component) {
        WCNS_REQUIRE_NEAR(closure(1, 1, 1, component), 0.0, 1.0e-15);
    }
    WCNS_REQUIRE_NEAR(closure(1, 1, 1, les_dynamic_coefficient), 0.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(closure(1, 1, 1, les_filter_width), 0.25, 1.0e-15);

    for (int k = -1; k <= 4; ++k) {
        for (int j = -1; j <= 4; ++j) {
            for (int i = -1; i <= 4; ++i) {
                resolved(i, j, k, 1) = static_cast<Real>(i);
                resolved(i, j, k, 4) = 1.0;
            }
        }
    }
    compute_les_dynamic_moments(moments, resolved, config);
    complete_les_fixed_stencil_ghosts(moments);
    compute_les_closure_field(closure, resolved, moments, config);
    complete_les_fixed_stencil_ghosts(closure);
    WCNS_REQUIRE_NEAR(closure(2, 2, 2, les_leonard_xx), 1.0 / 3.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(closure(2, 2, 2, les_leonard_yy), -1.0 / 6.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(closure(2, 2, 2, les_leonard_zz), -1.0 / 6.0, 1.0e-14);
    WCNS_REQUIRE(std::isfinite(closure(2, 2, 2, les_dynamic_coefficient)));
    WCNS_REQUIRE(closure(2, 2, 2, les_dynamic_coefficient)
                 >= config.les_dynamic_coefficient_minimum);
    WCNS_REQUIRE(closure(2, 2, 2, les_dynamic_coefficient)
                 <= config.les_dynamic_coefficient_maximum);
}

void test_statistics()
{
    using namespace wcns;
    WeightedMomentState moments;
    moments.add(1.0, 1.0);
    moments.add(3.0, 3.0);
    WCNS_REQUIRE_NEAR(moments.mean, 2.5, 1.0e-15);
    WCNS_REQUIRE_NEAR(moments.variance(), 0.75, 1.0e-15);

    WeightedMomentState left;
    WeightedMomentState right;
    left.add(1.0, 1.0);
    right.add(3.0, 3.0);
    left.merge(right);
    WCNS_REQUIRE_NEAR(left.mean, moments.mean, 1.0e-15);
    WCNS_REQUIRE_NEAR(left.second_central, moments.second_central, 1.0e-15);

    WeightedCovarianceState covariance;
    covariance.add(1.0, 2.0, 1.0);
    covariance.add(3.0, 6.0, 1.0);
    WCNS_REQUIRE_NEAR(covariance.covariance(), 2.0, 1.0e-15);

    AcceptedTimeStatistics statistics("stage-z-test", 2);
    WCNS_REQUIRE(!statistics.sample(0, 0.0, 0.1, {100.0, 100.0}, 1.0, false));
    WCNS_REQUIRE(statistics.sample(1, 0.1, 0.1, {1.0, 4.0}, 1.0, true));
    WCNS_REQUIRE(!statistics.sample(1, 0.1, 0.1, {9.0, 9.0}, 9.0, true));
    WCNS_REQUIRE(statistics.sample(2, 0.3, 0.2, {3.0, 8.0}, 2.0, true));
    WCNS_REQUIRE(statistics.state().accepted_events == 2);
    WCNS_REQUIRE_NEAR(statistics.state().reynolds[0].mean, 7.0 / 3.0, 1.0e-15);
    WCNS_REQUIRE_NEAR(statistics.state().favre[0].mean, 13.0 / 5.0, 1.0e-15);

    auto restored = AcceptedTimeStatistics::deserialize(statistics.serialize());
    WCNS_REQUIRE(restored.state().identity == "stage-z-test");
    WCNS_REQUIRE(restored.state().accepted_events == 2);
    WCNS_REQUIRE_NEAR(restored.state().favre[1].mean,
                      statistics.state().favre[1].mean,
                      1.0e-15);
    WCNS_REQUIRE(restored.sample(3, 0.4, 0.1, {5.0, 12.0}, 1.0, true));
    WCNS_REQUIRE_THROWS(std::invalid_argument,
                        restored.sample(2, 0.2, 0.1, {0.0, 0.0}, 1.0, true));
}

} // namespace

void test_stage_z_les()
{
    test_filter();
    test_models();
    test_closure();
    test_statistics();
}
