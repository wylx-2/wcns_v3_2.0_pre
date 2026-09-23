#include "test_support.hpp"

#include <wcns/runtime/stop_controller.hpp>

namespace {

wcns::ResidualNorms uniform_norm(wcns::Real value)
{
    wcns::ResidualNorms result;
    result.l2.fill(value);
    result.linf.fill(2.0 * value);
    return result;
}

} // namespace

// 验收定常参考残差冻结、连续检查计数和收敛优先于最大步数。
void test_stop_controller()
{
    {
        wcns::CaseRunConfig config;
        config.mode = wcns::RunMode::Steady;
        config.max_steps = 3;
        config.steady.min_steps = 2;
        config.steady.consecutive_checks = 2;
        config.steady.l2_absolute = 1.0e-20;
        config.steady.l2_relative = 1.0e-6;
        config.steady.linf_absolute = 1.0e-20;
        config.steady.linf_relative = 1.0e-6;
        wcns::StopController controller(config);

        wcns::SimulationProgress progress;
        progress.step = 1;
        progress.time_step = 0.1;
        progress.time = 0.1;
        progress.residuals = uniform_norm(1.0);
        WCNS_REQUIRE(controller.evaluate(progress).reason == wcns::StopReason::Running);
        WCNS_REQUIRE(controller.steady_state().reference_initialized);
        WCNS_REQUIRE_NEAR(controller.steady_state().reference_l2[0], 1.0, 1.0e-15);

        progress.step = 2;
        progress.time = 0.2;
        progress.residuals = uniform_norm(1.0e-8);
        WCNS_REQUIRE(controller.evaluate(progress).reason == wcns::StopReason::Running);
        WCNS_REQUIRE(controller.steady_state().consecutive_passes == 1);

        progress.step = 3;
        progress.time = 0.3;
        WCNS_REQUIRE(controller.evaluate(progress).reason == wcns::StopReason::SteadyConverged);
        WCNS_REQUIRE(wcns::stop_reason_exit_code(wcns::StopReason::SteadyConverged) == 0);
    }
    {
        wcns::CaseRunConfig config;
        config.mode = wcns::RunMode::Steady;
        config.max_steps = 10;
        config.steady.min_steps = 1;
        config.steady.consecutive_checks = 1;
        config.steady.l2_relative = 1.0e-6;
        config.steady.linf_relative = 1.0e-6;
        wcns::StopController controller(config);

        wcns::SimulationProgress progress;
        progress.step = 1;
        progress.time = 0.1;
        progress.time_step = 0.1;
        progress.residuals = uniform_norm(1.0);
        progress.residuals.model_l2 = {1.0};
        progress.residuals.model_linf = {2.0};
        WCNS_REQUIRE(controller.evaluate(progress).reason == wcns::StopReason::Running);
        WCNS_REQUIRE(controller.steady_state().model_reference_l2.size() == 1);

        progress.step = 2;
        progress.time = 0.2;
        progress.residuals = uniform_norm(1.0e-8);
        progress.residuals.model_l2 = {1.0};
        progress.residuals.model_linf = {2.0};
        WCNS_REQUIRE(controller.evaluate(progress).reason == wcns::StopReason::Running);

        progress.step = 3;
        progress.time = 0.3;
        progress.residuals.model_l2 = {1.0e-8};
        progress.residuals.model_linf = {2.0e-8};
        WCNS_REQUIRE(controller.evaluate(progress).reason == wcns::StopReason::SteadyConverged);
    }
    {
        wcns::CaseRunConfig config;
        config.mode = wcns::RunMode::Unsteady;
        config.end_time = 0.5;
        config.max_steps = 2;
        wcns::StopController controller(config);
        wcns::SimulationProgress progress;
        progress.step = 2;
        progress.time = 0.5;
        progress.time_step = 0.25;
        progress.wall_time = 100.0;
        progress.user_signal = true;
        progress.residuals = uniform_norm(1.0);
        WCNS_REQUIRE(controller.evaluate(progress).reason == wcns::StopReason::PhysicalTimeReached);
    }
    {
        wcns::CaseRunConfig config;
        config.mode = wcns::RunMode::Unsteady;
        config.end_time = 1.0;
        config.max_steps = 10;
        wcns::StopController controller(config);
        wcns::SimulationProgress progress;
        progress.step = 1;
        progress.time = 0.1;
        progress.time_step = 0.1;
        progress.residuals = uniform_norm(1.0);
        progress.numerical_failure = true;
        WCNS_REQUIRE(controller.evaluate(progress).reason == wcns::StopReason::NumericalFailure);
        WCNS_REQUIRE(wcns::stop_reason_exit_code(wcns::StopReason::NumericalFailure) != 0);
    }
}
