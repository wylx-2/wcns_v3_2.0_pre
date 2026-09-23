#pragma once

#include <wcns/runtime/stop_controller.hpp>
#include <wcns/solver/viscous_wcns_solver.hpp>

#include <functional>
#include <array>
#include <limits>
#include <string>
#include <vector>

namespace wcns {

struct SolverDiagnostics {
    std::size_t reconstruction_fallbacks = 0;
    std::size_t riemann_faces = 0;
    std::size_t riemann_fallbacks = 0;
    std::array<std::size_t, 4> robustness_face_levels {};
    std::size_t robustness_troubled_cells = 0;
    std::size_t robustness_local_recomputations = 0;
    std::size_t robustness_step_retries = 0;
    Real robustness_proposed_time_step = 0.0;
    Real robustness_accepted_time_step = 0.0;
    Real minimum_density = std::numeric_limits<Real>::infinity();
    Real minimum_pressure = std::numeric_limits<Real>::infinity();
    Real minimum_temperature = std::numeric_limits<Real>::infinity();
    Real minimum_internal_energy = std::numeric_limits<Real>::infinity();
    std::size_t wall_function_faces = 0;
    std::size_t wall_function_out_of_range = 0;
    Real wall_function_minimum_y_plus = std::numeric_limits<Real>::infinity();
    Real wall_function_maximum_y_plus = 0.0;
    std::vector<std::size_t> turbulence_floor_repairs;
};

class ISimulationSolver {
public:
    virtual ~ISimulationSolver() = default;

    [[nodiscard]] virtual Real global_time_step(Real cfl) = 0;
    [[nodiscard]] virtual Real advance(Real time_step, Real initial_time) = 0;
    virtual void refresh_residuals(Real time) = 0;
    [[nodiscard]] virtual ResidualNorms residual_norms() const = 0;
    [[nodiscard]] virtual SolverDiagnostics diagnostics() const = 0;
};

class InviscidSimulationSolver final : public ISimulationSolver {
public:
    InviscidSimulationSolver(InviscidWcnsSolver& solver,
                             const MpiRuntime& mpi,
                             const LocalBlockSet& local_blocks,
                             const BlockMetricMap& metrics,
                             const StructuredPartitionPlan& partition,
                             AlgorithmProfile profile,
                             TimeAlgorithmConfig time_algorithm = {},
                             RunMode run_mode = RunMode::Steady);

    [[nodiscard]] Real global_time_step(Real cfl) override;
    [[nodiscard]] Real advance(Real time_step, Real initial_time) override;
    void refresh_residuals(Real time) override;
    [[nodiscard]] ResidualNorms residual_norms() const override;
    [[nodiscard]] SolverDiagnostics diagnostics() const override;

private:
    InviscidWcnsSolver& solver_;
    const MpiRuntime& mpi_;
    const LocalBlockSet& local_blocks_;
    const BlockMetricMap& metrics_;
    const StructuredPartitionPlan& partition_;
    AlgorithmProfile profile_;
    TimeAlgorithmConfig time_algorithm_;
    RunMode run_mode_ = RunMode::Steady;
    Real pending_pseudo_cfl_ = 0.0;
};

class ViscousSimulationSolver final : public ISimulationSolver {
public:
    ViscousSimulationSolver(ViscousWcnsSolver& solver,
                            const MpiRuntime& mpi,
                            const LocalBlockSet& local_blocks,
                            const BlockMetricMap& metrics,
                            const StructuredPartitionPlan& partition,
                            AlgorithmProfile profile,
                            TimeAlgorithmConfig time_algorithm = {},
                            RunMode run_mode = RunMode::Steady);

    [[nodiscard]] Real global_time_step(Real cfl) override;
    [[nodiscard]] Real advance(Real time_step, Real initial_time) override;
    void refresh_residuals(Real time) override;
    [[nodiscard]] ResidualNorms residual_norms() const override;
    [[nodiscard]] SolverDiagnostics diagnostics() const override;

private:
    ViscousWcnsSolver& solver_;
    const MpiRuntime& mpi_;
    const LocalBlockSet& local_blocks_;
    const BlockMetricMap& metrics_;
    const StructuredPartitionPlan& partition_;
    AlgorithmProfile profile_;
    TimeAlgorithmConfig time_algorithm_;
    RunMode run_mode_ = RunMode::Steady;
    Real pending_pseudo_cfl_ = 0.0;
};

struct SimulationState {
    std::size_t step = 0;
    Real time = 0.0;
    Real time_step = 0.0;
    Real wall_time = 0.0;
    ResidualNorms residuals;
    SolverDiagnostics diagnostics;
    SteadyConvergenceState steady;
    StopReason stop_reason = StopReason::Running;
};

class ISimulationObserver {
public:
    virtual ~ISimulationObserver() = default;

    [[nodiscard]] virtual Real next_time_event(const SimulationState&) const
    {
        return std::numeric_limits<Real>::infinity();
    }
    virtual void on_initial(const SimulationState&) { }
    virtual void on_step(const SimulationState&, bool) { }
    virtual void on_final(const SimulationState&) { }
};

class NullSimulationObserver final : public ISimulationObserver { };

class CompositeSimulationObserver final : public ISimulationObserver {
public:
    void add(ISimulationObserver& observer);
    [[nodiscard]] Real next_time_event(const SimulationState& state) const override;
    void on_initial(const SimulationState& state) override;
    void on_step(const SimulationState& state, bool residual_checked) override;
    void on_final(const SimulationState& state) override;

private:
    std::vector<ISimulationObserver*> observers_;
};

struct SimulationInitialState {
    std::size_t step = 0;
    Real time = 0.0;
    SteadyConvergenceState steady;
};

class SimulationDriver {
public:
    using StopRequest = std::function<bool()>;
    using WallClock = std::function<Real()>;

    SimulationDriver(const MpiRuntime& mpi,
                     ISimulationSolver& solver,
                     CaseRunConfig config,
                     ISimulationObserver& observer,
                     StopRequest stop_requested = {},
                     WallClock wall_clock = {});

    [[nodiscard]] SimulationState run(SimulationInitialState initial = {});

private:
    [[nodiscard]] Real limited_time_step(Real proposed, const SimulationState& state) const;
    [[nodiscard]] bool all_ranks_succeeded(bool local_success) const;

    const MpiRuntime& mpi_;
    ISimulationSolver& solver_;
    CaseRunConfig config_;
    ISimulationObserver& observer_;
    StopRequest stop_requested_;
    WallClock wall_clock_;
};

} // namespace wcns
