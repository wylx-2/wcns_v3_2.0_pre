#pragma once

#include <wcns/runtime/quantity_registry.hpp>
#include <wcns/runtime/simulation_driver.hpp>
#include <wcns/runtime/weighted_statistics.hpp>

#include <fstream>
#include <functional>
#include <string>
#include <vector>

namespace wcns {

enum class OutputCategory {
    Field,
    Statistics,
    Boundary,
    Checkpoint,
};

class OutputSchedule {
public:
    explicit OutputSchedule(OutputScheduleConfig config = {});

    [[nodiscard]] Real next_time(Real current_time) const;
    [[nodiscard]] bool consume(const SimulationState& state, bool initial, bool final);

private:
    [[nodiscard]] bool scheduled(const SimulationState& state) const;
    [[nodiscard]] bool already_emitted(const SimulationState& state) const;

    OutputScheduleConfig config_;
    bool has_emitted_ = false;
    std::size_t last_step_ = 0;
    Real last_time_ = 0.0;
};

class RuntimeOutputManager final : public ISimulationObserver {
public:
    using EventWriter = std::function<std::vector<std::string>(
        OutputCategory, const SimulationState&, bool, bool)>;

    RuntimeOutputManager(const MpiRuntime& mpi,
                         const CaseConfig& config,
                         const StructuredPartitionPlan& partition,
                         std::string mesh_signature = {},
                         const StatisticContext* statistic_context = nullptr,
                         EventWriter event_writer = {},
                         StatisticRegistry statistic_registry
                         = StatisticRegistry::create_builtin(),
                         AcceptedTimeStatistics* time_statistics = nullptr);
    ~RuntimeOutputManager() override;

    [[nodiscard]] Real next_time_event(const SimulationState& state) const override;
    void on_initial(const SimulationState& state) override;
    void on_step(const SimulationState& state, bool residual_checked) override;
    void on_final(const SimulationState& state) override;

    void record_file(std::string path);
    [[nodiscard]] const std::vector<std::string>& files() const noexcept { return files_; }

private:
    void prepare_directory();
    void dispatch(OutputCategory category,
                  OutputSchedule& schedule,
                  const SimulationState& state,
                  bool initial,
                  bool final,
                  bool enabled);
    void write_history(const SimulationState& state, bool residual_checked);
    void finish_history();
    void write_statistics(const SimulationState& state);
    void sample_time_statistics(const SimulationState& state);
    void finish_statistics();
    void finish_time_statistics();
    void write_manifest(const SimulationState& state);

    const MpiRuntime& mpi_;
    const CaseConfig& config_;
    const StructuredPartitionPlan& partition_;
    std::string mesh_signature_;
    const StatisticContext* statistic_context_ = nullptr;
    StatisticRegistry statistic_registry_;
    AcceptedTimeStatistics* time_statistics_ = nullptr;
    EventWriter event_writer_;
    OutputSchedule field_schedule_;
    OutputSchedule history_schedule_;
    OutputSchedule statistics_schedule_;
    OutputSchedule boundary_schedule_;
    OutputSchedule checkpoint_schedule_;
    bool prepared_ = false;
    bool finalized_ = false;
    std::string output_directory_;
    std::string history_temporary_path_;
    std::string history_final_path_;
    std::ofstream history_stream_;
    std::string statistics_temporary_path_;
    std::string statistics_final_path_;
    std::ofstream statistics_stream_;
    std::string time_statistics_temporary_path_;
    std::string time_statistics_final_path_;
    std::vector<std::string> files_;
    std::vector<std::size_t> turbulence_floor_repairs_total_;
    std::size_t last_diagnostic_step_ = 0;
};

} // namespace wcns
