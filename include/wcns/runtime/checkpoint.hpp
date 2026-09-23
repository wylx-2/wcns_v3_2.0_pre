#pragma once

#include <wcns/runtime/field_output.hpp>
#include <wcns/runtime/weighted_statistics.hpp>

#include <string>
#include <vector>

namespace wcns {

struct CheckpointRestoreResult {
    SimulationInitialState initial;
    Real previous_time_step = 0.0;
};

class CheckpointService {
public:
    CheckpointService(const MpiRuntime& mpi,
                      const CaseConfig& config,
                      const StructuredPartitionPlan& partition,
                      LocalBlockSet& local_blocks,
                      const BlockMetricMap& metrics,
                      QuantityContext quantity_context,
                      std::string mesh_path,
                      AcceptedTimeStatistics* time_statistics = nullptr);

    [[nodiscard]] std::vector<std::string> write(const SimulationState& state) const;
    [[nodiscard]] CheckpointRestoreResult restore(const std::string& path) const;
    [[nodiscard]] const std::string& mesh_signature() const noexcept { return mesh_signature_; }

private:
    const MpiRuntime& mpi_;
    const CaseConfig& config_;
    const StructuredPartitionPlan& partition_;
    LocalBlockSet& local_blocks_;
    const BlockMetricMap& metrics_;
    QuantityContext quantity_context_;
    std::string mesh_path_;
    std::string mesh_signature_;
    FieldQuantityRegistry registry_;
    AcceptedTimeStatistics* time_statistics_ = nullptr;
};

} // namespace wcns
