#pragma once

#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/halo_exchanger.hpp>
#include <wcns/solver/inviscid_flux.hpp>
#include <wcns/solver/implicit_time_integrator.hpp>
#include <wcns/solver/robustness.hpp>
#include <wcns/solver/source_operator.hpp>
#include <wcns/solver/time_integrator.hpp>

#include <unordered_map>

namespace wcns {

using BlockMetricMap = std::unordered_map<BlockId, MetricField>;
using BlockBoundaryDataMap = std::unordered_map<BlockId, BoundaryDataMap>;

struct InviscidWcnsConfig {
    ReconstructionConfig reconstruction {};
    RiemannConfig riemann {};
    FluxDifferenceMode flux_difference = FluxDifferenceMode::Profile;
    InviscidBoundaryOptions boundary {};
    SourceTermConfig source_terms {};
    RobustnessConfig robustness {};

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

class InviscidWcnsSolver {
public:
    InviscidWcnsSolver(const MpiRuntime& mpi,
                       LocalBlockSet& local_blocks,
                       const StructuredMesh& global_mesh,
                       const DistributedTopology& topology,
                       int distribution_rank_count,
                       BlockMetricMap& metrics,
                       const BlockBoundaryDataMap& boundary_data,
                       AlgorithmProfile profile,
                       GasModel gas,
                       ReferenceScales reference,
                       NumericalFloors floors,
                       InviscidWcnsConfig config = {});

    void compute_residuals(Real stage_time, int rk_stage = 0);
    [[nodiscard]] Real advance(Real time_step, Real initial_time);
    [[nodiscard]] Real advance_lu_sgs(Real pseudo_cfl,
                                      Real time,
                                      const LuSgsIterationConfig& config = {});
    [[nodiscard]] Real advance_dual_time(Real physical_time_step,
                                         Real initial_time,
                                         const DualTimeIterationConfig& config);
    [[nodiscard]] Real global_time_step(Real cfl);

    [[nodiscard]] Real global_residual_l2() const;
    [[nodiscard]] const ReconstructionDiagnostics& reconstruction_diagnostics() const noexcept
    {
        return reconstruction_diagnostics_;
    }
    [[nodiscard]] const RiemannDiagnostics& riemann_diagnostics() const noexcept
    {
        return riemann_diagnostics_;
    }
    [[nodiscard]] std::size_t global_reconstruction_fallback_count() const;
    [[nodiscard]] std::size_t global_riemann_face_count() const;
    [[nodiscard]] std::size_t global_riemann_fallback_count() const;
    [[nodiscard]] RobustnessDiagnostics global_robustness_diagnostics() const;

private:
    void compute_residuals_impl(Real stage_time,
                                int rk_stage,
                                const BlockFaceRobustnessMap* robustness_levels);

    const MpiRuntime& mpi_;
    LocalBlockSet& local_blocks_;
    const StructuredMesh& global_mesh_;
    const DistributedTopology& topology_;
    HaloExchanger state_exchanger_;
    BlockMetricMap& metrics_;
    const BlockBoundaryDataMap& boundary_data_;
    AlgorithmProfile profile_;
    GasModel gas_;
    ReferenceScales reference_;
    NumericalFloors floors_;
    InviscidWcnsConfig config_;
    SourceTermRegistry source_registry_;
    RiemannSolver riemann_;
    RiemannSolver robust_riemann_;
    RobustnessLadder robustness_ladder_ {};
    std::vector<StructuredBlock*> block_workspace_;
    std::unordered_map<BlockId, InviscidFaceFluxField> face_flux_workspace_;
    FaceFluxFieldRegistry face_flux_registry_;
    FaceFluxHaloPlan face_flux_plan_;
    FaceFluxHaloExchanger face_flux_exchanger_;
    SsprkWorkspace time_workspace_;
    StateSnapshot previous_physical_state_;
    bool has_previous_physical_state_ = false;
    Real previous_physical_time_step_ = 0.0;
    std::uint64_t version_ = 0;
    ReconstructionDiagnostics reconstruction_diagnostics_ {};
    RiemannDiagnostics riemann_diagnostics_ {};
    RobustnessDiagnostics robustness_diagnostics_ {};
};

} // namespace wcns
