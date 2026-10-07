#pragma once

#include <wcns/solver/inviscid_wcns_solver.hpp>
#include <wcns/solver/les_closure.hpp>
#include <wcns/solver/rans_two_equation_transport.hpp>
#include <wcns/solver/sa_negative_transport.hpp>
#include <wcns/solver/viscous_operator.hpp>
#include <functional>

namespace wcns {

struct ViscousStabilityCoefficients {
    Real phenglei_2d_ssprk3 = 4.0;
    Real phenglei_3d_ssprk3 = 4.0;
    Real scmm6_2d_ssprk3 = 4.0;
    Real scmm6_3d_ssprk3 = 4.0;

    void validate() const;
    [[nodiscard]] Real for_ssprk3(AlgorithmProfileKind profile, int dimension) const;
    [[nodiscard]] std::string summary() const;
};

struct ViscousWcnsConfig {
    InviscidWcnsConfig inviscid {};
    TransportConfig transport {};
    TurbulenceModelConfig turbulence {};
    ViscousStabilityCoefficients stability {};

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

class ViscousWcnsSolver {
public:
    ViscousWcnsSolver(const MpiRuntime& mpi,
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
                      ViscousWcnsConfig config = {});

    void set_pressure_gradient_x(Real force) { source_registry_.set_pressure_gradient_x(force); }
    void set_accepted_step_transform(std::function<void(Real)> transform) { accepted_step_transform_ = std::move(transform); }
    // Collective stage callback, recomputed on retries. It may only add to the RHS.
    void set_collective_source(std::function<void(Real)> source) { collective_source_ = std::move(source); }
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
    [[nodiscard]] WallFunctionDiagnostics global_wall_function_diagnostics() const;
    [[nodiscard]] std::vector<std::size_t> global_turbulence_floor_repairs() const;
    [[nodiscard]] const std::unordered_map<BlockId, Field<Real>>&
    turbulence_residuals() const noexcept
    {
        return turbulence_residual_workspace_;
    }

private:
    void compute_residuals_impl(Real stage_time,
                                int rk_stage,
                                const BlockFaceRobustnessMap* robustness_levels);
    [[nodiscard]] bool turbulence_active() const noexcept;
    [[nodiscard]] bool les_active() const noexcept;
    void update_les_closure_fields();
    void capture_turbulence_stage_state();
    void update_turbulence_stage(Real initial_weight,
                                 Real stage_weight,
                                 Real residual_weight);

    const MpiRuntime& mpi_;
    LocalBlockSet& local_blocks_;
    const StructuredMesh& global_mesh_;
    const DistributedTopology& topology_;
    HaloExchanger state_exchanger_;
    HaloExchanger turbulence_exchanger_;
    HaloExchanger les_resolved_exchanger_;
    HaloExchanger les_dynamic_moment_exchanger_;
    HaloExchanger les_closure_exchanger_;
    BlockMetricMap& metrics_;
    const BlockBoundaryDataMap& boundary_data_;
    AlgorithmProfile profile_;
    GasModel gas_;
    ReferenceScales reference_;
    NumericalFloors floors_;
    ViscousWcnsConfig config_;
    SourceTermRegistry source_registry_;
    std::function<void(Real)> collective_source_;
    std::function<void(Real)> accepted_step_transform_;
    TransportModel transport_;
    std::unique_ptr<ITurbulenceModel> turbulence_model_;
    RiemannSolver riemann_;
    RiemannSolver robust_riemann_;
    RobustnessLadder robustness_ladder_ {};
    std::vector<StructuredBlock*> block_workspace_;
    std::unordered_map<BlockId, InviscidFaceFluxField> inviscid_flux_workspace_;
    FaceFluxFieldRegistry inviscid_flux_registry_;
    FaceFluxHaloPlan inviscid_flux_plan_;
    FaceFluxHaloExchanger inviscid_flux_exchanger_;
    std::unordered_map<BlockId, GradientOperandFaceField> operand_workspace_;
    GradientOperandFieldRegistry operand_registry_;
    GradientOperandFaceHaloPlan operand_plan_;
    GradientOperandFaceHaloExchanger operand_exchanger_;
    std::unordered_map<BlockId, PrimitiveGradientField> gradient_workspace_;
    GradientFieldRegistry gradient_registry_;
    GradientHaloPlan gradient_plan_;
    GradientHaloExchanger gradient_exchanger_;
    std::unordered_map<BlockId, PrimitiveGradientField> turbulence_gradient_workspace_;
    GradientFieldRegistry turbulence_gradient_registry_;
    std::unordered_map<BlockId, Field<Real>> two_equation_face_workspace_;
    BlockFieldRegistry two_equation_face_registry_ {two_equation_face_workspace_components};
    std::unordered_map<BlockId, ViscousFaceFluxField> viscous_flux_workspace_;
    ViscousFaceFluxFieldRegistry viscous_flux_registry_;
    ViscousFaceFluxHaloPlan viscous_flux_plan_;
    ViscousFaceFluxHaloExchanger viscous_flux_exchanger_;
    std::unordered_map<BlockId, ViscousFaceFluxField> turbulence_flux_workspace_;
    ViscousFaceFluxFieldRegistry turbulence_flux_registry_;
    ViscousFaceFluxHaloPlan turbulence_flux_plan_;
    ViscousFaceFluxHaloExchanger turbulence_flux_exchanger_;
    std::unordered_map<BlockId, Field<Real>> turbulence_residual_workspace_;
    std::unordered_map<BlockId, Field<Real>> turbulence_source_jacobian_workspace_;
    std::unordered_map<BlockId, Field<Real>> les_resolved_workspace_;
    std::unordered_map<BlockId, Field<Real>> les_dynamic_moment_workspace_;
    std::unordered_map<BlockId, Field<Real>> les_closure_workspace_;
    BlockFieldRegistry les_resolved_registry_ {les_resolved_components};
    BlockFieldRegistry les_dynamic_moment_registry_ {les_dynamic_moment_components};
    BlockFieldRegistry les_closure_registry_ {les_closure_components};
    std::unordered_map<BlockId, std::vector<Real>> turbulence_initial_state_;
    std::unordered_map<BlockId, std::vector<Real>> turbulence_stage_state_;
    SsprkWorkspace time_workspace_;
    StateSnapshot previous_physical_state_;
    std::unordered_map<BlockId, std::vector<Real>> previous_turbulence_physical_state_;
    bool has_previous_physical_state_ = false;
    Real previous_physical_time_step_ = 0.0;
    std::uint64_t version_ = 0;
    ReconstructionDiagnostics reconstruction_diagnostics_ {};
    RiemannDiagnostics riemann_diagnostics_ {};
    RobustnessDiagnostics robustness_diagnostics_ {};
    std::vector<std::size_t> turbulence_floor_repairs_;
};

} // namespace wcns
