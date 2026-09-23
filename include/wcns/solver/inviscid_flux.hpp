#pragma once

#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/mesh/linear_operators.hpp>
#include <wcns/mesh/structured_mesh.hpp>
#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/solver/physical_boundary.hpp>
#include <wcns/solver/riemann_solver.hpp>
#include <wcns/solver/wcns_reconstruction.hpp>

#include <cstdint>
#include <unordered_map>
#include <vector>

namespace wcns {

class FaceRobustnessField;
class RobustnessLadder;

enum class FluxDifferenceMode {
    Profile,
    ConservativeTwoPoint,
};

[[nodiscard]] const char* flux_difference_mode_name(FluxDifferenceMode mode);

class InviscidFaceFluxField {
public:
    InviscidFaceFluxField(Extent3 cell_extent,
                          int dimension,
                          AlgorithmProfileKind profile,
                          std::uint64_t version);

    [[nodiscard]] AlgorithmProfileKind profile() const noexcept { return profile_; }
    [[nodiscard]] std::uint64_t version() const noexcept { return version_; }
    [[nodiscard]] int halo_layers() const noexcept { return halo_layers_; }
    [[nodiscard]] int dimension() const noexcept { return dimension_; }

    [[nodiscard]] Field<Real>& field(Axis axis);
    [[nodiscard]] const Field<Real>& field(Axis axis) const;
    void reset(std::uint64_t version);

private:
    AlgorithmProfileKind profile_;
    std::uint64_t version_;
    int halo_layers_;
    int dimension_;
    Field<Real> i_;
    Field<Real> j_;
    Field<Real> k_;
};

struct FaceFluxPair {
    Index3 receiver;
    Index3 donor;
    int layer = 0;
};

struct FaceFluxExchangeDescriptor {
    ConnectionId connection = invalid_connection_id;
    BlockId receiver_block = invalid_block_id;
    BlockId donor_block = invalid_block_id;
    RankId receiver_rank = invalid_rank_id;
    RankId donor_rank = invalid_rank_id;
    BlockId shared_face_owner = invalid_block_id;
    int direction = 0;
    Axis receiver_axis = Axis::I;
    Axis donor_axis = Axis::I;
    Real orientation = 1.0;
    PeriodicTransform periodic {};
    AlgorithmProfileKind profile = AlgorithmProfileKind::PhengleiWcns;
    std::uint64_t version = 0;
    std::vector<FaceFluxPair> pairs;

    [[nodiscard]] int message_tag(int tag_base = 12288) const;
};

[[nodiscard]] ConservativeState
transform_inviscid_face_flux_for_receiver(const ConservativeState& donor,
                                          const FaceFluxExchangeDescriptor& descriptor);

[[nodiscard]] bool
is_non_owned_connection_face(const StructuredBlock& block, Axis axis, Index3 face);

[[nodiscard]] StencilRow inviscid_residual_stencil(const StructuredBlock& block,
                                                   const AlgorithmProfile& profile,
                                                   FluxDifferenceMode mode,
                                                   Axis axis,
                                                   Index3 cell);

class FaceFluxHaloPlan {
public:
    [[nodiscard]] static FaceFluxHaloPlan
    build(const StructuredMesh& mesh, const AlgorithmProfile& profile, std::uint64_t version);

    [[nodiscard]] const std::vector<FaceFluxExchangeDescriptor>& exchanges() const noexcept
    {
        return exchanges_;
    }
    void set_version(std::uint64_t version);

private:
    std::vector<FaceFluxExchangeDescriptor> exchanges_;
};

class FaceFluxFieldRegistry {
public:
    void add(BlockId block, InviscidFaceFluxField& field);
    [[nodiscard]] bool contains(BlockId block) const noexcept;
    [[nodiscard]] InviscidFaceFluxField& field(BlockId block) const;

private:
    std::unordered_map<BlockId, InviscidFaceFluxField*> fields_;
};

class FaceFluxHaloExchanger {
public:
    FaceFluxHaloExchanger(const MpiRuntime& mpi, const FaceFluxHaloPlan& plan)
        : mpi_(mpi)
        , plan_(plan)
    {
        prepare();
    }

    void exchange(const FaceFluxFieldRegistry& fields) const;

private:
    struct Pending {
        const FaceFluxExchangeDescriptor* descriptor = nullptr;
        std::vector<Real> values;
    };

    void prepare();

    const MpiRuntime& mpi_;
    const FaceFluxHaloPlan& plan_;
    mutable std::vector<Pending> receives_;
    mutable std::vector<Pending> sends_;
#if WCNS_HAS_MPI
    mutable std::vector<MPI_Request> requests_;
#endif
};

[[nodiscard]] InviscidFaceFluxField
compute_inviscid_face_fluxes(const StructuredBlock& block,
                             const MetricField& metric,
                             const AlgorithmProfile& profile,
                             const ReconstructionConfig& reconstruction,
                             const RiemannSolver& riemann,
                             const GasModel& gas,
                             const ReferenceScales& reference,
                             const NumericalFloors& floors,
                             const BoundaryDataMap& boundary_data,
                             const InviscidBoundaryOptions& boundary_options,
                             std::uint64_t version,
                             ReconstructionDiagnostics& diagnostics,
                             RiemannDiagnostics* riemann_diagnostics = nullptr,
                             int rk_stage = 0,
                             Real stage_time = 0.0,
                             const FaceRobustnessField* robustness_levels = nullptr,
                             const RobustnessLadder* robustness_ladder = nullptr,
                             const RiemannSolver* robust_riemann = nullptr,
                             Real viscous_preconditioner_scale = 0.0);

void compute_inviscid_face_fluxes_into(InviscidFaceFluxField& result,
                                       const StructuredBlock& block,
                                       const MetricField& metric,
                                       const AlgorithmProfile& profile,
                                       const ReconstructionConfig& reconstruction,
                                       const RiemannSolver& riemann,
                                       const GasModel& gas,
                                       const ReferenceScales& reference,
                                       const NumericalFloors& floors,
                                       const BoundaryDataMap& boundary_data,
                                       const InviscidBoundaryOptions& boundary_options,
                                       std::uint64_t version,
                                       ReconstructionDiagnostics& diagnostics,
                                       RiemannDiagnostics* riemann_diagnostics = nullptr,
                                       int rk_stage = 0,
                                       Real stage_time = 0.0,
                                       const FaceRobustnessField* robustness_levels = nullptr,
                                       const RobustnessLadder* robustness_ladder = nullptr,
                                       const RiemannSolver* robust_riemann = nullptr,
                                       Real viscous_preconditioner_scale = 0.0);

void compute_wcns_inviscid_residual(StructuredBlock& block,
                                    const MetricField& metric,
                                    const InviscidFaceFluxField& flux,
                                    const AlgorithmProfile& profile,
                                    FluxDifferenceMode mode = FluxDifferenceMode::Profile);

} // namespace wcns
