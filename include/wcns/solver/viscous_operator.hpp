#pragma once

#include <wcns/solver/les_closure.hpp>
#include <wcns/solver/turbulence_model.hpp>
#include <wcns/solver/viscous_boundary.hpp>
#include <wcns/solver/viscous_halo.hpp>

#include <stdexcept>
#include <unordered_map>

namespace wcns {

class ViscousFaceFluxField {
public:
    ViscousFaceFluxField(Extent3 cell_extent,
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

[[nodiscard]] ConservativeState
transform_viscous_face_flux_for_receiver(const ConservativeState& donor,
                                         const FaceFluxExchangeDescriptor& descriptor);

class ViscousFaceFluxHaloPlan {
public:
    [[nodiscard]] static ViscousFaceFluxHaloPlan
    build(const StructuredMesh& mesh, const AlgorithmProfile& profile, std::uint64_t version);

    [[nodiscard]] const std::vector<FaceFluxExchangeDescriptor>& exchanges() const noexcept
    {
        return exchanges_;
    }
    void set_version(std::uint64_t version);

private:
    std::vector<FaceFluxExchangeDescriptor> exchanges_;
};

class ViscousFaceFluxFieldRegistry {
public:
    void add(BlockId block, ViscousFaceFluxField& field);
    [[nodiscard]] bool contains(BlockId block) const noexcept;
    [[nodiscard]] ViscousFaceFluxField& field(BlockId block) const;

private:
    std::unordered_map<BlockId, ViscousFaceFluxField*> fields_;
};

class ViscousFaceFluxHaloExchanger {
public:
    ViscousFaceFluxHaloExchanger(const MpiRuntime& mpi,
                                 const ViscousFaceFluxHaloPlan& plan,
                                 int tag_base = 24576)
        : mpi_(mpi)
        , plan_(plan)
        , tag_base_(tag_base)
    {
        if (tag_base_ < 0) throw std::invalid_argument("face-flux MPI tag base is negative");
        prepare();
    }

    void exchange(const ViscousFaceFluxFieldRegistry& fields) const;

private:
    struct Pending {
        const FaceFluxExchangeDescriptor* descriptor = nullptr;
        std::vector<Real> values;
    };
    void prepare();

    const MpiRuntime& mpi_;
    const ViscousFaceFluxHaloPlan& plan_;
    int tag_base_ = 24576;
    mutable std::vector<Pending> receives_;
    mutable std::vector<Pending> sends_;
#if WCNS_HAS_MPI
    mutable std::vector<MPI_Request> requests_;
#endif
};

[[nodiscard]] TurbulenceViscousContribution evaluate_turbulence_viscous_face(
    const StructuredBlock& block,
    const ViscousFaceTrace& trace,
    const AlgorithmProfile& profile,
    const TransportModel& transport,
    const GasModel& gas,
    const ReferenceScales& reference,
    Axis axis,
    Index3 face,
    const BoundaryPatch* physical_boundary,
    const ITurbulenceModel* turbulence_model,
    const Field<Real>* les_closure = nullptr);

[[nodiscard]] ViscousFaceFluxField
compute_viscous_face_fluxes(const StructuredBlock& block,
                            const MetricField& metric,
                            const PrimitiveGradientField& gradients,
                            const AlgorithmProfile& profile,
                            const TransportModel& transport,
                            const BoundaryDataMap& boundary_data,
                            const GasModel& gas,
                            const ReferenceScales& reference,
                            const NumericalFloors& floors,
                            std::uint64_t version,
                            const ITurbulenceModel* turbulence_model = nullptr,
                            const Field<Real>* les_closure = nullptr);

void compute_viscous_face_fluxes_into(ViscousFaceFluxField& result,
                                      const StructuredBlock& block,
                                      const MetricField& metric,
                                      const PrimitiveGradientField& gradients,
                                      const AlgorithmProfile& profile,
                                      const TransportModel& transport,
                                      const BoundaryDataMap& boundary_data,
                                      const GasModel& gas,
                                      const ReferenceScales& reference,
                                      const NumericalFloors& floors,
                                      std::uint64_t version,
                                      const ITurbulenceModel* turbulence_model = nullptr,
                                      const Field<Real>* les_closure = nullptr);

void add_wcns_viscous_residual(StructuredBlock& block,
                               const MetricField& metric,
                               const ViscousFaceFluxField& flux,
                               const AlgorithmProfile& profile,
                               Real reynolds);

} // namespace wcns
