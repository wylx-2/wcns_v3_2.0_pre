#pragma once

#include <wcns/core/field.hpp>
#include <wcns/core/types.hpp>
#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/mesh/structured_block.hpp>
#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/low_mach_preconditioner.hpp>

namespace wcns {

struct LuSgsIterationConfig {
    int sweeps = 1;
    Real relaxation = 1.0;

    void validate() const;
};

struct DualTimeIterationConfig {
    std::size_t max_iterations = 100;
    Real absolute_tolerance = 1.0e-10;
    Real relative_tolerance = 1.0e-8;
    Real cfl = 5.0;
    LuSgsIterationConfig lu_sgs {};

    void validate() const;
};

struct ScalarSpectralSystem {
    Field<Real> diagonal;
    Field<Real> coupling;
};

[[nodiscard]] ScalarSpectralSystem build_scalar_spectral_system(
    const StructuredBlock& block,
    const MetricField& metric,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    Real pseudo_cfl,
    Real physical_time_diagonal = 0.0,
    const WeissSmithParameters* preconditioner = nullptr,
    Real viscous_preconditioner_scale = 0.0,
    const Field<Real>* additional_face_coupling = nullptr);

// Adds the physical BDF term alpha*Gamma_0*Gamma_p^{-1} to a scalar-spectral
// pseudo-time diagonal.  The scalar input must not already contain alpha.
[[nodiscard]] Field<Real> build_preconditioned_time_diagonal_blocks(
    const StructuredBlock& block,
    const MetricField& metric,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    const Field<Real>& scalar_diagonal,
    Real physical_time_diagonal,
    const WeissSmithParameters& preconditioner,
    Real viscous_preconditioner_scale = 0.0);

// Coupling components are lower/upper I, lower/upper J, lower/upper K.
// The diagonal may have one component (shared scalar spectral radius) or one
// component per right-hand-side variable.
[[nodiscard]] Field<Real> solve_scalar_lu_sgs(const Field<Real>& right_hand_side,
                                              const Field<Real>& diagonal,
                                              const Field<Real>& coupling,
                                              int dimension,
                                              const LuSgsIterationConfig& config = {});

// diagonal_blocks stores a dense row-major block for every cell. Spatial
// neighbour coupling remains scalar-spectral and therefore multiplies the
// identity; source-term cross coupling is retained inside each dense block.
[[nodiscard]] Field<Real> solve_block_lu_sgs(const Field<Real>& right_hand_side,
                                             const Field<Real>& diagonal_blocks,
                                             const Field<Real>& coupling,
                                             int dimension,
                                             const LuSgsIterationConfig& config = {});

// Builds dense Euler A+/(-A-) neighbour blocks while retaining the same
// scalar spectral diagonal used by the baseline LU-SGS approximation.
[[nodiscard]] Field<Real> build_euler_face_coupling_blocks(
    const StructuredBlock& block,
    const MetricField& metric,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    const Field<Real>* additional_face_coupling = nullptr);

// Face blocks are row-major matrices ordered lower/upper I,J,K.  The diagonal
// remains scalar, so this captures Euler characteristic coupling without a
// per-cell dense factorization.
[[nodiscard]] Field<Real> solve_face_block_lu_sgs(
    const Field<Real>& right_hand_side,
    const Field<Real>& diagonal,
    const Field<Real>& face_blocks,
    int dimension,
    const LuSgsIterationConfig& config = {});

enum class BdfOrder {
    First = 1,
    Second = 2,
};

[[nodiscard]] Real bdf_time_diagonal(BdfOrder order, Real physical_time_step);
[[nodiscard]] Real bdf_time_residual(BdfOrder order,
                                     Real candidate,
                                     Real current,
                                     Real previous,
                                     Real physical_time_step);

} // namespace wcns
