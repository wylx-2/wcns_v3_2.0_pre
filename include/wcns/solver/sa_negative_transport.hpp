#pragma once

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/parallel/halo_exchanger.hpp>
#include <wcns/solver/inviscid_flux.hpp>
#include <wcns/solver/turbulence_model.hpp>
#include <wcns/solver/viscous_gradient.hpp>
#include <wcns/solver/viscous_operator.hpp>

#include <unordered_map>

namespace wcns {

// Allocates the descriptor-defined SA fields, initializes the transported
// working variable from the configured farfield ratio, and computes d from
// the global set of resolved no-slip wall primitives.
void initialize_sa_negative_fields(const MpiRuntime& mpi,
                                   LocalBlockSet& local_blocks,
                                   const TurbulenceModelConfig& config,
                                   Real reference_reynolds);

// Applies connectivity exchange followed by model-specific physical boundary
// values. Only physical face slabs are written; edge and corner ghosts remain
// outside the model transport contract.
void synchronize_sa_negative_fields(const HaloExchanger& exchanger,
                                    LocalBlockSet& local_blocks,
                                    Real farfield_nu_tilde);

// High-order conservative-gradient construction for the SA working variable.
// Temperature slots of PrimitiveGradientField carry grad(nu_tilde); the remaining
// slots are zero so the existing vector/tensor-aware gradient halo path can be
// reused without adding a second MPI protocol.
void compute_sa_negative_gradient(PrimitiveGradientField& result,
                                  const StructuredBlock& block,
                                  const MetricField& metric,
                                  const AlgorithmProfile& profile,
                                  std::uint64_t version);

// Forms the conservative SA face flux in component zero of a standard viscous
// face-flux workspace and updates cell source diagnostics/Jacobian.
void compute_sa_negative_flux_and_source(
    ViscousFaceFluxField& flux,
    Field<Real>& residual,
    Field<Real>& source_jacobian,
    StructuredBlock& block,
    const MetricField& metric,
    const InviscidFaceFluxField& mean_flux,
    const PrimitiveGradientField& mean_gradients,
    const PrimitiveGradientField& model_gradients,
    const AlgorithmProfile& profile,
    const ITurbulenceModel& model,
    const TransportModel& transport,
    const GasModel& gas,
    const ReferenceScales& reference,
    std::uint64_t version);

// Adds -div(F_sa)+S_sa to the conservative rho*nu_tilde residual.
void assemble_sa_negative_residual(Field<Real>& residual,
                                   const StructuredBlock& block,
                                   const MetricField& metric,
                                   const ViscousFaceFluxField& flux,
                                   const AlgorithmProfile& profile);

} // namespace wcns
