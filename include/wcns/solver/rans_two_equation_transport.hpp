#pragma once

#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/parallel/halo_exchanger.hpp>
#include <wcns/solver/inviscid_flux.hpp>
#include <wcns/solver/transport_model.hpp>
#include <wcns/solver/turbulence_model.hpp>
#include <wcns/solver/viscous_gradient.hpp>

#include <cstddef>
#include <limits>

namespace wcns {

struct WallFunctionDiagnostics {
    std::size_t face_count = 0;
    std::size_t out_of_range_count = 0;
    Real minimum_y_plus = std::numeric_limits<Real>::infinity();
    Real maximum_y_plus = 0.0;
};

void initialize_two_equation_fields(const MpiRuntime& mpi,
                                    LocalBlockSet& local_blocks,
                                    const TurbulenceModelConfig& config);

void synchronize_two_equation_fields(const HaloExchanger& exchanger,
                                     LocalBlockSet& local_blocks,
                                     const TurbulenceModelConfig& config,
                                     const TransportModel& transport,
                                     const ReferenceScales& reference);

[[nodiscard]] WallFunctionDiagnostics
global_wall_function_diagnostics(const MpiRuntime& mpi,
                                 const LocalBlockSet& local_blocks,
                                 const TurbulenceModelConfig& config);

// The Y-stage two-equation path uses one conservative face mass flux and a
// cell-centred second-order diffusive flux for every transported scalar.
void compute_two_equation_residual_and_source(
    Field<Real>& residual,
    Field<Real>& source_jacobian,
    Field<Real>& model_gradients,
    StructuredBlock& block,
    const MetricField& metric,
    const InviscidFaceFluxField& mean_flux,
    const PrimitiveGradientField& mean_gradients,
    const ITurbulenceModel& model,
    const TransportModel& transport,
    const GasModel& gas,
    const ReferenceScales& reference);

} // namespace wcns
