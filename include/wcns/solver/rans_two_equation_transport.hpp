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

// Six Cartesian model-gradient components followed by two cell-centred
// diffusion coefficients.  The complete workspace is exchanged before a
// partition-interface flux is assembled.
inline constexpr int two_equation_face_workspace_components = 8;

void compute_two_equation_gradients_and_source(
    Field<Real>& residual,
    Field<Real>& source_jacobian,
    Field<Real>& face_workspace,
    StructuredBlock& block,
    const MetricField& metric,
    const PrimitiveGradientField& mean_gradients,
    const ITurbulenceModel& model,
    const TransportModel& transport,
    const GasModel& gas,
    const ReferenceScales& reference);

void fill_two_equation_workspace_physical_ghosts(const StructuredBlock& block,
                                                 Field<Real>& face_workspace);

// The Y-stage two-equation path uses the already synchronized mean-flow mass
// flux and one synchronized second-order scalar flux on every interface.
void assemble_two_equation_flux_residual(Field<Real>& residual,
                                         const Field<Real>& face_workspace,
                                         const StructuredBlock& block,
                                         const MetricField& metric,
                                         const InviscidFaceFluxField& mean_flux,
                                         const ITurbulenceModel& model);

} // namespace wcns
