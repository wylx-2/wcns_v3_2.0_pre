#pragma once

#include <wcns/runtime/quantity_registry.hpp>
#include <wcns/runtime/simulation_driver.hpp>
#include <wcns/solver/viscous_halo.hpp>

#include <array>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace wcns {

struct BoundaryFacePhysics {
    Real pressure = 0.0;
    Real temperature = 0.0;
    Real viscosity = 0.0;
    Real pressure_coefficient = 0.0;
    Real skin_friction_coefficient = 0.0;
    Real heat_flux_into_wall = 0.0;
    Real wall_distance = std::numeric_limits<Real>::quiet_NaN();
    Real friction_velocity = std::numeric_limits<Real>::quiet_NaN();
    Real wall_y_plus = std::numeric_limits<Real>::quiet_NaN();
    Real wall_y_plus_class = std::numeric_limits<Real>::quiet_NaN();
    Real model_viscosity_ratio = 0.0;
    std::array<Real, 3> pressure_traction {};
    std::array<Real, 3> viscous_traction {};
    std::array<Real, 3> total_traction {};
};

struct BoundaryOutputScales {
    Real coordinate = 1.0;
    Real area = 1.0;
    Real pressure = 1.0;
    Real temperature = 1.0;
    Real viscosity = 1.0;
    Real velocity = 1.0;
    Real force = 1.0;
    Real moment = 1.0;
};

[[nodiscard]] BoundaryOutputScales boundary_output_scales(const QuantityContext& context,
                                                          int dimension);

// stress_normal is tau*n before the equation-level 1/Re factor.  The
// temperature gradient uses the same outward normal supplied here.
[[nodiscard]] BoundaryFacePhysics
evaluate_boundary_face_physics(Real pressure,
                               Real temperature,
                               Real viscosity,
                               const std::array<Real, 3>& outward_normal,
                               const std::array<Real, 3>& stress_normal,
                               Real thermal_coefficient,
                               const std::array<Real, 3>& temperature_gradient,
                               Real reynolds,
                               bool viscous,
                               const BoundaryOutputConfig& config);

// The optional prescribed y+ is the wall-function value used by the model.
// Without it, wall units are reconstructed from the resolved tangential stress.
void populate_boundary_wall_units(BoundaryFacePhysics& physics,
                                  const std::array<Real, 3>& outward_normal,
                                  Real density,
                                  Real wall_distance,
                                  Real reynolds,
                                  std::optional<Real> prescribed_y_plus = std::nullopt);

class BoundaryOutputWriter {
public:
    BoundaryOutputWriter(const MpiRuntime& mpi,
                         const CaseConfig& config,
                         const StructuredPartitionPlan& partition,
                         const LocalBlockSet& local_blocks,
                         const StructuredMesh& global_mesh,
                         const DistributedTopology& topology,
                         const BlockMetricMap& metrics,
                         const BlockBoundaryDataMap& boundary_data,
                         const GlobalConservationWeights& conservation_weights,
                         AlgorithmProfile profile,
                         QuantityContext quantities);

    [[nodiscard]] std::vector<std::string> write(const SimulationState& state);

private:
    const MpiRuntime& mpi_;
    const CaseConfig& config_;
    const StructuredPartitionPlan& partition_;
    const LocalBlockSet& local_blocks_;
    const StructuredMesh& global_mesh_;
    const DistributedTopology& topology_;
    const BlockMetricMap& metrics_;
    const BlockBoundaryDataMap& boundary_data_;
    const GlobalConservationWeights& conservation_weights_;
    AlgorithmProfile profile_;
    QuantityContext quantities_;
    std::uint64_t version_ = 0;
    bool load_history_created_ = false;
    std::vector<std::vector<Real>> load_history_;
};

} // namespace wcns
