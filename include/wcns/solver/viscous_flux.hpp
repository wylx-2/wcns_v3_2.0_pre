#pragma once

#include <wcns/solver/euler.hpp>
#include <wcns/solver/transport_model.hpp>
#include <wcns/solver/turbulence_model.hpp>

#include <array>

namespace wcns {

inline constexpr int viscous_primitive_components = 4;

enum class ViscousPrimitive : int {
    VelocityX = 0,
    VelocityY = 1,
    VelocityZ = 2,
    Temperature = 3,
};

using CartesianGradient = std::array<Real, 3>;
using PrimitiveGradients = std::array<CartesianGradient, viscous_primitive_components>;

struct ViscousFaceTrace {
    TemperaturePrimitiveState state {};
    PrimitiveGradients gradients {};
};

struct ViscousCartesianFlux {
    ConservativeState x {};
    ConservativeState y {};
    ConservativeState z {};
    Real viscosity = 0.0;
    Real thermal_coefficient = 0.0;
};

[[nodiscard]] ViscousCartesianFlux compute_viscous_cartesian_flux(const ViscousFaceTrace& trace,
                                                                  const TransportModel& transport,
                                                                  const GasModel& gas,
                                                                  const ReferenceScales& reference,
                                                                  const NumericalFloors& floors,
                                                                  int dimension);

[[nodiscard]] ViscousCartesianFlux compute_viscous_cartesian_flux(
    const ViscousFaceTrace& trace,
    const TransportModel& transport,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    const TurbulenceViscousContribution& turbulence,
    int dimension);

} // namespace wcns
