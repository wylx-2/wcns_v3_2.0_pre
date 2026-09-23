#pragma once

#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/euler.hpp>

#include <array>

namespace wcns {

using Matrix5 = std::array<std::array<Real, 5>, 5>;

struct WeissSmithParameters {
    Real mach_cutoff = 1.0e-3;
    Real viscous_cutoff = 1.0;

    void validate() const;
};

struct WeissSmithState {
    Matrix5 gamma {};
    Matrix5 inverse {};
    std::array<Real, 5> eigenvalues {};
    Real reference_speed = 0.0;
    Real beta = 1.0;
};

[[nodiscard]] Real weiss_smith_reference_speed(Real velocity_magnitude,
                                                Real sound_speed,
                                                Real viscous_speed,
                                                const WeissSmithParameters& parameters);

[[nodiscard]] std::array<Real, 5> weiss_smith_eigenvalues(Real normal_velocity,
                                                          Real sound_speed,
                                                          Real beta);

[[nodiscard]] WeissSmithState weiss_smith_state(const TemperaturePrimitiveState& state,
                                                Normal3 unit_normal,
                                                const GasModel& gas,
                                                const ReferenceScales& reference,
                                                const NumericalFloors& floors,
                                                const WeissSmithParameters& parameters,
                                                Real viscous_speed = 0.0,
                                                int dimension = 3);

// Weiss-Smith matrix expressed in W=(rho,u,v,w,p).  Unlike the temperature-
// primitive form above, this form only needs the pressure state used by a
// Riemann solver and is therefore suitable for preconditioned Roe dissipation.
[[nodiscard]] WeissSmithState
weiss_smith_pressure_state(const PressurePrimitiveState& state,
                           Normal3 unit_normal,
                           const GasModel& gas,
                           const NumericalFloors& floors,
                           const WeissSmithParameters& parameters,
                           Real viscous_speed = 0.0,
                           int dimension = 3);

// Left preconditioning for a conservative residual.  Gamma_0 is the physical
// dQ/dW matrix and Gamma_p is the Weiss-Smith pseudo-time matrix, so the
// conservative update is driven by Gamma_0 Gamma_p^{-1} R(Q).
[[nodiscard]] Matrix5
weiss_smith_conservative_preconditioner(const PressurePrimitiveState& state,
                                        const GasModel& gas,
                                        const NumericalFloors& floors,
                                        const WeissSmithParameters& parameters,
                                        Real viscous_speed = 0.0,
                                        int dimension = 3);

[[nodiscard]] ConservativeState
weiss_smith_precondition_residual(const PressurePrimitiveState& state,
                                  const ConservativeState& residual,
                                  const GasModel& gas,
                                  const NumericalFloors& floors,
                                  const WeissSmithParameters& parameters,
                                  Real viscous_speed = 0.0,
                                  int dimension = 3);

[[nodiscard]] Matrix5 multiply(const Matrix5& left, const Matrix5& right);

} // namespace wcns
