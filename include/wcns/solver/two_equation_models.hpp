#pragma once

#include <wcns/solver/turbulence_model.hpp>

#include <array>

namespace wcns {

struct SstConstants {
    Real sigma_k1 = 0.85;
    Real sigma_k2 = 1.0;
    Real sigma_omega1 = 0.5;
    Real sigma_omega2 = 0.856;
    Real beta1 = 0.075;
    Real beta2 = 0.0828;
    Real beta_star = 0.09;
    Real kappa = 0.41;
    Real a1 = 0.31;
    Real gamma1 = 5.0 / 9.0;
    Real gamma2 = 0.44;
    Real production_limiter = 10.0;
    Real cross_diffusion_floor = 1.0e-10;

    void validate() const;
};

struct StandardKEpsilonConstants {
    Real c_mu = 0.09;
    Real c_epsilon1 = 1.44;
    Real c_epsilon2 = 1.92;
    Real sigma_k = 1.0;
    Real sigma_epsilon = 1.3;
    Real kappa = 0.41;

    void validate() const;
};

struct TwoEquationEvaluation {
    Real blending_f1 = 0.0;
    Real blending_f2 = 0.0;
    Real production = 0.0;
    Real destruction_k = 0.0;
    Real destruction_second = 0.0;
    Real cross_diffusion = 0.0;
    Real eddy_kinematic_viscosity = 0.0;
    std::array<Real, 2> diffusion_coefficients {{0.0, 0.0}};
    std::array<Real, 2> source {{0.0, 0.0}};
    std::array<Real, 4> source_jacobian {{0.0, 0.0, 0.0, 0.0}};
};

[[nodiscard]] TwoEquationEvaluation
evaluate_k_omega_sst(const TurbulenceCellContext& context,
                     const SstConstants& constants = {});

[[nodiscard]] TwoEquationEvaluation
evaluate_standard_k_epsilon(const TurbulenceCellContext& context,
                            const StandardKEpsilonConstants& constants = {});

// Values are nondimensionalized by U_ref^2, U_ref/L_ref, and U_ref^3/L_ref.
[[nodiscard]] std::array<Real, 2>
two_equation_farfield_values(const TurbulenceModelConfig& config);

} // namespace wcns
