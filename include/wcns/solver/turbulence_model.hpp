#pragma once

#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/euler.hpp>
#include <wcns/solver/turbulence_fields.hpp>

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace wcns {

enum class TurbulenceModelKind {
    None,
    SaNegative,
    KOmegaSst,
    KEpsilon,
    Smagorinsky,
    ScaleSimilarity,
    MixedSmagorinskySimilarity,
    DynamicSmagorinsky,
    Wale,
};

enum class TurbulenceModelFamily {
    None,
    RansTransport,
    LesAlgebraic,
};

enum class WallTreatment {
    Resolved,
    WallFunction,
};

enum class TurbulenceSourceTreatment {
    Explicit,
    LocalImplicit,
};

struct TurbulenceModelConfig {
    TurbulenceModelKind kind = TurbulenceModelKind::None;
    WallTreatment wall_treatment = WallTreatment::Resolved;
    bool experimental = false;
    Real turbulent_prandtl = 0.9;
    Real sa_farfield_nu_tilde_ratio = 3.0;
    Real freestream_turbulence_intensity = 0.01;
    Real freestream_length_scale = 0.1;
    Real model_floor = 1.0e-12;
    Real wall_function_y_plus_min = 30.0;
    Real wall_function_y_plus_max = 300.0;
    TurbulenceSourceTreatment source_treatment = TurbulenceSourceTreatment::Explicit;

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

struct SymmetricStress {
    Real xx = 0.0;
    Real yy = 0.0;
    Real zz = 0.0;
    Real xy = 0.0;
    Real xz = 0.0;
    Real yz = 0.0;
};

// Stress uses the same positive diffusive-flux convention as viscous_flux.
// energy_heat_flux is added directly to the three energy-flux components.
struct TurbulenceViscousContribution {
    SymmetricStress stress;
    std::array<Real, 3> energy_heat_flux {{0.0, 0.0, 0.0}};
    Real eddy_viscosity = 0.0;

    void validate(int dimension) const;
};

struct TurbulenceSourceLinearization {
    std::vector<Real> source;
    std::vector<Real> jacobian;

    void validate(std::size_t variable_count) const;
    [[nodiscard]] std::vector<Real> implicit_diagonal_block(Real time_diagonal) const;
};

struct TurbulenceCellContext {
    TemperaturePrimitiveState mean_state {};
    std::array<std::array<Real, 3>, 4> primitive_gradients {};
    std::vector<Real> model_values;
    std::vector<std::array<Real, 3>> model_gradients;
    Real molecular_kinematic_viscosity = 0.0;
    Real wall_distance = 0.0;
    Real filter_width = 0.0;
    Real reference_reynolds = 0.0;
    Real reference_mach = 0.0;
    Real heat_capacity_ratio = 0.0;
    int dimension = 3;
};

struct SaNegativeConstants {
    Real cb1 = 0.1355;
    Real sigma = 2.0 / 3.0;
    Real cb2 = 0.622;
    Real kappa = 0.41;
    Real cw2 = 0.3;
    Real cw3 = 2.0;
    Real cv1 = 7.1;
    Real ct3 = 1.2;
    Real ct4 = 0.5;
    Real cn1 = 16.0;
    Real c2 = 0.7;
    Real c3 = 0.9;

    [[nodiscard]] Real cw1() const noexcept;
    void validate() const;
};

struct SaNegativeEvaluation {
    Real chi = 0.0;
    Real fv1 = 0.0;
    Real fv2 = 0.0;
    Real ft2 = 0.0;
    Real fn = 1.0;
    Real vorticity_magnitude = 0.0;
    Real modified_vorticity = 0.0;
    Real r = 0.0;
    Real fw = 0.0;
    Real diffusion_coefficient = 0.0;
    Real eddy_kinematic_viscosity = 0.0;
    Real production_source = 0.0;
    Real destruction_source = 0.0;
    Real cross_diffusion_source = 0.0;
    Real source = 0.0;
    Real source_derivative = 0.0;
    bool negative_branch = false;
};

[[nodiscard]] SaNegativeEvaluation
evaluate_sa_negative(const TurbulenceCellContext& context,
                     const SaNegativeConstants& constants = {});

[[nodiscard]] Real sa_negative_farfield_value(const TurbulenceModelConfig& config,
                                              Real reference_reynolds);

[[nodiscard]] Real local_implicit_turbulence_increment(Real explicit_residual,
                                                       Real source_jacobian,
                                                       Real time_step);

class ITurbulenceModel {
public:
    virtual ~ITurbulenceModel() = default;
    [[nodiscard]] virtual const TurbulenceModelConfig& config() const noexcept = 0;
    [[nodiscard]] virtual TurbulenceModelFamily family() const noexcept = 0;
    [[nodiscard]] virtual std::vector<TurbulenceFieldDescriptor> fields() const = 0;
    [[nodiscard]] virtual TurbulenceViscousContribution
    viscous_contribution(const TurbulenceCellContext& context) const = 0;
    [[nodiscard]] virtual std::vector<Real>
    diffusion_coefficients(const TurbulenceCellContext& context) const = 0;
    [[nodiscard]] virtual TurbulenceSourceLinearization
    source_linearization(const TurbulenceCellContext& context) const = 0;
};

class TurbulenceModelRegistry {
public:
    using Factory
        = std::function<std::unique_ptr<ITurbulenceModel>(const TurbulenceModelConfig&)>;

    [[nodiscard]] static TurbulenceModelRegistry create_builtin();
    void register_model(TurbulenceModelKind kind, Factory factory);
    [[nodiscard]] bool contains(TurbulenceModelKind kind) const noexcept;
    [[nodiscard]] std::unique_ptr<ITurbulenceModel>
    create(const TurbulenceModelConfig& config) const;

private:
    std::unordered_map<int, Factory> factories_;
};

[[nodiscard]] TurbulenceModelKind turbulence_model_kind(const std::string& name);
[[nodiscard]] const char* turbulence_model_name(TurbulenceModelKind kind);
[[nodiscard]] TurbulenceModelFamily turbulence_model_family(TurbulenceModelKind kind);
[[nodiscard]] WallTreatment wall_treatment(const std::string& name);
[[nodiscard]] const char* wall_treatment_name(WallTreatment treatment);
[[nodiscard]] TurbulenceSourceTreatment turbulence_source_treatment(const std::string& name);
[[nodiscard]] const char* turbulence_source_treatment_name(TurbulenceSourceTreatment treatment);

} // namespace wcns
