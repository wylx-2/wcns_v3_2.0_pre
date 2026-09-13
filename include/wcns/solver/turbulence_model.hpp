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

struct TurbulenceModelConfig {
    TurbulenceModelKind kind = TurbulenceModelKind::None;
    WallTreatment wall_treatment = WallTreatment::Resolved;
    Real turbulent_prandtl = 0.9;

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
    Real wall_distance = 0.0;
    Real filter_width = 0.0;
    int dimension = 3;
};

class ITurbulenceModel {
public:
    virtual ~ITurbulenceModel() = default;
    [[nodiscard]] virtual const TurbulenceModelConfig& config() const noexcept = 0;
    [[nodiscard]] virtual TurbulenceModelFamily family() const noexcept = 0;
    [[nodiscard]] virtual std::vector<TurbulenceFieldDescriptor> fields() const = 0;
    [[nodiscard]] virtual TurbulenceViscousContribution
    viscous_contribution(const TurbulenceCellContext& context) const = 0;
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

} // namespace wcns
