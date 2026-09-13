#include <wcns/solver/turbulence_model.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace wcns {
namespace {

class NoneTurbulenceModel final : public ITurbulenceModel {
public:
    explicit NoneTurbulenceModel(TurbulenceModelConfig config)
        : config_(std::move(config))
    {
        config_.validate();
        if (config_.kind != TurbulenceModelKind::None) {
            throw std::invalid_argument("none turbulence factory received a non-none model");
        }
    }

    const TurbulenceModelConfig& config() const noexcept override { return config_; }
    TurbulenceModelFamily family() const noexcept override { return TurbulenceModelFamily::None; }
    std::vector<TurbulenceFieldDescriptor> fields() const override { return {}; }

    TurbulenceViscousContribution
    viscous_contribution(const TurbulenceCellContext&) const override
    {
        return {};
    }

    TurbulenceSourceLinearization
    source_linearization(const TurbulenceCellContext&) const override
    {
        return {};
    }

private:
    TurbulenceModelConfig config_;
};

bool finite_stress(const SymmetricStress& stress)
{
    return std::isfinite(stress.xx) && std::isfinite(stress.yy) && std::isfinite(stress.zz)
        && std::isfinite(stress.xy) && std::isfinite(stress.xz) && std::isfinite(stress.yz);
}

} // namespace

TurbulenceModelKind turbulence_model_kind(const std::string& name)
{
    if (name == "none") return TurbulenceModelKind::None;
    if (name == "sa_neg") return TurbulenceModelKind::SaNegative;
    if (name == "k_omega_sst") return TurbulenceModelKind::KOmegaSst;
    if (name == "k_epsilon") return TurbulenceModelKind::KEpsilon;
    if (name == "smagorinsky") return TurbulenceModelKind::Smagorinsky;
    if (name == "scale_similarity") return TurbulenceModelKind::ScaleSimilarity;
    if (name == "mixed_smagorinsky_similarity") {
        return TurbulenceModelKind::MixedSmagorinskySimilarity;
    }
    if (name == "dynamic_smagorinsky") return TurbulenceModelKind::DynamicSmagorinsky;
    if (name == "wale") return TurbulenceModelKind::Wale;
    throw std::invalid_argument("unknown turbulence model: " + name);
}

const char* turbulence_model_name(TurbulenceModelKind kind)
{
    switch (kind) {
    case TurbulenceModelKind::None: return "none";
    case TurbulenceModelKind::SaNegative: return "sa_neg";
    case TurbulenceModelKind::KOmegaSst: return "k_omega_sst";
    case TurbulenceModelKind::KEpsilon: return "k_epsilon";
    case TurbulenceModelKind::Smagorinsky: return "smagorinsky";
    case TurbulenceModelKind::ScaleSimilarity: return "scale_similarity";
    case TurbulenceModelKind::MixedSmagorinskySimilarity:
        return "mixed_smagorinsky_similarity";
    case TurbulenceModelKind::DynamicSmagorinsky: return "dynamic_smagorinsky";
    case TurbulenceModelKind::Wale: return "wale";
    }
    throw std::invalid_argument("invalid turbulence model kind");
}

TurbulenceModelFamily turbulence_model_family(TurbulenceModelKind kind)
{
    switch (kind) {
    case TurbulenceModelKind::None: return TurbulenceModelFamily::None;
    case TurbulenceModelKind::SaNegative:
    case TurbulenceModelKind::KOmegaSst:
    case TurbulenceModelKind::KEpsilon: return TurbulenceModelFamily::RansTransport;
    case TurbulenceModelKind::Smagorinsky:
    case TurbulenceModelKind::ScaleSimilarity:
    case TurbulenceModelKind::MixedSmagorinskySimilarity:
    case TurbulenceModelKind::DynamicSmagorinsky:
    case TurbulenceModelKind::Wale: return TurbulenceModelFamily::LesAlgebraic;
    }
    throw std::invalid_argument("invalid turbulence model kind");
}

WallTreatment wall_treatment(const std::string& name)
{
    if (name == "resolved") return WallTreatment::Resolved;
    if (name == "wall_function") return WallTreatment::WallFunction;
    throw std::invalid_argument("unknown turbulence wall treatment: " + name);
}

const char* wall_treatment_name(WallTreatment treatment)
{
    switch (treatment) {
    case WallTreatment::Resolved: return "resolved";
    case WallTreatment::WallFunction: return "wall_function";
    }
    throw std::invalid_argument("invalid turbulence wall treatment");
}

void TurbulenceModelConfig::validate() const
{
    static_cast<void>(turbulence_model_name(kind));
    static_cast<void>(wall_treatment_name(wall_treatment));
    if (!std::isfinite(turbulent_prandtl) || turbulent_prandtl <= 0.0) {
        throw std::invalid_argument("turbulent Prandtl number must be positive and finite");
    }
    if (kind == TurbulenceModelKind::KEpsilon
        && wall_treatment != WallTreatment::WallFunction) {
        throw std::invalid_argument("k_epsilon requires wall_function treatment");
    }
}

std::string TurbulenceModelConfig::summary() const
{
    validate();
    std::ostringstream result;
    result << "turbulence(model=" << turbulence_model_name(kind)
           << ",family=";
    switch (turbulence_model_family(kind)) {
    case TurbulenceModelFamily::None: result << "none"; break;
    case TurbulenceModelFamily::RansTransport: result << "rans_transport"; break;
    case TurbulenceModelFamily::LesAlgebraic: result << "les_algebraic"; break;
    }
    result << ",wall_treatment=" << wall_treatment_name(wall_treatment)
           << ",Pr_t=" << std::setprecision(17) << turbulent_prandtl << ')';
    return result.str();
}

std::string TurbulenceModelConfig::restart_signature() const
{
    return "turbulence_v1;" + summary();
}

void TurbulenceViscousContribution::validate(int dimension) const
{
    if ((dimension != 2 && dimension != 3) || !finite_stress(stress)
        || !std::all_of(energy_heat_flux.begin(), energy_heat_flux.end(), [](Real value) {
               return std::isfinite(value);
           })
        || !std::isfinite(eddy_viscosity) || eddy_viscosity < 0.0) {
        throw std::invalid_argument("invalid turbulence viscous contribution");
    }
    if (dimension == 2
        && (stress.xz != 0.0 || stress.yz != 0.0 || energy_heat_flux[2] != 0.0)) {
        throw std::invalid_argument("two-dimensional turbulence contribution contains z flux");
    }
}

void TurbulenceSourceLinearization::validate(std::size_t variable_count) const
{
    if (source.size() != variable_count || jacobian.size() != variable_count * variable_count) {
        throw std::invalid_argument("turbulence source/Jacobian size mismatch");
    }
    if (!std::all_of(source.begin(), source.end(), [](Real value) { return std::isfinite(value); })
        || !std::all_of(jacobian.begin(), jacobian.end(), [](Real value) {
               return std::isfinite(value);
           })) {
        throw std::invalid_argument("turbulence source/Jacobian contains non-finite values");
    }
}

std::vector<Real>
TurbulenceSourceLinearization::implicit_diagonal_block(Real time_diagonal) const
{
    if (!std::isfinite(time_diagonal) || time_diagonal <= 0.0) {
        throw std::invalid_argument("implicit turbulence time diagonal must be positive and finite");
    }
    const auto variable_count = source.size();
    validate(variable_count);
    auto result = jacobian;
    for (auto& value : result)
        value = -value;
    for (std::size_t variable = 0; variable < variable_count; ++variable) {
        result[variable * variable_count + variable] += time_diagonal;
    }
    return result;
}

TurbulenceModelRegistry TurbulenceModelRegistry::create_builtin()
{
    TurbulenceModelRegistry result;
    result.register_model(TurbulenceModelKind::None, [](const TurbulenceModelConfig& config) {
        return std::make_unique<NoneTurbulenceModel>(config);
    });
    return result;
}

void TurbulenceModelRegistry::register_model(TurbulenceModelKind kind, Factory factory)
{
    if (!factory) throw std::invalid_argument("turbulence model factory is empty");
    if (!factories_.emplace(static_cast<int>(kind), std::move(factory)).second) {
        throw std::invalid_argument(std::string("duplicate turbulence model factory: ")
                                    + turbulence_model_name(kind));
    }
}

bool TurbulenceModelRegistry::contains(TurbulenceModelKind kind) const noexcept
{
    return factories_.find(static_cast<int>(kind)) != factories_.end();
}

std::unique_ptr<ITurbulenceModel>
TurbulenceModelRegistry::create(const TurbulenceModelConfig& config) const
{
    config.validate();
    const auto iterator = factories_.find(static_cast<int>(config.kind));
    if (iterator == factories_.end()) {
        throw std::invalid_argument(std::string("turbulence model is not registered: ")
                                    + turbulence_model_name(config.kind));
    }
    auto result = iterator->second(config);
    if (!result || result->config().kind != config.kind
        || result->family() != turbulence_model_family(config.kind)) {
        throw std::runtime_error("turbulence model factory returned inconsistent metadata");
    }
    return result;
}

} // namespace wcns
