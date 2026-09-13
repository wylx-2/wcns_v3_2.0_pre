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

    std::vector<Real> diffusion_coefficients(const TurbulenceCellContext&) const override
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

class SaNegativeTurbulenceModel final : public ITurbulenceModel {
public:
    explicit SaNegativeTurbulenceModel(TurbulenceModelConfig config)
        : config_(std::move(config))
    {
        config_.validate();
        if (config_.kind != TurbulenceModelKind::SaNegative) {
            throw std::invalid_argument("SA-neg factory received a different model");
        }
    }

    const TurbulenceModelConfig& config() const noexcept override { return config_; }
    TurbulenceModelFamily family() const noexcept override
    {
        return TurbulenceModelFamily::RansTransport;
    }
    std::vector<TurbulenceFieldDescriptor> fields() const override
    {
        const Real any_finite = std::numeric_limits<Real>::lowest();
        return {
            {"nu_tilde",
             TurbulenceFieldRole::Transported,
             TurbulenceFieldScale::KinematicViscosity,
             false,
             any_finite},
            {"mu_t_over_mu",
             TurbulenceFieldRole::Diagnostic,
             TurbulenceFieldScale::Dimensionless,
             false,
             0.0},
            {"sa_production",
             TurbulenceFieldRole::Diagnostic,
             TurbulenceFieldScale::Dissipation,
             false,
             any_finite},
            {"sa_destruction",
             TurbulenceFieldRole::Diagnostic,
             TurbulenceFieldScale::Dissipation,
             false,
             any_finite},
            {"wall_distance",
             TurbulenceFieldRole::Auxiliary,
             TurbulenceFieldScale::Length,
             true,
             0.0},
            {"sa_negative_branch",
             TurbulenceFieldRole::Diagnostic,
             TurbulenceFieldScale::Dimensionless,
             false,
             0.0},
        };
    }

    TurbulenceViscousContribution
    viscous_contribution(const TurbulenceCellContext& context) const override
    {
        if ((context.dimension != 2 && context.dimension != 3)
            || context.model_values.size() != 1
            || !std::isfinite(context.model_values.front())
            || !std::isfinite(context.molecular_kinematic_viscosity)
            || context.molecular_kinematic_viscosity <= 0.0
            || !std::isfinite(context.mean_state[temperature_density])
            || context.mean_state[temperature_density] <= 0.0) {
            throw std::invalid_argument("SA-neg viscous context is invalid");
        }
        const Real nu_tilde = context.model_values.front();
        Real eddy_kinematic = 0.0;
        if (nu_tilde >= 0.0) {
            const Real chi = nu_tilde / context.molecular_kinematic_viscosity;
            const Real chi3 = chi * chi * chi;
            const Real cv13 = 7.1 * 7.1 * 7.1;
            eddy_kinematic = nu_tilde * chi3 / (chi3 + cv13);
        }
        if (!std::isfinite(eddy_kinematic) || eddy_kinematic < 0.0
            || !std::isfinite(context.reference_reynolds)
            || context.reference_reynolds <= 0.0) {
            throw PhysicsError("SA-neg eddy viscosity is invalid");
        }
        const Real mu_t = context.mean_state[temperature_density] * context.reference_reynolds
            * eddy_kinematic;
        const auto& du = context.primitive_gradients[0];
        const auto& dv = context.primitive_gradients[1];
        const auto& dw = context.primitive_gradients[2];
        const Real divergence = du[0] + dv[1] + dw[2];
        TurbulenceViscousContribution result;
        result.eddy_viscosity = mu_t;
        result.stress.xx = 2.0 * mu_t * (du[0] - divergence / 3.0);
        result.stress.yy = 2.0 * mu_t * (dv[1] - divergence / 3.0);
        result.stress.zz = 2.0 * mu_t * (dw[2] - divergence / 3.0);
        result.stress.xy = mu_t * (du[1] + dv[0]);
        result.stress.xz = mu_t * (du[2] + dw[0]);
        result.stress.yz = mu_t * (dv[2] + dw[1]);
        if (!std::isfinite(context.reference_mach) || context.reference_mach <= 0.0
            || !std::isfinite(context.heat_capacity_ratio)
            || context.heat_capacity_ratio <= 1.0) {
            throw std::invalid_argument("SA-neg thermal context is invalid");
        }
        const Real heat = mu_t
            / ((context.heat_capacity_ratio - 1.0) * context.reference_mach
               * context.reference_mach * config_.turbulent_prandtl);
        const auto& temperature_gradient = context.primitive_gradients[3];
        for (int direction = 0; direction < context.dimension; ++direction) {
            result.energy_heat_flux[static_cast<std::size_t>(direction)]
                = heat * temperature_gradient[static_cast<std::size_t>(direction)];
        }
        result.validate(context.dimension);
        return result;
    }

    std::vector<Real> diffusion_coefficients(const TurbulenceCellContext& context) const override
    {
        const auto evaluation = evaluate_sa_negative(context);
        return {context.mean_state[temperature_density] * evaluation.diffusion_coefficient};
    }

    TurbulenceSourceLinearization
    source_linearization(const TurbulenceCellContext& context) const override
    {
        const auto evaluation = evaluate_sa_negative(context);
        const Real density_value = context.mean_state[temperature_density];
        TurbulenceSourceLinearization result {
            {density_value * evaluation.source},
            {evaluation.source_derivative},
        };
        result.validate(1);
        return result;
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

TurbulenceSourceTreatment turbulence_source_treatment(const std::string& name)
{
    if (name == "explicit") return TurbulenceSourceTreatment::Explicit;
    if (name == "local_implicit") return TurbulenceSourceTreatment::LocalImplicit;
    throw std::invalid_argument("unknown turbulence source treatment: " + name);
}

const char* turbulence_source_treatment_name(TurbulenceSourceTreatment treatment)
{
    switch (treatment) {
    case TurbulenceSourceTreatment::Explicit: return "explicit";
    case TurbulenceSourceTreatment::LocalImplicit: return "local_implicit";
    }
    throw std::invalid_argument("invalid turbulence source treatment");
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
    static_cast<void>(turbulence_source_treatment_name(source_treatment));
    if (!std::isfinite(sa_farfield_nu_tilde_ratio)
        || sa_farfield_nu_tilde_ratio < 3.0 || sa_farfield_nu_tilde_ratio > 5.0) {
        throw std::invalid_argument("SA-neg farfield nu-tilde ratio must lie in [3,5]");
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
           << ",Pr_t=" << std::setprecision(17) << turbulent_prandtl;
    if (kind == TurbulenceModelKind::SaNegative) {
        result << ",farfield_nu_tilde_over_nu=" << sa_farfield_nu_tilde_ratio
               << ",source=" << turbulence_source_treatment_name(source_treatment);
    }
    result << ')';
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
    result.register_model(TurbulenceModelKind::SaNegative,
                          [](const TurbulenceModelConfig& config) {
                              return std::make_unique<SaNegativeTurbulenceModel>(config);
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
