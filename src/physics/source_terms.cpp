#include <wcns/physics/source_terms.hpp>
#include <wcns/physics/chapter5.hpp>

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <unordered_set>

namespace wcns {
namespace {

void validate_kind(SourceModelKind kind)
{
    switch (kind) {
    case SourceModelKind::UniformConservative:
    case SourceModelKind::BodyForce:
    case SourceModelKind::PressureGradient:
    case SourceModelKind::RampTrip:
    case SourceModelKind::ManufacturedSolution: return;
    }
    throw std::invalid_argument("source-term configuration contains an unknown model");
}

template <std::size_t N> bool all_zero(const std::array<Real, N>& values)
{
    return std::all_of(values.begin(), values.end(), [](Real value) { return value == 0.0; });
}

template <std::size_t N> void validate_finite(const std::array<Real, N>& values, const char* label)
{
    if (!std::all_of(
            values.begin(), values.end(), [](Real value) { return std::isfinite(value); })) {
        throw std::invalid_argument(std::string(label) + " must be finite");
    }
}

bool contains_model(const SourceTermConfig& config, SourceModelKind kind)
{
    return std::find(config.models.begin(), config.models.end(), kind) != config.models.end();
}

template <std::size_t N> std::string values_string(const std::array<Real, N>& values)
{
    std::ostringstream stream;
    stream << std::setprecision(std::numeric_limits<Real>::max_digits10);
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) stream << ',';
        stream << values[index];
    }
    return stream.str();
}

} // namespace

const char* source_model_name(SourceModelKind kind)
{
    validate_kind(kind);
    switch (kind) {
    case SourceModelKind::UniformConservative: return "uniform_conservative";
    case SourceModelKind::BodyForce: return "body_force";
    case SourceModelKind::PressureGradient: return "pressure_gradient";
    case SourceModelKind::ManufacturedSolution: return "manufactured_solution";
    case SourceModelKind::RampTrip: return "ramp_trip";
    }
    throw std::logic_error("unreachable source model");
}

void SourceTermConfig::validate() const
{
    if(!std::isfinite(ramp_trip_amplitude) || ramp_trip_amplitude<0
        || !std::isfinite(ramp_trip_span) || ramp_trip_span<=0
        || (ramp_trip_amplitude!=0 && (!enable_source_terms || !contains_model(*this,SourceModelKind::RampTrip))))
        throw std::invalid_argument("invalid ramp trip configuration");
    validate_finite(uniform_conservative, "uniform conservative source");
    validate_finite(body_acceleration, "body acceleration");
    validate_finite(pressure_gradient, "pressure-gradient force");
    validate_finite(manufactured_amplitude, "manufactured amplitude");
    if (!enable_source_terms) {
        if (!models.empty() || !all_zero(uniform_conservative) || !all_zero(body_acceleration)
            || !all_zero(pressure_gradient) || !all_zero(manufactured_amplitude)) {
            throw std::invalid_argument(
                "disabled source terms require empty models and zero parameters");
        }
        return;
    }
    if (models.empty()) {
        throw std::invalid_argument("enabled source terms require at least one model");
    }

    std::unordered_set<int> seen;
    for (const auto model : models) {
        validate_kind(model);
        if (!seen.insert(static_cast<int>(model)).second) {
            throw std::invalid_argument("source-term configuration contains a duplicate model");
        }
    }
    if (!contains_model(*this, SourceModelKind::UniformConservative)
        && !all_zero(uniform_conservative)) {
        throw std::invalid_argument("uniform source parameters require the uniform model");
    }
    if (!contains_model(*this, SourceModelKind::BodyForce) && !all_zero(body_acceleration)) {
        throw std::invalid_argument("body acceleration requires the body-force model");
    }
    if (!contains_model(*this, SourceModelKind::PressureGradient) && !all_zero(pressure_gradient)) {
        throw std::invalid_argument("pressure-gradient force requires the pressure-gradient model");
    }
    if (!contains_model(*this, SourceModelKind::ManufacturedSolution)
        && !all_zero(manufactured_amplitude)) {
        throw std::invalid_argument("manufactured parameters require the manufactured model");
    }
}

std::string SourceTermConfig::summary() const
{
    validate();
    std::string result = "enable_source_terms=";
    result += enable_source_terms ? "true;models=" : "false;models=";
    for (std::size_t index = 0; index < models.size(); ++index) {
        if (index != 0) {
            result += ',';
        }
        result += source_model_name(models[index]);
    }
    if (enable_source_terms) {
        result += ";uniform=" + values_string(uniform_conservative);
        result += ";body_acceleration=" + values_string(body_acceleration);
        result += ";pressure_gradient=" + values_string(pressure_gradient);
        result += ";manufactured_amplitude=" + values_string(manufactured_amplitude);
    }
    if(contains_model(*this,SourceModelKind::RampTrip))
        result += ";ramp_trip="+values_string(std::array<Real,2>{{ramp_trip_amplitude,ramp_trip_span}});
    return result;
}

std::string SourceTermConfig::restart_signature() const
{
    return "source_terms_v3;" + summary();
}

SourceTermRegistry SourceTermRegistry::create_stage_h(const SourceTermConfig& config)
{
    config.validate();
    if (config.enable_source_terms) {
        throw std::logic_error("source-term models cannot be instantiated before stage J");
    }
    return {};
}

SourceTermRegistry SourceTermRegistry::create_stage_j(const SourceTermConfig& config)
{
    config.validate();
    SourceTermRegistry result;
    if (config.enable_source_terms) {
        result.model_count_ = config.models.size();
        result.config_ = config;
    }
    return result;
}

void SourceTermRegistry::set_pressure_gradient_x(Real force)
{
    if (!std::isfinite(force) || !contains_model(config_, SourceModelKind::PressureGradient))
        throw std::invalid_argument("dynamic pressure gradient requires an enabled pressure-gradient model and finite force");
    config_.pressure_gradient[0] = force;
}

std::array<Real, 5> SourceTermRegistry::evaluate(const std::array<Real, 5>& conservative,
                                                 const std::array<Real, 3>& coordinates,
                                                 Real time,
                                                 int dimension) const
{
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("source evaluation dimension must be two or three");
    }
    validate_finite(conservative, "source conservative state");
    validate_finite(coordinates, "source coordinates");
    if (!std::isfinite(time) || conservative[0] <= 0.0) {
        throw std::invalid_argument("source time and density must be valid");
    }
    std::array<Real, 5> result {};
    for (const auto model : config_.models) {
        switch (model) {
        case SourceModelKind::UniformConservative:
            for (std::size_t component = 0; component < result.size(); ++component) {
                result[component] += config_.uniform_conservative[component];
            }
            break;
        case SourceModelKind::BodyForce:
            result[1] += conservative[0] * config_.body_acceleration[0];
            result[2] += conservative[0] * config_.body_acceleration[1];
            result[3] += conservative[0] * config_.body_acceleration[2];
            result[4] += conservative[1] * config_.body_acceleration[0]
                + conservative[2] * config_.body_acceleration[1]
                + conservative[3] * config_.body_acceleration[2];
            break;
        case SourceModelKind::PressureGradient:
            result[1] += config_.pressure_gradient[0];
            result[2] += config_.pressure_gradient[1];
            result[3] += config_.pressure_gradient[2];
            result[4] += (conservative[1] * config_.pressure_gradient[0]
                          + conservative[2] * config_.pressure_gradient[1]
                          + conservative[3] * config_.pressure_gradient[2])
                / conservative[0];
            break;
        case SourceModelKind::RampTrip: {
            if(dimension!=3) throw std::invalid_argument("ramp trip requires 3D");
            const Real force=ramp_trip(coordinates[0],coordinates[1],coordinates[2],time,
                config_.ramp_trip_amplitude,config_.ramp_trip_span);
            result[1]+=force; result[4]+=force*conservative[1]/conservative[0];
            break;
        }
        case SourceModelKind::ManufacturedSolution: {
            const Real shape = 1.0 + coordinates[0] + coordinates[1]
                + (dimension == 3 ? coordinates[2] : 0.0) + time;
            for (std::size_t component = 0; component < result.size(); ++component) {
                result[component] += config_.manufactured_amplitude[component] * shape;
            }
            break;
        }
        }
    }
    if (dimension == 2 && result[3] != 0.0) {
        throw std::invalid_argument("two-dimensional source must have zero z momentum");
    }
    validate_finite(result, "evaluated source");
    return result;
}

} // namespace wcns
