#include <wcns/solver/turbulence_fields.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace wcns {
namespace {

const char* role_name(TurbulenceFieldRole role)
{
    switch (role) {
    case TurbulenceFieldRole::Transported: return "transported";
    case TurbulenceFieldRole::Auxiliary: return "auxiliary";
    case TurbulenceFieldRole::Diagnostic: return "diagnostic";
    }
    throw std::invalid_argument("invalid turbulence field role");
}

const char* scale_name(TurbulenceFieldScale scale)
{
    switch (scale) {
    case TurbulenceFieldScale::Dimensionless: return "dimensionless";
    case TurbulenceFieldScale::VelocitySquared: return "velocity_squared";
    case TurbulenceFieldScale::KinematicViscosity: return "kinematic_viscosity";
    case TurbulenceFieldScale::DynamicViscosity: return "dynamic_viscosity";
    case TurbulenceFieldScale::Pressure: return "pressure";
    case TurbulenceFieldScale::InverseTime: return "inverse_time";
    case TurbulenceFieldScale::Dissipation: return "dissipation";
    case TurbulenceFieldScale::Length: return "length";
    }
    throw std::invalid_argument("invalid turbulence field scale");
}

} // namespace

void TurbulenceFieldDescriptor::validate() const
{
    if (name.empty()) {
        throw std::invalid_argument("turbulence field name must not be empty");
    }
    if (!std::all_of(name.begin(), name.end(), [](unsigned char character) {
            return (character >= 'a' && character <= 'z')
                || (character >= '0' && character <= '9') || character == '_';
        })) {
        throw std::invalid_argument(
            "turbulence field name must contain only lowercase ASCII, digits, or underscore");
    }
    if (!std::isfinite(lower_bound)) {
        throw std::invalid_argument("turbulence field lower bound must be finite");
    }
}

std::string TurbulenceFieldDescriptor::signature() const
{
    validate();
    std::ostringstream result;
    result << name << ':' << role_name(role) << ':' << scale_name(scale) << ':'
           << (strictly_positive ? "positive" : "bounded") << ':' << std::setprecision(17)
           << lower_bound;
    return result.str();
}

TurbulenceFloorProjection
project_positive_turbulence_state(Real conservative,
                                  Real density,
                                  const TurbulenceFieldDescriptor& descriptor)
{
    descriptor.validate();
    if (!descriptor.strictly_positive) {
        throw std::invalid_argument(
            "positive turbulence projection requires a strictly-positive descriptor");
    }
    if (!std::isfinite(conservative) || !std::isfinite(density) || density <= 0.0) {
        throw std::invalid_argument(
            "positive turbulence projection requires finite conservative state and density");
    }
    const Real specific = conservative / density;
    if (!std::isfinite(specific)) {
        throw std::invalid_argument(
            "positive turbulence projection produced a non-finite specific state");
    }
    if (specific > descriptor.lower_bound) {
        return {conservative, specific, false};
    }
    const Real projected
        = std::nextafter(descriptor.lower_bound, std::numeric_limits<Real>::infinity());
    return {density * projected, projected, true};
}

TurbulenceFieldSet::TurbulenceFieldSet(
    Extent3 extent,
    int ghost_width,
    std::vector<TurbulenceFieldDescriptor> descriptors)
{
    reset(extent, ghost_width, std::move(descriptors));
}

void TurbulenceFieldSet::reset(Extent3 extent,
                               int ghost_width,
                               std::vector<TurbulenceFieldDescriptor> descriptors)
{
    if (descriptors.empty()) {
        clear();
        return;
    }
    std::unordered_map<std::string, int> components;
    components.reserve(descriptors.size());
    for (std::size_t index = 0; index < descriptors.size(); ++index) {
        descriptors[index].validate();
        if (!components.emplace(descriptors[index].name, static_cast<int>(index)).second) {
            throw std::invalid_argument("duplicate turbulence field descriptor: "
                                        + descriptors[index].name);
        }
    }
    Field<Real> storage(extent,
                        static_cast<int>(descriptors.size()),
                        ghost_width,
                        std::numeric_limits<Real>::quiet_NaN());
    descriptors_ = std::move(descriptors);
    components_ = std::move(components);
    storage_ = std::move(storage);
}

void TurbulenceFieldSet::clear() noexcept
{
    descriptors_.clear();
    components_.clear();
    storage_ = Field<Real> {};
}

bool TurbulenceFieldSet::contains(const std::string& name) const noexcept
{
    return components_.find(name) != components_.end();
}

int TurbulenceFieldSet::component(const std::string& name) const
{
    const auto iterator = components_.find(name);
    if (iterator == components_.end()) {
        throw std::out_of_range("unknown turbulence field: " + name);
    }
    return iterator->second;
}

const TurbulenceFieldDescriptor&
TurbulenceFieldSet::descriptor(const std::string& name) const
{
    return descriptors_[static_cast<std::size_t>(component(name))];
}

Field<Real>& TurbulenceFieldSet::storage()
{
    if (empty()) throw std::logic_error("empty turbulence field set has no storage");
    return storage_;
}

const Field<Real>& TurbulenceFieldSet::storage() const
{
    if (empty()) throw std::logic_error("empty turbulence field set has no storage");
    return storage_;
}

Real& TurbulenceFieldSet::at(Index3 index, const std::string& name)
{
    return storage()(index.i, index.j, index.k, component(name));
}

const Real& TurbulenceFieldSet::at(Index3 index, const std::string& name) const
{
    return storage()(index.i, index.j, index.k, component(name));
}

void TurbulenceFieldSet::fill(Real value)
{
    if (!std::isfinite(value)) {
        throw std::invalid_argument("turbulence field fill value must be finite");
    }
    storage().fill(value);
}

void TurbulenceFieldSet::validate_interior() const
{
    if (empty()) return;
    const auto extent = storage_.interior_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                for (std::size_t component_index = 0; component_index < descriptors_.size();
                     ++component_index) {
                    const Real value
                        = storage_(i, j, k, static_cast<int>(component_index));
                    const auto& descriptor = descriptors_[component_index];
                    if (!std::isfinite(value)
                        || (descriptor.strictly_positive && !(value > descriptor.lower_bound))
                        || (!descriptor.strictly_positive && value < descriptor.lower_bound)) {
                        throw std::runtime_error("inadmissible turbulence field value: "
                                                 + descriptor.name);
                    }
                }
            }
        }
    }
}

std::vector<Real> TurbulenceFieldSet::pack_interior() const
{
    if (empty()) return {};
    const auto extent = storage_.interior_extent();
    std::vector<Real> result;
    result.reserve(extent.size() * descriptors_.size());
    for (std::size_t component_index = 0; component_index < descriptors_.size();
         ++component_index) {
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    result.push_back(
                        storage_(i, j, k, static_cast<int>(component_index)));
                }
            }
        }
    }
    return result;
}

void TurbulenceFieldSet::unpack_interior(const std::vector<Real>& values)
{
    if (empty()) {
        if (!values.empty()) {
            throw std::invalid_argument("cannot unpack model values into an empty field set");
        }
        return;
    }
    const auto extent = storage_.interior_extent();
    if (values.size() != extent.size() * descriptors_.size()) {
        throw std::invalid_argument("turbulence restart payload has the wrong size");
    }
    std::size_t offset = 0;
    for (std::size_t component_index = 0; component_index < descriptors_.size();
         ++component_index) {
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    storage_(i, j, k, static_cast<int>(component_index)) = values[offset++];
                }
            }
        }
    }
    validate_interior();
}

std::string TurbulenceFieldSet::descriptor_signature() const
{
    std::ostringstream result;
    result << "turbulence_fields_v1";
    for (const auto& descriptor : descriptors_) {
        result << ';' << descriptor.signature();
    }
    return result.str();
}

TurbulenceRestartPayload capture_turbulence_restart(const TurbulenceFieldSet& fields)
{
    fields.validate_interior();
    return {fields.descriptor_signature(), fields.pack_interior()};
}

void restore_turbulence_restart(const TurbulenceRestartPayload& payload,
                                TurbulenceFieldSet& fields)
{
    if (payload.descriptor_signature != fields.descriptor_signature()) {
        throw std::runtime_error("turbulence restart descriptor signature differs");
    }
    fields.unpack_interior(payload.interior_values);
}

} // namespace wcns
