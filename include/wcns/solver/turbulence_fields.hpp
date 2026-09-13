#pragma once

#include <wcns/core/field.hpp>
#include <wcns/core/types.hpp>

#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace wcns {

enum class TurbulenceFieldRole {
    Transported,
    Auxiliary,
    Diagnostic,
};

enum class TurbulenceFieldScale {
    Dimensionless,
    VelocitySquared,
    KinematicViscosity,
    DynamicViscosity,
    InverseTime,
    Dissipation,
    Length,
};

struct TurbulenceFieldDescriptor {
    std::string name;
    TurbulenceFieldRole role = TurbulenceFieldRole::Transported;
    TurbulenceFieldScale scale = TurbulenceFieldScale::Dimensionless;
    bool strictly_positive = false;
    Real lower_bound = 0.0;

    void validate() const;
    [[nodiscard]] std::string signature() const;
};

// Model scalars remain separate from the five-component mean-flow layout.
// Components are cell-major and share one allocation so the existing generic
// BlockFieldRegistry can exchange any registered model in one halo message.
class TurbulenceFieldSet {
public:
    TurbulenceFieldSet() = default;
    TurbulenceFieldSet(Extent3 extent,
                       int ghost_width,
                       std::vector<TurbulenceFieldDescriptor> descriptors);

    void reset(Extent3 extent,
               int ghost_width,
               std::vector<TurbulenceFieldDescriptor> descriptors);
    void clear() noexcept;

    [[nodiscard]] bool empty() const noexcept { return descriptors_.empty(); }
    [[nodiscard]] int components() const noexcept
    {
        return static_cast<int>(descriptors_.size());
    }
    [[nodiscard]] const std::vector<TurbulenceFieldDescriptor>& descriptors() const noexcept
    {
        return descriptors_;
    }
    [[nodiscard]] bool contains(const std::string& name) const noexcept;
    [[nodiscard]] int component(const std::string& name) const;
    [[nodiscard]] const TurbulenceFieldDescriptor& descriptor(const std::string& name) const;

    [[nodiscard]] Field<Real>& storage();
    [[nodiscard]] const Field<Real>& storage() const;
    [[nodiscard]] Real& at(Index3 index, const std::string& name);
    [[nodiscard]] const Real& at(Index3 index, const std::string& name) const;

    void fill(Real value);
    void validate_interior() const;
    [[nodiscard]] std::vector<Real> pack_interior() const;
    void unpack_interior(const std::vector<Real>& values);
    [[nodiscard]] std::string descriptor_signature() const;

private:
    std::vector<TurbulenceFieldDescriptor> descriptors_;
    std::unordered_map<std::string, int> components_;
    Field<Real> storage_;
};

struct TurbulenceRestartPayload {
    std::string descriptor_signature;
    std::vector<Real> interior_values;
};

[[nodiscard]] TurbulenceRestartPayload
capture_turbulence_restart(const TurbulenceFieldSet& fields);
void restore_turbulence_restart(const TurbulenceRestartPayload& payload,
                                TurbulenceFieldSet& fields);

} // namespace wcns
