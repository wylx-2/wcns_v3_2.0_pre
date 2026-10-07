#pragma once

#include <wcns/core/types.hpp>

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace wcns {

enum class SourceModelKind {
    UniformConservative,
    BodyForce,
    PressureGradient,
    ManufacturedSolution,
    RampTrip,
};

[[nodiscard]] const char* source_model_name(SourceModelKind kind);

struct SourceTermConfig {
    bool enable_source_terms = false;
    std::vector<SourceModelKind> models;
    std::array<Real, 5> uniform_conservative {{0.0, 0.0, 0.0, 0.0, 0.0}};
    std::array<Real, 3> body_acceleration {{0.0, 0.0, 0.0}};
    // Constant force per unit volume G=-grad(p); energy source is u dot G.
    std::array<Real, 3> pressure_gradient {{0.0, 0.0, 0.0}};
    std::array<Real, 5> manufactured_amplitude {{0.0, 0.0, 0.0, 0.0, 0.0}};
    Real ramp_trip_amplitude = 0.0;
    Real ramp_trip_span = 6.0;

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

// Stage H deliberately owns no concrete source model. The factory validates
// the complete configuration and rejects every enabled configuration until
// actual source evaluation is implemented in stage J.
class SourceTermRegistry {
public:
    [[nodiscard]] static SourceTermRegistry create_stage_h(const SourceTermConfig& config);
    [[nodiscard]] static SourceTermRegistry create_stage_j(const SourceTermConfig& config);

    [[nodiscard]] std::array<Real, 5> evaluate(const std::array<Real, 5>& conservative,
                                               const std::array<Real, 3>& coordinates,
                                               Real time,
                                               int dimension) const;

    void set_pressure_gradient_x(Real force);

    [[nodiscard]] constexpr bool empty() const noexcept { return model_count_ == 0; }
    [[nodiscard]] constexpr std::size_t size() const noexcept { return model_count_; }

private:
    std::size_t model_count_ = 0;
    SourceTermConfig config_ {};
};

} // namespace wcns
