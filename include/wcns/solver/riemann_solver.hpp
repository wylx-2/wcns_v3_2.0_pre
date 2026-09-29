#pragma once

#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/euler.hpp>
#include <wcns/solver/low_mach_preconditioner.hpp>
#include <wcns/solver/numerical_diagnostics.hpp>

#include <array>
#include <functional>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace wcns {

enum class RiemannSolverKind {
    Rusanov,
    Hll,
    Hllc,
    Roe,
    AllSpeedRoe,
    RieperAllSpeedRoe,
    RotatedRoe,
};

enum class RiemannFallbackReason {
    None,
    InvalidWaveSpeed,
    InvalidIntermediateState,
    InvalidRoeAverage,
    NonFiniteFlux,
};

struct LiGuAllSpeedRoeParameters {
    Real reference_mach = 0.1;
    Real dissipation_scale = 0.02;
    Real pressure_coefficient = 0.05;

    void validate() const;
};

struct RiemannSolverParameters {
    Real entropy_fix_coefficient = 0.1;
    Real denominator_tolerance = 1.0e-12;
    bool weiss_smith = false;
    WeissSmithParameters preconditioner {};
    LiGuAllSpeedRoeParameters all_speed {};

    void validate() const;
};

struct RiemannFaceContext {
    Real viscous_speed = 0.0;
};

struct RiemannConfig {
    std::string scheme = "rusanov";
    RiemannSolverParameters parameters {};

    void validate() const;
    void validate(const class RiemannSolverRegistry& registry) const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

struct RiemannFallbackStep {
    std::string from_solver;
    std::string to_solver;
    RiemannFallbackReason reason = RiemannFallbackReason::None;
};

struct RiemannResult {
    ConservativeState flux_per_unit_area {};
    Real spectral_radius = 0.0;
    std::string requested_solver;
    std::string used_solver;
    RiemannFallbackReason fallback_reason = RiemannFallbackReason::None;
    std::vector<RiemannFallbackStep> fallback_path;
};

struct RiemannFallbackEvent {
    FaceDiagnosticLocation location {};
    std::string from_solver;
    std::string to_solver;
    RiemannFallbackReason reason = RiemannFallbackReason::None;
};

struct RiemannDiagnostics {
    std::size_t total_faces = 0;
    std::unordered_map<std::string, std::size_t> requested_faces;
    std::unordered_map<std::string, std::size_t> used_faces;
    std::array<std::size_t, 5> fallback_reasons {};
    std::vector<RiemannFallbackEvent> fallback_events;

    void record(const RiemannResult& result, FaceDiagnosticLocation location);
    [[nodiscard]] std::size_t fallback_count() const noexcept;
    [[nodiscard]] std::size_t fallback_count(RiemannFallbackReason reason) const;
};

[[nodiscard]] std::string_view riemann_fallback_reason_name(RiemannFallbackReason reason);

class IRiemannSolver {
public:
    virtual ~IRiemannSolver() = default;
    [[nodiscard]] virtual std::string_view name() const noexcept = 0;
    [[nodiscard]] virtual RiemannResult solve(const PressurePrimitiveState& left,
                                              const PressurePrimitiveState& right,
                                              Normal3 unit_normal,
                                              const GasModel& gas,
                                              const NumericalFloors& floors) const = 0;
};

class RiemannSolverRegistry {
public:
    using Factory = std::function<std::unique_ptr<IRiemannSolver>()>;

    void register_solver(std::string name, Factory factory);
    [[nodiscard]] bool contains(std::string_view name) const noexcept;
    [[nodiscard]] std::unique_ptr<IRiemannSolver> create(std::string_view name) const;
    [[nodiscard]] std::vector<std::string> names() const;

    [[nodiscard]] static RiemannSolverRegistry
    with_builtins(const RiemannSolverParameters& parameters = {});

private:
    std::unordered_map<std::string, Factory> factories_;
};

[[nodiscard]] std::string_view riemann_solver_name(RiemannSolverKind kind);

class RiemannSolver {
public:
    explicit RiemannSolver(RiemannSolverKind kind = RiemannSolverKind::Rusanov,
                           RiemannSolverParameters parameters = {});
    explicit RiemannSolver(std::string_view name,
                           const RiemannSolverRegistry& registry,
                           RiemannSolverParameters parameters = {});

    [[nodiscard]] std::string_view name() const noexcept;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;

    [[nodiscard]] RiemannResult solve(const PressurePrimitiveState& left,
                                      const PressurePrimitiveState& right,
                                      Normal3 unit_normal,
                                      const GasModel& gas,
                                      const NumericalFloors& floors,
                                      RiemannFaceContext context = {}) const;

    [[nodiscard]] ConservativeState flux(const PressurePrimitiveState& left,
                                         const PressurePrimitiveState& right,
                                         Normal3 unit_normal,
                                         const GasModel& gas,
                                         const NumericalFloors& floors,
                                         RiemannFaceContext context = {}) const;

private:
    std::unique_ptr<IRiemannSolver> implementation_;
    RiemannSolverParameters parameters_ {};
};

} // namespace wcns
