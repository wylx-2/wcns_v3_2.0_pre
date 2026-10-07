#pragma once

#include <wcns/runtime/quantity_registry.hpp>
#include <wcns/runtime/simulation_driver.hpp>
#include <functional>
#include <string>
#include <vector>

namespace wcns {

// Single extruded structured source zone. Original i,j indices survive MPI
// repartitioning; all persisted moments are globally reduced and rank independent.
class PeriodicHillRuntime final : public ISimulationObserver {
public:
    PeriodicHillRuntime(const CaseConfig& config, const StatisticContext& context);
    void bind_force(std::function<void(Real)> setter);
    void on_initial(const SimulationState& state) override;
    void on_step(const SimulationState& state, bool) override;
    void on_final(const SimulationState& state) override;
    [[nodiscard]] std::string serialize() const;
    void restore(const std::string& text);
    [[nodiscard]] std::vector<std::string> output_paths() const;
    [[nodiscard]] Real force() const noexcept { return force_; }
    [[nodiscard]] Real weight() const noexcept { return weight_; }

private:
    void sample_flow();
    void write_history(const SimulationState& state);
    void write_averages() const;
    const CaseConfig& config_;
    const StatisticContext& context_;
    int nx_ = 0, ny_ = 0;
    // Geometry: x,y,span weight, y lower/upper wall, lower-wall slope,
    // crest strip area per unit span. Moments: rho,u,v,w,p,T,uu,vv,ww,uv,uw,vw,rho*u.
    std::vector<Real> geometry_, instant_, integral_, wall_, wall_integral_;
    Real force_ = 0, previous_error_ = 0, weight_ = 0;
    Real mass_flux_ = 0, volume_flux_ = 0;
    std::size_t samples_ = 0, last_step_ = 0;
    std::size_t last_history_step_ = static_cast<std::size_t>(-1);
    Real last_time_ = 0;
    bool restored_ = false;
    std::function<void(Real)> setter_;
};
} // namespace wcns
