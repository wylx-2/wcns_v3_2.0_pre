#pragma once
#include <wcns/runtime/quantity_registry.hpp>
#include <wcns/runtime/simulation_driver.hpp>
#include <vector>
#include <string>

namespace wcns {
// Time/span statistics on conforming extrusions with j=0 the wall in every
// original zone. Runtime partitioning in any direction is allowed.
class Chapter5Runtime final : public ISimulationObserver {
public:
    Chapter5Runtime(const CaseConfig&,const StatisticContext&);
    void on_initial(const SimulationState&) override;
    void on_step(const SimulationState&,bool) override;
    void on_final(const SimulationState&) override;
    std::string serialize() const;
    void restore(const std::string&);
    std::vector<std::string> output_paths() const;
private:
    struct Zone { BlockId id; int nx,ny,nz; std::size_t offset,wall_offset; };
    struct Probe { std::size_t cell; int k; Real x,y,z; };
    const Zone& zone(BlockId) const;
    void sample();
    void write_history(const SimulationState&);
    void write_means() const;
    const CaseConfig& config_;
    const StatisticContext& context_;
    std::vector<Zone> zones_;
    std::vector<Probe> probes_;
    std::vector<Real> geometry_,walls_,instant_,integral_,wall_,wall_integral_,probe_values_;
    std::array<Real,9> load_integral_{};
    std::array<Real,3> loads_{};
    Real weight_=0,last_time_=0;
    std::size_t samples_=0,last_step_=0,last_history_=static_cast<std::size_t>(-1);
    bool restored_=false;
};
} // namespace wcns
