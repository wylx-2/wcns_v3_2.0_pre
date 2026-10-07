#pragma once
#include <wcns/parallel/slab_fft.hpp>
#include <wcns/runtime/quantity_registry.hpp>
#include <wcns/runtime/simulation_driver.hpp>

namespace wcns {
// Repartition-independent deterministic initial spectrum, collective forcing and
// accepted-state diagnostics. Only O(N^3/P) flow/FFT storage is used per rank.
class HitRuntime final : public ISimulationObserver {
public:
    HitRuntime(const CaseConfig&, const StatisticContext&, LocalBlockSet&, std::string spectrum_path);
    void initialize();
    bool preparing() const { return !preparation_complete_; }
    void finish_preparation(const SimulationState&);
    void add_stage_source(Real time);
    void accepted_step_transform(Real dt);
    void on_initial(const SimulationState&) override;
    void on_step(const SimulationState&,bool) override;
    void on_final(const SimulationState&) override;
    Real next_time_event(const SimulationState&) const override;
    std::string serialize() const;
    void restore(const std::string&);
    std::vector<std::string> output_paths() const;
private:
    struct Cell { StructuredBlock* block; Index3 index; std::size_t global; };
    std::vector<Real> to_slab(const std::vector<Real>&,int components) const;
    std::vector<Real> from_slab(const std::vector<Real>&,int components) const;
    std::vector<Real> velocities(bool residual=false) const;
    std::array<FourierField,3> velocity_spectrum() const;
    void install_velocity(std::array<FourierField,3>, bool fresh);
    void sample(const SimulationState&,bool write);
    void write_means() const;
    void root_io(const std::function<void()>&) const;
    const CaseConfig& config_;
    const StatisticContext& context_;
    LocalBlockSet& blocks_;
    SlabFft fft_;
    std::vector<Cell> cells_;
    std::vector<std::size_t> send_order_, received_indices_;
    std::vector<int> send_counts_, receive_counts_;
    std::vector<std::array<Real,2>> target_spectrum_;
    std::string spectrum_identity_;
    std::vector<Real> instant_, integral_, previous_, spectrum_integral_;
    Real weight_=0, last_sample_time_=0, last_step_time_=0;
    Real force_power_=0, cooling_=0, force_alpha_=0;
    std::size_t samples_=0,last_sample_step_=0,last_step_=0;
    bool restored_=false, sampled_=false;
    bool preparation_complete_=true;
    std::size_t preparation_steps_=0;
};
} // namespace wcns
