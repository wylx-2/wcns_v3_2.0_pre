#pragma once

#include <wcns/core/types.hpp>

#include <cstddef>
#include <string>
#include <vector>

namespace wcns {

struct WeightedMomentState {
    std::size_t sample_count = 0;
    Real weight = 0.0;
    Real mean = 0.0;
    Real second_central = 0.0;

    void add(Real value, Real sample_weight);
    void merge(const WeightedMomentState& other);
    [[nodiscard]] Real variance() const;
    [[nodiscard]] Real rms() const;
    void validate() const;
};

struct WeightedCovarianceState {
    std::size_t sample_count = 0;
    Real weight = 0.0;
    Real mean_x = 0.0;
    Real mean_y = 0.0;
    Real co_moment = 0.0;

    void add(Real x, Real y, Real sample_weight);
    void merge(const WeightedCovarianceState& other);
    [[nodiscard]] Real covariance() const;
    void validate() const;
};

struct AcceptedStatisticsState {
    std::string identity;
    std::size_t accepted_events = 0;
    std::size_t first_step = 0;
    std::size_t last_step = 0;
    Real first_time = 0.0;
    Real last_time = 0.0;
    std::vector<WeightedMomentState> reynolds;
    std::vector<WeightedMomentState> favre;
    // Unique quantity pairs in lexicographic (first, second) order with
    // first < second.  This keeps checkpoint and text-output ordering stable.
    std::vector<WeightedCovarianceState> reynolds_covariances;
    std::vector<WeightedCovarianceState> favre_covariances;
};

class AcceptedTimeStatistics {
public:
    AcceptedTimeStatistics(std::string identity, std::size_t variable_count);

    // Returns false for an exact replay of the last accepted event. Pseudo-time
    // iterations and rejected attempts must call with accepted_physical_step=false.
    bool sample(std::size_t step,
                Real time,
                Real physical_time_step,
                const std::vector<Real>& values,
                Real density,
                bool accepted_physical_step);

    [[nodiscard]] const AcceptedStatisticsState& state() const noexcept { return state_; }
    [[nodiscard]] std::string serialize() const;
    [[nodiscard]] static AcceptedTimeStatistics deserialize(const std::string& text);

private:
    explicit AcceptedTimeStatistics(AcceptedStatisticsState state);
    AcceptedStatisticsState state_;
};

} // namespace wcns
