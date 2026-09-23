#include <wcns/runtime/weighted_statistics.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace wcns {
namespace {

bool same_time(Real lhs, Real rhs)
{
    const Real tolerance = 64.0 * std::numeric_limits<Real>::epsilon()
        * std::max({Real {1.0}, std::abs(lhs), std::abs(rhs)});
    return std::abs(lhs - rhs) <= tolerance;
}

} // namespace

void WeightedMomentState::validate() const
{
    if (!std::isfinite(weight) || weight < 0.0 || !std::isfinite(mean)
        || !std::isfinite(second_central) || second_central < -1.0e-13
        || (sample_count == 0) != (weight == 0.0)) {
        throw std::invalid_argument("weighted moment state is invalid");
    }
}

void WeightedMomentState::add(Real value, Real sample_weight)
{
    validate();
    if (!std::isfinite(value) || !std::isfinite(sample_weight) || sample_weight <= 0.0) {
        throw std::invalid_argument("weighted moment sample is invalid");
    }
    const Real next_weight = weight + sample_weight;
    const Real delta = value - mean;
    const Real next_mean = mean + sample_weight * delta / next_weight;
    second_central += sample_weight * delta * (value - next_mean);
    weight = next_weight;
    mean = next_mean;
    ++sample_count;
    if (second_central < 0.0 && second_central > -1.0e-13) second_central = 0.0;
    validate();
}

void WeightedMomentState::merge(const WeightedMomentState& other)
{
    validate();
    other.validate();
    if (other.sample_count == 0) return;
    if (sample_count == 0) {
        *this = other;
        return;
    }
    const Real combined_weight = weight + other.weight;
    const Real delta = other.mean - mean;
    second_central += other.second_central
        + delta * delta * weight * other.weight / combined_weight;
    mean += delta * other.weight / combined_weight;
    weight = combined_weight;
    sample_count += other.sample_count;
    validate();
}

Real WeightedMomentState::variance() const
{
    validate();
    return weight == 0.0 ? 0.0 : std::max(Real {0.0}, second_central / weight);
}

Real WeightedMomentState::rms() const
{
    return std::sqrt(variance());
}

void WeightedCovarianceState::validate() const
{
    if (!std::isfinite(weight) || weight < 0.0 || !std::isfinite(mean_x)
        || !std::isfinite(mean_y) || !std::isfinite(co_moment)
        || (sample_count == 0) != (weight == 0.0)) {
        throw std::invalid_argument("weighted covariance state is invalid");
    }
}

void WeightedCovarianceState::add(Real x, Real y, Real sample_weight)
{
    validate();
    if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(sample_weight)
        || sample_weight <= 0.0) {
        throw std::invalid_argument("weighted covariance sample is invalid");
    }
    const Real next_weight = weight + sample_weight;
    const Real dx = x - mean_x;
    const Real dy = y - mean_y;
    const Real next_mean_x = mean_x + sample_weight * dx / next_weight;
    const Real next_mean_y = mean_y + sample_weight * dy / next_weight;
    co_moment += sample_weight * dx * (y - next_mean_y);
    weight = next_weight;
    mean_x = next_mean_x;
    mean_y = next_mean_y;
    ++sample_count;
    validate();
}

void WeightedCovarianceState::merge(const WeightedCovarianceState& other)
{
    validate();
    other.validate();
    if (other.sample_count == 0) return;
    if (sample_count == 0) {
        *this = other;
        return;
    }
    const Real combined_weight = weight + other.weight;
    const Real dx = other.mean_x - mean_x;
    const Real dy = other.mean_y - mean_y;
    co_moment += other.co_moment + dx * dy * weight * other.weight / combined_weight;
    mean_x += dx * other.weight / combined_weight;
    mean_y += dy * other.weight / combined_weight;
    weight = combined_weight;
    sample_count += other.sample_count;
    validate();
}

Real WeightedCovarianceState::covariance() const
{
    validate();
    return weight == 0.0 ? 0.0 : co_moment / weight;
}

AcceptedTimeStatistics::AcceptedTimeStatistics(std::string identity,
                                               std::size_t variable_count)
{
    if (identity.empty() || variable_count == 0) {
        throw std::invalid_argument("accepted statistics identity or size is invalid");
    }
    state_.identity = std::move(identity);
    state_.reynolds.resize(variable_count);
    state_.favre.resize(variable_count);
}

AcceptedTimeStatistics::AcceptedTimeStatistics(AcceptedStatisticsState state)
    : state_(std::move(state))
{
    if (state_.identity.empty() || state_.reynolds.empty()
        || state_.reynolds.size() != state_.favre.size()) {
        throw std::invalid_argument("restored accepted statistics identity is invalid");
    }
    for (const auto& entry : state_.reynolds) entry.validate();
    for (const auto& entry : state_.favre) entry.validate();
    if (state_.accepted_events != state_.reynolds.front().sample_count
        || state_.accepted_events != state_.favre.front().sample_count) {
        throw std::invalid_argument("restored accepted statistics sample counts differ");
    }
}

bool AcceptedTimeStatistics::sample(std::size_t step,
                                    Real time,
                                    Real physical_time_step,
                                    const std::vector<Real>& values,
                                    Real density,
                                    bool accepted_physical_step)
{
    if (!accepted_physical_step) return false;
    if (!std::isfinite(time) || time < 0.0 || !std::isfinite(physical_time_step)
        || physical_time_step <= 0.0 || !std::isfinite(density) || density <= 0.0
        || values.size() != state_.reynolds.size()
        || !std::all_of(values.begin(), values.end(), [](Real value) {
               return std::isfinite(value);
           })) {
        throw std::invalid_argument("accepted statistics event is invalid");
    }
    if (state_.accepted_events > 0) {
        if (step == state_.last_step && same_time(time, state_.last_time)) return false;
        if (step <= state_.last_step || time <= state_.last_time) {
            throw std::invalid_argument("accepted statistics events are not strictly ordered");
        }
    }
    for (std::size_t variable = 0; variable < values.size(); ++variable) {
        state_.reynolds[variable].add(values[variable], physical_time_step);
        state_.favre[variable].add(values[variable], density * physical_time_step);
    }
    if (state_.accepted_events == 0) {
        state_.first_step = step;
        state_.first_time = time;
    }
    ++state_.accepted_events;
    state_.last_step = step;
    state_.last_time = time;
    return true;
}

std::string AcceptedTimeStatistics::serialize() const
{
    std::ostringstream output;
    output << "wcns_weighted_statistics_v1\n" << std::quoted(state_.identity) << '\n'
           << state_.accepted_events << ' ' << state_.first_step << ' ' << state_.last_step << ' '
           << std::setprecision(17) << state_.first_time << ' ' << state_.last_time << ' '
           << state_.reynolds.size() << '\n';
    for (std::size_t variable = 0; variable < state_.reynolds.size(); ++variable) {
        const auto& r = state_.reynolds[variable];
        const auto& f = state_.favre[variable];
        output << r.sample_count << ' ' << r.weight << ' ' << r.mean << ' '
               << r.second_central << ' ' << f.sample_count << ' ' << f.weight << ' '
               << f.mean << ' ' << f.second_central << '\n';
    }
    return output.str();
}

AcceptedTimeStatistics AcceptedTimeStatistics::deserialize(const std::string& text)
{
    std::istringstream input(text);
    std::string version;
    std::getline(input, version);
    if (version != "wcns_weighted_statistics_v1") {
        throw std::invalid_argument("unsupported weighted statistics state version");
    }
    AcceptedStatisticsState state;
    std::size_t variables = 0;
    if (!(input >> std::quoted(state.identity) >> state.accepted_events >> state.first_step
          >> state.last_step >> state.first_time >> state.last_time >> variables)
        || variables == 0) {
        throw std::invalid_argument("weighted statistics state header is invalid");
    }
    state.reynolds.resize(variables);
    state.favre.resize(variables);
    for (std::size_t variable = 0; variable < variables; ++variable) {
        auto& r = state.reynolds[variable];
        auto& f = state.favre[variable];
        if (!(input >> r.sample_count >> r.weight >> r.mean >> r.second_central
              >> f.sample_count >> f.weight >> f.mean >> f.second_central)) {
            throw std::invalid_argument("weighted statistics state is truncated");
        }
    }
    std::string trailing;
    if (input >> trailing) {
        throw std::invalid_argument("weighted statistics state has trailing data");
    }
    return AcceptedTimeStatistics(std::move(state));
}

} // namespace wcns
