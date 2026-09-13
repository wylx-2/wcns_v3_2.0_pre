#include <wcns/runtime/output_manager.hpp>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <utility>

#ifndef WCNS_GIT_COMMIT
#define WCNS_GIT_COMMIT "unknown"
#endif
#ifndef WCNS_PROGRAM_VERSION
#define WCNS_PROGRAM_VERSION "unknown"
#endif

#if defined(_WIN32)
#include <direct.h>
#define WCNS_MKDIR(path) _mkdir(path)
#else
#include <sys/types.h>
#define WCNS_MKDIR(path) mkdir(path, 0777)
#endif

namespace wcns {
namespace {

const std::array<const char*, euler_components> residual_names {{
    "rho",
    "rhou",
    "rhov",
    "rhow",
    "rhoE",
}};

Real time_tolerance(Real lhs, Real rhs)
{
    return 64.0 * std::numeric_limits<Real>::epsilon()
        * std::max({Real {1.0}, std::abs(lhs), std::abs(rhs)});
}

std::string join_path(const std::string& directory, const std::string& name)
{
    if (directory.empty()) return name;
    const char last = directory.back();
    return directory + ((last == '/' || last == '\\') ? "" : "/") + name;
}

bool path_exists(const std::string& path)
{
    struct stat information {};
    return stat(path.c_str(), &information) == 0;
}

bool path_is_directory(const std::string& path)
{
    struct stat information {};
    if (stat(path.c_str(), &information) != 0) return false;
#if defined(_WIN32)
    return (information.st_mode & _S_IFDIR) != 0;
#else
    return S_ISDIR(information.st_mode);
#endif
}

void make_directory_tree(const std::string& input)
{
    if (input.empty()) throw std::runtime_error("output directory is empty");
    std::string path = input;
    std::replace(path.begin(), path.end(), '\\', '/');
    std::size_t begin = 0;
    if (path.size() >= 3 && path[1] == ':' && path[2] == '/')
        begin = 3;
    else if (path.front() == '/')
        begin = 1;
    for (std::size_t index = begin; index <= path.size(); ++index) {
        if (index != path.size() && path[index] != '/') continue;
        const auto prefix = path.substr(0, index);
        if (prefix.empty() || (prefix.size() == 2 && prefix[1] == ':')) continue;
        if (path_exists(prefix)) {
            if (!path_is_directory(prefix)) {
                throw std::runtime_error("output path component is not a directory: " + prefix);
            }
            continue;
        }
        errno = 0;
        if (WCNS_MKDIR(prefix.c_str()) != 0 && errno != EEXIST) {
            throw std::runtime_error("cannot create output directory: " + prefix);
        }
    }
}

std::string safe_name(std::string name)
{
    for (char& character : name) {
        const bool safe = (character >= 'a' && character <= 'z')
            || (character >= 'A' && character <= 'Z') || (character >= '0' && character <= '9')
            || character == '-' || character == '_';
        if (!safe) character = '_';
    }
    return name.empty() ? "case" : name;
}

void atomic_replace(const std::string& temporary, const std::string& target, bool allow_existing)
{
    if (path_exists(target)) {
        if (!allow_existing) {
            throw std::runtime_error("output file already exists: " + target);
        }
        if (std::remove(target.c_str()) != 0) {
            throw std::runtime_error("cannot replace output file: " + target);
        }
    }
    if (std::rename(temporary.c_str(), target.c_str()) != 0) {
        throw std::runtime_error("cannot commit output file: " + target);
    }
}

} // namespace

OutputSchedule::OutputSchedule(OutputScheduleConfig config)
    : config_(std::move(config))
{
    config_.validate();
}

Real OutputSchedule::next_time(Real current_time) const
{
    Real result = std::numeric_limits<Real>::infinity();
    if (config_.every_time > 0.0) {
        const Real count = std::floor((current_time + time_tolerance(current_time, current_time))
                                      / config_.every_time)
            + 1.0;
        result = count * config_.every_time;
    }
    for (const Real time : config_.explicit_times) {
        if (time > current_time + time_tolerance(time, current_time)) {
            result = std::min(result, time);
            break;
        }
    }
    return result;
}

bool OutputSchedule::scheduled(const SimulationState& state) const
{
    if (config_.every_steps > 0 && state.step > 0 && state.step % config_.every_steps == 0) {
        return true;
    }
    if (config_.every_time > 0.0 && state.time > 0.0) {
        const Real multiple = std::round(state.time / config_.every_time) * config_.every_time;
        if (std::abs(state.time - multiple) <= time_tolerance(state.time, multiple)) {
            return true;
        }
    }
    for (const Real time : config_.explicit_times) {
        if (std::abs(state.time - time) <= time_tolerance(state.time, time)) {
            return true;
        }
    }
    return false;
}

bool OutputSchedule::already_emitted(const SimulationState& state) const
{
    return has_emitted_ && state.step == last_step_
        && std::abs(state.time - last_time_) <= time_tolerance(state.time, last_time_);
}

bool OutputSchedule::consume(const SimulationState& state, bool initial, bool final)
{
    const bool due
        = (initial && config_.write_initial) || (final && config_.write_final) || scheduled(state);
    if (!due || already_emitted(state)) return false;
    has_emitted_ = true;
    last_step_ = state.step;
    last_time_ = state.time;
    return true;
}

RuntimeOutputManager::RuntimeOutputManager(const MpiRuntime& mpi,
                                           const CaseConfig& config,
                                           const StructuredPartitionPlan& partition,
                                           std::string mesh_signature,
                                           const StatisticContext* statistic_context,
                                           EventWriter event_writer,
                                           StatisticRegistry statistic_registry)
    : mpi_(mpi)
    , config_(config)
    , partition_(partition)
    , mesh_signature_(std::move(mesh_signature))
    , statistic_context_(statistic_context)
    , statistic_registry_(std::move(statistic_registry))
    , event_writer_(std::move(event_writer))
    , field_schedule_(config.output.field.schedule)
    , history_schedule_(config.output.history.schedule)
    , statistics_schedule_(config.output.statistics.schedule)
    , boundary_schedule_(config.output.boundary.schedule)
    , checkpoint_schedule_(config.output.checkpoint.schedule)
    , output_directory_(config.output.directory)
{ }

RuntimeOutputManager::~RuntimeOutputManager()
{
    if (history_stream_.is_open()) history_stream_.close();
    if (statistics_stream_.is_open()) statistics_stream_.close();
}

void RuntimeOutputManager::prepare_directory()
{
    if (prepared_) return;
    std::string status;
    if (mpi_.rank() == 0) {
        try {
            const bool exists = path_exists(output_directory_);
            if (exists && !path_is_directory(output_directory_)) {
                throw std::runtime_error("output root exists and is not a directory: "
                                         + output_directory_);
            }
            if (exists && !config_.output.allow_existing) {
                throw std::runtime_error("output directory already exists: " + output_directory_);
            }
            if (!exists) make_directory_tree(output_directory_);
            status = "OK\n";
        } catch (const std::exception& error) {
            status = "ERROR\n" + std::string(error.what());
        }
    }
    status = mpi_.broadcast_string(std::move(status));
    if (status.rfind("ERROR\n", 0) == 0) {
        throw std::runtime_error(status.substr(6));
    }
    if (status != "OK\n") {
        throw std::runtime_error("invalid output directory status broadcast");
    }
    prepared_ = true;

    if (config_.output.history.enabled && mpi_.rank() == 0) {
        const auto base = safe_name(config_.case_name) + ".history.r" + std::to_string(mpi_.size());
        const char* extension
            = config_.output.history.format == SeriesOutputFormat::Text ? ".txt" : ".dat";
        history_final_path_ = join_path(output_directory_, base + extension);
        history_temporary_path_ = history_final_path_ + ".tmp";
        history_stream_.open(history_temporary_path_, std::ios::out | std::ios::trunc);
        if (!history_stream_) {
            throw std::runtime_error("cannot open residual history temporary file: "
                                     + history_temporary_path_);
        }
        if (config_.output.history.format == SeriesOutputFormat::Tecplot) {
            history_stream_ << "TITLE=\"WCNS residual history\"\n"
                            << "VARIABLES=\"step\",\"time\",\"dt\",\"cfl\","
                               "\"wall_time\",\"total_l2\"";
            for (const char* suffix : {"l2",
                                       "reference_l2",
                                       "normalized_l2",
                                       "linf",
                                       "reference_linf",
                                       "normalized_linf"}) {
                for (const char* component : residual_names) {
                    history_stream_ << ",\"" << component << '_' << suffix << "\"";
                }
                if (config_.turbulence.kind == TurbulenceModelKind::SaNegative) {
                    history_stream_ << ",\"nu_tilde_" << suffix << "\"";
                }
            }
            history_stream_ << ",\"consecutive\",\"reconstruction_fallbacks\","
                               "\"riemann_fallbacks\",\"robustness_level0_faces\","
                               "\"robustness_level1_faces\",\"robustness_level2_faces\","
                               "\"robustness_level3_faces\",\"troubled_cells\","
                               "\"local_recomputations\",\"step_retries\","
                               "\"proposed_dt\",\"accepted_dt\",\"minimum_rho\","
                               "\"minimum_p\",\"minimum_T\",\"minimum_e\","
                               "\"residual_checked\","
                               "\"stop_reason_code\"\n"
                               "AUXDATA STOP_REASON_CODES=\"0=running;1=steady_converged;"
                               "2=physical_time_reached;3=maximum_steps;"
                               "4=wall_time_checkpoint;5=user_signal_checkpoint;"
                               "6=numerical_failure\"\n"
                               "ZONE T=\"history\"\n";
        } else {
            history_stream_ << "# step time dt cfl wall_time total_l2";
            for (const char* suffix : {"l2",
                                       "reference_l2",
                                       "normalized_l2",
                                       "linf",
                                       "reference_linf",
                                       "normalized_linf"}) {
                for (const char* component : residual_names) {
                    history_stream_ << ' ' << component << '_' << suffix;
                }
                if (config_.turbulence.kind == TurbulenceModelKind::SaNegative) {
                    history_stream_ << " nu_tilde_" << suffix;
                }
            }
            history_stream_ << " consecutive reconstruction_fallbacks riemann_fallbacks "
                               "robustness_level0_faces robustness_level1_faces "
                               "robustness_level2_faces robustness_level3_faces "
                               "troubled_cells local_recomputations step_retries "
                               "proposed_dt accepted_dt minimum_rho minimum_p minimum_T minimum_e "
                               "residual_checked stop_reason\n";
        }
    }
    if (config_.output.statistics.enabled) {
        if (statistic_context_ == nullptr) {
            throw std::runtime_error("statistics output is enabled without a statistic context");
        }
        statistic_registry_.validate_selection(config_.output.statistics.quantities);
        if (mpi_.rank() == 0) {
            const auto base
                = safe_name(config_.case_name) + ".statistics.r" + std::to_string(mpi_.size());
            const char* extension
                = config_.output.statistics.format == SeriesOutputFormat::Text ? ".txt" : ".dat";
            statistics_final_path_ = join_path(output_directory_, base + extension);
            statistics_temporary_path_ = statistics_final_path_ + ".tmp";
            statistics_stream_.open(statistics_temporary_path_, std::ios::out | std::ios::trunc);
            if (!statistics_stream_) {
                throw std::runtime_error("cannot open statistics temporary file: "
                                         + statistics_temporary_path_);
            }
            if (config_.output.statistics.format == SeriesOutputFormat::Tecplot) {
                statistics_stream_ << "TITLE=\"WCNS statistics\"\n"
                                   << "VARIABLES=\"step\",\"time\"";
                for (const auto& name : config_.output.statistics.quantities) {
                    statistics_stream_ << ",\"" << name << "\"";
                }
                statistics_stream_ << "\nZONE T=\"statistics\"\n";
            } else {
                statistics_stream_ << "# step time";
                for (const auto& name : config_.output.statistics.quantities) {
                    statistics_stream_ << ' ' << name;
                }
                statistics_stream_ << '\n';
            }
        }
    }
}

Real RuntimeOutputManager::next_time_event(const SimulationState& state) const
{
    Real result = std::numeric_limits<Real>::infinity();
    if (config_.output.field.enabled) {
        result = std::min(result, field_schedule_.next_time(state.time));
    }
    if (config_.output.history.enabled) {
        result = std::min(result, history_schedule_.next_time(state.time));
    }
    if (config_.output.statistics.enabled) {
        result = std::min(result, statistics_schedule_.next_time(state.time));
    }
    if (config_.output.boundary.enabled) {
        result = std::min(result, boundary_schedule_.next_time(state.time));
    }
    if (config_.output.checkpoint.enabled) {
        result = std::min(result, checkpoint_schedule_.next_time(state.time));
    }
    return result;
}

void RuntimeOutputManager::dispatch(OutputCategory category,
                                    OutputSchedule& schedule,
                                    const SimulationState& state,
                                    bool initial,
                                    bool final,
                                    bool enabled)
{
    if (!enabled || !schedule.consume(state, initial, final)) return;
    if (!event_writer_) {
        throw std::runtime_error("configured output category has no production writer");
    }
    const auto paths = event_writer_(category, state, initial, final);
    for (const auto& path : paths)
        record_file(path);
}

void RuntimeOutputManager::write_history(const SimulationState& state, bool residual_checked)
{
    if (mpi_.rank() != 0) return;
    if (!history_stream_) {
        throw std::runtime_error("residual history stream is not writable");
    }
    const std::size_t expected_model_components
        = config_.turbulence.kind == TurbulenceModelKind::SaNegative ? 1u : 0u;
    if (state.residuals.model_l2.size() != expected_model_components
        || state.residuals.model_linf.size() != expected_model_components
        || (state.steady.reference_initialized
            && (state.steady.model_reference_l2.size() != expected_model_components
                || state.steady.model_reference_linf.size() != expected_model_components))) {
        throw std::runtime_error("history model residual identity differs from configuration");
    }
    history_stream_ << std::setprecision(17) << state.step << ' ' << state.time << ' '
                    << state.time_step << ' ' << config_.run.cfl << ' ' << state.wall_time << ' '
                    << state.residuals.total_l2();
    for (const Real value : state.residuals.l2)
        history_stream_ << ' ' << value;
    for (const Real value : state.residuals.model_l2)
        history_stream_ << ' ' << value;
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    for (const Real value : state.steady.reference_l2) {
        history_stream_ << ' ' << (state.steady.reference_initialized ? value : nan);
    }
    for (std::size_t component = 0; component < expected_model_components; ++component) {
        history_stream_ << ' '
                        << (state.steady.reference_initialized
                                ? state.steady.model_reference_l2[component]
                                : nan);
    }
    for (std::size_t component = 0; component < euler_components; ++component) {
        history_stream_ << ' '
                        << (state.steady.reference_initialized ? state.residuals.l2[component]
                                    / state.steady.reference_l2[component]
                                                               : nan);
    }
    for (std::size_t component = 0; component < state.residuals.model_l2.size(); ++component) {
        history_stream_ << ' '
                        << (state.steady.reference_initialized
                                ? state.residuals.model_l2[component]
                                    / state.steady.model_reference_l2[component]
                                : nan);
    }
    for (const Real value : state.residuals.linf)
        history_stream_ << ' ' << value;
    for (const Real value : state.residuals.model_linf)
        history_stream_ << ' ' << value;
    for (const Real value : state.steady.reference_linf) {
        history_stream_ << ' ' << (state.steady.reference_initialized ? value : nan);
    }
    for (std::size_t component = 0; component < expected_model_components; ++component) {
        history_stream_ << ' '
                        << (state.steady.reference_initialized
                                ? state.steady.model_reference_linf[component]
                                : nan);
    }
    for (std::size_t component = 0; component < euler_components; ++component) {
        history_stream_ << ' '
                        << (state.steady.reference_initialized ? state.residuals.linf[component]
                                    / state.steady.reference_linf[component]
                                                               : nan);
    }
    for (std::size_t component = 0; component < state.residuals.model_linf.size(); ++component) {
        history_stream_ << ' '
                        << (state.steady.reference_initialized
                                ? state.residuals.model_linf[component]
                                    / state.steady.model_reference_linf[component]
                                : nan);
    }
    history_stream_ << ' ' << state.steady.consecutive_passes << ' '
                    << state.diagnostics.reconstruction_fallbacks << ' '
                    << state.diagnostics.riemann_fallbacks;
    for (const auto count : state.diagnostics.robustness_face_levels) {
        history_stream_ << ' ' << count;
    }
    const auto finite_or_nan = [nan](Real value) { return std::isfinite(value) ? value : nan; };
    history_stream_ << ' ' << state.diagnostics.robustness_troubled_cells << ' '
                    << state.diagnostics.robustness_local_recomputations << ' '
                    << state.diagnostics.robustness_step_retries << ' '
                    << state.diagnostics.robustness_proposed_time_step << ' '
                    << state.diagnostics.robustness_accepted_time_step << ' '
                    << finite_or_nan(state.diagnostics.minimum_density) << ' '
                    << finite_or_nan(state.diagnostics.minimum_pressure) << ' '
                    << finite_or_nan(state.diagnostics.minimum_temperature) << ' '
                    << finite_or_nan(state.diagnostics.minimum_internal_energy) << ' '
                    << (residual_checked ? 1 : 0) << ' ';
    if (config_.output.history.format == SeriesOutputFormat::Tecplot) {
        history_stream_ << static_cast<int>(state.stop_reason);
    } else {
        history_stream_ << stop_reason_name(state.stop_reason);
    }
    history_stream_ << '\n';
    if (!history_stream_) {
        throw std::runtime_error("failed to write residual history");
    }
}

void RuntimeOutputManager::on_initial(const SimulationState& state)
{
    prepare_directory();
    if (config_.output.history.enabled && history_schedule_.consume(state, true, false)) {
        write_history(state, false);
    }
    dispatch(
        OutputCategory::Field, field_schedule_, state, true, false, config_.output.field.enabled);
    if (config_.output.statistics.enabled && statistics_schedule_.consume(state, true, false)) {
        write_statistics(state);
    }
    dispatch(OutputCategory::Boundary,
             boundary_schedule_,
             state,
             true,
             false,
             config_.output.boundary.enabled);
    dispatch(OutputCategory::Checkpoint,
             checkpoint_schedule_,
             state,
             true,
             false,
             config_.output.checkpoint.enabled);
}

void RuntimeOutputManager::on_step(const SimulationState& state, bool residual_checked)
{
    if (config_.output.history.enabled && history_schedule_.consume(state, false, false)) {
        write_history(state, residual_checked);
    }
    if (state.stop_reason == StopReason::NumericalFailure) return;
    dispatch(
        OutputCategory::Field, field_schedule_, state, false, false, config_.output.field.enabled);
    if (config_.output.statistics.enabled && statistics_schedule_.consume(state, false, false)) {
        write_statistics(state);
    }
    dispatch(OutputCategory::Boundary,
             boundary_schedule_,
             state,
             false,
             false,
             config_.output.boundary.enabled);
    dispatch(OutputCategory::Checkpoint,
             checkpoint_schedule_,
             state,
             false,
             false,
             config_.output.checkpoint.enabled);
}

void RuntimeOutputManager::finish_history()
{
    if (!config_.output.history.enabled || mpi_.rank() != 0 || !history_stream_.is_open()) {
        return;
    }
    history_stream_.flush();
    if (!history_stream_) {
        throw std::runtime_error("failed to flush residual history");
    }
    history_stream_.close();
    atomic_replace(history_temporary_path_, history_final_path_, config_.output.allow_existing);
    record_file(history_final_path_);
}

void RuntimeOutputManager::write_statistics(const SimulationState& state)
{
    if (statistic_context_ == nullptr) {
        throw std::runtime_error("statistics context is not configured");
    }
    std::vector<Real> values;
    values.reserve(config_.output.statistics.quantities.size());
    for (const auto& name : config_.output.statistics.quantities) {
        values.push_back(statistic_registry_.evaluate(name, *statistic_context_));
    }
    if (mpi_.rank() != 0) return;
    if (!statistics_stream_) {
        throw std::runtime_error("statistics stream is not writable");
    }
    statistics_stream_ << state.step << ' ' << std::setprecision(17) << state.time;
    for (const Real value : values)
        statistics_stream_ << ' ' << value;
    statistics_stream_ << '\n';
    if (!statistics_stream_) throw std::runtime_error("failed to write statistics");
}

void RuntimeOutputManager::finish_statistics()
{
    if (!config_.output.statistics.enabled || mpi_.rank() != 0 || !statistics_stream_.is_open()) {
        return;
    }
    statistics_stream_.flush();
    if (!statistics_stream_) throw std::runtime_error("failed to flush statistics");
    statistics_stream_.close();
    atomic_replace(
        statistics_temporary_path_, statistics_final_path_, config_.output.allow_existing);
    record_file(statistics_final_path_);
}

void RuntimeOutputManager::write_manifest(const SimulationState& state)
{
    if (mpi_.rank() != 0) return;
    const auto target = join_path(output_directory_,
                                  safe_name(config_.case_name) + ".manifest.r"
                                      + std::to_string(mpi_.size()) + ".txt");
    const auto temporary = target + ".tmp";
    std::ofstream output(temporary, std::ios::out | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("cannot open manifest temporary file: " + temporary);
    }
    output << "manifest_version=1\n"
           << "program_version=" << WCNS_PROGRAM_VERSION << '\n'
           << "git_commit=" << WCNS_GIT_COMMIT << '\n'
           << "compiler=" << __VERSION__ << '\n'
#ifdef NDEBUG
           << "build_type=Release\n"
#else
           << "build_type=Debug\n"
#endif
           << "case=" << config_.case_name << '\n'
           << "mpi_ranks=" << mpi_.size() << '\n'
           << "config_digest=" << config_.digest() << '\n'
           << "partition_digest=" << partition_.digest() << '\n'
           << "mesh_signature=" << mesh_signature_ << '\n'
           << "restart_signature=" << config_.restart_signature() << '\n'
           << "step=" << state.step << '\n'
           << std::setprecision(17) << "time=" << state.time << '\n'
           << "time_step=" << state.time_step << '\n'
           << "wall_time=" << state.wall_time << '\n'
           << "stop_reason=" << stop_reason_name(state.stop_reason) << '\n'
           << "config_summary=" << config_.summary() << '\n'
           << "partition_summary=" << partition_.summary() << '\n';
    for (std::size_t level = 0; level < state.diagnostics.robustness_face_levels.size(); ++level) {
        output << "robustness_level" << level
               << "_faces=" << state.diagnostics.robustness_face_levels[level] << '\n';
    }
    output << "robustness_troubled_cells=" << state.diagnostics.robustness_troubled_cells << '\n'
           << "robustness_local_recomputations="
           << state.diagnostics.robustness_local_recomputations << '\n'
           << "robustness_step_retries=" << state.diagnostics.robustness_step_retries << '\n'
           << "robustness_proposed_dt=" << state.diagnostics.robustness_proposed_time_step << '\n'
           << "robustness_accepted_dt=" << state.diagnostics.robustness_accepted_time_step << '\n'
           << "robustness_minimum_rho=" << state.diagnostics.minimum_density << '\n'
           << "robustness_minimum_p=" << state.diagnostics.minimum_pressure << '\n'
           << "robustness_minimum_T=" << state.diagnostics.minimum_temperature << '\n'
           << "robustness_minimum_e=" << state.diagnostics.minimum_internal_energy << '\n';
    for (const auto& file : files_)
        output << "file=" << file << '\n';
    output.close();
    if (!output) throw std::runtime_error("failed to write manifest: " + temporary);
    atomic_replace(temporary, target, config_.output.allow_existing);
    files_.push_back(target);
}

void RuntimeOutputManager::on_final(const SimulationState& state)
{
    if (finalized_) return;
    // Initialization can fail before on_initial() is reached.  The failure
    // manifest is still mandatory evidence, so prepare its directory here too.
    prepare_directory();
    if (config_.output.history.enabled && history_schedule_.consume(state, false, true)) {
        write_history(state, true);
    }
    if (state.stop_reason != StopReason::NumericalFailure) {
        dispatch(OutputCategory::Field,
                 field_schedule_,
                 state,
                 false,
                 true,
                 config_.output.field.enabled);
        if (config_.output.statistics.enabled && statistics_schedule_.consume(state, false, true)) {
            write_statistics(state);
        }
        dispatch(OutputCategory::Boundary,
                 boundary_schedule_,
                 state,
                 false,
                 true,
                 config_.output.boundary.enabled);
        const bool safety_checkpoint = state.stop_reason == StopReason::WallTimeCheckpoint
            || state.stop_reason == StopReason::UserSignalCheckpoint;
        dispatch(OutputCategory::Checkpoint,
                 checkpoint_schedule_,
                 state,
                 false,
                 true,
                 config_.output.checkpoint.enabled || safety_checkpoint);
    }
    finish_history();
    finish_statistics();
    write_manifest(state);
    finalized_ = true;
}

void RuntimeOutputManager::record_file(std::string path)
{
    if (mpi_.rank() == 0) files_.push_back(std::move(path));
}

} // namespace wcns
