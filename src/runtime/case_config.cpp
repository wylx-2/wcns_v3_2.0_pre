#include <wcns/runtime/case_config.hpp>
#include <wcns/runtime/quantity_registry.hpp>
#include <wcns/physics/double_mach_reflection.hpp>
#include <wcns/solver/low_mach_preconditioner.hpp>

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <string_view>
#include <utility>
#include <vector>

namespace wcns {
namespace {

using EntryMap = std::unordered_map<std::string, std::string>;

std::string trim(std::string_view text)
{
    std::size_t first = 0;
    while (first < text.size() && std::isspace(static_cast<unsigned char>(text[first])) != 0) {
        ++first;
    }
    std::size_t last = text.size();
    while (last > first && std::isspace(static_cast<unsigned char>(text[last - 1])) != 0) {
        --last;
    }
    return std::string(text.substr(first, last - first));
}

EntryMap parse_entries(const std::string& text)
{
    EntryMap entries;
    std::istringstream stream(text);
    std::string line;
    std::size_t line_number = 0;
    while (std::getline(stream, line)) {
        ++line_number;
        const auto clean = trim(line);
        if (clean.empty() || clean.front() == '#') {
            continue;
        }
        const auto equal = clean.find('=');
        if (equal == std::string::npos || clean.find('=', equal + 1) != std::string::npos) {
            throw CaseConfigurationError("configuration line " + std::to_string(line_number)
                                         + " must contain exactly one '='");
        }
        const auto key = trim(std::string_view(clean).substr(0, equal));
        const auto value = trim(std::string_view(clean).substr(equal + 1));
        if (key.empty() || value.empty()) {
            throw CaseConfigurationError("configuration line " + std::to_string(line_number)
                                         + " has an empty key or value");
        }
        if (!entries.emplace(key, value).second) {
            throw CaseConfigurationError("duplicate configuration key: " + key);
        }
    }
    return entries;
}

const std::string& require(const EntryMap& entries, const std::string& key)
{
    const auto iterator = entries.find(key);
    if (iterator == entries.end()) {
        throw CaseConfigurationError("missing required configuration key: " + key);
    }
    return iterator->second;
}

Real parse_real(const std::string& value, const std::string& key)
{
    std::size_t consumed = 0;
    Real result = 0.0;
    try {
        result = std::stod(value, &consumed);
    } catch (const std::exception&) {
        throw CaseConfigurationError("configuration key is not a real number: " + key);
    }
    if (consumed != value.size() || !std::isfinite(result)) {
        throw CaseConfigurationError("configuration key requires a finite real number: " + key);
    }
    return result;
}

long long parse_integer(const std::string& value, const std::string& key)
{
    std::size_t consumed = 0;
    long long result = 0;
    try {
        result = std::stoll(value, &consumed);
    } catch (const std::exception&) {
        throw CaseConfigurationError("configuration key is not an integer: " + key);
    }
    if (consumed != value.size()) {
        throw CaseConfigurationError("configuration key is not an integer: " + key);
    }
    return result;
}

bool parse_bool(const std::string& value, const std::string& key)
{
    if (value == "true") return true;
    if (value == "false") return false;
    throw CaseConfigurationError("configuration key requires true or false: " + key);
}

PartitionMode parse_partition_mode(const std::string& value)
{
    if (value == "zones_only") return PartitionMode::ZonesOnly;
    if (value == "auto_split") return PartitionMode::AutoSplit;
    if (value == "force_split") return PartitionMode::ForceSplit;
    throw CaseConfigurationError("unknown partition mode: " + value);
}

RunMode parse_run_mode(const std::string& value)
{
    if (value == "steady") return RunMode::Steady;
    if (value == "unsteady") return RunMode::Unsteady;
    throw CaseConfigurationError("unknown run mode: " + value);
}

RestartMode parse_restart_mode(const std::string& value)
{
    if (value == "strict") return RestartMode::Strict;
    if (value == "algorithm_change") return RestartMode::AlgorithmChange;
    throw CaseConfigurationError("unknown restart mode: " + value);
}

TimeIntegratorKind parse_time_integrator(const std::string& value)
{
    if (value == "ssprk3") return TimeIntegratorKind::SspRk3;
    if (value == "lu_sgs") return TimeIntegratorKind::LuSgs;
    throw CaseConfigurationError("unknown time integrator: " + value);
}

PreconditionerKind parse_preconditioner(const std::string& value)
{
    if (value == "none") return PreconditionerKind::None;
    if (value == "weiss_smith") return PreconditionerKind::WeissSmith;
    throw CaseConfigurationError("unknown preconditioner: " + value);
}

TurbulenceModelKind parse_turbulence_model(const std::string& value)
{
    try {
        return turbulence_model_kind(value);
    } catch (const std::invalid_argument&) {
        throw CaseConfigurationError("unknown turbulence model: " + value);
    }
}

TurbulenceSourceTreatment parse_turbulence_source_treatment(const std::string& value)
{
    try {
        return turbulence_source_treatment(value);
    } catch (const std::invalid_argument&) {
        throw CaseConfigurationError("unknown turbulence source treatment: " + value);
    }
}

WallTreatment parse_wall_treatment(const std::string& value)
{
    try {
        return wall_treatment(value);
    } catch (const std::invalid_argument&) {
        throw CaseConfigurationError("unknown turbulence wall treatment: " + value);
    }
}

LesWallDamping parse_les_wall_damping(const std::string& value)
{
    try {
        return les_wall_damping(value);
    } catch (const std::invalid_argument&) {
        throw CaseConfigurationError("unknown LES wall damping: " + value);
    }
}

FieldOutputFormat parse_field_output_format(const std::string& value)
{
    if (value == "cgns") return FieldOutputFormat::Cgns;
    if (value == "tecplot") return FieldOutputFormat::Tecplot;
    if (value == "both") return FieldOutputFormat::Both;
    throw CaseConfigurationError("unknown field output format: " + value);
}

SeriesOutputFormat parse_series_output_format(const std::string& value)
{
    if (value == "txt") return SeriesOutputFormat::Text;
    if (value == "tecplot") return SeriesOutputFormat::Tecplot;
    throw CaseConfigurationError("unknown series output format: " + value);
}

ReconstructionVariables parse_reconstruction_variables(const std::string& value)
{
    if (value == "conservative") return ReconstructionVariables::Conservative;
    if (value == "primitive") return ReconstructionVariables::Primitive;
    if (value == "characteristic") return ReconstructionVariables::Characteristic;
    throw CaseConfigurationError("unknown reconstruction variable space: " + value);
}

FluxDifferenceMode parse_flux_difference_mode(const std::string& value)
{
    if (value == "profile") return FluxDifferenceMode::Profile;
    if (value == "conservative_two_point") {
        return FluxDifferenceMode::ConservativeTwoPoint;
    }
    throw CaseConfigurationError("unknown flux-difference mode: " + value);
}

BoundaryType parse_boundary_type(const std::string& value)
{
    if (value == "farfield") return BoundaryType::Farfield;
    if (value == "inflow") return BoundaryType::Inflow;
    if (value == "outflow") return BoundaryType::Outflow;
    if (value == "slip_wall") return BoundaryType::SlipWall;
    if (value == "no_slip_adiabatic_wall") return BoundaryType::NoSlipAdiabaticWall;
    if (value == "no_slip_isothermal_wall") return BoundaryType::NoSlipIsothermalWall;
    if (value == "symmetry") return BoundaryType::Symmetry;
    if (value == "periodic") return BoundaryType::Periodic;
    if (value == "double_mach_reflection") {
        return BoundaryType::DoubleMachReflection;
    }
    throw CaseConfigurationError("unknown boundary type: " + value);
}

std::vector<std::string> split_list(const std::string& value)
{
    std::vector<std::string> result;
    std::size_t begin = 0;
    while (begin <= value.size()) {
        const auto end = value.find(',', begin);
        const auto item = trim(std::string_view(value).substr(
            begin, end == std::string::npos ? value.size() - begin : end - begin));
        if (item.empty()) {
            throw CaseConfigurationError("list contains an empty item");
        }
        result.push_back(item);
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return result;
}

SourceModelKind parse_source_model(const std::string& value)
{
    if (value == "uniform_conservative") return SourceModelKind::UniformConservative;
    if (value == "body_force") return SourceModelKind::BodyForce;
    if (value == "ramp_trip") return SourceModelKind::RampTrip;
    if (value == "pressure_gradient") return SourceModelKind::PressureGradient;
    if (value == "manufactured") return SourceModelKind::ManufacturedSolution;
    throw CaseConfigurationError("unknown source model: " + value);
}

std::uint64_t fnv1a(const std::string& text)
{
    std::uint64_t result = 14695981039346656037ull;
    for (const unsigned char byte : text) {
        result ^= static_cast<std::uint64_t>(byte);
        result *= 1099511628211ull;
    }
    return result;
}

std::string canonical_entries(const EntryMap& entries)
{
    std::vector<std::pair<std::string, std::string>> sorted(entries.begin(), entries.end());
    std::sort(sorted.begin(), sorted.end());
    std::ostringstream result;
    for (const auto& [key, value] : sorted) {
        result << key << '=' << value << '\n';
    }
    return result.str();
}

const std::set<std::string>& fixed_keys()
{
    static const std::set<std::string> keys {
        "hill.enabled",
        "hill.bulk_velocity",
        "hill.density",
        "hill.controller_time",
        "hill.force_limit",
        "hill.statistics.start",
        "hill.statistics.end",
        "hill.write_every_steps",
        "schema_version",
        "case.name",
        "mesh.path",
        "algorithm.profile",
        "algorithm.reconstruction",
        "algorithm.reconstruction_variables",
        "algorithm.riemann",
        "algorithm.roe_all_speed.reference_mach",
        "algorithm.roe_all_speed.dissipation_scale",
        "algorithm.roe_all_speed.pressure_coefficient",
        "algorithm.flux_difference",
        "algorithm.mdcd.disp",
        "algorithm.mdcd.diss",
        "robustness.enabled",
        "robustness.max_local_recomputations",
        "robustness.max_step_retries",
        "robustness.time_step_reduction",
        "robustness.minimum_time_step",
        "transport.model",
        "transport.prandtl",
        "transport.constant.viscosity_ratio",
        "transport.sutherland.reference_viscosity_ratio",
        "transport.sutherland.temperature",
        "transport.sutherland.temperature_ratio",
        "turbulence.model",
        "turbulence.experimental",
        "turbulence.prandtl",
        "turbulence.wall_treatment",
        "turbulence.sa.farfield_nu_tilde_ratio",
        "turbulence.sa.source_treatment",
        "turbulence.freestream.intensity",
        "turbulence.freestream.length_scale",
        "turbulence.model_floor",
        "turbulence.two_equation.source_treatment",
        "turbulence.wall_function.y_plus_min",
        "turbulence.wall_function.y_plus_max",
        "les.filter.type",
        "les.filter.width_ratio",
        "les.test_filter.ratio",
        "les.smagorinsky.cs",
        "les.smagorinsky.wall_damping",
        "les.similarity.cb",
        "les.dynamic.average",
        "les.dynamic.clipping",
        "les.dynamic.denominator_floor",
        "les.dynamic.coefficient_minimum",
        "les.dynamic.coefficient_maximum",
        "les.wale.cw",
        "les.sgs_prandtl",
        "time.integrator",
        "time.physical.scheme",
        "time.physical.step",
        "time.dual_time.max_iterations",
        "time.dual_time.absolute_tolerance",
        "time.dual_time.relative_tolerance",
        "time.dual_time.cfl",
        "lu_sgs.sweeps",
        "lu_sgs.jacobian",
        "lu_sgs.relaxation",
        "preconditioner.type",
        "preconditioner.mach_cutoff",
        "preconditioner.viscous_cutoff",
        "gas.gamma",
        "gas.molar_mass",
        "gas.specific_gas_constant",
        "reference.velocity",
        "reference.density",
        "reference.temperature",
        "reference.length",
        "reference.viscosity",
        "partition.mode",
        "partition.allow_idle_ranks",
        "partition.max_load_ratio",
        "partition.min_cells_per_active_direction",
        "geometry.metric.fallback",
        "geometry.metric.maximum_reference_relative_difference",
        "initial.type",
        "hit.type",
        "hit.n",
        "hit.length",
        "hit.seed",
        "hit.spectrum_file",
        "hit.initial_energy",
        "hit.initialization",
        "hit.spectrum_amplitude",
        "hit.preparation.time",
        "hit.preparation.max_steps",
        "hit.peak_wave",
        "hit.cutoff",
        "hit.forcing",
        "hit.forcing_power",
        "hit.forcing_kmax",
        "hit.remove_mean_acceleration",
        "hit.thermostat",
        "hit.statistics.start",
        "hit.statistics.end",
        "hit.sample_every_steps",
        "hit.write_every_samples",
        "hit.sample_times",
        "benchmark.type",
        "benchmark.statistics.start",
        "benchmark.statistics.end",
        "benchmark.write_every_steps",
        "benchmark.history_every_steps",
        "benchmark.probes",
        "source.ramp_trip.amplitude",
        "source.ramp_trip.span",
        "initial.alpha_degrees",

        "initial.rho",
        "initial.u",
        "initial.v",
        "initial.w",
        "initial.pressure",
        "initial.temperature",
        "initial.x0",
        "initial.y0",
        "initial.y1",
        "initial.lower_velocity",
        "initial.upper_velocity",
        "initial.centerline_velocity",
        "initial.lower_temperature",
        "initial.upper_temperature",
        "initial.temperature_curvature",
        "initial.velocity_curvature",
        "initial.beta",
        "initial.background_u",
        "initial.background_v",
        "initial.period_x",
        "initial.period_y",
        "initial.period_z",
        "initial.z0",
        "initial.re_tau",
        "initial.bulk_velocity",
        "initial.bulk_velocity_plus",
        "initial.perturbation_amplitude",
        "initial.left_rho",
        "initial.left_u",
        "initial.left_v",
        "initial.left_p",
        "initial.right_rho",
        "initial.right_u",
        "initial.right_v",
        "initial.right_p",
        "initial.ne_rho",
        "initial.ne_u",
        "initial.ne_v",
        "initial.ne_p",
        "initial.nw_rho",
        "initial.nw_u",
        "initial.nw_v",
        "initial.nw_p",
        "initial.sw_rho",
        "initial.sw_u",
        "initial.sw_v",
        "initial.sw_p",
        "initial.se_rho",
        "initial.se_u",
        "initial.se_v",
        "initial.se_p",
        "boundary.default",
        "boundary.point_vortex.enabled",
        "boundary.point_vortex.lift_coefficient",
        "boundary.point_vortex.center_x",
        "boundary.point_vortex.center_y",
        "boundary.point_vortex.chord",
        "source.enabled",
        "source.models",
        "source.uniform.rho",
        "source.uniform.momentum_x",
        "source.uniform.momentum_y",
        "source.uniform.momentum_z",
        "source.uniform.energy",
        "source.body.ax",
        "source.body.ay",
        "source.body.az",
        "source.pressure_gradient.x",
        "source.pressure_gradient.y",
        "source.pressure_gradient.z",
        "source.manufactured.rho",
        "source.manufactured.momentum_x",
        "source.manufactured.momentum_y",
        "source.manufactured.momentum_z",
        "source.manufactured.energy",
        "run.mode",
        "run.viscous",
        "run.cfl",
        "run.max_steps",
        "run.t_end",
        "run.max_wall_time",
        "steady.min_steps",
        "steady.check_interval_steps",
        "steady.consecutive_checks",
        "steady.reference_floor",
        "steady.l2_absolute",
        "steady.l2_relative",
        "steady.linf_enabled",
        "steady.linf_absolute",
        "steady.linf_relative",
        "output.directory",
        "output.allow_existing",
        "output.dimensional",
        "output.field.enabled",
        "output.field.format",
        "output.field.every_steps",
        "output.field.every_time",
        "output.field.explicit_times",
        "output.field.write_initial",
        "output.field.write_final",
        "output.field.quantities",
        "output.history.enabled",
        "output.history.format",
        "output.history.every_steps",
        "output.history.every_time",
        "output.history.explicit_times",
        "output.history.write_initial",
        "output.history.write_final",
        "output.history.quantities",
        "output.statistics.enabled",
        "output.statistics.format",
        "output.statistics.every_steps",
        "output.statistics.every_time",
        "output.statistics.explicit_times",
        "output.statistics.write_initial",
        "output.statistics.write_final",
        "output.statistics.quantities",
        "statistics.time.enabled",
        "statistics.time.start",
        "statistics.time.end",
        "statistics.time.every_steps",
        "statistics.time.weight",
        "output.boundary.enabled",
        "output.boundary.format",
        "output.boundary.every_steps",
        "output.boundary.every_time",
        "output.boundary.explicit_times",
        "output.boundary.write_initial",
        "output.boundary.write_final",
        "output.boundary.patches",
        "output.boundary.quantities",
        "output.boundary.reference_pressure",
        "output.boundary.reference_density",
        "output.boundary.reference_velocity_x",
        "output.boundary.reference_velocity_y",
        "output.boundary.reference_velocity_z",
        "output.boundary.reference_area",
        "output.boundary.reference_length",
        "output.boundary.moment_center_x",
        "output.boundary.moment_center_y",
        "output.boundary.moment_center_z",
        "output.boundary.drag_direction_x",
        "output.boundary.drag_direction_y",
        "output.boundary.drag_direction_z",
        "output.boundary.lift_direction_x",
        "output.boundary.lift_direction_y",
        "output.boundary.lift_direction_z",
        "output.boundary.tangent_direction_x",
        "output.boundary.tangent_direction_y",
        "output.boundary.tangent_direction_z",
        "output.boundary.span_direction_x",
        "output.boundary.span_direction_y",
        "output.boundary.span_direction_z",
        "output.boundary.span_bin_edges",
        "output.statistics.xz_planes.enabled",
        "output.statistics.xz_planes.cell_j_indices",
        "output.statistics.yz_planes.enabled",
        "output.statistics.yz_planes.target_x_coordinates",
        "output.statistics.channel_walls.enabled",
        "output.statistics.channel_walls.lower_patch",
        "output.statistics.channel_walls.upper_patch",
        "output.statistics.channel_walls.half_height",
        "output.checkpoint.enabled",
        "output.checkpoint.every_steps",
        "output.checkpoint.every_time",
        "output.checkpoint.explicit_times",
        "output.checkpoint.write_initial",
        "output.checkpoint.write_final",
        "restart.path",
        "restart.mode",
    };
    return keys;
}

bool dynamic_boundary_key(const std::string& key)
{
    constexpr std::string_view prefix = "boundary.";
    if (key.size() <= prefix.size()
        || key.compare(0, prefix.size(), prefix.data(), prefix.size()) != 0) {
        return false;
    }
    const auto separator = key.rfind('.');
    if (separator == std::string::npos || separator <= prefix.size()
        || separator + 1 >= key.size()) {
        return false;
    }
    static const std::set<std::string> properties {
        "type",
        "wall_velocity_x",
        "wall_velocity_y",
        "wall_velocity_z",
        "wall_temperature",
        "rho",
        "u",
        "v",
        "w",
        "temperature",
        "pressure",
    };
    return properties.find(key.substr(separator + 1)) != properties.end();
}

std::pair<std::string, std::string> split_boundary_key(const std::string& key)
{
    constexpr std::size_t prefix_size = std::string_view("boundary.").size();
    const auto separator = key.rfind('.');
    return {
        key.substr(prefix_size, separator - prefix_size),
        key.substr(separator + 1),
    };
}

void reject_unknown_keys(const EntryMap& entries)
{
    for (const auto& [key, value] : entries) {
        static_cast<void>(value);
        if (key == "Re" || key == "Ma" || key == "reference.reynolds" || key == "reference.mach") {
            throw CaseConfigurationError("Re and Ma are derived values and cannot be configured: "
                                         + key);
        }
        if (fixed_keys().find(key) == fixed_keys().end() && !dynamic_boundary_key(key)) {
            throw CaseConfigurationError("unknown configuration key: " + key);
        }
    }
}

Real optional_real(const EntryMap& entries, const std::string& key, Real default_value)
{
    const auto iterator = entries.find(key);
    return iterator == entries.end() ? default_value : parse_real(iterator->second, key);
}

bool optional_bool(const EntryMap& entries, const std::string& key, bool default_value)
{
    const auto iterator = entries.find(key);
    return iterator == entries.end() ? default_value : parse_bool(iterator->second, key);
}

std::size_t
optional_size(const EntryMap& entries, const std::string& key, std::size_t default_value)
{
    const auto iterator = entries.find(key);
    if (iterator == entries.end()) return default_value;
    const auto value = parse_integer(iterator->second, key);
    if (value < 0
        || static_cast<unsigned long long>(value) > std::numeric_limits<std::size_t>::max()) {
        throw CaseConfigurationError("configuration key is outside size_t range: " + key);
    }
    return static_cast<std::size_t>(value);
}

std::vector<Real> optional_real_list(const EntryMap& entries, const std::string& key)
{
    const auto iterator = entries.find(key);
    if (iterator == entries.end()) return {};
    std::vector<Real> result;
    for (const auto& item : split_list(iterator->second)) {
        result.push_back(parse_real(item, key));
    }
    return result;
}

std::vector<std::string> optional_string_list(const EntryMap& entries, const std::string& key)
{
    const auto iterator = entries.find(key);
    return iterator == entries.end() ? std::vector<std::string> {} : split_list(iterator->second);
}

std::string
optional_string(const EntryMap& entries, const std::string& key, std::string_view default_value)
{
    const auto iterator = entries.find(key);
    return iterator == entries.end() ? std::string(default_value) : iterator->second;
}

std::vector<int> optional_integer_list(const EntryMap& entries, const std::string& key)
{
    const auto iterator = entries.find(key);
    if (iterator == entries.end()) return {};
    std::vector<int> result;
    for (const auto& item : split_list(iterator->second)) {
        const auto value = parse_integer(item, key);
        if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max()) {
            throw CaseConfigurationError("configuration list item is outside int range: " + key);
        }
        result.push_back(static_cast<int>(value));
    }
    return result;
}

OutputScheduleConfig
parse_schedule(const EntryMap& entries, const std::string& prefix, bool default_final)
{
    OutputScheduleConfig result;
    result.every_steps = optional_size(entries, prefix + ".every_steps", 0);
    result.every_time = optional_real(entries, prefix + ".every_time", 0.0);
    result.explicit_times = optional_real_list(entries, prefix + ".explicit_times");
    result.write_initial = optional_bool(entries, prefix + ".write_initial", false);
    result.write_final = optional_bool(entries, prefix + ".write_final", default_final);
    return result;
}

} // namespace

const char* partition_mode_name(PartitionMode mode)
{
    switch (mode) {
    case PartitionMode::ZonesOnly: return "zones_only";
    case PartitionMode::AutoSplit: return "auto_split";
    case PartitionMode::ForceSplit: return "force_split";
    }
    throw CaseConfigurationError("invalid partition mode");
}

const char* boundary_type_name(BoundaryType type)
{
    switch (type) {
    case BoundaryType::Farfield: return "farfield";
    case BoundaryType::Inflow: return "inflow";
    case BoundaryType::Outflow: return "outflow";
    case BoundaryType::SlipWall: return "slip_wall";
    case BoundaryType::NoSlipAdiabaticWall: return "no_slip_adiabatic_wall";
    case BoundaryType::NoSlipIsothermalWall: return "no_slip_isothermal_wall";
    case BoundaryType::Symmetry: return "symmetry";
    case BoundaryType::Periodic: return "periodic";
    case BoundaryType::DoubleMachReflection: return "double_mach_reflection";
    case BoundaryType::Undefined: return "undefined";
    }
    throw CaseConfigurationError("invalid boundary type");
}

const char* run_mode_name(RunMode mode)
{
    switch (mode) {
    case RunMode::Steady: return "steady";
    case RunMode::Unsteady: return "unsteady";
    }
    throw CaseConfigurationError("invalid run mode");
}

const char* restart_mode_name(RestartMode mode)
{
    switch (mode) {
    case RestartMode::Strict: return "strict";
    case RestartMode::AlgorithmChange: return "algorithm_change";
    }
    throw std::invalid_argument("unknown restart mode");
}

const char* time_integrator_name(TimeIntegratorKind integrator)
{
    switch (integrator) {
    case TimeIntegratorKind::SspRk3: return "ssprk3";
    case TimeIntegratorKind::LuSgs: return "lu_sgs";
    }
    throw CaseConfigurationError("invalid time integrator");
}

const char* preconditioner_name(PreconditionerKind preconditioner)
{
    switch (preconditioner) {
    case PreconditionerKind::None: return "none";
    case PreconditionerKind::WeissSmith: return "weiss_smith";
    }
    throw CaseConfigurationError("invalid preconditioner");
}

const char* field_output_format_name(FieldOutputFormat format)
{
    switch (format) {
    case FieldOutputFormat::Cgns: return "cgns";
    case FieldOutputFormat::Tecplot: return "tecplot";
    case FieldOutputFormat::Both: return "both";
    }
    throw CaseConfigurationError("invalid field output format");
}

const char* series_output_format_name(SeriesOutputFormat format)
{
    switch (format) {
    case SeriesOutputFormat::Text: return "txt";
    case SeriesOutputFormat::Tecplot: return "tecplot";
    }
    throw CaseConfigurationError("invalid series output format");
}

void PartitionConfig::validate(AlgorithmProfileKind profile) const
{
    if (!std::isfinite(max_load_ratio) || max_load_ratio < 1.0) {
        throw CaseConfigurationError("partition max load ratio must be finite and >= 1");
    }
    const int strict_minimum = profile == AlgorithmProfileKind::PhengleiWcns ? 4 : 5;
    if (min_cells_per_active_direction < strict_minimum) {
        throw CaseConfigurationError(
            "partition minimum cells are below the selected profile limit");
    }
}

std::string PartitionConfig::summary() const
{
    std::ostringstream result;
    result << "partition(mode=" << partition_mode_name(mode)
           << ",allow_idle=" << (allow_idle_ranks ? "true" : "false")
           << ",max_load_ratio=" << std::setprecision(17) << max_load_ratio
           << ",min_cells=" << min_cells_per_active_direction << ')';
    return result.str();
}

Real InitialConditionConfig::parameter(const std::string& name, Real default_value) const
{
    const auto iterator = parameters.find(name);
    return iterator == parameters.end() ? default_value : iterator->second;
}

void InitialConditionConfig::validate(int dimension) const
{
    static const std::set<std::string> valid_types {
        "uniform",
        "quadrant_riemann",
        "sod_x",
        "isentropic_vortex",
        "couette",
        "poiseuille",
        "linear_conduction",
        "manufactured_periodic",
        "taylor_green_vortex",
        "double_mach_reflection",
        "hit",
        "turbulent_channel",
        "periodic_hill",
        "sd7003",
        "compression_ramp",
    };
    if (valid_types.find(type) == valid_types.end()) {
        throw CaseConfigurationError("unknown initial condition type: " + type);
    }
    if (dimension != 2 && dimension != 3) {
        throw CaseConfigurationError("initial condition dimension must be 2 or 3");
    }
    for (const auto& [name, value] : parameters) {
        if (!std::isfinite(value)) {
            throw CaseConfigurationError("initial condition parameter is not finite: " + name);
        }
    }
    if (parameter("rho", 1.0) <= 0.0
        || (parameters.find("temperature") != parameters.end()
            && parameter("temperature", 0.0) <= 0.0)
        || (parameters.find("pressure") != parameters.end() && parameter("pressure", 0.0) <= 0.0)) {
        throw CaseConfigurationError("initial density, temperature and pressure must be positive");
    }
    if (type == "isentropic_vortex"
        && (parameter("period_x", 0.0) < 0.0 || parameter("period_y", 0.0) < 0.0)) {
        throw CaseConfigurationError("isentropic-vortex periods must be zero or positive");
    }
    if (type == "periodic_hill") {
        if (dimension != 3 || parameter("bulk_velocity", 1.0) <= 0
            || parameter("perturbation_amplitude", 0.05) < 0
            || parameter("perturbation_amplitude", 0.05) > 0.2)
            throw CaseConfigurationError("periodic_hill requires 3D, positive bulk velocity and perturbation in [0,0.2]");
    }
    if (type == "turbulent_channel") {
        const Real y0 = parameter("y0", -1.0);
        const Real y1 = parameter("y1", 1.0);
        if (dimension != 3 || !(y1 > y0) || parameter("re_tau", 0.0) <= 0.0
            || parameter("bulk_velocity", 0.0) <= 0.0 || parameter("bulk_velocity_plus", 0.0) <= 0.0
            || parameter("period_x", 2.0 * 3.14159265358979323846) <= 0.0
            || parameter("period_z", 3.14159265358979323846) <= 0.0
            || parameter("perturbation_amplitude", 0.05) < 0.0
            || parameter("perturbation_amplitude", 0.05) > 0.5) {
            throw CaseConfigurationError(
                "turbulent-channel initial data require 3D, y1>y0, positive "
                "Re_tau/bulk scales/periods and perturbation amplitude in [0,0.5]");
        }
    }
    if (type == "taylor_green_vortex"
        && (dimension != 3 || parameter("beta", 0.1) <= 0.0
            || parameter("beta", 0.1) > 1.0)) {
        throw CaseConfigurationError(
            "Taylor-Green initial data require 3D and beta in (0,1]");
    }
    if (type == "couette" || type == "poiseuille" || type == "linear_conduction") {
        const Real y0 = parameter("y0", 0.0);
        const Real y1 = parameter("y1", 1.0);
        if (!(y1 > y0) || parameter("pressure", 1.0) <= 0.0) {
            throw CaseConfigurationError(
                "analytic viscous initial bounds and pressure are invalid");
        }
        if (type == "poiseuille") {
            const Real wall_temperature = parameter("temperature", 1.0);
            const Real thermal_amplitude = parameter("temperature_curvature", 0.0);
            if (wall_temperature <= 0.0
                || wall_temperature + std::min(Real {0.0}, thermal_amplitude) / 48.0 <= 0.0) {
                throw CaseConfigurationError(
                    "Poiseuille initial temperature profile is not positive");
            }
            return;
        }
        const Real lower_temperature
            = parameter("lower_temperature", parameter("temperature", 1.0));
        const Real upper_temperature
            = parameter("upper_temperature", parameter("temperature", 1.0));
        if (lower_temperature <= 0.0 || upper_temperature <= 0.0) {
            throw CaseConfigurationError("analytic viscous initial temperatures are invalid");
        }
        if (type == "couette") {
            const Real curvature = parameter("temperature_curvature", 0.0);
            // T(eta)=T0+(T1-T0)eta+c*eta*(1-eta), eta in [0,1].
            if (curvature < 0.0) {
                const Real vertex
                    = 0.5 * (1.0 + (upper_temperature - lower_temperature) / curvature);
                if (vertex > 0.0 && vertex < 1.0) {
                    const Real temperature = lower_temperature
                        + (upper_temperature - lower_temperature) * vertex
                        + curvature * vertex * (1.0 - vertex);
                    if (temperature <= 0.0) {
                        throw CaseConfigurationError(
                            "Couette analytic temperature is not positive");
                    }
                }
            }
        }
    }
}

std::string InitialConditionConfig::summary() const
{
    std::vector<std::pair<std::string, Real>> sorted(parameters.begin(), parameters.end());
    std::sort(sorted.begin(), sorted.end());
    std::ostringstream result;
    result << "initial(type=" << type;
    for (const auto& [name, value] : sorted) {
        result << ',' << name << '=' << std::setprecision(17) << value;
    }
    result << ')';
    return result.str();
}

bool BoundaryPhysicalDataConfig::has_wall_velocity() const noexcept
{
    return std::any_of(wall_velocity.begin(), wall_velocity.end(), [](const auto& value) {
        return value.has_value();
    });
}

bool BoundaryPhysicalDataConfig::has_target_state() const noexcept
{
    return rho || u || v || w || temperature || pressure;
}

void BoundaryPhysicalDataConfig::validate() const
{
    const auto require_finite = [](const std::optional<Real>& value, const char* name) {
        if (value && !std::isfinite(*value)) {
            throw CaseConfigurationError(std::string("boundary physical value is not finite: ")
                                         + name);
        }
    };
    for (std::size_t component = 0; component < wall_velocity.size(); ++component) {
        require_finite(wall_velocity[component], "wall_velocity");
    }
    require_finite(wall_temperature, "wall_temperature");
    require_finite(rho, "rho");
    require_finite(u, "u");
    require_finite(v, "v");
    require_finite(w, "w");
    require_finite(temperature, "temperature");
    require_finite(pressure, "pressure");
    if (wall_temperature && *wall_temperature <= 0.0) {
        throw CaseConfigurationError("boundary wall_temperature must be positive");
    }
    if (rho && *rho <= 0.0) {
        throw CaseConfigurationError("boundary target rho must be positive");
    }
    if (temperature && *temperature <= 0.0) {
        throw CaseConfigurationError("boundary target temperature must be positive");
    }
    if (pressure && *pressure <= 0.0) {
        throw CaseConfigurationError("boundary target pressure must be positive");
    }
    if (temperature && pressure) {
        throw CaseConfigurationError(
            "boundary target must specify temperature or pressure, not both");
    }
    if (has_target_state() && (!rho || (!temperature && !pressure))) {
        throw CaseConfigurationError(
            "boundary target requires rho and exactly one of temperature or pressure");
    }
}

std::string BoundaryPhysicalDataConfig::summary() const
{
    std::ostringstream result;
    result << std::setprecision(17);
    const auto append = [&](const char* name, const std::optional<Real>& value) {
        if (value) result << ',' << name << '=' << *value;
    };
    append("wall_velocity_x", wall_velocity[0]);
    append("wall_velocity_y", wall_velocity[1]);
    append("wall_velocity_z", wall_velocity[2]);
    append("wall_temperature", wall_temperature);
    append("rho", rho);
    append("u", u);
    append("v", v);
    append("w", w);
    append("temperature", temperature);
    append("pressure", pressure);
    auto text = result.str();
    return text.empty() ? text : text.substr(1);
}

void FarfieldPointVortexConfig::validate(int schema_version) const
{
    if (enabled && schema_version != 2) {
        throw CaseConfigurationError("farfield point vortex requires schema_version=2");
    }
    if (!std::isfinite(lift_coefficient) || !std::isfinite(center[0])
        || !std::isfinite(center[1]) || !std::isfinite(chord) || chord <= 0.0) {
        throw CaseConfigurationError("farfield point vortex parameters are invalid");
    }
}

std::string FarfieldPointVortexConfig::summary() const
{
    std::ostringstream result;
    result << std::setprecision(17) << "point_vortex(enabled=" << (enabled ? "true" : "false")
           << ",Cl=" << lift_coefficient << ",center=" << center[0] << ':' << center[1]
           << ",chord=" << chord << ')';
    return result.str();
}

std::string FarfieldPointVortexConfig::restart_signature() const
{
    return "farfield_" + summary();
}

void SteadyConvergenceConfig::validate() const
{
    if (min_steps == 0 || check_interval_steps == 0 || consecutive_checks == 0) {
        throw CaseConfigurationError("steady step/check counts must be positive");
    }
    const std::array<std::pair<const char*, Real>, 5> values {{
        {"reference_floor", reference_floor},
        {"l2_absolute", l2_absolute},
        {"l2_relative", l2_relative},
        {"linf_absolute", linf_absolute},
        {"linf_relative", linf_relative},
    }};
    for (const auto& value : values) {
        if (!std::isfinite(value.second) || value.second <= 0.0) {
            throw CaseConfigurationError(std::string("steady ") + value.first
                                         + " must be finite and positive");
        }
    }
}

std::string SteadyConvergenceConfig::summary() const
{
    std::ostringstream result;
    result << "steady(min_steps=" << min_steps << ",check_interval=" << check_interval_steps
           << ",consecutive=" << consecutive_checks << ",reference_floor=" << std::setprecision(17)
           << reference_floor << ",l2_abs=" << l2_absolute << ",l2_rel=" << l2_relative
           << ",linf_enabled=" << (linf_enabled ? "true" : "false") << ",linf_abs=" << linf_absolute
           << ",linf_rel=" << linf_relative << ')';
    return result.str();
}

void OutputScheduleConfig::validate() const
{
    if (!std::isfinite(every_time) || every_time < 0.0) {
        throw CaseConfigurationError("output schedule every_time must be finite and non-negative");
    }
    Real previous = -1.0;
    for (const Real time : explicit_times) {
        if (!std::isfinite(time) || time < 0.0 || time <= previous) {
            throw CaseConfigurationError(
                "output explicit times must be finite, non-negative and strictly increasing");
        }
        previous = time;
    }
}

std::string OutputScheduleConfig::summary() const
{
    std::ostringstream result;
    result << "schedule(every_steps=" << every_steps << ",every_time=" << std::setprecision(17)
           << every_time << ",explicit_times=";
    for (std::size_t index = 0; index < explicit_times.size(); ++index) {
        if (index != 0) result << ':';
        result << explicit_times[index];
    }
    result << ",initial=" << (write_initial ? "true" : "false")
           << ",final=" << (write_final ? "true" : "false") << ')';
    return result.str();
}

void FieldOutputConfig::validate() const
{
    schedule.validate();
    if (enabled && quantities.empty()) {
        throw CaseConfigurationError("enabled field output requires at least one quantity");
    }
}

std::string FieldOutputConfig::summary() const
{
    std::ostringstream result;
    result << "field(enabled=" << (enabled ? "true" : "false")
           << ",format=" << field_output_format_name(format) << ',' << schedule.summary()
           << ",quantities=";
    for (std::size_t index = 0; index < quantities.size(); ++index) {
        if (index != 0) result << ':';
        result << quantities[index];
    }
    result << ')';
    return result.str();
}

void SeriesOutputConfig::validate(const char* label) const
{
    schedule.validate();
    if (std::string(label) == "history" && !quantities.empty()) {
        throw CaseConfigurationError("residual history uses a fixed schema; "
                                     "output.history.quantities must be omitted");
    }
    if (enabled && quantities.empty() && std::string(label) == "statistics") {
        throw CaseConfigurationError("enabled statistics output requires at least one quantity");
    }
}

std::string SeriesOutputConfig::summary(const char* label) const
{
    std::ostringstream result;
    result << label << "(enabled=" << (enabled ? "true" : "false")
           << ",format=" << series_output_format_name(format) << ',' << schedule.summary()
           << ",quantities=";
    for (std::size_t index = 0; index < quantities.size(); ++index) {
        if (index != 0) result << ':';
        result << quantities[index];
    }
    result << ')';
    return result.str();
}

void BoundaryOutputConfig::validate(bool viscous) const
{
    schedule.validate();
    if (!enabled) return;
    const std::set<std::string> supported {
        "p_w",
        "T_w",
        "mu_w",
        "Cp",
        "Cf",
        "q_wall",
        "wall_distance",
        "friction_velocity",
        "wall_y_plus",
        "wall_y_plus_class",
        "mu_model_over_mu",
        "pressure_traction_x",
        "pressure_traction_y",
        "pressure_traction_z",
        "viscous_traction_x",
        "viscous_traction_y",
        "viscous_traction_z",
        "traction_x",
        "traction_y",
        "traction_z",
    };
    const std::set<std::string> viscous_only {
        "Cf",
        "q_wall",
        "wall_distance",
        "friction_velocity",
        "wall_y_plus",
        "wall_y_plus_class",
        "mu_model_over_mu",
        "viscous_traction_x",
        "viscous_traction_y",
        "viscous_traction_z",
    };
    std::set<std::string> unique_patches;
    for (const auto& patch : patches) {
        if (patch.empty() || !unique_patches.insert(patch).second) {
            throw CaseConfigurationError("boundary output patch names must be nonempty and unique");
        }
    }
    std::set<std::string> unique_quantities;
    for (const auto& quantity : quantities) {
        if (supported.find(quantity) == supported.end()
            || !unique_quantities.insert(quantity).second) {
            throw CaseConfigurationError("boundary output quantity is unknown or duplicated: "
                                         + quantity);
        }
        if (!viscous && viscous_only.find(quantity) != viscous_only.end()) {
            throw CaseConfigurationError(
                "inviscid boundary output cannot request viscous quantity: " + quantity);
        }
    }
    if (patches.empty() || quantities.empty()) {
        throw CaseConfigurationError("enabled boundary output requires patches and quantities");
    }
    const auto positive = [](Real value) { return std::isfinite(value) && value > 0.0; };
    const auto finite_vector = [](const std::array<Real, 3>& vector) {
        return std::all_of(
            vector.begin(), vector.end(), [](Real value) { return std::isfinite(value); });
    };
    const auto dot = [](const std::array<Real, 3>& lhs, const std::array<Real, 3>& rhs) {
        return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
    };
    const auto unit = [&](const std::array<Real, 3>& vector) {
        return finite_vector(vector) && std::abs(dot(vector, vector) - 1.0) <= 1.0e-12;
    };
    if (!positive(reference_pressure) || !positive(reference_density) || !positive(reference_area)
        || !positive(reference_length) || !finite_vector(reference_velocity)
        || !finite_vector(moment_center) || !unit(drag_direction) || !unit(lift_direction)
        || !unit(tangent_direction) || std::abs(dot(drag_direction, lift_direction)) > 1.0e-12) {
        throw CaseConfigurationError("boundary output reference data are invalid");
    }
    if (!span_bin_edges.empty()) {
        if (span_bin_edges.size() < 2 || !unit(span_direction)) {
            throw CaseConfigurationError(
                "boundary span bins require a unit span direction and at least two edges");
        }
        for (std::size_t edge = 0; edge < span_bin_edges.size(); ++edge) {
            if (!std::isfinite(span_bin_edges[edge])
                || (edge > 0 && span_bin_edges[edge] <= span_bin_edges[edge - 1])) {
                throw CaseConfigurationError(
                    "boundary span-bin edges must be finite and strictly increasing");
            }
        }
        const std::array<Real, 3> right_handed {{
            drag_direction[1] * lift_direction[2] - drag_direction[2] * lift_direction[1],
            drag_direction[2] * lift_direction[0] - drag_direction[0] * lift_direction[2],
            drag_direction[0] * lift_direction[1] - drag_direction[1] * lift_direction[0],
        }};
        if (dot(right_handed, span_direction) < 1.0 - 1.0e-12) {
            throw CaseConfigurationError(
                "boundary drag, lift and span directions must form a right-handed basis");
        }
    }
    const Real speed_squared = dot(reference_velocity, reference_velocity);
    if (!std::isfinite(speed_squared) || 0.5 * reference_density * speed_squared <= 1.0e-12) {
        throw CaseConfigurationError("boundary output reference dynamic pressure is too small");
    }
}

std::string BoundaryOutputConfig::summary() const
{
    if (!enabled) return "boundary(enabled=false)";
    std::ostringstream result;
    result << "boundary(enabled=true,format=" << series_output_format_name(format) << ','
           << schedule.summary() << ",patches=";
    for (std::size_t index = 0; index < patches.size(); ++index) {
        if (index != 0) result << ':';
        result << patches[index];
    }
    result << ",quantities=";
    for (std::size_t index = 0; index < quantities.size(); ++index) {
        if (index != 0) result << ':';
        result << quantities[index];
    }
    const auto vector = [&result](const std::array<Real, 3>& value) {
        result << value[0] << ':' << value[1] << ':' << value[2];
    };
    result << std::setprecision(17) << ",p_ref=" << reference_pressure
           << ",rho_ref=" << reference_density << ",u_ref=";
    vector(reference_velocity);
    result << ",A_ref=" << reference_area << ",L_ref=" << reference_length << ",moment_center=";
    vector(moment_center);
    result << ",drag=";
    vector(drag_direction);
    result << ",lift=";
    vector(lift_direction);
    result << ",tangent=";
    vector(tangent_direction);
    if (!span_bin_edges.empty()) {
        result << ",span=";
        vector(span_direction);
        result << ",span_bins=";
        for (std::size_t edge = 0; edge < span_bin_edges.size(); ++edge) {
            if (edge != 0) result << ':';
            result << span_bin_edges[edge];
        }
    }
    result << ')';
    return result.str();
}

void XzPlaneStatisticsConfig::validate(bool statistics_enabled) const
{
    std::set<int> unique;
    for (const int index : cell_j_indices) {
        if (index < 0 || !unique.insert(index).second) {
            throw CaseConfigurationError(
                "x-z plane cell-j indices must be non-negative and unique");
        }
    }
    if (enabled && !statistics_enabled) {
        throw CaseConfigurationError(
            "x-z plane monitoring requires output.statistics.enabled=true");
    }
    if (enabled && cell_j_indices.empty()) {
        throw CaseConfigurationError(
            "enabled x-z plane monitoring requires at least one cell-j index");
    }
}

std::string XzPlaneStatisticsConfig::summary() const
{
    std::ostringstream result;
    result << "xz_planes(enabled=" << (enabled ? "true" : "false") << ",cell_j_indices=";
    for (std::size_t index = 0; index < cell_j_indices.size(); ++index) {
        if (index != 0) result << ':';
        result << cell_j_indices[index];
    }
    result << ')';
    return result.str();
}

void YzPlaneStatisticsConfig::validate(bool statistics_enabled) const
{
    std::set<Real> unique;
    for (const Real coordinate : target_x_coordinates) {
        if (!std::isfinite(coordinate) || !unique.insert(coordinate).second) {
            throw CaseConfigurationError(
                "y-z plane target x coordinates must be finite and unique");
        }
    }
    if (enabled && !statistics_enabled) {
        throw CaseConfigurationError(
            "y-z plane monitoring requires output.statistics.enabled=true");
    }
    if (enabled && target_x_coordinates.empty()) {
        throw CaseConfigurationError(
            "enabled y-z plane monitoring requires at least one target x coordinate");
    }
}

std::string YzPlaneStatisticsConfig::summary() const
{
    std::ostringstream result;
    result << "yz_planes(enabled=" << (enabled ? "true" : "false")
           << ",target_x_coordinates=" << std::setprecision(17);
    for (std::size_t index = 0; index < target_x_coordinates.size(); ++index) {
        if (index != 0) result << ':';
        result << target_x_coordinates[index];
    }
    result << ')';
    return result.str();
}

void ChannelWallStatisticsConfig::validate(bool statistics_enabled) const
{
    if (lower_patch.empty() || upper_patch.empty() || lower_patch == upper_patch
        || !std::isfinite(half_height) || half_height <= 0.0) {
        throw CaseConfigurationError(
            "channel-wall monitoring requires two distinct patch names and "
            "a positive finite half height");
    }
    if (enabled && !statistics_enabled) {
        throw CaseConfigurationError(
            "channel-wall monitoring requires output.statistics.enabled=true");
    }
}

std::string ChannelWallStatisticsConfig::summary() const
{
    std::ostringstream result;
    result << "channel_walls(enabled=" << (enabled ? "true" : "false")
           << ",lower_patch=" << lower_patch << ",upper_patch=" << upper_patch
           << ",half_height=" << std::setprecision(17) << half_height << ')';
    return result.str();
}

void CheckpointOutputConfig::validate() const
{
    schedule.validate();
}

std::string CheckpointOutputConfig::summary() const
{
    return std::string("checkpoint(enabled=") + (enabled ? "true," : "false,") + schedule.summary()
        + ')';
}

void TimeStatisticsConfig::validate(RunMode mode,
                                    bool instantaneous_statistics_enabled) const
{
    if (!enabled) return;
    if (mode != RunMode::Unsteady) {
        throw CaseConfigurationError("time statistics require run.mode=unsteady");
    }
    if (!instantaneous_statistics_enabled) {
        throw CaseConfigurationError(
            "time statistics require output.statistics.enabled=true");
    }
    if (!std::isfinite(start_time) || !std::isfinite(end_time) || start_time < 0.0
        || end_time <= start_time || every_steps == 0) {
        throw CaseConfigurationError(
            "time statistics require 0 <= start < end and every_steps > 0");
    }
}

std::string TimeStatisticsConfig::summary() const
{
    std::ostringstream result;
    result << "time_statistics(enabled=" << (enabled ? "true" : "false");
    if (enabled) {
        result << std::setprecision(17) << ",start=" << start_time << ",end=" << end_time
               << ",every_steps=" << every_steps << ",weight=accepted_dt";
    }
    result << ')';
    return result.str();
}

std::string
TimeStatisticsConfig::restart_signature(const std::vector<std::string>& quantities) const
{
    if (!enabled) return "time_statistics=disabled";
    std::ostringstream result;
    result << std::setprecision(17) << "time_statistics=v1,start=" << start_time
           << ",end=" << end_time << ",every_steps=" << every_steps
           << ",weight=accepted_dt,quantities=";
    for (std::size_t index = 0; index < quantities.size(); ++index) {
        if (index != 0) result << ':';
        result << quantities[index];
    }
    return result.str();
}

void OutputConfig::validate(bool viscous) const
{
    if (directory.empty()) {
        throw CaseConfigurationError("output directory must not be empty");
    }
    field.validate();
    history.validate("history");
    statistics.validate("statistics");
    boundary.validate(viscous);
    xz_planes.validate(statistics.enabled);
    yz_planes.validate(statistics.enabled);
    channel_walls.validate(statistics.enabled);
    checkpoint.validate();
}

std::string OutputConfig::summary() const
{
    std::ostringstream result;
    result << "output(directory=" << directory
           << ",allow_existing=" << (allow_existing ? "true" : "false")
           << ",dimensional=" << (dimensional ? "true" : "false") << ',' << field.summary() << ','
           << history.summary("history") << ',' << statistics.summary("statistics") << ','
           << boundary.summary() << ',' << xz_planes.summary() << ',' << yz_planes.summary() << ','
           << channel_walls.summary() << ',' << checkpoint.summary() << ')';
    return result.str();
}

void TimeAlgorithmConfig::validate() const
{
    static_cast<void>(time_integrator_name(integrator));
    if (!std::isfinite(physical_time_step) || physical_time_step < 0.0
        || dual_time_max_iterations == 0
        || !std::isfinite(dual_time_absolute_tolerance)
        || dual_time_absolute_tolerance <= 0.0
        || !std::isfinite(dual_time_relative_tolerance)
        || dual_time_relative_tolerance <= 0.0 || !std::isfinite(dual_time_cfl)
        || dual_time_cfl <= 0.0 || lu_sgs_sweeps < 1 || lu_sgs_sweeps > 4
        || !std::isfinite(lu_sgs_relaxation) || lu_sgs_relaxation <= 0.0
        || lu_sgs_relaxation > 1.0) {
        throw std::invalid_argument("invalid LU-SGS or dual-time configuration");
    }
}

std::string TimeAlgorithmConfig::summary() const
{
    validate();
    std::ostringstream result;
    result << "time(integrator=" << time_integrator_name(integrator);
    if (integrator == TimeIntegratorKind::LuSgs) {
        result << ",physical_scheme=bdf2,physical_step=" << std::setprecision(17)
               << physical_time_step << ",dual_max=" << dual_time_max_iterations
               << ",dual_abs=" << dual_time_absolute_tolerance
               << ",dual_rel=" << dual_time_relative_tolerance
               << ",dual_cfl=" << dual_time_cfl << ",lu_sgs_sweeps=" << lu_sgs_sweeps
               << ",jacobian=scalar_spectral,relaxation=" << lu_sgs_relaxation;
    }
    result << ')';
    return result.str();
}

std::string TimeAlgorithmConfig::restart_signature() const
{
    return "time_v3;" + summary();
}

void PreconditionerConfig::validate() const
{
    static_cast<void>(preconditioner_name(kind));
    WeissSmithParameters {mach_cutoff, viscous_cutoff}.validate();
}

std::string PreconditionerConfig::summary() const
{
    validate();
    std::ostringstream result;
    result << "preconditioner(type=" << preconditioner_name(kind);
    if (kind == PreconditionerKind::WeissSmith) {
        result << ",mach_cutoff=" << std::setprecision(17) << mach_cutoff
               << ",viscous_cutoff=" << viscous_cutoff;
    }
    result << ')';
    return result.str();
}

std::string PreconditionerConfig::restart_signature() const
{
    return "preconditioner_v2;" + summary();
}

void CaseRunConfig::validate() const
{
    if (!std::isfinite(cfl) || cfl <= 0.0) {
        throw CaseConfigurationError("run CFL must be finite and positive");
    }
    if (max_steps == 0) {
        throw CaseConfigurationError("run max_steps must be positive");
    }
    if (!std::isfinite(end_time) || end_time < 0.0 || !std::isfinite(max_wall_time)
        || max_wall_time < 0.0) {
        throw CaseConfigurationError("run t_end and max_wall_time must be finite and non-negative");
    }
    if (mode == RunMode::Unsteady && end_time <= 0.0) {
        throw CaseConfigurationError("unsteady run requires a positive t_end");
    }
    steady.validate();
}

std::string CaseRunConfig::summary() const
{
    std::ostringstream result;
    result << "run(mode=" << run_mode_name(mode) << ",viscous=" << (viscous ? "true" : "false")
           << ",cfl=" << std::setprecision(17) << cfl << ",max_steps=" << max_steps
           << ",t_end=" << end_time << ",max_wall_time=" << max_wall_time << ',' << steady.summary()
           << ')';
    return result.str();
}

CaseConfig CaseConfig::from_text(const std::string& text)
{
    const auto entries = parse_entries(text);
    reject_unknown_keys(entries);

    CaseConfig result;
    const auto version = parse_integer(require(entries, "schema_version"), "schema_version");
    if (version < minimum_schema_version || version > supported_schema_version) {
        throw CaseConfigurationError("unsupported configuration schema version");
    }
    result.schema_version = static_cast<int>(version);
    const std::array<const char*, 46> v2_keys {{
        "turbulence.model",
        "turbulence.experimental",
        "turbulence.prandtl",
        "turbulence.wall_treatment",
        "turbulence.sa.farfield_nu_tilde_ratio",
        "turbulence.sa.source_treatment",
        "turbulence.freestream.intensity",
        "turbulence.freestream.length_scale",
        "turbulence.model_floor",
        "turbulence.two_equation.source_treatment",
        "turbulence.wall_function.y_plus_min",
        "turbulence.wall_function.y_plus_max",
        "les.filter.type",
        "les.filter.width_ratio",
        "les.test_filter.ratio",
        "les.smagorinsky.cs",
        "les.smagorinsky.wall_damping",
        "les.similarity.cb",
        "les.dynamic.average",
        "les.dynamic.clipping",
        "les.dynamic.denominator_floor",
        "les.dynamic.coefficient_minimum",
        "les.dynamic.coefficient_maximum",
        "les.wale.cw",
        "les.sgs_prandtl",
        "time.integrator",
        "time.physical.scheme",
        "time.physical.step",
        "time.dual_time.max_iterations",
        "time.dual_time.absolute_tolerance",
        "time.dual_time.relative_tolerance",
        "time.dual_time.cfl",
        "lu_sgs.sweeps",
        "lu_sgs.jacobian",
        "lu_sgs.relaxation",
        "preconditioner.type",
        "preconditioner.mach_cutoff",
        "preconditioner.viscous_cutoff",
        "algorithm.roe_all_speed.reference_mach",
        "algorithm.roe_all_speed.dissipation_scale",
        "algorithm.roe_all_speed.pressure_coefficient",
        "statistics.time.enabled",
        "statistics.time.start",
        "statistics.time.end",
        "statistics.time.every_steps",
        "statistics.time.weight",
    }};
    if (version == 1) {
        for (const auto* key : v2_keys) {
            if (entries.find(key) != entries.end()) {
                throw CaseConfigurationError(std::string("schema 1 does not accept v2 key: ")
                                             + key);
            }
        }
    } else {
        result.turbulence.kind
            = parse_turbulence_model(require(entries, "turbulence.model"));
        const auto has_key = [&](const char* key) { return entries.find(key) != entries.end(); };
        const bool has_turbulence_prandtl = entries.find("turbulence.prandtl") != entries.end();
        const bool has_turbulence_experimental
            = entries.find("turbulence.experimental") != entries.end();
        const bool has_wall_treatment
            = entries.find("turbulence.wall_treatment") != entries.end();
        const bool has_sa_farfield
            = entries.find("turbulence.sa.farfield_nu_tilde_ratio") != entries.end();
        const bool has_sa_source
            = entries.find("turbulence.sa.source_treatment") != entries.end();
        const bool has_two_equation_freestream
            = entries.find("turbulence.freestream.intensity") != entries.end()
            || entries.find("turbulence.freestream.length_scale") != entries.end();
        const bool has_model_floor = entries.find("turbulence.model_floor") != entries.end();
        const bool has_two_equation_source
            = entries.find("turbulence.two_equation.source_treatment") != entries.end();
        const bool has_wall_function_range
            = entries.find("turbulence.wall_function.y_plus_min") != entries.end()
            || entries.find("turbulence.wall_function.y_plus_max") != entries.end();
        const std::array<const char*, 13> les_keys {{
            "les.filter.type",
            "les.filter.width_ratio",
            "les.test_filter.ratio",
            "les.smagorinsky.cs",
            "les.smagorinsky.wall_damping",
            "les.similarity.cb",
            "les.dynamic.average",
            "les.dynamic.clipping",
            "les.dynamic.denominator_floor",
            "les.dynamic.coefficient_minimum",
            "les.dynamic.coefficient_maximum",
            "les.wale.cw",
            "les.sgs_prandtl",
        }};
        const bool has_les_key = std::any_of(les_keys.begin(), les_keys.end(), [&](const char* key) {
            return has_key(key);
        });
        const bool two_equation = result.turbulence.kind == TurbulenceModelKind::KOmegaSst
            || result.turbulence.kind == TurbulenceModelKind::KEpsilon;
        const bool les = turbulence_model_family(result.turbulence.kind)
            == TurbulenceModelFamily::LesAlgebraic;
        if (result.turbulence.kind == TurbulenceModelKind::None
            && (has_turbulence_prandtl || has_turbulence_experimental
                || has_wall_treatment || has_sa_farfield
                || has_sa_source || has_two_equation_freestream || has_model_floor
                || has_two_equation_source || has_wall_function_range || has_les_key)) {
            throw CaseConfigurationError(
                "turbulence.model=none does not accept model-specific turbulence keys");
        }
        if (result.turbulence.kind != TurbulenceModelKind::SaNegative
            && (has_sa_farfield || has_sa_source)) {
            throw CaseConfigurationError(
                "turbulence.sa.* keys require turbulence.model=sa_neg");
        }
        if (result.turbulence.kind != TurbulenceModelKind::KEpsilon
            && has_turbulence_experimental) {
            throw CaseConfigurationError(
                "turbulence.experimental is reserved for k_epsilon");
        }
        result.turbulence.experimental = optional_bool(
            entries, "turbulence.experimental", result.turbulence.experimental);
        if (!two_equation
            && (has_two_equation_freestream || has_model_floor || has_two_equation_source
                || has_wall_function_range)) {
            throw CaseConfigurationError(
                "two-equation turbulence keys require k_omega_sst or k_epsilon");
        }
        if (!les && has_les_key) {
            throw CaseConfigurationError("les.* keys require an LES turbulence model");
        }
        if (les) {
            if (has_turbulence_prandtl || has_turbulence_experimental || has_wall_treatment
                || has_sa_farfield || has_sa_source || has_two_equation_freestream
                || has_model_floor || has_two_equation_source || has_wall_function_range) {
                throw CaseConfigurationError(
                    "LES does not accept RANS-specific turbulence keys");
            }
            if (require(entries, "les.filter.type") != "box3_tensor") {
                throw CaseConfigurationError("LES filter type must be box3_tensor");
            }
            result.turbulence.les_filter_width_ratio = parse_real(
                require(entries, "les.filter.width_ratio"), "les.filter.width_ratio");
            result.turbulence.les_sgs_prandtl
                = parse_real(require(entries, "les.sgs_prandtl"), "les.sgs_prandtl");
            const bool similarity = result.turbulence.kind == TurbulenceModelKind::ScaleSimilarity
                || result.turbulence.kind
                    == TurbulenceModelKind::MixedSmagorinskySimilarity;
            const bool smag = result.turbulence.kind == TurbulenceModelKind::Smagorinsky
                || result.turbulence.kind
                    == TurbulenceModelKind::MixedSmagorinskySimilarity;
            const bool dynamic
                = result.turbulence.kind == TurbulenceModelKind::DynamicSmagorinsky;
            const bool uses_test_filter = similarity || dynamic;
            if (uses_test_filter) {
                result.turbulence.les_test_filter_ratio = parse_real(
                    require(entries, "les.test_filter.ratio"), "les.test_filter.ratio");
                if (result.turbulence.les_test_filter_ratio != 2.0) {
                    throw CaseConfigurationError(
                        "stage Z fixes les.test_filter.ratio at 2");
                }
            } else if (has_key("les.test_filter.ratio")) {
                throw CaseConfigurationError(
                    "les.test_filter.ratio requires a similarity or dynamic model");
            }
            if (smag) {
                result.turbulence.les_smagorinsky_coefficient = parse_real(
                    require(entries, "les.smagorinsky.cs"), "les.smagorinsky.cs");
                result.turbulence.les_wall_damping = parse_les_wall_damping(
                    require(entries, "les.smagorinsky.wall_damping"));
            } else if (has_key("les.smagorinsky.cs")
                       || has_key("les.smagorinsky.wall_damping")) {
                throw CaseConfigurationError(
                    "les.smagorinsky.* requires smagorinsky or mixed model");
            }
            if (similarity) {
                result.turbulence.les_similarity_coefficient = parse_real(
                    require(entries, "les.similarity.cb"), "les.similarity.cb");
            } else if (has_key("les.similarity.cb")) {
                throw CaseConfigurationError(
                    "les.similarity.cb requires scale_similarity or mixed model");
            }
            if (dynamic) {
                if (require(entries, "les.dynamic.average") != "local_box3_tensor") {
                    throw CaseConfigurationError(
                        "LES dynamic average must be local_box3_tensor");
                }
                if (require(entries, "les.dynamic.clipping") != "bounded") {
                    throw CaseConfigurationError("LES dynamic clipping must be bounded");
                }
                result.turbulence.les_dynamic_denominator_floor = parse_real(
                    require(entries, "les.dynamic.denominator_floor"),
                    "les.dynamic.denominator_floor");
                result.turbulence.les_dynamic_coefficient_minimum = parse_real(
                    require(entries, "les.dynamic.coefficient_minimum"),
                    "les.dynamic.coefficient_minimum");
                result.turbulence.les_dynamic_coefficient_maximum = parse_real(
                    require(entries, "les.dynamic.coefficient_maximum"),
                    "les.dynamic.coefficient_maximum");
            } else if (has_key("les.dynamic.average") || has_key("les.dynamic.clipping")
                       || has_key("les.dynamic.denominator_floor")
                       || has_key("les.dynamic.coefficient_minimum")
                       || has_key("les.dynamic.coefficient_maximum")) {
                throw CaseConfigurationError("les.dynamic.* requires dynamic_smagorinsky");
            }
            if (result.turbulence.kind == TurbulenceModelKind::Wale) {
                result.turbulence.les_wale_coefficient
                    = parse_real(require(entries, "les.wale.cw"), "les.wale.cw");
            } else if (has_key("les.wale.cw")) {
                throw CaseConfigurationError("les.wale.cw requires turbulence.model=wale");
            }
        }
        result.turbulence.turbulent_prandtl = optional_real(
            entries, "turbulence.prandtl", result.turbulence.turbulent_prandtl);
        if (has_wall_treatment) {
            result.turbulence.wall_treatment
                = parse_wall_treatment(require(entries, "turbulence.wall_treatment"));
        }
        result.turbulence.sa_farfield_nu_tilde_ratio = optional_real(
            entries,
            "turbulence.sa.farfield_nu_tilde_ratio",
            result.turbulence.sa_farfield_nu_tilde_ratio);
        if (has_sa_source) {
            result.turbulence.source_treatment = parse_turbulence_source_treatment(
                require(entries, "turbulence.sa.source_treatment"));
        }
        result.turbulence.freestream_turbulence_intensity = optional_real(
            entries,
            "turbulence.freestream.intensity",
            result.turbulence.freestream_turbulence_intensity);
        result.turbulence.freestream_length_scale = optional_real(
            entries,
            "turbulence.freestream.length_scale",
            result.turbulence.freestream_length_scale);
        result.turbulence.model_floor = optional_real(
            entries, "turbulence.model_floor", result.turbulence.model_floor);
        result.turbulence.wall_function_y_plus_min = optional_real(
            entries,
            "turbulence.wall_function.y_plus_min",
            result.turbulence.wall_function_y_plus_min);
        result.turbulence.wall_function_y_plus_max = optional_real(
            entries,
            "turbulence.wall_function.y_plus_max",
            result.turbulence.wall_function_y_plus_max);
        if (has_two_equation_source) {
            result.turbulence.source_treatment = parse_turbulence_source_treatment(
                require(entries, "turbulence.two_equation.source_treatment"));
        }
        result.time_algorithm.integrator
            = parse_time_integrator(require(entries, "time.integrator"));
        result.preconditioner.kind
            = parse_preconditioner(require(entries, "preconditioner.type"));
        const auto has = [&](const char* key) { return entries.find(key) != entries.end(); };
        const bool has_implicit_key = has("time.physical.scheme") || has("time.physical.step")
            || has("time.dual_time.max_iterations")
            || has("time.dual_time.absolute_tolerance")
            || has("time.dual_time.relative_tolerance") || has("time.dual_time.cfl")
            || has("lu_sgs.sweeps") || has("lu_sgs.jacobian")
            || has("lu_sgs.relaxation");
        if (result.time_algorithm.integrator == TimeIntegratorKind::SspRk3
            && has_implicit_key) {
            throw CaseConfigurationError(
                "time.integrator=ssprk3 does not accept LU-SGS or dual-time keys");
        }
        if (result.time_algorithm.integrator == TimeIntegratorKind::LuSgs) {
            if (has("time.physical.scheme")
                && require(entries, "time.physical.scheme") != "bdf2") {
                throw CaseConfigurationError("time.physical.scheme must be bdf2");
            }
            if (has("lu_sgs.jacobian")
                && require(entries, "lu_sgs.jacobian") != "scalar_spectral") {
                throw CaseConfigurationError("lu_sgs.jacobian must be scalar_spectral");
            }
            result.time_algorithm.physical_time_step
                = optional_real(entries, "time.physical.step", 0.0);
            result.time_algorithm.dual_time_max_iterations
                = optional_size(entries, "time.dual_time.max_iterations", 100);
            result.time_algorithm.dual_time_absolute_tolerance = optional_real(
                entries, "time.dual_time.absolute_tolerance", 1.0e-10);
            result.time_algorithm.dual_time_relative_tolerance = optional_real(
                entries, "time.dual_time.relative_tolerance", 1.0e-8);
            result.time_algorithm.dual_time_cfl
                = optional_real(entries, "time.dual_time.cfl", 5.0);
            if (has("lu_sgs.sweeps")) {
                const auto sweeps = parse_integer(require(entries, "lu_sgs.sweeps"),
                                                  "lu_sgs.sweeps");
                if (sweeps < std::numeric_limits<int>::min()
                    || sweeps > std::numeric_limits<int>::max()) {
                    throw CaseConfigurationError("lu_sgs.sweeps exceeds int range");
                }
                result.time_algorithm.lu_sgs_sweeps = static_cast<int>(sweeps);
            }
            result.time_algorithm.lu_sgs_relaxation
                = optional_real(entries, "lu_sgs.relaxation", 1.0);
        }
        const bool has_cutoff
            = has("preconditioner.mach_cutoff") || has("preconditioner.viscous_cutoff");
        if (result.preconditioner.kind == PreconditionerKind::None && has_cutoff) {
            throw CaseConfigurationError(
                "preconditioner.type=none does not accept cutoff keys");
        }
        if (result.preconditioner.kind == PreconditionerKind::WeissSmith) {
            result.preconditioner.mach_cutoff
                = optional_real(entries, "preconditioner.mach_cutoff", 1.0e-3);
            result.preconditioner.viscous_cutoff
                = optional_real(entries, "preconditioner.viscous_cutoff", 1.0);
        }
    }
    result.case_name = require(entries, "case.name");
    result.mesh_path = require(entries, "mesh.path");
    result.profile = ProfileFactory::from_string(require(entries, "algorithm.profile")).kind();
    if (const auto iterator = entries.find("algorithm.flux_difference");
        iterator != entries.end()) {
        result.flux_difference = parse_flux_difference_mode(iterator->second);
    }
    result.reconstruction.scheme = require(entries, "algorithm.reconstruction");
    result.reconstruction.variables
        = parse_reconstruction_variables(require(entries, "algorithm.reconstruction_variables"));
    result.reconstruction.nonlinear.mdcd_dispersion = optional_real(
        entries, "algorithm.mdcd.disp", result.reconstruction.nonlinear.mdcd_dispersion);
    result.reconstruction.nonlinear.mdcd_dissipation = optional_real(
        entries, "algorithm.mdcd.diss", result.reconstruction.nonlinear.mdcd_dissipation);
    result.riemann.scheme = require(entries, "algorithm.riemann");
    const bool has_all_speed_parameters
        = entries.find("algorithm.roe_all_speed.reference_mach") != entries.end()
        || entries.find("algorithm.roe_all_speed.dissipation_scale") != entries.end()
        || entries.find("algorithm.roe_all_speed.pressure_coefficient") != entries.end();
    if (has_all_speed_parameters && result.riemann.scheme != "roe_all_speed") {
        throw CaseConfigurationError(
            "algorithm.roe_all_speed parameters require algorithm.riemann=roe_all_speed");
    }
    if (result.riemann.scheme == "roe_all_speed") {
        result.riemann.parameters.all_speed.reference_mach = optional_real(
            entries,
            "algorithm.roe_all_speed.reference_mach",
            result.riemann.parameters.all_speed.reference_mach);
        result.riemann.parameters.all_speed.dissipation_scale = optional_real(
            entries,
            "algorithm.roe_all_speed.dissipation_scale",
            result.riemann.parameters.all_speed.dissipation_scale);
        result.riemann.parameters.all_speed.pressure_coefficient = optional_real(
            entries,
            "algorithm.roe_all_speed.pressure_coefficient",
            result.riemann.parameters.all_speed.pressure_coefficient);
    }
    result.robustness.enabled
        = optional_bool(entries, "robustness.enabled", result.robustness.enabled);
    if (const auto iterator = entries.find("robustness.max_local_recomputations");
        iterator != entries.end()) {
        const auto value = parse_integer(iterator->second, iterator->first);
        if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max()) {
            throw CaseConfigurationError("robustness max_local_recomputations exceeds int range");
        }
        result.robustness.max_local_recomputations = static_cast<int>(value);
    }
    if (const auto iterator = entries.find("robustness.max_step_retries");
        iterator != entries.end()) {
        const auto value = parse_integer(iterator->second, iterator->first);
        if (value < std::numeric_limits<int>::min() || value > std::numeric_limits<int>::max()) {
            throw CaseConfigurationError("robustness max_step_retries exceeds int range");
        }
        result.robustness.max_step_retries = static_cast<int>(value);
    }
    result.robustness.time_step_reduction = optional_real(
        entries, "robustness.time_step_reduction", result.robustness.time_step_reduction);
    result.robustness.minimum_time_step = optional_real(
        entries, "robustness.minimum_time_step", result.robustness.minimum_time_step);

    result.gas.gamma = parse_real(require(entries, "gas.gamma"), "gas.gamma");
    if (const auto iterator = entries.find("gas.molar_mass"); iterator != entries.end()) {
        result.gas.molar_mass = parse_real(iterator->second, iterator->first);
    }
    if (const auto iterator = entries.find("gas.specific_gas_constant");
        iterator != entries.end()) {
        result.gas.specific_gas_constant = parse_real(iterator->second, iterator->first);
    }

    result.reference.velocity
        = parse_real(require(entries, "reference.velocity"), "reference.velocity");
    result.reference.density
        = parse_real(require(entries, "reference.density"), "reference.density");
    result.reference.temperature
        = parse_real(require(entries, "reference.temperature"), "reference.temperature");
    result.reference.length = parse_real(require(entries, "reference.length"), "reference.length");
    result.reference.viscosity
        = parse_real(require(entries, "reference.viscosity"), "reference.viscosity");

    const auto transport_model = optional_string(entries, "transport.model", "constant");
    result.transport.prandtl
        = optional_real(entries, "transport.prandtl", result.transport.prandtl);
    const bool has_constant = entries.find("transport.constant.viscosity_ratio") != entries.end();
    const bool has_sutherland_mu
        = entries.find("transport.sutherland.reference_viscosity_ratio") != entries.end();
    const bool has_sutherland_temperature
        = entries.find("transport.sutherland.temperature") != entries.end();
    const bool has_sutherland_ratio
        = entries.find("transport.sutherland.temperature_ratio") != entries.end();
    if (transport_model == "constant") {
        if (has_sutherland_mu || has_sutherland_temperature || has_sutherland_ratio) {
            throw CaseConfigurationError(
                "Sutherland transport keys require transport.model=sutherland");
        }
        result.transport.viscosity
            = ConstantViscosity {optional_real(entries, "transport.constant.viscosity_ratio", 1.0)};
    } else if (transport_model == "sutherland") {
        if (has_constant) {
            throw CaseConfigurationError(
                "constant viscosity key requires transport.model=constant");
        }
        if (!has_sutherland_mu) {
            throw CaseConfigurationError("Sutherland transport requires reference_viscosity_ratio");
        }
        if (has_sutherland_temperature == has_sutherland_ratio) {
            throw CaseConfigurationError("Sutherland transport requires exactly one temperature or "
                                         "temperature_ratio key");
        }
        const Real ratio = has_sutherland_ratio
            ? parse_real(entries.at("transport.sutherland.temperature_ratio"),
                         "transport.sutherland.temperature_ratio")
            : parse_real(entries.at("transport.sutherland.temperature"),
                         "transport.sutherland.temperature")
                / result.reference.temperature;
        result.transport.viscosity = SutherlandViscosity {
            parse_real(entries.at("transport.sutherland.reference_viscosity_ratio"),
                       "transport.sutherland.reference_viscosity_ratio"),
            ratio,
        };
    } else {
        throw CaseConfigurationError("unknown transport model: " + transport_model);
    }

    result.partition.mode = parse_partition_mode(require(entries, "partition.mode"));
    result.partition.allow_idle_ranks
        = parse_bool(require(entries, "partition.allow_idle_ranks"), "partition.allow_idle_ranks");
    result.partition.max_load_ratio
        = parse_real(require(entries, "partition.max_load_ratio"), "partition.max_load_ratio");
    const auto min_cells
        = parse_integer(require(entries, "partition.min_cells_per_active_direction"),
                        "partition.min_cells_per_active_direction");
    if (min_cells < std::numeric_limits<int>::min()
        || min_cells > std::numeric_limits<int>::max()) {
        throw CaseConfigurationError("partition minimum cells exceed int range");
    }
    result.partition.min_cells_per_active_direction = static_cast<int>(min_cells);
    if (const auto iterator = entries.find("geometry.metric.fallback");
        iterator != entries.end()) {
        if (iterator->second == "strict") {
            result.metric_options.fallback = MetricFallback::Strict;
        } else if (iterator->second == "phenglei_finite_volume") {
            result.metric_options.fallback = MetricFallback::PhengleiFiniteVolume;
        } else {
            throw CaseConfigurationError("unknown geometry metric fallback: "
                                         + iterator->second);
        }
    }
    result.metric_options.maximum_reference_relative_difference = optional_real(
        entries,
        "geometry.metric.maximum_reference_relative_difference",
        result.metric_options.maximum_reference_relative_difference);

    result.initial.type = require(entries, "initial.type");
    constexpr std::string_view initial_prefix = "initial.";
    for (const auto& [key, value] : entries) {
        if (key.size() >= initial_prefix.size()
            && key.compare(0, initial_prefix.size(), initial_prefix.data(), initial_prefix.size())
                == 0
            && key != "initial.type") {
            result.initial.parameters.emplace(key.substr(initial_prefix.size()),
                                              parse_real(value, key));
        }
    }

    result.default_boundary = parse_boundary_type(require(entries, "boundary.default"));
    result.farfield_point_vortex.enabled
        = optional_bool(entries, "boundary.point_vortex.enabled", false);
    result.farfield_point_vortex.lift_coefficient = optional_real(
        entries, "boundary.point_vortex.lift_coefficient", 0.0);
    result.farfield_point_vortex.center = {{
        optional_real(entries, "boundary.point_vortex.center_x", 0.25),
        optional_real(entries, "boundary.point_vortex.center_y", 0.0),
    }};
    result.farfield_point_vortex.chord
        = optional_real(entries, "boundary.point_vortex.chord", 1.0);
    for (const auto& [key, value] : entries) {
        if (!dynamic_boundary_key(key)) continue;
        const auto [name, property] = split_boundary_key(key);
        if (property == "type") {
            result.boundary_overrides.emplace(name, parse_boundary_type(value));
            continue;
        }
        auto& data = result.boundary_data[name];
        const Real parsed = parse_real(value, key);
        if (property == "wall_velocity_x")
            data.wall_velocity[0] = parsed;
        else if (property == "wall_velocity_y")
            data.wall_velocity[1] = parsed;
        else if (property == "wall_velocity_z")
            data.wall_velocity[2] = parsed;
        else if (property == "wall_temperature")
            data.wall_temperature = parsed;
        else if (property == "rho")
            data.rho = parsed;
        else if (property == "u")
            data.u = parsed;
        else if (property == "v")
            data.v = parsed;
        else if (property == "w")
            data.w = parsed;
        else if (property == "temperature")
            data.temperature = parsed;
        else if (property == "pressure")
            data.pressure = parsed;
    }

    result.source_terms.enable_source_terms
        = parse_bool(require(entries, "source.enabled"), "source.enabled");
    if (const auto iterator = entries.find("source.models"); iterator != entries.end()) {
        for (const auto& name : split_list(iterator->second)) {
            result.source_terms.models.push_back(parse_source_model(name));
        }
    }
    const std::array<std::string, 5> uniform_keys {{
        "source.uniform.rho",
        "source.uniform.momentum_x",
        "source.uniform.momentum_y",
        "source.uniform.momentum_z",
        "source.uniform.energy",
    }};
    const std::array<std::string, 5> manufactured_keys {{
        "source.manufactured.rho",
        "source.manufactured.momentum_x",
        "source.manufactured.momentum_y",
        "source.manufactured.momentum_z",
        "source.manufactured.energy",
    }};
    for (std::size_t component = 0; component < 5; ++component) {
        result.source_terms.uniform_conservative[component]
            = optional_real(entries, uniform_keys[component], 0.0);
        result.source_terms.manufactured_amplitude[component]
            = optional_real(entries, manufactured_keys[component], 0.0);
    }
    result.source_terms.body_acceleration = {{
        optional_real(entries, "source.body.ax", 0.0),
        optional_real(entries, "source.body.ay", 0.0),
        optional_real(entries, "source.body.az", 0.0),
    }};
    result.source_terms.pressure_gradient = {{
        optional_real(entries, "source.pressure_gradient.x", 0.0),
        optional_real(entries, "source.pressure_gradient.y", 0.0),
        optional_real(entries, "source.pressure_gradient.z", 0.0),
    }};

    result.run.mode = parse_run_mode(require(entries, "run.mode"));
    result.run.viscous = parse_bool(require(entries, "run.viscous"), "run.viscous");
    result.run.cfl = parse_real(require(entries, "run.cfl"), "run.cfl");
    const auto max_steps = parse_integer(require(entries, "run.max_steps"), "run.max_steps");
    if (max_steps <= 0
        || static_cast<unsigned long long>(max_steps) > std::numeric_limits<std::size_t>::max()) {
        throw CaseConfigurationError("run max_steps is outside size_t range");
    }
    result.run.max_steps = static_cast<std::size_t>(max_steps);
    result.run.end_time = optional_real(entries, "run.t_end", 0.0);
    result.run.max_wall_time = optional_real(entries, "run.max_wall_time", 0.0);
    result.run.steady.min_steps = optional_size(entries, "steady.min_steps", 1);
    result.run.steady.check_interval_steps
        = optional_size(entries, "steady.check_interval_steps", 1);
    result.run.steady.consecutive_checks = optional_size(entries, "steady.consecutive_checks", 1);
    result.run.steady.reference_floor = optional_real(entries, "steady.reference_floor", 1.0e-30);
    result.run.steady.l2_absolute = optional_real(entries, "steady.l2_absolute", 1.0e-12);
    result.run.steady.l2_relative = optional_real(entries, "steady.l2_relative", 1.0e-8);
    result.run.steady.linf_enabled = optional_bool(entries, "steady.linf_enabled", true);
    result.run.steady.linf_absolute = optional_real(entries, "steady.linf_absolute", 1.0e-11);
    result.run.steady.linf_relative = optional_real(entries, "steady.linf_relative", 1.0e-8);

    result.output.directory = require(entries, "output.directory");
    result.output.allow_existing
        = parse_bool(require(entries, "output.allow_existing"), "output.allow_existing");
    result.output.dimensional
        = parse_bool(require(entries, "output.dimensional"), "output.dimensional");

    result.output.field.enabled
        = parse_bool(require(entries, "output.field.enabled"), "output.field.enabled");
    if (const auto iterator = entries.find("output.field.format"); iterator != entries.end()) {
        result.output.field.format = parse_field_output_format(iterator->second);
    }
    result.output.field.schedule = parse_schedule(entries, "output.field", true);
    result.output.field.quantities = optional_string_list(entries, "output.field.quantities");

    result.output.history.enabled
        = parse_bool(require(entries, "output.history.enabled"), "output.history.enabled");
    if (const auto iterator = entries.find("output.history.format"); iterator != entries.end()) {
        result.output.history.format = parse_series_output_format(iterator->second);
    }
    result.output.history.schedule = parse_schedule(entries, "output.history", true);
    result.output.history.quantities = optional_string_list(entries, "output.history.quantities");

    result.output.statistics.enabled
        = parse_bool(require(entries, "output.statistics.enabled"), "output.statistics.enabled");
    if (const auto iterator = entries.find("output.statistics.format"); iterator != entries.end()) {
        result.output.statistics.format = parse_series_output_format(iterator->second);
    }
    result.output.statistics.schedule = parse_schedule(entries, "output.statistics", true);
    result.output.statistics.quantities
        = optional_string_list(entries, "output.statistics.quantities");
    result.time_statistics.enabled
        = optional_bool(entries, "statistics.time.enabled", false);
    result.time_statistics.start_time
        = optional_real(entries, "statistics.time.start", 0.0);
    result.time_statistics.end_time
        = optional_real(entries, "statistics.time.end", result.run.end_time);
    result.time_statistics.every_steps
        = optional_size(entries, "statistics.time.every_steps", 1);
    if (const auto iterator = entries.find("statistics.time.weight");
        iterator != entries.end() && iterator->second != "accepted_dt") {
        throw CaseConfigurationError(
            "statistics.time.weight must be accepted_dt");
    }
    result.output.boundary.enabled = optional_bool(entries, "output.boundary.enabled", false);
    if (const auto iterator = entries.find("output.boundary.format"); iterator != entries.end()) {
        result.output.boundary.format = parse_series_output_format(iterator->second);
    }
    result.output.boundary.schedule = parse_schedule(entries, "output.boundary", true);
    result.output.boundary.patches = optional_string_list(entries, "output.boundary.patches");
    result.output.boundary.quantities = optional_string_list(entries, "output.boundary.quantities");
    result.output.boundary.reference_pressure = optional_real(
        entries, "output.boundary.reference_pressure", result.output.boundary.reference_pressure);
    result.output.boundary.reference_density = optional_real(
        entries, "output.boundary.reference_density", result.output.boundary.reference_density);
    result.output.boundary.reference_velocity = {{
        optional_real(entries,
                      "output.boundary.reference_velocity_x",
                      result.output.boundary.reference_velocity[0]),
        optional_real(entries,
                      "output.boundary.reference_velocity_y",
                      result.output.boundary.reference_velocity[1]),
        optional_real(entries,
                      "output.boundary.reference_velocity_z",
                      result.output.boundary.reference_velocity[2]),
    }};
    result.output.boundary.reference_area = optional_real(
        entries, "output.boundary.reference_area", result.output.boundary.reference_area);
    result.output.boundary.reference_length = optional_real(
        entries, "output.boundary.reference_length", result.output.boundary.reference_length);
    result.output.boundary.moment_center = {{
        optional_real(entries, "output.boundary.moment_center_x", 0.0),
        optional_real(entries, "output.boundary.moment_center_y", 0.0),
        optional_real(entries, "output.boundary.moment_center_z", 0.0),
    }};
    result.output.boundary.drag_direction = {{
        optional_real(
            entries, "output.boundary.drag_direction_x", result.output.boundary.drag_direction[0]),
        optional_real(
            entries, "output.boundary.drag_direction_y", result.output.boundary.drag_direction[1]),
        optional_real(
            entries, "output.boundary.drag_direction_z", result.output.boundary.drag_direction[2]),
    }};
    result.output.boundary.lift_direction = {{
        optional_real(
            entries, "output.boundary.lift_direction_x", result.output.boundary.lift_direction[0]),
        optional_real(
            entries, "output.boundary.lift_direction_y", result.output.boundary.lift_direction[1]),
        optional_real(
            entries, "output.boundary.lift_direction_z", result.output.boundary.lift_direction[2]),
    }};
    result.output.boundary.tangent_direction = {{
        optional_real(entries,
                      "output.boundary.tangent_direction_x",
                      result.output.boundary.tangent_direction[0]),
        optional_real(entries,
                      "output.boundary.tangent_direction_y",
                      result.output.boundary.tangent_direction[1]),
        optional_real(entries,
                      "output.boundary.tangent_direction_z",
                      result.output.boundary.tangent_direction[2]),
    }};
    result.output.boundary.span_direction = {{
        optional_real(entries,
                      "output.boundary.span_direction_x",
                      result.output.boundary.span_direction[0]),
        optional_real(entries,
                      "output.boundary.span_direction_y",
                      result.output.boundary.span_direction[1]),
        optional_real(entries,
                      "output.boundary.span_direction_z",
                      result.output.boundary.span_direction[2]),
    }};
    result.output.boundary.span_bin_edges
        = optional_real_list(entries, "output.boundary.span_bin_edges");
    result.output.xz_planes.enabled
        = optional_bool(entries, "output.statistics.xz_planes.enabled", false);
    result.output.xz_planes.cell_j_indices
        = optional_integer_list(entries, "output.statistics.xz_planes.cell_j_indices");
    if (result.output.xz_planes.enabled) {
        const auto names = xz_plane_statistic_names(result.output.xz_planes.cell_j_indices);
        result.output.statistics.quantities.insert(
            result.output.statistics.quantities.end(), names.begin(), names.end());
    }
    result.output.yz_planes.enabled
        = optional_bool(entries, "output.statistics.yz_planes.enabled", false);
    result.output.yz_planes.target_x_coordinates
        = optional_real_list(entries, "output.statistics.yz_planes.target_x_coordinates");
    if (result.output.yz_planes.enabled) {
        const auto names
            = yz_plane_statistic_names(result.output.yz_planes.target_x_coordinates.size());
        result.output.statistics.quantities.insert(
            result.output.statistics.quantities.end(), names.begin(), names.end());
    }
    result.output.channel_walls.enabled
        = optional_bool(entries, "output.statistics.channel_walls.enabled", false);
    result.output.channel_walls.lower_patch
        = optional_string(entries, "output.statistics.channel_walls.lower_patch", "bottom");
    result.output.channel_walls.upper_patch
        = optional_string(entries, "output.statistics.channel_walls.upper_patch", "top");
    result.output.channel_walls.half_height
        = optional_real(entries, "output.statistics.channel_walls.half_height", 1.0);
    if (result.output.channel_walls.enabled) {
        const auto names = channel_wall_statistic_names();
        result.output.statistics.quantities.insert(
            result.output.statistics.quantities.end(), names.begin(), names.end());
    }

    result.output.checkpoint.enabled
        = parse_bool(require(entries, "output.checkpoint.enabled"), "output.checkpoint.enabled");
    result.output.checkpoint.schedule = parse_schedule(entries, "output.checkpoint", true);
    if (const auto iterator = entries.find("restart.path"); iterator != entries.end()) {
        result.restart_path = iterator->second;
    }
    if (const auto iterator = entries.find("restart.mode"); iterator != entries.end()) {
        result.restart_mode = parse_restart_mode(iterator->second);
    }
    result.digest_ = fnv1a(canonical_entries(entries));
    result.hit.type = optional_string(entries, "hit.type", "none");
    const auto hit_n = optional_size(entries, "hit.n", 32);
    if(hit_n>512) throw CaseConfigurationError("hit.n exceeds 512");
    result.hit.n = static_cast<int>(hit_n);
    result.hit.seed = optional_size(entries, "hit.seed", 20261003);
    result.hit.length = optional_real(entries, "hit.length", 6.2831853071795864769);
    result.hit.spectrum_file = optional_string(entries, "hit.spectrum_file", "");
    result.hit.initial_energy = optional_real(entries, "hit.initial_energy", .705);
    result.hit.initialization = optional_string(entries, "hit.initialization", "shell_spectrum");
    result.hit.spectrum_amplitude = optional_real(entries, "hit.spectrum_amplitude", .00013);
    result.hit.preparation_time = optional_real(entries, "hit.preparation.time", 0);
    result.hit.preparation_max_steps = optional_size(entries, "hit.preparation.max_steps", 10000000);
    if(result.hit.initialization=="analytic_random_phase" && entries.count("hit.initial_energy"))
        throw CaseConfigurationError("analytic_random_phase derives energy from A and peak_wave; remove hit.initial_energy");
    result.hit.peak_wave = optional_real(entries, "hit.peak_wave", 4);
    result.hit.cutoff = optional_real(entries, "hit.cutoff", 0);
    result.hit.forcing = optional_string(entries, "hit.forcing", "constant_power");
    result.hit.forcing_power = optional_real(entries, "hit.forcing_power", .103);
    result.hit.forcing_kmax = optional_real(entries, "hit.forcing_kmax", 2);
    result.hit.remove_mean_acceleration = optional_bool(entries, "hit.remove_mean_acceleration", true);
    result.hit.thermostat = optional_bool(entries, "hit.thermostat", true);
    result.hit.statistics_start = optional_real(entries, "hit.statistics.start", 0);
    result.hit.statistics_end = optional_real(entries, "hit.statistics.end", 100);
    result.hit.sample_every_steps = optional_size(entries, "hit.sample_every_steps", 10);
    result.hit.write_every_samples = optional_size(entries, "hit.write_every_samples", 10);
    result.hit.sample_times = optional_real_list(entries, "hit.sample_times");
    if (!result.hit.enabled()) for (const auto& entry : entries)
        if (entry.first.rfind("hit.",0)==0 && entry.first!="hit.type")
            throw CaseConfigurationError("disabled HIT has configured parameters");

    result.chapter5.type = optional_string(entries, "benchmark.type", "none");
    result.chapter5.statistics_start = optional_real(entries, "benchmark.statistics.start", 0);
    result.chapter5.statistics_end = optional_real(entries, "benchmark.statistics.end", 100);
    result.chapter5.write_every_steps = optional_size(entries, "benchmark.write_every_steps", 1000);
    result.chapter5.history_every_steps = optional_size(entries, "benchmark.history_every_steps", 20);
    result.chapter5.probes = optional_real_list(entries, "benchmark.probes");
    result.source_terms.ramp_trip_amplitude = optional_real(entries, "source.ramp_trip.amplitude", 0);
    result.source_terms.ramp_trip_span = optional_real(entries, "source.ramp_trip.span", 6);
    result.periodic_hill.enabled = optional_bool(entries, "hill.enabled", false);
    result.periodic_hill.bulk_velocity = optional_real(entries, "hill.bulk_velocity", 1.0);
    result.periodic_hill.density = optional_real(entries, "hill.density", 1.0);
    result.periodic_hill.controller_time = optional_real(entries, "hill.controller_time", 1.0);
    result.periodic_hill.force_limit = optional_real(entries, "hill.force_limit", 0.5);
    result.periodic_hill.statistics_start = optional_real(entries, "hill.statistics.start", 100.0);
    result.periodic_hill.statistics_end = optional_real(entries, "hill.statistics.end", 1000.0);
    if (entries.count("hill.write_every_steps")) {
        const auto n = parse_integer(entries.at("hill.write_every_steps"), "hill.write_every_steps");
        if (n <= 0) throw CaseConfigurationError("hill.write_every_steps must be positive");
        result.periodic_hill.write_every_steps = static_cast<std::size_t>(n);
    }
    if (!result.periodic_hill.enabled) {
        for (const auto& entry : entries)
            if (entry.first.rfind("hill.", 0) == 0 && entry.first != "hill.enabled")
                throw CaseConfigurationError("hill settings require hill.enabled=true");
    }
    result.validate();
    return result;
}

CaseConfig CaseConfig::from_file(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) {
        throw CaseConfigurationError("cannot open case configuration: " + path);
    }
    std::ostringstream text;
    text << input.rdbuf();
    if (!input.good() && !input.eof()) {
        throw CaseConfigurationError("failed to read case configuration: " + path);
    }
    return from_text(text.str());
}

void CaseConfig::validate() const
{
    if (schema_version < minimum_schema_version || schema_version > supported_schema_version) {
        throw CaseConfigurationError("unsupported configuration schema version");
    }
    if (case_name.empty() || mesh_path.empty()) {
        throw CaseConfigurationError("case name and mesh path must not be empty");
    }
    if (restart_mode != RestartMode::Strict && restart_path.empty()) {
        throw CaseConfigurationError("restart.mode=algorithm_change requires restart.path");
    }
    const auto gas_model = make_gas_model();
    static_cast<void>(make_reference_scales(gas_model));
    transport.validate();
    try {
        turbulence.validate();
        time_algorithm.validate();
        preconditioner.validate();
    } catch (const std::invalid_argument& error) {
        throw CaseConfigurationError(error.what());
    }
    if (schema_version == 2) {
        const bool les = turbulence_model_family(turbulence.kind)
            == TurbulenceModelFamily::LesAlgebraic;
        if (turbulence.kind == TurbulenceModelKind::SaNegative) {
            if (!run.viscous) {
                throw CaseConfigurationError("SA-neg requires run.viscous=true");
            }
            if (turbulence.wall_treatment != WallTreatment::Resolved) {
                throw CaseConfigurationError("stage X SA-neg requires resolved wall treatment");
            }
            if (robustness.enabled) {
                throw CaseConfigurationError(
                    "stage X SA-neg does not yet support mean-flow step retry transactions");
            }
        }
        if (turbulence.kind == TurbulenceModelKind::KOmegaSst
            || turbulence.kind == TurbulenceModelKind::KEpsilon) {
            if (!run.viscous) {
                throw CaseConfigurationError("two-equation RANS requires run.viscous=true");
            }
            if (robustness.enabled) {
                throw CaseConfigurationError(
                    "two-equation RANS does not support mean-flow step retry transactions");
            }
        }
        if (les) {
            if (!run.viscous) {
                throw CaseConfigurationError("LES requires run.viscous=true");
            }
            if (run.mode != RunMode::Unsteady) {
                throw CaseConfigurationError("LES requires run.mode=unsteady");
            }
            if (time_algorithm.integrator != TimeIntegratorKind::LuSgs) {
                throw CaseConfigurationError("LES requires time.integrator=lu_sgs");
            }
            if (time_algorithm.physical_time_step <= 0.0) {
                throw CaseConfigurationError(
                    "LES requires BDF2 and a positive physical time step");
            }
            if (robustness.enabled) {
                throw CaseConfigurationError(
                    "LES does not support explicit robustness retries");
            }
            if (turbulence.les_wall_damping == LesWallDamping::VanDriest) {
                throw CaseConfigurationError(
                    "LES van-Driest damping requires a production wall-y+ field; "
                    "the current stage Z runtime cannot provide it");
            }
        }
        if (time_algorithm.integrator == TimeIntegratorKind::LuSgs) {
            if (robustness.enabled) {
                throw CaseConfigurationError(
                    "stage AA LU-SGS does not support explicit robustness retries");
            }
            if (run.mode == RunMode::Unsteady && time_algorithm.physical_time_step <= 0.0) {
                throw CaseConfigurationError(
                    "unsteady LU-SGS requires a positive time.physical.step");
            }
            if (run.mode == RunMode::Steady && time_algorithm.physical_time_step != 0.0) {
                throw CaseConfigurationError(
                    "steady LU-SGS does not accept time.physical.step");
            }
        }
        if (preconditioner.kind == PreconditionerKind::WeissSmith) {
            if (time_algorithm.integrator != TimeIntegratorKind::LuSgs) {
                throw CaseConfigurationError("Weiss-Smith requires time.integrator=lu_sgs");
            }
            if (riemann.scheme != "roe") {
                throw CaseConfigurationError("Weiss-Smith requires algorithm.riemann=roe");
            }
        }
    }
    hit.validate();
    if (hit.enabled()) {
        if(hit.preparation_time>0 && (time_statistics.enabled || restart_mode==RestartMode::AlgorithmChange))
            throw CaseConfigurationError("prepared HIT uses its own statistics and requires strict restart (no generic time statistics or algorithm_change)");
        if(schema_version!=2 || initial.type!="hit" || run.mode!=RunMode::Unsteady || !run.viscous
            || turbulence.kind!=TurbulenceModelKind::None || chapter5.enabled() || periodic_hill.enabled
            || time_algorithm.integrator!=TimeIntegratorKind::SspRk3 || preconditioner.kind!=PreconditionerKind::None
            || source_terms.enable_source_terms || output.dimensional
            || !std::holds_alternative<ConstantViscosity>(transport.viscosity))
            throw CaseConfigurationError("HIT requires schema 2, hit initial, viscous unsteady ILES/SSPRK3, constant viscosity, nondimensional output, no other source/benchmark");
    } else if(initial.type=="hit") throw CaseConfigurationError("hit initial requires hit.type=decay or forced");
    if (chapter5.enabled()) {
        if (chapter5.type != "sd7003" && chapter5.type != "compression_ramp")
            throw CaseConfigurationError("unknown benchmark.type");
        if (schema_version!=2 || !run.viscous || run.mode!=RunMode::Unsteady
            || turbulence.kind!=TurbulenceModelKind::None || periodic_hill.enabled
            || time_algorithm.integrator!=TimeIntegratorKind::SspRk3
            || preconditioner.kind!=PreconditionerKind::None || initial.type!=chapter5.type)
            throw CaseConfigurationError("chapter5 requires matching initial, schema 2, unsteady viscous ILES/SSPRK3, no preconditioner");
        if (!(chapter5.statistics_start>=0 && chapter5.statistics_end>chapter5.statistics_start)
            || !chapter5.write_every_steps || !chapter5.history_every_steps || chapter5.probes.size()%3)
            throw CaseConfigurationError("invalid benchmark statistics/probes");
        const auto ref=make_reference_scales(make_gas_model());
        const bool ramp=chapter5.type=="compression_ramp";
        if(std::abs(ref.mach()-(ramp?2.25:.2))>1e-10 || std::abs(ref.reynolds()-(ramp?15800:60000))>1e-6)
            throw CaseConfigurationError("chapter5 reference Ma/Re differs from benchmark");
        if(!ramp && std::abs(initial.parameter("alpha_degrees",4)-4)>1e-12)
            throw CaseConfigurationError("SD7003 chapter5 benchmark requires alpha=4 degrees");
        if (ramp && (std::abs(make_gas_model().gamma()-1.4)>1e-12
            || std::abs(transport.prandtl-.72)>1e-12
            || !std::holds_alternative<SutherlandViscosity>(transport.viscosity)
            || std::abs(std::get<SutherlandViscosity>(transport.viscosity).reference_viscosity_ratio-1)>1e-12
            || std::abs(std::get<SutherlandViscosity>(transport.viscosity).constant_temperature_ratio-110.4/170)>1e-12))
            throw CaseConfigurationError("fixed ramp similarity table requires gamma=1.4, Pr=.72, S/Tinf=110.4/170");
        if(ramp && std::abs(source_terms.ramp_trip_span-6)>1e-12)
            throw CaseConfigurationError("compression ramp trip requires span=6");
    }
    if (periodic_hill.enabled) {
        const auto& h = periodic_hill;
        if (schema_version != 2 || !run.viscous || run.mode != RunMode::Unsteady
            || turbulence.kind != TurbulenceModelKind::None
            || time_algorithm.integrator != TimeIntegratorKind::SspRk3
            || preconditioner.kind != PreconditionerKind::None || robustness.enabled)
            throw CaseConfigurationError("hill runtime requires schema 2 viscous unsteady ILES, SSPRK3, no preconditioner/retries");
        if (!(h.bulk_velocity > 0 && h.density > 0 && h.controller_time > 0
              && h.force_limit > 0 && h.statistics_start >= 0
              && h.statistics_end > h.statistics_start && h.write_every_steps > 0))
            throw CaseConfigurationError("invalid hill control/statistics parameters");
        if (!source_terms.enable_source_terms || source_terms.models.size() != 1
            || source_terms.models.front() != SourceModelKind::PressureGradient
            || source_terms.pressure_gradient[1] != 0 || source_terms.pressure_gradient[2] != 0
            || std::abs(source_terms.pressure_gradient[0]) > h.force_limit)
            throw CaseConfigurationError("hill runtime requires only an x pressure-gradient source within force_limit");
        if (!std::holds_alternative<ConstantViscosity>(transport.viscosity))
            throw CaseConfigurationError("hill wall statistics currently require constant viscosity");
    }
    static_cast<void>(make_profile());
    partition.validate(profile);
    metric_options.validate(make_profile());
    initial.validate();
    run.validate();
    output.validate(run.viscous);
    time_statistics.validate(run.mode, output.statistics.enabled);
    if (time_statistics.enabled && output.statistics.quantities.empty()) {
        throw CaseConfigurationError(
            "time statistics require at least one output.statistics quantity");
    }
    const std::set<std::string> modeled_wall_quantities {
        "wall_distance", "friction_velocity", "wall_y_plus", "wall_y_plus_class"};
    const bool rans_transport = turbulence.kind == TurbulenceModelKind::SaNegative
        || turbulence.kind == TurbulenceModelKind::KOmegaSst
        || turbulence.kind == TurbulenceModelKind::KEpsilon;
    const bool les = turbulence_model_family(turbulence.kind)
        == TurbulenceModelFamily::LesAlgebraic;
    if (!rans_transport && !les
        && std::any_of(output.boundary.quantities.begin(),
                       output.boundary.quantities.end(),
                       [&](const std::string& quantity) {
                           return modeled_wall_quantities.find(quantity)
                               != modeled_wall_quantities.end();
                       })) {
        throw CaseConfigurationError(
            "wall-unit boundary output requires an active RANS or LES model");
    }
    if (output.channel_walls.enabled && !run.viscous) {
        throw CaseConfigurationError("channel-wall friction monitoring requires run.viscous=true");
    }
    if (run.max_wall_time > 0.0 && !output.checkpoint.enabled) {
        throw CaseConfigurationError("positive max_wall_time requires checkpoint output");
    }
    if (default_boundary == BoundaryType::Undefined) {
        throw CaseConfigurationError("default boundary type must be defined");
    }
    farfield_point_vortex.validate(schema_version);
    if (farfield_point_vortex.enabled) {
        bool has_farfield = default_boundary == BoundaryType::Farfield;
        for (const auto& [name, type] : boundary_overrides) {
            static_cast<void>(name);
            has_farfield = has_farfield || type == BoundaryType::Farfield;
        }
        if (!has_farfield) {
            throw CaseConfigurationError(
                "enabled farfield point vortex requires a farfield boundary");
        }
    }
    for (const auto& [name, type] : boundary_overrides) {
        if (name.empty() || type == BoundaryType::Undefined) {
            throw CaseConfigurationError("boundary override is invalid");
        }
    }
    for (const auto& [name, data] : boundary_data) {
        if (name.empty()) {
            throw CaseConfigurationError("boundary data name must not be empty");
        }
        data.validate();
        const auto type_iterator = boundary_overrides.find(name);
        const BoundaryType type
            = type_iterator == boundary_overrides.end() ? default_boundary : type_iterator->second;
        const bool wall = type == BoundaryType::SlipWall
            || type == BoundaryType::NoSlipAdiabaticWall
            || type == BoundaryType::NoSlipIsothermalWall;
        const bool target = type == BoundaryType::Farfield || type == BoundaryType::Inflow
            || type == BoundaryType::Outflow;
        if (data.has_wall_velocity() && !wall) {
            throw CaseConfigurationError("wall velocity is only valid for wall boundary: " + name);
        }
        if (data.wall_temperature && type != BoundaryType::NoSlipIsothermalWall) {
            throw CaseConfigurationError("wall temperature requires an isothermal wall: " + name);
        }
        if (data.has_target_state() && !target) {
            throw CaseConfigurationError(
                "target state is only valid for inflow, farfield or outflow: " + name);
        }
    }
    const auto is_double_mach_boundary
        = [](BoundaryType type) { return type == BoundaryType::DoubleMachReflection; };
    bool has_double_mach_boundary = is_double_mach_boundary(default_boundary);
    for (const auto& [name, type] : boundary_overrides) {
        static_cast<void>(name);
        has_double_mach_boundary = has_double_mach_boundary || is_double_mach_boundary(type);
    }
    if (initial.type == "double_mach_reflection") {
        DoubleMachReflection(initial.parameter("x0", 1.0 / 6.0)).validate(gas_model.gamma(), 2);
        if (!has_double_mach_boundary) {
            throw CaseConfigurationError("double-Mach-reflection initial data require at least one "
                                         "double_mach_reflection boundary");
        }
        if (run.viscous || source_terms.enable_source_terms) {
            throw CaseConfigurationError(
                "classical double-Mach reflection must be inviscid and source-free");
        }
    } else if (has_double_mach_boundary) {
        throw CaseConfigurationError(
            "double_mach_reflection boundary requires matching initial.type");
    }
    auto inviscid = make_inviscid_config();
    inviscid.validate();
}

GasModel CaseConfig::make_gas_model() const
{
    return GasModel::from_input(gas);
}

ReferenceScales CaseConfig::make_reference_scales(const GasModel& gas_model) const
{
    return ReferenceScales::derive(reference, gas_model);
}

AlgorithmProfile CaseConfig::make_profile() const
{
    return ProfileFactory::create(profile);
}

InviscidWcnsConfig CaseConfig::make_inviscid_config() const
{
    InviscidWcnsConfig result;
    result.reconstruction = reconstruction;
    result.riemann = riemann;
    result.riemann.parameters.weiss_smith
        = preconditioner.kind == PreconditionerKind::WeissSmith;
    result.riemann.parameters.preconditioner
        = {preconditioner.mach_cutoff, preconditioner.viscous_cutoff};
    result.flux_difference = flux_difference;
    result.source_terms = source_terms;
    result.robustness = robustness;
    return result;
}

TransportConfig CaseConfig::make_transport_config() const
{
    transport.validate();
    return transport;
}

std::string CaseConfig::summary() const
{
    const auto gas_model = make_gas_model();
    const auto reference_scales = make_reference_scales(gas_model);
    std::vector<std::pair<std::string, BoundaryType>> boundaries(boundary_overrides.begin(),
                                                                 boundary_overrides.end());
    std::sort(boundaries.begin(), boundaries.end());
    std::vector<std::pair<std::string, BoundaryPhysicalDataConfig>> physical_data(
        boundary_data.begin(), boundary_data.end());
    std::sort(physical_data.begin(), physical_data.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });
    std::ostringstream result;
    result << "case(schema=" << schema_version << ",name=" << case_name << ",mesh=" << mesh_path
           << ",profile=" << make_profile().name()
           << ",flux_difference=" << flux_difference_mode_name(flux_difference) << ","
           << reconstruction.summary() << ',' << riemann.summary() << ',' << robustness.summary()
           << ',' << transport.summary() << ',' << gas_model.summary() << ','
           << reference_scales.summary() << ',' << partition.summary() << ',' << initial.summary()
           << ",metric(fallback="
           << (metric_options.fallback == MetricFallback::Strict ? "strict"
                                                                  : "phenglei_finite_volume")
           << ",max_reference_relative_difference=" << std::setprecision(17)
           << metric_options.maximum_reference_relative_difference << ')'
           << ",boundary.default=" << boundary_type_name(default_boundary) << ','
           << farfield_point_vortex.summary();
    for (const auto& [name, type] : boundaries) {
        result << ",boundary." << name << '=' << boundary_type_name(type);
    }
    for (const auto& [name, data] : physical_data) {
        result << ",boundary_data." << name << '(' << data.summary() << ')';
    }
    if (schema_version == 2) {
        result << ',' << turbulence.summary() << ',' << time_algorithm.summary() << ','
               << preconditioner.summary();
    }
    if (hit.enabled()) result << ',' << hit.signature();
    if (chapter5.enabled()) result << ',' << chapter5.signature();
    if (periodic_hill.enabled) result << ',' << periodic_hill.signature();
    result << ',' << source_terms.summary() << ',' << run.summary() << ',' << output.summary()
           << ',' << time_statistics.summary()
           << ",restart.path=" << (restart_path.empty() ? "<none>" : restart_path)
           << ",restart.mode=" << restart_mode_name(restart_mode) << ",digest=0x"
           << std::hex << digest_ << ')';
    return result.str();
}

std::string CaseConfig::legacy_v1_restart_signature() const
{
    std::vector<std::pair<std::string, BoundaryType>> boundaries(boundary_overrides.begin(),
                                                                 boundary_overrides.end());
    std::sort(boundaries.begin(), boundaries.end());
    std::vector<std::pair<std::string, BoundaryPhysicalDataConfig>> physical_data(
        boundary_data.begin(), boundary_data.end());
    std::sort(physical_data.begin(), physical_data.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.first < rhs.first;
    });
    std::ostringstream result;
    // This signature is the schema-1 compatibility identity even when the
    // active parser input is schema 2 with the exact none/SSPRK3/none mapping.
    result << "schema=1;profile=" << make_profile().restart_signature()
           << ";flux_difference=" << flux_difference_mode_name(flux_difference)
           << ";reconstruction=" << reconstruction.restart_signature()
           << ";riemann=" << riemann.restart_signature()
           << ";robustness=" << robustness.restart_signature()
           << ";gas=" << make_gas_model().restart_signature()
           << ";reference=" << make_reference_scales(make_gas_model()).restart_signature()
           << ";boundary.default=" << boundary_type_name(default_boundary);
    for (const auto& boundary : boundaries) {
        result << ";boundary." << boundary.first << '=' << boundary_type_name(boundary.second);
    }
    for (const auto& [name, data] : physical_data) {
        result << ";boundary_data." << name << '=' << data.summary();
    }
    result << ";source=" << source_terms.restart_signature()
           << ";viscous=" << (run.viscous ? "true" : "false");
    if (initial.type == "double_mach_reflection") {
        result << ';'
               << DoubleMachReflection(initial.parameter("x0", 1.0 / 6.0)).restart_signature();
    }
    return result.str();
}

std::string Chapter5Config::signature() const {
    std::ostringstream s; s<<std::setprecision(17)<<"chapter5_v1;type="<<type
        <<";start="<<statistics_start<<";end="<<statistics_end;
    for(Real p:probes) s<<';'<<p;
    return s.str();
}

std::string PeriodicHillConfig::signature() const
{
    std::ostringstream out;
    out << std::setprecision(17) << "hill_v1;Ub=" << bulk_velocity << ";rho=" << density
        << ";tau=" << controller_time << ";limit=" << force_limit
        << ";start=" << statistics_start << ";end=" << statistics_end;
    return out.str();
}

std::string CaseConfig::restart_signature() const
{
    std::string result
        = legacy_v1_restart_signature() + ";transport=" + transport.restart_signature();
    if (schema_version == 2) {
        result += ";" + turbulence.restart_signature() + ";"
            + time_algorithm.restart_signature() + ";" + preconditioner.restart_signature();
        if (farfield_point_vortex.enabled) {
            result += ";" + farfield_point_vortex.restart_signature();
        }
        const MetricBuildOptions defaults;
        if (metric_options.fallback != defaults.fallback
            || metric_options.maximum_reference_relative_difference
                != defaults.maximum_reference_relative_difference) {
            std::ostringstream metric_signature;
            metric_signature << std::setprecision(17) << ";metric_fallback="
                             << (metric_options.fallback == MetricFallback::Strict
                                     ? "strict"
                                     : "phenglei_finite_volume")
                             << ";metric_reference_tolerance="
                             << metric_options.maximum_reference_relative_difference;
            result += metric_signature.str();
        }
        if (time_statistics.enabled) {
            result += ";" + time_statistics_identity();
        }
    }
    if (hit.enabled()) result += ";" + hit.signature();
    if (chapter5.enabled()) result += ";" + chapter5.signature();
    if (periodic_hill.enabled) result += ";" + periodic_hill.signature();
    return result;
}

std::string CaseConfig::time_statistics_identity() const
{
    return time_statistics.restart_signature(output.statistics.quantities);
}

} // namespace wcns
