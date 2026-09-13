#pragma once

#include <wcns/mesh/algorithm_profile.hpp>
#include <wcns/physics/source_terms.hpp>
#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/inviscid_wcns_solver.hpp>
#include <wcns/solver/transport_model.hpp>
#include <wcns/solver/turbulence_model.hpp>

#include <cstddef>
#include <cstdint>
#include <array>
#include <optional>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace wcns {

class CaseConfigurationError : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

enum class PartitionMode {
    ZonesOnly,
    AutoSplit,
    ForceSplit,
};

struct PartitionConfig {
    PartitionMode mode = PartitionMode::AutoSplit;
    bool allow_idle_ranks = false;
    Real max_load_ratio = 1.2;
    int min_cells_per_active_direction = 8;

    void validate(AlgorithmProfileKind profile) const;
    [[nodiscard]] std::string summary() const;
};

struct InitialConditionConfig {
    std::string type = "uniform";
    std::unordered_map<std::string, Real> parameters;

    void validate(int dimension = 3) const;
    [[nodiscard]] Real parameter(const std::string& name, Real default_value) const;
    [[nodiscard]] std::string summary() const;
};

struct BoundaryPhysicalDataConfig {
    std::array<std::optional<Real>, 3> wall_velocity;
    std::optional<Real> wall_temperature;
    std::optional<Real> rho;
    std::optional<Real> u;
    std::optional<Real> v;
    std::optional<Real> w;
    std::optional<Real> temperature;
    std::optional<Real> pressure;

    [[nodiscard]] bool has_wall_velocity() const noexcept;
    [[nodiscard]] bool has_target_state() const noexcept;
    void validate() const;
    [[nodiscard]] std::string summary() const;
};

enum class RunMode {
    Steady,
    Unsteady,
};

enum class TimeIntegratorKind {
    SspRk3,
    LuSgs,
};

enum class PreconditionerKind {
    None,
    WeissSmith,
};

struct TimeAlgorithmConfig {
    TimeIntegratorKind integrator = TimeIntegratorKind::SspRk3;
    Real physical_time_step = 0.0;
    std::size_t dual_time_max_iterations = 100;
    Real dual_time_absolute_tolerance = 1.0e-10;
    Real dual_time_relative_tolerance = 1.0e-8;
    Real dual_time_cfl = 5.0;
    int lu_sgs_sweeps = 1;
    Real lu_sgs_relaxation = 1.0;

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

struct PreconditionerConfig {
    PreconditionerKind kind = PreconditionerKind::None;
    Real mach_cutoff = 1.0e-3;
    Real viscous_cutoff = 1.0;

    void validate() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

struct SteadyConvergenceConfig {
    std::size_t min_steps = 1;
    std::size_t check_interval_steps = 1;
    std::size_t consecutive_checks = 1;
    Real reference_floor = 1.0e-30;
    Real l2_absolute = 1.0e-12;
    Real l2_relative = 1.0e-8;
    bool linf_enabled = true;
    Real linf_absolute = 1.0e-11;
    Real linf_relative = 1.0e-8;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

struct OutputScheduleConfig {
    std::size_t every_steps = 0;
    Real every_time = 0.0;
    std::vector<Real> explicit_times;
    bool write_initial = false;
    bool write_final = true;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

enum class FieldOutputFormat {
    Cgns,
    Tecplot,
    Both,
};

enum class SeriesOutputFormat {
    Text,
    Tecplot,
};

struct FieldOutputConfig {
    bool enabled = false;
    FieldOutputFormat format = FieldOutputFormat::Cgns;
    OutputScheduleConfig schedule;
    std::vector<std::string> quantities;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

struct SeriesOutputConfig {
    bool enabled = false;
    SeriesOutputFormat format = SeriesOutputFormat::Text;
    OutputScheduleConfig schedule;
    std::vector<std::string> quantities;

    void validate(const char* label) const;
    [[nodiscard]] std::string summary(const char* label) const;
};

struct BoundaryOutputConfig {
    bool enabled = false;
    SeriesOutputFormat format = SeriesOutputFormat::Text;
    OutputScheduleConfig schedule;
    std::vector<std::string> patches;
    std::vector<std::string> quantities;
    Real reference_pressure = std::numeric_limits<Real>::quiet_NaN();
    Real reference_density = std::numeric_limits<Real>::quiet_NaN();
    std::array<Real, 3> reference_velocity {{std::numeric_limits<Real>::quiet_NaN(),
                                             std::numeric_limits<Real>::quiet_NaN(),
                                             std::numeric_limits<Real>::quiet_NaN()}};
    Real reference_area = std::numeric_limits<Real>::quiet_NaN();
    Real reference_length = std::numeric_limits<Real>::quiet_NaN();
    std::array<Real, 3> moment_center {{0.0, 0.0, 0.0}};
    std::array<Real, 3> drag_direction {{std::numeric_limits<Real>::quiet_NaN(),
                                         std::numeric_limits<Real>::quiet_NaN(),
                                         std::numeric_limits<Real>::quiet_NaN()}};
    std::array<Real, 3> lift_direction {{std::numeric_limits<Real>::quiet_NaN(),
                                         std::numeric_limits<Real>::quiet_NaN(),
                                         std::numeric_limits<Real>::quiet_NaN()}};
    std::array<Real, 3> tangent_direction {{std::numeric_limits<Real>::quiet_NaN(),
                                            std::numeric_limits<Real>::quiet_NaN(),
                                            std::numeric_limits<Real>::quiet_NaN()}};

    void validate(bool viscous) const;
    [[nodiscard]] std::string summary() const;
};

struct XzPlaneStatisticsConfig {
    bool enabled = false;
    std::vector<int> cell_j_indices;

    void validate(bool statistics_enabled) const;
    [[nodiscard]] std::string summary() const;
};

struct YzPlaneStatisticsConfig {
    bool enabled = false;
    // Each target selects the nearest geometrically planar cell-centre section
    // on its positive-x side (an exactly coincident section is selected).
    std::vector<Real> target_x_coordinates;

    void validate(bool statistics_enabled) const;
    [[nodiscard]] std::string summary() const;
};

struct ChannelWallStatisticsConfig {
    bool enabled = false;
    std::string lower_patch = "bottom";
    std::string upper_patch = "top";
    Real half_height = 1.0;

    void validate(bool statistics_enabled) const;
    [[nodiscard]] std::string summary() const;
};

struct CheckpointOutputConfig {
    bool enabled = false;
    OutputScheduleConfig schedule;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

struct OutputConfig {
    std::string directory = "output";
    bool allow_existing = false;
    bool dimensional = false;
    FieldOutputConfig field;
    SeriesOutputConfig history;
    SeriesOutputConfig statistics;
    BoundaryOutputConfig boundary;
    XzPlaneStatisticsConfig xz_planes;
    YzPlaneStatisticsConfig yz_planes;
    ChannelWallStatisticsConfig channel_walls;
    CheckpointOutputConfig checkpoint;

    void validate(bool viscous) const;
    [[nodiscard]] std::string summary() const;
};

struct CaseRunConfig {
    RunMode mode = RunMode::Steady;
    bool viscous = false;
    Real cfl = 0.2;
    std::size_t max_steps = 1;
    Real end_time = 0.0;
    Real max_wall_time = 0.0;
    SteadyConvergenceConfig steady;

    void validate() const;
    [[nodiscard]] std::string summary() const;
};

struct CaseConfig {
    static constexpr int minimum_schema_version = 1;
    static constexpr int supported_schema_version = 2;

    int schema_version = minimum_schema_version;
    std::string case_name;
    std::string mesh_path;
    AlgorithmProfileKind profile = AlgorithmProfileKind::PhengleiWcns;
    FluxDifferenceMode flux_difference = FluxDifferenceMode::Profile;
    ReconstructionConfig reconstruction {};
    RiemannConfig riemann {};
    RobustnessConfig robustness {};
    TransportConfig transport {};
    TurbulenceModelConfig turbulence {};
    TimeAlgorithmConfig time_algorithm {};
    PreconditionerConfig preconditioner {};
    GasModelInput gas;
    ReferenceInput reference;
    PartitionConfig partition;
    InitialConditionConfig initial;
    BoundaryType default_boundary = BoundaryType::Farfield;
    std::unordered_map<std::string, BoundaryType> boundary_overrides;
    std::unordered_map<std::string, BoundaryPhysicalDataConfig> boundary_data;
    SourceTermConfig source_terms;
    CaseRunConfig run;
    OutputConfig output;
    std::string restart_path;

    [[nodiscard]] static CaseConfig from_text(const std::string& text);
    [[nodiscard]] static CaseConfig from_file(const std::string& path);

    void validate() const;
    [[nodiscard]] GasModel make_gas_model() const;
    [[nodiscard]] ReferenceScales make_reference_scales(const GasModel& gas_model) const;
    [[nodiscard]] AlgorithmProfile make_profile() const;
    [[nodiscard]] InviscidWcnsConfig make_inviscid_config() const;
    [[nodiscard]] TransportConfig make_transport_config() const;
    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
    [[nodiscard]] std::string legacy_v1_restart_signature() const;
    [[nodiscard]] std::uint64_t digest() const noexcept { return digest_; }

private:
    std::uint64_t digest_ = 0;
};

[[nodiscard]] const char* partition_mode_name(PartitionMode mode);
[[nodiscard]] const char* boundary_type_name(BoundaryType type);
[[nodiscard]] const char* run_mode_name(RunMode mode);
[[nodiscard]] const char* time_integrator_name(TimeIntegratorKind integrator);
[[nodiscard]] const char* preconditioner_name(PreconditionerKind preconditioner);
[[nodiscard]] const char* field_output_format_name(FieldOutputFormat format);
[[nodiscard]] const char* series_output_format_name(SeriesOutputFormat format);

} // namespace wcns
