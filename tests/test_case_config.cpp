#include "test_support.hpp"

#include <wcns/runtime/case_config.hpp>
#include <wcns/runtime/quantity_registry.hpp>

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>

namespace {

std::string valid_config()
{
    return R"(
schema_version = 1
case.name = parser-test
mesh.path = mesh.cgns
algorithm.profile = phenglei_wcns
algorithm.reconstruction = weno_z
algorithm.reconstruction_variables = characteristic
algorithm.riemann = hllc
algorithm.mdcd.disp = 0.04
algorithm.mdcd.diss = 0.005
robustness.enabled = true
robustness.max_local_recomputations = 2
robustness.max_step_retries = 5
robustness.time_step_reduction = 0.4
robustness.minimum_time_step = 1e-10
gas.gamma = 1.4
gas.molar_mass = 0.029
reference.velocity = 340
reference.density = 1.225
reference.temperature = 288.15
reference.length = 1
reference.viscosity = 1.7894e-5
partition.mode = auto_split
partition.allow_idle_ranks = false
partition.max_load_ratio = 1.2
partition.min_cells_per_active_direction = 8
initial.type = uniform
initial.rho = 1
initial.u = 0.2
initial.v = 0
initial.w = 0
initial.temperature = 1
boundary.default = farfield
boundary.wall.type = no_slip_isothermal_wall
boundary.wall.wall_velocity_x = 0.75
boundary.wall.wall_temperature = 1.1
boundary.inlet.type = inflow
boundary.inlet.rho = 1.05
boundary.inlet.u = 0.3
boundary.inlet.temperature = 0.95
source.enabled = false
run.mode = steady
run.viscous = false
run.cfl = 0.2
run.max_steps = 10
output.directory = output
output.allow_existing = false
output.dimensional = false
output.field.enabled = false
output.history.enabled = false
output.statistics.enabled = false
output.checkpoint.enabled = false
)";
}

std::string valid_v2_config()
{
    auto result = valid_config();
    const auto version = result.find("schema_version = 1");
    result.replace(version, std::string("schema_version = 1").size(), "schema_version = 2");
    const auto line_end = result.find('\n', version);
    result.insert(line_end + 1,
                  "turbulence.model = none\n"
                  "time.integrator = ssprk3\n"
                  "preconditioner.type = none\n");
    return result;
}

std::string valid_sa_config()
{
    auto result = valid_v2_config();
    auto replace = [&](const std::string& from, const std::string& to) {
        const auto position = result.find(from);
        if (position == std::string::npos) throw std::logic_error("SA test fixture is incomplete");
        result.replace(position, from.size(), to);
    };
    replace("turbulence.model = none", "turbulence.model = sa_neg");
    replace("robustness.enabled = true", "robustness.enabled = false");
    replace("run.viscous = false", "run.viscous = true");
    const auto model = result.find("turbulence.model = sa_neg");
    const auto line_end = result.find('\n', model);
    result.insert(line_end + 1,
                  "turbulence.prandtl = 0.9\n"
                  "turbulence.wall_treatment = resolved\n"
                  "turbulence.sa.farfield_nu_tilde_ratio = 3\n"
                  "turbulence.sa.source_treatment = local_implicit\n");
    return result;
}

std::string valid_two_equation_config(const std::string& model,
                                      const std::string& wall_treatment)
{
    auto result = valid_v2_config();
    auto replace = [&](const std::string& from, const std::string& to) {
        const auto position = result.find(from);
        if (position == std::string::npos) {
            throw std::logic_error("two-equation test fixture is incomplete");
        }
        result.replace(position, from.size(), to);
    };
    replace("turbulence.model = none", "turbulence.model = " + model);
    replace("robustness.enabled = true", "robustness.enabled = false");
    replace("run.viscous = false", "run.viscous = true");
    const auto model_key = result.find("turbulence.model = " + model);
    const auto line_end = result.find('\n', model_key);
    result.insert(line_end + 1,
                  (model == "k_epsilon" ? std::string("turbulence.experimental = true\n")
                                        : std::string())
                      + "turbulence.prandtl = 0.9\n"
                  "turbulence.wall_treatment = "
                      + wall_treatment
                      + "\n"
                        "turbulence.freestream.intensity = 0.02\n"
                        "turbulence.freestream.length_scale = 0.15\n"
                        "turbulence.model_floor = 1e-11\n"
                        "turbulence.two_equation.source_treatment = local_implicit\n"
                      + (wall_treatment == "wall_function"
                             ? "turbulence.wall_function.y_plus_min = 30\n"
                               "turbulence.wall_function.y_plus_max = 300\n"
                             : ""));
    return result;
}

std::string valid_lu_sgs_config(bool unsteady = false)
{
    auto result = valid_v2_config();
    auto replace = [&](const std::string& from, const std::string& to) {
        const auto position = result.find(from);
        if (position == std::string::npos) throw std::logic_error("LU-SGS test fixture is incomplete");
        result.replace(position, from.size(), to);
    };
    replace("time.integrator = ssprk3", "time.integrator = lu_sgs");
    replace("robustness.enabled = true", "robustness.enabled = false");
    if (unsteady) {
        replace("run.mode = steady", "run.mode = unsteady\nrun.t_end = 0.25");
    }
    const auto integrator = result.find("time.integrator = lu_sgs");
    const auto line_end = result.find('\n', integrator);
    result.insert(line_end + 1,
                  "time.physical.scheme = bdf2\n"
                  + std::string(unsteady ? "time.physical.step = 0.01\n" : "")
                  + "time.dual_time.max_iterations = 40\n"
                    "time.dual_time.absolute_tolerance = 1e-11\n"
                    "time.dual_time.relative_tolerance = 1e-7\n"
                    "time.dual_time.cfl = 4\n"
                    "lu_sgs.sweeps = 1\n"
                    "lu_sgs.jacobian = scalar_spectral\n"
                    "lu_sgs.relaxation = 0.9\n");
    return result;
}

} // namespace

void test_case_config()
{
    {
        const auto config = wcns::CaseConfig::from_text(valid_v2_config());
        WCNS_REQUIRE(config.schema_version == 2);
        WCNS_REQUIRE(config.turbulence.kind == wcns::TurbulenceModelKind::None);
        WCNS_REQUIRE(config.time_algorithm.integrator == wcns::TimeIntegratorKind::SspRk3);
        WCNS_REQUIRE(config.preconditioner.kind == wcns::PreconditionerKind::None);
        WCNS_REQUIRE(config.summary().find("turbulence(model=none") != std::string::npos);
        WCNS_REQUIRE(config.restart_signature().find("turbulence_v1") != std::string::npos);
        const auto schema_one = wcns::CaseConfig::from_text(valid_config());
        WCNS_REQUIRE(config.legacy_v1_restart_signature()
                         + ";transport=" + config.transport.restart_signature()
                     == schema_one.restart_signature());

        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(
                                valid_config() + "turbulence.model = none\n"));
        auto missing = valid_v2_config();
        const auto preconditioner = missing.find("preconditioner.type = none\n");
        missing.erase(preconditioner, std::string("preconditioner.type = none\n").size());
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError, wcns::CaseConfig::from_text(missing));

        auto model_key = valid_v2_config();
        const auto model = model_key.find("turbulence.model = none");
        model_key.replace(model,
                          std::string("turbulence.model = none").size(),
                          "turbulence.model = sa_neg");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(model_key));

        const auto sa = wcns::CaseConfig::from_text(valid_sa_config());
        WCNS_REQUIRE(sa.turbulence.kind == wcns::TurbulenceModelKind::SaNegative);
        WCNS_REQUIRE(sa.turbulence.wall_treatment == wcns::WallTreatment::Resolved);
        WCNS_REQUIRE(sa.turbulence.source_treatment
                     == wcns::TurbulenceSourceTreatment::LocalImplicit);
        WCNS_REQUIRE_NEAR(sa.turbulence.sa_farfield_nu_tilde_ratio, 3.0, 0.0);
        WCNS_REQUIRE(sa.summary().find("model=sa_neg") != std::string::npos);
        WCNS_REQUIRE(sa.restart_signature().find("source=local_implicit")
                     != std::string::npos);

        auto invalid_ratio = valid_sa_config();
        const auto ratio = invalid_ratio.find("farfield_nu_tilde_ratio = 3");
        invalid_ratio.replace(ratio,
                              std::string("farfield_nu_tilde_ratio = 3").size(),
                              "farfield_nu_tilde_ratio = 2.99");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(invalid_ratio));

        auto invalid_source = valid_sa_config();
        const auto source = invalid_source.find("source_treatment = local_implicit");
        invalid_source.replace(source,
                               std::string("source_treatment = local_implicit").size(),
                               "source_treatment = hidden_clipping");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(invalid_source));

        const auto sst
            = wcns::CaseConfig::from_text(valid_two_equation_config("k_omega_sst", "resolved"));
        WCNS_REQUIRE(sst.turbulence.kind == wcns::TurbulenceModelKind::KOmegaSst);
        WCNS_REQUIRE(sst.turbulence.wall_treatment == wcns::WallTreatment::Resolved);
        WCNS_REQUIRE_NEAR(sst.turbulence.freestream_turbulence_intensity, 0.02, 0.0);
        WCNS_REQUIRE_NEAR(sst.turbulence.freestream_length_scale, 0.15, 0.0);
        WCNS_REQUIRE_NEAR(sst.turbulence.model_floor, 1.0e-11, 0.0);
        WCNS_REQUIRE(sst.restart_signature().find("model=k_omega_sst")
                     != std::string::npos);

        const auto k_epsilon = wcns::CaseConfig::from_text(
            valid_two_equation_config("k_epsilon", "wall_function"));
        WCNS_REQUIRE(k_epsilon.turbulence.kind == wcns::TurbulenceModelKind::KEpsilon);
        WCNS_REQUIRE(k_epsilon.turbulence.wall_treatment
                     == wcns::WallTreatment::WallFunction);
        WCNS_REQUIRE(k_epsilon.turbulence.experimental);
        WCNS_REQUIRE_NEAR(k_epsilon.turbulence.wall_function_y_plus_min, 30.0, 0.0);
        WCNS_REQUIRE_NEAR(k_epsilon.turbulence.wall_function_y_plus_max, 300.0, 0.0);

        auto invalid_k_epsilon = valid_two_equation_config("k_epsilon", "wall_function");
        const auto wall = invalid_k_epsilon.find("wall_function");
        invalid_k_epsilon.replace(wall, std::string("wall_function").size(), "resolved");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(invalid_k_epsilon));

        auto missing_experimental
            = valid_two_equation_config("k_epsilon", "wall_function");
        const auto experimental_line
            = missing_experimental.find("turbulence.experimental = true\n");
        missing_experimental.erase(experimental_line,
                                   std::string("turbulence.experimental = true\n").size());
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(missing_experimental));

        auto invalid_intensity = valid_two_equation_config("k_omega_sst", "resolved");
        const auto intensity = invalid_intensity.find("intensity = 0.02");
        invalid_intensity.replace(intensity,
                                  std::string("intensity = 0.02").size(),
                                  "intensity = 0");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(invalid_intensity));

        auto none_parameter = valid_v2_config();
        none_parameter += "turbulence.prandtl = 0.9\n";
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(none_parameter));

        const auto implicit = wcns::CaseConfig::from_text(valid_lu_sgs_config());
        WCNS_REQUIRE(implicit.time_algorithm.integrator == wcns::TimeIntegratorKind::LuSgs);
        WCNS_REQUIRE(implicit.time_algorithm.dual_time_max_iterations == 40);
        WCNS_REQUIRE_NEAR(implicit.time_algorithm.dual_time_cfl, 4.0, 0.0);
        WCNS_REQUIRE_NEAR(implicit.time_algorithm.lu_sgs_relaxation, 0.9, 0.0);
        WCNS_REQUIRE(implicit.summary().find("jacobian=scalar_spectral")
                     != std::string::npos);

        auto three_sweeps_text = valid_lu_sgs_config();
        const auto sweep_value = three_sweeps_text.find("lu_sgs.sweeps = 1");
        three_sweeps_text.replace(sweep_value,
                                  std::string("lu_sgs.sweeps = 1").size(),
                                  "lu_sgs.sweeps = 3");
        const auto three_sweeps = wcns::CaseConfig::from_text(three_sweeps_text);
        WCNS_REQUIRE(three_sweeps.time_algorithm.lu_sgs_sweeps == 3);

        auto invalid_sweeps = valid_lu_sgs_config();
        const auto invalid_sweep_value = invalid_sweeps.find("lu_sgs.sweeps = 1");
        invalid_sweeps.replace(invalid_sweep_value,
                               std::string("lu_sgs.sweeps = 1").size(),
                               "lu_sgs.sweeps = 5");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(invalid_sweeps));

        const auto dual = wcns::CaseConfig::from_text(valid_lu_sgs_config(true));
        WCNS_REQUIRE(dual.run.mode == wcns::RunMode::Unsteady);
        WCNS_REQUIRE_NEAR(dual.time_algorithm.physical_time_step, 0.01, 0.0);

        auto explicit_with_lu_key = valid_v2_config();
        explicit_with_lu_key += "lu_sgs.sweeps = 1\n";
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(explicit_with_lu_key));

        auto missing_physical_step = valid_lu_sgs_config(true);
        const auto step = missing_physical_step.find("time.physical.step = 0.01\n");
        missing_physical_step.erase(step, std::string("time.physical.step = 0.01\n").size());
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(missing_physical_step));

        auto preconditioned = valid_lu_sgs_config();
        const auto preconditioner_key = preconditioned.find("preconditioner.type = none");
        preconditioned.replace(preconditioner_key,
                               std::string("preconditioner.type = none").size(),
                               "preconditioner.type = weiss_smith\n"
                               "preconditioner.mach_cutoff = 0.002\n"
                               "preconditioner.viscous_cutoff = 1.5");
        const auto riemann = preconditioned.find("algorithm.riemann = hllc");
        preconditioned.replace(
            riemann, std::string("algorithm.riemann = hllc").size(), "algorithm.riemann = roe");
        const auto low_mach = wcns::CaseConfig::from_text(preconditioned);
        WCNS_REQUIRE(low_mach.preconditioner.kind == wcns::PreconditionerKind::WeissSmith);
        WCNS_REQUIRE_NEAR(low_mach.preconditioner.mach_cutoff, 0.002, 0.0);
    }
    {
        const auto config = wcns::CaseConfig::from_text(valid_config());
        WCNS_REQUIRE(config.schema_version == 1);
        WCNS_REQUIRE(config.case_name == "parser-test");
        WCNS_REQUIRE(config.mesh_path == "mesh.cgns");
        WCNS_REQUIRE(config.profile == wcns::AlgorithmProfileKind::PhengleiWcns);
        WCNS_REQUIRE(config.flux_difference == wcns::FluxDifferenceMode::Profile);
        WCNS_REQUIRE(config.reconstruction.scheme == "weno_z");
        WCNS_REQUIRE(config.reconstruction.variables
                     == wcns::ReconstructionVariables::Characteristic);
        WCNS_REQUIRE(config.riemann.scheme == "hllc");
        WCNS_REQUIRE(config.robustness.enabled);
        WCNS_REQUIRE(config.robustness.max_local_recomputations == 2);
        WCNS_REQUIRE(config.robustness.max_step_retries == 5);
        WCNS_REQUIRE_NEAR(config.robustness.time_step_reduction, 0.4, 1.0e-15);
        WCNS_REQUIRE_NEAR(config.robustness.minimum_time_step, 1.0e-10, 1.0e-20);
        WCNS_REQUIRE_NEAR(config.reconstruction.nonlinear.mdcd_dispersion, 0.04, 1.0e-15);
        WCNS_REQUIRE_NEAR(config.reconstruction.nonlinear.mdcd_dissipation, 0.005, 1.0e-15);
        WCNS_REQUIRE(config.restart_signature().find("mdcd_dispersion=0.04") != std::string::npos);
        WCNS_REQUIRE(config.restart_signature().find("robustness_v1") != std::string::npos);
        WCNS_REQUIRE(config.partition.mode == wcns::PartitionMode::AutoSplit);
        WCNS_REQUIRE(config.boundary_overrides.at("wall")
                     == wcns::BoundaryType::NoSlipIsothermalWall);
        WCNS_REQUIRE_NEAR(*config.boundary_data.at("wall").wall_velocity[0], 0.75, 1.0e-15);
        WCNS_REQUIRE_NEAR(*config.boundary_data.at("wall").wall_temperature, 1.1, 1.0e-15);
        WCNS_REQUIRE_NEAR(*config.boundary_data.at("inlet").rho, 1.05, 1.0e-15);
        WCNS_REQUIRE(config.restart_signature().find("wall_velocity_x=0.75") != std::string::npos);
        WCNS_REQUIRE(config.run.max_steps == 10);
        WCNS_REQUIRE(config.run.mode == wcns::RunMode::Steady);
        WCNS_REQUIRE(config.output.directory == "output");
        WCNS_REQUIRE(config.digest() != 0);
        WCNS_REQUIRE(config.summary().find("Re=") != std::string::npos);
        WCNS_REQUIRE(config.summary().find("Ma=") != std::string::npos);
    }
    {
        auto invalid = valid_config();
        const auto reduction = invalid.find("robustness.time_step_reduction = 0.4");
        invalid.replace(reduction,
                        std::string("robustness.time_step_reduction = 0.4").size(),
                        "robustness.time_step_reduction = 1.0");
        WCNS_REQUIRE_THROWS(std::invalid_argument, wcns::CaseConfig::from_text(invalid));
    }
    {
        auto two_point = valid_config();
        const auto profile = two_point.find("algorithm.profile = phenglei_wcns");
        two_point.insert(profile + std::string("algorithm.profile = phenglei_wcns").size(),
                         "\nalgorithm.flux_difference = conservative_two_point");
        const auto config = wcns::CaseConfig::from_text(two_point);
        WCNS_REQUIRE(config.flux_difference == wcns::FluxDifferenceMode::ConservativeTwoPoint);
        WCNS_REQUIRE(config.summary().find("flux_difference=conservative_two_point")
                     != std::string::npos);
        WCNS_REQUIRE(config.restart_signature().find("flux_difference=conservative_two_point")
                     != std::string::npos);

        auto invalid = valid_config();
        const auto invalid_profile = invalid.find("algorithm.profile = phenglei_wcns");
        invalid.insert(invalid_profile + std::string("algorithm.profile = phenglei_wcns").size(),
                       "\nalgorithm.flux_difference = unknown");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError, wcns::CaseConfig::from_text(invalid));
    }
    {
        auto monitored = valid_config();
        const auto enabled = monitored.find("output.statistics.enabled = false");
        monitored.replace(enabled,
                          std::string("output.statistics.enabled = false").size(),
                          "output.statistics.enabled = true\n"
                          "output.statistics.every_steps = 2\n"
                          "output.statistics.quantities = total_mass\n"
                          "output.statistics.xz_planes.enabled = true\n"
                          "output.statistics.xz_planes.cell_j_indices = 0, 23, 47");
        const auto config = wcns::CaseConfig::from_text(monitored);
        WCNS_REQUIRE(config.output.xz_planes.enabled);
        WCNS_REQUIRE(config.output.xz_planes.cell_j_indices == std::vector<int>({0, 23, 47}));
        WCNS_REQUIRE(config.output.statistics.quantities
                     == std::vector<std::string>({"total_mass",
                                                  "xz_mean_u_j0",
                                                  "xz_mass_flow_x_j0",
                                                  "xz_mean_u_j23",
                                                  "xz_mass_flow_x_j23",
                                                  "xz_mean_u_j47",
                                                  "xz_mass_flow_x_j47"}));

        auto duplicate = monitored;
        const auto indices
            = duplicate.find("output.statistics.xz_planes.cell_j_indices = 0, 23, 47");
        duplicate.replace(
            indices,
            std::string("output.statistics.xz_planes.cell_j_indices = 0, 23, 47").size(),
            "output.statistics.xz_planes.cell_j_indices = 0, 23, 23");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError, wcns::CaseConfig::from_text(duplicate));
    }
    {
        auto unsteady = valid_config();
        const auto mode = unsteady.find("run.mode = steady");
        unsteady.replace(
            mode, std::string("run.mode = steady").size(), "run.mode = unsteady\nrun.t_end = 0.25");
        const auto config = wcns::CaseConfig::from_text(unsteady);
        WCNS_REQUIRE(config.run.mode == wcns::RunMode::Unsteady);
        WCNS_REQUIRE_NEAR(config.run.end_time, 0.25, 1.0e-15);

        auto missing_end = valid_config();
        const auto missing_mode = missing_end.find("run.mode = steady");
        missing_end.replace(
            missing_mode, std::string("run.mode = steady").size(), "run.mode = unsteady");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError, wcns::CaseConfig::from_text(missing_end));
    }
    // 验证 y-z 截面目标坐标按配置顺序生成平均速度和真实质量流量列。
    {
        auto monitored = valid_config();
        const auto statistics = monitored.find("output.statistics.enabled = false");
        monitored.replace(
            statistics,
            std::string("output.statistics.enabled = false").size(),
            "output.statistics.enabled = true\n"
            "output.statistics.quantities = total_mass\n"
            "output.statistics.yz_planes.enabled = true\n"
            "output.statistics.yz_planes.target_x_coordinates = 0, 3.141592653589793");
        const auto config = wcns::CaseConfig::from_text(monitored);
        WCNS_REQUIRE(config.output.yz_planes.enabled);
        WCNS_REQUIRE(config.output.yz_planes.target_x_coordinates.size() == 2);
        WCNS_REQUIRE_NEAR(
            config.output.yz_planes.target_x_coordinates[1], 3.141592653589793, 1.0e-15);
        WCNS_REQUIRE(config.output.statistics.quantities
                     == std::vector<std::string>({"total_mass",
                                                  "yz_mean_u_plane0",
                                                  "yz_mass_flow_x_plane0",
                                                  "yz_mean_u_plane1",
                                                  "yz_mass_flow_x_plane1"}));

        auto duplicate = monitored;
        const auto targets = duplicate.find(
            "output.statistics.yz_planes.target_x_coordinates = 0, 3.141592653589793");
        duplicate.replace(
            targets,
            std::string("output.statistics.yz_planes.target_x_coordinates = 0, 3.141592653589793")
                .size(),
            "output.statistics.yz_planes.target_x_coordinates = 0, 0");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError, wcns::CaseConfig::from_text(duplicate));
    }
    // 验证槽道壁面统计会自动注册全部派生列，并拒绝无粘或非法几何配置。
    {
        auto channel = valid_config();
        const auto viscous = channel.find("run.viscous = false");
        channel.replace(viscous, std::string("run.viscous = false").size(), "run.viscous = true");
        const auto statistics = channel.find("output.statistics.enabled = false");
        channel.replace(statistics,
                        std::string("output.statistics.enabled = false").size(),
                        "output.statistics.enabled = true\n"
                        "output.statistics.quantities = total_mass\n"
                        "output.statistics.channel_walls.enabled = true\n"
                        "output.statistics.channel_walls.lower_patch = bottom\n"
                        "output.statistics.channel_walls.upper_patch = top\n"
                        "output.statistics.channel_walls.half_height = 1.0");
        const auto config = wcns::CaseConfig::from_text(channel);
        WCNS_REQUIRE(config.output.channel_walls.enabled);
        WCNS_REQUIRE(config.output.channel_walls.lower_patch == "bottom");
        WCNS_REQUIRE(config.output.channel_walls.upper_patch == "top");
        WCNS_REQUIRE_NEAR(config.output.channel_walls.half_height, 1.0, 1.0e-15);
        const auto expected_names = wcns::channel_wall_statistic_names();
        WCNS_REQUIRE(std::equal(expected_names.begin(),
                                expected_names.end(),
                                config.output.statistics.quantities.end()
                                    - static_cast<std::ptrdiff_t>(expected_names.size())));

        auto inviscid = channel;
        const auto enabled_viscous = inviscid.find("run.viscous = true");
        inviscid.replace(
            enabled_viscous, std::string("run.viscous = true").size(), "run.viscous = false");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError, wcns::CaseConfig::from_text(inviscid));
        auto invalid_height = channel;
        const auto height
            = invalid_height.find("output.statistics.channel_walls.half_height = 1.0");
        invalid_height.replace(
            height,
            std::string("output.statistics.channel_walls.half_height = 1.0").size(),
            "output.statistics.channel_walls.half_height = 0.0");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(invalid_height));
    }
    {
        auto boundary = valid_config();
        const auto viscous = boundary.find("run.viscous = false");
        boundary.replace(viscous, std::string("run.viscous = false").size(), "run.viscous = true");
        boundary += R"(
output.boundary.enabled = true
output.boundary.format = tecplot
output.boundary.every_steps = 2
output.boundary.patches = wall
output.boundary.quantities = p_w,T_w,mu_w,Cp,Cf,q_wall,pressure_traction_x,viscous_traction_x,traction_x
output.boundary.reference_pressure = 0.02
output.boundary.reference_density = 1
output.boundary.reference_velocity_x = 1
output.boundary.reference_velocity_y = 0
output.boundary.reference_velocity_z = 0
output.boundary.reference_area = 1
output.boundary.reference_length = 1
output.boundary.drag_direction_x = 1
output.boundary.drag_direction_y = 0
output.boundary.drag_direction_z = 0
output.boundary.lift_direction_x = 0
output.boundary.lift_direction_y = 1
output.boundary.lift_direction_z = 0
output.boundary.tangent_direction_x = 1
output.boundary.tangent_direction_y = 0
output.boundary.tangent_direction_z = 0
)";
        const auto config = wcns::CaseConfig::from_text(boundary);
        WCNS_REQUIRE(config.output.boundary.enabled);
        WCNS_REQUIRE(config.output.boundary.patches == std::vector<std::string>({"wall"}));
        WCNS_REQUIRE_NEAR(config.output.boundary.reference_pressure, 0.02, 1.0e-15);
        WCNS_REQUIRE(config.output.boundary.summary().find("A_ref=1") != std::string::npos);
        WCNS_REQUIRE(config.restart_signature().find("output.boundary") == std::string::npos);

        auto inviscid = boundary;
        const auto enabled_viscous = inviscid.find("run.viscous = true");
        inviscid.replace(
            enabled_viscous, std::string("run.viscous = true").size(), "run.viscous = false");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError, wcns::CaseConfig::from_text(inviscid));

        auto duplicate = boundary;
        const auto patches = duplicate.find("output.boundary.patches = wall");
        duplicate.replace(patches,
                          std::string("output.boundary.patches = wall").size(),
                          "output.boundary.patches = wall,wall");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError, wcns::CaseConfig::from_text(duplicate));
    }
    {
        const auto legacy = wcns::CaseConfig::from_text(valid_config());
        WCNS_REQUIRE(std::holds_alternative<wcns::ConstantViscosity>(legacy.transport.viscosity));
        WCNS_REQUIRE_NEAR(legacy.transport.prandtl, 0.72, 0.0);
        WCNS_REQUIRE_NEAR(
            std::get<wcns::ConstantViscosity>(legacy.transport.viscosity).viscosity_ratio,
            1.0,
            0.0);
        WCNS_REQUIRE(legacy.summary().find("law=constant") != std::string::npos);
        WCNS_REQUIRE(legacy.restart_signature().find("transport_v1") != std::string::npos);
        WCNS_REQUIRE(legacy.legacy_v1_restart_signature().find("transport_v1")
                     == std::string::npos);

        auto in_kelvin = valid_config()
            + "transport.model = sutherland\n"
              "transport.prandtl = 0.71\n"
              "transport.sutherland.reference_viscosity_ratio = 1.25\n"
              "transport.sutherland.temperature = 110.4\n";
        auto as_ratio = valid_config()
            + "transport.model = sutherland\n"
              "transport.prandtl = 0.71\n"
              "transport.sutherland.reference_viscosity_ratio = 1.25\n"
              "transport.sutherland.temperature_ratio = "
              "0.38313378448724628\n";
        const auto kelvin_config = wcns::CaseConfig::from_text(in_kelvin);
        const auto ratio_config = wcns::CaseConfig::from_text(as_ratio);
        const auto& law = std::get<wcns::SutherlandViscosity>(kelvin_config.transport.viscosity);
        WCNS_REQUIRE_NEAR(law.reference_viscosity_ratio, 1.25, 0.0);
        WCNS_REQUIRE_NEAR(law.constant_temperature_ratio, 110.4 / 288.15, 1.0e-15);
        WCNS_REQUIRE(kelvin_config.restart_signature() == ratio_config.restart_signature());
        WCNS_REQUIRE(kelvin_config.digest() != ratio_config.digest());

        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(
                                in_kelvin + "transport.sutherland.temperature_ratio = 0.383\n"));
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(valid_config()
                                        + "transport.model = constant\n"
                                          "transport.sutherland.temperature = 110.4\n"));
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(valid_config()
                                        + "transport.model = sutherland\n"
                                          "transport.sutherland.reference_viscosity_ratio = 1\n"));
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(valid_config() + "transport.model = power_law\n"));
        WCNS_REQUIRE_THROWS(
            std::invalid_argument,
            wcns::CaseConfig::from_text(valid_config() + "transport.prandtl = 0\n"));
    }
    {
        auto reordered = valid_config();
        reordered = "# ignored\n" + reordered;
        const auto first = wcns::CaseConfig::from_text(valid_config());
        const auto second = wcns::CaseConfig::from_text(reordered);
        WCNS_REQUIRE(first.digest() == second.digest());
    }
    {
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(valid_config() + "run.cfl = 0.1\n"));
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(valid_config() + "unknown.field = 1\n"));
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(valid_config() + "Re = 1000\n"));
        auto invalid_mdcd = valid_config();
        const auto diss_position = invalid_mdcd.find("algorithm.mdcd.diss = 0.005");
        invalid_mdcd.replace(diss_position,
                             std::string("algorithm.mdcd.diss = 0.005").size(),
                             "algorithm.mdcd.diss = 0.04");
        WCNS_REQUIRE_THROWS(std::invalid_argument, wcns::CaseConfig::from_text(invalid_mdcd));
        WCNS_REQUIRE_THROWS(
            wcns::CaseConfigurationError,
            wcns::CaseConfig::from_text(valid_config() + "output.history.quantities = rho\n"));
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(valid_config() + "boundary.bad.rho = 1\n"));
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(valid_config()
                                                        + "boundary.bad.type = inflow\n"
                                                          "boundary.bad.rho = 1\n"
                                                          "boundary.bad.temperature = 1\n"
                                                          "boundary.bad.pressure = 1\n"));
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(valid_config()
                                                        + "boundary.bad.type = farfield\n"
                                                          "boundary.bad.wall_temperature = 1\n"));

        auto invalid_vortex_period = valid_config();
        const auto initial_type = invalid_vortex_period.find("initial.type = uniform");
        invalid_vortex_period.replace(initial_type,
                                      std::string("initial.type = uniform").size(),
                                      "initial.type = isentropic_vortex\ninitial.period_x = -10");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError,
                            wcns::CaseConfig::from_text(invalid_vortex_period));
    }
    {
        auto two_gas_inputs = valid_config() + "gas.specific_gas_constant = 287\n";
        WCNS_REQUIRE_THROWS(wcns::PhysicsConfigurationError,
                            wcns::CaseConfig::from_text(two_gas_inputs));

        auto scmm = valid_config();
        const auto profile_position = scmm.find("phenglei_wcns");
        scmm.replace(profile_position, std::string("phenglei_wcns").size(), "scmm6_wcns");
        const auto minimum_position = scmm.find("partition.min_cells_per_active_direction = 8");
        scmm.replace(minimum_position,
                     std::string("partition.min_cells_per_active_direction = 8").size(),
                     "partition.min_cells_per_active_direction = 4");
        WCNS_REQUIRE_THROWS(wcns::CaseConfigurationError, wcns::CaseConfig::from_text(scmm));
    }
    {
        auto double_mach = valid_config();
        const auto initial = double_mach.find("initial.type = uniform");
        double_mach.replace(
            initial,
            std::string("initial.type = uniform").size(),
            "initial.type = double_mach_reflection\ninitial.x0 = 0.16666666666666667");
        double_mach += "boundary.double_mach.type = double_mach_reflection\n";
        const auto config = wcns::CaseConfig::from_text(double_mach);
        WCNS_REQUIRE(config.boundary_overrides.at("double_mach")
                     == wcns::BoundaryType::DoubleMachReflection);
        WCNS_REQUIRE(config.restart_signature().find("double_mach_reflection_v1")
                     != std::string::npos);

        auto wrong_gamma = double_mach;
        const auto gamma = wrong_gamma.find("gas.gamma = 1.4");
        wrong_gamma.replace(gamma, std::string("gas.gamma = 1.4").size(), "gas.gamma = 1.3");
        WCNS_REQUIRE_THROWS(wcns::PhysicsConfigurationError,
                            wcns::CaseConfig::from_text(wrong_gamma));
    }
}
