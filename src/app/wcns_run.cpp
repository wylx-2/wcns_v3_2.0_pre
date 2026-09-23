#include <wcns/io/cgns_reader.hpp>
#include <wcns/mesh/high_order_metrics.hpp>
#include <wcns/mesh/metrics.hpp>
#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/mpi_runtime.hpp>
#include <wcns/runtime/case_config.hpp>
#include <wcns/runtime/boundary_output.hpp>
#include <wcns/runtime/checkpoint.hpp>
#include <wcns/runtime/field_output.hpp>
#include <wcns/runtime/flow_initializer.hpp>
#include <wcns/runtime/output_manager.hpp>
#include <wcns/runtime/simulation_driver.hpp>
#include <wcns/runtime/structured_partition.hpp>
#include <wcns/solver/inviscid_wcns_solver.hpp>
#include <wcns/solver/sa_negative_transport.hpp>
#include <wcns/solver/rans_two_equation_transport.hpp>
#include <wcns/solver/viscous_wcns_solver.hpp>

#include <cctype>
#include <algorithm>
#include <csignal>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

namespace {

volatile std::sig_atomic_t stop_requested = 0;

void request_stop(int)
{
    stop_requested = 1;
}

struct CommandLine {
    std::string config_path;
    bool dry_run = false;
};

CommandLine parse_command_line(int argc, char** argv)
{
    CommandLine result;
    for (int index = 1; index < argc; ++index) {
        const std::string argument = argv[index];
        if (argument == "--dry-run") {
            result.dry_run = true;
        } else if (argument == "--config" && index + 1 < argc) {
            result.config_path = argv[++index];
        } else {
            throw std::invalid_argument("usage: wcns_run --config <case.wcns> [--dry-run]");
        }
    }
    if (result.config_path.empty()) {
        throw std::invalid_argument("usage: wcns_run --config <case.wcns> [--dry-run]");
    }
    return result;
}

std::string read_text(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    if (!input) throw std::runtime_error("cannot open configuration: " + path);
    std::ostringstream result;
    result << input.rdbuf();
    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed to read configuration: " + path);
    }
    return result.str();
}

bool is_absolute_path(const std::string& path)
{
    if (path.empty()) return false;
    if (path.front() == '/' || path.front() == '\\') return true;
    return path.size() >= 3 && std::isalpha(static_cast<unsigned char>(path[0])) != 0
        && path[1] == ':' && (path[2] == '/' || path[2] == '\\');
}

std::string resolve_mesh_path(const std::string& config_path, const std::string& mesh_path)
{
    if (is_absolute_path(mesh_path)) return mesh_path;
    const auto separator = config_path.find_last_of("/\\");
    if (separator == std::string::npos) return mesh_path;
    return config_path.substr(0, separator + 1) + mesh_path;
}

std::string broadcast_configuration(const wcns::MpiRuntime& mpi, const std::string& path)
{
    std::string payload;
    if (mpi.rank() == 0) {
        try {
            payload = "OK\n" + read_text(path);
        } catch (const std::exception& error) {
            payload = "ERROR\n" + std::string(error.what());
        }
    }
    payload = mpi.broadcast_string(std::move(payload));
    if (payload.rfind("ERROR\n", 0) == 0) {
        throw std::runtime_error(payload.substr(6));
    }
    if (payload.rfind("OK\n", 0) != 0) {
        throw std::runtime_error("invalid broadcast configuration envelope");
    }
    return payload.substr(3);
}

std::vector<wcns::PartitionZone> partition_zones(const wcns::CgnsMeshMetadata& metadata)
{
    std::vector<wcns::PartitionZone> result;
    result.reserve(metadata.zones.size());
    for (const auto& zone : metadata.zones) {
        result.push_back({
            zone.block_id,
            zone.name,
            zone.cell_dimension,
            zone.cell_extent,
        });
    }
    return result;
}

std::vector<wcns::CgnsPartitionLeaf> cgns_leaves(const wcns::StructuredPartitionPlan& plan)
{
    std::vector<wcns::CgnsPartitionLeaf> result;
    result.reserve(plan.leaves().size());
    for (const auto& leaf : plan.leaves()) {
        result.push_back({
            leaf.block,
            leaf.source_zone,
            leaf.cells.begin,
            leaf.cells.end,
            leaf.owner,
        });
    }
    return result;
}

std::tuple<wcns::Real, wcns::Real, wcns::Real, wcns::Real>
transport_range(const wcns::MpiRuntime& mpi,
                wcns::LocalBlockSet& local_blocks,
                const wcns::GasModel& gas,
                const wcns::ReferenceScales& reference,
                const wcns::NumericalFloors& floors,
                const wcns::TransportModel& transport)
{
    wcns::Real local_minimum_temperature = std::numeric_limits<wcns::Real>::infinity();
    wcns::Real local_maximum_temperature = -std::numeric_limits<wcns::Real>::infinity();
    for (auto& block : local_blocks.blocks()) {
        wcns::update_temperature_primitive_interior(block, gas, reference, floors);
        const auto cells = block.cell_extent();
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    const auto temperature
                        = block.flow.temperature_primitive(i, j, k, wcns::temperature_value);
                    // Evaluate every runtime temperature before any normal output.
                    static_cast<void>(transport.viscosity(temperature));
                    local_minimum_temperature = std::min(local_minimum_temperature, temperature);
                    local_maximum_temperature = std::max(local_maximum_temperature, temperature);
                }
            }
        }
    }
    const auto minimum_temperature = mpi.min(local_minimum_temperature);
    const auto maximum_temperature = mpi.max(local_maximum_temperature);
    return {
        minimum_temperature,
        maximum_temperature,
        transport.viscosity(minimum_temperature),
        transport.viscosity(maximum_temperature),
    };
}

wcns::BoundaryType configured_boundary_type(const wcns::CaseConfig& config, const std::string& name)
{
    const auto iterator = config.boundary_overrides.find(name);
    return iterator == config.boundary_overrides.end() ? config.default_boundary : iterator->second;
}

void configure_boundaries(wcns::StructuredMesh& global_mesh,
                          std::vector<wcns::StructuredBlock>& local_blocks,
                          const wcns::CaseConfig& config)
{
    for (const auto& descriptor : global_mesh.blocks()) {
        auto& block = global_mesh.block(descriptor.id());
        for (auto& patch : block.boundaries) {
            patch.type = configured_boundary_type(config, patch.name);
        }
    }
    for (auto& block : local_blocks) {
        for (auto& patch : block.boundaries) {
            patch.type = configured_boundary_type(config, patch.name);
        }
    }
}

wcns::BlockBoundaryDataMap make_boundary_data(const wcns::LocalBlockSet& local_blocks,
                                              const wcns::CaseConfig& config,
                                              const wcns::GasModel& gas,
                                              const wcns::ReferenceScales& reference,
                                              const wcns::NumericalFloors& floors)
{
    wcns::BlockBoundaryDataMap result;
    for (const auto& block : local_blocks.blocks()) {
        const auto target = wcns::FlowInitializer::evaluate(
            config.initial, {0.0, 0.0, 0.0}, gas, reference, floors, block.cell_dimension());
        wcns::BoundaryDataMap data;
        for (const auto& patch : block.boundaries) {
            wcns::BoundaryData patch_data;
            const auto configured = config.boundary_data.find(patch.name);
            const wcns::BoundaryPhysicalDataConfig* physical
                = configured == config.boundary_data.end() ? nullptr : &configured->second;
            if (physical != nullptr) {
                for (std::size_t component = 0; component < 3; ++component) {
                    patch_data.wall_velocity[component]
                        = physical->wall_velocity[component].value_or(0.0);
                }
            }
            if (patch.type == wcns::BoundaryType::Farfield
                || patch.type == wcns::BoundaryType::Inflow) {
                patch_data.target_state = target;
            }
            if (patch.type == wcns::BoundaryType::Farfield
                && config.farfield_point_vortex.enabled) {
                patch_data.farfield_point_vortex = wcns::FarfieldPointVortex {
                    config.farfield_point_vortex.lift_coefficient,
                    config.farfield_point_vortex.center,
                    config.farfield_point_vortex.chord,
                };
            }
            if (physical != nullptr && physical->has_target_state()) {
                const wcns::Real rho = *physical->rho;
                const wcns::Real u = physical->u.value_or(0.0);
                const wcns::Real v = physical->v.value_or(0.0);
                const wcns::Real w = physical->w.value_or(0.0);
                if (physical->temperature) {
                    const wcns::TemperaturePrimitiveState state {
                        rho,
                        u,
                        v,
                        w,
                        *physical->temperature,
                    };
                    static_cast<void>(wcns::pressure_primitive(
                        state, gas, reference, floors, block.cell_dimension()));
                    patch_data.target_state = state;
                } else {
                    patch_data.target_state
                        = wcns::temperature_primitive({rho, u, v, w, *physical->pressure},
                                                      gas,
                                                      reference,
                                                      floors,
                                                      block.cell_dimension());
                }
            }
            if (patch.type == wcns::BoundaryType::NoSlipIsothermalWall) {
                patch_data.wall_temperature = physical != nullptr && physical->wall_temperature
                    ? physical->wall_temperature
                    : std::optional<wcns::Real>(config.initial.parameter("temperature", 1.0));
            }
            if (patch.type == wcns::BoundaryType::DoubleMachReflection) {
                patch_data.double_mach_reflection.emplace(
                    config.initial.parameter("x0", 1.0 / 6.0));
            }
            patch_data.validate(patch.type, block.cell_dimension());
            data.emplace(patch.name, patch_data);
        }
        result.emplace(block.id(), std::move(data));
    }
    return result;
}

const wcns::CgnsZoneMetadata& source_zone_metadata(const wcns::CgnsMeshMetadata& metadata,
                                                   wcns::BlockId source_zone)
{
    const auto iterator
        = std::find_if(metadata.zones.begin(),
                       metadata.zones.end(),
                       [source_zone](const auto& zone) { return zone.block_id == source_zone; });
    if (iterator == metadata.zones.end()) {
        throw std::runtime_error("partition metric source zone is missing");
    }
    return *iterator;
}

wcns::BlockMetricMap initialize_partitioned_metrics(const wcns::MpiRuntime& mpi,
                                                    const wcns::CgnsReader& reader,
                                                    const std::string& mesh_name,
                                                    const wcns::CgnsMeshMetadata& metadata,
                                                    const wcns::StructuredPartitionPlan& plan,
                                                    wcns::LocalBlockSet& local_blocks,
                                                    const wcns::AlgorithmProfile& profile,
                                                    const wcns::MetricBuildOptions& metric_options)
{
    // Low-order physical-boundary normals remain block-local. High-order metric
    // operands are evaluated once on each original CGNS zone and then sliced,
    // so an artificial MPI cut cannot become a one-sided geometry boundary.
    for (auto& block : local_blocks.blocks())
        wcns::compute_metrics(block);

    wcns::BlockMetricMap result;
    for (const auto& zone : plan.zones()) {
        wcns::RankId metric_owner = std::numeric_limits<wcns::RankId>::max();
        for (const auto& leaf : plan.leaves()) {
            if (leaf.source_zone == zone.source_zone) {
                metric_owner = std::min(metric_owner, leaf.owner);
            }
        }
        if (metric_owner == std::numeric_limits<wcns::RankId>::max()) {
            throw std::runtime_error("partition source zone has no leaves");
        }

        std::vector<std::size_t> counts;
        std::vector<wcns::Real> owner_payload;
        if (mpi.rank() == metric_owner) {
            counts.assign(static_cast<std::size_t>(mpi.size()), 0);
            auto source_block = reader.read_block(
                mesh_name, source_zone_metadata(metadata, zone.source_zone), metric_owner, 0);
            auto initialized
                = wcns::initialize_metric_field(source_block, profile, metric_options);
            std::cout << "metric zone=" << zone.source_zone
                      << " max_jacobian_relative_difference="
                      << initialized.diagnostics.maximum_jacobian_relative_difference
                      << " fallback_cells=" << initialized.diagnostics.fallback_cell_count << '\n';
            const auto source_metric = std::move(initialized.metric);
            for (int rank = 0; rank < mpi.size(); ++rank) {
                for (const auto& leaf : plan.leaves()) {
                    if (leaf.source_zone != zone.source_zone || leaf.owner != rank) {
                        continue;
                    }
                    const auto metric = wcns::extract_metric_field(
                        source_metric, leaf.cells.begin, leaf.cell_extent());
                    const auto packed = wcns::pack_metric_field(metric);
                    counts[static_cast<std::size_t>(rank)] += packed.size();
                    owner_payload.insert(owner_payload.end(), packed.begin(), packed.end());
                }
            }
        }
        const auto local_payload = mpi.scatter_reals(owner_payload, counts, metric_owner);
        std::size_t offset = 0;
        for (const auto& leaf : plan.leaves()) {
            if (leaf.source_zone != zone.source_zone || leaf.owner != mpi.rank()) {
                continue;
            }
            const auto count
                = wcns::metric_field_payload_size(leaf.cell_extent(), leaf.cell_dimension);
            if (offset + count > local_payload.size()) {
                throw std::runtime_error("partition metric payload is truncated");
            }
            const std::vector<wcns::Real> packed(
                local_payload.begin() + static_cast<std::ptrdiff_t>(offset),
                local_payload.begin() + static_cast<std::ptrdiff_t>(offset + count));
            result.emplace(leaf.block,
                           wcns::unpack_metric_field(
                               profile.kind(), leaf.cell_extent(), leaf.cell_dimension, packed));
            offset += count;
        }
        if (offset != local_payload.size()) {
            throw std::runtime_error("partition metric payload has trailing values");
        }
    }
    if (result.size() != local_blocks.blocks().size()) {
        throw std::runtime_error("partition metric map does not cover local blocks");
    }
    return result;
}

class ConsoleObserver final : public wcns::ISimulationObserver {
public:
    explicit ConsoleObserver(const wcns::MpiRuntime& mpi)
        : mpi_(mpi)
    { }

    void on_step(const wcns::SimulationState& state, bool residual_checked) override
    {
        if (mpi_.rank() != 0) return;
        std::cout << "step=" << state.step << " time=" << std::setprecision(17) << state.time
                  << " dt=" << state.time_step << " residual=" << state.residuals.total_l2()
                  << " checked=" << (residual_checked ? "true" : "false")
                  << " stop=" << wcns::stop_reason_name(state.stop_reason) << '\n';
    }

private:
    const wcns::MpiRuntime& mpi_;
};

} // namespace

int main(int argc, char** argv)
{
    try {
        const auto command = parse_command_line(argc, argv);
        wcns::MpiRuntime mpi(argc, argv);
        const auto text = broadcast_configuration(mpi, command.config_path);
        auto config = wcns::CaseConfig::from_text(text);
        if (!mpi.all_equal(config.digest())) {
            throw std::runtime_error("case configuration digest differs across MPI ranks");
        }

        const auto mesh_name = resolve_mesh_path(command.config_path, config.mesh_path);
        wcns::CgnsReader reader;
        const auto metadata = reader.read_metadata(mesh_name);
        const auto plan = wcns::StructuredPartitionPlan::build(
            partition_zones(metadata), mpi.size(), config.partition);
        if (config.output.xz_planes.enabled) {
            wcns::validate_xz_plane_statistics(config.output.xz_planes.cell_j_indices, plan);
        }
        if (config.output.yz_planes.enabled) {
            wcns::validate_yz_plane_statistics(config.output.yz_planes.target_x_coordinates, plan);
        }
        if (!mpi.all_equal(plan.digest())) {
            throw std::runtime_error("structured partition digest differs across MPI ranks");
        }
        auto partitioned
            = reader.read_partitioned_mesh(mesh_name, cgns_leaves(plan), mpi.rank(), 3);
        configure_boundaries(partitioned.global_mesh, partitioned.local_blocks, config);
        const auto topology
            = wcns::DistributedTopology::build(partitioned.global_mesh, plan.distribution(), false);
        wcns::LocalBlockSet local_blocks(
            mpi.rank(), std::move(partitioned.local_blocks), plan.distribution());

        const auto gas = config.make_gas_model();
        const auto reference = config.make_reference_scales(gas);
        const auto profile = config.make_profile();
        const auto transport_config = config.make_transport_config();
        const wcns::TransportModel transport(transport_config);
        const wcns::NumericalFloors floors;
        auto metrics = initialize_partitioned_metrics(
            mpi, reader, mesh_name, metadata, plan, local_blocks, profile, config.metric_options);
        const auto boundary_data = make_boundary_data(local_blocks, config, gas, reference, floors);

        const auto turbulence_model
            = wcns::TurbulenceModelRegistry::create_builtin().create(config.turbulence);
        if (config.turbulence.kind == wcns::TurbulenceModelKind::SaNegative) {
            wcns::initialize_sa_negative_fields(
                mpi, local_blocks, config.turbulence, reference.reynolds());
        } else if (config.turbulence.kind == wcns::TurbulenceModelKind::KOmegaSst
                   || config.turbulence.kind == wcns::TurbulenceModelKind::KEpsilon) {
            wcns::initialize_two_equation_fields(mpi, local_blocks, config.turbulence);
        }

        wcns::QuantityContext quantity_context {
            gas,
            reference,
            floors,
            transport,
            config.output.dimensional,
        };
        const auto conservation_weights
            = wcns::GlobalConservationWeights::build(partitioned.global_mesh, profile);
        wcns::CheckpointService checkpoint(
            mpi, config, plan, local_blocks, metrics, quantity_context, mesh_name);
        wcns::SimulationInitialState simulation_initial;
        if (config.restart_path.empty()) {
            wcns::FlowInitializer::initialize_local_blocks(
                local_blocks, metrics, config.initial, gas, reference, floors);
        } else {
            const auto restart_name = resolve_mesh_path(command.config_path, config.restart_path);
            simulation_initial = checkpoint.restore(restart_name).initial;
        }

        const auto [minimum_temperature,
                    maximum_temperature,
                    minimum_temperature_viscosity,
                    maximum_temperature_viscosity]
            = transport_range(mpi, local_blocks, gas, reference, floors, transport);

        if (mpi.rank() == 0) {
            std::cout << config.summary() << '\n'
                      << plan.summary() << '\n'
                      << "mesh_signature=" << checkpoint.mesh_signature() << '\n'
                      << "derived Re=" << std::setprecision(17) << reference.reynolds()
                      << " Ma=" << reference.mach() << '\n'
                      << "transport_reference mu_Tref_over_mu_ref=" << transport.viscosity(1.0)
                      << " Pr=" << transport_config.prandtl;
            if (const auto* sutherland
                = std::get_if<wcns::SutherlandViscosity>(&transport_config.viscosity)) {
                std::cout << " S_over_Tref=" << sutherland->constant_temperature_ratio;
            } else {
                std::cout << " S_over_Tref=not_applicable";
            }
            std::cout << " T_range=[" << minimum_temperature << ',' << maximum_temperature
                      << "] mu_range=[" << minimum_temperature_viscosity << ','
                      << maximum_temperature_viscosity << "]\n";
        }
        if (command.dry_run) {
            if (mpi.rank() == 0) {
                std::cout << "WCNS dry-run completed\n";
            }
            return EXIT_SUCCESS;
        }

        std::signal(SIGINT, request_stop);
        std::signal(SIGTERM, request_stop);
        ConsoleObserver console(mpi);
        wcns::StatisticContext statistic_context {
            mpi,
            local_blocks,
            metrics,
            plan,
            conservation_weights,
            profile,
            quantity_context,
            &boundary_data,
            config.run.viscous,
        };
        auto field_registry = wcns::FieldQuantityRegistry::create_builtin();
        for (const auto& descriptor : turbulence_model->fields()) {
            field_registry.register_turbulence_field(descriptor);
        }
        wcns::ProductionFieldWriter field_writer(mpi,
                                                 config,
                                                 plan,
                                                 local_blocks,
                                                 metrics,
                                                 quantity_context,
                                                 mesh_name,
                                                 std::move(field_registry));
        wcns::BoundaryOutputWriter boundary_writer(mpi,
                                                   config,
                                                   plan,
                                                   local_blocks,
                                                   partitioned.global_mesh,
                                                   topology,
                                                   metrics,
                                                   boundary_data,
                                                   conservation_weights,
                                                   profile,
                                                   quantity_context);
        auto statistic_registry = wcns::StatisticRegistry::create_builtin();
        if (config.output.xz_planes.enabled) {
            wcns::register_xz_plane_statistics(statistic_registry,
                                               config.output.xz_planes.cell_j_indices);
        }
        if (config.output.yz_planes.enabled) {
            wcns::register_yz_plane_statistics(statistic_registry,
                                               config.output.yz_planes.target_x_coordinates);
        }
        if (config.output.channel_walls.enabled) {
            wcns::register_channel_wall_statistics(statistic_registry,
                                                   config.output.channel_walls.lower_patch,
                                                   config.output.channel_walls.upper_patch,
                                                   config.output.channel_walls.half_height);
        }
        wcns::RuntimeOutputManager output(
            mpi,
            config,
            plan,
            checkpoint.mesh_signature(),
            &statistic_context,
            [&](wcns::OutputCategory category, const wcns::SimulationState& state, bool, bool)
                -> std::vector<std::string> {
                if (category == wcns::OutputCategory::Field) {
                    return field_writer.write(state);
                }
                if (category == wcns::OutputCategory::Checkpoint) {
                    return checkpoint.write(state);
                }
                if (category == wcns::OutputCategory::Boundary) {
                    return boundary_writer.write(state);
                }
                return {};
            },
            std::move(statistic_registry));
        wcns::CompositeSimulationObserver observer;
        observer.add(console);
        observer.add(output);
        wcns::SimulationState final_state;
        if (config.run.viscous) {
            wcns::ViscousWcnsConfig solver_config;
            solver_config.inviscid = config.make_inviscid_config();
            solver_config.transport = transport_config;
            solver_config.turbulence = config.turbulence;
            wcns::ViscousWcnsSolver solver(mpi,
                                           local_blocks,
                                           partitioned.global_mesh,
                                           topology,
                                           plan.distribution().rank_count(),
                                           metrics,
                                           boundary_data,
                                           profile,
                                           gas,
                                           reference,
                                           floors,
                                           solver_config);
            wcns::ViscousSimulationSolver adapter(
                solver,
                mpi,
                local_blocks,
                metrics,
                plan,
                profile,
                config.time_algorithm,
                config.run.mode);
            wcns::SimulationDriver driver(
                mpi, adapter, config.run, observer, [] { return stop_requested != 0; });
            final_state = driver.run(simulation_initial);
        } else {
            auto solver_config = config.make_inviscid_config();
            wcns::InviscidWcnsSolver solver(mpi,
                                            local_blocks,
                                            partitioned.global_mesh,
                                            topology,
                                            plan.distribution().rank_count(),
                                            metrics,
                                            boundary_data,
                                            profile,
                                            gas,
                                            reference,
                                            floors,
                                            solver_config);
            wcns::InviscidSimulationSolver adapter(
                solver,
                mpi,
                local_blocks,
                metrics,
                plan,
                profile,
                config.time_algorithm,
                config.run.mode);
            wcns::SimulationDriver driver(
                mpi, adapter, config.run, observer, [] { return stop_requested != 0; });
            final_state = driver.run(simulation_initial);
        }
        if (mpi.rank() == 0) {
            std::cout << "WCNS run stopped: reason="
                      << wcns::stop_reason_name(final_state.stop_reason)
                      << " step=" << final_state.step << " time=" << std::setprecision(17)
                      << final_state.time << '\n';
        }
        return wcns::stop_reason_exit_code(final_state.stop_reason);
    } catch (const std::exception& error) {
        std::cerr << "wcns_run: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
