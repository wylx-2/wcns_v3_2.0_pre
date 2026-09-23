#include <wcns/runtime/checkpoint.hpp>

#include <wcns/io/cgns_reader.hpp>
#include <wcns/solver/euler.hpp>

#include <cgnslib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace wcns {
namespace {

constexpr int laminar_checkpoint_version = 1;
constexpr int turbulence_checkpoint_version = 2;
constexpr int implicit_checkpoint_version = 3;
const std::array<std::string, 5> checkpoint_quantities {{
    "rho",
    "rho_u",
    "rho_v",
    "rho_w",
    "rho_E",
}};
const std::array<const char*, 5> checkpoint_fields {{
    "Density",
    "MomentumX",
    "MomentumY",
    "MomentumZ",
    "EnergyStagnationDensity",
}};
const std::array<std::string, 5> previous_mean_quantities {{
    "checkpoint_previous_rho",
    "checkpoint_previous_rho_u",
    "checkpoint_previous_rho_v",
    "checkpoint_previous_rho_w",
    "checkpoint_previous_rho_E",
}};
const std::array<const char*, 5> previous_mean_fields {{
    "PreviousDensity",
    "PreviousMomentumX",
    "PreviousMomentumY",
    "PreviousMomentumZ",
    "PreviousEnergyStagnationDensity",
}};

int checkpoint_version(const CaseConfig& config)
{
    if (config.time_algorithm.integrator == TimeIntegratorKind::LuSgs
        && config.run.mode == RunMode::Unsteady) {
        return implicit_checkpoint_version;
    }
    return config.turbulence.kind == TurbulenceModelKind::None ? laminar_checkpoint_version
                                                               : turbulence_checkpoint_version;
}

class ImplicitHistoryQuantity final : public IFieldQuantity {
public:
    ImplicitHistoryQuantity(std::string name, int component, bool model)
        : descriptor_ {std::move(name),
                       TopologyLocation::Cell,
                       {},
                       "1",
                       "1",
                       QuantityScale::Dimensionless,
                       0,
                       -1}
        , component_(component)
        , model_(model)
    {
        descriptor_.validate();
    }

    const QuantityDescriptor& descriptor() const override { return descriptor_; }

    Real evaluate_cell(const StructuredBlock& block,
                       const MetricField&,
                       Index3 index,
                       const QuantityContext&) const override
    {
        const auto& history = block.flow.implicit_history;
        if (!history.valid) throw std::runtime_error("implicit checkpoint history is unavailable");
        const auto& field = model_ ? history.previous_model : history.previous_mean;
        return field(index.i, index.j, index.k, component_);
    }

private:
    QuantityDescriptor descriptor_;
    int component_ = 0;
    bool model_ = false;
};

struct ModelCheckpointField {
    TurbulenceFieldDescriptor descriptor;
    std::string field_name;
    std::string previous_quantity;
    std::string previous_field;
};

std::vector<ModelCheckpointField> model_checkpoint_fields(const CaseConfig& config)
{
    if (config.turbulence.kind == TurbulenceModelKind::None) return {};
    const auto model = TurbulenceModelRegistry::create_builtin().create(config.turbulence);
    std::vector<ModelCheckpointField> result;
    for (const auto& descriptor : model->fields()) {
        if (descriptor.role != TurbulenceFieldRole::Transported) continue;
        ModelCheckpointField field;
        field.descriptor = descriptor;
        if (descriptor.name == "nu_tilde") {
            field.field_name = "NuTilde";
            field.previous_quantity = "checkpoint_previous_rho_nu_tilde";
            field.previous_field = "PreviousNuTildeConservative";
        } else {
            field.field_name = "Model_" + descriptor.name;
            field.previous_quantity = "checkpoint_previous_model_" + descriptor.name;
            field.previous_field = "PreviousModel_" + descriptor.name;
        }
        if (field.field_name.size() > 32 || field.previous_field.size() > 32) {
            throw std::invalid_argument("turbulence checkpoint field name exceeds CGNS limit");
        }
        result.push_back(std::move(field));
    }
    return result;
}

std::string turbulence_descriptor_signature(const CaseConfig& config)
{
    if (config.turbulence.kind == TurbulenceModelKind::None) return "none";
    const auto model = TurbulenceModelRegistry::create_builtin().create(config.turbulence);
    return TurbulenceFieldSet({1, 1, 1}, 0, model->fields()).descriptor_signature();
}

void check_cgns(int status, const char* operation)
{
    if (status != CG_OK) {
        throw std::runtime_error(std::string(operation) + ": " + cg_get_error());
    }
}

class CgnsFile {
public:
    CgnsFile(const std::string& path, int mode)
    {
        check_cgns(cg_open(path.c_str(), mode, &file_), "cg_open checkpoint");
    }
    ~CgnsFile()
    {
        if (file_ >= 0) cg_close(file_);
    }
    int id() const noexcept { return file_; }
    void close()
    {
        if (file_ < 0) return;
        check_cgns(cg_close(file_), "cg_close checkpoint");
        file_ = -1;
    }

private:
    int file_ = -1;
};

std::uint64_t fnv_byte(std::uint64_t hash, unsigned char byte)
{
    hash ^= static_cast<std::uint64_t>(byte);
    return hash * 1099511628211ull;
}

void hash_text(std::uint64_t& hash, const std::string& text)
{
    for (const unsigned char byte : text)
        hash = fnv_byte(hash, byte);
}

template <class T> void hash_value(std::uint64_t& hash, const T& value)
{
    const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
    for (std::size_t index = 0; index < sizeof(T); ++index) {
        hash = fnv_byte(hash, bytes[index]);
    }
}

std::string compute_mesh_signature(const std::string& path)
{
    CgnsReader reader;
    const auto metadata = reader.read_metadata(path);
    std::uint64_t hash = 14695981039346656037ull;
    for (const auto& zone : metadata.zones) {
        hash_text(hash, zone.base_name);
        hash_text(hash, zone.name);
        hash_value(hash, zone.cell_dimension);
        hash_value(hash, zone.physical_dimension);
        hash_value(hash, zone.vertex_extent.ni);
        hash_value(hash, zone.vertex_extent.nj);
        hash_value(hash, zone.vertex_extent.nk);
        auto block = reader.read_block(path, zone, 0, 0);
        const std::array<const Array3D<Real>*, 3> coordinates {{
            &block.coordinates.x,
            &block.coordinates.y,
            &block.coordinates.z,
        }};
        for (int axis = 0; axis < zone.physical_dimension; ++axis) {
            for (int k = 0; k < zone.vertex_extent.nk; ++k) {
                for (int j = 0; j < zone.vertex_extent.nj; ++j) {
                    for (int i = 0; i < zone.vertex_extent.ni; ++i) {
                        hash_value(hash, (*coordinates[static_cast<std::size_t>(axis)])(i, j, k));
                    }
                }
            }
        }
    }
    return std::to_string(hash);
}

std::string join_path(const std::string& directory, const std::string& name)
{
    if (directory.empty()) return name;
    const char last = directory.back();
    return directory + ((last == '/' || last == '\\') ? "" : "/") + name;
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

std::string time_tag(Real time)
{
    std::ostringstream stream;
    stream << std::scientific << std::setprecision(9) << time;
    auto result = stream.str();
    for (char& character : result) {
        if (character == '.')
            character = 'p';
        else if (character == '+')
            character = 'P';
        else if (character == '-')
            character = 'M';
    }
    return result;
}

bool file_exists(const std::string& path)
{
    std::ifstream input(path, std::ios::binary);
    return static_cast<bool>(input);
}

void commit_file(const std::string& temporary, const std::string& target, bool allow_existing)
{
    if (file_exists(target)) {
        if (!allow_existing) {
            throw std::runtime_error("checkpoint already exists: " + target);
        }
        if (std::remove(target.c_str()) != 0) {
            throw std::runtime_error("cannot replace checkpoint: " + target);
        }
    }
    if (std::rename(temporary.c_str(), target.c_str()) != 0) {
        throw std::runtime_error("cannot commit checkpoint: " + target);
    }
}

void copy_file(const std::string& source, const std::string& target)
{
    std::ifstream input(source, std::ios::binary);
    std::ofstream output(target, std::ios::binary | std::ios::trunc);
    if (!input || !output) {
        throw std::runtime_error("cannot create latest checkpoint copy");
    }
    output << input.rdbuf();
    output.close();
    if (!input.good() && !input.eof()) {
        throw std::runtime_error("failed to read committed checkpoint");
    }
    if (!output) throw std::runtime_error("failed to write latest checkpoint copy");
}

std::string real_list(const std::array<Real, euler_components>& values)
{
    std::ostringstream result;
    result << std::setprecision(17);
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) result << ',';
        result << values[index];
    }
    return result.str();
}

std::string real_list(const std::vector<Real>& values)
{
    if (values.empty()) return "none";
    std::ostringstream result;
    result << std::setprecision(17);
    for (std::size_t index = 0; index < values.size(); ++index) {
        if (index != 0) result << ',';
        result << values[index];
    }
    return result.str();
}

std::array<Real, euler_components> parse_real_list(const std::string& text)
{
    std::array<Real, euler_components> result {{}};
    std::istringstream stream(text);
    std::string item;
    for (std::size_t index = 0; index < result.size(); ++index) {
        if (!std::getline(stream, item, ',')) {
            throw std::runtime_error("checkpoint residual reference is truncated");
        }
        std::size_t consumed = 0;
        result[index] = std::stod(item, &consumed);
        if (consumed != item.size() || !std::isfinite(result[index])) {
            throw std::runtime_error("checkpoint residual reference is invalid");
        }
    }
    if (std::getline(stream, item, ',')) {
        throw std::runtime_error("checkpoint residual reference has extra values");
    }
    return result;
}

std::vector<Real> parse_real_vector(const std::string& text, std::size_t expected)
{
    if (expected == 0) {
        if (text != "none") {
            throw std::runtime_error("checkpoint model residual reference has extra values");
        }
        return {};
    }
    std::vector<Real> result;
    std::istringstream stream(text);
    std::string item;
    while (std::getline(stream, item, ',')) {
        std::size_t consumed = 0;
        const Real value = std::stod(item, &consumed);
        if (consumed != item.size() || !std::isfinite(value)) {
            throw std::runtime_error("checkpoint model residual reference is invalid");
        }
        result.push_back(value);
    }
    if (result.size() != expected) {
        throw std::runtime_error("checkpoint model residual reference size differs");
    }
    return result;
}

void write_descriptor(int, const char* name, const std::string& value)
{
    check_cgns(cg_descriptor_write(name, value.c_str()), "cg_descriptor_write checkpoint");
}

std::map<std::string, std::string> read_descriptors(int file)
{
    check_cgns(cg_goto(file, 1, "end"), "cg_goto checkpoint base");
    int count = 0;
    check_cgns(cg_ndescriptors(&count), "cg_ndescriptors checkpoint");
    std::map<std::string, std::string> result;
    for (int descriptor = 1; descriptor <= count; ++descriptor) {
        char name[33] = {};
        char* value = nullptr;
        check_cgns(cg_descriptor_read(descriptor, name, &value), "cg_descriptor_read checkpoint");
        result.emplace(name, value == nullptr ? "" : value);
        if (value != nullptr) cg_free(value);
    }
    return result;
}

const std::string& required(const std::map<std::string, std::string>& values,
                            const std::string& name)
{
    const auto iterator = values.find(name);
    if (iterator == values.end()) {
        throw std::runtime_error("checkpoint is missing descriptor " + name);
    }
    return iterator->second;
}

std::size_t parse_size(const std::string& text, const char* label)
{
    std::size_t consumed = 0;
    const auto value = std::stoull(text, &consumed);
    if (consumed != text.size() || value > std::numeric_limits<std::size_t>::max()) {
        throw std::runtime_error(std::string("invalid checkpoint ") + label);
    }
    return static_cast<std::size_t>(value);
}

Real parse_real(const std::string& text, const char* label)
{
    std::size_t consumed = 0;
    const Real value = std::stod(text, &consumed);
    if (consumed != text.size() || !std::isfinite(value)) {
        throw std::runtime_error(std::string("invalid checkpoint ") + label);
    }
    return value;
}

std::size_t flat_index(Extent3 extent, int i, int j, int k)
{
    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(extent.nj)
            + static_cast<std::size_t>(j))
        * static_cast<std::size_t>(extent.ni)
        + static_cast<std::size_t>(i);
}

struct RootCheckpointData {
    CheckpointRestoreResult restored;
    bool implicit_history_valid = false;
    Real implicit_history_time_step = 0.0;
    std::string time_statistics;
    std::vector<Real> rank_payload;
    std::vector<std::size_t> rank_counts;
};

struct ZoneCheckpointFields {
    std::array<std::vector<Real>, euler_components> mean;
    std::vector<std::vector<Real>> model;
    std::array<std::vector<Real>, euler_components> previous_mean;
    std::vector<std::vector<Real>> previous_model;
};

void collective_checkpoint_action(const MpiRuntime& mpi, const std::function<void()>& action)
{
    std::string status;
    if (mpi.rank() == 0) {
        try {
            action();
            status = "OK";
        } catch (const std::exception& error) {
            status = "ERROR\n" + std::string(error.what());
        }
    }
    status = mpi.broadcast_string(std::move(status));
    if (status.rfind("ERROR\n", 0) == 0) {
        throw std::runtime_error(status.substr(6));
    }
    if (status != "OK") {
        throw std::runtime_error("invalid collective checkpoint status");
    }
}

} // namespace

CheckpointService::CheckpointService(const MpiRuntime& mpi,
                                     const CaseConfig& config,
                                     const StructuredPartitionPlan& partition,
                                     LocalBlockSet& local_blocks,
                                     const BlockMetricMap& metrics,
                                     QuantityContext quantity_context,
                                     std::string mesh_path,
                                     AcceptedTimeStatistics* time_statistics)
    : mpi_(mpi)
    , config_(config)
    , partition_(partition)
    , local_blocks_(local_blocks)
    , metrics_(metrics)
    , quantity_context_(std::move(quantity_context))
    , mesh_path_(std::move(mesh_path))
    , registry_(FieldQuantityRegistry::create_builtin())
    , time_statistics_(time_statistics)
{
    if ((time_statistics_ != nullptr) != config_.time_statistics.enabled) {
        throw std::invalid_argument(
            "checkpoint time-statistics state differs from configuration");
    }
    if (time_statistics_ != nullptr
        && time_statistics_->state().identity != config_.time_statistics_identity()) {
        throw std::invalid_argument("checkpoint time-statistics identity differs");
    }
    quantity_context_.dimensional = false;
    const auto model = TurbulenceModelRegistry::create_builtin().create(config_.turbulence);
    for (const auto& descriptor : model->fields()) {
        registry_.register_turbulence_field(descriptor);
    }
    if (checkpoint_version(config_) == implicit_checkpoint_version) {
        for (int component = 0; component < euler_components; ++component) {
            registry_.register_quantity(std::make_shared<ImplicitHistoryQuantity>(
                previous_mean_quantities[static_cast<std::size_t>(component)],
                component,
                false));
        }
        const auto model_fields = model_checkpoint_fields(config_);
        for (std::size_t component = 0; component < model_fields.size(); ++component) {
            registry_.register_quantity(std::make_shared<ImplicitHistoryQuantity>(
                model_fields[component].previous_quantity,
                static_cast<int>(component),
                true));
        }
    }
    std::string status;
    if (mpi_.rank() == 0) {
        try {
            status = "OK\n" + compute_mesh_signature(mesh_path_);
        } catch (const std::exception& error) {
            status = "ERROR\n" + std::string(error.what());
        }
    }
    status = mpi_.broadcast_string(std::move(status));
    if (status.rfind("ERROR\n", 0) == 0) {
        throw std::runtime_error(status.substr(6));
    }
    if (status.rfind("OK\n", 0) != 0) {
        throw std::runtime_error("invalid mesh signature broadcast");
    }
    mesh_signature_ = status.substr(3);
}

std::vector<std::string> CheckpointService::write(const SimulationState& state) const
{
    const auto model_fields = model_checkpoint_fields(config_);
    const std::size_t expected_model_references
        = state.steady.reference_initialized ? model_fields.size() : 0u;
    if (state.steady.model_reference_l2.size() != expected_model_references
        || state.steady.model_reference_linf.size() != expected_model_references) {
        throw std::runtime_error("checkpoint model residual reference identity differs");
    }
    const bool implicit_checkpoint
        = checkpoint_version(config_) == implicit_checkpoint_version;
    bool local_any_history = false;
    bool local_all_history = true;
    Real local_history_step = 0.0;
    for (const auto& block : local_blocks_.blocks()) {
        local_any_history = local_any_history || block.flow.implicit_history.valid;
        local_all_history = local_all_history && block.flow.implicit_history.valid;
        if (block.flow.implicit_history.valid) {
            local_history_step
                = std::max(local_history_step, block.flow.implicit_history.physical_time_step);
        }
    }
    const bool history_valid = !mpi_.all_true(!local_any_history);
    if (history_valid && !mpi_.all_true(local_all_history)) {
        throw std::runtime_error("implicit history is valid on only some checkpoint blocks");
    }
    const Real history_time_step = mpi_.max(local_history_step);
    if (implicit_checkpoint && history_valid
        && (!std::isfinite(history_time_step) || history_time_step <= 0.0)) {
        throw std::runtime_error("implicit checkpoint history time step is invalid");
    }
    std::ostringstream name;
    name << safe_name(config_.case_name) << ".checkpoint.step" << std::setw(8) << std::setfill('0')
         << state.step << ".time" << time_tag(state.time) << ".cgns";
    const auto path = join_path(config_.output.directory, name.str());
    const auto temporary = path + ".tmp";
    if (partition_.zones().empty()) {
        throw std::runtime_error("cannot write checkpoint without zones");
    }
    const int dimension = partition_.zones().front().cell_dimension;
    CgnsReader mesh_reader;
    std::unique_ptr<CgnsMeshMetadata> mesh_metadata;
    std::unique_ptr<CgnsFile> file;
    std::unordered_map<BlockId, const CgnsZoneMetadata*> mesh_zones;
    int base = 0;
    int physical_dimension = 0;
    collective_checkpoint_action(mpi_, [&] {
        mesh_metadata = std::make_unique<CgnsMeshMetadata>(mesh_reader.read_metadata(mesh_path_));
        if (mesh_metadata->zones.size() != partition_.zones().size()) {
            throw std::runtime_error("checkpoint mesh/source-zone count differs");
        }
        for (const auto& zone : mesh_metadata->zones) {
            mesh_zones.emplace(zone.block_id, &zone);
        }
        physical_dimension = mesh_metadata->zones.front().physical_dimension;
        file = std::make_unique<CgnsFile>(temporary, CG_MODE_WRITE);
        check_cgns(
            cg_base_write(file->id(), "WCNSCheckpoint", dimension, physical_dimension, &base),
            "cg_base_write checkpoint");
        check_cgns(cg_goto(file->id(), base, "end"), "cg_goto checkpoint output base");
        check_cgns(cg_dataclass_write(NormalizedByDimensional), "cg_dataclass_write checkpoint");
        write_descriptor(file->id(),
                         "WCNS_Version",
                         std::to_string(checkpoint_version(config_)));
        if (config_.turbulence.kind != TurbulenceModelKind::None) {
            write_descriptor(file->id(),
                             "WCNS_TurbulenceDescriptor",
                             turbulence_descriptor_signature(config_));
        }
        write_descriptor(file->id(), "WCNS_Step", std::to_string(state.step));
        write_descriptor(file->id(), "WCNS_Time", [&] {
            std::ostringstream value;
            value << std::setprecision(17) << state.time;
            return value.str();
        }());
        write_descriptor(file->id(), "WCNS_TimeStep", [&] {
            std::ostringstream value;
            value << std::setprecision(17) << state.time_step;
            return value.str();
        }());
        write_descriptor(file->id(), "WCNS_MeshSignature", mesh_signature_);
        write_descriptor(file->id(), "WCNS_RestartSignature", config_.restart_signature());
        if (time_statistics_ != nullptr) {
            write_descriptor(
                file->id(), "WCNS_TimeStatistics", time_statistics_->serialize());
        }
        if (implicit_checkpoint) {
            write_descriptor(file->id(), "WCNS_ImplicitHistoryValid", history_valid ? "1" : "0");
            std::ostringstream history_step;
            history_step << std::setprecision(17) << history_time_step;
            write_descriptor(
                file->id(), "WCNS_ImplicitHistoryTimeStep", history_step.str());
        }
        write_descriptor(
            file->id(), "WCNS_SteadyInitialized", state.steady.reference_initialized ? "1" : "0");
        write_descriptor(
            file->id(), "WCNS_Consecutive", std::to_string(state.steady.consecutive_passes));
        write_descriptor(file->id(), "WCNS_ReferenceL2", real_list(state.steady.reference_l2));
        write_descriptor(file->id(), "WCNS_ReferenceLinf", real_list(state.steady.reference_linf));
        if (!model_fields.empty()) {
            write_descriptor(file->id(),
                             "WCNS_ModelReferenceL2",
                             real_list(state.steady.model_reference_l2));
            write_descriptor(file->id(),
                             "WCNS_ModelReferenceLinf",
                             real_list(state.steady.model_reference_linf));
        }
    });

    for (const auto& zone : partition_.zones()) {
        if (zone.cell_dimension != dimension) {
            throw std::runtime_error("checkpoint zones use mixed dimensions");
        }
        int output_zone = 0;
        int solution = 0;
        collective_checkpoint_action(mpi_, [&] {
            const auto& mesh_zone = *mesh_zones.at(zone.source_zone);
            auto mesh_block = mesh_reader.read_block(mesh_path_, mesh_zone, 0, 0);
            std::array<cgsize_t, 9> size {{}};
            for (int axis = 0; axis < dimension; ++axis) {
                size[static_cast<std::size_t>(axis)]
                    = static_cast<cgsize_t>(mesh_zone.vertex_extent[axis]);
                size[static_cast<std::size_t>(dimension + axis)]
                    = static_cast<cgsize_t>(zone.cell_extent[axis]);
            }
            check_cgns(
                cg_zone_write(
                    file->id(), base, zone.name.c_str(), size.data(), Structured, &output_zone),
                "cg_zone_write checkpoint");
            const std::array<std::pair<const char*, const Array3D<Real>*>, 3> coordinates {{
                {"CoordinateX", &mesh_block.coordinates.x},
                {"CoordinateY", &mesh_block.coordinates.y},
                {"CoordinateZ", &mesh_block.coordinates.z},
            }};
            std::vector<Real> coordinate(mesh_zone.vertex_extent.size());
            for (int axis = 0; axis < physical_dimension; ++axis) {
                std::size_t offset = 0;
                for (int k = 0; k < mesh_zone.vertex_extent.nk; ++k) {
                    for (int j = 0; j < mesh_zone.vertex_extent.nj; ++j) {
                        for (int i = 0; i < mesh_zone.vertex_extent.ni; ++i) {
                            coordinate[offset++]
                                = (*coordinates[static_cast<std::size_t>(axis)].second)(i, j, k);
                        }
                    }
                }
                int coordinate_index = 0;
                check_cgns(cg_coord_write(file->id(),
                                          base,
                                          output_zone,
                                          RealDouble,
                                          coordinates[static_cast<std::size_t>(axis)].first,
                                          coordinate.data(),
                                          &coordinate_index),
                           "cg_coord_write checkpoint");
            }
            check_cgns(
                cg_sol_write(
                    file->id(), base, output_zone, "ConservativeState", CellCenter, &solution),
                "cg_sol_write checkpoint");
        });
        for (std::size_t component = 0; component < checkpoint_fields.size(); ++component) {
            auto values = gather_original_zone_quantity(mpi_,
                                                        local_blocks_,
                                                        metrics_,
                                                        partition_,
                                                        zone,
                                                        registry_,
                                                        checkpoint_quantities[component],
                                                        quantity_context_);
            collective_checkpoint_action(mpi_, [&] {
                int field = 0;
                check_cgns(cg_field_write(file->id(),
                                          base,
                                          output_zone,
                                          solution,
                                          RealDouble,
                                          checkpoint_fields[component],
                                          values.data(),
                                          &field),
                           "cg_field_write checkpoint");
            });
        }
        if (implicit_checkpoint && history_valid) {
            for (std::size_t component = 0; component < previous_mean_fields.size();
                 ++component) {
                auto values = gather_original_zone_quantity(mpi_,
                                                            local_blocks_,
                                                            metrics_,
                                                            partition_,
                                                            zone,
                                                            registry_,
                                                            previous_mean_quantities[component],
                                                            quantity_context_);
                collective_checkpoint_action(mpi_, [&] {
                    int field = 0;
                    check_cgns(cg_field_write(file->id(),
                                              base,
                                              output_zone,
                                              solution,
                                              RealDouble,
                                              previous_mean_fields[component],
                                              values.data(),
                                              &field),
                               "cg_field_write implicit mean history");
                });
            }
            for (const auto& model_field : model_fields) {
                auto values = gather_original_zone_quantity(mpi_,
                                                            local_blocks_,
                                                            metrics_,
                                                            partition_,
                                                            zone,
                                                            registry_,
                                                            model_field.previous_quantity,
                                                            quantity_context_);
                collective_checkpoint_action(mpi_, [&] {
                    int field = 0;
                    check_cgns(cg_field_write(file->id(),
                                              base,
                                              output_zone,
                                              solution,
                                              RealDouble,
                                              model_field.previous_field.c_str(),
                                              values.data(),
                                              &field),
                               "cg_field_write implicit model history");
                });
            }
        }
        for (const auto& model_field : model_fields) {
            auto values = gather_original_zone_quantity(mpi_,
                                                        local_blocks_,
                                                        metrics_,
                                                        partition_,
                                                        zone,
                                                        registry_,
                                                        model_field.descriptor.name,
                                                        quantity_context_);
            collective_checkpoint_action(mpi_, [&] {
                int field = 0;
                check_cgns(cg_field_write(file->id(),
                                          base,
                                          output_zone,
                                          solution,
                                          RealDouble,
                                          model_field.field_name.c_str(),
                                          values.data(),
                                          &field),
                           "cg_field_write turbulence checkpoint");
            });
        }
    }
    const auto latest = join_path(config_.output.directory,
                                  safe_name(config_.case_name) + ".checkpoint.latest.cgns");
    const auto latest_temporary = latest + ".tmp";
    collective_checkpoint_action(mpi_, [&] {
        file->close();
        commit_file(temporary, path, config_.output.allow_existing);
        copy_file(path, latest_temporary);
        // The latest file is a rolling alias owned by this run. The output
        // root policy already rejected a pre-existing directory when needed.
        commit_file(latest_temporary, latest, true);
    });
    return mpi_.rank() == 0 ? std::vector<std::string> {path, latest} : std::vector<std::string> {};
}

CheckpointRestoreResult CheckpointService::restore(const std::string& path) const
{
    const auto model_fields = model_checkpoint_fields(config_);
    RootCheckpointData root;
    std::string status;
    if (mpi_.rank() == 0) {
        try {
            CgnsFile file(path, CG_MODE_READ);
            const auto descriptors = read_descriptors(file.id());
            if (required(descriptors, "WCNS_Version")
                != std::to_string(checkpoint_version(config_))) {
                throw std::runtime_error("unsupported checkpoint version");
            }
            if (config_.turbulence.kind != TurbulenceModelKind::None
                && required(descriptors, "WCNS_TurbulenceDescriptor")
                    != turbulence_descriptor_signature(config_)) {
                throw std::runtime_error("checkpoint turbulence descriptor differs");
            }
            if (required(descriptors, "WCNS_MeshSignature") != mesh_signature_) {
                throw std::runtime_error("checkpoint mesh signature differs");
            }
            const auto stored_signature = required(descriptors, "WCNS_RestartSignature");
            const bool legacy_default_transport
                = config_.transport.restart_signature() == TransportConfig {}.restart_signature();
            const auto schema_one_signature = config_.legacy_v1_restart_signature()
                + ";transport=" + config_.transport.restart_signature();
            const bool schema_two_none_compatibility = config_.schema_version == 2
                && config_.turbulence.kind == TurbulenceModelKind::None
                && config_.time_algorithm.integrator == TimeIntegratorKind::SspRk3
                && config_.preconditioner.kind == PreconditionerKind::None
                && stored_signature == schema_one_signature;
            if (stored_signature != config_.restart_signature()
                && !schema_two_none_compatibility
                && !(legacy_default_transport
                     && stored_signature == config_.legacy_v1_restart_signature())) {
                throw std::runtime_error(
                    "checkpoint numerical signature differs (including v2 algorithms)");
            }
            root.restored.initial.step = parse_size(required(descriptors, "WCNS_Step"), "step");
            root.restored.initial.time = parse_real(required(descriptors, "WCNS_Time"), "time");
            root.restored.previous_time_step
                = parse_real(required(descriptors, "WCNS_TimeStep"), "time step");
            if (time_statistics_ != nullptr) {
                root.time_statistics = required(descriptors, "WCNS_TimeStatistics");
                const auto restored_statistics
                    = AcceptedTimeStatistics::deserialize(root.time_statistics);
                if (restored_statistics.state().identity
                    != config_.time_statistics_identity()) {
                    throw std::runtime_error(
                        "checkpoint time-statistics identity differs");
                }
            }
            if (checkpoint_version(config_) == implicit_checkpoint_version) {
                root.implicit_history_valid
                    = required(descriptors, "WCNS_ImplicitHistoryValid") == "1";
                root.implicit_history_time_step = parse_real(
                    required(descriptors, "WCNS_ImplicitHistoryTimeStep"),
                    "implicit history time step");
                if (root.implicit_history_valid
                    && root.implicit_history_time_step <= 0.0) {
                    throw std::runtime_error("checkpoint implicit history time step is invalid");
                }
            }
            root.restored.initial.steady.reference_initialized
                = required(descriptors, "WCNS_SteadyInitialized") == "1";
            root.restored.initial.steady.consecutive_passes
                = parse_size(required(descriptors, "WCNS_Consecutive"), "consecutive count");
            root.restored.initial.steady.reference_l2
                = parse_real_list(required(descriptors, "WCNS_ReferenceL2"));
            root.restored.initial.steady.reference_linf
                = parse_real_list(required(descriptors, "WCNS_ReferenceLinf"));
            if (!model_fields.empty()) {
                const std::size_t expected_model_references
                    = root.restored.initial.steady.reference_initialized ? model_fields.size() : 0u;
                root.restored.initial.steady.model_reference_l2 = parse_real_vector(
                    required(descriptors, "WCNS_ModelReferenceL2"), expected_model_references);
                root.restored.initial.steady.model_reference_linf = parse_real_vector(
                    required(descriptors, "WCNS_ModelReferenceLinf"), expected_model_references);
            }

            int zones = 0;
            check_cgns(cg_nzones(file.id(), 1, &zones), "cg_nzones checkpoint");
            if (zones != static_cast<int>(partition_.zones().size())) {
                throw std::runtime_error("checkpoint source-zone count differs");
            }
            std::unordered_map<BlockId, ZoneCheckpointFields> fields;
            for (int zone_index = 1; zone_index <= zones; ++zone_index) {
                const auto& expected = partition_.zones()[static_cast<std::size_t>(zone_index - 1)];
                char name[33] = {};
                std::array<cgsize_t, 9> size {{}};
                check_cgns(cg_zone_read(file.id(), 1, zone_index, name, size.data()),
                           "cg_zone_read checkpoint");
                if (expected.name != name) {
                    throw std::runtime_error("checkpoint source-zone name differs");
                }
                for (int axis = 0; axis < expected.cell_dimension; ++axis) {
                    if (size[static_cast<std::size_t>(expected.cell_dimension + axis)]
                        != static_cast<cgsize_t>(expected.cell_extent[axis])) {
                        throw std::runtime_error("checkpoint source-zone extent differs");
                    }
                }
                ZoneCheckpointFields zone_fields;
                const auto count = expected.cell_extent.size();
                std::array<cgsize_t, 3> lower {{1, 1, 1}};
                std::array<cgsize_t, 3> upper {{
                    static_cast<cgsize_t>(expected.cell_extent.ni),
                    static_cast<cgsize_t>(expected.cell_extent.nj),
                    static_cast<cgsize_t>(expected.cell_extent.nk),
                }};
                for (std::size_t component = 0; component < 5; ++component) {
                    zone_fields.mean[component].resize(count);
                    check_cgns(cg_field_read(file.id(),
                                             1,
                                             zone_index,
                                             1,
                                             checkpoint_fields[component],
                                             RealDouble,
                                             lower.data(),
                                             upper.data(),
                                             zone_fields.mean[component].data()),
                               "cg_field_read checkpoint");
                }
                zone_fields.model.resize(model_fields.size());
                for (std::size_t component = 0; component < model_fields.size(); ++component) {
                    zone_fields.model[component].resize(count);
                    check_cgns(cg_field_read(file.id(),
                                             1,
                                             zone_index,
                                             1,
                                             model_fields[component].field_name.c_str(),
                                             RealDouble,
                                             lower.data(),
                                             upper.data(),
                                             zone_fields.model[component].data()),
                               "cg_field_read turbulence checkpoint");
                }
                if (root.implicit_history_valid) {
                    for (std::size_t component = 0; component < euler_components; ++component) {
                        zone_fields.previous_mean[component].resize(count);
                        check_cgns(cg_field_read(file.id(),
                                                 1,
                                                 zone_index,
                                                 1,
                                                 previous_mean_fields[component],
                                                 RealDouble,
                                                 lower.data(),
                                                 upper.data(),
                                                 zone_fields.previous_mean[component].data()),
                                   "cg_field_read implicit mean history");
                    }
                    zone_fields.previous_model.resize(model_fields.size());
                    for (std::size_t component = 0; component < model_fields.size(); ++component) {
                        zone_fields.previous_model[component].resize(count);
                        check_cgns(cg_field_read(file.id(),
                                                 1,
                                                 zone_index,
                                                 1,
                                                 model_fields[component].previous_field.c_str(),
                                                 RealDouble,
                                                 lower.data(),
                                                 upper.data(),
                                                 zone_fields.previous_model[component].data()),
                                   "cg_field_read implicit model history");
                    }
                }
                fields.emplace(expected.source_zone, std::move(zone_fields));
            }

            root.rank_counts.resize(static_cast<std::size_t>(mpi_.size()));
            for (int rank = 0; rank < mpi_.size(); ++rank) {
                const auto begin = root.rank_payload.size();
                for (const auto& leaf : partition_.leaves()) {
                    if (leaf.owner != rank) continue;
                    const auto& zone = fields.at(leaf.source_zone);
                    const auto& source
                        = *std::find_if(partition_.zones().begin(),
                                        partition_.zones().end(),
                                        [&](const PartitionZone& candidate) {
                                            return candidate.source_zone == leaf.source_zone;
                                        });
                    const auto extent = leaf.cell_extent();
                    for (int component = 0; component < euler_components; ++component) {
                        for (int k = 0; k < extent.nk; ++k) {
                            for (int j = 0; j < extent.nj; ++j) {
                                for (int i = 0; i < extent.ni; ++i) {
                                    const auto global = flat_index(source.cell_extent,
                                                                   leaf.cells.begin.i + i,
                                                                   leaf.cells.begin.j + j,
                                                                   leaf.cells.begin.k + k);
                                    root.rank_payload.push_back(
                                        zone.mean[static_cast<std::size_t>(component)][global]);
                                }
                            }
                        }
                    }
                    for (std::size_t component = 0; component < model_fields.size(); ++component) {
                        for (int k = 0; k < extent.nk; ++k) {
                            for (int j = 0; j < extent.nj; ++j) {
                                for (int i = 0; i < extent.ni; ++i) {
                                    const auto global = flat_index(source.cell_extent,
                                                                   leaf.cells.begin.i + i,
                                                                   leaf.cells.begin.j + j,
                                                                   leaf.cells.begin.k + k);
                                    root.rank_payload.push_back(zone.model[component][global]);
                                }
                            }
                        }
                    }
                    if (root.implicit_history_valid) {
                        for (int component = 0; component < euler_components; ++component) {
                            for (int k = 0; k < extent.nk; ++k) {
                                for (int j = 0; j < extent.nj; ++j) {
                                    for (int i = 0; i < extent.ni; ++i) {
                                        const auto global = flat_index(source.cell_extent,
                                                                       leaf.cells.begin.i + i,
                                                                       leaf.cells.begin.j + j,
                                                                       leaf.cells.begin.k + k);
                                        root.rank_payload.push_back(
                                            zone.previous_mean[static_cast<std::size_t>(component)]
                                                              [global]);
                                    }
                                }
                            }
                        }
                        for (std::size_t component = 0; component < model_fields.size(); ++component) {
                            for (int k = 0; k < extent.nk; ++k) {
                                for (int j = 0; j < extent.nj; ++j) {
                                    for (int i = 0; i < extent.ni; ++i) {
                                        const auto global = flat_index(source.cell_extent,
                                                                       leaf.cells.begin.i + i,
                                                                       leaf.cells.begin.j + j,
                                                                       leaf.cells.begin.k + k);
                                        root.rank_payload.push_back(
                                            zone.previous_model[component][global]);
                                    }
                                }
                            }
                        }
                    }
                }
                root.rank_counts[static_cast<std::size_t>(rank)] = root.rank_payload.size() - begin;
            }
            std::ostringstream header;
            header << "OK\n"
                   << root.restored.initial.step << '\n'
                   << std::setprecision(17) << root.restored.initial.time << '\n'
                   << root.restored.previous_time_step << '\n'
                   << (root.implicit_history_valid ? 1 : 0) << '\n'
                   << root.implicit_history_time_step << '\n'
                   << (root.restored.initial.steady.reference_initialized ? 1 : 0) << '\n'
                   << root.restored.initial.steady.consecutive_passes << '\n'
                   << real_list(root.restored.initial.steady.reference_l2) << '\n'
                   << real_list(root.restored.initial.steady.reference_linf) << '\n';
            if (!model_fields.empty()) {
                header << real_list(root.restored.initial.steady.model_reference_l2) << '\n'
                       << real_list(root.restored.initial.steady.model_reference_linf) << '\n';
            }
            status = header.str();
        } catch (const std::exception& error) {
            status = "ERROR\n" + std::string(error.what());
        }
    }
    status = mpi_.broadcast_string(std::move(status));
    if (status.rfind("ERROR\n", 0) == 0) {
        throw std::runtime_error(status.substr(6));
    }
    if (status.rfind("OK\n", 0) != 0) {
        throw std::runtime_error("invalid checkpoint restore broadcast");
    }
    if (time_statistics_ != nullptr) {
        auto serialized_statistics = mpi_.broadcast_string(
            mpi_.rank() == 0 ? std::move(root.time_statistics) : std::string {});
        auto restored_statistics
            = AcceptedTimeStatistics::deserialize(serialized_statistics);
        if (restored_statistics.state().identity
            != config_.time_statistics_identity()) {
            throw std::runtime_error("checkpoint time-statistics identity differs");
        }
        *time_statistics_ = std::move(restored_statistics);
    }
    std::istringstream header(status.substr(3));
    std::string line;
    CheckpointRestoreResult result;
    std::getline(header, line);
    result.initial.step = parse_size(line, "step");
    std::getline(header, line);
    result.initial.time = parse_real(line, "time");
    std::getline(header, line);
    result.previous_time_step = parse_real(line, "time step");
    std::getline(header, line);
    const bool implicit_history_valid = line == "1";
    std::getline(header, line);
    const Real implicit_history_time_step = parse_real(line, "implicit history time step");
    std::getline(header, line);
    result.initial.steady.reference_initialized = line == "1";
    std::getline(header, line);
    result.initial.steady.consecutive_passes = parse_size(line, "consecutive count");
    std::getline(header, line);
    result.initial.steady.reference_l2 = parse_real_list(line);
    std::getline(header, line);
    result.initial.steady.reference_linf = parse_real_list(line);
    if (!model_fields.empty()) {
        const std::size_t expected_model_references
            = result.initial.steady.reference_initialized ? model_fields.size() : 0u;
        std::getline(header, line);
        result.initial.steady.model_reference_l2
            = parse_real_vector(line, expected_model_references);
        std::getline(header, line);
        result.initial.steady.model_reference_linf
            = parse_real_vector(line, expected_model_references);
    }

    const auto payload = mpi_.scatter_reals(root.rank_payload, root.rank_counts);
    std::unordered_map<BlockId, StructuredBlock*> blocks;
    for (auto& block : local_blocks_.blocks())
        blocks.emplace(block.id(), &block);
    std::size_t offset = 0;
    for (const auto& leaf : partition_.leaves()) {
        if (leaf.owner != mpi_.rank()) continue;
        auto& block = *blocks.at(leaf.block);
        const auto extent = block.cell_extent();
        const auto count = extent.size();
        if (offset + euler_components * count > payload.size()) {
            throw std::runtime_error("checkpoint rank payload is truncated");
        }
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const auto local = flat_index(extent, i, j, k);
                    ConservativeState state {{}};
                    for (int component = 0; component < euler_components; ++component) {
                        state[static_cast<std::size_t>(component)]
                            = payload[offset + static_cast<std::size_t>(component) * count + local];
                    }
                    static_cast<void>(
                        temperature_primitive_from_conservative(state,
                                                                quantity_context_.gas,
                                                                quantity_context_.reference,
                                                                quantity_context_.floors,
                                                                block.cell_dimension()));
                    store_state(block.flow.conservative, {i, j, k}, state);
                }
            }
        }
        offset += euler_components * count;
        for (const auto& model_field : model_fields) {
            if (!block.turbulence.contains(model_field.descriptor.name)
                || offset + count > payload.size()) {
                throw std::runtime_error("turbulence checkpoint rank payload is truncated");
            }
            for (int k = 0; k < extent.nk; ++k) {
                for (int j = 0; j < extent.nj; ++j) {
                    for (int i = 0; i < extent.ni; ++i) {
                        const auto local = flat_index(extent, i, j, k);
                        const Real value = payload[offset + local];
                        if (!std::isfinite(value)) {
                            throw std::runtime_error("turbulence checkpoint contains non-finite value");
                        }
                        block.turbulence.at({i, j, k}, model_field.descriptor.name) = value;
                    }
                }
            }
            offset += count;
        }
        if (implicit_history_valid) {
            const std::size_t model_count = model_fields.size();
            if (offset + (euler_components + model_count) * count > payload.size()) {
                throw std::runtime_error("implicit checkpoint history payload is truncated");
            }
            auto& history = block.flow.implicit_history;
            history.prepare(extent, static_cast<int>(model_count));
            history.valid = true;
            history.physical_time_step = implicit_history_time_step;
            for (int k = 0; k < extent.nk; ++k) {
                for (int j = 0; j < extent.nj; ++j) {
                    for (int i = 0; i < extent.ni; ++i) {
                        const auto local = flat_index(extent, i, j, k);
                        for (int component = 0; component < euler_components; ++component) {
                            history.previous_mean(i, j, k, component)
                                = payload[offset
                                          + static_cast<std::size_t>(component) * count + local];
                        }
                        for (std::size_t component = 0; component < model_count; ++component) {
                            history.previous_model(i, j, k, static_cast<int>(component))
                                = payload[offset + (euler_components + component) * count + local];
                        }
                    }
                }
            }
            offset += (euler_components + model_count) * count;
        }
    }
    if (offset != payload.size()) {
        throw std::runtime_error("checkpoint rank payload has trailing values");
    }
    return result;
}

} // namespace wcns
