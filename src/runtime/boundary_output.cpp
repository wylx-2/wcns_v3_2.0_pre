#include <wcns/runtime/boundary_output.hpp>

#include <wcns/solver/viscous_boundary.hpp>
#include <wcns/solver/viscous_flux.hpp>

#include <algorithm>
#include <array>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <sys/stat.h>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace wcns {
namespace {

constexpr std::size_t sample_width = 38;

enum SampleOffset : std::size_t {
    sample_patch = 0,
    sample_rank,
    sample_block,
    sample_zone,
    sample_axis,
    sample_side,
    sample_i,
    sample_j,
    sample_k,
    sample_x,
    sample_y,
    sample_z,
    sample_area,
    sample_nx,
    sample_ny,
    sample_nz,
    sample_pressure,
    sample_temperature,
    sample_viscosity,
    sample_cp,
    sample_cf,
    sample_heat,
    sample_pressure_tx,
    sample_pressure_ty,
    sample_pressure_tz,
    sample_viscous_tx,
    sample_viscous_ty,
    sample_viscous_tz,
    sample_total_tx,
    sample_total_ty,
    sample_total_tz,
    sample_wall_distance,
    sample_friction_velocity,
    sample_wall_y_plus,
    sample_wall_y_plus_class,
    sample_dimension,
    sample_local_ordinal,
    sample_reserved,
};

struct Sample {
    std::size_t patch = 0;
    RankId rank = 0;
    BlockId block = invalid_block_id;
    BlockId zone = invalid_block_id;
    Axis axis = Axis::I;
    Side side = Side::Lower;
    Index3 global {};
    std::array<Real, 3> center {};
    Real area = 0.0;
    std::array<Real, 3> normal {};
    BoundaryFacePhysics physics;
    int dimension = 0;
    std::size_t local_ordinal = 0;
};

using Vector3 = std::array<Real, 3>;

Real dot(const Vector3& lhs, const Vector3& rhs)
{
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

Vector3 cross(const Vector3& lhs, const Vector3& rhs)
{
    return {{
        lhs[1] * rhs[2] - lhs[2] * rhs[1],
        lhs[2] * rhs[0] - lhs[0] * rhs[2],
        lhs[0] * rhs[1] - lhs[1] * rhs[0],
    }};
}

Vector3 add(const Vector3& lhs, const Vector3& rhs)
{
    return {{lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]}};
}

Vector3 subtract(const Vector3& lhs, const Vector3& rhs)
{
    return {{lhs[0] - rhs[0], lhs[1] - rhs[1], lhs[2] - rhs[2]}};
}

Vector3 multiply(const Vector3& value, Real factor)
{
    return {{factor * value[0], factor * value[1], factor * value[2]}};
}

bool finite_vector(const Vector3& value)
{
    return std::all_of(
        value.begin(), value.end(), [](Real component) { return std::isfinite(component); });
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

void atomic_replace(const std::string& temporary, const std::string& target, bool allow_existing)
{
    if (path_exists(target)) {
        if (!allow_existing) {
            throw std::runtime_error("boundary output file already exists: " + target);
        }
        if (std::remove(target.c_str()) != 0) {
            throw std::runtime_error("cannot replace boundary output file: " + target);
        }
    }
    if (std::rename(temporary.c_str(), target.c_str()) != 0) {
        throw std::runtime_error("cannot commit boundary output file: " + target);
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

std::string time_token(Real time)
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

std::string step_token(std::size_t step)
{
    std::ostringstream stream;
    stream << std::setw(8) << std::setfill('0') << step;
    return stream.str();
}

const FaceAreaVectors& face_metrics(const MetricField& metric, Axis axis)
{
    if (axis == Axis::I) return metric.i_faces();
    if (axis == Axis::J) return metric.j_faces();
    return metric.k_faces();
}

const Array3D<Real>& face_weights(const BlockConservationWeights& weights, Axis axis)
{
    if (axis == Axis::I) return weights.physical_i_face;
    if (axis == Axis::J) return weights.physical_j_face;
    return weights.physical_k_face;
}

bool no_slip_wall(BoundaryType type)
{
    return type == BoundaryType::NoSlipAdiabaticWall || type == BoundaryType::NoSlipIsothermalWall;
}

bool requests_viscous_quantity(const BoundaryOutputConfig& config)
{
    const std::set<std::string> viscous {
        "Cf",
        "q_wall",
        "wall_distance",
        "friction_velocity",
        "wall_y_plus",
        "wall_y_plus_class",
        "viscous_traction_x",
        "viscous_traction_y",
        "viscous_traction_z",
    };
    return std::any_of(
        config.quantities.begin(), config.quantities.end(), [&](const std::string& quantity) {
            return viscous.find(quantity) != viscous.end();
        });
}

bool requests_thermal_quantity(const BoundaryOutputConfig& config)
{
    const std::set<std::string> thermal {"T_w", "mu_w"};
    return std::any_of(
        config.quantities.begin(), config.quantities.end(), [&](const std::string& quantity) {
            return thermal.find(quantity) != thermal.end();
        });
}

Index3 nearest_interior_cell(const StructuredBlock& block, Axis axis, Side side, Index3 face)
{
    face[static_cast<std::size_t>(axis)]
        = side == Side::Lower ? 0 : block.cell_extent()[static_cast<std::size_t>(axis)] - 1;
    return face;
}

Real positive_output_trace(Real high_order, Real nearest_cell, Real floor, const char* quantity)
{
    if (std::isfinite(high_order) && high_order > floor) return high_order;
    if (std::isfinite(nearest_cell) && nearest_cell > floor) return nearest_cell;
    throw PhysicsError(std::string("boundary ") + quantity
                       + " trace and nearest interior value are invalid");
}

std::size_t exact_index(Real value, const char* label)
{
    if (!std::isfinite(value) || value < 0.0
        || value > static_cast<Real>(std::numeric_limits<std::int32_t>::max())
        || std::floor(value) != value) {
        throw std::runtime_error(std::string("invalid integer in boundary payload: ") + label);
    }
    return static_cast<std::size_t>(value);
}

int exact_signed_index(Real value, const char* label)
{
    if (!std::isfinite(value) || value < static_cast<Real>(std::numeric_limits<std::int32_t>::min())
        || value > static_cast<Real>(std::numeric_limits<std::int32_t>::max())
        || std::floor(value) != value) {
        throw std::runtime_error(std::string("invalid signed integer in boundary payload: ")
                                 + label);
    }
    return static_cast<int>(value);
}

Sample decode_sample(const std::vector<Real>& values, std::size_t begin)
{
    const auto at
        = [&](SampleOffset offset) { return values[begin + static_cast<std::size_t>(offset)]; };
    Sample result;
    result.patch = exact_index(at(sample_patch), "patch");
    result.rank = static_cast<RankId>(exact_signed_index(at(sample_rank), "rank"));
    result.block = static_cast<BlockId>(exact_signed_index(at(sample_block), "block"));
    result.zone = static_cast<BlockId>(exact_signed_index(at(sample_zone), "zone"));
    result.axis = static_cast<Axis>(exact_signed_index(at(sample_axis), "axis"));
    result.side = static_cast<Side>(exact_signed_index(at(sample_side), "side"));
    result.global = {
        exact_signed_index(at(sample_i), "i"),
        exact_signed_index(at(sample_j), "j"),
        exact_signed_index(at(sample_k), "k"),
    };
    result.center = {{at(sample_x), at(sample_y), at(sample_z)}};
    result.area = at(sample_area);
    result.normal = {{at(sample_nx), at(sample_ny), at(sample_nz)}};
    result.physics.pressure = at(sample_pressure);
    result.physics.temperature = at(sample_temperature);
    result.physics.viscosity = at(sample_viscosity);
    result.physics.pressure_coefficient = at(sample_cp);
    result.physics.skin_friction_coefficient = at(sample_cf);
    result.physics.heat_flux_into_wall = at(sample_heat);
    result.physics.pressure_traction
        = {{at(sample_pressure_tx), at(sample_pressure_ty), at(sample_pressure_tz)}};
    result.physics.viscous_traction
        = {{at(sample_viscous_tx), at(sample_viscous_ty), at(sample_viscous_tz)}};
    result.physics.total_traction
        = {{at(sample_total_tx), at(sample_total_ty), at(sample_total_tz)}};
    result.physics.wall_distance = at(sample_wall_distance);
    result.physics.friction_velocity = at(sample_friction_velocity);
    result.physics.wall_y_plus = at(sample_wall_y_plus);
    result.physics.wall_y_plus_class = at(sample_wall_y_plus_class);
    result.dimension = exact_signed_index(at(sample_dimension), "dimension");
    result.local_ordinal = exact_index(at(sample_local_ordinal), "ordinal");
    return result;
}

void append_vector(std::vector<Real>& values, const Vector3& vector)
{
    values.insert(values.end(), vector.begin(), vector.end());
}

void append_sample(std::vector<Real>& values, const Sample& sample)
{
    const auto old_size = values.size();
    values.push_back(static_cast<Real>(sample.patch));
    values.push_back(static_cast<Real>(sample.rank));
    values.push_back(static_cast<Real>(sample.block));
    values.push_back(static_cast<Real>(sample.zone));
    values.push_back(static_cast<Real>(sample.axis));
    values.push_back(static_cast<Real>(sample.side));
    values.push_back(static_cast<Real>(sample.global.i));
    values.push_back(static_cast<Real>(sample.global.j));
    values.push_back(static_cast<Real>(sample.global.k));
    append_vector(values, sample.center);
    values.push_back(sample.area);
    append_vector(values, sample.normal);
    values.push_back(sample.physics.pressure);
    values.push_back(sample.physics.temperature);
    values.push_back(sample.physics.viscosity);
    values.push_back(sample.physics.pressure_coefficient);
    values.push_back(sample.physics.skin_friction_coefficient);
    values.push_back(sample.physics.heat_flux_into_wall);
    append_vector(values, sample.physics.pressure_traction);
    append_vector(values, sample.physics.viscous_traction);
    append_vector(values, sample.physics.total_traction);
    values.push_back(sample.physics.wall_distance);
    values.push_back(sample.physics.friction_velocity);
    values.push_back(sample.physics.wall_y_plus);
    values.push_back(sample.physics.wall_y_plus_class);
    values.push_back(static_cast<Real>(sample.dimension));
    values.push_back(static_cast<Real>(sample.local_ordinal));
    values.push_back(0.0);
    if (values.size() - old_size != sample_width) {
        throw std::logic_error("boundary sample payload width mismatch");
    }
}

Real selected_quantity(const Sample& sample, const std::string& name)
{
    if (name == "p_w") return sample.physics.pressure;
    if (name == "T_w") return sample.physics.temperature;
    if (name == "mu_w") return sample.physics.viscosity;
    if (name == "Cp") return sample.physics.pressure_coefficient;
    if (name == "Cf") return sample.physics.skin_friction_coefficient;
    if (name == "q_wall") return sample.physics.heat_flux_into_wall;
    if (name == "wall_distance") return sample.physics.wall_distance;
    if (name == "friction_velocity") return sample.physics.friction_velocity;
    if (name == "wall_y_plus") return sample.physics.wall_y_plus;
    if (name == "wall_y_plus_class") return sample.physics.wall_y_plus_class;
    const auto component = [&](const char* prefix, const Vector3& value) -> Real {
        const std::string base(prefix);
        if (name == base + "x") return value[0];
        if (name == base + "y") return value[1];
        if (name == base + "z") return value[2];
        return std::numeric_limits<Real>::quiet_NaN();
    };
    for (const auto entry :
         {std::pair<const char*, const Vector3&>("pressure_traction_",
                                                 sample.physics.pressure_traction),
          std::pair<const char*, const Vector3&>("viscous_traction_",
                                                 sample.physics.viscous_traction),
          std::pair<const char*, const Vector3&>("traction_", sample.physics.total_traction)}) {
        const Real value = component(entry.first, entry.second);
        if (std::isfinite(value)) return value;
    }
    throw std::logic_error("unsupported boundary quantity reached writer: " + name);
}

Real quantity_scale(const std::string& name, const QuantityContext& context)
{
    if (!context.dimensional || name == "Cp" || name == "Cf"
        || name == "wall_y_plus" || name == "wall_y_plus_class") {
        return 1.0;
    }
    const auto scales = boundary_output_scales(context, 2);
    if (name == "T_w") return scales.temperature;
    if (name == "mu_w") return scales.viscosity;
    if (name == "wall_distance") return scales.coordinate;
    if (name == "friction_velocity") return scales.velocity;
    return scales.pressure;
}

void write_metadata(std::ostream& output,
                    const CaseConfig& config,
                    const SimulationState& state,
                    int dimension,
                    const char* prefix)
{
    const auto& boundary = config.output.boundary;
    output << prefix << " case=" << config.case_name << '\n'
           << prefix << " step=" << state.step << '\n'
           << prefix << " time=" << std::setprecision(17) << state.time << '\n'
           << prefix << " dimensional=" << (config.output.dimensional ? "true" : "false") << '\n'
           << prefix << " force_per_unit_span=" << (dimension == 2 ? "true" : "false") << '\n'
           << prefix << " p_ref=" << boundary.reference_pressure << '\n'
           << prefix << " rho_ref=" << boundary.reference_density << '\n'
           << prefix << " u_ref=" << boundary.reference_velocity[0] << ','
           << boundary.reference_velocity[1] << ',' << boundary.reference_velocity[2] << '\n'
           << prefix << " A_ref=" << boundary.reference_area << '\n'
           << prefix << " L_ref_load=" << boundary.reference_length << '\n'
           << prefix << " moment_center=" << boundary.moment_center[0] << ','
           << boundary.moment_center[1] << ',' << boundary.moment_center[2] << '\n';
}

struct Loads {
    Vector3 pressure_force {};
    Vector3 viscous_force {};
    Vector3 pressure_moment {};
    Vector3 viscous_moment {};
};

std::vector<Real> load_row(const std::vector<Sample>& samples,
                           const CaseConfig& config,
                           const QuantityContext& quantities,
                           const SimulationState& state)
{
    Loads loads;
    for (const auto& sample : samples) {
        const auto pressure_force = multiply(sample.physics.pressure_traction, sample.area);
        const auto viscous_force = multiply(sample.physics.viscous_traction, sample.area);
        loads.pressure_force = add(loads.pressure_force, pressure_force);
        loads.viscous_force = add(loads.viscous_force, viscous_force);
        const auto arm = subtract(sample.center, config.output.boundary.moment_center);
        loads.pressure_moment = add(loads.pressure_moment, cross(arm, pressure_force));
        loads.viscous_moment = add(loads.viscous_moment, cross(arm, viscous_force));
    }
    const auto total_force = add(loads.pressure_force, loads.viscous_force);
    const auto total_moment = add(loads.pressure_moment, loads.viscous_moment);
    const auto& boundary = config.output.boundary;
    const Real qref = 0.5 * boundary.reference_density
        * dot(boundary.reference_velocity, boundary.reference_velocity);
    const Real force_denominator = qref * boundary.reference_area;
    const Real moment_denominator = force_denominator * boundary.reference_length;
    const int dimension = samples.front().dimension;
    const auto scales = boundary_output_scales(quantities, dimension);
    std::vector<Real> result {
        static_cast<Real>(state.step),
        state.time,
    };
    for (const auto& group : {loads.pressure_force, loads.viscous_force, total_force}) {
        append_vector(result, multiply(group, scales.force));
    }
    for (const auto& group : {loads.pressure_moment, loads.viscous_moment, total_moment}) {
        append_vector(result, multiply(group, scales.moment));
    }
    for (const auto& force : {loads.pressure_force, loads.viscous_force, total_force}) {
        result.push_back(dot(force, boundary.drag_direction) / force_denominator);
    }
    for (const auto& force : {loads.pressure_force, loads.viscous_force, total_force}) {
        result.push_back(dot(force, boundary.lift_direction) / force_denominator);
    }
    append_vector(result, multiply(total_moment, 1.0 / moment_denominator));
    return result;
}

const char* const load_columns[] = {
    "step",        "time",        "pressure_Fx", "pressure_Fy", "pressure_Fz", "viscous_Fx",
    "viscous_Fy",  "viscous_Fz",  "total_Fx",    "total_Fy",    "total_Fz",    "pressure_Mx",
    "pressure_My", "pressure_Mz", "viscous_Mx",  "viscous_My",  "viscous_Mz",  "total_Mx",
    "total_My",    "total_Mz",    "Cd_pressure", "Cd_viscous",  "Cd_total",    "Cl_pressure",
    "Cl_viscous",  "Cl_total",    "Cm_x",        "Cm_y",        "Cm_z",
};

} // namespace

BoundaryOutputScales boundary_output_scales(const QuantityContext& context, int dimension)
{
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("boundary output scales require dimension 2 or 3");
    }
    if (!context.dimensional) return {};
    const Real length = context.reference.length();
    const Real pressure = context.reference.dynamic_pressure();
    BoundaryOutputScales result;
    result.coordinate = length;
    result.area = std::pow(length, static_cast<Real>(dimension - 1));
    result.pressure = pressure;
    result.temperature = context.reference.temperature();
    result.viscosity = context.reference.viscosity();
    result.velocity = context.reference.velocity();
    result.force = pressure * result.area;
    result.moment = result.force * length;
    return result;
}

BoundaryFacePhysics evaluate_boundary_face_physics(Real pressure,
                                                   Real temperature,
                                                   Real viscosity,
                                                   const Vector3& outward_normal,
                                                   const Vector3& stress_normal,
                                                   Real thermal_coefficient,
                                                   const Vector3& temperature_gradient,
                                                   Real reynolds,
                                                   bool viscous,
                                                   const BoundaryOutputConfig& config)
{
    const Real normal_norm = dot(outward_normal, outward_normal);
    if (!std::isfinite(pressure) || !std::isfinite(temperature) || !std::isfinite(viscosity)
        || pressure <= 0.0 || temperature <= 0.0 || viscosity <= 0.0
        || !finite_vector(outward_normal) || std::abs(normal_norm - 1.0) > 1.0e-12
        || !finite_vector(stress_normal) || !finite_vector(temperature_gradient)
        || !std::isfinite(thermal_coefficient) || thermal_coefficient <= 0.0
        || !std::isfinite(reynolds) || reynolds <= 0.0) {
        throw PhysicsError("boundary face physics input is invalid");
    }
    const Real qref = 0.5 * config.reference_density
        * dot(config.reference_velocity, config.reference_velocity);
    if (!std::isfinite(qref) || qref <= 1.0e-12) {
        throw PhysicsConfigurationError("boundary reference dynamic pressure is too small");
    }
    BoundaryFacePhysics result;
    result.pressure = pressure;
    result.temperature = temperature;
    result.viscosity = viscosity;
    result.pressure_coefficient = (pressure - config.reference_pressure) / qref;
    result.pressure_traction = multiply(outward_normal, pressure);
    if (viscous) {
        result.viscous_traction = multiply(stress_normal, -1.0 / reynolds);
        const auto tangential
            = subtract(result.viscous_traction,
                       multiply(outward_normal, dot(result.viscous_traction, outward_normal)));
        result.skin_friction_coefficient = dot(tangential, config.tangent_direction) / qref;
        result.heat_flux_into_wall
            = -thermal_coefficient * dot(temperature_gradient, outward_normal) / reynolds;
    }
    result.total_traction = add(result.pressure_traction, result.viscous_traction);
    if (!finite_vector(result.pressure_traction) || !finite_vector(result.viscous_traction)
        || !finite_vector(result.total_traction) || !std::isfinite(result.pressure_coefficient)
        || !std::isfinite(result.skin_friction_coefficient)
        || !std::isfinite(result.heat_flux_into_wall)) {
        throw PhysicsError("boundary face physics result is non-finite");
    }
    return result;
}

void populate_boundary_wall_units(BoundaryFacePhysics& physics,
                                  const Vector3& outward_normal,
                                  Real density,
                                  Real wall_distance,
                                  Real reynolds,
                                  std::optional<Real> prescribed_y_plus)
{
    const Real normal_norm = dot(outward_normal, outward_normal);
    if (!std::isfinite(density) || density <= 0.0 || !std::isfinite(wall_distance)
        || wall_distance <= 0.0 || !std::isfinite(reynolds) || reynolds <= 0.0
        || !std::isfinite(physics.viscosity) || physics.viscosity <= 0.0
        || !finite_vector(physics.viscous_traction) || !finite_vector(outward_normal)
        || std::abs(normal_norm - 1.0) > 1.0e-12
        || (prescribed_y_plus
            && (!std::isfinite(*prescribed_y_plus) || *prescribed_y_plus <= 0.0))) {
        throw PhysicsError("boundary wall-unit input is invalid");
    }
    physics.wall_distance = wall_distance;
    if (prescribed_y_plus) {
        physics.wall_y_plus = *prescribed_y_plus;
        physics.friction_velocity
            = physics.wall_y_plus * physics.viscosity
            / (density * reynolds * wall_distance);
    } else {
        const auto tangential = subtract(
            physics.viscous_traction,
            multiply(outward_normal, dot(physics.viscous_traction, outward_normal)));
        const Real shear = std::sqrt(dot(tangential, tangential));
        physics.friction_velocity = std::sqrt(shear / density);
        physics.wall_y_plus = physics.friction_velocity * wall_distance * density * reynolds
            / physics.viscosity;
    }
    if (!std::isfinite(physics.friction_velocity) || physics.friction_velocity < 0.0
        || !std::isfinite(physics.wall_y_plus) || physics.wall_y_plus < 0.0) {
        throw PhysicsError("boundary wall-unit result is invalid");
    }
    // 0: viscous sublayer, 1: buffer layer, 2: log layer, 3: outer/high-y+.
    physics.wall_y_plus_class = physics.wall_y_plus <= 5.0   ? 0.0
        : physics.wall_y_plus < 30.0                         ? 1.0
        : physics.wall_y_plus <= 300.0                       ? 2.0
                                                              : 3.0;
}

BoundaryOutputWriter::BoundaryOutputWriter(const MpiRuntime& mpi,
                                           const CaseConfig& config,
                                           const StructuredPartitionPlan& partition,
                                           const LocalBlockSet& local_blocks,
                                           const StructuredMesh& global_mesh,
                                           const DistributedTopology& topology,
                                           const BlockMetricMap& metrics,
                                           const BlockBoundaryDataMap& boundary_data,
                                           const GlobalConservationWeights& conservation_weights,
                                           AlgorithmProfile profile,
                                           QuantityContext quantities)
    : mpi_(mpi)
    , config_(config)
    , partition_(partition)
    , local_blocks_(local_blocks)
    , global_mesh_(global_mesh)
    , topology_(topology)
    , metrics_(metrics)
    , boundary_data_(boundary_data)
    , conservation_weights_(conservation_weights)
    , profile_(std::move(profile))
    , quantities_(std::move(quantities))
{
    config_.output.boundary.validate(config_.run.viscous);
}

std::vector<std::string> BoundaryOutputWriter::write(const SimulationState& state)
{
    if (!config_.output.boundary.enabled) return {};
    ++version_;
    std::unordered_map<BlockId, const PartitionLeaf*> leaves;
    for (const auto& leaf : partition_.leaves()) {
        leaves.emplace(leaf.block, &leaf);
    }

    std::unordered_map<BlockId, GradientOperandFaceField> operands;
    std::unordered_map<BlockId, PrimitiveGradientField> gradients;
    if (config_.run.viscous) {
        GradientOperandFieldRegistry operand_registry;
        for (const auto& block : local_blocks_.blocks()) {
            auto [entry, inserted]
                = operands.emplace(std::piecewise_construct,
                                   std::forward_as_tuple(block.id()),
                                   std::forward_as_tuple(compute_gradient_face_operands(
                                       block, metrics_.at(block.id()), profile_, version_)));
            if (!inserted) throw std::logic_error("duplicate boundary operand block");
            operand_registry.add(block.id(), entry->second);
        }
        GradientOperandFaceHaloExchanger(
            mpi_, GradientOperandFaceHaloPlan::build(global_mesh_, profile_, version_))
            .exchange(operand_registry);
        GradientFieldRegistry gradient_registry;
        for (const auto& block : local_blocks_.blocks()) {
            auto [entry, inserted] = gradients.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(block.id()),
                std::forward_as_tuple(compute_primitive_gradients(
                    block, metrics_.at(block.id()), operands.at(block.id()), profile_)));
            if (!inserted) throw std::logic_error("duplicate boundary gradient block");
            gradient_registry.add(block.id(), entry->second);
        }
        GradientHaloExchanger(mpi_,
                              GradientHaloPlan::build(global_mesh_, topology_, profile_, version_))
            .exchange(gradient_registry);
    }

    const bool needs_viscous = requests_viscous_quantity(config_.output.boundary);
    const bool needs_thermal
        = config_.run.viscous || requests_thermal_quantity(config_.output.boundary);
    const auto turbulence_model
        = TurbulenceModelRegistry::create_builtin().create(config_.turbulence);
    std::vector<Real> local_payload;
    std::size_t ordinal = 0;
    for (const auto& block : local_blocks_.blocks()) {
        const auto leaf_iterator = leaves.find(block.id());
        if (leaf_iterator == leaves.end()) {
            throw std::logic_error("local boundary block is missing partition leaf");
        }
        const auto& leaf = *leaf_iterator->second;
        const auto& metric = metrics_.at(block.id());
        const auto& block_weights = conservation_weights_.block(block.id());
        const auto data_iterator = boundary_data_.find(block.id());
        if (data_iterator == boundary_data_.end()) {
            throw PhysicsConfigurationError("boundary data are missing for local block");
        }
        for (const auto& patch : block.boundaries) {
            const auto selected = std::find(config_.output.boundary.patches.begin(),
                                            config_.output.boundary.patches.end(),
                                            patch.name);
            if (selected == config_.output.boundary.patches.end()) continue;
            const auto patch_index
                = static_cast<std::size_t>(selected - config_.output.boundary.patches.begin());
            if (needs_viscous && !no_slip_wall(patch.type)) {
                throw PhysicsConfigurationError(
                    "viscous boundary quantities require a no-slip wall patch: " + patch.name);
            }
            const auto patch_data = data_iterator->second.find(patch.name);
            if (patch_data == data_iterator->second.end()) {
                throw PhysicsConfigurationError("boundary output data are missing for patch "
                                                + patch.name);
            }
            const auto counts = patch.boundary_face_range.counts();
            for (int ok = 0; ok < counts.nk; ++ok) {
                for (int oj = 0; oj < counts.nj; ++oj) {
                    for (int oi = 0; oi < counts.ni; ++oi) {
                        const auto face = patch.boundary_face_range.at({oi, oj, ok});
                        const auto& faces = face_metrics(metric, patch.face.axis);
                        const Real metric_area = faces.area(face.i, face.j, face.k);
                        if (!std::isfinite(metric_area) || metric_area <= 0.0) {
                            throw PhysicsError("boundary face area is invalid");
                        }
                        const Real sign = patch.face.side == Side::Lower ? -1.0 : 1.0;
                        const Vector3 normal {{
                            sign * faces.x(face.i, face.j, face.k) / metric_area,
                            sign * faces.y(face.i, face.j, face.k) / metric_area,
                            sign * faces.z(face.i, face.j, face.k) / metric_area,
                        }};
                        const Real weight
                            = face_weights(block_weights, patch.face.axis)(face.i, face.j, face.k);
                        if (!std::isfinite(weight) || weight <= 0.0) {
                            throw PhysicsError(
                                "boundary conservation quadrature weight is invalid");
                        }
                        const auto nearest
                            = nearest_interior_cell(block, patch.face.axis, patch.face.side, face);
                        Real high_order_pressure = std::numeric_limits<Real>::quiet_NaN();
                        try {
                            high_order_pressure = interpolate_internal_pressure_trace(
                                block, profile_, patch.face.axis, face);
                        } catch (const PhysicsError&) {
                            // A finite, positive nearest-cell value remains a
                            // valid output-only fallback for an invalid trace.
                        }
                        const Real pressure = positive_output_trace(
                            high_order_pressure,
                            block.flow.primitive(nearest.i, nearest.j, nearest.k, wcns::pressure),
                            quantities_.floors.pressure,
                            "pressure");
                        // Pressure-only inviscid output must not depend on a
                        // high-order temperature trace that is neither requested
                        // nor used by the Euler traction.  The neutral values are
                        // carried only in the internal gather payload and are not
                        // emitted unless a thermal quantity was requested.
                        Real temperature = 1.0;
                        Real viscosity = 1.0;
                        Vector3 stress_normal {};
                        Vector3 temperature_gradient {};
                        Real thermal_coefficient = 1.0;
                        if (needs_thermal) {
                            Real high_order_temperature = std::numeric_limits<Real>::quiet_NaN();
                            try {
                                const auto raw_temperature = interpolate_temperature_face(
                                    block, profile_, patch.face.axis, face);
                                high_order_temperature = raw_temperature[temperature_value];
                            } catch (const PhysicsError&) {
                                // Apply the same output-only fallback as pressure.
                            }
                            temperature = positive_output_trace(
                                high_order_temperature,
                                block.flow.temperature_primitive(
                                    nearest.i, nearest.j, nearest.k, temperature_value),
                                quantities_.floors.temperature,
                                "temperature");
                            viscosity = quantities_.transport.viscosity(temperature);
                            thermal_coefficient = quantities_.transport.thermal_coefficient(
                                temperature, quantities_.gas, quantities_.reference);
                        }
                        if (config_.run.viscous) {
                            auto trace = interpolate_viscous_face_trace(
                                block, gradients.at(block.id()), profile_, patch.face.axis, face);
                            trace = apply_viscous_boundary_trace(block,
                                                                 metric,
                                                                 patch,
                                                                 face,
                                                                 trace,
                                                                 pressure,
                                                                 {normal[0], normal[1], normal[2]},
                                                                 patch_data->second,
                                                                 profile_,
                                                                 quantities_.gas,
                                                                 quantities_.reference,
                                                                 quantities_.floors);
                            temperature = trace.state[temperature_value];
                            temperature_gradient
                                = trace.gradients[static_cast<int>(ViscousPrimitive::Temperature)];
                            const auto turbulence = evaluate_turbulence_viscous_face(
                                block,
                                trace,
                                profile_,
                                quantities_.transport,
                                quantities_.gas,
                                quantities_.reference,
                                patch.face.axis,
                                face,
                                &patch,
                                turbulence_model.get());
                            const auto cartesian = compute_viscous_cartesian_flux(
                                trace,
                                quantities_.transport,
                                quantities_.gas,
                                quantities_.reference,
                                quantities_.floors,
                                turbulence,
                                block.cell_dimension());
                            viscosity = cartesian.viscosity;
                            thermal_coefficient = cartesian.thermal_coefficient;
                            stress_normal = {{
                                cartesian.x[momentum_x] * normal[0]
                                    + cartesian.y[momentum_x] * normal[1]
                                    + cartesian.z[momentum_x] * normal[2],
                                cartesian.x[momentum_y] * normal[0]
                                    + cartesian.y[momentum_y] * normal[1]
                                    + cartesian.z[momentum_y] * normal[2],
                                cartesian.x[momentum_z] * normal[0]
                                    + cartesian.y[momentum_z] * normal[1]
                                    + cartesian.z[momentum_z] * normal[2],
                            }};
                        }
                        Sample sample;
                        sample.patch = patch_index;
                        sample.rank = mpi_.rank();
                        sample.block = block.id();
                        sample.zone = leaf.source_zone;
                        sample.axis = patch.face.axis;
                        sample.side = patch.face.side;
                        sample.global = {
                            face.i + leaf.cells.begin.i,
                            face.j + leaf.cells.begin.j,
                            face.k + leaf.cells.begin.k,
                        };
                        sample.center = boundary_face_coordinates(block, patch, face);
                        sample.area = weight * metric_area;
                        sample.normal = normal;
                        sample.physics
                            = evaluate_boundary_face_physics(pressure,
                                                             temperature,
                                                             viscosity,
                                                             normal,
                                                             stress_normal,
                                                             thermal_coefficient,
                                                             temperature_gradient,
                                                             quantities_.reference.reynolds(),
                                                             config_.run.viscous,
                                                             config_.output.boundary);
                        const bool rans_transport
                            = config_.turbulence.kind == TurbulenceModelKind::SaNegative
                            || config_.turbulence.kind == TurbulenceModelKind::KOmegaSst
                            || config_.turbulence.kind == TurbulenceModelKind::KEpsilon;
                        if (rans_transport && no_slip_wall(patch.type)) {
                            const Real density = block.flow.temperature_primitive(
                                nearest.i, nearest.j, nearest.k, temperature_density);
                            const Real wall_distance
                                = block.turbulence.at(nearest, "wall_distance");
                            std::optional<Real> prescribed_y_plus;
                            if (config_.turbulence.wall_treatment
                                == WallTreatment::WallFunction) {
                                prescribed_y_plus
                                    = block.turbulence.at(nearest, "wall_y_plus");
                            }
                            populate_boundary_wall_units(sample.physics,
                                                         normal,
                                                         density,
                                                         wall_distance,
                                                         quantities_.reference.reynolds(),
                                                         prescribed_y_plus);
                        }
                        sample.dimension = block.cell_dimension();
                        sample.local_ordinal = ordinal++;
                        append_sample(local_payload, sample);
                    }
                }
            }
        }
    }

    const auto gathered = mpi_.gather_reals(std::move(local_payload), 0);
    if (mpi_.rank() != 0) return {};
    if (gathered.size() % sample_width != 0) {
        throw std::runtime_error("boundary gather payload is truncated");
    }
    std::vector<Sample> samples;
    samples.reserve(gathered.size() / sample_width);
    for (std::size_t begin = 0; begin < gathered.size(); begin += sample_width) {
        samples.push_back(decode_sample(gathered, begin));
    }
    std::sort(samples.begin(), samples.end(), [](const Sample& lhs, const Sample& rhs) {
        return std::tie(lhs.patch,
                        lhs.zone,
                        lhs.global.k,
                        lhs.global.j,
                        lhs.global.i,
                        lhs.axis,
                        lhs.side)
            < std::tie(rhs.patch,
                       rhs.zone,
                       rhs.global.k,
                       rhs.global.j,
                       rhs.global.i,
                       rhs.axis,
                       rhs.side);
    });
    if (samples.empty()) {
        throw PhysicsConfigurationError("selected boundary patches have no global physical faces");
    }
    std::vector<bool> patch_seen(config_.output.boundary.patches.size(), false);
    for (std::size_t index = 0; index < samples.size(); ++index) {
        const auto& sample = samples[index];
        if (sample.patch >= patch_seen.size()) {
            throw std::runtime_error("boundary payload patch index is out of range");
        }
        patch_seen[sample.patch] = true;
        if (sample.dimension != samples.front().dimension) {
            throw std::runtime_error("boundary payload mixes cell dimensions");
        }
        if (index > 0) {
            const auto& previous = samples[index - 1];
            const bool duplicate = std::tie(sample.patch,
                                            sample.zone,
                                            sample.global.k,
                                            sample.global.j,
                                            sample.global.i,
                                            sample.axis,
                                            sample.side)
                == std::tie(previous.patch,
                            previous.zone,
                            previous.global.k,
                            previous.global.j,
                            previous.global.i,
                            previous.axis,
                            previous.side);
            if (duplicate) {
                throw std::runtime_error("duplicate physical face in boundary output gather");
            }
        }
    }
    for (std::size_t index = 0; index < patch_seen.size(); ++index) {
        if (!patch_seen[index]) {
            throw PhysicsConfigurationError("selected boundary patch was not found: "
                                            + config_.output.boundary.patches[index]);
        }
    }

    const auto basename = safe_name(config_.case_name);
    const auto extension
        = config_.output.boundary.format == SeriesOutputFormat::Tecplot ? ".dat" : ".txt";
    const auto face_path
        = join_path(config_.output.directory,
                    basename + ".boundary.r" + std::to_string(mpi_.size()) + ".step"
                        + step_token(state.step) + ".time" + time_token(state.time) + extension);
    const auto temporary = face_path + ".tmp";
    std::ofstream output(temporary, std::ios::out | std::ios::trunc);
    if (!output) {
        throw std::runtime_error("cannot open boundary output temporary file: " + temporary);
    }
    output << std::setprecision(17);
    if (config_.output.boundary.format == SeriesOutputFormat::Tecplot) {
        output << "TITLE=\"WCNS boundary faces\"\n";
        write_metadata(output, config_, state, samples.front().dimension, "#");
        output << "VARIABLES=\"patch_index\",\"source_zone\",\"global_i\","
                  "\"global_j\",\"global_k\",\"axis\",\"side\",\"x\","
                  "\"y\",\"z\",\"area\",\"nx\",\"ny\",\"nz\"";
        for (const auto& name : config_.output.boundary.quantities) {
            output << ",\"" << name << "\"";
        }
        output << "\nZONE T=\"boundary\", I=" << samples.size() << ", F=POINT\n";
    } else {
        write_metadata(output, config_, state, samples.front().dimension, "#");
        output << "# patch_index source_zone global_i global_j global_k axis side "
                  "x y z area nx ny nz";
        for (const auto& name : config_.output.boundary.quantities) {
            output << ' ' << name;
        }
        output << '\n';
    }
    const auto scales = boundary_output_scales(quantities_, samples.front().dimension);
    for (const auto& sample : samples) {
        output << sample.patch << ' ' << sample.zone << ' ' << sample.global.i << ' '
               << sample.global.j << ' ' << sample.global.k << ' ' << static_cast<int>(sample.axis)
               << ' ' << static_cast<int>(sample.side) << ' '
               << sample.center[0] * scales.coordinate << ' '
               << sample.center[1] * scales.coordinate << ' '
               << sample.center[2] * scales.coordinate << ' ' << sample.area * scales.area << ' '
               << sample.normal[0] << ' ' << sample.normal[1] << ' ' << sample.normal[2];
        for (const auto& name : config_.output.boundary.quantities) {
            output << ' ' << selected_quantity(sample, name) * quantity_scale(name, quantities_);
        }
        output << '\n';
    }
    output.close();
    if (!output) {
        throw std::runtime_error("failed to write boundary face output");
    }
    atomic_replace(temporary, face_path, config_.output.allow_existing);

    load_history_.push_back(load_row(samples, config_, quantities_, state));
    const auto load_path = join_path(config_.output.directory,
                                     basename + ".loads.r" + std::to_string(mpi_.size()) + ".txt");
    const auto load_temporary = load_path + ".tmp";
    std::ofstream loads(load_temporary, std::ios::out | std::ios::trunc);
    if (!loads) {
        throw std::runtime_error("cannot open load history temporary file: " + load_temporary);
    }
    loads << std::setprecision(17);
    write_metadata(loads, config_, state, samples.front().dimension, "#");
    loads << "#";
    for (const auto* column : load_columns)
        loads << ' ' << column;
    loads << '\n';
    for (const auto& row : load_history_) {
        for (std::size_t index = 0; index < row.size(); ++index) {
            if (index != 0) loads << ' ';
            loads << row[index];
        }
        loads << '\n';
    }
    loads.close();
    if (!loads) throw std::runtime_error("failed to write load history");
    atomic_replace(
        load_temporary, load_path, config_.output.allow_existing || load_history_created_);
    const bool first_load = !load_history_created_;
    load_history_created_ = true;
    return first_load ? std::vector<std::string> {face_path, load_path}
                      : std::vector<std::string> {face_path};
}

} // namespace wcns
