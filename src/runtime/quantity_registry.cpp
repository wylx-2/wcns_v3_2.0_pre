#include <wcns/runtime/quantity_registry.hpp>

#include <wcns/mesh/conservation_weights.hpp>
#include <wcns/solver/euler.hpp>
#include <wcns/solver/viscous_boundary.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace wcns {
namespace {

using FieldEvaluator = std::function<Real(
    const StructuredBlock&, const MetricField&, Index3, const QuantityContext&)>;

class FunctionalFieldQuantity final : public IFieldQuantity {
public:
    FunctionalFieldQuantity(QuantityDescriptor descriptor, FieldEvaluator evaluator)
        : descriptor_(std::move(descriptor))
        , evaluator_(std::move(evaluator))
    {
        descriptor_.validate();
        if (!evaluator_) throw std::invalid_argument("field evaluator is empty");
    }

    const QuantityDescriptor& descriptor() const override { return descriptor_; }
    Real evaluate_cell(const StructuredBlock& block,
                       const MetricField& metric,
                       Index3 index,
                       const QuantityContext& context) const override
    {
        return evaluator_(block, metric, index, context);
    }

private:
    QuantityDescriptor descriptor_;
    FieldEvaluator evaluator_;
};

using StatisticEvaluator = std::function<Real(const StatisticContext&)>;

class FunctionalStatisticQuantity final : public IStatisticQuantity {
public:
    FunctionalStatisticQuantity(QuantityDescriptor descriptor, StatisticEvaluator evaluator)
        : descriptor_(std::move(descriptor))
        , evaluator_(std::move(evaluator))
    {
        descriptor_.validate();
        if (!evaluator_) throw std::invalid_argument("statistic evaluator is empty");
    }

    const QuantityDescriptor& descriptor() const override { return descriptor_; }
    Real evaluate(const StatisticContext& context) const override { return evaluator_(context); }

private:
    QuantityDescriptor descriptor_;
    StatisticEvaluator evaluator_;
};

TemperaturePrimitiveState
cell_temperature_state(const StructuredBlock& block, Index3 index, const QuantityContext& context)
{
    return temperature_primitive_from_conservative(
        load_conservative(block.flow.conservative, index),
        context.gas,
        context.reference,
        context.floors,
        block.cell_dimension());
}

PressurePrimitiveState
cell_pressure_state(const StructuredBlock& block, Index3 index, const QuantityContext& context)
{
    return pressure_primitive(cell_temperature_state(block, index, context),
                              context.gas,
                              context.reference,
                              context.floors,
                              block.cell_dimension());
}

QuantityDescriptor
make_descriptor(std::string name, std::string dimensional_unit, QuantityScale scale)
{
    QuantityDescriptor result;
    result.name = std::move(name);
    result.dimensional_unit = std::move(dimensional_unit);
    result.scale = scale;
    return result;
}

std::shared_ptr<const IFieldQuantity> field(QuantityDescriptor information,
                                            FieldEvaluator evaluator)
{
    return std::make_shared<FunctionalFieldQuantity>(std::move(information), std::move(evaluator));
}

std::shared_ptr<const IStatisticQuantity> statistic(QuantityDescriptor information,
                                                    StatisticEvaluator evaluator)
{
    return std::make_shared<FunctionalStatisticQuantity>(std::move(information),
                                                         std::move(evaluator));
}

Real weighted_conservative_integral(const StatisticContext& context, int component)
{
    Real local = 0.0;
    for (const auto& block : context.local_blocks.blocks()) {
        const auto metric_iterator = context.metrics.find(block.id());
        if (metric_iterator == context.metrics.end()) {
            throw std::invalid_argument("statistic is missing metric data");
        }
        const auto& weights = context.conservation_weights.block(block.id()).cell;
        const auto& jacobian = metric_iterator->second.jacobian();
        const auto extent = block.cell_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Real weight = weights(i, j, k) * jacobian(i, j, k);
                    const Real value = block.flow.conservative(i, j, k, component);
                    if (!std::isfinite(weight) || weight <= 0.0 || !std::isfinite(value)) {
                        throw PhysicsError("statistic encountered an invalid cell value");
                    }
                    local += weight * value;
                }
            }
        }
    }
    return context.mpi.sum(local);
}

struct ChannelWallAverage {
    Real area = 0.0;
    Real shear = 0.0;
    Real density = 0.0;
    Real viscosity_ratio = 0.0;
};

ChannelWallAverage channel_wall_average(const StatisticContext& context,
                                        const std::string& patch_name,
                                        Side expected_side)
{
    if (!context.viscous || context.boundary_data == nullptr) {
        throw PhysicsError("channel-wall statistics require a viscous run and boundary data");
    }
    Real local_area = 0.0;
    Real local_shear = 0.0;
    Real local_density = 0.0;
    Real local_viscosity = 0.0;
    for (const auto& block : context.local_blocks.blocks()) {
        if (block.cell_dimension() != 3) {
            throw PhysicsError("channel-wall statistics require a 3D mesh");
        }
        const auto metric_iterator = context.metrics.find(block.id());
        const auto data_block_iterator = context.boundary_data->find(block.id());
        if (metric_iterator == context.metrics.end()
            || data_block_iterator == context.boundary_data->end()) {
            throw PhysicsError("channel-wall statistics are missing metric or boundary data");
        }
        const auto& metric = metric_iterator->second;
        for (const auto& patch : block.boundaries) {
            if (patch.name != patch_name) continue;
            if (patch.type != BoundaryType::NoSlipIsothermalWall || patch.face.axis != Axis::J
                || patch.face.side != expected_side) {
                throw PhysicsError(
                    "channel-wall statistic patch is not the requested isothermal J wall");
            }
            const auto data_iterator = data_block_iterator->second.find(patch.name);
            if (data_iterator == data_block_iterator->second.end()
                || !data_iterator->second.wall_temperature) {
                throw PhysicsError("channel-wall statistic is missing isothermal wall data");
            }
            const auto& wall_data = data_iterator->second;
            const int stencil_count
                = context.profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 4 : 6;
            if (block.cell_extent().nj < stencil_count) {
                throw ProfileError("channel-wall derivative stencil does not fit the block");
            }
            const auto counts = patch.boundary_face_range.counts();
            for (int ok = 0; ok < counts.nk; ++ok) {
                for (int oj = 0; oj < counts.nj; ++oj) {
                    for (int oi = 0; oi < counts.ni; ++oi) {
                        const auto face = patch.boundary_face_range.at({oi, oj, ok});
                        const auto& faces = metric.j_faces();
                        const Real area = faces.area(face.i, face.j, face.k);
                        const Real inward_sign = expected_side == Side::Lower ? 1.0 : -1.0;
                        const std::array<Real, 3> inward {{
                            inward_sign * faces.x(face.i, face.j, face.k) / area,
                            inward_sign * faces.y(face.i, face.j, face.k) / area,
                            inward_sign * faces.z(face.i, face.j, face.k) / area,
                        }};
                        if (std::abs(inward[0]) > 1.0e-10 || std::abs(inward[2]) > 1.0e-10
                            || std::abs(std::abs(inward[1]) - 1.0) > 1.0e-10) {
                            throw PhysicsError("channel-wall statistics require planar x-z walls");
                        }
                        auto first = face;
                        first.j = expected_side == Side::Lower ? 0 : block.cell_extent().nj - 1;
                        const auto wall = boundary_face_coordinates(block, patch, face);
                        const auto& centers = metric.cell_coordinates();
                        const std::array<Real, 3> delta {{
                            centers.x(first.i, first.j, first.k) - wall[0],
                            centers.y(first.i, first.j, first.k) - wall[1],
                            centers.z(first.i, first.j, first.k) - wall[2],
                        }};
                        const Real distance
                            = delta[0] * inward[0] + delta[1] * inward[1] + delta[2] * inward[2];
                        if (!std::isfinite(area) || area <= 0.0 || !std::isfinite(distance)
                            || distance <= 0.0) {
                            throw PhysicsError(
                                "channel-wall statistic encountered invalid geometry");
                        }
                        std::vector<Real> velocities;
                        velocities.reserve(static_cast<std::size_t>(stencil_count));
                        for (int ordinal = 0; ordinal < stencil_count; ++ordinal) {
                            auto cell = face;
                            cell.j = expected_side == Side::Lower
                                ? ordinal
                                : block.cell_extent().nj - 1 - ordinal;
                            velocities.push_back(block.flow.temperature_primitive(
                                cell.i, cell.j, cell.k, temperature_velocity_x));
                        }
                        const Real inward_derivative = 0.5 / distance
                            * wall_dirichlet_computational_derivative(wall_data.wall_velocity[0],
                                                                      velocities,
                                                                      context.profile);
                        const Real wall_temperature = *wall_data.wall_temperature;
                        const Real wall_pressure = interpolate_internal_pressure_trace(
                            block, context.profile, Axis::J, face);
                        const Real wall_density = context.quantities.gas.gamma()
                            * context.quantities.reference.mach()
                            * context.quantities.reference.mach() * wall_pressure
                            / wall_temperature;
                        const Real viscosity_ratio
                            = context.quantities.transport.viscosity(wall_temperature);
                        const Real shear = viscosity_ratio / context.quantities.reference.reynolds()
                            * inward_derivative;
                        if (!std::isfinite(wall_density) || wall_density <= 0.0
                            || !std::isfinite(shear)) {
                            throw PhysicsError(
                                "channel-wall statistic encountered invalid physics");
                        }
                        local_area += area;
                        local_shear += area * std::abs(shear);
                        local_density += area * wall_density;
                        local_viscosity += area * viscosity_ratio;
                    }
                }
            }
        }
    }
    ChannelWallAverage result;
    result.area = context.mpi.sum(local_area);
    result.shear = context.mpi.sum(local_shear);
    result.density = context.mpi.sum(local_density);
    result.viscosity_ratio = context.mpi.sum(local_viscosity);
    if (!std::isfinite(result.area) || result.area <= 0.0) {
        throw PhysicsError("channel-wall statistic patch has no finite global support: "
                           + patch_name);
    }
    result.shear /= result.area;
    result.density /= result.area;
    result.viscosity_ratio /= result.area;
    return result;
}

ChannelWallAverage combined_channel_walls(const ChannelWallAverage& lower,
                                          const ChannelWallAverage& upper)
{
    const Real area = lower.area + upper.area;
    return {
        area,
        (lower.area * lower.shear + upper.area * upper.shear) / area,
        (lower.area * lower.density + upper.area * upper.density) / area,
        (lower.area * lower.viscosity_ratio + upper.area * upper.viscosity_ratio) / area,
    };
}

} // namespace

void QuantityDescriptor::validate() const
{
    if (name.empty()) throw std::invalid_argument("quantity name is empty");
    if (location != TopologyLocation::Cell) {
        throw std::invalid_argument(
            "stage N quantity registry currently accepts cell quantities only");
    }
    if (scale == QuantityScale::LengthPower && length_power != -1 && length_power <= 0) {
        throw std::invalid_argument(
            "length-power quantity requires a positive power or dimension sentinel");
    }
    if (integration_length_power < -1) {
        throw std::invalid_argument(
            "statistic integration length power must be -1 or non-negative");
    }
    std::set<std::string> unique;
    for (const auto& dependency : dependencies) {
        if (dependency.empty() || !unique.insert(dependency).second) {
            throw std::invalid_argument("quantity dependency list is invalid");
        }
    }
}

Real quantity_scale_factor(const QuantityDescriptor& information,
                           const QuantityContext& context,
                           int dimension)
{
    if (!context.dimensional) return 1.0;
    switch (information.scale) {
    case QuantityScale::Dimensionless: return 1.0;
    case QuantityScale::Density: return context.reference.density();
    case QuantityScale::Velocity: return context.reference.velocity();
    case QuantityScale::Pressure:
    case QuantityScale::Energy: return context.reference.dynamic_pressure();
    case QuantityScale::SpecificEnergy:
        return context.reference.velocity() * context.reference.velocity();
    case QuantityScale::Temperature: return context.reference.temperature();
    case QuantityScale::Momentum: return context.reference.density() * context.reference.velocity();
    case QuantityScale::Viscosity: return context.reference.viscosity();
    case QuantityScale::KinematicViscosity:
        return context.reference.velocity() * context.reference.length();
    case QuantityScale::InverseTime:
        return context.reference.velocity() / context.reference.length();
    case QuantityScale::Dissipation:
        return std::pow(context.reference.velocity(), 3) / context.reference.length();
    case QuantityScale::LengthPower:
        return std::pow(context.reference.length(),
                        information.length_power == -1 ? dimension : information.length_power);
    }
    throw std::invalid_argument("invalid quantity scale");
}

FieldQuantityRegistry FieldQuantityRegistry::create_builtin()
{
    FieldQuantityRegistry result;
    const auto conservative = [](int component) {
        return [component](const StructuredBlock& block,
                           const MetricField&,
                           Index3 index,
                           const QuantityContext&) {
            return block.flow.conservative(index.i, index.j, index.k, component);
        };
    };
    result.register_quantity(field(make_descriptor("rho", "kg/m^3", QuantityScale::Density),
                                   [](const StructuredBlock& block,
                                      const MetricField&,
                                      Index3 index,
                                      const QuantityContext& context) {
                                       return cell_temperature_state(
                                           block, index, context)[temperature_density];
                                   }));
    const std::array<std::pair<const char*, int>, 3> velocities {{
        {"u", temperature_velocity_x},
        {"v", temperature_velocity_y},
        {"w", temperature_velocity_z},
    }};
    for (const auto& item : velocities) {
        result.register_quantity(
            field(make_descriptor(item.first, "m/s", QuantityScale::Velocity),
                  [component = item.second](const StructuredBlock& block,
                                            const MetricField&,
                                            Index3 index,
                                            const QuantityContext& context) {
                      return cell_temperature_state(
                          block, index, context)[static_cast<std::size_t>(component)];
                  }));
    }
    result.register_quantity(field(make_descriptor("p", "Pa", QuantityScale::Pressure),
                                   [](const StructuredBlock& block,
                                      const MetricField&,
                                      Index3 index,
                                      const QuantityContext& context) {
                                       return cell_pressure_state(block, index, context)[pressure];
                                   }));
    result.register_quantity(field(make_descriptor("T", "K", QuantityScale::Temperature),
                                   [](const StructuredBlock& block,
                                      const MetricField&,
                                      Index3 index,
                                      const QuantityContext& context) {
                                       return cell_temperature_state(
                                           block, index, context)[temperature_value];
                                   }));
    result.register_quantity(field(make_descriptor("rho_u", "kg/(m^2*s)", QuantityScale::Momentum),
                                   conservative(momentum_x)));
    result.register_quantity(field(make_descriptor("rho_v", "kg/(m^2*s)", QuantityScale::Momentum),
                                   conservative(momentum_y)));
    result.register_quantity(field(make_descriptor("rho_w", "kg/(m^2*s)", QuantityScale::Momentum),
                                   conservative(momentum_z)));
    result.register_quantity(
        field(make_descriptor("rho_E", "Pa", QuantityScale::Energy), conservative(total_energy)));
    result.register_quantity(field(make_descriptor("sound_speed", "m/s", QuantityScale::Velocity),
                                   [](const StructuredBlock& block,
                                      const MetricField&,
                                      Index3 index,
                                      const QuantityContext& context) {
                                       return thermodynamic_sound_speed(
                                           cell_temperature_state(block, index, context),
                                           context.gas,
                                           context.reference,
                                           context.floors,
                                           block.cell_dimension());
                                   }));
    result.register_quantity(
        field(make_descriptor("mach", "1", QuantityScale::Dimensionless),
              [](const StructuredBlock& block,
                 const MetricField&,
                 Index3 index,
                 const QuantityContext& context) {
                  const auto state = cell_temperature_state(block, index, context);
                  const Real speed
                      = std::sqrt(state[temperature_velocity_x] * state[temperature_velocity_x]
                                  + state[temperature_velocity_y] * state[temperature_velocity_y]
                                  + state[temperature_velocity_z] * state[temperature_velocity_z]);
                  return speed
                      / thermodynamic_sound_speed(state,
                                                  context.gas,
                                                  context.reference,
                                                  context.floors,
                                                  block.cell_dimension());
              }));
    result.register_quantity(
        field(make_descriptor("total_enthalpy", "m^2/s^2", QuantityScale::SpecificEnergy),
              [](const StructuredBlock& block,
                 const MetricField&,
                 Index3 index,
                 const QuantityContext& context) {
                  return thermodynamic_total_enthalpy(cell_temperature_state(block, index, context),
                                                      context.gas,
                                                      context.reference,
                                                      context.floors,
                                                      block.cell_dimension());
              }));
    result.register_quantity(
        field(make_descriptor("entropy_proxy", "1", QuantityScale::Dimensionless),
              [](const StructuredBlock& block,
                 const MetricField&,
                 Index3 index,
                 const QuantityContext& context) {
                  const auto state = cell_pressure_state(block, index, context);
                  return state[pressure] / std::pow(state[primitive_density], context.gas.gamma());
              }));
    result.register_quantity(field(make_descriptor("viscosity", "Pa*s", QuantityScale::Viscosity),
                                   [](const StructuredBlock& block,
                                      const MetricField&,
                                      Index3 index,
                                      const QuantityContext& context) {
                                       return context.transport.viscosity(cell_temperature_state(
                                           block, index, context)[temperature_value]);
                                   }));
    auto jacobian = make_descriptor("jacobian", "m^d", QuantityScale::LengthPower);
    jacobian.length_power = -1;
    result.register_quantity(
        field(std::move(jacobian),
              [](const StructuredBlock&,
                 const MetricField& metric,
                 Index3 index,
                 const QuantityContext&) { return metric.jacobian()(index.i, index.j, index.k); }));
    return result;
}

void FieldQuantityRegistry::register_quantity(std::shared_ptr<const IFieldQuantity> quantity)
{
    if (!quantity) throw std::invalid_argument("field quantity is null");
    quantity->descriptor().validate();
    const auto name = quantity->descriptor().name;
    if (!quantities_.emplace(name, std::move(quantity)).second) {
        throw std::invalid_argument("duplicate field quantity: " + name);
    }
}

void FieldQuantityRegistry::register_turbulence_field(
    const TurbulenceFieldDescriptor& model_descriptor)
{
    model_descriptor.validate();
    QuantityDescriptor descriptor;
    descriptor.name = model_descriptor.name;
    switch (model_descriptor.scale) {
    case TurbulenceFieldScale::Dimensionless:
        descriptor.scale = QuantityScale::Dimensionless;
        break;
    case TurbulenceFieldScale::VelocitySquared:
        descriptor.scale = QuantityScale::SpecificEnergy;
        descriptor.dimensional_unit = "m2/s2";
        break;
    case TurbulenceFieldScale::KinematicViscosity:
        descriptor.scale = QuantityScale::KinematicViscosity;
        descriptor.dimensional_unit = "m2/s";
        break;
    case TurbulenceFieldScale::DynamicViscosity:
        descriptor.scale = QuantityScale::Viscosity;
        descriptor.dimensional_unit = "Pa s";
        break;
    case TurbulenceFieldScale::InverseTime:
        descriptor.scale = QuantityScale::InverseTime;
        descriptor.dimensional_unit = "1/s";
        break;
    case TurbulenceFieldScale::Dissipation:
        descriptor.scale = QuantityScale::Dissipation;
        descriptor.dimensional_unit = "m2/s3";
        break;
    case TurbulenceFieldScale::Length:
        descriptor.scale = QuantityScale::LengthPower;
        descriptor.length_power = 1;
        descriptor.dimensional_unit = "m";
        break;
    }
    const auto name = model_descriptor.name;
    register_quantity(field(std::move(descriptor),
                            [name](const StructuredBlock& block,
                                   const MetricField&,
                                   Index3 index,
                                   const QuantityContext&) {
                                if (!block.turbulence.contains(name)) {
                                    throw PhysicsError("block is missing turbulence field: "
                                                       + name);
                                }
                                return block.turbulence.at(index, name);
                            }));
}

bool FieldQuantityRegistry::contains(const std::string& name) const noexcept
{
    return quantities_.find(name) != quantities_.end();
}

const QuantityDescriptor& FieldQuantityRegistry::descriptor(const std::string& name) const
{
    const auto iterator = quantities_.find(name);
    if (iterator == quantities_.end()) {
        throw std::invalid_argument("unknown field quantity: " + name);
    }
    return iterator->second->descriptor();
}

void FieldQuantityRegistry::validate_dependencies(const std::string& root) const
{
    std::set<std::string> active;
    std::set<std::string> finished;
    std::function<void(const std::string&)> visit = [&](const std::string& name) {
        if (finished.find(name) != finished.end()) return;
        if (!active.insert(name).second) {
            throw std::invalid_argument("cyclic field quantity dependency at " + name);
        }
        const auto iterator = quantities_.find(name);
        if (iterator == quantities_.end()) {
            throw std::invalid_argument("unknown field quantity dependency: " + name);
        }
        for (const auto& dependency : iterator->second->descriptor().dependencies) {
            visit(dependency);
        }
        active.erase(name);
        finished.insert(name);
    };
    visit(root);
}

void FieldQuantityRegistry::validate_selection(const std::vector<std::string>& names) const
{
    std::set<std::string> unique;
    for (const auto& name : names) {
        if (!unique.insert(name).second) {
            throw std::invalid_argument("field quantity was selected twice: " + name);
        }
        validate_dependencies(name);
    }
}

QuantityField FieldQuantityRegistry::evaluate(const std::string& name,
                                              const StructuredBlock& block,
                                              const MetricField& metric,
                                              const QuantityContext& context) const
{
    validate_dependencies(name);
    const auto& quantity = *quantities_.at(name);
    QuantityField result;
    result.descriptor = quantity.descriptor();
    result.extent = block.cell_extent();
    result.values.resize(result.extent.size());
    const Real scale = quantity_scale_factor(result.descriptor, context, block.cell_dimension());
    std::size_t offset = 0;
    for (int k = 0; k < result.extent.nk; ++k) {
        for (int j = 0; j < result.extent.nj; ++j) {
            for (int i = 0; i < result.extent.ni; ++i) {
                const Real value
                    = quantity.evaluate_cell(block, metric, {i, j, k}, context) * scale;
                if (!std::isfinite(value)) {
                    throw PhysicsError("field quantity produced a non-finite value: " + name);
                }
                result.values[offset++] = value;
            }
        }
    }
    return result;
}

StatisticRegistry StatisticRegistry::create_builtin()
{
    StatisticRegistry result;
    const std::array<std::tuple<const char*, int, QuantityScale>, 5> items {{
        {"total_mass", density, QuantityScale::Density},
        {"total_momentum_x", momentum_x, QuantityScale::Momentum},
        {"total_momentum_y", momentum_y, QuantityScale::Momentum},
        {"total_momentum_z", momentum_z, QuantityScale::Momentum},
        {"total_energy", total_energy, QuantityScale::Energy},
    }};
    for (const auto& item : items) {
        auto information = make_descriptor(std::get<0>(item), "integral", std::get<2>(item));
        result.register_quantity(
            statistic(std::move(information),
                      [component = std::get<1>(item)](const StatisticContext& context) {
                          return weighted_conservative_integral(context, component);
                      }));
    }
    return result;
}

void StatisticRegistry::register_quantity(std::shared_ptr<const IStatisticQuantity> quantity)
{
    if (!quantity) throw std::invalid_argument("statistic quantity is null");
    quantity->descriptor().validate();
    const auto name = quantity->descriptor().name;
    if (!quantities_.emplace(name, std::move(quantity)).second) {
        throw std::invalid_argument("duplicate statistic quantity: " + name);
    }
}

bool StatisticRegistry::contains(const std::string& name) const noexcept
{
    return quantities_.find(name) != quantities_.end();
}

Real StatisticRegistry::evaluate(const std::string& name, const StatisticContext& context) const
{
    const auto iterator = quantities_.find(name);
    if (iterator == quantities_.end()) {
        throw std::invalid_argument("unknown statistic quantity: " + name);
    }
    const Real value = iterator->second->evaluate(context);
    if (!std::isfinite(value)) {
        throw PhysicsError("statistic produced a non-finite value: " + name);
    }
    if (!context.quantities.dimensional) return value;
    if (context.partition.zones().empty()) {
        throw std::invalid_argument("statistic partition has no source zones");
    }
    const int dimension = context.partition.zones().front().cell_dimension;
    for (const auto& zone : context.partition.zones()) {
        if (zone.cell_dimension != dimension) {
            throw std::invalid_argument("statistic source zones use mixed cell dimensions");
        }
    }
    const Real field_scale
        = quantity_scale_factor(iterator->second->descriptor(), context.quantities, dimension);
    return value * field_scale
        * std::pow(context.quantities.reference.length(),
                   iterator->second->descriptor().integration_length_power < 0
                       ? dimension
                       : iterator->second->descriptor().integration_length_power);
}

void StatisticRegistry::validate_selection(const std::vector<std::string>& names) const
{
    std::set<std::string> unique;
    for (const auto& name : names) {
        if (!unique.insert(name).second) {
            throw std::invalid_argument("statistic was selected twice: " + name);
        }
        if (!contains(name)) {
            throw std::invalid_argument("unknown statistic quantity: " + name);
        }
    }
}

std::vector<std::string> xz_plane_statistic_names(const std::vector<int>& cell_j_indices)
{
    std::vector<std::string> result;
    result.reserve(2 * cell_j_indices.size());
    for (const int index : cell_j_indices) {
        if (index < 0) {
            throw std::invalid_argument("x-z plane cell-j index must be non-negative");
        }
        result.push_back("xz_mean_u_j" + std::to_string(index));
        result.push_back("xz_mass_flow_x_j" + std::to_string(index));
    }
    return result;
}

void validate_xz_plane_statistics(const std::vector<int>& cell_j_indices,
                                  const StructuredPartitionPlan& partition)
{
    if (partition.zones().empty()) {
        throw std::invalid_argument("x-z plane statistic partition has no zones");
    }
    std::set<int> unique;
    for (const int target_j : cell_j_indices) {
        if (target_j < 0 || !unique.insert(target_j).second) {
            throw std::invalid_argument("x-z plane cell-j indices must be non-negative and unique");
        }
        for (const auto& zone : partition.zones()) {
            if (zone.cell_dimension != 3) {
                throw std::invalid_argument(
                    "x-z plane statistics require three-dimensional source zones");
            }
            if (target_j >= zone.cell_extent.nj) {
                throw std::invalid_argument("x-z plane cell-j index lies outside source zone "
                                            + zone.name);
            }
        }
    }
}

void register_xz_plane_statistics(StatisticRegistry& registry,
                                  const std::vector<int>& cell_j_indices)
{
    const auto accumulate = [](const StatisticContext& context, int target_j, bool mass_flow) {
        validate_xz_plane_statistics({target_j}, context.partition);
        std::unordered_map<BlockId, const PartitionLeaf*> leaves;
        for (const auto& leaf : context.partition.leaves()) {
            leaves.emplace(leaf.block, &leaf);
        }
        Real local_area = 0.0;
        Real local_integral = 0.0;
        Real local_y_min = std::numeric_limits<Real>::infinity();
        Real local_y_max = -std::numeric_limits<Real>::infinity();
        for (const auto& block : context.local_blocks.blocks()) {
            const auto leaf_iterator = leaves.find(block.id());
            const auto metric_iterator = context.metrics.find(block.id());
            if (leaf_iterator == leaves.end() || metric_iterator == context.metrics.end()) {
                throw std::invalid_argument(
                    "x-z plane statistic is missing partition or metric data");
            }
            const auto& leaf = *leaf_iterator->second;
            if (target_j < leaf.cells.begin.j || target_j >= leaf.cells.end.j) {
                continue;
            }
            const int local_j = target_j - leaf.cells.begin.j;
            const auto& metric = metric_iterator->second;
            const auto extent = block.cell_extent();
            for (int k = 0; k < extent.nk; ++k) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Real area = 0.5
                        * (metric.j_faces().area(i, local_j, k)
                           + metric.j_faces().area(i, local_j + 1, k));
                    const Real rho = block.flow.conservative(i, local_j, k, density);
                    const Real rho_u = block.flow.conservative(i, local_j, k, momentum_x);
                    const Real value = mass_flow ? rho_u : rho_u / rho;
                    const Real y = metric.cell_coordinates().y(i, local_j, k);
                    if (!std::isfinite(area) || area <= 0.0 || !std::isfinite(rho)
                        || rho <= context.quantities.floors.density || !std::isfinite(value)
                        || !std::isfinite(y)) {
                        throw PhysicsError("x-z plane statistic encountered invalid data");
                    }
                    local_area += area;
                    local_integral += area * value;
                    local_y_min = std::min(local_y_min, y);
                    local_y_max = std::max(local_y_max, y);
                }
            }
        }
        const Real area = context.mpi.sum(local_area);
        const Real integral = context.mpi.sum(local_integral);
        const Real y_min = context.mpi.min(local_y_min);
        const Real y_max = context.mpi.max(local_y_max);
        if (!std::isfinite(area) || area <= 0.0 || !std::isfinite(integral) || !std::isfinite(y_min)
            || !std::isfinite(y_max)) {
            throw PhysicsError("x-z plane statistic has no finite global support");
        }
        const Real tolerance = 1.0e-10 * (1.0 + std::max(std::abs(y_min), std::abs(y_max)));
        if (y_max - y_min > tolerance) {
            throw PhysicsError("requested cell-j layer is not a geometrically planar x-z section");
        }
        return mass_flow ? integral : integral / area;
    };

    for (const int target_j : cell_j_indices) {
        const auto names = xz_plane_statistic_names({target_j});
        auto mean = make_descriptor(names[0], "m/s", QuantityScale::Velocity);
        mean.integration_length_power = 0;
        registry.register_quantity(
            statistic(std::move(mean), [accumulate, target_j](const StatisticContext& context) {
                return accumulate(context, target_j, false);
            }));
        auto flow = make_descriptor(names[1], "kg/s", QuantityScale::Momentum);
        flow.integration_length_power = 2;
        registry.register_quantity(
            statistic(std::move(flow), [accumulate, target_j](const StatisticContext& context) {
                return accumulate(context, target_j, true);
            }));
    }
}

std::vector<std::string> yz_plane_statistic_names(std::size_t plane_count)
{
    std::vector<std::string> result;
    result.reserve(2 * plane_count);
    for (std::size_t plane = 0; plane < plane_count; ++plane) {
        const auto suffix = "_plane" + std::to_string(plane);
        result.push_back("yz_mean_u" + suffix);
        result.push_back("yz_mass_flow_x" + suffix);
    }
    return result;
}

void validate_yz_plane_statistics(const std::vector<Real>& target_x_coordinates,
                                  const StructuredPartitionPlan& partition)
{
    if (target_x_coordinates.empty()) {
        throw std::invalid_argument(
            "y-z plane statistics require at least one target x coordinate");
    }
    std::set<Real> unique;
    for (const Real target : target_x_coordinates) {
        if (!std::isfinite(target) || !unique.insert(target).second) {
            throw std::invalid_argument("y-z plane target x coordinates must be finite and unique");
        }
    }
    if (partition.zones().empty()) {
        throw std::invalid_argument("y-z plane statistic partition has no zones");
    }
    for (const auto& zone : partition.zones()) {
        if (zone.cell_dimension != 3) {
            throw std::invalid_argument(
                "y-z plane statistics require three-dimensional source zones");
        }
    }
}

void register_yz_plane_statistics(StatisticRegistry& registry,
                                  const std::vector<Real>& target_x_coordinates)
{
    const auto planar_x = [](const StructuredBlock& block, const MetricField& metric, int cell_i) {
        const auto extent = block.cell_extent();
        Real minimum = std::numeric_limits<Real>::infinity();
        Real maximum = -std::numeric_limits<Real>::infinity();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                const Real x = metric.cell_coordinates().x(cell_i, j, k);
                if (!std::isfinite(x)) {
                    throw PhysicsError("y-z plane statistic encountered a non-finite coordinate");
                }
                minimum = std::min(minimum, x);
                maximum = std::max(maximum, x);
            }
        }
        const Real tolerance = 1.0e-10 * (1.0 + std::max(std::abs(minimum), std::abs(maximum)));
        if (maximum - minimum > tolerance) {
            throw PhysicsError(
                "y-z plane statistics require geometrically planar constant-x sections");
        }
        return 0.5 * (minimum + maximum);
    };

    const auto accumulate = [planar_x](
                                const StatisticContext& context, Real target_x, bool mass_flow) {
        Real local_selected_x = std::numeric_limits<Real>::infinity();
        for (const auto& block : context.local_blocks.blocks()) {
            if (block.cell_dimension() != 3) {
                throw PhysicsError("y-z plane statistics require a 3D mesh");
            }
            const auto metric_iterator = context.metrics.find(block.id());
            if (metric_iterator == context.metrics.end()) {
                throw std::invalid_argument("y-z plane statistic is missing metric data");
            }
            const auto extent = block.cell_extent();
            for (int i = 0; i < extent.ni; ++i) {
                const Real x = planar_x(block, metric_iterator->second, i);
                const Real tolerance = 1.0e-12 * (1.0 + std::max(std::abs(x), std::abs(target_x)));
                if (x + tolerance >= target_x) {
                    local_selected_x = std::min(local_selected_x, x);
                }
            }
        }
        const Real selected_x = context.mpi.min(local_selected_x);
        if (!std::isfinite(selected_x)) {
            throw PhysicsError(
                "y-z plane target has no cell-centre section on its positive-x side");
        }

        Real local_area = 0.0;
        Real local_integral = 0.0;
        const Real selection_tolerance = 1.0e-10 * (1.0 + std::abs(selected_x));
        for (const auto& block : context.local_blocks.blocks()) {
            const auto& metric = context.metrics.at(block.id());
            const auto extent = block.cell_extent();
            for (int i = 0; i < extent.ni; ++i) {
                const Real x = planar_x(block, metric, i);
                if (std::abs(x - selected_x) > selection_tolerance) continue;
                for (int k = 0; k < extent.nk; ++k) {
                    for (int j = 0; j < extent.nj; ++j) {
                        const Real area = 0.5
                            * (metric.i_faces().area(i, j, k) + metric.i_faces().area(i + 1, j, k));
                        const Real rho = block.flow.conservative(i, j, k, density);
                        const Real rho_u = block.flow.conservative(i, j, k, momentum_x);
                        const Real value = mass_flow ? rho_u : rho_u / rho;
                        if (!std::isfinite(area) || area <= 0.0 || !std::isfinite(rho)
                            || rho <= context.quantities.floors.density || !std::isfinite(value)) {
                            throw PhysicsError("y-z plane statistic encountered invalid flow data");
                        }
                        local_area += area;
                        local_integral += area * value;
                    }
                }
            }
        }
        const Real area = context.mpi.sum(local_area);
        const Real integral = context.mpi.sum(local_integral);
        if (!std::isfinite(area) || area <= 0.0 || !std::isfinite(integral)) {
            throw PhysicsError("y-z plane statistic has no finite global support");
        }
        return mass_flow ? integral : integral / area;
    };

    const auto names = yz_plane_statistic_names(target_x_coordinates.size());
    for (std::size_t plane = 0; plane < target_x_coordinates.size(); ++plane) {
        auto mean = make_descriptor(names[2 * plane], "m/s", QuantityScale::Velocity);
        mean.integration_length_power = 0;
        registry.register_quantity(statistic(
            std::move(mean),
            [accumulate, target_x = target_x_coordinates[plane]](const StatisticContext& context) {
                return accumulate(context, target_x, false);
            }));
        auto flow = make_descriptor(names[2 * plane + 1], "kg/s", QuantityScale::Momentum);
        flow.integration_length_power = 2;
        registry.register_quantity(statistic(
            std::move(flow),
            [accumulate, target_x = target_x_coordinates[plane]](const StatisticContext& context) {
                return accumulate(context, target_x, true);
            }));
    }
}

std::vector<std::string> channel_wall_statistic_names()
{
    return {
        "channel_wall_shear_lower",
        "channel_wall_shear_upper",
        "channel_wall_shear_mean",
        "channel_friction_velocity",
        "channel_re_tau",
    };
}

void register_channel_wall_statistics(StatisticRegistry& registry,
                                      const std::string& lower_patch,
                                      const std::string& upper_patch,
                                      Real half_height)
{
    if (lower_patch.empty() || upper_patch.empty() || lower_patch == upper_patch
        || !std::isfinite(half_height) || half_height <= 0.0) {
        throw std::invalid_argument("channel-wall statistic definition is invalid");
    }
    const auto average = [lower_patch, upper_patch](const StatisticContext& context) {
        return combined_channel_walls(channel_wall_average(context, lower_patch, Side::Lower),
                                      channel_wall_average(context, upper_patch, Side::Upper));
    };
    const auto names = channel_wall_statistic_names();
    auto lower = make_descriptor(names[0], "Pa", QuantityScale::Pressure);
    lower.integration_length_power = 0;
    registry.register_quantity(
        statistic(std::move(lower), [lower_patch](const StatisticContext& context) {
            return channel_wall_average(context, lower_patch, Side::Lower).shear;
        }));
    auto upper = make_descriptor(names[1], "Pa", QuantityScale::Pressure);
    upper.integration_length_power = 0;
    registry.register_quantity(
        statistic(std::move(upper), [upper_patch](const StatisticContext& context) {
            return channel_wall_average(context, upper_patch, Side::Upper).shear;
        }));
    auto mean = make_descriptor(names[2], "Pa", QuantityScale::Pressure);
    mean.integration_length_power = 0;
    registry.register_quantity(
        statistic(std::move(mean),
                  [average](const StatisticContext& context) { return average(context).shear; }));
    auto friction = make_descriptor(names[3], "m/s", QuantityScale::Velocity);
    friction.integration_length_power = 0;
    registry.register_quantity(
        statistic(std::move(friction), [average](const StatisticContext& context) {
            const auto wall = average(context);
            return std::sqrt(wall.shear / wall.density);
        }));
    auto re_tau = make_descriptor(names[4], "1", QuantityScale::Dimensionless);
    re_tau.integration_length_power = 0;
    registry.register_quantity(
        statistic(std::move(re_tau), [average, half_height](const StatisticContext& context) {
            const auto wall = average(context);
            const Real friction_velocity = std::sqrt(wall.shear / wall.density);
            return wall.density * friction_velocity * half_height
                * context.quantities.reference.reynolds() / wall.viscosity_ratio;
        }));
}

} // namespace wcns
