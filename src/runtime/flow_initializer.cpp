#include <wcns/runtime/flow_initializer.hpp>
#include <wcns/physics/double_mach_reflection.hpp>

#include <algorithm>
#include <cmath>
#include <string>

namespace wcns {
namespace {

constexpr Real pi = 3.141592653589793238462643383279502884;

TemperaturePrimitiveState from_temperature(Real rho,
                                           Real u,
                                           Real v,
                                           Real w,
                                           Real temperature,
                                           const GasModel& gas,
                                           const ReferenceScales& reference,
                                           const NumericalFloors& floors,
                                           int dimension)
{
    TemperaturePrimitiveState state {rho, u, v, w, temperature};
    static_cast<void>(pressure_primitive(state, gas, reference, floors, dimension));
    return state;
}

TemperaturePrimitiveState from_pressure(Real rho,
                                        Real u,
                                        Real v,
                                        Real w,
                                        Real pressure,
                                        const GasModel& gas,
                                        const ReferenceScales& reference,
                                        const NumericalFloors& floors,
                                        int dimension)
{
    return temperature_primitive({rho, u, v, w, pressure}, gas, reference, floors, dimension);
}

TemperaturePrimitiveState constant_pressure_temperature(Real pressure,
                                                        Real u,
                                                        Real v,
                                                        Real w,
                                                        Real temperature,
                                                        const GasModel& gas,
                                                        const ReferenceScales& reference,
                                                        const NumericalFloors& floors,
                                                        int dimension)
{
    const Real rho = gas.gamma() * reference.mach() * reference.mach() * pressure / temperature;
    return from_temperature(rho, u, v, w, temperature, gas, reference, floors, dimension);
}

Real normalized_wall_coordinate(const InitialConditionConfig& config, Real y)
{
    const Real y0 = config.parameter("y0", 0.0);
    const Real y1 = config.parameter("y1", 1.0);
    return (y - y0) / (y1 - y0);
}

Real analytic_pressure(const InitialConditionConfig& config,
                       const GasModel& gas,
                       const ReferenceScales& reference)
{
    return config.parameter("pressure", 1.0 / (gas.gamma() * reference.mach() * reference.mach()));
}

Real reichardt_velocity_plus(Real y_plus)
{
    constexpr Real kappa = 0.41;
    constexpr Real wake_amplitude = 7.8;
    constexpr Real outer_scale = 11.0;
    return std::log1p(kappa * y_plus) / kappa
        + wake_amplitude
        * (1.0 - std::exp(-y_plus / outer_scale)
           - (y_plus / outer_scale) * std::exp(-y_plus / 3.0));
}

Real reichardt_velocity_plus_derivative(Real y_plus)
{
    constexpr Real kappa = 0.41;
    constexpr Real wake_amplitude = 7.8;
    constexpr Real outer_scale = 11.0;
    return 1.0 / (1.0 + kappa * y_plus)
        + wake_amplitude
        * (std::exp(-y_plus / outer_scale) / outer_scale - std::exp(-y_plus / 3.0) / outer_scale
           + y_plus * std::exp(-y_plus / 3.0) / (3.0 * outer_scale));
}

TemperaturePrimitiveState turbulent_channel_state(const InitialConditionConfig& config,
                                                  std::array<Real, 3> coordinates,
                                                  const GasModel& gas,
                                                  const ReferenceScales& reference,
                                                  const NumericalFloors& floors,
                                                  int dimension)
{
    if (dimension != 3) {
        throw FlowInitializationError(
            "turbulent-channel initial condition requires a three-dimensional mesh");
    }
    const Real y0 = config.parameter("y0", -1.0);
    const Real y1 = config.parameter("y1", 1.0);
    const Real half_height = 0.5 * (y1 - y0);
    const Real center = 0.5 * (y0 + y1);
    const Real centered_y = (coordinates[1] - center) / half_height;
    Real wall_fraction = 1.0 - std::abs(centered_y);
    constexpr Real coordinate_tolerance = 1.0e-12;
    if (wall_fraction < -coordinate_tolerance || wall_fraction > 1.0 + coordinate_tolerance) {
        throw FlowInitializationError(
            "turbulent-channel cell lies outside the configured wall interval");
    }
    wall_fraction = std::max(Real {0.0}, std::min(Real {1.0}, wall_fraction));

    const Real re_tau = config.parameter("re_tau", 180.0);
    const Real bulk_velocity = config.parameter("bulk_velocity", 1.0);
    const Real bulk_velocity_plus = config.parameter("bulk_velocity_plus", 15.5);
    const Real y_plus = re_tau * wall_fraction;
    // The small eta^8 correction makes the mirrored wall-law profile have zero
    // first derivative at the channel centre without disturbing the near-wall law.
    const Real center_slope = re_tau * reichardt_velocity_plus_derivative(re_tau);
    const Real mean_plus
        = reichardt_velocity_plus(y_plus) - center_slope * std::pow(wall_fraction, 8) / 8.0;
    const Real mean_u = bulk_velocity * mean_plus / bulk_velocity_plus;

    const Real period_x = config.parameter("period_x", 2.0 * pi);
    const Real period_z = config.parameter("period_z", pi);
    const Real phase_x = 2.0 * pi * (coordinates[0] - config.parameter("x0", 0.0)) / period_x;
    const Real phase_z = 2.0 * pi * (coordinates[2] - config.parameter("z0", 0.0)) / period_z;
    const Real envelope = std::pow(1.0 - centered_y * centered_y, 2);
    const Real amplitude
        = config.parameter("perturbation_amplitude", 0.05) * bulk_velocity * envelope;

    // Integer low modes are periodic on the requested box, resolved by the
    // 36x36 x-z grid, have zero plane mean and vanish at both no-slip walls.
    const Real u = mean_u
        + amplitude
            * (std::sin(2.0 * phase_x) * std::cos(phase_z)
               + 0.5 * std::sin(3.0 * phase_x + 2.0 * phase_z));
    const Real v = amplitude
        * (std::sin(phase_x) * std::sin(phase_z) + 0.5 * std::cos(2.0 * phase_x - phase_z));
    const Real w = amplitude
        * (std::cos(phase_x) * std::sin(2.0 * phase_z) - 0.5 * std::sin(3.0 * phase_x - phase_z));
    return from_temperature(config.parameter("rho", 1.0),
                            u,
                            v,
                            w,
                            config.parameter("temperature", 1.0),
                            gas,
                            reference,
                            floors,
                            dimension);
}

TemperaturePrimitiveState couette_state(const InitialConditionConfig& config,
                                        Real y,
                                        const GasModel& gas,
                                        const ReferenceScales& reference,
                                        const NumericalFloors& floors,
                                        int dimension)
{
    const Real eta = normalized_wall_coordinate(config, y);
    const Real lower_u = config.parameter("lower_velocity", 0.0);
    const Real upper_u = config.parameter("upper_velocity", 1.0);
    const Real lower_temperature
        = config.parameter("lower_temperature", config.parameter("temperature", 1.0));
    const Real upper_temperature
        = config.parameter("upper_temperature", config.parameter("temperature", 1.0));
    const Real temperature = lower_temperature + (upper_temperature - lower_temperature) * eta
        + config.parameter("temperature_curvature", 0.0) * eta * (1.0 - eta);
    const Real velocity = lower_u + (upper_u - lower_u) * eta
        + config.parameter("velocity_curvature", 0.0) * eta * (1.0 - eta);
    return constant_pressure_temperature(analytic_pressure(config, gas, reference),
                                         velocity,
                                         0.0,
                                         0.0,
                                         temperature,
                                         gas,
                                         reference,
                                         floors,
                                         dimension);
}

TemperaturePrimitiveState poiseuille_state(const InitialConditionConfig& config,
                                           Real y,
                                           const GasModel& gas,
                                           const ReferenceScales& reference,
                                           const NumericalFloors& floors,
                                           int dimension)
{
    const Real eta = normalized_wall_coordinate(config, y);
    const Real centerline_velocity = config.parameter("centerline_velocity", 1.0);
    const Real velocity = 4.0 * centerline_velocity * eta * (1.0 - eta);
    // For u=A*eta*(1-eta), the polynomial below solves
    // T''=-(gamma-1)*Ma^2*Pr*(du/deta)^2 when its configured amplitude is
    // (gamma-1)*Ma^2*Pr*A^2.  A zero amplitude gives the usual isothermal
    // incompressible/low-Mach analytic initial profile.
    const Real thermal_shape
        = eta / 6.0 - eta * eta / 2.0 + 2.0 * eta * eta * eta / 3.0 - eta * eta * eta * eta / 3.0;
    const Real temperature = config.parameter("temperature", 1.0)
        + config.parameter("temperature_curvature", 0.0) * thermal_shape;
    return constant_pressure_temperature(analytic_pressure(config, gas, reference),
                                         velocity,
                                         0.0,
                                         0.0,
                                         temperature,
                                         gas,
                                         reference,
                                         floors,
                                         dimension);
}

TemperaturePrimitiveState linear_conduction_state(const InitialConditionConfig& config,
                                                  Real y,
                                                  const GasModel& gas,
                                                  const ReferenceScales& reference,
                                                  const NumericalFloors& floors,
                                                  int dimension)
{
    const Real eta = normalized_wall_coordinate(config, y);
    const Real lower_temperature = config.parameter("lower_temperature", 1.0);
    const Real upper_temperature = config.parameter("upper_temperature", 2.0);
    const Real temperature = lower_temperature + (upper_temperature - lower_temperature) * eta;
    return constant_pressure_temperature(analytic_pressure(config, gas, reference),
                                         0.0,
                                         0.0,
                                         0.0,
                                         temperature,
                                         gas,
                                         reference,
                                         floors,
                                         dimension);
}

TemperaturePrimitiveState uniform_state(const InitialConditionConfig& config,
                                        const GasModel& gas,
                                        const ReferenceScales& reference,
                                        const NumericalFloors& floors,
                                        int dimension)
{
    const Real rho = config.parameter("rho", 1.0);
    const Real u = config.parameter("u", 0.0);
    const Real v = config.parameter("v", 0.0);
    const Real w = config.parameter("w", 0.0);
    const auto pressure = config.parameters.find("pressure");
    if (pressure != config.parameters.end()) {
        return from_pressure(rho, u, v, w, pressure->second, gas, reference, floors, dimension);
    }
    return from_temperature(
        rho, u, v, w, config.parameter("temperature", 1.0), gas, reference, floors, dimension);
}

TemperaturePrimitiveState sod_state(const InitialConditionConfig& config,
                                    Real x,
                                    const GasModel& gas,
                                    const ReferenceScales& reference,
                                    const NumericalFloors& floors,
                                    int dimension)
{
    const bool left = x < config.parameter("x0", 0.5);
    const std::string prefix = left ? "left_" : "right_";
    return from_pressure(config.parameter(prefix + "rho", left ? 1.0 : 0.125),
                         config.parameter(prefix + "u", 0.0),
                         config.parameter(prefix + "v", 0.0),
                         0.0,
                         config.parameter(prefix + "p", left ? 1.0 : 0.1),
                         gas,
                         reference,
                         floors,
                         dimension);
}

TemperaturePrimitiveState double_mach_reflection_state(const InitialConditionConfig& config,
                                                       Real x,
                                                       Real y,
                                                       const GasModel& gas,
                                                       const ReferenceScales& reference,
                                                       const NumericalFloors& floors,
                                                       int dimension)
{
    const DoubleMachReflection model(config.parameter("x0", 1.0 / 6.0));
    model.validate(gas.gamma(), dimension);
    return temperature_primitive(model.exact_state(x, y, 0.0), gas, reference, floors, dimension);
}

TemperaturePrimitiveState quadrant_state(const InitialConditionConfig& config,
                                         Real x,
                                         Real y,
                                         const GasModel& gas,
                                         const ReferenceScales& reference,
                                         const NumericalFloors& floors,
                                         int dimension)
{
    const bool east = x >= config.parameter("x0", 0.5);
    const bool north = y >= config.parameter("y0", 0.5);
    const char* prefix = east ? (north ? "ne_" : "se_") : (north ? "nw_" : "sw_");
    const bool ne = std::string(prefix) == "ne_";
    const bool nw = std::string(prefix) == "nw_";
    const bool sw = std::string(prefix) == "sw_";
    const Real default_rho = ne ? 1.5 : (nw || !sw ? 0.5323 : 0.138);
    const Real default_u = (nw || sw) ? 1.206 : 0.0;
    const Real default_v = (sw || (!ne && !nw)) ? 1.206 : 0.0;
    const Real default_p = ne ? 1.5 : (sw ? 0.029 : 0.3);
    return from_pressure(config.parameter(std::string(prefix) + "rho", default_rho),
                         config.parameter(std::string(prefix) + "u", default_u),
                         config.parameter(std::string(prefix) + "v", default_v),
                         0.0,
                         config.parameter(std::string(prefix) + "p", default_p),
                         gas,
                         reference,
                         floors,
                         dimension);
}

TemperaturePrimitiveState isentropic_vortex_state(const InitialConditionConfig& config,
                                                  Real x,
                                                  Real y,
                                                  const GasModel& gas,
                                                  const ReferenceScales& reference,
                                                  const NumericalFloors& floors,
                                                  int dimension)
{
    const Real x0 = config.parameter("x0", 5.0);
    const Real y0 = config.parameter("y0", 5.0);
    const Real beta = config.parameter("beta", 5.0);
    const Real period_x = config.parameter("period_x", 0.0);
    const Real period_y = config.parameter("period_y", 0.0);
    const Real dx = period_x > 0.0 ? std::remainder(x - x0, period_x) : x - x0;
    const Real dy = period_y > 0.0 ? std::remainder(y - y0, period_y) : y - y0;
    const Real radius_squared = dx * dx + dy * dy;
    const Real exponential = std::exp(0.5 * (1.0 - radius_squared));
    const Real u = config.parameter("background_u", 1.0) - beta * exponential * dy / (2.0 * pi);
    const Real v = config.parameter("background_v", 1.0) + beta * exponential * dx / (2.0 * pi);
    const Real temperature = 1.0
        - (gas.gamma() - 1.0) * reference.mach() * reference.mach() * beta * beta
            * std::exp(1.0 - radius_squared) / (8.0 * pi * pi);
    if (!std::isfinite(temperature) || temperature <= floors.temperature) {
        throw FlowInitializationError("isentropic vortex produces a non-positive temperature");
    }
    const Real rho = std::pow(temperature, 1.0 / (gas.gamma() - 1.0));
    return from_temperature(rho, u, v, 0.0, temperature, gas, reference, floors, dimension);
}

TemperaturePrimitiveState manufactured_state(const InitialConditionConfig& config,
                                             Real x,
                                             Real y,
                                             Real z,
                                             const GasModel& gas,
                                             const ReferenceScales& reference,
                                             const NumericalFloors& floors,
                                             int dimension)
{
    const Real amplitude = config.parameter("beta", 0.01);
    const Real sx = std::sin(2.0 * pi * x);
    const Real sy = std::sin(2.0 * pi * y);
    const Real sz = dimension == 3 ? std::sin(2.0 * pi * z) : 1.0;
    const Real cx = std::cos(2.0 * pi * x);
    const Real cy = std::cos(2.0 * pi * y);
    const Real cz = dimension == 3 ? std::cos(2.0 * pi * z) : 0.0;
    return from_temperature(1.0 + amplitude * sx * sy * sz,
                            config.parameter("background_u", 0.2) + amplitude * cx * sy,
                            config.parameter("background_v", 0.1) - amplitude * sx * cy,
                            dimension == 3 ? amplitude * sx * sy * cz : 0.0,
                            1.0 + amplitude * cx * cy,
                            gas,
                            reference,
                            floors,
                            dimension);
}

TemperaturePrimitiveState taylor_green_state(const InitialConditionConfig& config,
                                              Real x,
                                              Real y,
                                              Real z,
                                              const GasModel& gas,
                                              const ReferenceScales& reference,
                                              const NumericalFloors& floors,
                                              int dimension)
{
    if (dimension != 3) {
        throw FlowInitializationError("Taylor-Green vortex requires three dimensions");
    }
    const Real amplitude = config.parameter("beta", 0.1);
    const Real phase_x = 2.0 * pi * x;
    const Real phase_y = 2.0 * pi * y;
    const Real phase_z = 2.0 * pi * z;
    const Real u = amplitude * std::sin(phase_x) * std::cos(phase_y) * std::cos(phase_z);
    const Real v = -amplitude * std::cos(phase_x) * std::sin(phase_y) * std::cos(phase_z);
    const Real pressure_pattern
        = (std::cos(2.0 * phase_x) + std::cos(2.0 * phase_y))
        * (std::cos(2.0 * phase_z) + 2.0);
    const Real temperature = config.parameter("temperature", 1.0)
        + gas.gamma() * reference.mach() * reference.mach() * amplitude * amplitude
            * pressure_pattern / 16.0;
    return from_temperature(config.parameter("rho", 1.0),
                            u,
                            v,
                            0.0,
                            temperature,
                            gas,
                            reference,
                            floors,
                            dimension);
}

} // namespace

TemperaturePrimitiveState FlowInitializer::evaluate(const InitialConditionConfig& config,
                                                    std::array<Real, 3> coordinates,
                                                    const GasModel& gas,
                                                    const ReferenceScales& reference,
                                                    const NumericalFloors& floors,
                                                    int dimension)
{
    config.validate(dimension);
    if (config.type == "uniform") {
        return uniform_state(config, gas, reference, floors, dimension);
    }
    if (config.type == "sod_x") {
        return sod_state(config, coordinates[0], gas, reference, floors, dimension);
    }
    if (config.type == "double_mach_reflection") {
        return double_mach_reflection_state(
            config, coordinates[0], coordinates[1], gas, reference, floors, dimension);
    }
    if (config.type == "quadrant_riemann") {
        return quadrant_state(
            config, coordinates[0], coordinates[1], gas, reference, floors, dimension);
    }
    if (config.type == "isentropic_vortex") {
        return isentropic_vortex_state(
            config, coordinates[0], coordinates[1], gas, reference, floors, dimension);
    }
    if (config.type == "couette") {
        return couette_state(config, coordinates[1], gas, reference, floors, dimension);
    }
    if (config.type == "poiseuille") {
        return poiseuille_state(config, coordinates[1], gas, reference, floors, dimension);
    }
    if (config.type == "turbulent_channel") {
        return turbulent_channel_state(config, coordinates, gas, reference, floors, dimension);
    }
    if (config.type == "linear_conduction") {
        return linear_conduction_state(config, coordinates[1], gas, reference, floors, dimension);
    }
    if (config.type == "manufactured_periodic") {
        return manufactured_state(config,
                                  coordinates[0],
                                  coordinates[1],
                                  coordinates[2],
                                  gas,
                                  reference,
                                  floors,
                                  dimension);
    }
    if (config.type == "taylor_green_vortex") {
        return taylor_green_state(config,
                                  coordinates[0],
                                  coordinates[1],
                                  coordinates[2],
                                  gas,
                                  reference,
                                  floors,
                                  dimension);
    }
    throw FlowInitializationError("initial condition registry has no implementation for "
                                  + config.type);
}

void FlowInitializer::initialize_block(StructuredBlock& block,
                                       const MetricField& metrics,
                                       const InitialConditionConfig& config,
                                       const GasModel& gas,
                                       const ReferenceScales& reference,
                                       const NumericalFloors& floors)
{
    if (metrics.dimension() != block.cell_dimension()
        || metrics.cell_coordinates().x.interior_extent() != block.cell_extent()) {
        throw FlowInitializationError("initial condition metric/block metadata mismatch");
    }
    const auto extent = block.cell_extent();
    const auto& coordinates = metrics.cell_coordinates();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const auto primitive = evaluate(config,
                                                {
                                                    coordinates.x(i, j, k),
                                                    coordinates.y(i, j, k),
                                                    coordinates.z(i, j, k),
                                                },
                                                gas,
                                                reference,
                                                floors,
                                                block.cell_dimension());
                store_state(block.flow.conservative,
                            {i, j, k},
                            thermodynamic_conservative(
                                primitive, gas, reference, floors, block.cell_dimension()));
            }
        }
    }
}

void FlowInitializer::initialize_local_blocks(LocalBlockSet& local_blocks,
                                              const BlockMetricMap& metrics,
                                              const InitialConditionConfig& config,
                                              const GasModel& gas,
                                              const ReferenceScales& reference,
                                              const NumericalFloors& floors)
{
    for (auto& block : local_blocks.blocks()) {
        const auto iterator = metrics.find(block.id());
        if (iterator == metrics.end()) {
            throw FlowInitializationError("initial condition is missing metrics for block "
                                          + std::to_string(block.id()));
        }
        initialize_block(block, iterator->second, config, gas, reference, floors);
    }
}

} // namespace wcns
