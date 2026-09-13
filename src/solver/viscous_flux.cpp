#include <wcns/solver/viscous_flux.hpp>

#include <cmath>
#include <stdexcept>

namespace wcns {
namespace {

void require_finite_gradients(const PrimitiveGradients& gradients, int dimension)
{
    for (int variable = 0; variable < viscous_primitive_components; ++variable) {
        for (int direction = 0; direction < 3; ++direction) {
            const Real value = gradients[static_cast<std::size_t>(variable)]
                                        [static_cast<std::size_t>(direction)];
            if (!std::isfinite(value)) {
                throw PhysicsError("viscous primitive gradient is non-finite");
            }
            if (dimension == 2
                && (variable == static_cast<int>(ViscousPrimitive::VelocityZ) || direction == 2)
                && value != 0.0) {
                throw PhysicsError("two-dimensional viscous z gradients must be zero");
            }
        }
    }
}

} // namespace

ViscousCartesianFlux compute_viscous_cartesian_flux(
    const ViscousFaceTrace& trace,
    const TransportModel& transport,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    const TurbulenceViscousContribution& turbulence,
    int dimension)
{
    // The conversion is also the authoritative finite/positivity/dimension check.
    (void)pressure_primitive(trace.state, gas, reference, floors, dimension);
    require_finite_gradients(trace.gradients, dimension);
    turbulence.validate(dimension);

    const auto& gu = trace.gradients[static_cast<int>(ViscousPrimitive::VelocityX)];
    const auto& gv = trace.gradients[static_cast<int>(ViscousPrimitive::VelocityY)];
    const auto& gw = trace.gradients[static_cast<int>(ViscousPrimitive::VelocityZ)];
    const auto& gt = trace.gradients[static_cast<int>(ViscousPrimitive::Temperature)];
    const Real mu = transport.viscosity(trace.state[temperature_value]);
    const Real chi = transport.thermal_coefficient(trace.state[temperature_value], gas, reference);
    const Real theta = gu[0] + gv[1] + gw[2];
    const Real isotropic = (2.0 / 3.0) * mu * theta;
    const Real tau_xx = 2.0 * mu * gu[0] - isotropic + turbulence.stress.xx;
    const Real tau_yy = 2.0 * mu * gv[1] - isotropic + turbulence.stress.yy;
    const Real tau_zz = 2.0 * mu * gw[2] - isotropic + turbulence.stress.zz;
    const Real tau_xy = mu * (gu[1] + gv[0]) + turbulence.stress.xy;
    const Real tau_xz = mu * (gu[2] + gw[0]) + turbulence.stress.xz;
    const Real tau_yz = mu * (gv[2] + gw[1]) + turbulence.stress.yz;
    const Real u = trace.state[temperature_velocity_x];
    const Real v = trace.state[temperature_velocity_y];
    const Real w = trace.state[temperature_velocity_z];

    ViscousCartesianFlux result;
    result.viscosity = mu + turbulence.eddy_viscosity;
    result.thermal_coefficient = chi;
    result.x = {{0.0,
                 tau_xx,
                 tau_xy,
                 tau_xz,
                 u * tau_xx + v * tau_xy + w * tau_xz + chi * gt[0]
                     + turbulence.energy_heat_flux[0]}};
    result.y = {{0.0,
                 tau_xy,
                 tau_yy,
                 tau_yz,
                 u * tau_xy + v * tau_yy + w * tau_yz + chi * gt[1]
                     + turbulence.energy_heat_flux[1]}};
    result.z = {{0.0,
                 tau_xz,
                 tau_yz,
                 tau_zz,
                 u * tau_xz + v * tau_yz + w * tau_zz + chi * gt[2]
                     + turbulence.energy_heat_flux[2]}};
    for (const auto* flux : {&result.x, &result.y, &result.z}) {
        for (const Real value : *flux) {
            if (!std::isfinite(value)) {
                throw PhysicsError("viscous Cartesian flux is non-finite");
            }
        }
    }
    if (dimension == 2) {
        result.x[momentum_z] = 0.0;
        result.y[momentum_z] = 0.0;
        result.z.fill(0.0);
    }
    return result;
}

ViscousCartesianFlux compute_viscous_cartesian_flux(const ViscousFaceTrace& trace,
                                                    const TransportModel& transport,
                                                    const GasModel& gas,
                                                    const ReferenceScales& reference,
                                                    const NumericalFloors& floors,
                                                    int dimension)
{
    return compute_viscous_cartesian_flux(
        trace, transport, gas, reference, floors, TurbulenceViscousContribution {}, dimension);
}

} // namespace wcns
