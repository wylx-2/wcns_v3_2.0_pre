#pragma once

#include <wcns/mesh/structured_block.hpp>
#include <wcns/physics/double_mach_reflection.hpp>
#include <wcns/physics/thermodynamics.hpp>
#include <wcns/solver/euler.hpp>

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_map>

namespace wcns {

struct FarfieldPointVortex {
    Real lift_coefficient = 0.0;
    std::array<Real, 2> center {{0.25, 0.0}};
    Real chord = 1.0;

    void validate() const;
};

struct BoundaryData {
    std::optional<TemperaturePrimitiveState> target_state;
    std::array<Real, 3> wall_velocity {{0.0, 0.0, 0.0}};
    std::optional<Real> wall_temperature;
    std::optional<DoubleMachReflection> double_mach_reflection;
    std::optional<FarfieldPointVortex> farfield_point_vortex;
    bool compression_ramp_inlet = false;

    void validate(BoundaryType type, int dimension) const;
};

[[nodiscard]] TemperaturePrimitiveState
farfield_target_at(const BoundaryData& data,
                   std::array<Real, 3> face_coordinates,
                   int dimension);

using BoundaryDataMap = std::unordered_map<std::string, BoundaryData>;

struct PhysicalGhostFillResult {
    std::uint64_t version = 0;
    std::size_t state_count = 0;
};

class PhysicalGhostStateOperator {
public:
    [[nodiscard]] static PhysicalGhostFillResult fill(StructuredBlock& block,
                                                      const BoundaryDataMap& boundary_data,
                                                      const GasModel& gas,
                                                      const ReferenceScales& reference,
                                                      const NumericalFloors& floors,
                                                      std::uint64_t version,
                                                      Real time = 0.0);
};

[[nodiscard]] std::array<Real, 3>
boundary_face_coordinates(const StructuredBlock& block, const BoundaryPatch& patch, Index3 face);

struct InviscidBoundaryOptions {
    bool strong_boundary_face_state = true;

    [[nodiscard]] std::string summary() const;
    [[nodiscard]] std::string restart_signature() const;
};

[[nodiscard]] PressurePrimitiveState
apply_inviscid_boundary_face_state(const BoundaryPatch& patch,
                                   const PressurePrimitiveState& interior_trace,
                                   const PressurePrimitiveState& reconstructed_exterior_trace,
                                   Normal3 outward_unit_normal,
                                   const BoundaryData& data,
                                   const InviscidBoundaryOptions& options,
                                   const GasModel& gas,
                                   const ReferenceScales& reference,
                                   const NumericalFloors& floors,
                                   int dimension,
                                   std::array<Real, 3> face_coordinates = {{0.0, 0.0, 0.0}},
                                   Real time = 0.0);

void update_temperature_primitive_interior(StructuredBlock& block,
                                           const GasModel& gas,
                                           const ReferenceScales& reference,
                                           const NumericalFloors& floors);

void update_temperature_primitive_cell(StructuredBlock& block,
                                       Index3 index,
                                       const GasModel& gas,
                                       const ReferenceScales& reference,
                                       const NumericalFloors& floors);

} // namespace wcns
