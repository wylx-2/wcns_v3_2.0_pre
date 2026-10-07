#pragma once

#include <wcns/core/topology_field.hpp>
#include <wcns/runtime/structured_partition.hpp>
#include <wcns/mesh/conservation_weights.hpp>
#include <wcns/solver/inviscid_wcns_solver.hpp>
#include <wcns/solver/physical_boundary.hpp>
#include <wcns/solver/transport_model.hpp>
#include <wcns/solver/turbulence_fields.hpp>

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace wcns {

enum class QuantityScale {
    Dimensionless,
    Density,
    Velocity,
    Pressure,
    Temperature,
    Momentum,
    Energy,
    SpecificEnergy,
    Viscosity,
    KinematicViscosity,
    InverseTime,
    InverseTimeSquared,
    Dissipation,
    LengthPower,
};

struct QuantityDescriptor {
    std::string name;
    TopologyLocation location = TopologyLocation::Cell;
    std::vector<std::string> dependencies;
    std::string nondimensional_unit = "1";
    std::string dimensional_unit = "1";
    QuantityScale scale = QuantityScale::Dimensionless;
    int length_power = 0;
    // Statistics normally integrate over the full domain dimension.  A
    // non-negative value overrides that exponent for means and plane integrals.
    int integration_length_power = -1;

    void validate() const;
};

struct QuantityContext {
    GasModel gas;
    ReferenceScales reference;
    NumericalFloors floors;
    TransportModel transport;
    bool dimensional = false;
};

struct QuantityField {
    QuantityDescriptor descriptor;
    Extent3 extent {};
    std::vector<Real> values;
};

class IFieldQuantity {
public:
    virtual ~IFieldQuantity() = default;
    [[nodiscard]] virtual const QuantityDescriptor& descriptor() const = 0;
    [[nodiscard]] virtual Real evaluate_cell(const StructuredBlock& block,
                                             const MetricField& metric,
                                             Index3 index,
                                             const QuantityContext& context) const = 0;
};

class FieldQuantityRegistry {
public:
    FieldQuantityRegistry() = default;
    [[nodiscard]] static FieldQuantityRegistry create_builtin();

    void register_quantity(std::shared_ptr<const IFieldQuantity> quantity);
    void register_turbulence_field(const TurbulenceFieldDescriptor& descriptor);
    [[nodiscard]] bool contains(const std::string& name) const noexcept;
    [[nodiscard]] const QuantityDescriptor& descriptor(const std::string& name) const;
    [[nodiscard]] QuantityField evaluate(const std::string& name,
                                         const StructuredBlock& block,
                                         const MetricField& metric,
                                         const QuantityContext& context) const;
    void validate_selection(const std::vector<std::string>& names) const;

private:
    void validate_dependencies(const std::string& root) const;

    std::unordered_map<std::string, std::shared_ptr<const IFieldQuantity>> quantities_;
};

struct StatisticContext {
    const MpiRuntime& mpi;
    const LocalBlockSet& local_blocks;
    const BlockMetricMap& metrics;
    const StructuredPartitionPlan& partition;
    const GlobalConservationWeights& conservation_weights;
    AlgorithmProfile profile;
    QuantityContext quantities;
    const BlockBoundaryDataMap* boundary_data = nullptr;
    bool viscous = false;
};

class IStatisticQuantity {
public:
    virtual ~IStatisticQuantity() = default;
    [[nodiscard]] virtual const QuantityDescriptor& descriptor() const = 0;
    [[nodiscard]] virtual Real evaluate(const StatisticContext& context) const = 0;
};

class StatisticRegistry {
public:
    StatisticRegistry() = default;
    [[nodiscard]] static StatisticRegistry create_builtin();

    void register_quantity(std::shared_ptr<const IStatisticQuantity> quantity);
    [[nodiscard]] bool contains(const std::string& name) const noexcept;
    [[nodiscard]] Real evaluate(const std::string& name, const StatisticContext& context) const;
    void validate_selection(const std::vector<std::string>& names) const;

private:
    std::unordered_map<std::string, std::shared_ptr<const IStatisticQuantity>> quantities_;
};

[[nodiscard]] std::vector<std::string>
xz_plane_statistic_names(const std::vector<int>& cell_j_indices);
void validate_xz_plane_statistics(const std::vector<int>& cell_j_indices,
                                  const StructuredPartitionPlan& partition);
void register_xz_plane_statistics(StatisticRegistry& registry,
                                  const std::vector<int>& cell_j_indices);

[[nodiscard]] std::vector<std::string> yz_plane_statistic_names(std::size_t plane_count);
void validate_yz_plane_statistics(const std::vector<Real>& target_x_coordinates,
                                  const StructuredPartitionPlan& partition);
void register_yz_plane_statistics(StatisticRegistry& registry,
                                  const std::vector<Real>& target_x_coordinates);

[[nodiscard]] std::vector<std::string> channel_wall_statistic_names();
void register_channel_wall_statistics(StatisticRegistry& registry,
                                      const std::string& lower_patch,
                                      const std::string& upper_patch,
                                      Real half_height);

[[nodiscard]] Real quantity_scale_factor(const QuantityDescriptor& descriptor,
                                         const QuantityContext& context,
                                         int dimension);

} // namespace wcns
