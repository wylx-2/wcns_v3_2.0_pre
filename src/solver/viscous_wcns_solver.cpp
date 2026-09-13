#include <wcns/solver/viscous_wcns_solver.hpp>

#include <wcns/solver/implicit_time_integrator.hpp>
#include <wcns/solver/time_integrator.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

namespace wcns {
namespace {

constexpr std::size_t maximum_exact_diagnostic_count
    = static_cast<std::size_t>(9007199254740992ULL);

bool same_floors(const NumericalFloors& lhs, const NumericalFloors& rhs)
{
    return lhs.density == rhs.density && lhs.pressure == rhs.pressure
        && lhs.temperature == rhs.temperature && lhs.jacobian_absolute == rhs.jacobian_absolute
        && lhs.jacobian_relative == rhs.jacobian_relative
        && lhs.face_area_absolute == rhs.face_area_absolute
        && lhs.face_area_relative == rhs.face_area_relative
        && lhs.reconstruction_scale == rhs.reconstruction_scale
        && lhs.reconstruction_epsilon == rhs.reconstruction_epsilon;
}

bool positive_finite(Real value)
{
    return std::isfinite(value) && value > 0.0;
}

const FaceAreaVectors& faces(const MetricField& metric, Axis axis)
{
    if (axis == Axis::I) return metric.i_faces();
    if (axis == Axis::J) return metric.j_faces();
    return metric.k_faces();
}

std::array<Real, 3> area_vector(const MetricField& metric, Axis axis, Index3 face)
{
    const auto& field = faces(metric, axis);
    return {{
        field.x(face.i, face.j, face.k),
        field.y(face.i, face.j, face.k),
        field.z(face.i, face.j, face.k),
    }};
}

Real area_squared(const std::array<Real, 3>& area)
{
    return area[0] * area[0] + area[1] * area[1] + area[2] * area[2];
}

std::size_t global_diagnostic_count(const MpiRuntime& mpi, std::size_t local, const char* label)
{
    if (local > maximum_exact_diagnostic_count) {
        throw std::overflow_error(std::string(label) + " exceeds exact MPI reduction range");
    }
    const Real global = mpi.sum(static_cast<Real>(local));
    if (!std::isfinite(global) || global < 0.0
        || global > static_cast<Real>(maximum_exact_diagnostic_count)) {
        throw std::overflow_error(std::string(label) + " global reduction is invalid");
    }
    return static_cast<std::size_t>(global);
}

} // namespace

void ViscousStabilityCoefficients::validate() const
{
    if (!positive_finite(phenglei_2d_ssprk3) || !positive_finite(phenglei_3d_ssprk3)
        || !positive_finite(scmm6_2d_ssprk3) || !positive_finite(scmm6_3d_ssprk3)) {
        throw std::invalid_argument("all profile/dimension SSPRK3 viscous stability coefficients "
                                    "must be positive and finite");
    }
}

Real ViscousStabilityCoefficients::for_ssprk3(AlgorithmProfileKind profile, int dimension) const
{
    validate();
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument(
            "viscous SSPRK3 stability coefficient requires dimension 2 or 3");
    }
    if (profile == AlgorithmProfileKind::PhengleiWcns) {
        return dimension == 2 ? phenglei_2d_ssprk3 : phenglei_3d_ssprk3;
    }
    if (profile == AlgorithmProfileKind::Scmm6Wcns) {
        return dimension == 2 ? scmm6_2d_ssprk3 : scmm6_3d_ssprk3;
    }
    throw std::invalid_argument("unknown viscous SSPRK3 algorithm profile");
}

std::string ViscousStabilityCoefficients::summary() const
{
    validate();
    return "Cv_phenglei_2d_ssprk3=" + std::to_string(phenglei_2d_ssprk3)
        + " Cv_phenglei_3d_ssprk3=" + std::to_string(phenglei_3d_ssprk3)
        + " Cv_scmm6_2d_ssprk3=" + std::to_string(scmm6_2d_ssprk3)
        + " Cv_scmm6_3d_ssprk3=" + std::to_string(scmm6_3d_ssprk3);
}

void ViscousWcnsConfig::validate() const
{
    inviscid.validate();
    transport.validate();
    turbulence.validate();
    stability.validate();
}

std::string ViscousWcnsConfig::summary() const
{
    validate();
    auto result = inviscid.summary() + ';' + transport.summary() + "; " + stability.summary();
    if (turbulence.kind != TurbulenceModelKind::None) {
        result += "; " + turbulence.summary();
    }
    return result;
}

std::string ViscousWcnsConfig::restart_signature() const
{
    return std::string(turbulence.kind == TurbulenceModelKind::None ? "viscous_wcns_v1;"
                                                                    : "viscous_wcns_v2;")
        + summary();
}

ViscousWcnsSolver::ViscousWcnsSolver(const MpiRuntime& mpi,
                                     LocalBlockSet& local_blocks,
                                     const StructuredMesh& global_mesh,
                                     const DistributedTopology& topology,
                                     int distribution_rank_count,
                                     BlockMetricMap& metrics,
                                     const BlockBoundaryDataMap& boundary_data,
                                     AlgorithmProfile profile,
                                     GasModel gas,
                                     ReferenceScales reference,
                                     NumericalFloors floors,
                                     ViscousWcnsConfig config)
    : mpi_(mpi)
    , local_blocks_(local_blocks)
    , global_mesh_(global_mesh)
    , topology_(topology)
    , state_exchanger_(mpi, topology, distribution_rank_count, euler_components)
    , turbulence_exchanger_(mpi, topology, distribution_rank_count)
    , metrics_(metrics)
    , boundary_data_(boundary_data)
    , profile_(std::move(profile))
    , gas_(std::move(gas))
    , reference_(std::move(reference))
    , floors_(floors)
    , config_(std::move(config))
    , source_registry_(SourceTermRegistry::create_stage_j(config_.inviscid.source_terms))
    , transport_(config_.transport)
    , inviscid_flux_plan_(FaceFluxHaloPlan::build(global_mesh_, profile_, 1))
    , inviscid_flux_exchanger_(mpi_, inviscid_flux_plan_)
    , operand_plan_(GradientOperandFaceHaloPlan::build(global_mesh_, profile_, 1))
    , operand_exchanger_(mpi_, operand_plan_)
    , gradient_plan_(GradientHaloPlan::build(global_mesh_, topology_, profile_, 1))
    , gradient_exchanger_(mpi_, gradient_plan_)
    , viscous_flux_plan_(ViscousFaceFluxHaloPlan::build(global_mesh_, profile_, 1))
    , viscous_flux_exchanger_(mpi_, viscous_flux_plan_)
    , turbulence_flux_plan_(ViscousFaceFluxHaloPlan::build(global_mesh_, profile_, 1))
    , turbulence_flux_exchanger_(mpi_, turbulence_flux_plan_, 28672)
{
    config_.validate();
    turbulence_model_ = TurbulenceModelRegistry::create_builtin().create(config_.turbulence);
    const auto riemann_registry
        = RiemannSolverRegistry::with_builtins(config_.inviscid.riemann.parameters);
    config_.inviscid.riemann.validate(riemann_registry);
    riemann_ = RiemannSolver(
        config_.inviscid.riemann.scheme, riemann_registry, config_.inviscid.riemann.parameters);
    auto robust_parameters = config_.inviscid.riemann.parameters;
    robust_parameters.weiss_smith = false;
    robust_riemann_ = RiemannSolver("rusanov", riemann_registry, robust_parameters);
    robustness_ladder_
        = RobustnessLadder::build(config_.inviscid.reconstruction, config_.inviscid.riemann);
    floors_.validate();
    if (!same_floors(config_.inviscid.reconstruction.floors, floors_)) {
        throw std::invalid_argument("viscous WCNS reconstruction and solver floors differ");
    }
    ProfileFactory::validate_bundle(profile_.components());
    if (local_blocks_.rank() != mpi_.rank()) {
        throw std::invalid_argument("viscous WCNS local block rank differs from MPI rank");
    }
    for (const auto& block : local_blocks_.blocks()) {
        const auto metric = metrics_.find(block.id());
        if (metric == metrics_.end() || metric->second.profile() != profile_.kind()) {
            throw std::invalid_argument("viscous WCNS metric is missing or incompatible");
        }
        if (boundary_data_.find(block.id()) == boundary_data_.end()) {
            throw std::invalid_argument("viscous WCNS boundary data is missing");
        }
        for (int logical = 0; logical < block.cell_dimension(); ++logical) {
            static_cast<void>(cached_line_operators(
                profile_, block.cell_extent()[static_cast<std::size_t>(logical)]));
        }
    }
    const auto local_count = local_blocks_.blocks().size();
    block_workspace_.reserve(local_count);
    inviscid_flux_workspace_.reserve(local_count);
    operand_workspace_.reserve(local_count);
    gradient_workspace_.reserve(local_count);
    turbulence_gradient_workspace_.reserve(local_count);
    viscous_flux_workspace_.reserve(local_count);
    turbulence_flux_workspace_.reserve(local_count);
    turbulence_residual_workspace_.reserve(local_count);
    turbulence_source_jacobian_workspace_.reserve(local_count);
    for (auto& block : local_blocks_.blocks()) {
        block_workspace_.push_back(&block);
        const auto inviscid = inviscid_flux_workspace_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(block.cell_extent(), block.cell_dimension(), profile_.kind(), 1));
        const auto operand = operand_workspace_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(block.cell_extent(), block.cell_dimension(), profile_.kind(), 1));
        const auto gradient = gradient_workspace_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(block.cell_extent(), block.cell_dimension(), profile_.kind(), 1));
        const auto viscous = viscous_flux_workspace_.emplace(
            std::piecewise_construct,
            std::forward_as_tuple(block.id()),
            std::forward_as_tuple(block.cell_extent(), block.cell_dimension(), profile_.kind(), 1));
        if (!inviscid.second || !operand.second || !gradient.second || !viscous.second) {
            throw std::logic_error("duplicate viscous workspace block");
        }
        inviscid_flux_registry_.add(block.id(), inviscid.first->second);
        operand_registry_.add(block.id(), operand.first->second);
        gradient_registry_.add(block.id(), gradient.first->second);
        viscous_flux_registry_.add(block.id(), viscous.first->second);
        if (turbulence_active()) {
            TurbulenceFieldSet expected(
                block.cell_extent(), block.ghost_width(), turbulence_model_->fields());
            if (block.turbulence.empty()
                || block.turbulence.descriptor_signature() != expected.descriptor_signature()) {
                throw std::invalid_argument(
                    "active SA-neg fields must be initialized before solver construction");
            }
            const auto model_gradient = turbulence_gradient_workspace_.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(block.id()),
                std::forward_as_tuple(
                    block.cell_extent(), block.cell_dimension(), profile_.kind(), 1));
            const auto model_flux = turbulence_flux_workspace_.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(block.id()),
                std::forward_as_tuple(
                    block.cell_extent(), block.cell_dimension(), profile_.kind(), 1));
            const auto model_residual = turbulence_residual_workspace_.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(block.id()),
                std::forward_as_tuple(block.cell_extent(), 1, 0, 0.0));
            const auto source_jacobian = turbulence_source_jacobian_workspace_.emplace(
                std::piecewise_construct,
                std::forward_as_tuple(block.id()),
                std::forward_as_tuple(block.cell_extent(), 1, 0, 0.0));
            if (!model_gradient.second || !model_flux.second || !model_residual.second
                || !source_jacobian.second) {
                throw std::logic_error("duplicate SA-neg workspace block");
            }
            turbulence_gradient_registry_.add(block.id(), model_gradient.first->second);
            turbulence_flux_registry_.add(block.id(), model_flux.first->second);
        }
    }
    const bool any_history = std::any_of(local_blocks_.blocks().begin(),
                                         local_blocks_.blocks().end(),
                                         [](const auto& block) {
                                             return block.flow.implicit_history.valid;
                                         });
    const bool all_history = std::all_of(local_blocks_.blocks().begin(),
                                         local_blocks_.blocks().end(),
                                         [](const auto& block) {
                                             return block.flow.implicit_history.valid;
                                         });
    if (any_history != all_history) {
        throw std::invalid_argument("implicit history is valid on only some viscous blocks");
    }
    if (all_history) {
        previous_physical_time_step_
            = local_blocks_.blocks().front().flow.implicit_history.physical_time_step;
        previous_physical_state_.reserve(local_blocks_.blocks().size());
        for (const auto& block : local_blocks_.blocks()) {
            const auto& history = block.flow.implicit_history;
            if (history.physical_time_step != previous_physical_time_step_
                || history.previous_mean.interior_extent() != block.cell_extent()
                || history.previous_mean.components() != euler_components
                || history.previous_mean.ghost_width() != 0) {
                throw std::invalid_argument("viscous implicit mean history metadata differ");
            }
            previous_physical_state_.push_back(
                {block.id(),
                 block.cell_extent(),
                 std::vector<Real>(history.previous_mean.data(),
                                   history.previous_mean.data() + history.previous_mean.size())});
            if (turbulence_active()) {
                if (history.previous_model.interior_extent() != block.cell_extent()
                    || history.previous_model.components() != 1
                    || history.previous_model.ghost_width() != 0) {
                    throw std::invalid_argument("viscous implicit model history metadata differ");
                }
                previous_turbulence_physical_state_.emplace(
                    block.id(),
                    std::vector<Real>(history.previous_model.data(),
                                      history.previous_model.data()
                                          + history.previous_model.size()));
            }
        }
        has_previous_physical_state_ = true;
    }
}

bool ViscousWcnsSolver::turbulence_active() const noexcept
{
    return config_.turbulence.kind == TurbulenceModelKind::SaNegative;
}

void ViscousWcnsSolver::compute_residuals(Real stage_time, int rk_stage)
{
    compute_residuals_impl(stage_time, rk_stage, nullptr);
}

void ViscousWcnsSolver::compute_residuals_impl(Real stage_time,
                                               int rk_stage,
                                               const BlockFaceRobustnessMap* robustness_levels)
{
    if (!std::isfinite(stage_time)) {
        throw std::invalid_argument("viscous WCNS stage time must be finite");
    }
    if (rk_stage < 0 || rk_stage > 3) {
        throw std::invalid_argument("viscous WCNS RK stage must lie in [0,3]");
    }
    if (version_ == std::numeric_limits<std::uint64_t>::max()) {
        throw std::overflow_error("viscous WCNS state version overflow");
    }
    ++version_;
    BlockFieldRegistry conservative_fields(euler_components);
    for (auto& block : local_blocks_.blocks()) {
        update_temperature_primitive_interior(block, gas_, reference_, floors_);
        conservative_fields.add(block.id(), block.flow.conservative);
    }
    state_exchanger_.exchange(conservative_fields);
    for (const auto& exchange : topology_.exchanges()) {
        if (exchange.receiver_rank != mpi_.rank()) continue;
        auto& receiver = local_blocks_.block(exchange.halo.receiver_block);
        for (const auto& pair : exchange.halo.cell_pairs) {
            update_temperature_primitive_cell(
                receiver, pair.receiver_ghost, gas_, reference_, floors_);
        }
    }

    if (turbulence_active()) {
        synchronize_sa_negative_fields(turbulence_exchanger_,
                                       local_blocks_,
                                       sa_negative_farfield_value(config_.turbulence,
                                                                  reference_.reynolds()));
    }

    reconstruction_diagnostics_ = {};
    riemann_diagnostics_ = {};
    for (auto& block : local_blocks_.blocks()) {
        const auto& data = boundary_data_.at(block.id());
        const auto ghost = PhysicalGhostStateOperator::fill(
            block, data, gas_, reference_, floors_, version_, stage_time);
        if (ghost.version != version_) {
            throw std::logic_error("viscous WCNS physical ghost version mismatch");
        }
        compute_inviscid_face_fluxes_into(
            inviscid_flux_workspace_.at(block.id()),
            block,
            metrics_.at(block.id()),
            profile_,
            config_.inviscid.reconstruction,
            riemann_,
            gas_,
            reference_,
            floors_,
            data,
            config_.inviscid.boundary,
            version_,
            reconstruction_diagnostics_,
            &riemann_diagnostics_,
            rk_stage,
            stage_time,
            robustness_levels == nullptr ? nullptr : &robustness_levels->at(block.id()),
            robustness_levels == nullptr ? nullptr : &robustness_ladder_,
            robustness_levels == nullptr ? nullptr : &robust_riemann_,
            1.0 / reference_.reynolds());
        compute_gradient_face_operands_into(
            operand_workspace_.at(block.id()), block, metrics_.at(block.id()), profile_, version_);
    }
    inviscid_flux_plan_.set_version(version_);
    operand_plan_.set_version(version_);
    inviscid_flux_exchanger_.exchange(inviscid_flux_registry_);
    operand_exchanger_.exchange(operand_registry_);

    for (auto& block : local_blocks_.blocks()) {
        compute_primitive_gradients_into(gradient_workspace_.at(block.id()),
                                         block,
                                         metrics_.at(block.id()),
                                         operand_workspace_.at(block.id()),
                                         profile_);
    }
    gradient_plan_.set_version(version_);
    gradient_exchanger_.exchange(gradient_registry_);

    if (turbulence_active()) {
        for (auto& block : local_blocks_.blocks()) {
            compute_sa_negative_gradient(turbulence_gradient_workspace_.at(block.id()),
                                         block,
                                         metrics_.at(block.id()),
                                         profile_,
                                         version_);
        }
        gradient_exchanger_.exchange(turbulence_gradient_registry_);
    }

    for (auto& block : local_blocks_.blocks()) {
        compute_viscous_face_fluxes_into(viscous_flux_workspace_.at(block.id()),
                                         block,
                                         metrics_.at(block.id()),
                                         gradient_workspace_.at(block.id()),
                                         profile_,
                                         transport_,
                                         boundary_data_.at(block.id()),
                                         gas_,
                                         reference_,
                                         floors_,
                                         version_,
                                         turbulence_active() ? turbulence_model_.get() : nullptr);
        if (turbulence_active()) {
            compute_sa_negative_flux_and_source(
                turbulence_flux_workspace_.at(block.id()),
                turbulence_residual_workspace_.at(block.id()),
                turbulence_source_jacobian_workspace_.at(block.id()),
                block,
                metrics_.at(block.id()),
                inviscid_flux_workspace_.at(block.id()),
                gradient_workspace_.at(block.id()),
                turbulence_gradient_workspace_.at(block.id()),
                profile_,
                *turbulence_model_,
                transport_,
                gas_,
                reference_,
                version_);
        }
    }
    viscous_flux_plan_.set_version(version_);
    viscous_flux_exchanger_.exchange(viscous_flux_registry_);
    if (turbulence_active()) {
        turbulence_flux_plan_.set_version(version_);
        turbulence_flux_exchanger_.exchange(turbulence_flux_registry_);
    }

    for (auto& block : local_blocks_.blocks()) {
        compute_wcns_inviscid_residual(block,
                                       metrics_.at(block.id()),
                                       inviscid_flux_workspace_.at(block.id()),
                                       profile_,
                                       config_.inviscid.flux_difference);
        add_wcns_viscous_residual(block,
                                  metrics_.at(block.id()),
                                  viscous_flux_workspace_.at(block.id()),
                                  profile_,
                                  reference_.reynolds());
        add_source_terms(block, metrics_.at(block.id()), source_registry_, stage_time);
        if (turbulence_active()) {
            assemble_sa_negative_residual(turbulence_residual_workspace_.at(block.id()),
                                          block,
                                          metrics_.at(block.id()),
                                          turbulence_flux_workspace_.at(block.id()),
                                          profile_);
        }
    }
}

void ViscousWcnsSolver::capture_turbulence_stage_state()
{
    turbulence_initial_state_.clear();
    turbulence_stage_state_.clear();
    for (const auto& block : local_blocks_.blocks()) {
        const auto cells = block.cell_extent();
        std::vector<Real> values(cells.size());
        std::size_t offset = 0;
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    values[offset++] = block.flow.conservative(i, j, k, density)
                        * block.turbulence.at({i, j, k}, "nu_tilde");
                }
            }
        }
        turbulence_initial_state_.emplace(block.id(), values);
        turbulence_stage_state_.emplace(block.id(), std::move(values));
    }
}

void ViscousWcnsSolver::update_turbulence_stage(Real initial_weight,
                                                Real stage_weight,
                                                Real residual_weight)
{
    if (!std::isfinite(initial_weight) || !std::isfinite(stage_weight)
        || !positive_finite(residual_weight)) {
        throw std::invalid_argument("SA-neg SSPRK stage weights are invalid");
    }
    for (auto& block : local_blocks_.blocks()) {
        const auto cells = block.cell_extent();
        const auto& initial = turbulence_initial_state_.at(block.id());
        auto& current = turbulence_stage_state_.at(block.id());
        const auto& residual = turbulence_residual_workspace_.at(block.id());
        const auto& jacobian = turbulence_source_jacobian_workspace_.at(block.id());
        std::size_t offset = 0;
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    Real increment = residual_weight * residual(i, j, k, 0);
                    if (config_.turbulence.source_treatment
                        == TurbulenceSourceTreatment::LocalImplicit) {
                        increment = local_implicit_turbulence_increment(
                            residual(i, j, k, 0), jacobian(i, j, k, 0), residual_weight);
                    }
                    const Real candidate = initial_weight * initial[offset]
                        + stage_weight * current[offset] + increment;
                    const Real rho = block.flow.conservative(i, j, k, density);
                    const Real specific = candidate / rho;
                    if (!std::isfinite(candidate) || !positive_finite(rho)
                        || !std::isfinite(specific)) {
                        throw PhysicsError("SA-neg SSPRK update is non-finite");
                    }
                    current[offset] = candidate;
                    block.turbulence.at({i, j, k}, "nu_tilde") = specific;
                    ++offset;
                }
            }
        }
    }
}

namespace {

using ModelStateMap = std::unordered_map<BlockId, std::vector<Real>>;

std::vector<TurbulenceFieldDescriptor>
transported_descriptors(const ITurbulenceModel* model)
{
    std::vector<TurbulenceFieldDescriptor> result;
    if (model == nullptr) return result;
    for (const auto& descriptor : model->fields()) {
        if (descriptor.role == TurbulenceFieldRole::Transported) result.push_back(descriptor);
    }
    return result;
}

const BlockStateBuffer& find_snapshot(const StateSnapshot& snapshot, BlockId block)
{
    const auto iterator = std::find_if(snapshot.begin(), snapshot.end(), [&](const auto& value) {
        return value.block == block;
    });
    if (iterator == snapshot.end()) throw std::logic_error("viscous implicit snapshot is missing");
    return *iterator;
}

ModelStateMap capture_model_state(const LocalBlockSet& local_blocks,
                                  const std::vector<TurbulenceFieldDescriptor>& descriptors)
{
    ModelStateMap result;
    for (const auto& block : local_blocks.blocks()) {
        const auto extent = block.cell_extent();
        std::vector<Real> values(
            extent.size() * static_cast<std::size_t>(descriptors.size()));
        std::size_t offset = 0;
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Real rho = block.flow.conservative(i, j, k, density);
                    for (const auto& descriptor : descriptors) {
                        values[offset++]
                            = rho * block.turbulence.at({i, j, k}, descriptor.name);
                    }
                }
            }
        }
        result.emplace(block.id(), std::move(values));
    }
    return result;
}

struct ViscousImplicitIncrementSet {
    std::unordered_map<BlockId, Field<Real>> mean;
    std::unordered_map<BlockId, Field<Real>> model;
    Real global_residual_l2 = 0.0;
};

ViscousImplicitIncrementSet form_viscous_implicit_increments(
    const MpiRuntime& mpi,
    const LocalBlockSet& local_blocks,
    const BlockMetricMap& metrics,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    const TransportModel& transport,
    const ViscousStabilityCoefficients& stability,
    AlgorithmProfileKind profile,
    const RiemannSolverParameters& riemann_parameters,
    const std::vector<TurbulenceFieldDescriptor>& descriptors,
    const std::unordered_map<BlockId, Field<Real>>& model_residuals,
    const std::unordered_map<BlockId, Field<Real>>& source_jacobians,
    Real pseudo_cfl,
    const LuSgsIterationConfig& iteration,
    Real physical_time_step = 0.0,
    BdfOrder order = BdfOrder::First,
    const StateSnapshot* current_mean = nullptr,
    const StateSnapshot* previous_mean = nullptr,
    const ModelStateMap* current_model = nullptr,
    const ModelStateMap* previous_model = nullptr)
{
    const bool physical = physical_time_step > 0.0;
    if (physical
        != (current_mean != nullptr && previous_mean != nullptr && current_model != nullptr
            && previous_model != nullptr)) {
        throw std::invalid_argument("viscous implicit BDF history inputs are incomplete");
    }
    const Real physical_diagonal
        = physical ? bdf_time_diagonal(order, physical_time_step) : 0.0;
    ViscousImplicitIncrementSet result;
    Real local_square_sum = 0.0;
    Real local_count = 0.0;
    for (const auto& block : local_blocks.blocks()) {
        const auto extent = block.cell_extent();
        const auto& metric = metrics.at(block.id());
        Field<Real> additional(extent, 1, 0, 0.0);
        const Real coefficient = stability.for_ssprk3(profile, block.cell_dimension());
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 cell {i, j, k};
                    Real area_square_sum = 0.0;
                    for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                        const auto axis = static_cast<Axis>(logical);
                        for (int side = 0; side <= 1; ++side) {
                            auto face = cell;
                            face[static_cast<std::size_t>(axis)] += side;
                            area_square_sum += area_squared(area_vector(metric, axis, face));
                        }
                    }
                    const Real rho = block.flow.temperature_primitive(
                        i, j, k, temperature_density);
                    const Real temperature = block.flow.temperature_primitive(
                        i, j, k, temperature_value);
                    const Real molecular = transport.viscosity(temperature)
                        / (rho * reference.reynolds());
                    const Real volume = metric.jacobian()(i, j, k);
                    additional(i, j, k, 0)
                        = coefficient * molecular * area_square_sum
                        / (2.0 * volume * volume);
                }
            }
        }
        const auto* preconditioner = riemann_parameters.weiss_smith
            ? &riemann_parameters.preconditioner
            : nullptr;
        const auto system = build_scalar_spectral_system(block,
                                                         metric,
                                                         gas,
                                                         reference,
                                                         floors,
                                                         pseudo_cfl,
                                                         physical_diagonal,
                                                         preconditioner,
                                                         1.0 / reference.reynolds(),
                                                         &additional);
        Field<Real> mean_rhs(extent, euler_components, 0, 0.0);
        const BlockStateBuffer* mean_n = nullptr;
        const BlockStateBuffer* mean_nm1 = nullptr;
        if (physical) {
            mean_n = &find_snapshot(*current_mean, block.id());
            mean_nm1 = &find_snapshot(*previous_mean, block.id());
        }
        std::size_t mean_offset = 0;
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    for (int component = 0; component < euler_components; ++component) {
                        Real value = block.flow.residual(i, j, k, component);
                        if (physical) {
                            value -= bdf_time_residual(
                                order,
                                block.flow.conservative(i, j, k, component),
                                mean_n->values[mean_offset],
                                mean_nm1->values[mean_offset],
                                physical_time_step);
                        }
                        mean_rhs(i, j, k, component) = value;
                        local_square_sum += value * value;
                        local_count += 1.0;
                        ++mean_offset;
                    }
                }
            }
        }
        result.mean.emplace(block.id(),
                            solve_scalar_lu_sgs(mean_rhs,
                                                system.diagonal,
                                                system.coupling,
                                                block.cell_dimension(),
                                                iteration));
        if (!descriptors.empty()) {
            const int variables = static_cast<int>(descriptors.size());
            const auto& residual = model_residuals.at(block.id());
            const auto& source_jacobian = source_jacobians.at(block.id());
            Field<Real> model_rhs(extent, variables, 0, 0.0);
            Field<Real> model_diagonal(extent, variables, 0, 0.0);
            const auto* model_n = physical ? &current_model->at(block.id()) : nullptr;
            const auto* model_nm1 = physical ? &previous_model->at(block.id()) : nullptr;
            std::size_t model_offset = 0;
            for (int k = 0; k < extent.nk; ++k) {
                for (int j = 0; j < extent.nj; ++j) {
                    for (int i = 0; i < extent.ni; ++i) {
                        const Real rho = block.flow.conservative(i, j, k, density);
                        for (int variable = 0; variable < variables; ++variable) {
                            const Real conservative = rho
                                * block.turbulence.at(
                                    {i, j, k}, descriptors[static_cast<std::size_t>(variable)].name);
                            Real value = residual(i, j, k, variable);
                            if (physical) {
                                value -= bdf_time_residual(order,
                                                           conservative,
                                                           (*model_n)[model_offset],
                                                           (*model_nm1)[model_offset],
                                                           physical_time_step);
                            }
                            model_rhs(i, j, k, variable) = value;
                            model_diagonal(i, j, k, variable)
                                = system.diagonal(i, j, k, 0)
                                + std::max(0.0, -source_jacobian(i, j, k, variable));
                            local_square_sum += value * value;
                            local_count += 1.0;
                            ++model_offset;
                        }
                    }
                }
            }
            result.model.emplace(block.id(),
                                 solve_scalar_lu_sgs(model_rhs,
                                                     model_diagonal,
                                                     system.coupling,
                                                     block.cell_dimension(),
                                                     iteration));
        }
    }
    const Real global_sum = mpi.sum(local_square_sum);
    const Real global_count = mpi.sum(local_count);
    if (!std::isfinite(global_sum) || global_count <= 0.0) {
        throw PhysicsError("viscous implicit residual norm is invalid");
    }
    result.global_residual_l2 = std::sqrt(global_sum / global_count);
    return result;
}

Real commit_viscous_implicit_increment(
    const MpiRuntime& mpi,
    const std::vector<StructuredBlock*>& blocks,
    const std::unordered_map<BlockId, Field<Real>>& mean_increments,
    const std::unordered_map<BlockId, Field<Real>>& model_increments,
    const std::vector<TurbulenceFieldDescriptor>& descriptors,
    const GasModel& gas,
    const ReferenceScales& reference,
    const NumericalFloors& floors,
    Real time)
{
    const auto mean_baseline = capture_conservative_state(blocks);
    ModelStateMap model_baseline;
    for (const auto* block : blocks) {
        const auto extent = block->cell_extent();
        std::vector<Real> values(
            extent.size() * static_cast<std::size_t>(descriptors.size()));
        std::size_t offset = 0;
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Real rho = block->flow.conservative(i, j, k, density);
                    for (const auto& descriptor : descriptors) {
                        values[offset++]
                            = rho * block->turbulence.at({i, j, k}, descriptor.name);
                    }
                }
            }
        }
        model_baseline.emplace(block->id(), std::move(values));
    }
    for (int backtrack = 0; backtrack <= 7; ++backtrack) {
        const Real relaxation = std::ldexp(1.0, -backtrack);
        auto mean_candidate = mean_baseline;
        auto model_candidate = model_baseline;
        bool local_model_valid = true;
        for (std::size_t block_index = 0; block_index < blocks.size(); ++block_index) {
            const auto* block = blocks[block_index];
            const auto extent = block->cell_extent();
            const auto& mean_increment = mean_increments.at(block->id());
            const auto model_iterator = model_increments.find(block->id());
            std::size_t mean_offset = 0;
            std::size_t model_offset = 0;
            for (int k = 0; k < extent.nk; ++k) {
                for (int j = 0; j < extent.nj; ++j) {
                    for (int i = 0; i < extent.ni; ++i) {
                        const std::size_t density_offset = mean_offset;
                        for (int component = 0; component < euler_components; ++component) {
                            mean_candidate[block_index].values[mean_offset++]
                                += relaxation * mean_increment(i, j, k, component);
                        }
                        if (model_iterator != model_increments.end()) {
                            const Real rho
                                = mean_candidate[block_index].values[density_offset];
                            for (std::size_t variable = 0; variable < descriptors.size();
                                 ++variable) {
                                auto& conservative
                                    = model_candidate.at(block->id())[model_offset];
                                conservative += relaxation
                                    * model_iterator->second(
                                        i, j, k, static_cast<int>(variable));
                                const Real specific = conservative / rho;
                                if (!std::isfinite(specific)
                                    || (descriptors[variable].strictly_positive
                                        && specific <= descriptors[variable].lower_bound)
                                    || (!descriptors[variable].strictly_positive
                                        && specific < descriptors[variable].lower_bound)) {
                                    local_model_valid = false;
                                }
                                ++model_offset;
                            }
                        }
                    }
                }
            }
        }
        const auto validation = validate_candidate_state(
            mean_candidate, blocks, gas, reference, floors, 1, time);
        if (mpi.all_true(validation.valid() && local_model_valid)) {
            commit_candidate_state(blocks, mean_candidate);
            for (auto* block : blocks) {
                const auto extent = block->cell_extent();
                const auto& values = model_candidate.at(block->id());
                std::size_t offset = 0;
                for (int k = 0; k < extent.nk; ++k) {
                    for (int j = 0; j < extent.nj; ++j) {
                        for (int i = 0; i < extent.ni; ++i) {
                            const Real rho = block->flow.conservative(i, j, k, density);
                            for (const auto& descriptor : descriptors) {
                                block->turbulence.at({i, j, k}, descriptor.name)
                                    = values[offset++] / rho;
                            }
                        }
                    }
                }
            }
            return relaxation;
        }
    }
    throw PhysicsError("viscous LU-SGS increment remained inadmissible after backtracking");
}

} // namespace

Real ViscousWcnsSolver::advance(Real time_step, Real initial_time)
{
    Real accepted_time_step = time_step;
    if (!config_.inviscid.robustness.enabled) {
        robustness_diagnostics_ = {};
        robustness_diagnostics_.proposed_time_step = time_step;
        robustness_diagnostics_.accepted_time_step = time_step;
        int rk_stage = 0;
        if (turbulence_active()) capture_turbulence_stage_state();
        advance_ssprk3(
            block_workspace_,
            time_workspace_,
            time_step,
            initial_time,
            [this, &rk_stage, time_step](Real stage_time) {
                ++rk_stage;
                if (turbulence_active() && rk_stage == 2) {
                    update_turbulence_stage(1.0, 0.0, time_step);
                } else if (turbulence_active() && rk_stage == 3) {
                    update_turbulence_stage(0.75, 0.25, 0.25 * time_step);
                }
                compute_residuals(stage_time, rk_stage);
            });
        if (turbulence_active()) {
            update_turbulence_stage(1.0 / 3.0, 2.0 / 3.0, 2.0 * time_step / 3.0);
        }
    } else {
        accepted_time_step = advance_ssprk3_with_robustness(
            mpi_,
            block_workspace_,
            global_mesh_,
            profile_,
            config_.inviscid.flux_difference,
            gas_,
            reference_,
            floors_,
            config_.inviscid.robustness,
            robustness_ladder_,
            time_step,
            initial_time,
            [this](Real stage_time, int rk_stage, const BlockFaceRobustnessMap& levels) {
                compute_residuals_impl(stage_time, rk_stage, &levels);
            },
            robustness_diagnostics_);
    }
    for (auto& block : local_blocks_.blocks()) {
        update_temperature_primitive_interior(block, gas_, reference_, floors_);
    }
    return accepted_time_step;
}

Real ViscousWcnsSolver::advance_lu_sgs(Real pseudo_cfl,
                                       Real time,
                                       const LuSgsIterationConfig& iteration)
{
    iteration.validate();
    compute_residuals(time);
    const auto descriptors
        = transported_descriptors(turbulence_active() ? turbulence_model_.get() : nullptr);
    const auto increments = form_viscous_implicit_increments(mpi_,
                                                             local_blocks_,
                                                             metrics_,
                                                             gas_,
                                                             reference_,
                                                             floors_,
                                                             transport_,
                                                             config_.stability,
                                                             profile_.kind(),
                                                             config_.inviscid.riemann.parameters,
                                                             descriptors,
                                                             turbulence_residual_workspace_,
                                                             turbulence_source_jacobian_workspace_,
                                                             pseudo_cfl,
                                                             iteration);
    static_cast<void>(commit_viscous_implicit_increment(mpi_,
                                                        block_workspace_,
                                                        increments.mean,
                                                        increments.model,
                                                        descriptors,
                                                        gas_,
                                                        reference_,
                                                        floors_,
                                                        time));
    for (auto& block : local_blocks_.blocks())
        update_temperature_primitive_interior(block, gas_, reference_, floors_);
    robustness_diagnostics_ = {};
    robustness_diagnostics_.proposed_time_step = 1.0;
    robustness_diagnostics_.accepted_time_step = 1.0;
    return 1.0;
}

Real ViscousWcnsSolver::advance_dual_time(Real physical_time_step,
                                          Real initial_time,
                                          const DualTimeIterationConfig& config)
{
    config.validate();
    if (!positive_finite(physical_time_step)) {
        throw std::invalid_argument("viscous dual-time physical step is invalid");
    }
    const auto descriptors
        = transported_descriptors(turbulence_active() ? turbulence_model_.get() : nullptr);
    const auto current_mean = capture_conservative_state(block_workspace_);
    const auto current_model = capture_model_state(local_blocks_, descriptors);
    const auto previous_mean
        = has_previous_physical_state_ ? previous_physical_state_ : current_mean;
    const auto previous_model = has_previous_physical_state_
        ? previous_turbulence_physical_state_
        : current_model;
    const Real step_scale
        = std::max({Real {1.0}, physical_time_step, previous_physical_time_step_});
    const bool uniform_history = has_previous_physical_state_
        && std::abs(previous_physical_time_step_ - physical_time_step)
            <= 64.0 * std::numeric_limits<Real>::epsilon() * step_scale;
    const BdfOrder order = uniform_history ? BdfOrder::Second : BdfOrder::First;
    Real initial_norm = 0.0;
    for (std::size_t iteration = 0; iteration <= config.max_iterations; ++iteration) {
        compute_residuals(initial_time + physical_time_step);
        const auto increments = form_viscous_implicit_increments(
            mpi_,
            local_blocks_,
            metrics_,
            gas_,
            reference_,
            floors_,
            transport_,
            config_.stability,
            profile_.kind(),
            config_.inviscid.riemann.parameters,
            descriptors,
            turbulence_residual_workspace_,
            turbulence_source_jacobian_workspace_,
            config.cfl,
            config.lu_sgs,
            physical_time_step,
            order,
            &current_mean,
            &previous_mean,
            &current_model,
            &previous_model);
        if (iteration == 0) initial_norm = increments.global_residual_l2;
        const bool converged = increments.global_residual_l2 <= config.absolute_tolerance
            || (initial_norm > 0.0
                && increments.global_residual_l2 / initial_norm <= config.relative_tolerance);
        if (converged) {
            previous_physical_state_ = current_mean;
            previous_turbulence_physical_state_ = current_model;
            has_previous_physical_state_ = true;
            previous_physical_time_step_ = physical_time_step;
            for (auto& block : local_blocks_.blocks()) {
                auto& history = block.flow.implicit_history;
                history.prepare(block.cell_extent(), static_cast<int>(descriptors.size()));
                history.valid = true;
                history.physical_time_step = physical_time_step;
                const auto& mean = find_snapshot(current_mean, block.id()).values;
                std::copy(mean.begin(), mean.end(), history.previous_mean.data());
                if (!descriptors.empty()) {
                    const auto& model = current_model.at(block.id());
                    std::copy(model.begin(), model.end(), history.previous_model.data());
                }
            }
            for (auto& block : local_blocks_.blocks())
                update_temperature_primitive_interior(block, gas_, reference_, floors_);
            return physical_time_step;
        }
        if (iteration == config.max_iterations) break;
        static_cast<void>(commit_viscous_implicit_increment(mpi_,
                                                            block_workspace_,
                                                            increments.mean,
                                                            increments.model,
                                                            descriptors,
                                                            gas_,
                                                            reference_,
                                                            floors_,
                                                            initial_time
                                                                + physical_time_step));
    }
    restore_conservative_state(block_workspace_, current_mean);
    for (auto& block : local_blocks_.blocks()) {
        const auto extent = block.cell_extent();
        const auto& values = current_model.at(block.id());
        std::size_t offset = 0;
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Real rho = block.flow.conservative(i, j, k, density);
                    for (const auto& descriptor : descriptors) {
                        block.turbulence.at({i, j, k}, descriptor.name)
                            = values[offset++] / rho;
                    }
                }
            }
        }
        update_temperature_primitive_interior(block, gas_, reference_, floors_);
    }
    throw PhysicsError("viscous dual-time LU-SGS failed; physical layer was restored");
}

std::size_t ViscousWcnsSolver::global_reconstruction_fallback_count() const
{
    return global_diagnostic_count(mpi_,
                                   reconstruction_diagnostics_.fallback_events.size(),
                                   "viscous reconstruction fallback count");
}

std::size_t ViscousWcnsSolver::global_riemann_fallback_count() const
{
    return global_diagnostic_count(
        mpi_, riemann_diagnostics_.fallback_count(), "viscous Riemann fallback count");
}

std::size_t ViscousWcnsSolver::global_riemann_face_count() const
{
    return global_diagnostic_count(
        mpi_, riemann_diagnostics_.total_faces, "viscous Riemann face count");
}

RobustnessDiagnostics ViscousWcnsSolver::global_robustness_diagnostics() const
{
    RobustnessDiagnostics result = robustness_diagnostics_;
    for (std::size_t level = 0; level < result.face_levels.size(); ++level) {
        result.face_levels[level]
            = global_diagnostic_count(mpi_,
                                      robustness_diagnostics_.face_levels[level],
                                      "viscous robustness face-level count");
    }
    result.troubled_cells = global_diagnostic_count(
        mpi_, robustness_diagnostics_.troubled_cells, "viscous robustness troubled-cell count");
    result.local_recomputations = static_cast<std::size_t>(
        mpi_.max(static_cast<Real>(robustness_diagnostics_.local_recomputations)));
    result.step_retries = static_cast<std::size_t>(
        mpi_.max(static_cast<Real>(robustness_diagnostics_.step_retries)));
    result.minimum_density = mpi_.min(robustness_diagnostics_.minimum_density);
    result.minimum_pressure = mpi_.min(robustness_diagnostics_.minimum_pressure);
    result.minimum_temperature = mpi_.min(robustness_diagnostics_.minimum_temperature);
    result.minimum_internal_energy = mpi_.min(robustness_diagnostics_.minimum_internal_energy);
    return result;
}

Real ViscousWcnsSolver::global_time_step(Real cfl)
{
    if (!positive_finite(cfl)) {
        throw std::invalid_argument("viscous WCNS CFL must be positive and finite");
    }
    Real local_minimum = std::numeric_limits<Real>::infinity();
    for (auto& block : local_blocks_.blocks()) {
        update_temperature_primitive_interior(block, gas_, reference_, floors_);
        const auto& metric = metrics_.at(block.id());
        const auto cells = block.cell_extent();
        const Real viscous_stability
            = config_.stability.for_ssprk3(profile_.kind(), block.cell_dimension());
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    const Index3 cell {i, j, k};
                    TemperaturePrimitiveState state {};
                    for (int component = 0; component < fluid_components; ++component) {
                        state[static_cast<std::size_t>(component)]
                            = block.flow.temperature_primitive(i, j, k, component);
                    }
                    const std::array<Real, 3> velocity {{
                        state[temperature_velocity_x],
                        state[temperature_velocity_y],
                        state[temperature_velocity_z],
                    }};
                    const Real sound = thermodynamic_sound_speed(
                        state, gas_, reference_, floors_, block.cell_dimension());
                    Real inviscid_sum = 0.0;
                    Real area_square_sum = 0.0;
                    for (int logical = 0; logical < block.cell_dimension(); ++logical) {
                        const auto axis = static_cast<Axis>(logical);
                        for (int side = 0; side <= 1; ++side) {
                            auto face = cell;
                            face[static_cast<std::size_t>(axis)] += side;
                            const auto area = area_vector(metric, axis, face);
                            const Real area2 = area_squared(area);
                            inviscid_sum += std::abs(velocity[0] * area[0] + velocity[1] * area[1]
                                                     + velocity[2] * area[2])
                                + sound * std::sqrt(area2);
                            area_square_sum += area2;
                        }
                    }
                    const Real jacobian = metric.jacobian()(i, j, k);
                    const Real mu = transport_.viscosity(state[temperature_value]);
                    const Real rho = state[temperature_density];
                    const Real nu_effective
                        = std::max(mu / rho, mu / (rho * config_.transport.prandtl));
                    Real mean_diffusivity = nu_effective / reference_.reynolds();
                    Real model_diffusivity = 0.0;
                    Real source_rate = 0.0;
                    if (turbulence_active()) {
                        TurbulenceCellContext context;
                        context.mean_state = state;
                        context.model_values
                            = {block.turbulence.at(cell, "nu_tilde")};
                        context.model_gradients = {{{0.0, 0.0, 0.0}}};
                        context.molecular_kinematic_viscosity
                            = mu / (rho * reference_.reynolds());
                        context.wall_distance
                            = block.turbulence.at(cell, "wall_distance");
                        context.reference_reynolds = reference_.reynolds();
                        context.reference_mach = reference_.mach();
                        context.heat_capacity_ratio = gas_.gamma();
                        context.dimension = block.cell_dimension();
                        const auto evaluation = evaluate_sa_negative(context);
                        mean_diffusivity = std::max(
                            mean_diffusivity,
                            mu / (rho * reference_.reynolds())
                                + evaluation.eddy_kinematic_viscosity);
                        model_diffusivity = evaluation.diffusion_coefficient;
                        if (config_.turbulence.source_treatment
                            == TurbulenceSourceTreatment::Explicit) {
                            source_rate = std::max(0.0, -evaluation.source_derivative);
                        }
                    }
                    const Real denominator = inviscid_sum / (2.0 * jacobian)
                        + viscous_stability * std::max(mean_diffusivity, model_diffusivity)
                            * area_square_sum / (2.0 * jacobian * jacobian)
                        + source_rate;
                    if (!positive_finite(denominator)) {
                        throw PhysicsError("viscous WCNS time-step denominator is invalid");
                    }
                    local_minimum = std::min(local_minimum, cfl / denominator);
                }
            }
        }
    }
    const Real result = mpi_.min(local_minimum);
    if (!positive_finite(result)) {
        throw PhysicsError("viscous WCNS global time step is invalid");
    }
    return result;
}

Real ViscousWcnsSolver::global_residual_l2() const
{
    Real local_sum = 0.0;
    Real local_count = 0.0;
    for (const auto& block : local_blocks_.blocks()) {
        const auto extent = block.cell_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    for (int component = 0; component < euler_components; ++component) {
                        const Real value = block.flow.residual(i, j, k, component);
                        local_sum += value * value;
                        local_count += 1.0;
                    }
                }
            }
        }
    }
    const Real count = mpi_.sum(local_count);
    if (count <= 0.0) throw PhysicsError("viscous residual norm has no cells");
    return std::sqrt(mpi_.sum(local_sum) / count);
}

} // namespace wcns
