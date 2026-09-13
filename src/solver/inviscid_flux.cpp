#include <wcns/solver/inviscid_flux.hpp>

#include <wcns/solver/robustness.hpp>

#include <wcns/mesh/linear_operators.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>
#include <tuple>

namespace wcns {
namespace {

constexpr std::uint64_t maximum_exact_message_version = 9007199254740992ULL;

int side_sign(Side side)
{
    return side == Side::Lower ? -1 : 1;
}

int inward_sign(Side side)
{
    return side == Side::Lower ? 1 : -1;
}

bool contains(const IndexRange3& range, Index3 index)
{
    for (int axis = 0; axis < 3; ++axis) {
        const auto a = static_cast<std::size_t>(axis);
        const int lower = std::min(range.begin[a], range.end[a]);
        const int upper = std::max(range.begin[a], range.end[a]);
        if (index[a] < lower || index[a] > upper) {
            return false;
        }
    }
    return true;
}

template <class LeftRange, class RightRange>
bool same_undirected_range(const LeftRange& lhs, const RightRange& rhs)
{
    return (lhs.begin == rhs.begin && lhs.end == rhs.end)
        || (lhs.begin == rhs.end && lhs.end == rhs.begin);
}

bool reciprocal_pair(const ConnectivityPatch& connection,
                     const ConnectivityPatch& candidate,
                     int dimension)
{
    return candidate.receiver_block == connection.donor_block
        && candidate.donor_block == connection.receiver_block
        && candidate.receiver_face == connection.donor_face
        && candidate.donor_face == connection.receiver_face
        && same_undirected_range(candidate.receiver_vertex_range,
                                 connection.donor_vertex_range)
        && same_undirected_range(candidate.donor_vertex_range,
                                 connection.receiver_vertex_range)
        && candidate.transform == connection.transform.inverse(dimension)
        && candidate.periodic == connection.periodic.inverse();
}

const ConnectivityPatch& reciprocal(const StructuredMesh& mesh, const ConnectivityPatch& connection)
{
    const auto& donor = mesh.block(connection.donor_block);
    const auto iterator = std::find_if(
        donor.connectivities.begin(),
        donor.connectivities.end(),
        [&](const ConnectivityPatch& candidate) {
            return reciprocal_pair(connection, candidate, donor.cell_dimension());
        });
    if (iterator == donor.connectivities.end()) {
        throw TopologyError("face-flux plan cannot find reciprocal connectivity");
    }
    return *iterator;
}

auto connection_key(const ConnectivityPatch& connection)
{
    return std::tuple {
        connection.receiver_block,
        connection.donor_block,
        static_cast<int>(connection.receiver_face.axis),
        static_cast<int>(connection.receiver_face.side),
        connection.receiver_vertex_range.begin.i,
        connection.receiver_vertex_range.begin.j,
        connection.receiver_vertex_range.begin.k,
        connection.receiver_vertex_range.end.i,
        connection.receiver_vertex_range.end.j,
        connection.receiver_vertex_range.end.k,
        connection.name,
    };
}

FaceFluxExchangeDescriptor make_descriptor(const StructuredMesh& mesh,
                                           const ConnectivityPatch& connection,
                                           ConnectionId id,
                                           BlockId owner,
                                           bool receiver_is_owner,
                                           const AlgorithmProfile& profile,
                                           std::uint64_t version)
{
    FaceFluxExchangeDescriptor descriptor;
    descriptor.connection = id;
    descriptor.receiver_block = connection.receiver_block;
    descriptor.donor_block = connection.donor_block;
    descriptor.receiver_rank = mesh.block(connection.receiver_block).owner_rank();
    descriptor.donor_rank = mesh.block(connection.donor_block).owner_rank();
    descriptor.shared_face_owner = owner;
    descriptor.direction = receiver_is_owner ? 0 : 1;
    descriptor.receiver_axis = connection.receiver_face.axis;
    descriptor.donor_axis = connection.donor_face.axis;
    descriptor.orientation = static_cast<Real>(-side_sign(connection.receiver_face.side)
                                               / side_sign(connection.donor_face.side));
    descriptor.periodic = connection.periodic;
    descriptor.profile = profile.kind();
    descriptor.version = version;
    const int maximum_layer = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 1 : 2;
    const int first_layer = receiver_is_owner ? 1 : 0;
    const auto& reverse = reciprocal(mesh, connection);
    const auto counts = connection.shared_face_range.counts();
    for (int layer = first_layer; layer <= maximum_layer; ++layer) {
        for (int k = 0; k < counts.nk; ++k) {
            for (int j = 0; j < counts.nj; ++j) {
                for (int i = 0; i < counts.ni; ++i) {
                    const Index3 receiver_ordinal {i, j, k};
                    Index3 donor_ordinal;
                    for (int receiver_axis = 0;
                         receiver_axis < mesh.block(connection.receiver_block).cell_dimension();
                         ++receiver_axis) {
                        const int donor_axis
                            = std::abs(
                                  connection.transform
                                      .receiver_to_donor[static_cast<std::size_t>(receiver_axis)])
                            - 1;
                        donor_ordinal[static_cast<std::size_t>(donor_axis)]
                            = receiver_ordinal[static_cast<std::size_t>(receiver_axis)];
                    }
                    auto receiver_index = connection.shared_face_range.at(receiver_ordinal);
                    auto donor_index = reverse.shared_face_range.at(donor_ordinal);
                    receiver_index[static_cast<std::size_t>(connection.receiver_face.axis)]
                        += side_sign(connection.receiver_face.side) * layer;
                    donor_index[static_cast<std::size_t>(connection.donor_face.axis)]
                        += inward_sign(connection.donor_face.side) * layer;
                    descriptor.pairs.push_back({receiver_index, donor_index, layer});
                }
            }
        }
    }
    return descriptor;
}

ConservativeState transform_flux_impl(const ConservativeState& donor,
                                      const FaceFluxExchangeDescriptor& descriptor)
{
    ConservativeState result = donor;
    const std::array<Real, 3> momentum {{donor[1], donor[2], donor[3]}};
    const auto rotated = descriptor.periodic.inverse().apply_vector(momentum);
    result[0] *= descriptor.orientation;
    result[1] = descriptor.orientation * rotated[0];
    result[2] = descriptor.orientation * rotated[1];
    result[3] = descriptor.orientation * rotated[2];
    result[4] *= descriptor.orientation;
    return result;
}

ConservativeState load_flux(const InviscidFaceFluxField& field, Axis axis, Index3 index)
{
    const auto& values = field.field(axis);
    ConservativeState result {};
    for (int component = 0; component < euler_components; ++component) {
        result[static_cast<std::size_t>(component)] = values(index.i, index.j, index.k, component);
    }
    return result;
}

void store_flux(InviscidFaceFluxField& field,
                Axis axis,
                Index3 index,
                const ConservativeState& state)
{
    auto& values = field.field(axis);
    for (int component = 0; component < euler_components; ++component) {
        values(index.i, index.j, index.k, component) = state[static_cast<std::size_t>(component)];
    }
}

void validate_field(const InviscidFaceFluxField& field,
                    const FaceFluxExchangeDescriptor& descriptor)
{
    if (field.profile() != descriptor.profile || field.version() != descriptor.version) {
        throw std::invalid_argument("face-flux field profile or version mismatch");
    }
}

const BoundaryPatch* physical_patch(const StructuredBlock& block, Axis axis, Index3 face)
{
    for (const auto& patch : block.boundaries) {
        if (patch.face.axis == axis && contains(patch.boundary_face_range.untyped(), face)) {
            return &patch;
        }
    }
    return nullptr;
}

const FaceAreaVectors& metric_faces(const MetricField& metric, Axis axis)
{
    switch (axis) {
    case Axis::I: return metric.i_faces();
    case Axis::J: return metric.j_faces();
    case Axis::K: return metric.k_faces();
    }
    throw std::invalid_argument("invalid metric face axis");
}

Normal3 unit_normal(const FaceAreaVectors& metric, Index3 face, Real& area)
{
    area = metric.area(face.i, face.j, face.k);
    if (!std::isfinite(area) || area <= 0.0) {
        throw PhysicsError("inviscid face has invalid area");
    }
    return {
        metric.x(face.i, face.j, face.k) / area,
        metric.y(face.i, face.j, face.k) / area,
        metric.z(face.i, face.j, face.k) / area,
    };
}

Normal3 outward(Normal3 positive, Side side)
{
    const Real sign = side == Side::Lower ? -1.0 : 1.0;
    return {sign * positive.x, sign * positive.y, sign * positive.z};
}

bool connection_covers(const StructuredBlock& block, Axis axis, Side side, Index3 face)
{
    for (const auto& connection : block.connectivities) {
        if (connection.receiver_face.axis == axis && connection.receiver_face.side == side
            && contains(connection.shared_face_range.untyped(), face)) {
            return true;
        }
    }
    return false;
}

bool non_owned_connection_face(const StructuredBlock& block, Axis axis, Index3 face)
{
    for (const auto& connection : block.connectivities) {
        if (connection.receiver_face.axis == axis
            && contains(connection.shared_face_range.untyped(), face)) {
            const BlockId owner = std::min(connection.receiver_block, connection.donor_block);
            return block.id() != owner;
        }
    }
    return false;
}

#if WCNS_HAS_MPI
int mpi_count(std::size_t count)
{
    if (count > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("face-flux MPI message exceeds int count");
    }
    return static_cast<int>(count);
}
#endif

} // namespace

const char* flux_difference_mode_name(FluxDifferenceMode mode)
{
    switch (mode) {
    case FluxDifferenceMode::Profile: return "profile";
    case FluxDifferenceMode::ConservativeTwoPoint: return "conservative_two_point";
    }
    throw std::invalid_argument("unknown flux-difference mode");
}

bool is_non_owned_connection_face(const StructuredBlock& block, Axis axis, Index3 face)
{
    return non_owned_connection_face(block, axis, face);
}

StencilRow inviscid_residual_stencil(const StructuredBlock& block,
                                     const AlgorithmProfile& profile,
                                     FluxDifferenceMode mode,
                                     Axis axis,
                                     Index3 cell)
{
    ProfileFactory::validate_bundle(profile.components());
    const auto cells = block.cell_extent();
    if (static_cast<int>(axis) >= block.cell_dimension() || cell.i < 0 || cell.i >= cells.ni
        || cell.j < 0 || cell.j >= cells.nj || cell.k < 0 || cell.k >= cells.nk) {
        throw std::out_of_range("inviscid residual stencil cell/axis is invalid");
    }
    const int normal = cell[static_cast<std::size_t>(axis)];
    if (mode == FluxDifferenceMode::ConservativeTwoPoint) {
        return {{normal, -1.0}, {normal + 1, 1.0}};
    }
    static_cast<void>(flux_difference_mode_name(mode));
    const int count = cells[static_cast<std::size_t>(axis)];
    Index3 lower_face = cell;
    lower_face[static_cast<std::size_t>(axis)] = 0;
    Index3 upper_face = cell;
    upper_face[static_cast<std::size_t>(axis)] = count;
    const bool lower_connection = connection_covers(block, axis, Side::Lower, lower_face);
    const bool upper_connection = connection_covers(block, axis, Side::Upper, upper_face);
    const int boundary_width = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 1 : 2;
    if ((lower_connection && normal < boundary_width)
        || (upper_connection && normal >= count - boundary_width)) {
        StencilRow result;
        if (profile.kind() == AlgorithmProfileKind::PhengleiWcns) {
            constexpr std::array<int, 4> offsets {{-1, 0, 1, 2}};
            constexpr std::array<Real, 4> coefficients {
                {1.0 / 24.0, -27.0 / 24.0, 27.0 / 24.0, -1.0 / 24.0}};
            for (std::size_t index = 0; index < offsets.size(); ++index) {
                result.emplace_back(normal + offsets[index], coefficients[index]);
            }
        } else {
            constexpr std::array<int, 6> offsets {{-2, -1, 0, 1, 2, 3}};
            constexpr std::array<Real, 6> coefficients {{-9.0 / 1920.0,
                                                         125.0 / 1920.0,
                                                         -2250.0 / 1920.0,
                                                         2250.0 / 1920.0,
                                                         -125.0 / 1920.0,
                                                         9.0 / 1920.0}};
            for (std::size_t index = 0; index < offsets.size(); ++index) {
                result.emplace_back(normal + offsets[index], coefficients[index]);
            }
        }
        return result;
    }
    return cached_line_operators(profile, count)
        .derivative_rows()[static_cast<std::size_t>(normal)];
}

ConservativeState
transform_inviscid_face_flux_for_receiver(const ConservativeState& donor,
                                          const FaceFluxExchangeDescriptor& descriptor)
{
    return transform_flux_impl(donor, descriptor);
}

InviscidFaceFluxField::InviscidFaceFluxField(Extent3 cells,
                                             int dimension,
                                             AlgorithmProfileKind profile,
                                             std::uint64_t version)
    : profile_(profile)
    , version_(version)
    , halo_layers_(profile == AlgorithmProfileKind::PhengleiWcns ? 1 : 2)
    , dimension_(dimension)
    , i_({cells.ni + 1, cells.nj, cells.nk}, euler_components, halo_layers_)
    , j_({cells.ni, cells.nj + 1, cells.nk}, euler_components, halo_layers_)
    , k_({cells.ni, cells.nj, cells.nk + 1}, euler_components, halo_layers_)
{
    if ((dimension != 2 && dimension != 3) || version == 0
        || version > maximum_exact_message_version) {
        throw std::invalid_argument("face-flux field has invalid dimension or version");
    }
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    i_.fill(nan);
    j_.fill(nan);
    k_.fill(nan);
}

Field<Real>& InviscidFaceFluxField::field(Axis axis)
{
    switch (axis) {
    case Axis::I: return i_;
    case Axis::J: return j_;
    case Axis::K:
        if (dimension_ != 3) throw std::out_of_range("2D face flux has no K field");
        return k_;
    }
    throw std::invalid_argument("invalid face-flux axis");
}

const Field<Real>& InviscidFaceFluxField::field(Axis axis) const
{
    return const_cast<InviscidFaceFluxField*>(this)->field(axis);
}

void InviscidFaceFluxField::reset(std::uint64_t version)
{
    if (version == 0 || version > maximum_exact_message_version) {
        throw std::invalid_argument("face-flux reset version is invalid");
    }
    version_ = version;
    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    i_.fill(nan);
    j_.fill(nan);
    k_.fill(nan);
}

int FaceFluxExchangeDescriptor::message_tag(int tag_base) const
{
    if (tag_base < 0 || connection < 0 || receiver_block < 0 || donor_block < 0) {
        throw TopologyError("face-flux message tag inputs are invalid");
    }
    const long long tag = static_cast<long long>(tag_base) + 8LL * connection
        + 2LL * static_cast<int>(profile) + direction;
    if (tag > std::numeric_limits<int>::max()) {
        throw TopologyError("face-flux message tag exceeds int range");
    }
    return static_cast<int>(tag);
}

FaceFluxHaloPlan FaceFluxHaloPlan::build(const StructuredMesh& mesh,
                                         const AlgorithmProfile& profile,
                                         std::uint64_t version)
{
    if (version == 0 || version > maximum_exact_message_version) {
        throw TopologyError("face-flux plan version must be non-zero");
    }
    // Distributed production meshes keep only topology for remote blocks;
    // CGNS coordinate continuity is validated during mesh ingestion.
    mesh.validate_connectivities(false);
    std::vector<const ConnectivityPatch*> canonical;
    for (const auto& block : mesh.blocks()) {
        for (const auto& connection : block.connectivities) {
            if (connection.receiver_block < connection.donor_block) {
                canonical.push_back(&connection);
            } else if (connection.receiver_block == connection.donor_block) {
                const auto reciprocal_iterator = std::find_if(
                    block.connectivities.begin(),
                    block.connectivities.end(),
                    [&](const ConnectivityPatch& candidate) {
                        return &candidate != &connection
                            && reciprocal_pair(
                                connection, candidate, block.cell_dimension());
                    });
                if (reciprocal_iterator != block.connectivities.end()
                    && connection_key(connection) < connection_key(*reciprocal_iterator)) {
                    canonical.push_back(&connection);
                } else if (reciprocal_pair(connection, connection, block.cell_dimension())) {
                    canonical.push_back(&connection);
                }
            }
        }
    }
    std::sort(canonical.begin(), canonical.end(), [](const auto* lhs, const auto* rhs) {
        return connection_key(*lhs) < connection_key(*rhs);
    });
    FaceFluxHaloPlan result;
    for (std::size_t index = 0; index < canonical.size(); ++index) {
        const auto id = static_cast<ConnectionId>(index);
        const auto& forward = *canonical[index];
        const auto& reverse = reciprocal(mesh, forward);
        const BlockId owner = std::min(forward.receiver_block, forward.donor_block);
        result.exchanges_.push_back(
            make_descriptor(mesh, forward, id, owner, true, profile, version));
        if (&forward != &reverse) {
            result.exchanges_.push_back(
                make_descriptor(mesh, reverse, id, owner, false, profile, version));
        }
    }
    std::set<int> tags;
    for (const auto& descriptor : result.exchanges_) {
        if (!tags.insert(descriptor.message_tag()).second) {
            throw TopologyError("face-flux plan generated duplicate message tags");
        }
    }
    return result;
}

void FaceFluxHaloPlan::set_version(std::uint64_t version)
{
    if (version == 0 || version > maximum_exact_message_version) {
        throw TopologyError("face-flux plan version must be non-zero");
    }
    for (auto& descriptor : exchanges_)
        descriptor.version = version;
}

void FaceFluxFieldRegistry::add(BlockId block, InviscidFaceFluxField& field)
{
    if (block < 0 || !fields_.emplace(block, &field).second) {
        throw std::invalid_argument("face-flux registry contains invalid or duplicate block");
    }
}

bool FaceFluxFieldRegistry::contains(BlockId block) const noexcept
{
    return fields_.find(block) != fields_.end();
}

InviscidFaceFluxField& FaceFluxFieldRegistry::field(BlockId block) const
{
    const auto iterator = fields_.find(block);
    if (iterator == fields_.end()) {
        throw std::out_of_range("face-flux field is not registered");
    }
    return *iterator->second;
}

void FaceFluxHaloExchanger::prepare()
{
    const RankId rank = mpi_.rank();
    for (const auto& descriptor : plan_.exchanges()) {
        const std::size_t count
            = 1 + descriptor.pairs.size() * static_cast<std::size_t>(euler_components);
        if (descriptor.receiver_rank == rank && descriptor.donor_rank != rank) {
            receives_.push_back({&descriptor, std::vector<Real>(count)});
        } else if (descriptor.donor_rank == rank && descriptor.receiver_rank != rank) {
            sends_.push_back({&descriptor, std::vector<Real>(count)});
        }
    }
#if WCNS_HAS_MPI
    requests_.resize(receives_.size() + sends_.size(), MPI_REQUEST_NULL);
#endif
}

void FaceFluxHaloExchanger::exchange(const FaceFluxFieldRegistry& fields) const
{
    const RankId rank = mpi_.rank();
    for (const auto& descriptor : plan_.exchanges()) {
        if (descriptor.receiver_rank == rank && descriptor.donor_rank == rank) {
            auto& receiver = fields.field(descriptor.receiver_block);
            const auto& donor = fields.field(descriptor.donor_block);
            validate_field(receiver, descriptor);
            validate_field(donor, descriptor);
            for (const auto& pair : descriptor.pairs) {
                store_flux(receiver,
                           descriptor.receiver_axis,
                           pair.receiver,
                           transform_inviscid_face_flux_for_receiver(
                               load_flux(donor, descriptor.donor_axis, pair.donor), descriptor));
            }
        }
    }
    for (const auto& pending : receives_) {
        validate_field(fields.field(pending.descriptor->receiver_block), *pending.descriptor);
    }
    for (auto& pending : sends_) {
        const auto& descriptor = *pending.descriptor;
        const auto& donor = fields.field(descriptor.donor_block);
        validate_field(donor, descriptor);
        pending.values[0] = static_cast<Real>(descriptor.version);
        std::size_t offset = 1;
        for (const auto& pair : descriptor.pairs) {
            const auto value = load_flux(donor, descriptor.donor_axis, pair.donor);
            for (const auto component : value)
                pending.values[offset++] = component;
        }
    }

#if WCNS_HAS_MPI
    std::fill(requests_.begin(), requests_.end(), MPI_REQUEST_NULL);
    std::size_t request = 0;
    for (auto& pending : receives_) {
        check_mpi(MPI_Irecv(pending.values.data(),
                            mpi_count(pending.values.size()),
                            MPI_DOUBLE,
                            pending.descriptor->donor_rank,
                            pending.descriptor->message_tag(),
                            mpi_.communicator(),
                            &requests_[request++]),
                  "MPI_Irecv face flux");
    }
    for (auto& pending : sends_) {
        check_mpi(MPI_Isend(pending.values.data(),
                            mpi_count(pending.values.size()),
                            MPI_DOUBLE,
                            pending.descriptor->receiver_rank,
                            pending.descriptor->message_tag(),
                            mpi_.communicator(),
                            &requests_[request++]),
                  "MPI_Isend face flux");
    }
    if (!requests_.empty()) {
        check_mpi(
            MPI_Waitall(static_cast<int>(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall face flux");
    }
#else
    if (!receives_.empty() || !sends_.empty()) {
        throw MpiError("remote face-flux exchange requires WCNS_ENABLE_MPI");
    }
#endif
    for (const auto& pending : receives_) {
        auto& receiver = fields.field(pending.descriptor->receiver_block);
        if (pending.values.empty()
            || pending.values[0] != static_cast<Real>(pending.descriptor->version)) {
            throw MpiError("face-flux message version mismatch");
        }
        std::size_t offset = 1;
        for (const auto& pair : pending.descriptor->pairs) {
            ConservativeState donor {};
            for (auto& component : donor)
                component = pending.values[offset++];
            store_flux(receiver,
                       pending.descriptor->receiver_axis,
                       pair.receiver,
                       transform_inviscid_face_flux_for_receiver(donor, *pending.descriptor));
        }
    }
}

void compute_inviscid_face_fluxes_into(InviscidFaceFluxField& result,
                                       const StructuredBlock& block,
                                       const MetricField& metric,
                                       const AlgorithmProfile& profile,
                                       const ReconstructionConfig& reconstruction,
                                       const RiemannSolver& riemann,
                                       const GasModel& gas,
                                       const ReferenceScales& reference,
                                       const NumericalFloors& floors,
                                       const BoundaryDataMap& boundary_data,
                                       const InviscidBoundaryOptions& boundary_options,
                                       std::uint64_t version,
                                       ReconstructionDiagnostics& diagnostics,
                                       RiemannDiagnostics* riemann_diagnostics,
                                       int rk_stage,
                                       Real stage_time,
                                       const FaceRobustnessField* robustness_levels,
                                       const RobustnessLadder* robustness_ladder,
                                       const RiemannSolver* robust_riemann)
{
    ProfileFactory::validate_bundle(profile.components());
    if (metric.profile() != profile.kind() || metric.dimension() != block.cell_dimension()) {
        throw ProfileError("inviscid flux metric belongs to another profile or dimension");
    }
    if (rk_stage < 0 || rk_stage > 3) {
        throw std::invalid_argument("inviscid flux RK stage must lie in [0,3]");
    }
    reconstruction.validate();
    const bool robustness_enabled = robustness_levels != nullptr;
    if (robustness_enabled != (robustness_ladder != nullptr)
        || robustness_enabled != (robust_riemann != nullptr)) {
        throw std::invalid_argument("inviscid flux robustness inputs must be supplied together");
    }
    if (robustness_enabled
        && (robustness_levels->profile() != profile.kind()
            || robustness_levels->dimension() != block.cell_dimension())) {
        throw ProfileError("inviscid flux robustness field has a mismatched profile");
    }
    const auto cells = block.cell_extent();
    if (result.profile() != profile.kind() || result.dimension() != block.cell_dimension()
        || result.field(Axis::I).interior_extent() != Extent3 {cells.ni + 1, cells.nj, cells.nk}
        || result.field(Axis::J).interior_extent() != Extent3 {cells.ni, cells.nj + 1, cells.nk}) {
        throw ProfileError("inviscid flux workspace metadata mismatch");
    }
    result.reset(version);
    const auto compute_axis = [&](Axis axis) {
        const auto& faces = metric_faces(metric, axis);
        const auto extent = faces.x.interior_extent();
        auto& output = result.field(axis);
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Index3 face {i, j, k};
                    // The non-owner copy is received from FaceFluxHaloExchanger.
                    // Skipping it here also keeps fallback diagnostics unique.
                    if (non_owned_connection_face(block, axis, face)) continue;
                    const FaceDiagnosticLocation diagnostic_location {
                        block.id(), block.owner_rank(), axis, face, version, rk_stage};
                    const int robustness_level
                        = robustness_enabled ? robustness_levels->level(axis, face) : 0;
                    const auto* face_reconstruction = &reconstruction;
                    const auto* face_riemann = &riemann;
                    if (robustness_enabled) {
                        const auto& strategy = robustness_ladder->strategy(robustness_level);
                        face_reconstruction = &strategy.reconstruction;
                        if (strategy.force_rusanov) face_riemann = robust_riemann;
                    }
                    Real area = 0.0;
                    const auto normal = unit_normal(faces, face, area);
                    auto states = reconstruct_thermodynamic_face(block.flow.conservative,
                                                                 block.flow.primitive,
                                                                 axis,
                                                                 face,
                                                                 *face_reconstruction,
                                                                 gas,
                                                                 reference,
                                                                 diagnostics,
                                                                 block.cell_dimension(),
                                                                 normal,
                                                                 diagnostic_location,
                                                                 !robustness_enabled);
                    if (const auto* patch = physical_patch(block, axis, face)) {
                        const auto data_iterator = boundary_data.find(patch->name);
                        if (data_iterator == boundary_data.end()) {
                            throw PhysicsConfigurationError(
                                "inviscid face boundary data is missing for patch " + patch->name);
                        }
                        if (patch->face.side == Side::Lower) {
                            const auto coordinates = boundary_face_coordinates(block, *patch, face);
                            states.left
                                = apply_inviscid_boundary_face_state(*patch,
                                                                     states.right,
                                                                     states.left,
                                                                     outward(normal, Side::Lower),
                                                                     data_iterator->second,
                                                                     boundary_options,
                                                                     gas,
                                                                     reference,
                                                                     floors,
                                                                     block.cell_dimension(),
                                                                     coordinates,
                                                                     stage_time);
                        } else {
                            const auto coordinates = boundary_face_coordinates(block, *patch, face);
                            states.right
                                = apply_inviscid_boundary_face_state(*patch,
                                                                     states.left,
                                                                     states.right,
                                                                     outward(normal, Side::Upper),
                                                                     data_iterator->second,
                                                                     boundary_options,
                                                                     gas,
                                                                     reference,
                                                                     floors,
                                                                     block.cell_dimension(),
                                                                     coordinates,
                                                                     stage_time);
                        }
                    }
                    const auto numerical
                        = face_riemann->solve(states.left, states.right, normal, gas, floors);
                    if (riemann_diagnostics != nullptr) {
                        riemann_diagnostics->record(numerical, diagnostic_location);
                    }
                    for (int component = 0; component < euler_components; ++component) {
                        output(i, j, k, component) = area
                            * numerical.flux_per_unit_area[static_cast<std::size_t>(component)];
                    }
                }
            }
        }
    };
    compute_axis(Axis::I);
    compute_axis(Axis::J);
    if (block.cell_dimension() == 3) compute_axis(Axis::K);
}

InviscidFaceFluxField compute_inviscid_face_fluxes(const StructuredBlock& block,
                                                   const MetricField& metric,
                                                   const AlgorithmProfile& profile,
                                                   const ReconstructionConfig& reconstruction,
                                                   const RiemannSolver& riemann,
                                                   const GasModel& gas,
                                                   const ReferenceScales& reference,
                                                   const NumericalFloors& floors,
                                                   const BoundaryDataMap& boundary_data,
                                                   const InviscidBoundaryOptions& boundary_options,
                                                   std::uint64_t version,
                                                   ReconstructionDiagnostics& diagnostics,
                                                   RiemannDiagnostics* riemann_diagnostics,
                                                   int rk_stage,
                                                   Real stage_time,
                                                   const FaceRobustnessField* robustness_levels,
                                                   const RobustnessLadder* robustness_ladder,
                                                   const RiemannSolver* robust_riemann)
{
    InviscidFaceFluxField result(
        block.cell_extent(), block.cell_dimension(), profile.kind(), version);
    compute_inviscid_face_fluxes_into(result,
                                      block,
                                      metric,
                                      profile,
                                      reconstruction,
                                      riemann,
                                      gas,
                                      reference,
                                      floors,
                                      boundary_data,
                                      boundary_options,
                                      version,
                                      diagnostics,
                                      riemann_diagnostics,
                                      rk_stage,
                                      stage_time,
                                      robustness_levels,
                                      robustness_ladder,
                                      robust_riemann);
    return result;
}

void compute_wcns_inviscid_residual(StructuredBlock& block,
                                    const MetricField& metric,
                                    const InviscidFaceFluxField& flux,
                                    const AlgorithmProfile& profile,
                                    FluxDifferenceMode mode)
{
    if (metric.profile() != profile.kind() || flux.profile() != profile.kind()
        || metric.dimension() != block.cell_dimension()
        || flux.dimension() != block.cell_dimension()) {
        throw ProfileError("WCNS residual inputs belong to different profiles or dimensions");
    }
    block.flow.residual.fill(0.0);
    const auto cells = block.cell_extent();
    static_cast<void>(flux_difference_mode_name(mode));
    const auto accumulate_axis = [&](Axis axis) {
        const auto& values = flux.field(axis);
        if (mode == FluxDifferenceMode::ConservativeTwoPoint) {
            for (int k = 0; k < cells.nk; ++k) {
                for (int j = 0; j < cells.nj; ++j) {
                    for (int i = 0; i < cells.ni; ++i) {
                        const Index3 cell {i, j, k};
                        const int normal = cell[static_cast<std::size_t>(axis)];
                        auto upper = cell;
                        upper[static_cast<std::size_t>(axis)] = normal + 1;
                        auto lower = cell;
                        lower[static_cast<std::size_t>(axis)] = normal;
                        const Real jacobian = metric.jacobian()(i, j, k);
                        if (!std::isfinite(jacobian) || jacobian <= 0.0) {
                            throw PhysicsError("WCNS flux divergence has an invalid Jacobian");
                        }
                        for (int component = 0; component < euler_components; ++component) {
                            const Real derivative = values(upper.i, upper.j, upper.k, component)
                                - values(lower.i, lower.j, lower.k, component);
                            if (!std::isfinite(derivative)) {
                                throw PhysicsError("WCNS flux divergence is non-finite");
                            }
                            block.flow.residual(i, j, k, component) -= derivative / jacobian;
                        }
                    }
                }
            }
            return;
        }

        const int count = cells[static_cast<std::size_t>(axis)];
        const auto& rows = cached_line_operators(profile, count).derivative_rows();
        const int boundary_width = profile.kind() == AlgorithmProfileKind::PhengleiWcns ? 1 : 2;
        constexpr std::array<int, 4> ph_offsets {{-1, 0, 1, 2}};
        constexpr std::array<Real, 4> ph_coefficients {
            {1.0 / 24.0, -27.0 / 24.0, 27.0 / 24.0, -1.0 / 24.0}};
        constexpr std::array<int, 6> scmm_offsets {{-2, -1, 0, 1, 2, 3}};
        constexpr std::array<Real, 6> scmm_coefficients {{-9.0 / 1920.0,
                                                          125.0 / 1920.0,
                                                          -2250.0 / 1920.0,
                                                          2250.0 / 1920.0,
                                                          -125.0 / 1920.0,
                                                          9.0 / 1920.0}};
        for (int k = 0; k < cells.nk; ++k) {
            for (int j = 0; j < cells.nj; ++j) {
                for (int i = 0; i < cells.ni; ++i) {
                    const Index3 cell {i, j, k};
                    const int normal = cell[static_cast<std::size_t>(axis)];
                    Index3 lower_face = cell;
                    lower_face[static_cast<std::size_t>(axis)] = 0;
                    Index3 upper_face = cell;
                    upper_face[static_cast<std::size_t>(axis)] = count;
                    const bool use_centered_connection_stencil
                        = (normal < boundary_width
                           && connection_covers(block, axis, Side::Lower, lower_face))
                        || (normal >= count - boundary_width
                            && connection_covers(block, axis, Side::Upper, upper_face));
                    for (int component = 0; component < euler_components; ++component) {
                        Real derivative = 0.0;
                        if (use_centered_connection_stencil) {
                            if (profile.kind() == AlgorithmProfileKind::PhengleiWcns) {
                                for (std::size_t entry = 0; entry < ph_offsets.size(); ++entry) {
                                    auto face = cell;
                                    face[static_cast<std::size_t>(axis)]
                                        = normal + ph_offsets[entry];
                                    derivative += ph_coefficients[entry]
                                        * values(face.i, face.j, face.k, component);
                                }
                            } else {
                                for (std::size_t entry = 0; entry < scmm_offsets.size(); ++entry) {
                                    auto face = cell;
                                    face[static_cast<std::size_t>(axis)]
                                        = normal + scmm_offsets[entry];
                                    derivative += scmm_coefficients[entry]
                                        * values(face.i, face.j, face.k, component);
                                }
                            }
                        } else {
                            const auto& row = rows[static_cast<std::size_t>(normal)];
                            for (const auto [face_index, coefficient] : row) {
                                auto face = cell;
                                face[static_cast<std::size_t>(axis)] = face_index;
                                derivative
                                    += coefficient * values(face.i, face.j, face.k, component);
                            }
                        }
                        const Real jacobian = metric.jacobian()(i, j, k);
                        if (!std::isfinite(derivative) || !std::isfinite(jacobian)
                            || jacobian <= 0.0) {
                            throw PhysicsError(
                                "WCNS flux divergence is non-finite or has invalid Jacobian");
                        }
                        block.flow.residual(i, j, k, component) -= derivative / jacobian;
                    }
                }
            }
        }
    };
    accumulate_axis(Axis::I);
    accumulate_axis(Axis::J);
    if (block.cell_dimension() == 3) accumulate_axis(Axis::K);
}

} // namespace wcns
