#include <wcns/parallel/halo_exchanger.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

namespace wcns {
namespace {

std::size_t value_count(const DirectedExchange& exchange, int components)
{
    const auto cells = exchange.halo.cell_pairs.size();
    const auto component_count = static_cast<std::size_t>(components);
    if (cells > std::numeric_limits<std::size_t>::max() / component_count) {
        throw std::overflow_error("halo message value count exceeds size_t range");
    }
    return cells * component_count;
}

void validate_field_for_exchange(const Field<Real>& field,
                                 const DirectedExchange& exchange,
                                 bool receiver)
{
    const auto extent = field.interior_extent();
    const int ghost = field.ghost_width();
    for (const auto& pair : exchange.halo.cell_pairs) {
        const auto index = receiver ? pair.receiver_ghost : pair.donor_interior;
        const bool valid = receiver
            ? index.i >= -ghost && index.i < extent.ni + ghost && index.j >= -ghost
                && index.j < extent.nj + ghost && index.k >= -ghost && index.k < extent.nk + ghost
            : index.i >= 0 && index.i < extent.ni && index.j >= 0 && index.j < extent.nj
                && index.k >= 0 && index.k < extent.nk;
        if (!valid) {
            throw std::invalid_argument(receiver ? "halo receiver index is outside field storage"
                                                 : "halo donor index is outside field interior");
        }
    }
}

void copy_local(const DirectedExchange& exchange, const BlockFieldRegistry& fields)
{
    auto& receiver = fields.field(exchange.halo.receiver_block);
    const auto& donor = fields.field(exchange.halo.donor_block);
    validate_field_for_exchange(receiver, exchange, true);
    validate_field_for_exchange(donor, exchange, false);
    for (const auto& pair : exchange.halo.cell_pairs) {
        for (int component = 0; component < fields.components(); ++component) {
            receiver(pair.receiver_ghost.i, pair.receiver_ghost.j, pair.receiver_ghost.k, component)
                = donor(
                    pair.donor_interior.i, pair.donor_interior.j, pair.donor_interior.k, component);
        }
    }
}

void pack_send(HaloMessageBuffer& pending, const BlockFieldRegistry& fields)
{
    const auto& donor = fields.field(pending.exchange->halo.donor_block);
    validate_field_for_exchange(donor, *pending.exchange, false);
    std::size_t output = 0;
    for (const auto& pair : pending.exchange->halo.cell_pairs) {
        for (int component = 0; component < fields.components(); ++component) {
            pending.values[output++] = donor(
                pair.donor_interior.i, pair.donor_interior.j, pair.donor_interior.k, component);
        }
    }
}

void unpack_receive(const HaloMessageBuffer& pending, const BlockFieldRegistry& fields)
{
    auto& receiver = fields.field(pending.exchange->halo.receiver_block);
    validate_field_for_exchange(receiver, *pending.exchange, true);
    std::size_t input = 0;
    for (const auto& pair : pending.exchange->halo.cell_pairs) {
        for (int component = 0; component < fields.components(); ++component) {
            receiver(pair.receiver_ghost.i, pair.receiver_ghost.j, pair.receiver_ghost.k, component)
                = pending.values[input++];
        }
    }
}

bool stored(Index3 index, Extent3 extent, int ghost, int dimension)
{
    for (int axis = 0; axis < dimension; ++axis) {
        const auto direction = static_cast<std::size_t>(axis);
        if (index[direction] < -ghost || index[direction] >= extent[direction] + ghost) {
            return false;
        }
    }
    return dimension == 3 || index.k == 0;
}

std::vector<HaloCellPair> tangential_pairs(const HaloExchangePlan& plan, int width)
{
    if (width <= 0 || plan.dimension < 2 || plan.dimension > 3
        || !plan.transform.valid(plan.dimension)) {
        throw std::invalid_argument("tangential halo metadata are invalid");
    }
    std::vector<HaloCellPair> result;
    for (const auto& base : plan.cell_pairs) {
        int receiver_normal = -1;
        for (int axis = 0; axis < plan.dimension; ++axis) {
            const auto direction = static_cast<std::size_t>(axis);
            if (base.receiver_ghost[direction] < 0
                || base.receiver_ghost[direction] >= plan.receiver_extent[direction]) {
                receiver_normal = axis;
                break;
            }
        }
        if (receiver_normal < 0) {
            throw std::logic_error("halo pair has no receiver-normal ghost coordinate");
        }
        for (int dk = plan.dimension == 3 ? -width : 0;
             dk <= (plan.dimension == 3 ? width : 0);
             ++dk) {
            for (int dj = -width; dj <= width; ++dj) {
                for (int di = -width; di <= width; ++di) {
                    const std::array<int, 3> offset {{di, dj, dk}};
                    if (offset[static_cast<std::size_t>(receiver_normal)] != 0) continue;
                    bool tangential_ghost = false;
                    Index3 receiver = base.receiver_ghost;
                    Index3 donor = base.donor_interior;
                    for (int receiver_axis = 0; receiver_axis < plan.dimension;
                         ++receiver_axis) {
                        const auto receiver_direction
                            = static_cast<std::size_t>(receiver_axis);
                        const int delta = offset[receiver_direction];
                        receiver[receiver_direction] += delta;
                        const int entry
                            = plan.transform.receiver_to_donor[receiver_direction];
                        const auto donor_direction
                            = static_cast<std::size_t>(std::abs(entry) - 1);
                        donor[donor_direction] += (entry < 0 ? -delta : delta);
                        if (receiver_axis != receiver_normal
                            && (receiver[receiver_direction] < 0
                                || receiver[receiver_direction]
                                    >= plan.receiver_extent[receiver_direction])) {
                            tangential_ghost = true;
                        }
                    }
                    if (!tangential_ghost
                        || !stored(receiver, plan.receiver_extent, width, plan.dimension)
                        || !stored(donor, plan.donor_extent, width, plan.dimension)) {
                        continue;
                    }
                    result.push_back({receiver, donor});
                }
            }
        }
    }
    std::sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs) {
        return std::tie(lhs.receiver_ghost.i,
                        lhs.receiver_ghost.j,
                        lhs.receiver_ghost.k,
                        lhs.donor_interior.i,
                        lhs.donor_interior.j,
                        lhs.donor_interior.k)
            < std::tie(rhs.receiver_ghost.i,
                       rhs.receiver_ghost.j,
                       rhs.receiver_ghost.k,
                       rhs.donor_interior.i,
                       rhs.donor_interior.j,
                       rhs.donor_interior.k);
    });
    result.erase(std::unique(result.begin(), result.end()), result.end());
    return result;
}

struct TangentialBuffer {
    const DirectedExchange* exchange = nullptr;
    std::vector<HaloCellPair> pairs;
    std::vector<Real> values;
};

void copy_tangential(const TangentialBuffer& pending, const BlockFieldRegistry& fields)
{
    auto& receiver = fields.field(pending.exchange->halo.receiver_block);
    const auto& donor = fields.field(pending.exchange->halo.donor_block);
    for (const auto& pair : pending.pairs) {
        for (int component = 0; component < fields.components(); ++component) {
            const Real value = donor(
                pair.donor_interior.i, pair.donor_interior.j, pair.donor_interior.k, component);
            if (!std::isfinite(value)) {
                throw std::invalid_argument("tangential halo donor value is non-finite");
            }
            receiver(pair.receiver_ghost.i, pair.receiver_ghost.j, pair.receiver_ghost.k, component)
                = value;
        }
    }
}

void pack_tangential(TangentialBuffer& pending, const BlockFieldRegistry& fields)
{
    const auto& donor = fields.field(pending.exchange->halo.donor_block);
    std::size_t output = 0;
    for (const auto& pair : pending.pairs) {
        for (int component = 0; component < fields.components(); ++component) {
            const Real value = donor(
                pair.donor_interior.i, pair.donor_interior.j, pair.donor_interior.k, component);
            if (!std::isfinite(value)) {
                throw std::invalid_argument("tangential halo donor value is non-finite");
            }
            pending.values[output++] = value;
        }
    }
}

void unpack_tangential(const TangentialBuffer& pending, const BlockFieldRegistry& fields)
{
    auto& receiver = fields.field(pending.exchange->halo.receiver_block);
    std::size_t input = 0;
    for (const auto& pair : pending.pairs) {
        for (int component = 0; component < fields.components(); ++component) {
            receiver(pair.receiver_ghost.i, pair.receiver_ghost.j, pair.receiver_ghost.k, component)
                = pending.values[input++];
        }
    }
}

#if WCNS_HAS_MPI
int mpi_count(std::size_t values)
{
    if (values > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::overflow_error("halo message exceeds MPI int count range");
    }
    return static_cast<int>(values);
}
#endif

} // namespace

BlockFieldRegistry::BlockFieldRegistry(int components)
    : components_(components)
{
    if (components <= 0) {
        throw std::invalid_argument("field registry component count must be positive");
    }
}

void BlockFieldRegistry::add(BlockId block, Field<Real>& field)
{
    if (block < 0 || field.components() != components_) {
        throw std::invalid_argument("registered block field has incompatible metadata");
    }
    if (!fields_.emplace(block, &field).second) {
        throw std::invalid_argument("a block field is already registered");
    }
}

bool BlockFieldRegistry::contains(BlockId block) const noexcept
{
    return fields_.find(block) != fields_.end();
}

Field<Real>& BlockFieldRegistry::field(BlockId block) const
{
    const auto iterator = fields_.find(block);
    if (iterator == fields_.end()) {
        throw std::out_of_range("block field is not registered on this rank");
    }
    return *iterator->second;
}

HaloExchanger::HaloExchanger(const MpiRuntime& mpi,
                             const DistributedTopology& topology,
                             int distribution_rank_count,
                             int prepared_components)
    : mpi_(mpi)
    , topology_(topology)
{
    if (distribution_rank_count != mpi.size()) {
        throw std::invalid_argument("MPI size differs from the block distribution rank count");
    }
    if (prepared_components < 0) {
        throw std::invalid_argument("prepared halo component count must not be negative");
    }
    local_ = topology_.local_copies(mpi_.rank());
    receive_descriptors_ = topology_.receives(mpi_.rank());
    send_descriptors_ = topology_.sends(mpi_.rank());
    if (prepared_components > 0) prepare_buffers(prepared_components);
}

void HaloExchanger::prepare_buffers(int components) const
{
    receive_buffers_.clear();
    send_buffers_.clear();
    receive_buffers_.reserve(receive_descriptors_.size());
    send_buffers_.reserve(send_descriptors_.size());
    for (const auto* exchange : receive_descriptors_) {
        receive_buffers_.push_back(
            {exchange, std::vector<Real>(value_count(*exchange, components))});
    }
    for (const auto* exchange : send_descriptors_) {
        send_buffers_.push_back({exchange, std::vector<Real>(value_count(*exchange, components))});
    }
    prepared_components_ = components;
#if WCNS_HAS_MPI
    requests_.resize(receive_buffers_.size() + send_buffers_.size(), MPI_REQUEST_NULL);
#endif
}

void HaloExchanger::exchange(const BlockFieldRegistry& fields) const
{
    if (prepared_components_ != fields.components()) {
        prepare_buffers(fields.components());
    }
    for (const auto* exchange : local_) {
        copy_local(*exchange, fields);
    }

    for (const auto* exchange : receive_descriptors_) {
        if (!fields.contains(exchange->halo.receiver_block)) {
            throw std::invalid_argument("receiver field is missing on its owner rank");
        }
    }
    for (const auto* exchange : send_descriptors_) {
        if (!fields.contains(exchange->halo.donor_block)) {
            throw std::invalid_argument("donor field is missing on its owner rank");
        }
    }
    for (auto& pending : send_buffers_) {
        pack_send(pending, fields);
    }

#if WCNS_HAS_MPI
    int* tag_upper_bound = nullptr;
    int has_tag_upper_bound = 0;
    check_mpi(
        MPI_Comm_get_attr(mpi_.communicator(), MPI_TAG_UB, &tag_upper_bound, &has_tag_upper_bound),
        "MPI_Comm_get_attr MPI_TAG_UB");
    if (has_tag_upper_bound == 0 || tag_upper_bound == nullptr) {
        throw MpiError("MPI_TAG_UB is unavailable");
    }

    std::fill(requests_.begin(), requests_.end(), MPI_REQUEST_NULL);
    std::size_t request_index = 0;
    for (auto& pending : receive_buffers_) {
        const int tag = pending.exchange->message_tag();
        if (tag > *tag_upper_bound) {
            throw MpiError("halo receive tag exceeds MPI_TAG_UB");
        }
        check_mpi(MPI_Irecv(pending.values.data(),
                            mpi_count(pending.values.size()),
                            MPI_DOUBLE,
                            pending.exchange->donor_rank,
                            tag,
                            mpi_.communicator(),
                            &requests_[request_index++]),
                  "MPI_Irecv halo");
    }
    for (auto& pending : send_buffers_) {
        const int tag = pending.exchange->message_tag();
        if (tag > *tag_upper_bound) {
            throw MpiError("halo send tag exceeds MPI_TAG_UB");
        }
        check_mpi(MPI_Isend(pending.values.data(),
                            mpi_count(pending.values.size()),
                            MPI_DOUBLE,
                            pending.exchange->receiver_rank,
                            tag,
                            mpi_.communicator(),
                            &requests_[request_index++]),
                  "MPI_Isend halo");
    }
    if (!requests_.empty()) {
        check_mpi(
            MPI_Waitall(static_cast<int>(requests_.size()), requests_.data(), MPI_STATUSES_IGNORE),
            "MPI_Waitall halo");
    }
#else
    if (!receive_buffers_.empty() || !send_buffers_.empty()) {
        throw MpiError("remote halo exchange requires WCNS_ENABLE_MPI");
    }
#endif

    for (const auto& pending : receive_buffers_) {
        unpack_receive(pending, fields);
    }
}

void HaloExchanger::exchange_tangential_ghosts(const BlockFieldRegistry& fields,
                                               int tangential_width) const
{
    if (tangential_width <= 0) {
        throw std::invalid_argument("tangential halo width must be positive");
    }
    auto make_buffer = [&](const DirectedExchange* exchange) {
        auto pairs = tangential_pairs(exchange->halo, tangential_width);
        const auto count = pairs.size() * static_cast<std::size_t>(fields.components());
        return TangentialBuffer {exchange, std::move(pairs), std::vector<Real>(count)};
    };
    for (const auto* exchange : local_) {
        auto pending = make_buffer(exchange);
        copy_tangential(pending, fields);
    }
    std::vector<TangentialBuffer> receives;
    std::vector<TangentialBuffer> sends;
    receives.reserve(receive_descriptors_.size());
    sends.reserve(send_descriptors_.size());
    for (const auto* exchange : receive_descriptors_) receives.push_back(make_buffer(exchange));
    for (const auto* exchange : send_descriptors_) sends.push_back(make_buffer(exchange));
    for (auto& pending : sends) pack_tangential(pending, fields);
#if WCNS_HAS_MPI
    std::vector<MPI_Request> requests(receives.size() + sends.size(), MPI_REQUEST_NULL);
    std::size_t request = 0;
    for (auto& pending : receives) {
        check_mpi(MPI_Irecv(pending.values.data(),
                            mpi_count(pending.values.size()),
                            MPI_DOUBLE,
                            pending.exchange->donor_rank,
                            pending.exchange->message_tag(),
                            mpi_.communicator(),
                            &requests[request++]),
                  "MPI_Irecv tangential halo");
    }
    for (auto& pending : sends) {
        check_mpi(MPI_Isend(pending.values.data(),
                            mpi_count(pending.values.size()),
                            MPI_DOUBLE,
                            pending.exchange->receiver_rank,
                            pending.exchange->message_tag(),
                            mpi_.communicator(),
                            &requests[request++]),
                  "MPI_Isend tangential halo");
    }
    if (!requests.empty()) {
        check_mpi(MPI_Waitall(
                      static_cast<int>(requests.size()), requests.data(), MPI_STATUSES_IGNORE),
                  "MPI_Waitall tangential halo");
    }
#else
    if (!receives.empty() || !sends.empty()) {
        throw MpiError("remote tangential halo exchange requires WCNS_ENABLE_MPI");
    }
#endif
    for (const auto& pending : receives) unpack_tangential(pending, fields);
}

} // namespace wcns
