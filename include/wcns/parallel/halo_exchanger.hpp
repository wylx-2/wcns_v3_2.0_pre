#pragma once

#include <wcns/core/field.hpp>
#include <wcns/core/types.hpp>
#include <wcns/parallel/distributed_topology.hpp>
#include <wcns/parallel/mpi_runtime.hpp>

#include <unordered_map>

namespace wcns {

class BlockFieldRegistry {
public:
    explicit BlockFieldRegistry(int components);

    void add(BlockId block, Field<Real>& field);
    [[nodiscard]] bool contains(BlockId block) const noexcept;
    [[nodiscard]] Field<Real>& field(BlockId block) const;
    [[nodiscard]] int components() const noexcept { return components_; }

private:
    int components_ = 0;
    std::unordered_map<BlockId, Field<Real>*> fields_;
};

struct HaloMessageBuffer {
    const DirectedExchange* exchange = nullptr;
    std::vector<Real> values;
};

class HaloExchanger {
public:
    HaloExchanger(const MpiRuntime& mpi,
                  const DistributedTopology& topology,
                  int distribution_rank_count,
                  int prepared_components = 0);

    void exchange(const BlockFieldRegistry& fields) const;

    // Propagate already-filled face ghosts into connected edge/corner ghosts.
    // LES tensor filters call this after the ordinary face exchange and
    // physical-boundary extension.
    void exchange_tangential_ghosts(const BlockFieldRegistry& fields,
                                    int tangential_width) const;

private:
    void prepare_buffers(int components) const;

    const MpiRuntime& mpi_;
    const DistributedTopology& topology_;
    std::vector<const DirectedExchange*> local_;
    std::vector<const DirectedExchange*> receive_descriptors_;
    std::vector<const DirectedExchange*> send_descriptors_;
    mutable int prepared_components_ = 0;
    mutable std::vector<HaloMessageBuffer> receive_buffers_;
    mutable std::vector<HaloMessageBuffer> send_buffers_;
#if WCNS_HAS_MPI
    mutable std::vector<MPI_Request> requests_;
#endif
};

} // namespace wcns
