#pragma once

#include <wcns/core/index.hpp>
#include <wcns/core/types.hpp>
#include <wcns/mesh/topology.hpp>

#include <string>
#include <vector>

namespace wcns {

struct HaloCellPair {
    Index3 receiver_ghost;
    Index3 donor_interior;

    friend bool operator==(const HaloCellPair& lhs, const HaloCellPair& rhs)
    {
        return lhs.receiver_ghost == rhs.receiver_ghost && lhs.donor_interior == rhs.donor_interior;
    }
};

struct HaloExchangePlan {
    std::string connectivity_name;
    BlockId receiver_block = invalid_block_id;
    BlockId donor_block = invalid_block_id;
    RankId donor_rank = invalid_rank_id;
    Extent3 receiver_extent {};
    Extent3 donor_extent {};
    int dimension = 0;
    IndexTransform transform {};
    std::vector<HaloCellPair> cell_pairs;
};

[[nodiscard]] HaloExchangePlan make_halo_exchange_plan(const ConnectivityPatch& connection,
                                                       Extent3 receiver_cell_extent,
                                                       Extent3 donor_cell_extent,
                                                       int dimension);

} // namespace wcns
