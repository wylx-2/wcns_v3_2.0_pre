#pragma once

#include <wcns/core/index.hpp>
#include <wcns/core/types.hpp>
#include <wcns/mesh/geometry.hpp>
#include <wcns/mesh/topology.hpp>
#include <wcns/solver/flow_fields.hpp>
#include <wcns/solver/turbulence_fields.hpp>

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace wcns {

class StructuredBlock {
private:
    BlockId id_;
    std::string name_;
    RankId owner_rank_;
    int cell_dimension_;
    int physical_dimension_;
    Extent3 vertex_extent_;
    Extent3 cell_extent_;
    int ghost_width_;

public:
    StructuredBlock(BlockId id,
                    std::string name,
                    RankId owner_rank,
                    int cell_dimension,
                    int physical_dimension,
                    Extent3 vertex_extent,
                    int ghost_width)
        : id_(checked_id(id))
        , name_(checked_name(std::move(name)))
        , owner_rank_(checked_rank(owner_rank))
        , cell_dimension_(checked_dimensions(cell_dimension, physical_dimension))
        , physical_dimension_(physical_dimension)
        , vertex_extent_(checked_vertex_extent(vertex_extent, cell_dimension_))
        , cell_extent_(make_cell_extent(vertex_extent_, cell_dimension_))
        , ghost_width_(checked_ghost_width(ghost_width))
        , coordinates(vertex_extent_)
        , cell_metrics(cell_extent_, ghost_width_)
        , face_metrics(cell_extent_)
        , flow(cell_extent_, ghost_width_)
    { }

    [[nodiscard]] BlockId id() const noexcept { return id_; }
    [[nodiscard]] const std::string& name() const noexcept { return name_; }
    [[nodiscard]] RankId owner_rank() const noexcept { return owner_rank_; }
    [[nodiscard]] int cell_dimension() const noexcept { return cell_dimension_; }
    [[nodiscard]] int physical_dimension() const noexcept { return physical_dimension_; }
    [[nodiscard]] const Extent3& vertex_extent() const noexcept { return vertex_extent_; }
    [[nodiscard]] const Extent3& cell_extent() const noexcept { return cell_extent_; }
    [[nodiscard]] int ghost_width() const noexcept { return ghost_width_; }

    void set_owner_rank(RankId owner_rank) { owner_rank_ = checked_rank(owner_rank); }

    NodeCoordinates coordinates;
    CellMetrics cell_metrics;
    FaceMetrics face_metrics;
    std::vector<BoundaryPatch> boundaries;
    std::vector<ConnectivityPatch> connectivities;
    FlowFields flow;
    TurbulenceFieldSet turbulence;

private:
    static BlockId checked_id(BlockId id)
    {
        if (id < 0) {
            throw std::invalid_argument("StructuredBlock id must be non-negative");
        }
        return id;
    }

    static std::string checked_name(std::string name)
    {
        if (name.empty()) {
            throw std::invalid_argument("StructuredBlock name must not be empty");
        }
        return name;
    }

    static RankId checked_rank(RankId rank)
    {
        if (rank < invalid_rank_id) {
            throw std::invalid_argument("StructuredBlock owner rank is invalid");
        }
        return rank;
    }

    static int checked_dimensions(int cell_dimension, int physical_dimension)
    {
        if ((cell_dimension != 2 && cell_dimension != 3)
            || (physical_dimension != 2 && physical_dimension != 3)
            || cell_dimension > physical_dimension) {
            throw std::invalid_argument(
                "StructuredBlock dimensions must satisfy 2 <= cell <= physical <= 3");
        }
        return cell_dimension;
    }

    static Extent3 checked_vertex_extent(Extent3 extent, int cell_dimension)
    {
        const bool valid_2d
            = cell_dimension == 2 && extent.ni >= 2 && extent.nj >= 2 && extent.nk == 1;
        const bool valid_3d
            = cell_dimension == 3 && extent.ni >= 2 && extent.nj >= 2 && extent.nk >= 2;
        if (!valid_2d && !valid_3d) {
            throw std::invalid_argument(
                "StructuredBlock vertex extent is incompatible with cell dimension");
        }
        return extent;
    }

    static Extent3 make_cell_extent(Extent3 vertex_extent, int cell_dimension)
    {
        return {
            vertex_extent.ni - 1,
            vertex_extent.nj - 1,
            cell_dimension == 2 ? 1 : vertex_extent.nk - 1,
        };
    }

    static int checked_ghost_width(int ghost_width)
    {
        if (ghost_width < 0) {
            throw std::invalid_argument("StructuredBlock ghost width must be non-negative");
        }
        return ghost_width;
    }
};

} // namespace wcns
