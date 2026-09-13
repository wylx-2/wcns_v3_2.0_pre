#pragma once

#include <wcns/core/field.hpp>
#include <wcns/mesh/structured_block.hpp>
#include <wcns/parallel/mpi_runtime.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace wcns {

using WallPoint = std::array<Real, 3>;

enum class WallPrimitiveKind {
    Segment2D,
    Triangle3D,
};

struct WallPrimitive {
    WallPrimitiveKind kind = WallPrimitiveKind::Segment2D;
    std::array<WallPoint, 3> vertices {};

    [[nodiscard]] static WallPrimitive segment(WallPoint first, WallPoint second);
    [[nodiscard]] static WallPrimitive
    triangle(WallPoint first, WallPoint second, WallPoint third);
    void validate() const;
};

struct WallDistanceResult {
    Real distance = 0.0;
    std::size_t primitive = 0;
};

class WallDistanceIndex {
public:
    explicit WallDistanceIndex(std::vector<WallPrimitive> primitives);

    [[nodiscard]] const std::vector<WallPrimitive>& primitives() const noexcept
    {
        return primitives_;
    }
    [[nodiscard]] WallDistanceResult query(WallPoint point) const;
    void fill_cell_field(const StructuredBlock& block,
                         Field<Real>& distance,
                         int component = 0) const;

private:
    struct Node {
        WallPoint lower {};
        WallPoint upper {};
        std::size_t begin = 0;
        std::size_t end = 0;
        int left = -1;
        int right = -1;
    };

    int build_node(std::size_t begin, std::size_t end);
    void query_node(int node,
                    WallPoint point,
                    Real& best_squared,
                    std::size_t& best_primitive) const;

    std::vector<WallPrimitive> primitives_;
    std::vector<std::size_t> order_;
    std::vector<Node> nodes_;
};

[[nodiscard]] std::vector<WallPrimitive>
extract_wall_primitives(const std::vector<StructuredBlock>& local_blocks);
[[nodiscard]] std::vector<WallPrimitive>
collect_global_wall_primitives(const MpiRuntime& mpi,
                               const std::vector<WallPrimitive>& local_primitives);

} // namespace wcns
