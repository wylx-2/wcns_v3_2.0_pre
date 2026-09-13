#include <wcns/solver/wall_distance.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>

namespace wcns {
namespace {

constexpr std::size_t serialized_primitive_size = 10;

WallPoint add(WallPoint lhs, WallPoint rhs)
{
    for (std::size_t axis = 0; axis < 3; ++axis)
        lhs[axis] += rhs[axis];
    return lhs;
}

WallPoint subtract(WallPoint lhs, WallPoint rhs)
{
    for (std::size_t axis = 0; axis < 3; ++axis)
        lhs[axis] -= rhs[axis];
    return lhs;
}

WallPoint multiply(Real scale, WallPoint value)
{
    for (auto& component : value)
        component *= scale;
    return value;
}

Real dot(WallPoint lhs, WallPoint rhs)
{
    Real result = 0.0;
    for (std::size_t axis = 0; axis < 3; ++axis)
        result += lhs[axis] * rhs[axis];
    return result;
}

WallPoint cross(WallPoint lhs, WallPoint rhs)
{
    return {{lhs[1] * rhs[2] - lhs[2] * rhs[1],
             lhs[2] * rhs[0] - lhs[0] * rhs[2],
             lhs[0] * rhs[1] - lhs[1] * rhs[0]}};
}

Real squared_norm(WallPoint value)
{
    return dot(value, value);
}

bool finite(WallPoint point)
{
    return std::all_of(point.begin(), point.end(), [](Real value) { return std::isfinite(value); });
}

void canonicalize(WallPrimitive& primitive)
{
    if (primitive.kind == WallPrimitiveKind::Segment2D) {
        if (primitive.vertices[1] < primitive.vertices[0]) {
            std::swap(primitive.vertices[0], primitive.vertices[1]);
        }
        primitive.vertices[2] = primitive.vertices[1];
    } else {
        std::sort(primitive.vertices.begin(), primitive.vertices.end());
    }
}

auto primitive_key(const WallPrimitive& primitive)
{
    return std::tuple {static_cast<int>(primitive.kind), primitive.vertices};
}

Real point_segment_squared(WallPoint point, const WallPrimitive& primitive)
{
    const auto a = primitive.vertices[0];
    const auto ab = subtract(primitive.vertices[1], a);
    const Real denominator = squared_norm(ab);
    const Real coordinate = std::clamp(dot(subtract(point, a), ab) / denominator, 0.0, 1.0);
    return squared_norm(subtract(point, add(a, multiply(coordinate, ab))));
}

// Closest-point region tests from the standard barycentric triangle algorithm.
Real point_triangle_squared(WallPoint point, const WallPrimitive& primitive)
{
    const auto a = primitive.vertices[0];
    const auto b = primitive.vertices[1];
    const auto c = primitive.vertices[2];
    const auto ab = subtract(b, a);
    const auto ac = subtract(c, a);
    const auto ap = subtract(point, a);
    const Real d1 = dot(ab, ap);
    const Real d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return squared_norm(ap);

    const auto bp = subtract(point, b);
    const Real d3 = dot(ab, bp);
    const Real d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return squared_norm(bp);

    const Real vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) {
        const Real coordinate = d1 / (d1 - d3);
        return squared_norm(subtract(point, add(a, multiply(coordinate, ab))));
    }

    const auto cp = subtract(point, c);
    const Real d5 = dot(ab, cp);
    const Real d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return squared_norm(cp);

    const Real vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) {
        const Real coordinate = d2 / (d2 - d6);
        return squared_norm(subtract(point, add(a, multiply(coordinate, ac))));
    }

    const Real va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && d4 - d3 >= 0.0 && d5 - d6 >= 0.0) {
        const auto bc = subtract(c, b);
        const Real coordinate = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return squared_norm(subtract(point, add(b, multiply(coordinate, bc))));
    }

    const Real denominator = 1.0 / (va + vb + vc);
    const Real v = vb * denominator;
    const Real w = vc * denominator;
    return squared_norm(subtract(point, add(a, add(multiply(v, ab), multiply(w, ac)))));
}

Real primitive_distance_squared(WallPoint point, const WallPrimitive& primitive)
{
    return primitive.kind == WallPrimitiveKind::Segment2D
        ? point_segment_squared(point, primitive)
        : point_triangle_squared(point, primitive);
}

Real box_distance_squared(WallPoint point, WallPoint lower, WallPoint upper)
{
    Real result = 0.0;
    for (std::size_t axis = 0; axis < 3; ++axis) {
        const Real delta = point[axis] < lower[axis]
            ? lower[axis] - point[axis]
            : (point[axis] > upper[axis] ? point[axis] - upper[axis] : 0.0);
        result += delta * delta;
    }
    return result;
}

WallPoint coordinate(const StructuredBlock& block, Index3 index)
{
    return {{block.coordinates.x(index.i, index.j, index.k),
             block.coordinates.y(index.i, index.j, index.k),
             block.coordinates.z(index.i, index.j, index.k)}};
}

bool is_viscous_wall(BoundaryType type)
{
    return type == BoundaryType::NoSlipAdiabaticWall
        || type == BoundaryType::NoSlipIsothermalWall;
}

std::vector<Real> serialize(const std::vector<WallPrimitive>& primitives)
{
    std::vector<Real> result;
    result.reserve(primitives.size() * serialized_primitive_size);
    for (const auto& primitive : primitives) {
        primitive.validate();
        result.push_back(primitive.kind == WallPrimitiveKind::Segment2D ? 2.0 : 3.0);
        for (const auto& vertex : primitive.vertices) {
            result.insert(result.end(), vertex.begin(), vertex.end());
        }
    }
    return result;
}

std::vector<WallPrimitive> deserialize(const std::vector<Real>& values)
{
    if (values.size() % serialized_primitive_size != 0) {
        throw std::runtime_error("global wall primitive payload is truncated");
    }
    std::vector<WallPrimitive> result;
    for (std::size_t offset = 0; offset < values.size(); offset += serialized_primitive_size) {
        WallPrimitive primitive;
        if (values[offset] == 2.0) {
            primitive.kind = WallPrimitiveKind::Segment2D;
        } else if (values[offset] == 3.0) {
            primitive.kind = WallPrimitiveKind::Triangle3D;
        } else {
            throw std::runtime_error("global wall primitive payload has an invalid kind");
        }
        std::size_t input = offset + 1;
        for (auto& vertex : primitive.vertices) {
            for (auto& component : vertex)
                component = values[input++];
        }
        canonicalize(primitive);
        primitive.validate();
        result.push_back(primitive);
    }
    std::sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs) {
        return primitive_key(lhs) < primitive_key(rhs);
    });
    result.erase(std::unique(result.begin(), result.end(), [](const auto& lhs, const auto& rhs) {
                     return primitive_key(lhs) == primitive_key(rhs);
                 }),
                 result.end());
    return result;
}

} // namespace

WallPrimitive WallPrimitive::segment(WallPoint first, WallPoint second)
{
    WallPrimitive result {WallPrimitiveKind::Segment2D, {first, second, second}};
    canonicalize(result);
    result.validate();
    return result;
}

WallPrimitive WallPrimitive::triangle(WallPoint first, WallPoint second, WallPoint third)
{
    WallPrimitive result {WallPrimitiveKind::Triangle3D, {first, second, third}};
    canonicalize(result);
    result.validate();
    return result;
}

void WallPrimitive::validate() const
{
    if (!std::all_of(vertices.begin(), vertices.end(), finite)) {
        throw std::invalid_argument("wall primitive contains non-finite coordinates");
    }
    if (kind == WallPrimitiveKind::Segment2D) {
        if (!(squared_norm(subtract(vertices[1], vertices[0])) > 0.0)
            || vertices[2] != vertices[1]) {
            throw std::invalid_argument("wall segment is degenerate or non-canonical");
        }
    } else if (!(squared_norm(cross(subtract(vertices[1], vertices[0]),
                                    subtract(vertices[2], vertices[0])))
                        > 0.0)) {
        throw std::invalid_argument("wall triangle is degenerate");
    }
}

WallDistanceIndex::WallDistanceIndex(std::vector<WallPrimitive> primitives)
    : primitives_(std::move(primitives))
{
    for (auto& primitive : primitives_) {
        canonicalize(primitive);
        primitive.validate();
    }
    std::sort(primitives_.begin(), primitives_.end(), [](const auto& lhs, const auto& rhs) {
        return primitive_key(lhs) < primitive_key(rhs);
    });
    primitives_.erase(
        std::unique(primitives_.begin(), primitives_.end(), [](const auto& lhs, const auto& rhs) {
            return primitive_key(lhs) == primitive_key(rhs);
        }),
        primitives_.end());
    if (primitives_.empty()) {
        throw std::invalid_argument("wall-distance index requires at least one wall primitive");
    }
    order_.resize(primitives_.size());
    std::iota(order_.begin(), order_.end(), std::size_t {0});
    nodes_.reserve(2 * primitives_.size());
    static_cast<void>(build_node(0, order_.size()));
}

int WallDistanceIndex::build_node(std::size_t begin, std::size_t end)
{
    Node node;
    node.begin = begin;
    node.end = end;
    node.lower.fill(std::numeric_limits<Real>::infinity());
    node.upper.fill(-std::numeric_limits<Real>::infinity());
    WallPoint centroid_lower = node.lower;
    WallPoint centroid_upper = node.upper;
    for (std::size_t position = begin; position < end; ++position) {
        const auto& primitive = primitives_[order_[position]];
        WallPoint centroid {{0.0, 0.0, 0.0}};
        const int count = primitive.kind == WallPrimitiveKind::Segment2D ? 2 : 3;
        for (int vertex = 0; vertex < count; ++vertex) {
            for (std::size_t axis = 0; axis < 3; ++axis) {
                const Real coordinate_value
                    = primitive.vertices[static_cast<std::size_t>(vertex)][axis];
                node.lower[axis] = std::min(node.lower[axis], coordinate_value);
                node.upper[axis] = std::max(node.upper[axis], coordinate_value);
                centroid[axis] += coordinate_value / static_cast<Real>(count);
            }
        }
        for (std::size_t axis = 0; axis < 3; ++axis) {
            centroid_lower[axis] = std::min(centroid_lower[axis], centroid[axis]);
            centroid_upper[axis] = std::max(centroid_upper[axis], centroid[axis]);
        }
    }
    const int node_index = static_cast<int>(nodes_.size());
    nodes_.push_back(node);
    if (end - begin <= 8) return node_index;

    std::size_t split_axis = 0;
    for (std::size_t axis = 1; axis < 3; ++axis) {
        if (centroid_upper[axis] - centroid_lower[axis]
            > centroid_upper[split_axis] - centroid_lower[split_axis]) {
            split_axis = axis;
        }
    }
    std::stable_sort(order_.begin() + static_cast<std::ptrdiff_t>(begin),
                     order_.begin() + static_cast<std::ptrdiff_t>(end),
                     [&](std::size_t lhs, std::size_t rhs) {
                         const auto centroid = [&](std::size_t primitive_index) {
                             const auto& primitive = primitives_[primitive_index];
                             const int count
                                 = primitive.kind == WallPrimitiveKind::Segment2D ? 2 : 3;
                             Real value = 0.0;
                             for (int vertex = 0; vertex < count; ++vertex) {
                                 value += primitive.vertices[static_cast<std::size_t>(vertex)]
                                                             [split_axis]
                                     / static_cast<Real>(count);
                             }
                             return value;
                         };
                         return std::tuple {centroid(lhs), lhs} < std::tuple {centroid(rhs), rhs};
                     });
    const auto middle = begin + (end - begin) / 2;
    nodes_[static_cast<std::size_t>(node_index)].left = build_node(begin, middle);
    nodes_[static_cast<std::size_t>(node_index)].right = build_node(middle, end);
    return node_index;
}

void WallDistanceIndex::query_node(int node_index,
                                   WallPoint point,
                                   Real& best_squared,
                                   std::size_t& best_primitive) const
{
    const auto& node = nodes_[static_cast<std::size_t>(node_index)];
    if (box_distance_squared(point, node.lower, node.upper) > best_squared) return;
    if (node.left < 0) {
        for (std::size_t position = node.begin; position < node.end; ++position) {
            const auto primitive_index = order_[position];
            const Real distance_squared
                = primitive_distance_squared(point, primitives_[primitive_index]);
            if (distance_squared < best_squared
                || (distance_squared == best_squared && primitive_index < best_primitive)) {
                best_squared = distance_squared;
                best_primitive = primitive_index;
            }
        }
        return;
    }
    const auto& left = nodes_[static_cast<std::size_t>(node.left)];
    const auto& right = nodes_[static_cast<std::size_t>(node.right)];
    const Real left_distance = box_distance_squared(point, left.lower, left.upper);
    const Real right_distance = box_distance_squared(point, right.lower, right.upper);
    const int first = left_distance <= right_distance ? node.left : node.right;
    const int second = left_distance <= right_distance ? node.right : node.left;
    query_node(first, point, best_squared, best_primitive);
    query_node(second, point, best_squared, best_primitive);
}

WallDistanceResult WallDistanceIndex::query(WallPoint point) const
{
    if (!finite(point)) throw std::invalid_argument("wall-distance query point is non-finite");
    Real best_squared = std::numeric_limits<Real>::infinity();
    std::size_t best_primitive = std::numeric_limits<std::size_t>::max();
    query_node(0, point, best_squared, best_primitive);
    if (!std::isfinite(best_squared) || best_primitive >= primitives_.size()) {
        throw std::runtime_error("wall-distance query did not find a primitive");
    }
    return {std::sqrt(std::max(0.0, best_squared)), best_primitive};
}

void WallDistanceIndex::fill_cell_field(const StructuredBlock& block,
                                        Field<Real>& distance,
                                        int component) const
{
    if (distance.interior_extent() != block.cell_extent() || component < 0
        || component >= distance.components()) {
        throw std::invalid_argument("wall-distance field metadata differs from the block");
    }
    const auto extent = block.cell_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                const WallPoint point {{block.cell_metrics.center_x(i, j, k),
                                        block.cell_metrics.center_y(i, j, k),
                                        block.cell_metrics.center_z(i, j, k)}};
                distance(i, j, k, component) = query(point).distance;
            }
        }
    }
}

std::vector<WallPrimitive>
extract_wall_primitives(const std::vector<StructuredBlock>& local_blocks)
{
    std::vector<WallPrimitive> result;
    for (const auto& block : local_blocks) {
        const int dimension = block.cell_dimension();
        for (const auto& patch : block.boundaries) {
            if (!is_viscous_wall(patch.type)) continue;
            const auto counts = patch.vertex_range.counts();
            if (dimension == 2) {
                const int normal = static_cast<int>(patch.face.axis);
                const int tangent = normal == 0 ? 1 : 0;
                const int count = counts[static_cast<std::size_t>(tangent)];
                for (int ordinal = 0; ordinal + 1 < count; ++ordinal) {
                    Index3 first_ordinal {};
                    Index3 second_ordinal {};
                    first_ordinal[static_cast<std::size_t>(tangent)] = ordinal;
                    second_ordinal[static_cast<std::size_t>(tangent)] = ordinal + 1;
                    result.push_back(WallPrimitive::segment(
                        coordinate(block, patch.vertex_range.at(first_ordinal)),
                        coordinate(block, patch.vertex_range.at(second_ordinal))));
                }
            } else {
                std::array<int, 2> tangents {};
                int output = 0;
                for (int axis = 0; axis < 3; ++axis) {
                    if (axis != static_cast<int>(patch.face.axis)) tangents[output++] = axis;
                }
                const int first_count = counts[static_cast<std::size_t>(tangents[0])];
                const int second_count = counts[static_cast<std::size_t>(tangents[1])];
                for (int second = 0; second + 1 < second_count; ++second) {
                    for (int first = 0; first + 1 < first_count; ++first) {
                        Index3 q00 {};
                        Index3 q10 {};
                        Index3 q11 {};
                        Index3 q01 {};
                        q00[static_cast<std::size_t>(tangents[0])] = first;
                        q00[static_cast<std::size_t>(tangents[1])] = second;
                        q10 = q00;
                        q10[static_cast<std::size_t>(tangents[0])] += 1;
                        q11 = q10;
                        q11[static_cast<std::size_t>(tangents[1])] += 1;
                        q01 = q00;
                        q01[static_cast<std::size_t>(tangents[1])] += 1;
                        const auto p00 = coordinate(block, patch.vertex_range.at(q00));
                        const auto p10 = coordinate(block, patch.vertex_range.at(q10));
                        const auto p11 = coordinate(block, patch.vertex_range.at(q11));
                        const auto p01 = coordinate(block, patch.vertex_range.at(q01));
                        result.push_back(WallPrimitive::triangle(p00, p10, p11));
                        result.push_back(WallPrimitive::triangle(p00, p11, p01));
                    }
                }
            }
        }
    }
    return deserialize(serialize(result));
}

std::vector<WallPrimitive>
collect_global_wall_primitives(const MpiRuntime& mpi,
                               const std::vector<WallPrimitive>& local_primitives)
{
    return deserialize(mpi.all_gather_reals(serialize(local_primitives)));
}

} // namespace wcns
