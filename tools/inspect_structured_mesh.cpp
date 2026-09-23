#include <wcns/io/cgns_reader.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using Point = std::array<double, 3>;

Point point(const wcns::StructuredBlock& block, wcns::Index3 index)
{
    return {
        block.coordinates.x(index.i, index.j, index.k),
        block.coordinates.y(index.i, index.j, index.k),
        block.coordinates.z(index.i, index.j, index.k),
    };
}

double distance(const Point& lhs, const Point& rhs)
{
    double squared = 0.0;
    for (std::size_t component = 0; component < lhs.size(); ++component) {
        const double delta = lhs[component] - rhs[component];
        squared += delta * delta;
    }
    return std::sqrt(squared);
}

const char* axis_name(wcns::Axis axis)
{
    if (axis == wcns::Axis::I) return "i";
    if (axis == wcns::Axis::J) return "j";
    return "k";
}

const char* side_name(wcns::Side side)
{
    return side == wcns::Side::Lower ? "lower" : "upper";
}

const char* boundary_name(wcns::BoundaryType type)
{
    using wcns::BoundaryType;
    switch (type) {
    case BoundaryType::Farfield: return "farfield";
    case BoundaryType::Inflow: return "inflow";
    case BoundaryType::Outflow: return "outflow";
    case BoundaryType::SlipWall: return "slip_wall";
    case BoundaryType::NoSlipAdiabaticWall: return "no_slip_adiabatic_wall";
    case BoundaryType::NoSlipIsothermalWall: return "no_slip_isothermal_wall";
    case BoundaryType::Symmetry: return "symmetry";
    case BoundaryType::Periodic: return "periodic";
    case BoundaryType::DoubleMachReflection: return "double_mach_reflection";
    case BoundaryType::Undefined: return "undefined";
    }
    return "undefined";
}

std::string index_text(wcns::Index3 index)
{
    return std::to_string(index.i) + ":" + std::to_string(index.j) + ":"
        + std::to_string(index.k);
}

void print_zone_geometry(const wcns::StructuredBlock& block)
{
    Point minimum {{std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity(),
                    std::numeric_limits<double>::infinity()}};
    Point maximum {{-std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity(),
                    -std::numeric_limits<double>::infinity()}};
    const auto vertices = block.vertex_extent();
    for (int k = 0; k < vertices.nk; ++k) {
        for (int j = 0; j < vertices.nj; ++j) {
            for (int i = 0; i < vertices.ni; ++i) {
                const auto value = point(block, {i, j, k});
                for (std::size_t component = 0; component < value.size(); ++component) {
                    minimum[component] = std::min(minimum[component], value[component]);
                    maximum[component] = std::max(maximum[component], value[component]);
                }
            }
        }
    }
    std::cout << "zone=" << block.name() << " id=" << block.id()
              << " cell_dimension=" << block.cell_dimension()
              << " physical_dimension=" << block.physical_dimension() << " vertices="
              << vertices.ni << ':' << vertices.nj << ':' << vertices.nk << " cells="
              << block.cell_extent().ni << ':' << block.cell_extent().nj << ':'
              << block.cell_extent().nk << '\n';
    std::cout << "coordinate_bounds x=" << minimum[0] << ':' << maximum[0]
              << " y=" << minimum[1] << ':' << maximum[1] << " z=" << minimum[2] << ':'
              << maximum[2] << '\n';

    if (block.cell_dimension() != 2) return;
    double minimum_area = std::numeric_limits<double>::infinity();
    double maximum_area = 0.0;
    double minimum_edge = std::numeric_limits<double>::infinity();
    double maximum_edge = 0.0;
    double maximum_edge_ratio = 0.0;
    std::size_t non_positive = 0;
    const auto cells = block.cell_extent();
    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            const std::array<Point, 4> corners {{point(block, {i, j, 0}),
                                                  point(block, {i + 1, j, 0}),
                                                  point(block, {i + 1, j + 1, 0}),
                                                  point(block, {i, j + 1, 0})}};
            double twice_area = 0.0;
            double local_minimum_edge = std::numeric_limits<double>::infinity();
            double local_maximum_edge = 0.0;
            for (std::size_t edge = 0; edge < corners.size(); ++edge) {
                const auto next = (edge + 1) % corners.size();
                twice_area += corners[edge][0] * corners[next][1]
                    - corners[next][0] * corners[edge][1];
                const double length = distance(corners[edge], corners[next]);
                local_minimum_edge = std::min(local_minimum_edge, length);
                local_maximum_edge = std::max(local_maximum_edge, length);
            }
            const double signed_area = 0.5 * twice_area;
            if (!(signed_area > 0.0)) ++non_positive;
            const double area = std::abs(signed_area);
            minimum_area = std::min(minimum_area, area);
            maximum_area = std::max(maximum_area, area);
            minimum_edge = std::min(minimum_edge, local_minimum_edge);
            maximum_edge = std::max(maximum_edge, local_maximum_edge);
            maximum_edge_ratio
                = std::max(maximum_edge_ratio, local_maximum_edge / local_minimum_edge);
        }
    }
    std::cout << "quad_quality area_min=" << minimum_area << " area_max=" << maximum_area
              << " edge_min=" << minimum_edge << " edge_max=" << maximum_edge
              << " edge_ratio_max=" << maximum_edge_ratio
              << " non_positive_cells=" << non_positive << '\n';
}

void print_wall_geometry(const wcns::StructuredBlock& block, const wcns::BoundaryPatch& wall)
{
    if (block.cell_dimension() != 2 || wall.face.axis == wcns::Axis::K) {
        throw std::runtime_error(
            "wall geometry inspection currently requires a two-dimensional face");
    }
    const auto normal = static_cast<std::size_t>(wall.face.axis);
    const auto tangent = normal == 0 ? std::size_t {1} : std::size_t {0};
    const auto begin = wall.vertex_range.begin;
    const auto end = wall.vertex_range.end;
    const int tangent_step = end[tangent] >= begin[tangent] ? 1 : -1;
    const int segment_count = std::abs(end[tangent] - begin[tangent]);
    if (segment_count <= 0) throw std::runtime_error("wall patch has no line segments");
    const int inward_step = wall.face.side == wcns::Side::Lower ? 1 : -1;

    double x_minimum = std::numeric_limits<double>::infinity();
    double x_maximum = -std::numeric_limits<double>::infinity();
    double surface_edge_minimum = std::numeric_limits<double>::infinity();
    double surface_edge_maximum = 0.0;
    double first_layer_minimum = std::numeric_limits<double>::infinity();
    double first_layer_maximum = 0.0;
    double first_layer_sum = 0.0;
    double normal_alignment_minimum = 1.0;
    Point leading_edge {};
    int leading_edge_offset = 0;
    const auto start_point = point(block, begin);
    const auto end_point = point(block, end);

    for (int offset = 0; offset <= segment_count; ++offset) {
        auto vertex = begin;
        vertex[tangent] += tangent_step * offset;
        const auto value = point(block, vertex);
        if (value[0] < x_minimum) {
            x_minimum = value[0];
            leading_edge = value;
            leading_edge_offset = offset;
        }
        x_maximum = std::max(x_maximum, value[0]);
        if (offset == segment_count) continue;

        auto next = vertex;
        next[tangent] += tangent_step;
        auto off_vertex = vertex;
        auto off_next = next;
        off_vertex[normal] += inward_step;
        off_next[normal] += inward_step;
        const auto current = value;
        const auto following = point(block, next);
        const auto current_off = point(block, off_vertex);
        const auto following_off = point(block, off_next);
        const Point tangent_vector {{following[0] - current[0],
                                     following[1] - current[1],
                                     following[2] - current[2]}};
        const Point normal_vector {{0.5 * (current_off[0] + following_off[0] - current[0]
                                           - following[0]),
                                    0.5 * (current_off[1] + following_off[1] - current[1]
                                           - following[1]),
                                    0.5 * (current_off[2] + following_off[2] - current[2]
                                           - following[2])}};
        const double surface_edge = distance(current, following);
        const double first_layer
            = std::sqrt(normal_vector[0] * normal_vector[0]
                        + normal_vector[1] * normal_vector[1]
                        + normal_vector[2] * normal_vector[2]);
        const double cross = std::abs(tangent_vector[0] * normal_vector[1]
                                      - tangent_vector[1] * normal_vector[0]);
        const double alignment = cross / (surface_edge * first_layer);
        surface_edge_minimum = std::min(surface_edge_minimum, surface_edge);
        surface_edge_maximum = std::max(surface_edge_maximum, surface_edge);
        first_layer_minimum = std::min(first_layer_minimum, first_layer);
        first_layer_maximum = std::max(first_layer_maximum, first_layer);
        first_layer_sum += first_layer;
        normal_alignment_minimum = std::min(normal_alignment_minimum, alignment);
    }

    std::cout << "wall_geometry name=" << wall.name << " segments=" << segment_count
              << " chord_x=" << x_maximum - x_minimum << " x_min=" << x_minimum
              << " x_max=" << x_maximum << " leading_edge_offset=" << leading_edge_offset
              << " leading_edge=" << leading_edge[0] << ':' << leading_edge[1]
              << " endpoint_gap=" << distance(start_point, end_point)
              << " surface_edge_min=" << surface_edge_minimum
              << " surface_edge_max=" << surface_edge_maximum
              << " first_layer_min=" << first_layer_minimum
              << " first_layer_mean=" << first_layer_sum / segment_count
              << " first_layer_max=" << first_layer_maximum
              << " normal_alignment_min=" << normal_alignment_minimum << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 3) {
        std::cerr << "usage: wcns_inspect_structured_mesh <mesh.cgns> [wall-patch]\n";
        return 1;
    }
    try {
        const std::string mesh_path = argv[1];
        const std::string requested_wall = argc == 3 ? argv[2] : "";
        wcns::CgnsReader reader;
        const auto metadata = reader.read_metadata(mesh_path);
        std::cout << std::setprecision(17) << "check=structured_mesh_intake bases="
                  << metadata.bases.size() << " zones=" << metadata.zones.size() << '\n';
        bool wall_found = requested_wall.empty();
        for (const auto& zone : metadata.zones) {
            const auto block = reader.read_block(mesh_path, zone, 0, 0);
            print_zone_geometry(block);
            for (const auto& boundary : block.boundaries) {
                std::cout << "boundary name=" << boundary.name
                          << " type=" << boundary_name(boundary.type)
                          << " face=" << axis_name(boundary.face.axis) << ':'
                          << side_name(boundary.face.side)
                          << " vertices=" << index_text(boundary.vertex_range.begin) << ':'
                          << index_text(boundary.vertex_range.end) << '\n';
                const bool selected = requested_wall.empty()
                    ? boundary.type == wcns::BoundaryType::SlipWall
                        || boundary.type == wcns::BoundaryType::NoSlipAdiabaticWall
                        || boundary.type == wcns::BoundaryType::NoSlipIsothermalWall
                    : boundary.name == requested_wall;
                if (selected) {
                    print_wall_geometry(block, boundary);
                    if (!requested_wall.empty()) wall_found = true;
                }
            }
            for (const auto& connection : block.connectivities) {
                std::cout << "connectivity name=" << connection.name
                          << " receiver=" << connection.receiver_block
                          << " donor=" << connection.donor_block
                          << " receiver_face=" << axis_name(connection.receiver_face.axis) << ':'
                          << side_name(connection.receiver_face.side)
                          << " donor_face=" << axis_name(connection.donor_face.axis) << ':'
                          << side_name(connection.donor_face.side)
                          << " receiver_vertices="
                          << index_text(connection.receiver_vertex_range.begin) << ':'
                          << index_text(connection.receiver_vertex_range.end)
                          << " donor_vertices=" << index_text(connection.donor_vertex_range.begin)
                          << ':' << index_text(connection.donor_vertex_range.end) << '\n';
            }
        }
        if (!wall_found) throw std::runtime_error("requested wall patch was not found");
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "structured mesh inspection failed: " << error.what() << '\n';
        return 1;
    }
}
