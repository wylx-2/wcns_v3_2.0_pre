#include <wcns/mesh/metrics.hpp>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>

namespace wcns {
namespace {

struct Vector3 {
    Real x = 0.0;
    Real y = 0.0;
    Real z = 0.0;
};

Vector3 operator+(Vector3 lhs, Vector3 rhs)
{
    return {lhs.x + rhs.x, lhs.y + rhs.y, lhs.z + rhs.z};
}

Vector3 operator-(Vector3 lhs, Vector3 rhs)
{
    return {lhs.x - rhs.x, lhs.y - rhs.y, lhs.z - rhs.z};
}

Vector3 operator*(Real scale, Vector3 value)
{
    return {scale * value.x, scale * value.y, scale * value.z};
}

Real dot(Vector3 lhs, Vector3 rhs)
{
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z;
}

Vector3 cross(Vector3 lhs, Vector3 rhs)
{
    return {
        lhs.y * rhs.z - lhs.z * rhs.y,
        lhs.z * rhs.x - lhs.x * rhs.z,
        lhs.x * rhs.y - lhs.y * rhs.x,
    };
}

Real magnitude(Vector3 value)
{
    return std::sqrt(dot(value, value));
}

Vector3 average(std::initializer_list<Vector3> values)
{
    Vector3 result;
    for (const auto value : values) {
        result = result + value;
    }
    return (1.0 / static_cast<Real>(values.size())) * result;
}

Vector3 node(const StructuredBlock& block, int i, int j, int k)
{
    return {
        block.coordinates.x(i, j, k),
        block.coordinates.y(i, j, k),
        block.coordinates.z(i, j, k),
    };
}

void store_face(FaceMetric& metric, int i, int j, int k, Vector3 area_vector)
{
    const Real area = magnitude(area_vector);
    if (!(area > 0.0)) {
        throw GeometryError("encountered a degenerate structured-grid face");
    }
    metric.normal_x(i, j, k) = area_vector.x / area;
    metric.normal_y(i, j, k) = area_vector.y / area;
    metric.normal_z(i, j, k) = area_vector.z / area;
    metric.area(i, j, k) = area;
}

Vector3 area_vector(const FaceMetric& metric, int i, int j, int k)
{
    const Real area = metric.area(i, j, k);
    return {
        area * metric.normal_x(i, j, k),
        area * metric.normal_y(i, j, k),
        area * metric.normal_z(i, j, k),
    };
}

Vector3 quadrilateral_area_vector(Vector3 p0, Vector3 p1, Vector3 p2, Vector3 p3)
{
    // Polygon area vector for the directed order p0 -> p1 -> p2 -> p3.
    return 0.5 * (cross(p0, p1) + cross(p1, p2) + cross(p2, p3) + cross(p3, p0));
}

void store_cell_center(StructuredBlock& block, int i, int j, int k, Vector3 center)
{
    block.cell_metrics.center_x(i, j, k) = center.x;
    block.cell_metrics.center_y(i, j, k) = center.y;
    block.cell_metrics.center_z(i, j, k) = center.z;
}

void compute_2d_metrics(StructuredBlock& block)
{
    if (block.physical_dimension() == 3) {
        // Some CGNS writers store planar 2D grids with XYZ coordinates.  The
        // existing 2D metric formulas remain valid only for a constant-Z plane.
        Real minimum_z = std::numeric_limits<Real>::infinity();
        Real maximum_z = -std::numeric_limits<Real>::infinity();
        const auto vertices = block.vertex_extent();
        for (int j = 0; j < vertices.nj; ++j) {
            for (int i = 0; i < vertices.ni; ++i) {
                const Real z = block.coordinates.z(i, j, 0);
                if (!std::isfinite(z)) {
                    throw GeometryError("2D embedded grid contains a non-finite Z coordinate");
                }
                minimum_z = std::min(minimum_z, z);
                maximum_z = std::max(maximum_z, z);
            }
        }
        const Real scale = std::max({Real {1}, std::abs(minimum_z), std::abs(maximum_z)});
        const Real tolerance = Real {512} * std::numeric_limits<Real>::epsilon() * scale;
        if (maximum_z - minimum_z > tolerance) {
            throw GeometryError(
                "2D physical-dimension-3 grid must lie in a constant-Z Cartesian plane");
        }
    }
    const auto cells = block.cell_extent();

    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i <= cells.ni; ++i) {
            const auto lower = node(block, i, j, 0);
            const auto upper = node(block, i, j + 1, 0);
            const auto tangent = upper - lower;
            store_face(block.face_metrics.i_faces, i, j, 0, {tangent.y, -tangent.x, 0.0});
        }
    }

    for (int j = 0; j <= cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            const auto left = node(block, i, j, 0);
            const auto right = node(block, i + 1, j, 0);
            const auto tangent = right - left;
            store_face(block.face_metrics.j_faces, i, j, 0, {-tangent.y, tangent.x, 0.0});
        }
    }

    for (int j = 0; j < cells.nj; ++j) {
        for (int i = 0; i < cells.ni; ++i) {
            const auto p00 = node(block, i, j, 0);
            const auto p10 = node(block, i + 1, j, 0);
            const auto p11 = node(block, i + 1, j + 1, 0);
            const auto p01 = node(block, i, j + 1, 0);
            store_cell_center(block, i, j, 0, average({p00, p10, p11, p01}));
            const Real area = 0.5
                * (p00.x * p10.y - p00.y * p10.x + p10.x * p11.y - p10.y * p11.x + p11.x * p01.y
                   - p11.y * p01.x + p01.x * p00.y - p01.y * p00.x);
            if (!(area > 0.0)) {
                throw GeometryError("2D structured cell has non-positive signed area");
            }
            block.cell_metrics.volume(i, j, 0) = area;
            block.cell_metrics.jacobian(i, j, 0) = area;
        }
    }
}

void compute_3d_metrics(StructuredBlock& block)
{
    const auto cells = block.cell_extent();

    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i <= cells.ni; ++i) {
                const auto p0 = node(block, i, j, k);
                const auto p1 = node(block, i, j + 1, k);
                const auto p2 = node(block, i, j + 1, k + 1);
                const auto p3 = node(block, i, j, k + 1);
                store_face(
                    block.face_metrics.i_faces, i, j, k, quadrilateral_area_vector(p0, p1, p2, p3));
            }
        }
    }

    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j <= cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const auto p0 = node(block, i, j, k);
                const auto p1 = node(block, i, j, k + 1);
                const auto p2 = node(block, i + 1, j, k + 1);
                const auto p3 = node(block, i + 1, j, k);
                store_face(
                    block.face_metrics.j_faces, i, j, k, quadrilateral_area_vector(p0, p1, p2, p3));
            }
        }
    }

    for (int k = 0; k <= cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const auto p0 = node(block, i, j, k);
                const auto p1 = node(block, i + 1, j, k);
                const auto p2 = node(block, i + 1, j + 1, k);
                const auto p3 = node(block, i, j + 1, k);
                store_face(
                    block.face_metrics.k_faces, i, j, k, quadrilateral_area_vector(p0, p1, p2, p3));
            }
        }
    }

    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const auto p000 = node(block, i, j, k);
                const auto p100 = node(block, i + 1, j, k);
                const auto p110 = node(block, i + 1, j + 1, k);
                const auto p010 = node(block, i, j + 1, k);
                const auto p001 = node(block, i, j, k + 1);
                const auto p101 = node(block, i + 1, j, k + 1);
                const auto p111 = node(block, i + 1, j + 1, k + 1);
                const auto p011 = node(block, i, j + 1, k + 1);
                store_cell_center(
                    block, i, j, k, average({p000, p100, p110, p010, p001, p101, p111, p011}));

                const auto ci_lower = average({p000, p010, p011, p001});
                const auto ci_upper = average({p100, p110, p111, p101});
                const auto cj_lower = average({p000, p001, p101, p100});
                const auto cj_upper = average({p010, p011, p111, p110});
                const auto ck_lower = average({p000, p100, p110, p010});
                const auto ck_upper = average({p001, p101, p111, p011});

                const Real volume = (1.0 / 3.0)
                    * (dot(ci_upper, area_vector(block.face_metrics.i_faces, i + 1, j, k))
                       - dot(ci_lower, area_vector(block.face_metrics.i_faces, i, j, k))
                       + dot(cj_upper, area_vector(block.face_metrics.j_faces, i, j + 1, k))
                       - dot(cj_lower, area_vector(block.face_metrics.j_faces, i, j, k))
                       + dot(ck_upper, area_vector(block.face_metrics.k_faces, i, j, k + 1))
                       - dot(ck_lower, area_vector(block.face_metrics.k_faces, i, j, k)));
                if (!(volume > 0.0)) {
                    throw GeometryError("3D structured cell has non-positive volume");
                }
                block.cell_metrics.volume(i, j, k) = volume;
                block.cell_metrics.jacobian(i, j, k) = volume;
            }
        }
    }
}

} // namespace

void compute_metrics(StructuredBlock& block)
{
    if (block.cell_dimension() == 2) {
        compute_2d_metrics(block);
    } else if (block.cell_dimension() == 3) {
        compute_3d_metrics(block);
    } else {
        throw GeometryError("unsupported structured-grid cell dimension");
    }
}

} // namespace wcns
