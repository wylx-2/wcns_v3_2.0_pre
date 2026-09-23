#include <wcns/mesh/high_order_metrics.hpp>

#include <wcns/mesh/linear_operators.hpp>
#include <wcns/mesh/metrics.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace wcns {

struct MetricFieldBuilderAccess {
    static MetricField create(AlgorithmProfileKind profile, Extent3 cell_extent, int dimension)
    {
        return MetricField(profile, cell_extent, dimension);
    }

    static CellCoordinates& cell_coordinates(MetricField& metric)
    {
        return metric.cell_coordinates_;
    }

    static Array3D<Real>& jacobian(MetricField& metric) { return metric.jacobian_; }

    static FaceAreaVectors& i_faces(MetricField& metric) { return metric.i_faces_; }
    static FaceAreaVectors& j_faces(MetricField& metric) { return metric.j_faces_; }
    static FaceAreaVectors& k_faces(MetricField& metric) { return metric.k_faces_; }
};

namespace {

using ScalarField = Array3D<Real>;

struct VectorFields {
    explicit VectorFields(Extent3 extent)
        : x(extent)
        , y(extent)
        , z(extent)
    { }

    ScalarField x;
    ScalarField y;
    ScalarField z;
};

Extent3 face_extent(Extent3 cells, int axis)
{
    ++cells[static_cast<std::size_t>(axis)];
    return cells;
}

Real finite_relative_difference(Real value, Real reference)
{
    return std::abs(value - reference) / reference;
}

ScalarField multiply(const ScalarField& lhs, const ScalarField& rhs)
{
    if (lhs.interior_extent() != rhs.interior_extent()) {
        throw GeometryError("cannot multiply geometry fields with different extents");
    }
    const auto extent = lhs.interior_extent();
    ScalarField result(extent);
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                result(i, j, k) = lhs(i, j, k) * rhs(i, j, k);
            }
        }
    }
    return result;
}

ScalarField combine(const ScalarField& a,
                    const ScalarField& b,
                    const ScalarField& c,
                    const ScalarField& d,
                    Real scale)
{
    const auto extent = a.interior_extent();
    ScalarField result(extent);
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                result(i, j, k) = scale * (a(i, j, k) + b(i, j, k) - c(i, j, k) - d(i, j, k));
            }
        }
    }
    return result;
}

template <class LineOperation>
ScalarField
apply_lines(const ScalarField& input, int axis, Extent3 output_extent, LineOperation operation)
{
    const auto input_extent = input.interior_extent();
    ScalarField output(output_extent);
    const int line_count = input_extent[static_cast<std::size_t>(axis)];

    const int outer_k = axis == 2 ? 1 : input_extent.nk;
    const int outer_j = axis == 1 ? 1 : input_extent.nj;
    const int outer_i = axis == 0 ? 1 : input_extent.ni;
    for (int k0 = 0; k0 < outer_k; ++k0) {
        for (int j0 = 0; j0 < outer_j; ++j0) {
            for (int i0 = 0; i0 < outer_i; ++i0) {
                std::vector<Real> line(static_cast<std::size_t>(line_count));
                for (int coordinate = 0; coordinate < line_count; ++coordinate) {
                    Index3 index {i0, j0, k0};
                    index[static_cast<std::size_t>(axis)] = coordinate;
                    line[static_cast<std::size_t>(coordinate)] = input(index.i, index.j, index.k);
                }
                const auto transformed = operation(line);
                const int output_count = output_extent[static_cast<std::size_t>(axis)];
                if (static_cast<int>(transformed.size()) != output_count) {
                    throw GeometryError("line operator returned an unexpected extent");
                }
                for (int coordinate = 0; coordinate < output_count; ++coordinate) {
                    Index3 index {i0, j0, k0};
                    index[static_cast<std::size_t>(axis)] = coordinate;
                    output(index.i, index.j, index.k)
                        = transformed[static_cast<std::size_t>(coordinate)];
                }
            }
        }
    }
    return output;
}

ScalarField refined_derivative(const ScalarField& field, int axis)
{
    return apply_lines(field, axis, field.interior_extent(), [](const std::vector<Real>& line) {
        return grid_delta(line);
    });
}

ScalarField
scmm_center_derivative(const ScalarField& field, int axis, const AlgorithmProfile& profile)
{
    const int count = field.interior_extent()[static_cast<std::size_t>(axis)];
    const auto operators = LineOperators::build(profile, count);
    return apply_lines(
        field, axis, field.interior_extent(), [&operators](const std::vector<Real>& line) {
            return operators.delta(line);
        });
}

ScalarField
scmm_center_to_faces(const ScalarField& field, int axis, const AlgorithmProfile& profile)
{
    const int count = field.interior_extent()[static_cast<std::size_t>(axis)];
    const auto operators = LineOperators::build(profile, count);
    auto output_extent = field.interior_extent();
    ++output_extent[static_cast<std::size_t>(axis)];
    return apply_lines(field, axis, output_extent, [&operators](const std::vector<Real>& line) {
        return operators.interpolate(line);
    });
}

ScalarField vertices_to_centers_axis(const ScalarField& field, int axis)
{
    auto output_extent = field.interior_extent();
    --output_extent[static_cast<std::size_t>(axis)];
    return apply_lines(field, axis, output_extent, [](const std::vector<Real>& line) {
        return interpolate_vertices_to_centers_i6(line);
    });
}

ScalarField
product_derivative(const ScalarField& first, const ScalarField& second_derivative, int outer_axis)
{
    return refined_derivative(multiply(first, second_derivative), outer_axis);
}

ScalarField symmetric_component(const ScalarField& first,
                                const ScalarField& second,
                                int first_tangent,
                                int second_tangent)
{
    const auto d_first_first = refined_derivative(first, first_tangent);
    const auto d_first_second = refined_derivative(first, second_tangent);
    const auto d_second_first = refined_derivative(second, first_tangent);
    const auto d_second_second = refined_derivative(second, second_tangent);
    const auto a = product_derivative(first, d_second_first, second_tangent);
    const auto b = product_derivative(second, d_first_second, first_tangent);
    const auto c = product_derivative(first, d_second_second, first_tangent);
    const auto d = product_derivative(second, d_first_first, second_tangent);
    return combine(a, b, c, d, 0.5);
}

ScalarField scmm_product_derivative(const ScalarField& first,
                                    const ScalarField& second_derivative,
                                    int outer_axis,
                                    const AlgorithmProfile& profile)
{
    return scmm_center_derivative(multiply(first, second_derivative), outer_axis, profile);
}

ScalarField scmm_symmetric_component(const ScalarField& first,
                                     const ScalarField& second,
                                     int first_tangent,
                                     int second_tangent,
                                     const AlgorithmProfile& profile)
{
    const auto d_first_first = scmm_center_derivative(first, first_tangent, profile);
    const auto d_first_second = scmm_center_derivative(first, second_tangent, profile);
    const auto d_second_first = scmm_center_derivative(second, first_tangent, profile);
    const auto d_second_second = scmm_center_derivative(second, second_tangent, profile);
    const auto a = scmm_product_derivative(first, d_second_first, second_tangent, profile);
    const auto b = scmm_product_derivative(second, d_first_second, first_tangent, profile);
    const auto c = scmm_product_derivative(first, d_second_second, first_tangent, profile);
    const auto d = scmm_product_derivative(second, d_first_first, second_tangent, profile);
    return combine(a, b, c, d, 0.5);
}

ScalarField dot_product(const VectorFields& coordinates, const VectorFields& vectors)
{
    const auto extent = coordinates.x.interior_extent();
    ScalarField result(extent);
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                result(i, j, k) = coordinates.x(i, j, k) * vectors.x(i, j, k)
                    + coordinates.y(i, j, k) * vectors.y(i, j, k)
                    + coordinates.z(i, j, k) * vectors.z(i, j, k);
            }
        }
    }
    return result;
}

VectorFields refined_coordinates(const StructuredBlock& block)
{
    const auto vertices = block.vertex_extent();
    const Extent3 refined {
        2 * vertices.ni - 1,
        2 * vertices.nj - 1,
        block.cell_dimension() == 3 ? 2 * vertices.nk - 1 : 1,
    };
    VectorFields result(refined);
    for (int rk = 0; rk < refined.nk; ++rk) {
        for (int rj = 0; rj < refined.nj; ++rj) {
            for (int ri = 0; ri < refined.ni; ++ri) {
                const int i0 = ri / 2;
                const int j0 = rj / 2;
                const int k0 = rk / 2;
                const int i_count = ri % 2 == 0 ? 1 : 2;
                const int j_count = rj % 2 == 0 ? 1 : 2;
                const int k_count = block.cell_dimension() == 2 || rk % 2 == 0 ? 1 : 2;
                Real x = 0.0;
                Real y = 0.0;
                Real z = 0.0;
                const Real weight = 1.0 / static_cast<Real>(i_count * j_count * k_count);
                for (int dk = 0; dk < k_count; ++dk) {
                    for (int dj = 0; dj < j_count; ++dj) {
                        for (int di = 0; di < i_count; ++di) {
                            x += weight * block.coordinates.x(i0 + di, j0 + dj, k0 + dk);
                            y += weight * block.coordinates.y(i0 + di, j0 + dj, k0 + dk);
                            z += weight * block.coordinates.z(i0 + di, j0 + dj, k0 + dk);
                        }
                    }
                }
                result.x(ri, rj, rk) = x;
                result.y(ri, rj, rk) = y;
                result.z(ri, rj, rk) = z;
            }
        }
    }
    return result;
}

void compute_refined_symmetric_metrics(const VectorFields& coordinates,
                                       int dimension,
                                       VectorFields& s_i,
                                       VectorFields& s_j,
                                       VectorFields& s_k,
                                       ScalarField& jacobian)
{
    if (dimension == 2) {
        s_i.x = refined_derivative(coordinates.y, 1);
        const auto dx_eta = refined_derivative(coordinates.x, 1);
        const auto dy_xi = refined_derivative(coordinates.y, 0);
        s_j.y = refined_derivative(coordinates.x, 0);
        const auto extent = coordinates.x.interior_extent();
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                s_i.y(i, j, 0) = -dx_eta(i, j, 0);
                s_i.z(i, j, 0) = 0.0;
                s_j.x(i, j, 0) = -dy_xi(i, j, 0);
                s_j.z(i, j, 0) = 0.0;
            }
        }
        const auto di = refined_derivative(dot_product(coordinates, s_i), 0);
        const auto dj = refined_derivative(dot_product(coordinates, s_j), 1);
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                jacobian(i, j, 0) = 0.5 * (di(i, j, 0) + dj(i, j, 0));
            }
        }
        return;
    }

    s_i.x = symmetric_component(coordinates.z, coordinates.y, 1, 2);
    s_i.y = symmetric_component(coordinates.x, coordinates.z, 1, 2);
    s_i.z = symmetric_component(coordinates.y, coordinates.x, 1, 2);
    s_j.x = symmetric_component(coordinates.z, coordinates.y, 2, 0);
    s_j.y = symmetric_component(coordinates.x, coordinates.z, 2, 0);
    s_j.z = symmetric_component(coordinates.y, coordinates.x, 2, 0);
    s_k.x = symmetric_component(coordinates.z, coordinates.y, 0, 1);
    s_k.y = symmetric_component(coordinates.x, coordinates.z, 0, 1);
    s_k.z = symmetric_component(coordinates.y, coordinates.x, 0, 1);

    const auto di = refined_derivative(dot_product(coordinates, s_i), 0);
    const auto dj = refined_derivative(dot_product(coordinates, s_j), 1);
    const auto dk = refined_derivative(dot_product(coordinates, s_k), 2);
    const auto extent = coordinates.x.interior_extent();
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                jacobian(i, j, k) = (di(i, j, k) + dj(i, j, k) + dk(i, j, k)) / 3.0;
            }
        }
    }
}

void copy_refined_to_metric(const VectorFields& coordinates,
                            const VectorFields& s_i,
                            const VectorFields& s_j,
                            const VectorFields& s_k,
                            const ScalarField& jacobian,
                            MetricField& metric)
{
    const auto cells = metric.jacobian().interior_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const Index3 refined {
                    2 * i + 1, 2 * j + 1, metric.dimension() == 3 ? 2 * k + 1 : 0};
                MetricFieldBuilderAccess::cell_coordinates(metric).x(i, j, k)
                    = coordinates.x(refined.i, refined.j, refined.k);
                MetricFieldBuilderAccess::cell_coordinates(metric).y(i, j, k)
                    = coordinates.y(refined.i, refined.j, refined.k);
                MetricFieldBuilderAccess::cell_coordinates(metric).z(i, j, k)
                    = coordinates.z(refined.i, refined.j, refined.k);
                MetricFieldBuilderAccess::jacobian(metric)(i, j, k)
                    = jacobian(refined.i, refined.j, refined.k);
            }
        }
    }
    const auto copy_faces = [&](const VectorFields& source,
                                FaceAreaVectors& target,
                                int normal_axis) {
        const auto extent = target.x.interior_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    Index3 refined {2 * i + 1, 2 * j + 1, metric.dimension() == 3 ? 2 * k + 1 : 0};
                    refined[static_cast<std::size_t>(normal_axis)]
                        = 2 * Index3 {i, j, k}[static_cast<std::size_t>(normal_axis)];
                    target.x(i, j, k) = source.x(refined.i, refined.j, refined.k);
                    target.y(i, j, k) = source.y(refined.i, refined.j, refined.k);
                    target.z(i, j, k) = source.z(refined.i, refined.j, refined.k);
                }
            }
        }
    };
    copy_faces(s_i, MetricFieldBuilderAccess::i_faces(metric), 0);
    copy_faces(s_j, MetricFieldBuilderAccess::j_faces(metric), 1);
    if (metric.dimension() == 3) {
        copy_faces(s_k, MetricFieldBuilderAccess::k_faces(metric), 2);
    }
}

MetricField build_phenglei_metric(const StructuredBlock& block, const AlgorithmProfile& profile)
{
    const auto vertices = block.vertex_extent();
    for (int axis = 0; axis < block.cell_dimension(); ++axis) {
        if (vertices[static_cast<std::size_t>(axis)] < 3) {
            throw GeometryError(
                "phenglei_wcns metrics require at least three vertices per active direction");
        }
    }
    auto coordinates = refined_coordinates(block);
    const auto refined_extent = coordinates.x.interior_extent();
    VectorFields s_i(refined_extent);
    VectorFields s_j(refined_extent);
    VectorFields s_k(refined_extent);
    ScalarField jacobian(refined_extent);
    compute_refined_symmetric_metrics(coordinates, block.cell_dimension(), s_i, s_j, s_k, jacobian);
    auto metric = MetricFieldBuilderAccess::create(
        profile.kind(), block.cell_extent(), block.cell_dimension());
    copy_refined_to_metric(coordinates, s_i, s_j, s_k, jacobian, metric);
    return metric;
}

VectorFields scmm_cell_coordinates(const StructuredBlock& block)
{
    const auto vertices = block.vertex_extent();
    VectorFields vertex_fields(vertices);
    const auto copy = [&](const Array3D<Real>& source, ScalarField& target) {
        for (int k = 0; k < vertices.nk; ++k) {
            for (int j = 0; j < vertices.nj; ++j) {
                for (int i = 0; i < vertices.ni; ++i) {
                    target(i, j, k) = source(i, j, k);
                }
            }
        }
    };
    copy(block.coordinates.x, vertex_fields.x);
    copy(block.coordinates.y, vertex_fields.y);
    copy(block.coordinates.z, vertex_fields.z);

    const auto interpolate = [&](ScalarField field) {
        field = vertices_to_centers_axis(field, 0);
        field = vertices_to_centers_axis(field, 1);
        if (block.cell_dimension() == 3) {
            field = vertices_to_centers_axis(field, 2);
        }
        return field;
    };
    VectorFields centers(block.cell_extent());
    centers.x = interpolate(std::move(vertex_fields.x));
    centers.y = interpolate(std::move(vertex_fields.y));
    centers.z = interpolate(std::move(vertex_fields.z));
    return centers;
}

void compute_scmm_symmetric_metrics(const VectorFields& coordinates,
                                    int dimension,
                                    const AlgorithmProfile& profile,
                                    VectorFields& s_i,
                                    VectorFields& s_j,
                                    VectorFields& s_k,
                                    ScalarField& jacobian)
{
    const auto extent = coordinates.x.interior_extent();
    if (dimension == 2) {
        // 二维退化公式的面积矢量严格位于 x-y 平面，不依赖容器的默认初始化。
        s_i.z.fill(0.0);
        s_j.z.fill(0.0);
        s_i.x = scmm_center_derivative(coordinates.y, 1, profile);
        const auto dx_eta = scmm_center_derivative(coordinates.x, 1, profile);
        const auto dy_xi = scmm_center_derivative(coordinates.y, 0, profile);
        s_j.y = scmm_center_derivative(coordinates.x, 0, profile);
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                s_i.y(i, j, 0) = -dx_eta(i, j, 0);
                s_j.x(i, j, 0) = -dy_xi(i, j, 0);
            }
        }
        const auto di = scmm_center_derivative(dot_product(coordinates, s_i), 0, profile);
        const auto dj = scmm_center_derivative(dot_product(coordinates, s_j), 1, profile);
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                jacobian(i, j, 0) = 0.5 * (di(i, j, 0) + dj(i, j, 0));
            }
        }
        return;
    }

    s_i.x = scmm_symmetric_component(coordinates.z, coordinates.y, 1, 2, profile);
    s_i.y = scmm_symmetric_component(coordinates.x, coordinates.z, 1, 2, profile);
    s_i.z = scmm_symmetric_component(coordinates.y, coordinates.x, 1, 2, profile);
    s_j.x = scmm_symmetric_component(coordinates.z, coordinates.y, 2, 0, profile);
    s_j.y = scmm_symmetric_component(coordinates.x, coordinates.z, 2, 0, profile);
    s_j.z = scmm_symmetric_component(coordinates.y, coordinates.x, 2, 0, profile);
    s_k.x = scmm_symmetric_component(coordinates.z, coordinates.y, 0, 1, profile);
    s_k.y = scmm_symmetric_component(coordinates.x, coordinates.z, 0, 1, profile);
    s_k.z = scmm_symmetric_component(coordinates.y, coordinates.x, 0, 1, profile);
    const auto di = scmm_center_derivative(dot_product(coordinates, s_i), 0, profile);
    const auto dj = scmm_center_derivative(dot_product(coordinates, s_j), 1, profile);
    const auto dk = scmm_center_derivative(dot_product(coordinates, s_k), 2, profile);
    for (int k = 0; k < extent.nk; ++k) {
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                jacobian(i, j, k) = (di(i, j, k) + dj(i, j, k) + dk(i, j, k)) / 3.0;
            }
        }
    }
}

void copy_scmm_to_metric(const VectorFields& coordinates,
                         const VectorFields& s_i,
                         const VectorFields& s_j,
                         const VectorFields& s_k,
                         const ScalarField& jacobian,
                         const AlgorithmProfile& profile,
                         MetricField& metric)
{
    auto& target_coordinates = MetricFieldBuilderAccess::cell_coordinates(metric);
    auto& target_jacobian = MetricFieldBuilderAccess::jacobian(metric);
    const auto cells = jacobian.interior_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                target_coordinates.x(i, j, k) = coordinates.x(i, j, k);
                target_coordinates.y(i, j, k) = coordinates.y(i, j, k);
                target_coordinates.z(i, j, k) = coordinates.z(i, j, k);
                target_jacobian(i, j, k) = jacobian(i, j, k);
            }
        }
    }
    const auto interpolate_vector
        = [&](const VectorFields& source, FaceAreaVectors& target, int axis) {
              target.x = scmm_center_to_faces(source.x, axis, profile);
              target.y = scmm_center_to_faces(source.y, axis, profile);
              target.z = scmm_center_to_faces(source.z, axis, profile);
          };
    interpolate_vector(s_i, MetricFieldBuilderAccess::i_faces(metric), 0);
    interpolate_vector(s_j, MetricFieldBuilderAccess::j_faces(metric), 1);
    if (metric.dimension() == 3) {
        interpolate_vector(s_k, MetricFieldBuilderAccess::k_faces(metric), 2);
    }
}

MetricField build_scmm_metric(const StructuredBlock& block, const AlgorithmProfile& profile)
{
    const auto vertices = block.vertex_extent();
    for (int axis = 0; axis < block.cell_dimension(); ++axis) {
        if (vertices[static_cast<std::size_t>(axis)] < 6) {
            throw GeometryError(
                "scmm6_wcns metrics require at least six vertices per active direction");
        }
    }
    auto coordinates = scmm_cell_coordinates(block);
    const auto cells = block.cell_extent();
    VectorFields s_i(cells);
    VectorFields s_j(cells);
    VectorFields s_k(cells);
    ScalarField jacobian(cells);
    compute_scmm_symmetric_metrics(
        coordinates, block.cell_dimension(), profile, s_i, s_j, s_k, jacobian);
    auto metric = MetricFieldBuilderAccess::create(profile.kind(), cells, block.cell_dimension());
    copy_scmm_to_metric(coordinates, s_i, s_j, s_k, jacobian, profile, metric);
    return metric;
}

GeometryDiagnostics capture_reference_geometry(const StructuredBlock& block)
{
    const auto cells = block.cell_extent();
    GeometryDiagnostics diagnostics(cells);
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                diagnostics.reference_volume(i, j, k) = block.cell_metrics.volume(i, j, k);
            }
        }
    }
    const auto copy_area = [](const FaceMetric& source, ScalarField& target) {
        const auto extent = target.interior_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    target(i, j, k) = source.area(i, j, k);
                }
            }
        }
    };
    copy_area(block.face_metrics.i_faces, diagnostics.reference_i_area);
    copy_area(block.face_metrics.j_faces, diagnostics.reference_j_area);
    if (block.cell_dimension() == 3) {
        copy_area(block.face_metrics.k_faces, diagnostics.reference_k_area);
    }
    return diagnostics;
}

void validate_metric(MetricField& metric,
                     GeometryDiagnostics& diagnostics,
                     const MetricBuildOptions& options)
{
    const auto cells = metric.jacobian().interior_extent();
    for (int k = 0; k < cells.nk; ++k) {
        for (int j = 0; j < cells.nj; ++j) {
            for (int i = 0; i < cells.ni; ++i) {
                const Real reference = diagnostics.reference_volume(i, j, k);
                Real& value = MetricFieldBuilderAccess::jacobian(metric)(i, j, k);
                const bool finite_positive
                    = std::isfinite(value) && value > options.floors.jacobian_floor(reference);
                const Real difference = finite_positive
                    ? finite_relative_difference(value, reference)
                    : std::numeric_limits<Real>::infinity();
                diagnostics.maximum_jacobian_relative_difference
                    = std::max(diagnostics.maximum_jacobian_relative_difference, difference);
                if (!finite_positive
                    || difference > options.maximum_reference_relative_difference) {
                    if (metric.profile() == AlgorithmProfileKind::PhengleiWcns
                        && options.fallback == MetricFallback::PhengleiFiniteVolume) {
                        value = reference;
                        ++diagnostics.fallback_cell_count;
                    } else {
                        throw GeometryError(
                            "high-order Jacobian failed strict reference validation at cell ("
                            + std::to_string(i) + ',' + std::to_string(j) + ',' + std::to_string(k)
                            + "): value=" + std::to_string(value)
                            + ", reference=" + std::to_string(reference)
                            + ", relative_difference=" + std::to_string(difference));
                    }
                }
            }
        }
    }

    const auto validate_faces = [&](const FaceAreaVectors& faces, const ScalarField& reference) {
        const auto extent = reference.interior_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Real area = faces.area(i, j, k);
                    if (!std::isfinite(area)
                        || area <= options.floors.face_area_floor(reference(i, j, k))) {
                        throw GeometryError("high-order metric has a degenerate real face");
                    }
                }
            }
        }
    };
    validate_faces(metric.i_faces(), diagnostics.reference_i_area);
    validate_faces(metric.j_faces(), diagnostics.reference_j_area);
    if (metric.dimension() == 3) {
        validate_faces(metric.k_faces(), diagnostics.reference_k_area);
    }
}

} // namespace

void MetricBuildOptions::validate(const AlgorithmProfile& profile) const
{
    floors.validate();
    if (!std::isfinite(maximum_reference_relative_difference)
        || maximum_reference_relative_difference < 0.0) {
        throw std::invalid_argument(
            "maximum metric/reference relative difference must be finite and non-negative");
    }
    if (fallback == MetricFallback::PhengleiFiniteVolume
        && profile.kind() != AlgorithmProfileKind::PhengleiWcns) {
        throw ProfileError("finite-volume metric fallback is only valid for phenglei_wcns");
    }
}

Real FaceAreaVectors::area(int i, int j, int k) const
{
    return std::sqrt(x(i, j, k) * x(i, j, k) + y(i, j, k) * y(i, j, k) + z(i, j, k) * z(i, j, k));
}

MetricField::MetricField(AlgorithmProfileKind profile, Extent3 cell_extent, int dimension)
    : profile_(profile)
    , dimension_(dimension)
    , cell_coordinates_(cell_extent)
    , jacobian_(cell_extent)
    , i_faces_(face_extent(cell_extent, 0))
    , j_faces_(face_extent(cell_extent, 1))
    , k_faces_(face_extent(cell_extent, 2))
{
    if (dimension != 2 && dimension != 3) {
        throw std::invalid_argument("metric field dimension must be 2 or 3");
    }
}

const FaceAreaVectors& MetricField::k_faces() const
{
    if (dimension_ != 3) {
        throw std::logic_error("a 2D metric field has no K-face family");
    }
    return k_faces_;
}

MetricInitializationResult initialize_metric_field(StructuredBlock& block,
                                                   const AlgorithmProfile& profile,
                                                   const MetricBuildOptions& options)
{
    options.validate(profile);
    compute_metrics(block);
    auto diagnostics = capture_reference_geometry(block);
    MetricField metric = [&] {
        if (profile.kind() == AlgorithmProfileKind::PhengleiWcns) {
            return build_phenglei_metric(block, profile);
        }
        return build_scmm_metric(block, profile);
    }();
    validate_metric(metric, diagnostics, options);
    return {std::move(metric), std::move(diagnostics)};
}

MetricField extract_metric_field(const MetricField& source, Index3 cell_begin, Extent3 cell_extent)
{
    if (!cell_extent.valid() || cell_extent.size() == 0) {
        throw GeometryError("metric extraction extent must be nonempty");
    }
    const auto source_cells = source.jacobian().interior_extent();
    for (int axis = 0; axis < source.dimension(); ++axis) {
        if (cell_begin[static_cast<std::size_t>(axis)] < 0
            || cell_begin[static_cast<std::size_t>(axis)]
                    + cell_extent[static_cast<std::size_t>(axis)]
                > source_cells[static_cast<std::size_t>(axis)]) {
            throw GeometryError("metric extraction range is outside source field");
        }
    }
    if (source.dimension() == 2 && (cell_begin.k != 0 || cell_extent.nk != 1)) {
        throw GeometryError("two-dimensional metric extraction has an invalid k range");
    }
    MetricField result
        = MetricFieldBuilderAccess::create(source.profile(), cell_extent, source.dimension());
    const auto copy_scalar = [](const Array3D<Real>& input, Array3D<Real>& output, Index3 begin) {
        const auto extent = output.interior_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    output(i, j, k) = input(begin.i + i, begin.j + j, begin.k + k);
                }
            }
        }
    };
    copy_scalar(source.cell_coordinates().x,
                MetricFieldBuilderAccess::cell_coordinates(result).x,
                cell_begin);
    copy_scalar(source.cell_coordinates().y,
                MetricFieldBuilderAccess::cell_coordinates(result).y,
                cell_begin);
    copy_scalar(source.cell_coordinates().z,
                MetricFieldBuilderAccess::cell_coordinates(result).z,
                cell_begin);
    copy_scalar(source.jacobian(), MetricFieldBuilderAccess::jacobian(result), cell_begin);
    const auto copy_faces = [&](const FaceAreaVectors& input, FaceAreaVectors& output) {
        copy_scalar(input.x, output.x, cell_begin);
        copy_scalar(input.y, output.y, cell_begin);
        copy_scalar(input.z, output.z, cell_begin);
    };
    copy_faces(source.i_faces(), MetricFieldBuilderAccess::i_faces(result));
    copy_faces(source.j_faces(), MetricFieldBuilderAccess::j_faces(result));
    if (source.dimension() == 3) {
        copy_faces(source.k_faces(), MetricFieldBuilderAccess::k_faces(result));
    }
    return result;
}

std::vector<Real> pack_metric_field(const MetricField& metric)
{
    std::vector<Real> result;
    const auto append = [&](const Array3D<Real>& values) {
        const auto extent = values.interior_extent();
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    result.push_back(values(i, j, k));
                }
            }
        }
    };
    append(metric.cell_coordinates().x);
    append(metric.cell_coordinates().y);
    append(metric.cell_coordinates().z);
    append(metric.jacobian());
    const auto append_faces = [&](const FaceAreaVectors& faces) {
        append(faces.x);
        append(faces.y);
        append(faces.z);
    };
    append_faces(metric.i_faces());
    append_faces(metric.j_faces());
    if (metric.dimension() == 3) append_faces(metric.k_faces());
    return result;
}

std::size_t metric_field_payload_size(Extent3 cell_extent, int dimension)
{
    if ((dimension != 2 && dimension != 3) || !cell_extent.valid() || cell_extent.size() == 0
        || (dimension == 2 && cell_extent.nk != 1)) {
        throw GeometryError("metric payload extent or dimension is invalid");
    }
    const auto face_count = [](Extent3 cells, int axis) {
        ++cells[static_cast<std::size_t>(axis)];
        return cells.size();
    };
    std::size_t result
        = 4 * cell_extent.size() + 3 * face_count(cell_extent, 0) + 3 * face_count(cell_extent, 1);
    if (dimension == 3) result += 3 * face_count(cell_extent, 2);
    return result;
}

MetricField unpack_metric_field(AlgorithmProfileKind profile,
                                Extent3 cell_extent,
                                int dimension,
                                const std::vector<Real>& payload)
{
    MetricField result = MetricFieldBuilderAccess::create(profile, cell_extent, dimension);
    std::size_t offset = 0;
    const auto extract = [&](Array3D<Real>& values) {
        const auto extent = values.interior_extent();
        if (offset + extent.size() > payload.size()) {
            throw GeometryError("packed metric payload is truncated");
        }
        for (int k = 0; k < extent.nk; ++k) {
            for (int j = 0; j < extent.nj; ++j) {
                for (int i = 0; i < extent.ni; ++i) {
                    const Real value = payload[offset++];
                    if (!std::isfinite(value)) {
                        throw GeometryError("packed metric payload is non-finite");
                    }
                    values(i, j, k) = value;
                }
            }
        }
    };
    extract(MetricFieldBuilderAccess::cell_coordinates(result).x);
    extract(MetricFieldBuilderAccess::cell_coordinates(result).y);
    extract(MetricFieldBuilderAccess::cell_coordinates(result).z);
    extract(MetricFieldBuilderAccess::jacobian(result));
    const auto extract_faces = [&](FaceAreaVectors& faces) {
        extract(faces.x);
        extract(faces.y);
        extract(faces.z);
    };
    extract_faces(MetricFieldBuilderAccess::i_faces(result));
    extract_faces(MetricFieldBuilderAccess::j_faces(result));
    if (dimension == 3) extract_faces(MetricFieldBuilderAccess::k_faces(result));
    if (offset != payload.size()) {
        throw GeometryError("packed metric payload has trailing values");
    }
    return result;
}

} // namespace wcns
