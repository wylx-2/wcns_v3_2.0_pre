#include "test_support.hpp"

#include <wcns/mesh/metrics.hpp>

// 验收有限体积参考中心、面积、体积和退化网格拒绝路径。
void test_metrics()
{
    using wcns::StructuredBlock;

    StructuredBlock trapezoid(0, "trapezoid", 0, 2, 2, {2, 2, 1}, 1);
    trapezoid.coordinates.x(0, 0, 0) = 0.0;
    trapezoid.coordinates.y(0, 0, 0) = 0.0;
    trapezoid.coordinates.x(1, 0, 0) = 2.0;
    trapezoid.coordinates.y(1, 0, 0) = 0.0;
    trapezoid.coordinates.x(1, 1, 0) = 1.5;
    trapezoid.coordinates.y(1, 1, 0) = 1.0;
    trapezoid.coordinates.x(0, 1, 0) = 0.0;
    trapezoid.coordinates.y(0, 1, 0) = 1.0;
    trapezoid.coordinates.z.fill(0.0);

    wcns::compute_metrics(trapezoid);
    WCNS_REQUIRE_NEAR(trapezoid.cell_metrics.volume(0, 0, 0), 1.75, 1.0e-14);
    WCNS_REQUIRE_NEAR(trapezoid.cell_metrics.center_x(0, 0, 0), 0.875, 1.0e-14);
    WCNS_REQUIRE_NEAR(trapezoid.cell_metrics.center_y(0, 0, 0), 0.5, 1.0e-14);
    WCNS_REQUIRE_NEAR(trapezoid.face_metrics.i_faces.normal_x(0, 0, 0), 1.0, 1.0e-14);
    WCNS_REQUIRE_NEAR(trapezoid.face_metrics.j_faces.normal_y(0, 0, 0), 1.0, 1.0e-14);

    StructuredBlock degenerate(1, "degenerate", 0, 2, 2, {2, 2, 1}, 1);
    degenerate.coordinates.x.fill(0.0);
    degenerate.coordinates.y.fill(0.0);
    degenerate.coordinates.z.fill(0.0);
    WCNS_REQUIRE_THROWS(wcns::GeometryError, wcns::compute_metrics(degenerate));

    StructuredBlock inverted(2, "inverted", 0, 2, 2, {2, 2, 1}, 1);
    inverted.coordinates.x(0, 0, 0) = 1.0;
    inverted.coordinates.y(0, 0, 0) = 0.0;
    inverted.coordinates.x(1, 0, 0) = 0.0;
    inverted.coordinates.y(1, 0, 0) = 0.0;
    inverted.coordinates.x(1, 1, 0) = 0.0;
    inverted.coordinates.y(1, 1, 0) = 1.0;
    inverted.coordinates.x(0, 1, 0) = 1.0;
    inverted.coordinates.y(0, 1, 0) = 1.0;
    inverted.coordinates.z.fill(0.0);
    WCNS_REQUIRE_THROWS(wcns::GeometryError, wcns::compute_metrics(inverted));

    StructuredBlock embedded(3, "embedded", 0, 2, 3, {2, 2, 1}, 1);
    embedded.coordinates.x(0, 0, 0) = 0.0;
    embedded.coordinates.y(0, 0, 0) = 0.0;
    embedded.coordinates.x(1, 0, 0) = 1.0;
    embedded.coordinates.y(1, 0, 0) = 0.0;
    embedded.coordinates.x(1, 1, 0) = 1.0;
    embedded.coordinates.y(1, 1, 0) = 1.0;
    embedded.coordinates.x(0, 1, 0) = 0.0;
    embedded.coordinates.y(0, 1, 0) = 1.0;
    embedded.coordinates.z.fill(4.0);
    wcns::compute_metrics(embedded);
    WCNS_REQUIRE_NEAR(embedded.cell_metrics.volume(0, 0, 0), 1.0, 1.0e-15);
    embedded.coordinates.z(1, 1, 0) = 4.1;
    WCNS_REQUIRE_THROWS(wcns::GeometryError, wcns::compute_metrics(embedded));
}
