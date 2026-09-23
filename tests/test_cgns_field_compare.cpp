#include "test_support.hpp"

#include <cgnslib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check_cgns(int status, const char* operation)
{
    if (status != CG_OK) {
        throw std::runtime_error(std::string(operation) + ": " + cg_get_error());
    }
}

std::vector<double> read_field(int file,
                               int zone,
                               int dimension,
                               const char* name,
                               const std::array<cgsize_t, 9>& size)
{
    std::size_t count = 1;
    std::array<cgsize_t, 3> lower {{1, 1, 1}};
    std::array<cgsize_t, 3> upper {{1, 1, 1}};
    for (int axis = 0; axis < dimension; ++axis) {
        upper[static_cast<std::size_t>(axis)]
            = size[static_cast<std::size_t>(dimension + axis)];
        count *= static_cast<std::size_t>(upper[static_cast<std::size_t>(axis)]);
    }
    std::vector<double> values(count);
    check_cgns(cg_field_read(file,
                             1,
                             zone,
                             1,
                             name,
                             RealDouble,
                             lower.data(),
                             upper.data(),
                             values.data()),
               "cg_field_read comparison");
    return values;
}

} // namespace

// 独立重读并比较连续计算与检查点续算的守恒场及可选模型场。
int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "usage: wcns_cgns_field_compare <expected.cgns> "
                     "<actual.cgns> <absolute-tolerance> [additional-field ...]\n";
        return EXIT_FAILURE;
    }
    int expected_file = 0;
    int actual_file = 0;
    try {
        const double tolerance = std::stod(argv[3]);
        check_cgns(cg_open(argv[1], CG_MODE_READ, &expected_file), "cg_open expected field");
        check_cgns(cg_open(argv[2], CG_MODE_READ, &actual_file), "cg_open actual field");
        int expected_zones = 0;
        int actual_zones = 0;
        int expected_dimension = 0;
        int expected_physical_dimension = 0;
        int actual_dimension = 0;
        int actual_physical_dimension = 0;
        char expected_base[33] = {};
        char actual_base[33] = {};
        check_cgns(cg_base_read(expected_file,
                                1,
                                expected_base,
                                &expected_dimension,
                                &expected_physical_dimension),
                   "cg_base_read expected");
        check_cgns(cg_base_read(actual_file,
                                1,
                                actual_base,
                                &actual_dimension,
                                &actual_physical_dimension),
                   "cg_base_read actual");
        WCNS_REQUIRE(expected_dimension == actual_dimension);
        WCNS_REQUIRE(expected_physical_dimension == actual_physical_dimension);
        check_cgns(cg_nzones(expected_file, 1, &expected_zones), "cg_nzones expected");
        check_cgns(cg_nzones(actual_file, 1, &actual_zones), "cg_nzones actual");
        WCNS_REQUIRE(actual_zones == expected_zones);
        std::vector<std::string> fields {
            "Density",
            "MomentumX",
            "MomentumY",
            "MomentumZ",
            "EnergyStagnationDensity",
        };
        for (int argument = 4; argument < argc; ++argument) {
            fields.emplace_back(argv[argument]);
        }
        double maximum_difference = 0.0;
        std::string maximum_field;
        int maximum_zone = 0;
        std::size_t maximum_index = 0;
        double maximum_expected = 0.0;
        double maximum_actual = 0.0;
        for (int zone = 1; zone <= expected_zones; ++zone) {
            std::array<cgsize_t, 9> expected_size {{}};
            std::array<cgsize_t, 9> actual_size {{}};
            char expected_zone[33] = {};
            char actual_zone[33] = {};
            check_cgns(cg_zone_read(
                           expected_file, 1, zone, expected_zone, expected_size.data()),
                       "cg_zone_read expected");
            check_cgns(cg_zone_read(actual_file, 1, zone, actual_zone, actual_size.data()),
                       "cg_zone_read actual");
            WCNS_REQUIRE(std::string(expected_zone) == actual_zone);
            for (int axis = 0; axis < expected_dimension; ++axis) {
                WCNS_REQUIRE(expected_size[static_cast<std::size_t>(expected_dimension + axis)]
                             == actual_size[static_cast<std::size_t>(actual_dimension + axis)]);
            }
            for (const auto& field : fields) {
                const auto expected = read_field(
                    expected_file, zone, expected_dimension, field.c_str(), expected_size);
                const auto actual = read_field(
                    actual_file, zone, actual_dimension, field.c_str(), actual_size);
                WCNS_REQUIRE(expected.size() == actual.size());
                for (std::size_t index = 0; index < expected.size(); ++index) {
                    WCNS_REQUIRE(std::isfinite(expected[index]));
                    WCNS_REQUIRE(std::isfinite(actual[index]));
                    const double difference = std::abs(expected[index] - actual[index]);
                    if (difference > maximum_difference) {
                        maximum_difference = difference;
                        maximum_field = field;
                        maximum_zone = zone;
                        maximum_index = index;
                        maximum_expected = expected[index];
                        maximum_actual = actual[index];
                    }
                }
            }
        }
        std::cout << "field comparison max_abs=" << maximum_difference
                  << ", tolerance=" << tolerance << ", field=" << maximum_field
                  << ", zone=" << maximum_zone << ", index=" << maximum_index
                  << ", expected=" << maximum_expected << ", actual=" << maximum_actual
                  << '\n';
        WCNS_REQUIRE(maximum_difference <= tolerance);
        check_cgns(cg_close(expected_file), "cg_close expected field");
        expected_file = 0;
        check_cgns(cg_close(actual_file), "cg_close actual field");
        actual_file = 0;
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        if (expected_file > 0) cg_close(expected_file);
        if (actual_file > 0) cg_close(actual_file);
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
