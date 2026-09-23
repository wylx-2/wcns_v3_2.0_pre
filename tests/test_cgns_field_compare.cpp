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

std::vector<double> read_field(int file, const char* name, const std::array<cgsize_t, 9>& size)
{
    const auto count = static_cast<std::size_t>(size[2] * size[3]);
    std::vector<double> values(count);
    cgsize_t lower[2] = {1, 1};
    cgsize_t upper[2] = {size[2], size[3]};
    check_cgns(cg_field_read(file, 1, 1, 1, name, RealDouble, lower, upper, values.data()),
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
        check_cgns(cg_nzones(expected_file, 1, &expected_zones), "cg_nzones expected");
        check_cgns(cg_nzones(actual_file, 1, &actual_zones), "cg_nzones actual");
        WCNS_REQUIRE(expected_zones == 1);
        WCNS_REQUIRE(actual_zones == expected_zones);
        std::array<cgsize_t, 9> expected_size {{}};
        std::array<cgsize_t, 9> actual_size {{}};
        char expected_zone[33] = {};
        char actual_zone[33] = {};
        check_cgns(cg_zone_read(expected_file, 1, 1, expected_zone, expected_size.data()),
                   "cg_zone_read expected");
        check_cgns(cg_zone_read(actual_file, 1, 1, actual_zone, actual_size.data()),
                   "cg_zone_read actual");
        WCNS_REQUIRE(std::string(expected_zone) == actual_zone);
        WCNS_REQUIRE(expected_size[2] == actual_size[2]);
        WCNS_REQUIRE(expected_size[3] == actual_size[3]);
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
        for (const auto& field : fields) {
            const auto expected = read_field(expected_file, field.c_str(), expected_size);
            const auto actual = read_field(actual_file, field.c_str(), actual_size);
            WCNS_REQUIRE(expected.size() == actual.size());
            for (std::size_t index = 0; index < expected.size(); ++index) {
                WCNS_REQUIRE(std::isfinite(expected[index]));
                WCNS_REQUIRE(std::isfinite(actual[index]));
                maximum_difference
                    = std::max(maximum_difference, std::abs(expected[index] - actual[index]));
            }
        }
        std::cout << "field comparison max_abs=" << maximum_difference
                  << ", tolerance=" << tolerance << '\n';
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
