#include <cgnslib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

void check(int status, const char* operation)
{
    if (status != CG_OK) {
        throw std::runtime_error(std::string(operation) + ": " + cg_get_error());
    }
}

std::size_t index3(cgsize_t ni, cgsize_t nj, int i, int j, int k)
{
    return (static_cast<std::size_t>(k) * static_cast<std::size_t>(nj)
        + static_cast<std::size_t>(j)) * static_cast<std::size_t>(ni)
        + static_cast<std::size_t>(i);
}

std::vector<double> read_field(
    int file, int zone, const char* name, const std::array<cgsize_t, 9>& size)
{
    std::array<cgsize_t, 3> lower {{1, 1, 1}};
    std::array<cgsize_t, 3> upper {{size[3], size[4], size[5]}};
    std::vector<double> values(static_cast<std::size_t>(
        size[3] * size[4] * size[5]));
    check(cg_field_read(
        file, 1, zone, 1, name, RealDouble,
        lower.data(), upper.data(), values.data()),
        "read field");
    return values;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: extract_channel_profile <field.cgns> <profile.txt>\n";
        return EXIT_FAILURE;
    }
    int file = 0;
    try {
        check(cg_open(argv[1], CG_MODE_READ, &file), "open field");
        int zones = 0;
        check(cg_nzones(file, 1, &zones), "count zones");
        if (zones <= 0) throw std::runtime_error("field has no zones");

        std::vector<double> y;
        std::vector<double> sum_u, sum_v, sum_w, sum_rho, sum_temperature;
        std::vector<double> sum_mach, sum_uu, sum_vv, sum_ww, sum_uv;
        std::vector<std::size_t> counts;
        for (int zone = 1; zone <= zones; ++zone) {
            char zone_name[33] = {};
            std::array<cgsize_t, 9> size {{}};
            check(cg_zone_read(file, 1, zone, zone_name, size.data()), "read zone");
            if (size[4] <= 0 || size[3] <= 0 || size[5] <= 0) {
                throw std::runtime_error("field zone is not three-dimensional");
            }
            if (y.empty()) {
                const auto nj = static_cast<std::size_t>(size[4]);
                y.assign(nj, 0.0);
                sum_u.assign(nj, 0.0); sum_v.assign(nj, 0.0);
                sum_w.assign(nj, 0.0); sum_rho.assign(nj, 0.0);
                sum_temperature.assign(nj, 0.0); sum_mach.assign(nj, 0.0);
                sum_uu.assign(nj, 0.0); sum_vv.assign(nj, 0.0);
                sum_ww.assign(nj, 0.0); sum_uv.assign(nj, 0.0);
                counts.assign(nj, 0);

                std::array<cgsize_t, 3> lower {{1, 1, 1}};
                std::array<cgsize_t, 3> upper {{size[0], size[1], size[2]}};
                std::vector<double> coordinate_y(static_cast<std::size_t>(
                    size[0] * size[1] * size[2]));
                check(cg_coord_read(
                    file, 1, zone, "CoordinateY", RealDouble,
                    lower.data(), upper.data(), coordinate_y.data()),
                    "read CoordinateY");
                for (int j = 0; j < size[4]; ++j) {
                    y[static_cast<std::size_t>(j)] = 0.5 * (
                        coordinate_y[index3(size[0], size[1], 0, j, 0)]
                        + coordinate_y[index3(size[0], size[1], 0, j + 1, 0)]);
                }
            } else if (static_cast<cgsize_t>(y.size()) != size[4]) {
                throw std::runtime_error("zones have different wall-normal extents");
            }

            const auto u = read_field(file, zone, "VelocityX", size);
            const auto v = read_field(file, zone, "VelocityY", size);
            const auto w = read_field(file, zone, "VelocityZ", size);
            const auto rho = read_field(file, zone, "Density", size);
            const auto temperature = read_field(file, zone, "Temperature", size);
            const auto mach = read_field(file, zone, "Mach", size);
            for (int k = 0; k < size[5]; ++k) {
                for (int j = 0; j < size[4]; ++j) {
                    const auto jj = static_cast<std::size_t>(j);
                    for (int i = 0; i < size[3]; ++i) {
                        const auto q = index3(size[3], size[4], i, j, k);
                        sum_u[jj] += u[q]; sum_v[jj] += v[q]; sum_w[jj] += w[q];
                        sum_rho[jj] += rho[q];
                        sum_temperature[jj] += temperature[q];
                        sum_mach[jj] += mach[q];
                        sum_uu[jj] += u[q] * u[q];
                        sum_vv[jj] += v[q] * v[q];
                        sum_ww[jj] += w[q] * w[q];
                        sum_uv[jj] += u[q] * v[q];
                        ++counts[jj];
                    }
                }
            }
        }
        check(cg_close(file), "close field");
        file = 0;

        std::ofstream output(argv[2]);
        if (!output) throw std::runtime_error("cannot open profile output");
        output << "# y mean_u mean_v mean_w rms_u rms_v rms_w reynolds_uv"
                  " mean_rho mean_temperature mean_mach\n";
        output << std::setprecision(17);
        for (std::size_t j = 0; j < y.size(); ++j) {
            const double count = static_cast<double>(counts[j]);
            const double mean_u = sum_u[j] / count;
            const double mean_v = sum_v[j] / count;
            const double mean_w = sum_w[j] / count;
            const double rms_u = std::sqrt(std::max(0.0, sum_uu[j] / count - mean_u * mean_u));
            const double rms_v = std::sqrt(std::max(0.0, sum_vv[j] / count - mean_v * mean_v));
            const double rms_w = std::sqrt(std::max(0.0, sum_ww[j] / count - mean_w * mean_w));
            output << y[j] << ' ' << mean_u << ' ' << mean_v << ' ' << mean_w
                   << ' ' << rms_u << ' ' << rms_v << ' ' << rms_w
                   << ' ' << (sum_uv[j] / count - mean_u * mean_v)
                   << ' ' << sum_rho[j] / count
                   << ' ' << sum_temperature[j] / count
                   << ' ' << sum_mach[j] / count << '\n';
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        if (file != 0) cg_close(file);
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
