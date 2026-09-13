#include <cgnslib.h>

#include <array>
#include <cmath>
#include <exception>
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

int parse_points(const char* text)
{
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != std::string(text).size()) {
        throw std::invalid_argument("points-j must be an integer");
    }
    return value;
}

double parse_stretching(const char* text)
{
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != std::string(text).size()) {
        throw std::invalid_argument("stretching must be a real number");
    }
    return value;
}

void boundary(int file,
              int base,
              int zone,
              const char* name,
              const std::array<cgsize_t, 4>& range)
{
    int index = 0;
    check(cg_boco_write(
              file, base, zone, name, BCWall, PointRange, 2, range.data(), &index),
          "cg_boco_write");
    check(cg_boco_gridlocation_write(file, base, zone, index, Vertex),
          "cg_boco_gridlocation_write");
}

int connection(int file,
               int base,
               int zone,
               const char* name,
               const char* donor,
               const std::array<cgsize_t, 4>& receiver,
               const std::array<cgsize_t, 4>& donor_range)
{
    const std::array<int, 2> transform {{1, 2}};
    int index = 0;
    check(cg_1to1_write(file,
                        base,
                        zone,
                        name,
                        donor,
                        receiver.data(),
                        donor_range.data(),
                        transform.data(),
                        &index),
          "cg_1to1_write");
    return index;
}

void periodic(int file, int base, int zone, int index, float translation)
{
    const std::array<float, 3> zero {{0.0F, 0.0F, 0.0F}};
    const std::array<float, 3> shift {{translation, 0.0F, 0.0F}};
    check(cg_1to1_periodic_write(
              file, base, zone, index, zero.data(), zero.data(), shift.data()),
          "cg_1to1_periodic_write");
}

void write_grid(const std::string& path, int nj, double stretching)
{
    constexpr int ni = 17;
    if (nj < 9 || nj > 4097 || !(stretching > 0.0) || stretching > 20.0
        || !std::isfinite(stretching)) {
        throw std::invalid_argument(
            "points-j must be in [9,4097] and stretching must be in (0,20]");
    }
    int file = -1;
    check(cg_open(path.c_str(), CG_MODE_WRITE, &file), "cg_open");
    try {
        int base = 0;
        int left = 0;
        int right = 0;
        int coordinate = 0;
        check(cg_base_write(file, "CouetteValidation", 2, 2, &base), "cg_base_write");
        cgsize_t size[6] = {ni, nj, ni - 1, nj - 1, 0, 0};
        check(cg_zone_write(file, base, "CouetteLeft", size, Structured, &left),
              "cg_zone_write left");
        check(cg_zone_write(file, base, "CouetteRight", size, Structured, &right),
              "cg_zone_write right");

        std::vector<double> x(static_cast<std::size_t>(ni * nj));
        std::vector<double> y(x.size());
        const auto coordinates = [&](int zone, double x0, double x1) {
            for (int j = 0; j < nj; ++j) {
                const double eta = static_cast<double>(j) / static_cast<double>(nj - 1);
                const double ordinate = 0.5
                    * (1.0 + std::tanh(stretching * (eta - 0.5))
                                   / std::tanh(0.5 * stretching));
                for (int i = 0; i < ni; ++i) {
                    const auto index = static_cast<std::size_t>(j * ni + i);
                    const double fraction
                        = static_cast<double>(i) / static_cast<double>(ni - 1);
                    x[index] = x0 + (x1 - x0) * fraction;
                    y[index] = ordinate;
                }
            }
            check(cg_coord_write(
                      file, base, zone, RealDouble, "CoordinateX", x.data(), &coordinate),
                  "cg_coord_write X");
            check(cg_coord_write(
                      file, base, zone, RealDouble, "CoordinateY", y.data(), &coordinate),
                  "cg_coord_write Y");
        };
        coordinates(left, 0.0, 0.5);
        coordinates(right, 0.5, 1.0);

        boundary(file, base, left, "left-jmin", {1, 1, ni, 1});
        boundary(file, base, left, "left-jmax", {1, nj, ni, nj});
        boundary(file, base, right, "right-jmin", {1, 1, ni, 1});
        boundary(file, base, right, "right-jmax", {1, nj, ni, nj});
        connection(file,
                   base,
                   left,
                   "left-to-right",
                   "CouetteRight",
                   {ni, 1, ni, nj},
                   {1, 1, 1, nj});
        connection(file,
                   base,
                   right,
                   "right-to-left",
                   "CouetteLeft",
                   {1, 1, 1, nj},
                   {ni, 1, ni, nj});
        const int forward = connection(file,
                                       base,
                                       left,
                                       "periodic-left-to-right",
                                       "CouetteRight",
                                       {1, 1, 1, nj},
                                       {ni, 1, ni, nj});
        const int reverse = connection(file,
                                       base,
                                       right,
                                       "periodic-right-to-left",
                                       "CouetteLeft",
                                       {ni, 1, ni, nj},
                                       {1, 1, 1, nj});
        periodic(file, base, left, forward, 1.0F);
        periodic(file, base, right, reverse, -1.0F);
        check(cg_close(file), "cg_close");
    } catch (...) {
        cg_close(file);
        throw;
    }
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2 || argc > 4) {
        std::cerr << "usage: wcns_generate_couette_validation_cgns <output.cgns> "
                     "[points-j=65] [stretching=6]\n";
        return 2;
    }
    try {
        const int nj = argc >= 3 ? parse_points(argv[2]) : 65;
        const double stretching = argc >= 4 ? parse_stretching(argv[3]) : 6.0;
        write_grid(argv[1], nj, stretching);
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "wcns_generate_couette_validation_cgns: " << error.what() << '\n';
        return 1;
    }
}
