#include <cgnslib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
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

class CgnsFile {
public:
    explicit CgnsFile(const char* path)
    {
        check_cgns(cg_open(path, CG_MODE_READ, &id_), "cg_open RANS field");
    }
    ~CgnsFile()
    {
        if (id_ >= 0) cg_close(id_);
    }
    int id() const noexcept { return id_; }

private:
    int id_ = -1;
};

double parse_real(const char* text)
{
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != std::string(text).size() || !std::isfinite(value)) {
        throw std::invalid_argument("target x must be finite");
    }
    return value;
}

std::size_t index(int ni, int i, int j)
{
    return static_cast<std::size_t>(i) + static_cast<std::size_t>(ni) * j;
}

std::vector<double> read_field(int file,
                               int base,
                               int zone,
                               int solution,
                               const char* name,
                               int ni,
                               int nj)
{
    std::vector<double> result(static_cast<std::size_t>(ni) * nj);
    const std::array<cgsize_t, 2> lower {1, 1};
    const std::array<cgsize_t, 2> upper {ni, nj};
    check_cgns(cg_field_read(file,
                             base,
                             zone,
                             solution,
                             name,
                             RealDouble,
                             lower.data(),
                             upper.data(),
                             result.data()),
               name);
    return result;
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 4) {
        std::cerr << "usage: wcns_extract_rans_profile <field.cgns> <target-x> <output.txt>\n";
        return EXIT_FAILURE;
    }
    try {
        const double target_x = parse_real(argv[2]);
        CgnsFile input(argv[1]);
        int bases = 0;
        check_cgns(cg_nbases(input.id(), &bases), "cg_nbases RANS field");
        if (bases != 1) throw std::runtime_error("RANS profile requires exactly one CGNS base");
        int cell_dimension = 0;
        int physical_dimension = 0;
        char base_name[33] {};
        check_cgns(cg_base_read(
                       input.id(), 1, base_name, &cell_dimension, &physical_dimension),
                   "cg_base_read RANS field");
        if (cell_dimension != 2 || physical_dimension != 2) {
            throw std::runtime_error("RANS profile currently requires a two-dimensional field");
        }
        int zones = 0;
        check_cgns(cg_nzones(input.id(), 1, &zones), "cg_nzones RANS field");
        if (zones != 1) throw std::runtime_error("RANS profile requires exactly one CGNS zone");
        std::array<cgsize_t, 6> size {};
        char zone_name[33] {};
        check_cgns(cg_zone_read(input.id(), 1, 1, zone_name, size.data()),
                   "cg_zone_read RANS field");
        const int vertex_ni = static_cast<int>(size[0]);
        const int vertex_nj = static_cast<int>(size[1]);
        const int cell_ni = static_cast<int>(size[2]);
        const int cell_nj = static_cast<int>(size[3]);
        if (cell_ni != vertex_ni - 1 || cell_nj != vertex_nj - 1) {
            throw std::runtime_error("RANS profile requires an ordinary structured zone");
        }

        const auto vertex_count = static_cast<std::size_t>(vertex_ni) * vertex_nj;
        std::vector<double> vertex_x(vertex_count);
        std::vector<double> vertex_y(vertex_count);
        const std::array<cgsize_t, 2> lower {1, 1};
        const std::array<cgsize_t, 2> upper {vertex_ni, vertex_nj};
        check_cgns(cg_coord_read(input.id(),
                                 1,
                                 1,
                                 "CoordinateX",
                                 RealDouble,
                                 lower.data(),
                                 upper.data(),
                                 vertex_x.data()),
                   "cg_coord_read CoordinateX");
        check_cgns(cg_coord_read(input.id(),
                                 1,
                                 1,
                                 "CoordinateY",
                                 RealDouble,
                                 lower.data(),
                                 upper.data(),
                                 vertex_y.data()),
                   "cg_coord_read CoordinateY");
        std::vector<double> cell_x(static_cast<std::size_t>(cell_ni) * cell_nj);
        std::vector<double> cell_y(cell_x.size());
        for (int j = 0; j < cell_nj; ++j) {
            for (int i = 0; i < cell_ni; ++i) {
                const auto cell = index(cell_ni, i, j);
                const std::array<std::size_t, 4> vertices {
                    index(vertex_ni, i, j),
                    index(vertex_ni, i + 1, j),
                    index(vertex_ni, i, j + 1),
                    index(vertex_ni, i + 1, j + 1),
                };
                for (const auto vertex : vertices) {
                    cell_x[cell] += 0.25 * vertex_x[vertex];
                    cell_y[cell] += 0.25 * vertex_y[vertex];
                }
            }
        }
        int selected_i = 0;
        double selected_error = std::numeric_limits<double>::infinity();
        for (int i = 0; i < cell_ni; ++i) {
            double mean_x = 0.0;
            for (int j = 0; j < cell_nj; ++j) mean_x += cell_x[index(cell_ni, i, j)];
            mean_x /= cell_nj;
            const double error = std::abs(mean_x - target_x);
            if (error < selected_error) {
                selected_error = error;
                selected_i = i;
            }
        }

        int solutions = 0;
        check_cgns(cg_nsols(input.id(), 1, 1, &solutions), "cg_nsols RANS field");
        if (solutions != 1) throw std::runtime_error("RANS profile requires one flow solution");
        GridLocation_t location = GridLocationNull;
        char solution_name[33] {};
        check_cgns(cg_sol_info(input.id(), 1, 1, 1, solution_name, &location),
                   "cg_sol_info RANS field");
        if (location != CellCenter) {
            throw std::runtime_error("RANS profile requires a cell-centered solution");
        }
        const auto density = read_field(input.id(), 1, 1, 1, "Density", cell_ni, cell_nj);
        const auto momentum = read_field(input.id(), 1, 1, 1, "MomentumX", cell_ni, cell_nj);
        const auto k = read_field(input.id(), 1, 1, 1, "k", cell_ni, cell_nj);
        const auto omega = read_field(input.id(), 1, 1, 1, "omega", cell_ni, cell_nj);
        const auto viscosity_ratio
            = read_field(input.id(), 1, 1, 1, "mu_t_over_mu", cell_ni, cell_nj);
        const auto y_plus = read_field(input.id(), 1, 1, 1, "wall_y_plus", cell_ni, cell_nj);

        std::ofstream output(argv[3]);
        if (!output) throw std::runtime_error("cannot open RANS profile output");
        const double lower_wall_y
            = 0.5 * (vertex_y[index(vertex_ni, selected_i, 0)]
                     + vertex_y[index(vertex_ni, selected_i + 1, 0)]);
        const double upper_wall_y
            = 0.5 * (vertex_y[index(vertex_ni, selected_i, vertex_nj - 1)]
                     + vertex_y[index(vertex_ni, selected_i + 1, vertex_nj - 1)]);
        output << std::setprecision(17);
        output << "# target_x " << target_x << "\n# selected_i " << selected_i
               << "\n# x y wall_distance_lower wall_distance_upper rho u k omega "
                  "mu_t_over_mu wall_y_plus\n";
        for (int j = 0; j < cell_nj; ++j) {
            const auto cell = index(cell_ni, selected_i, j);
            if (!std::isfinite(density[cell]) || density[cell] <= 0.0) {
                throw std::runtime_error("RANS profile contains invalid density");
            }
            output << cell_x[cell] << ' ' << cell_y[cell] << ' '
                   << std::abs(cell_y[cell] - lower_wall_y) << ' '
                   << std::abs(upper_wall_y - cell_y[cell]) << ' ' << density[cell] << ' '
                   << momentum[cell] / density[cell] << ' ' << k[cell] << ' ' << omega[cell]
                   << ' ' << viscosity_ratio[cell] << ' ' << y_plus[cell] << '\n';
        }
        if (!output) throw std::runtime_error("failed while writing RANS profile");
        std::cout << "extracted i=" << selected_i << " target_x=" << target_x
                  << " error=" << selected_error << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
