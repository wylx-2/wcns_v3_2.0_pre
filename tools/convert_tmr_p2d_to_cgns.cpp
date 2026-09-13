#include <cgnslib.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
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

class CgnsFile {
public:
    explicit CgnsFile(const std::string& path)
    {
        check_cgns(cg_open(path.c_str(), CG_MODE_WRITE, &id_), "cg_open converted grid");
    }
    ~CgnsFile()
    {
        if (id_ >= 0) cg_close(id_);
    }
    int id() const noexcept { return id_; }
    void close()
    {
        check_cgns(cg_close(id_), "cg_close converted grid");
        id_ = -1;
    }

private:
    int id_ = -1;
};

int parse_index(const char* text, const char* label)
{
    std::size_t consumed = 0;
    const int value = std::stoi(text, &consumed);
    if (consumed != std::string(text).size()) {
        throw std::invalid_argument(std::string("invalid ") + label);
    }
    return value;
}

void write_boundary(int file,
                    int base,
                    int zone,
                    const char* name,
                    BCType_t type,
                    const std::vector<cgsize_t>& range)
{
    int boundary = 0;
    check_cgns(cg_boco_write(
                   file, base, zone, name, type, PointRange, 2, range.data(), &boundary),
               "cg_boco_write converted grid");
    check_cgns(cg_boco_gridlocation_write(file, base, zone, boundary, Vertex),
               "cg_boco_gridlocation_write converted grid");
}

void write_connection(int file,
                      int base,
                      int zone,
                      const char* name,
                      const std::vector<cgsize_t>& receiver,
                      const std::vector<cgsize_t>& donor)
{
    const std::vector<int> transform {-1, 2};
    int connection = 0;
    check_cgns(cg_1to1_write(file,
                             base,
                             zone,
                             name,
                             "Zone",
                             receiver.data(),
                             donor.data(),
                             transform.data(),
                             &connection),
               "cg_1to1_write converted grid");
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 5) {
        std::cerr << "usage: wcns_convert_tmr_p2d_to_cgns <input.p2dfmt> <output.cgns> "
                     "<wall-first-i> <wall-last-i>\n";
        return EXIT_FAILURE;
    }
    try {
        const int wall_first = parse_index(argv[3], "wall-first-i");
        const int wall_last = parse_index(argv[4], "wall-last-i");
        std::ifstream input(argv[1]);
        if (!input) throw std::runtime_error("cannot open TMR PLOT3D input");
        int blocks = 0;
        int ni = 0;
        int nj = 0;
        if (!(input >> blocks >> ni >> nj) || blocks != 1 || ni < 3 || nj < 3
            || wall_first <= 1 || wall_first >= wall_last || wall_last >= ni) {
            throw std::runtime_error("unsupported TMR PLOT3D header or wall range");
        }
        const auto count = static_cast<std::size_t>(ni) * static_cast<std::size_t>(nj);
        std::vector<double> x(count);
        std::vector<double> y(count);
        for (double& value : x) {
            if (!(input >> value) || !std::isfinite(value)) {
                throw std::runtime_error("invalid CoordinateX in TMR PLOT3D input");
            }
        }
        for (double& value : y) {
            if (!(input >> value) || !std::isfinite(value)) {
                throw std::runtime_error("invalid CoordinateY in TMR PLOT3D input");
            }
        }
        std::string trailing;
        if (input >> trailing) {
            throw std::runtime_error("unexpected trailing TMR PLOT3D data");
        }

        CgnsFile output(argv[2]);
        int base = 0;
        int zone = 0;
        int coordinate = 0;
        check_cgns(cg_base_write(output.id(), "Base", 2, 2, &base),
                   "cg_base_write converted grid");
        const cgsize_t size[6] = {ni, nj, ni - 1, nj - 1, 0, 0};
        check_cgns(cg_zone_write(output.id(), base, "Zone", size, Structured, &zone),
                   "cg_zone_write converted grid");
        check_cgns(cg_coord_write(output.id(),
                                  base,
                                  zone,
                                  RealDouble,
                                  "CoordinateX",
                                  x.data(),
                                  &coordinate),
                   "cg_coord_write converted CoordinateX");
        check_cgns(cg_coord_write(output.id(),
                                  base,
                                  zone,
                                  RealDouble,
                                  "CoordinateY",
                                  y.data(),
                                  &coordinate),
                   "cg_coord_write converted CoordinateY");

        write_boundary(output.id(), base, zone, "FarfieldILower", BCFarfield, {1, 1, 1, nj});
        write_boundary(output.id(), base, zone, "FarfieldIUpper", BCFarfield, {ni, 1, ni, nj});
        write_boundary(output.id(), base, zone, "FarfieldJUpper", BCFarfield, {1, nj, ni, nj});
        write_boundary(output.id(),
                       base,
                       zone,
                       "Wall",
                       BCWallViscousHeatFlux,
                       {wall_first, 1, wall_last, 1});
        write_connection(output.id(),
                         base,
                         zone,
                         "WakeLowerToUpper",
                         {1, 1, wall_first, 1},
                         {ni, 1, wall_last, 1});
        write_connection(output.id(),
                         base,
                         zone,
                         "WakeUpperToLower",
                         {ni, 1, wall_last, 1},
                         {1, 1, wall_first, 1});
        output.close();
        std::cout << "converted TMR PLOT3D grid " << ni << 'x' << nj << " wall=["
                  << wall_first << ',' << wall_last << "]\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
