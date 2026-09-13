#include <cgnslib.h>

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
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

template <class T>
T read_value(std::ifstream& input, const char* label)
{
    T value {};
    input.read(reinterpret_cast<char*>(&value), sizeof(value));
    if (!input) throw std::runtime_error(std::string("truncated ") + label);
    return value;
}

void expect_record(std::ifstream& input, std::int32_t expected, const char* label)
{
    const auto actual = read_value<std::int32_t>(input, label);
    if (actual != expected) {
        throw std::runtime_error(std::string("unexpected Fortran record size for ") + label);
    }
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
               "cg_boco_write flat-plate grid");
    check_cgns(cg_boco_gridlocation_write(file, base, zone, boundary, Vertex),
               "cg_boco_gridlocation_write flat-plate grid");
}

} // namespace

int main(int argc, char** argv)
{
    if (argc != 3) {
        std::cerr << "usage: wcns_convert_tmr_flatplate_p3d_to_cgns "
                     "<grid.p3d> <output.cgns>\n";
        return EXIT_FAILURE;
    }
    try {
        std::ifstream input(argv[1], std::ios::binary);
        if (!input) throw std::runtime_error("cannot open TMR flat-plate PLOT3D input");

        expect_record(input, sizeof(std::int32_t), "block count");
        const auto blocks = read_value<std::int32_t>(input, "block count");
        expect_record(input, sizeof(std::int32_t), "block count trailer");
        expect_record(input, 3 * sizeof(std::int32_t), "dimensions");
        const auto ni = read_value<std::int32_t>(input, "ni");
        const auto nj = read_value<std::int32_t>(input, "nj");
        const auto nk = read_value<std::int32_t>(input, "nk");
        expect_record(input, 3 * sizeof(std::int32_t), "dimensions trailer");
        if (blocks != 1 || ni != 2 || nj < 9 || nk < 3 || (nj - 1) % 4 != 0) {
            throw std::runtime_error("unsupported TMR flat-plate dimensions");
        }

        const std::size_t count = static_cast<std::size_t>(ni)
            * static_cast<std::size_t>(nj) * static_cast<std::size_t>(nk);
        if (count > static_cast<std::size_t>(std::numeric_limits<std::int32_t>::max())
                / (3 * sizeof(double))) {
            throw std::runtime_error("TMR flat-plate coordinate record is too large");
        }
        const auto coordinate_bytes = static_cast<std::int32_t>(3 * count * sizeof(double));
        expect_record(input, coordinate_bytes, "coordinates");
        std::vector<double> x3(count);
        std::vector<double> y3(count);
        std::vector<double> z3(count);
        input.read(reinterpret_cast<char*>(x3.data()),
                   static_cast<std::streamsize>(x3.size() * sizeof(double)));
        input.read(reinterpret_cast<char*>(y3.data()),
                   static_cast<std::streamsize>(y3.size() * sizeof(double)));
        input.read(reinterpret_cast<char*>(z3.data()),
                   static_cast<std::streamsize>(z3.size() * sizeof(double)));
        if (!input) throw std::runtime_error("truncated TMR flat-plate coordinates");
        expect_record(input, coordinate_bytes, "coordinates trailer");
        if (input.peek() != std::char_traits<char>::eof()) {
            throw std::runtime_error("unexpected trailing TMR flat-plate data");
        }

        const auto source = [=](int i, int j, int k) {
            return static_cast<std::size_t>(i)
                + static_cast<std::size_t>(ni)
                    * (static_cast<std::size_t>(j)
                       + static_cast<std::size_t>(nj) * static_cast<std::size_t>(k));
        };
        const std::size_t count2 = static_cast<std::size_t>(nj) * static_cast<std::size_t>(nk);
        std::vector<double> x(count2);
        std::vector<double> y(count2);
        for (int k = 0; k < nk; ++k) {
            for (int j = 0; j < nj; ++j) {
                const auto target = static_cast<std::size_t>(j)
                    + static_cast<std::size_t>(nj) * static_cast<std::size_t>(k);
                x[target] = x3[source(0, j, k)];
                y[target] = z3[source(0, j, k)];
                if (!std::isfinite(x[target]) || !std::isfinite(y[target])
                    || x3[source(0, j, k)] != x3[source(1, j, k)]
                    || z3[source(0, j, k)] != z3[source(1, j, k)]) {
                    throw std::runtime_error("invalid or non-degenerate flat-plate coordinates");
                }
            }
        }

        CgnsFile output(argv[2]);
        int base = 0;
        int zone = 0;
        int coordinate = 0;
        check_cgns(cg_base_write(output.id(), "Base", 2, 2, &base),
                   "cg_base_write flat-plate grid");
        const cgsize_t size[6] = {nj, nk, nj - 1, nk - 1, 0, 0};
        check_cgns(cg_zone_write(output.id(), base, "Zone", size, Structured, &zone),
                   "cg_zone_write flat-plate grid");
        check_cgns(cg_coord_write(output.id(),
                                  base,
                                  zone,
                                  RealDouble,
                                  "CoordinateX",
                                  x.data(),
                                  &coordinate),
                   "cg_coord_write flat-plate CoordinateX");
        check_cgns(cg_coord_write(output.id(),
                                  base,
                                  zone,
                                  RealDouble,
                                  "CoordinateY",
                                  y.data(),
                                  &coordinate),
                   "cg_coord_write flat-plate CoordinateY");

        const cgsize_t plate_first = (nj - 1) / 4 + 1;
        const cgsize_t plate_last = 3 * (nj - 1) / 4 + 1;
        write_boundary(output.id(), base, zone, "Inflow", BCInflowSubsonic, {1, 1, 1, nk});
        write_boundary(output.id(), base, zone, "Outflow", BCOutflowSubsonic, {nj, 1, nj, nk});
        write_boundary(output.id(),
                       base,
                       zone,
                       "BottomUpstream",
                       BCSymmetryPlane,
                       {1, 1, plate_first, 1});
        write_boundary(output.id(),
                       base,
                       zone,
                       "Wall",
                       BCWallViscousHeatFlux,
                       {plate_first, 1, plate_last, 1});
        write_boundary(output.id(),
                       base,
                       zone,
                       "BottomDownstream",
                       BCSymmetryPlane,
                       {plate_last, 1, nj, 1});
        write_boundary(output.id(), base, zone, "Top", BCSymmetryPlane, {1, nk, nj, nk});
        output.close();
        std::cout << "converted TMR flat-plate grid " << nj << 'x' << nk << " wall=["
                  << plate_first << ',' << plate_last << "]\n";
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
