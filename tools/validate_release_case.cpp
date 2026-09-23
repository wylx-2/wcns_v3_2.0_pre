#include <cgnslib.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
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

double parse_real(const char* text, const char* name)
{
    std::size_t consumed = 0;
    const double value = std::stod(text, &consumed);
    if (consumed != std::string(text).size() || !std::isfinite(value)) {
        throw std::invalid_argument(std::string(name) + " must be finite");
    }
    return value;
}

struct ZoneFields {
    int dimension = 0;
    std::vector<cgsize_t> cells;
    std::vector<std::array<double, 3>> cell_centers;
    std::map<std::string, std::vector<double>> fields;
};

using FieldFile = std::map<std::string, ZoneFields>;

FieldFile read_fields(const std::string& path)
{
    int file = 0;
    check_cgns(cg_open(path.c_str(), CG_MODE_READ, &file), "cg_open release field");
    try {
        int bases = 0;
        check_cgns(cg_nbases(file, &bases), "cg_nbases release field");
        if (bases != 1) throw std::runtime_error("release field requires one CGNS base");
        char base_name[33] = {};
        int dimension = 0;
        int physical_dimension = 0;
        check_cgns(cg_base_read(file, 1, base_name, &dimension, &physical_dimension),
                   "cg_base_read release field");
        if (dimension != 2 && dimension != 3) {
            throw std::runtime_error("release field dimension must be two or three");
        }
        int zone_count = 0;
        check_cgns(cg_nzones(file, 1, &zone_count), "cg_nzones release field");
        FieldFile result;
        for (int zone = 1; zone <= zone_count; ++zone) {
            char zone_name[33] = {};
            cgsize_t size[9] = {};
            check_cgns(cg_zone_read(file, 1, zone, zone_name, size), "cg_zone_read release field");
            ZoneFields values;
            values.dimension = dimension;
            std::size_t count = 1;
            std::vector<cgsize_t> vertices;
            for (int axis = 0; axis < dimension; ++axis) {
                vertices.push_back(size[axis]);
                values.cells.push_back(size[dimension + axis]);
                count *= static_cast<std::size_t>(size[dimension + axis]);
            }
            std::size_t vertex_count = 1;
            for (const auto extent : vertices) {
                vertex_count *= static_cast<std::size_t>(extent);
            }
            std::array<std::vector<double>, 3> coordinates;
            for (auto& coordinate : coordinates)
                coordinate.assign(vertex_count, 0.0);
            const std::array<const char*, 3> coordinate_names {{
                "CoordinateX",
                "CoordinateY",
                "CoordinateZ",
            }};
            std::vector<cgsize_t> vertex_lower(static_cast<std::size_t>(dimension), 1);
            for (int axis = 0; axis < dimension; ++axis) {
                check_cgns(cg_coord_read(file,
                                         1,
                                         zone,
                                         coordinate_names[static_cast<std::size_t>(axis)],
                                         RealDouble,
                                         vertex_lower.data(),
                                         vertices.data(),
                                         coordinates[static_cast<std::size_t>(axis)].data()),
                           "cg_coord_read release field");
            }
            values.cell_centers.resize(count);
            const auto vertex_index = [&](int i, int j, int k) {
                const std::size_t ni = static_cast<std::size_t>(vertices[0]);
                const std::size_t nj = static_cast<std::size_t>(vertices[1]);
                return (static_cast<std::size_t>(k) * nj + static_cast<std::size_t>(j)) * ni
                    + static_cast<std::size_t>(i);
            };
            const int cell_nk = dimension == 3 ? static_cast<int>(values.cells[2]) : 1;
            std::size_t cell = 0;
            for (int k = 0; k < cell_nk; ++k) {
                for (int j = 0; j < static_cast<int>(values.cells[1]); ++j) {
                    for (int i = 0; i < static_cast<int>(values.cells[0]); ++i) {
                        std::array<double, 3> center {{0.0, 0.0, 0.0}};
                        const int corner_count = dimension == 3 ? 8 : 4;
                        for (int corner = 0; corner < corner_count; ++corner) {
                            const int di = corner & 1;
                            const int dj = (corner >> 1) & 1;
                            const int dk = dimension == 3 ? ((corner >> 2) & 1) : 0;
                            const auto vertex = vertex_index(i + di, j + dj, k + dk);
                            for (int axis = 0; axis < dimension; ++axis) {
                                center[static_cast<std::size_t>(axis)]
                                    += coordinates[static_cast<std::size_t>(axis)][vertex]
                                    / static_cast<double>(corner_count);
                            }
                        }
                        values.cell_centers[cell++] = center;
                    }
                }
            }
            int solutions = 0;
            check_cgns(cg_nsols(file, 1, zone, &solutions), "cg_nsols release field");
            if (solutions != 1) {
                throw std::runtime_error("release field zone requires one solution");
            }
            int field_count = 0;
            check_cgns(cg_nfields(file, 1, zone, 1, &field_count), "cg_nfields release field");
            std::vector<cgsize_t> lower(static_cast<std::size_t>(dimension), 1);
            std::vector<cgsize_t> upper = values.cells;
            for (int field = 1; field <= field_count; ++field) {
                DataType_t data_type = DataTypeNull;
                char field_name[33] = {};
                check_cgns(cg_field_info(file, 1, zone, 1, field, &data_type, field_name),
                           "cg_field_info release field");
                std::vector<double> field_values(count);
                check_cgns(cg_field_read(file,
                                         1,
                                         zone,
                                         1,
                                         field_name,
                                         RealDouble,
                                         lower.data(),
                                         upper.data(),
                                         field_values.data()),
                           "cg_field_read release field");
                values.fields.emplace(field_name, std::move(field_values));
            }
            if (!result.emplace(zone_name, std::move(values)).second) {
                throw std::runtime_error("duplicate zone name in release field");
            }
        }
        check_cgns(cg_close(file), "cg_close release field");
        return result;
    } catch (...) {
        cg_close(file);
        throw;
    }
}

struct TecplotZone {
    int dimension = 0;
    std::vector<std::string> variables;
    std::vector<std::vector<double>> rows;
};

using TecplotFile = std::map<std::string, TecplotZone>;

std::vector<std::string> quoted_values(const std::string& line)
{
    std::vector<std::string> result;
    std::size_t begin = 0;
    while ((begin = line.find('"', begin)) != std::string::npos) {
        const auto end = line.find('"', begin + 1);
        if (end == std::string::npos) {
            throw std::runtime_error("unterminated quoted Tecplot value");
        }
        result.push_back(line.substr(begin + 1, end - begin - 1));
        begin = end + 1;
    }
    return result;
}

std::size_t tecplot_extent(const std::string& line, const std::string& label, bool required)
{
    const auto position = line.find(label);
    if (position == std::string::npos) {
        if (!required) return 1;
        throw std::runtime_error("Tecplot zone is missing " + label);
    }
    const auto begin = position + label.size();
    auto end = begin;
    while (end < line.size() && line[end] >= '0' && line[end] <= '9')
        ++end;
    if (end == begin) throw std::runtime_error("invalid Tecplot zone extent");
    return static_cast<std::size_t>(std::stoull(line.substr(begin, end - begin)));
}

TecplotFile read_tecplot(const std::string& path)
{
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open Tecplot field " + path);
    std::string line;
    std::vector<std::string> variables;
    TecplotFile result;
    while (std::getline(input, line)) {
        if (line.rfind("VARIABLES=", 0) == 0) {
            variables = quoted_values(line);
            continue;
        }
        if (line.rfind("ZONE ", 0) != 0) continue;
        if (variables.size() < 3) {
            throw std::runtime_error("Tecplot field has no valid variable list");
        }
        const auto names = quoted_values(line);
        if (names.empty()) throw std::runtime_error("Tecplot zone has no name");
        TecplotZone zone;
        zone.dimension = variables.size() >= 3 && variables[2] == "Z" ? 3 : 2;
        zone.variables = variables;
        const auto count = tecplot_extent(line, "I=", true) * tecplot_extent(line, "J=", true)
            * tecplot_extent(line, "K=", false);
        zone.rows.reserve(count);
        if (line.find("DATAPACKING=BLOCK") != std::string::npos) {
            std::vector<double> packed(count * variables.size());
            for (auto& value : packed) {
                if (!(input >> value)) {
                    throw std::runtime_error("Tecplot block data are truncated");
                }
            }
            std::getline(input, line);
            for (std::size_t row = 0; row < count; ++row) {
                std::vector<double> record(variables.size());
                for (std::size_t column = 0; column < variables.size(); ++column) {
                    record[column] = packed[column * count + row];
                }
                zone.rows.push_back(std::move(record));
            }
        } else {
            for (std::size_t row = 0; row < count; ++row) {
                if (!std::getline(input, line)) {
                    throw std::runtime_error("Tecplot zone data are truncated");
                }
                std::istringstream values(line);
                std::vector<double> record;
                double value = 0.0;
                while (values >> value)
                    record.push_back(value);
                if (record.size() != variables.size()) {
                    throw std::runtime_error("Tecplot row has the wrong column count");
                }
                zone.rows.push_back(std::move(record));
            }
        }
        if (!result.emplace(names.front(), std::move(zone)).second) {
            throw std::runtime_error("duplicate Tecplot zone name");
        }
    }
    if (variables.empty() || result.empty()) {
        throw std::runtime_error("Tecplot field is empty");
    }
    return result;
}

std::string tecplot_cgns_field(const std::string& name)
{
    static const std::map<std::string, std::string> names {
        {"rho", "Density"},
        {"u", "VelocityX"},
        {"v", "VelocityY"},
        {"w", "VelocityZ"},
        {"p", "Pressure"},
        {"T", "Temperature"},
        {"rho_u", "MomentumX"},
        {"rho_v", "MomentumY"},
        {"rho_w", "MomentumZ"},
        {"rho_E", "EnergyStagnationDensity"},
        {"mach", "Mach"},
        {"jacobian", "Jacobian"},
    };
    const auto iterator = names.find(name);
    return iterator == names.end() ? name : iterator->second;
}

const std::vector<double>& require_field(const ZoneFields& zone, const std::string& name)
{
    const auto iterator = zone.fields.find(name);
    if (iterator == zone.fields.end()) {
        throw std::runtime_error("release output is missing field " + name);
    }
    return iterator->second;
}

void require_tolerance(double tolerance)
{
    if (!(tolerance >= 0.0) || !std::isfinite(tolerance)) {
        throw std::invalid_argument("tolerance must be finite and nonnegative");
    }
}

void validate_tecplot_consistency(const std::string& cgns_path,
                                  const std::string& tecplot_path,
                                  double tolerance)
{
    require_tolerance(tolerance);
    const auto cgns = read_fields(cgns_path);
    const auto tecplot = read_tecplot(tecplot_path);
    if (cgns.size() != tecplot.size()) {
        throw std::runtime_error("CGNS/Tecplot zone counts differ");
    }
    double maximum = 0.0;
    std::size_t samples = 0;
    for (const auto& [name, cgns_zone] : cgns) {
        const auto zone_iterator = tecplot.find(name);
        if (zone_iterator == tecplot.end()) {
            throw std::runtime_error("Tecplot field is missing zone " + name);
        }
        const auto& zone = zone_iterator->second;
        if (zone.dimension != cgns_zone.dimension
            || zone.rows.size() != cgns_zone.cell_centers.size()) {
            throw std::runtime_error("CGNS/Tecplot zone dimensions differ");
        }
        const auto coordinate_columns = static_cast<std::size_t>(zone.dimension);
        if (zone.variables.size() - coordinate_columns != cgns_zone.fields.size()) {
            throw std::runtime_error("CGNS/Tecplot field counts differ");
        }
        for (std::size_t cell = 0; cell < zone.rows.size(); ++cell) {
            for (int axis = 0; axis < zone.dimension; ++axis) {
                maximum = std::max(
                    maximum,
                    std::abs(zone.rows[cell][static_cast<std::size_t>(axis)]
                             - cgns_zone.cell_centers[cell][static_cast<std::size_t>(axis)]));
                ++samples;
            }
            for (std::size_t column = coordinate_columns; column < zone.variables.size();
                 ++column) {
                const auto field_name = tecplot_cgns_field(zone.variables[column]);
                const auto& field = require_field(cgns_zone, field_name);
                maximum = std::max(maximum, std::abs(zone.rows[cell][column] - field[cell]));
                ++samples;
            }
        }
    }
    if (maximum > tolerance) {
        throw std::runtime_error("CGNS/Tecplot values differ beyond tolerance");
    }
    std::cout << std::setprecision(17) << "check=tecplot_consistency samples=" << samples
              << " max_abs=" << maximum << " tolerance=" << tolerance << '\n';
}

void validate_derived_fields(const std::string& path,
                             double gamma,
                             double expected_viscosity,
                             double expected_jacobian,
                             double tolerance)
{
    if (!(gamma > 1.0) || !(expected_viscosity > 0.0) || !(expected_jacobian > 0.0)) {
        throw std::invalid_argument("derived-field references must be positive");
    }
    require_tolerance(tolerance);
    const auto file = read_fields(path);
    double maximum = 0.0;
    std::size_t samples = 0;
    const auto check = [&](double actual, double exact) {
        if (!std::isfinite(actual) || !std::isfinite(exact)) {
            throw std::runtime_error("derived-field comparison is non-finite");
        }
        maximum = std::max(maximum, std::abs(actual - exact));
        ++samples;
    };
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        const auto& rho = require_field(zone, "Density");
        const auto& u = require_field(zone, "VelocityX");
        const auto& v = require_field(zone, "VelocityY");
        const auto& w = require_field(zone, "VelocityZ");
        const auto& pressure = require_field(zone, "Pressure");
        const auto& temperature = require_field(zone, "Temperature");
        const auto& rho_u = require_field(zone, "MomentumX");
        const auto& rho_v = require_field(zone, "MomentumY");
        const auto& rho_w = require_field(zone, "MomentumZ");
        const auto& energy = require_field(zone, "EnergyStagnationDensity");
        const auto& sound_speed = require_field(zone, "sound_speed");
        const auto& mach = require_field(zone, "Mach");
        const auto& enthalpy = require_field(zone, "total_enthalpy");
        const auto& entropy = require_field(zone, "entropy_proxy");
        const auto& viscosity = require_field(zone, "viscosity");
        const auto& jacobian = require_field(zone, "Jacobian");
        for (std::size_t cell = 0; cell < rho.size(); ++cell) {
            const double speed_squared = u[cell] * u[cell] + v[cell] * v[cell] + w[cell] * w[cell];
            const double expected_energy
                = pressure[cell] / (gamma - 1.0) + 0.5 * rho[cell] * speed_squared;
            const double expected_sound = std::sqrt(gamma * pressure[cell] / rho[cell]);
            check(rho_u[cell], rho[cell] * u[cell]);
            check(rho_v[cell], rho[cell] * v[cell]);
            check(rho_w[cell], rho[cell] * w[cell]);
            check(energy[cell], expected_energy);
            check(temperature[cell], pressure[cell] / rho[cell]);
            check(sound_speed[cell], expected_sound);
            check(mach[cell], std::sqrt(speed_squared) / expected_sound);
            check(enthalpy[cell], (energy[cell] + pressure[cell]) / rho[cell]);
            check(entropy[cell], pressure[cell] / std::pow(rho[cell], gamma));
            check(viscosity[cell], expected_viscosity);
            check(jacobian[cell], expected_jacobian);
        }
    }
    if (maximum > tolerance) {
        throw std::runtime_error("derived field exceeds tolerance");
    }
    std::cout << std::setprecision(17) << "check=derived_fields samples=" << samples
              << " max_abs=" << maximum << " tolerance=" << tolerance << '\n';
}

void validate_nonzero_field(const std::string& path,
                            const std::string& field_name,
                            double threshold)
{
    require_tolerance(threshold);
    const auto file = read_fields(path);
    double maximum = 0.0;
    std::size_t samples = 0;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        for (const double value : require_field(zone, field_name)) {
            if (!std::isfinite(value)) {
                throw std::runtime_error("nonzero-field check encountered non-finite data");
            }
            maximum = std::max(maximum, std::abs(value));
            ++samples;
        }
    }
    if (samples == 0 || maximum < threshold) {
        throw std::runtime_error("field does not reach the required nonzero magnitude");
    }
    std::cout << std::setprecision(17) << "check=nonzero field=" << field_name
              << " samples=" << samples << " max_abs=" << maximum << " threshold=" << threshold
              << '\n';
}

void report_field_extrema(const std::string& path, const std::string& field_name)
{
    const auto file = read_fields(path);
    double minimum = std::numeric_limits<double>::infinity();
    double maximum = -std::numeric_limits<double>::infinity();
    double maximum_absolute = 0.0;
    std::size_t samples = 0;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        for (const double value : require_field(zone, field_name)) {
            if (!std::isfinite(value)) {
                throw std::runtime_error("field-extrema check encountered non-finite data");
            }
            minimum = std::min(minimum, value);
            maximum = std::max(maximum, value);
            maximum_absolute = std::max(maximum_absolute, std::abs(value));
            ++samples;
        }
    }
    if (samples == 0) throw std::runtime_error("field-extrema check found no samples");
    std::cout << std::setprecision(17) << "check=field_extrema field=" << field_name
              << " samples=" << samples << " min=" << minimum << " max=" << maximum
              << " max_abs=" << maximum_absolute << '\n';
}

void report_flatplate_profile(const std::string& path, double target_x)
{
    const auto file = read_fields(path);
    for (const auto& [zone_name, zone] : file) {
        if (zone.dimension != 2 || zone.cells.size() != 2) {
            throw std::runtime_error("flatplate-profile requires a two-dimensional field");
        }
        const auto ni = static_cast<std::size_t>(zone.cells[0]);
        const auto nj = static_cast<std::size_t>(zone.cells[1]);
        const auto& density = require_field(zone, "Density");
        const auto& momentum_x = require_field(zone, "MomentumX");
        const auto& viscosity_ratio = require_field(zone, "mu_t_over_mu");
        const auto& wall_distance = require_field(zone, "wall_distance");
        const bool two_equation = zone.fields.find("k") != zone.fields.end()
            && zone.fields.find("omega") != zone.fields.end();
        const bool sa = zone.fields.find("nu_tilde") != zone.fields.end();
        if (!two_equation && !sa) {
            throw std::runtime_error(
                "flatplate-profile requires either k/omega or nu_tilde");
        }
        std::cout << std::setprecision(17)
                  << (two_equation
                          ? "# zone i j x y velocity_x k omega mu_t_over_mu wall_distance\n"
                          : "# zone i j x y velocity_x nu_tilde mu_t_over_mu wall_distance\n");
        for (std::size_t j = 0; j < nj; ++j) {
            std::size_t selected = j * ni;
            double distance = std::numeric_limits<double>::infinity();
            for (std::size_t i = 0; i < ni; ++i) {
                const auto cell = j * ni + i;
                const double candidate
                    = std::abs(zone.cell_centers[cell][0] - target_x);
                if (candidate < distance) {
                    selected = cell;
                    distance = candidate;
                }
            }
            if (!(density[selected] > 0.0) || !std::isfinite(momentum_x[selected])) {
                throw std::runtime_error("flatplate-profile encountered an invalid mean state");
            }
            std::cout << zone_name << ' ' << selected % ni << ' ' << j << ' '
                      << zone.cell_centers[selected][0] << ' '
                      << zone.cell_centers[selected][1] << ' '
                      << momentum_x[selected] / density[selected] << ' ';
            if (two_equation) {
                std::cout << require_field(zone, "k")[selected] << ' '
                          << require_field(zone, "omega")[selected] << ' ';
            } else {
                std::cout << require_field(zone, "nu_tilde")[selected] << ' ';
            }
            std::cout << viscosity_ratio[selected] << ' ' << wall_distance[selected] << '\n';
        }
    }
}

void report_flatplate_y_plus(const std::string& path,
                             double target_x,
                             double skin_friction,
                             double reynolds)
{
    if (!std::isfinite(skin_friction) || skin_friction < 0.0
        || !std::isfinite(reynolds) || reynolds <= 0.0) {
        throw std::runtime_error("flatplate-yplus requires Cf >= 0 and Re > 0");
    }
    const auto file = read_fields(path);
    for (const auto& [zone_name, zone] : file) {
        if (zone.dimension != 2 || zone.cells.size() != 2) {
            throw std::runtime_error("flatplate-yplus requires a two-dimensional field");
        }
        const auto ni = static_cast<std::size_t>(zone.cells[0]);
        const auto& wall_distance = require_field(zone, "wall_distance");
        std::size_t selected = 0;
        double streamwise_distance = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < ni; ++i) {
            const double candidate = std::abs(zone.cell_centers[i][0] - target_x);
            if (candidate < streamwise_distance) {
                selected = i;
                streamwise_distance = candidate;
            }
        }
        const double distance = wall_distance[selected];
        if (!std::isfinite(distance) || distance < 0.0) {
            throw std::runtime_error("flatplate-yplus encountered an invalid wall distance");
        }
        const double y_plus = reynolds * distance * std::sqrt(0.5 * skin_friction);
        std::cout << std::setprecision(17) << "check=flatplate_y_plus zone=" << zone_name
                  << " i=" << selected << " x=" << zone.cell_centers[selected][0]
                  << " wall_distance=" << distance << " Cf=" << skin_friction
                  << " Re=" << reynolds << " y_plus=" << y_plus << '\n';
    }
}

double validate_uniform(const std::string& path,
                        const std::map<std::string, double>& expected,
                        double tolerance)
{
    require_tolerance(tolerance);
    const auto file = read_fields(path);
    double maximum = 0.0;
    std::size_t samples = 0;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        for (const auto& [field_name, reference] : expected) {
            for (const double value : require_field(zone, field_name)) {
                if (!std::isfinite(value)) {
                    throw std::runtime_error("uniform field contains a non-finite value");
                }
                maximum = std::max(maximum, std::abs(value - reference));
                ++samples;
            }
        }
    }
    if (maximum > tolerance) {
        throw std::runtime_error("uniform maximum error " + std::to_string(maximum)
                                 + " exceeds tolerance " + std::to_string(tolerance));
    }
    std::cout << std::setprecision(17) << "check=uniform samples=" << samples
              << " max_abs=" << maximum << " tolerance=" << tolerance << '\n';
    return maximum;
}

void validate_finite(const std::string& path)
{
    const auto file = read_fields(path);
    std::size_t samples = 0;
    double minimum_density = std::numeric_limits<double>::infinity();
    double minimum_pressure = std::numeric_limits<double>::infinity();
    double minimum_temperature = std::numeric_limits<double>::infinity();
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        for (const auto& [field_name, values] : zone.fields) {
            for (const double value : values) {
                if (!std::isfinite(value)) {
                    throw std::runtime_error("release field contains a non-finite value");
                }
                ++samples;
                if (field_name == "Density") minimum_density = std::min(minimum_density, value);
                if (field_name == "Pressure") minimum_pressure = std::min(minimum_pressure, value);
                if (field_name == "Temperature") {
                    minimum_temperature = std::min(minimum_temperature, value);
                }
            }
        }
    }
    for (const auto& value : {std::pair<const char*, double> {"density", minimum_density},
                              {"pressure", minimum_pressure},
                              {"temperature", minimum_temperature}}) {
        if (std::isfinite(value.second) && !(value.second > 0.0)) {
            throw std::runtime_error(std::string(value.first) + " is not positive");
        }
    }
    std::cout << std::setprecision(17) << "check=finite samples=" << samples
              << " min_rho=" << minimum_density << " min_p=" << minimum_pressure
              << " min_T=" << minimum_temperature << '\n';
}

struct FieldErrorStatistics {
    double weighted_absolute = 0.0;
    double weighted_squared = 0.0;
    double weighted_reference_squared = 0.0;
    double weight_sum = 0.0;
    double maximum = 0.0;
    std::size_t samples = 0;

    void add(double reference, double value, double weight)
    {
        if (!std::isfinite(reference) || !std::isfinite(value) || !std::isfinite(weight)
            || !(weight > 0.0)) {
            throw std::runtime_error("field-error comparison is non-finite");
        }
        const double error = std::abs(value - reference);
        weighted_absolute += weight * error;
        weighted_squared += weight * error * error;
        weighted_reference_squared += weight * reference * reference;
        weight_sum += weight;
        maximum = std::max(maximum, error);
        ++samples;
    }

    void print(const std::string& field_name) const
    {
        if (samples == 0 || !(weight_sum > 0.0)) {
            throw std::runtime_error("field-error comparison has no samples");
        }
        const double l2 = std::sqrt(weighted_squared / weight_sum);
        const double reference_l2 = std::sqrt(weighted_reference_squared / weight_sum);
        std::cout << "check=field_error field=" << field_name << " samples=" << samples
                  << " l1=" << weighted_absolute / weight_sum << " l2=" << l2 << " linf=" << maximum
                  << " relative_l2="
                  << l2 / std::max(reference_l2, std::numeric_limits<double>::min()) << '\n';
    }
};

void report_field_errors(const std::string& reference_path, const std::string& value_path)
{
    const auto reference = read_fields(reference_path);
    const auto value = read_fields(value_path);
    if (reference.size() != value.size() || reference.empty()) {
        throw std::runtime_error("field-error comparison zone counts differ");
    }
    std::map<std::string, FieldErrorStatistics> statistics;
    double maximum_coordinate_difference = 0.0;
    for (const auto& [zone_name, reference_zone] : reference) {
        const auto value_iterator = value.find(zone_name);
        if (value_iterator == value.end()) {
            throw std::runtime_error("field-error comparison zone names differ");
        }
        const auto& value_zone = value_iterator->second;
        if (reference_zone.dimension != value_zone.dimension
            || reference_zone.cells != value_zone.cells
            || reference_zone.cell_centers.size() != value_zone.cell_centers.size()
            || reference_zone.fields.size() != value_zone.fields.size()) {
            throw std::runtime_error("field-error comparison metadata differs");
        }
        const auto jacobian = reference_zone.fields.find("Jacobian");
        for (const auto& [field_name, reference_values] : reference_zone.fields) {
            const auto& values = require_field(value_zone, field_name);
            if (reference_values.size() != values.size()) {
                throw std::runtime_error("field-error comparison array sizes differ");
            }
            for (std::size_t cell = 0; cell < values.size(); ++cell) {
                const double weight
                    = jacobian == reference_zone.fields.end() ? 1.0 : jacobian->second[cell];
                statistics[field_name].add(reference_values[cell], values[cell], weight);
            }
        }
        for (std::size_t cell = 0; cell < reference_zone.cell_centers.size(); ++cell) {
            for (int axis = 0; axis < reference_zone.dimension; ++axis) {
                maximum_coordinate_difference = std::max(
                    maximum_coordinate_difference,
                    std::abs(reference_zone.cell_centers[cell][static_cast<std::size_t>(axis)]
                             - value_zone.cell_centers[cell][static_cast<std::size_t>(axis)]));
            }
        }
    }
    if (maximum_coordinate_difference > 1.0e-12) {
        throw std::runtime_error("field-error comparison coordinates differ");
    }
    std::cout << std::setprecision(17) << "check=field_error_header zones=" << reference.size()
              << " max_coordinate_difference=" << maximum_coordinate_difference << " weighting="
              << (statistics.find("Jacobian") == statistics.end() ? "uniform"
                                                                  : "reference_jacobian")
              << '\n';
    for (const auto& [field_name, field_statistics] : statistics) {
        field_statistics.print(field_name);
    }
}

void compare_fields(const std::string& lhs_path, const std::string& rhs_path, double tolerance)
{
    require_tolerance(tolerance);
    const auto lhs = read_fields(lhs_path);
    const auto rhs = read_fields(rhs_path);
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("release field comparison zone counts differ");
    }
    double maximum = 0.0;
    std::size_t samples = 0;
    for (const auto& [name, lhs_zone] : lhs) {
        const auto rhs_iterator = rhs.find(name);
        if (rhs_iterator == rhs.end()) {
            throw std::runtime_error("release field comparison zone names differ");
        }
        const auto& rhs_zone = rhs_iterator->second;
        if (lhs_zone.dimension != rhs_zone.dimension || lhs_zone.cells != rhs_zone.cells
            || lhs_zone.fields.size() != rhs_zone.fields.size()) {
            throw std::runtime_error("release field comparison metadata differs");
        }
        for (const auto& [field_name, lhs_values] : lhs_zone.fields) {
            const auto& rhs_values = require_field(rhs_zone, field_name);
            if (lhs_values.size() != rhs_values.size()) {
                throw std::runtime_error("release field comparison array sizes differ");
            }
            for (std::size_t index = 0; index < lhs_values.size(); ++index) {
                if (!std::isfinite(lhs_values[index]) || !std::isfinite(rhs_values[index])) {
                    throw std::runtime_error("release field comparison is non-finite");
                }
                maximum = std::max(maximum, std::abs(lhs_values[index] - rhs_values[index]));
                ++samples;
            }
        }
    }
    if (maximum > tolerance) {
        throw std::runtime_error("field comparison maximum difference " + std::to_string(maximum)
                                 + " exceeds tolerance " + std::to_string(tolerance));
    }
    std::cout << std::setprecision(17) << "check=compare samples=" << samples
              << " max_abs=" << maximum << " tolerance=" << tolerance << '\n';
}

void compare_scaled_primitive_fields(const std::string& lhs_path,
                                     const std::string& rhs_path,
                                     double gamma,
                                     double mach,
                                     double tolerance)
{
    require_tolerance(tolerance);
    if (!(gamma > 1.0) || !(mach > 0.0) || !std::isfinite(gamma)
        || !std::isfinite(mach)) {
        throw std::invalid_argument("scaled primitive comparison requires gamma>1 and Mach>0");
    }
    const auto lhs = read_fields(lhs_path);
    const auto rhs = read_fields(rhs_path);
    if (lhs.size() != rhs.size()) {
        throw std::runtime_error("scaled primitive comparison zone counts differ");
    }
    double maximum = 0.0;
    std::string maximum_field;
    std::size_t samples = 0;
    for (const auto& [name, lhs_zone] : lhs) {
        const auto rhs_iterator = rhs.find(name);
        if (rhs_iterator == rhs.end()) {
            throw std::runtime_error("scaled primitive comparison zone names differ");
        }
        const auto& rhs_zone = rhs_iterator->second;
        if (lhs_zone.dimension != rhs_zone.dimension || lhs_zone.cells != rhs_zone.cells
            || lhs_zone.fields.size() != rhs_zone.fields.size()) {
            throw std::runtime_error("scaled primitive comparison metadata differs");
        }
        for (const auto& [field_name, lhs_values] : lhs_zone.fields) {
            const auto& rhs_values = require_field(rhs_zone, field_name);
            if (lhs_values.size() != rhs_values.size()) {
                throw std::runtime_error("scaled primitive comparison array sizes differ");
            }
            const double scale = field_name == "Pressure" ? gamma * mach * mach : 1.0;
            for (std::size_t index = 0; index < lhs_values.size(); ++index) {
                if (!std::isfinite(lhs_values[index]) || !std::isfinite(rhs_values[index])) {
                    throw std::runtime_error("scaled primitive comparison is non-finite");
                }
                const double difference
                    = scale * std::abs(lhs_values[index] - rhs_values[index]);
                if (difference > maximum) {
                    maximum = difference;
                    maximum_field = field_name;
                }
                ++samples;
            }
        }
    }
    if (maximum > tolerance) {
        throw std::runtime_error("scaled primitive field " + maximum_field
                                 + " exceeds tolerance " + std::to_string(tolerance));
    }
    std::cout << std::setprecision(17) << "check=compare_scaled_primitive samples=" << samples
              << " max_abs=" << maximum << " maximum_field=" << maximum_field
              << " pressure_scale=" << gamma * mach * mach
              << " tolerance=" << tolerance << '\n';
}

struct SpatialValue {
    std::array<double, 3> coordinates;
    std::vector<double> fields;
};

std::vector<SpatialValue> spatial_values(const FieldFile& file,
                                         const std::vector<std::string>& field_names)
{
    std::vector<SpatialValue> result;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        for (std::size_t cell = 0; cell < zone.cell_centers.size(); ++cell) {
            SpatialValue value;
            value.coordinates = zone.cell_centers[cell];
            value.fields.reserve(field_names.size());
            for (const auto& name : field_names) {
                value.fields.push_back(require_field(zone, name).at(cell));
            }
            result.push_back(std::move(value));
        }
    }
    std::sort(result.begin(), result.end(), [](const auto& lhs, const auto& rhs) {
        return lhs.coordinates < rhs.coordinates;
    });
    return result;
}

void compare_spatial_fields(const std::string& lhs_path,
                            const std::string& rhs_path,
                            double field_tolerance,
                            double coordinate_tolerance)
{
    require_tolerance(field_tolerance);
    require_tolerance(coordinate_tolerance);
    const auto lhs = read_fields(lhs_path);
    const auto rhs = read_fields(rhs_path);
    if (lhs.empty() || rhs.empty()) {
        throw std::runtime_error("spatial comparison requires nonempty field files");
    }
    std::vector<std::string> fields;
    for (const auto& [name, values] : lhs.begin()->second.fields) {
        static_cast<void>(values);
        fields.push_back(name);
    }
    const auto check_field_set = [&](const FieldFile& file) {
        for (const auto& [zone_name, zone] : file) {
            static_cast<void>(zone_name);
            std::vector<std::string> names;
            for (const auto& [name, values] : zone.fields) {
                static_cast<void>(values);
                names.push_back(name);
            }
            if (names != fields) {
                throw std::runtime_error("spatial comparison field sets differ");
            }
        }
    };
    check_field_set(lhs);
    check_field_set(rhs);
    const auto lhs_values = spatial_values(lhs, fields);
    const auto rhs_values = spatial_values(rhs, fields);
    if (lhs_values.size() != rhs_values.size()) {
        throw std::runtime_error("spatial comparison cell counts differ");
    }
    double maximum_coordinate = 0.0;
    double maximum_field = 0.0;
    for (std::size_t cell = 0; cell < lhs_values.size(); ++cell) {
        for (std::size_t component = 0; component < 3; ++component) {
            maximum_coordinate = std::max(maximum_coordinate,
                                          std::abs(lhs_values[cell].coordinates[component]
                                                   - rhs_values[cell].coordinates[component]));
        }
        for (std::size_t field = 0; field < fields.size(); ++field) {
            if (!std::isfinite(lhs_values[cell].fields[field])
                || !std::isfinite(rhs_values[cell].fields[field])) {
                throw std::runtime_error("spatial comparison field is non-finite");
            }
            maximum_field = std::max(
                maximum_field,
                std::abs(lhs_values[cell].fields[field] - rhs_values[cell].fields[field]));
        }
    }
    if (maximum_coordinate > coordinate_tolerance || maximum_field > field_tolerance) {
        throw std::runtime_error("spatial field comparison exceeds tolerance");
    }
    std::cout << std::setprecision(17) << "check=compare_spatial cells=" << lhs_values.size()
              << " max_coordinate=" << maximum_coordinate
              << " coordinate_tolerance=" << coordinate_tolerance << " max_field=" << maximum_field
              << " field_tolerance=" << field_tolerance << '\n';
}

void validate_constant_series(const std::string& path, double tolerance)
{
    require_tolerance(tolerance);
    std::ifstream input(path);
    if (!input) throw std::runtime_error("cannot open release series: " + path);
    std::vector<std::vector<double>> rows;
    std::string line;
    while (std::getline(input, line)) {
        if (line.empty() || line.front() == '#') continue;
        std::istringstream values(line);
        std::vector<double> row;
        double value = 0.0;
        while (values >> value)
            row.push_back(value);
        if (!values.eof() || row.size() < 3) {
            throw std::runtime_error("release series contains an invalid data row");
        }
        rows.push_back(std::move(row));
    }
    if (rows.size() < 2) {
        throw std::runtime_error("release series requires at least two samples");
    }
    double maximum = 0.0;
    for (const auto& row : rows) {
        if (row.size() != rows.front().size()) {
            throw std::runtime_error("release series row widths differ");
        }
        for (std::size_t column = 2; column < row.size(); ++column) {
            if (!std::isfinite(row[column])) {
                throw std::runtime_error("release series contains a non-finite value");
            }
            maximum = std::max(maximum,
                               std::abs(row[column] - rows.front()[column])
                                   / std::max(1.0, std::abs(rows.front()[column])));
        }
    }
    if (maximum > tolerance) {
        throw std::runtime_error("series relative drift " + std::to_string(maximum)
                                 + " exceeds tolerance " + std::to_string(tolerance));
    }
    std::cout << std::setprecision(17) << "check=series_constant rows=" << rows.size()
              << " max_relative_drift=" << maximum << " tolerance=" << tolerance << '\n';
}

void validate_isentropic_vortex(const std::string& path,
                                double time,
                                double length,
                                double x0,
                                double y0,
                                double beta,
                                double background_u,
                                double background_v,
                                double gamma,
                                double mach,
                                double l1_tolerance)
{
    if (!(length > 0.0) || !(gamma > 1.0) || !(mach > 0.0)) {
        throw std::invalid_argument("vortex length, gamma and Mach are invalid");
    }
    require_tolerance(l1_tolerance);
    constexpr double pi = 3.141592653589793238462643383279502884;
    const double center_x = std::fmod(x0 + background_u * time, length);
    const double center_y = std::fmod(y0 + background_v * time, length);
    const auto file = read_fields(path);
    double density_l1 = 0.0;
    double density_l2 = 0.0;
    double density_linf = 0.0;
    double maximum_primitive = 0.0;
    std::size_t cells = 0;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        const auto& density = require_field(zone, "Density");
        const auto& velocity_x = require_field(zone, "VelocityX");
        const auto& velocity_y = require_field(zone, "VelocityY");
        const auto& temperature = require_field(zone, "Temperature");
        for (std::size_t cell = 0; cell < zone.cell_centers.size(); ++cell) {
            const double dx = std::remainder(zone.cell_centers[cell][0] - center_x, length);
            const double dy = std::remainder(zone.cell_centers[cell][1] - center_y, length);
            const double radius_squared = dx * dx + dy * dy;
            const double exponential = std::exp(0.5 * (1.0 - radius_squared));
            const double exact_u = background_u - beta * exponential * dy / (2.0 * pi);
            const double exact_v = background_v + beta * exponential * dx / (2.0 * pi);
            const double exact_temperature = 1.0
                - (gamma - 1.0) * mach * mach * beta * beta * std::exp(1.0 - radius_squared)
                    / (8.0 * pi * pi);
            const double exact_density = std::pow(exact_temperature, 1.0 / (gamma - 1.0));
            const double density_error = std::abs(density[cell] - exact_density);
            density_l1 += density_error;
            density_l2 += density_error * density_error;
            density_linf = std::max(density_linf, density_error);
            maximum_primitive = std::max({
                maximum_primitive,
                density_error,
                std::abs(velocity_x[cell] - exact_u),
                std::abs(velocity_y[cell] - exact_v),
                std::abs(temperature[cell] - exact_temperature),
            });
            ++cells;
        }
    }
    if (cells == 0) throw std::runtime_error("vortex output has no cells");
    density_l1 /= static_cast<double>(cells);
    density_l2 = std::sqrt(density_l2 / static_cast<double>(cells));
    if (density_l1 > l1_tolerance) {
        throw std::runtime_error("vortex density L1 error " + std::to_string(density_l1)
                                 + " exceeds tolerance " + std::to_string(l1_tolerance));
    }
    std::cout << std::setprecision(17) << "check=isentropic_vortex cells=" << cells
              << " rho_l1=" << density_l1 << " rho_l2=" << density_l2
              << " rho_linf=" << density_linf << " primitive_linf=" << maximum_primitive
              << " tolerance=" << l1_tolerance << '\n';
}

struct SodState {
    double density = 0.0;
    double velocity = 0.0;
    double pressure = 0.0;
};

std::pair<double, double> pressure_wave(double pressure, const SodState& state, double gamma)
{
    const double sound = std::sqrt(gamma * state.pressure / state.density);
    if (pressure > state.pressure) {
        const double a = 2.0 / ((gamma + 1.0) * state.density);
        const double b = (gamma - 1.0) * state.pressure / (gamma + 1.0);
        const double root = std::sqrt(a / (pressure + b));
        return {
            (pressure - state.pressure) * root,
            root * (1.0 - 0.5 * (pressure - state.pressure) / (pressure + b)),
        };
    }
    const double exponent = (gamma - 1.0) / (2.0 * gamma);
    const double ratio = pressure / state.pressure;
    return {
        2.0 * sound / (gamma - 1.0) * (std::pow(ratio, exponent) - 1.0),
        std::pow(ratio, -(gamma + 1.0) / (2.0 * gamma)) / (state.density * sound),
    };
}

std::pair<double, double> sod_star(const SodState& left, const SodState& right, double gamma)
{
    const double left_sound = std::sqrt(gamma * left.pressure / left.density);
    const double right_sound = std::sqrt(gamma * right.pressure / right.density);
    double pressure
        = std::max(1.0e-14,
                   0.5 * (left.pressure + right.pressure)
                       - 0.125 * (right.velocity - left.velocity) * (left.density + right.density)
                           * (left_sound + right_sound));
    for (int iteration = 0; iteration < 50; ++iteration) {
        const auto left_wave = pressure_wave(pressure, left, gamma);
        const auto right_wave = pressure_wave(pressure, right, gamma);
        const double update = (left_wave.first + right_wave.first + right.velocity - left.velocity)
            / (left_wave.second + right_wave.second);
        const double next = std::max(1.0e-14, pressure - update);
        if (std::abs(next - pressure) <= 1.0e-13 * std::max(1.0, pressure)) {
            pressure = next;
            break;
        }
        pressure = next;
    }
    const auto left_wave = pressure_wave(pressure, left, gamma);
    const auto right_wave = pressure_wave(pressure, right, gamma);
    const double velocity
        = 0.5 * (left.velocity + right.velocity + right_wave.first - left_wave.first);
    return {pressure, velocity};
}

SodState sample_sod(double similarity,
                    const SodState& left,
                    const SodState& right,
                    double gamma,
                    double star_pressure,
                    double star_velocity)
{
    const double left_sound = std::sqrt(gamma * left.pressure / left.density);
    const double right_sound = std::sqrt(gamma * right.pressure / right.density);
    const auto star_density = [&](const SodState& state) {
        const double ratio = star_pressure / state.pressure;
        if (ratio > 1.0) {
            const double g = (gamma - 1.0) / (gamma + 1.0);
            return state.density * (ratio + g) / (g * ratio + 1.0);
        }
        return state.density * std::pow(ratio, 1.0 / gamma);
    };

    if (similarity <= star_velocity) {
        const double density = star_density(left);
        if (star_pressure > left.pressure) {
            const double speed = left.velocity
                - left_sound
                    * std::sqrt((gamma + 1.0) * star_pressure / (2.0 * gamma * left.pressure)
                                + (gamma - 1.0) / (2.0 * gamma));
            return similarity <= speed ? left : SodState {density, star_velocity, star_pressure};
        }
        const double head = left.velocity - left_sound;
        const double star_sound
            = left_sound * std::pow(star_pressure / left.pressure, (gamma - 1.0) / (2.0 * gamma));
        const double tail = star_velocity - star_sound;
        if (similarity <= head) return left;
        if (similarity >= tail) return {density, star_velocity, star_pressure};
        const double velocity
            = 2.0 / (gamma + 1.0) * (left_sound + 0.5 * (gamma - 1.0) * left.velocity + similarity);
        const double sound = 2.0 / (gamma + 1.0)
            * (left_sound + 0.5 * (gamma - 1.0) * (left.velocity - similarity));
        return {
            left.density * std::pow(sound / left_sound, 2.0 / (gamma - 1.0)),
            velocity,
            left.pressure * std::pow(sound / left_sound, 2.0 * gamma / (gamma - 1.0)),
        };
    }

    const double density = star_density(right);
    if (star_pressure > right.pressure) {
        const double speed = right.velocity
            + right_sound
                * std::sqrt((gamma + 1.0) * star_pressure / (2.0 * gamma * right.pressure)
                            + (gamma - 1.0) / (2.0 * gamma));
        return similarity >= speed ? right : SodState {density, star_velocity, star_pressure};
    }
    const double head = right.velocity + right_sound;
    const double star_sound
        = right_sound * std::pow(star_pressure / right.pressure, (gamma - 1.0) / (2.0 * gamma));
    const double tail = star_velocity + star_sound;
    if (similarity >= head) return right;
    if (similarity <= tail) return {density, star_velocity, star_pressure};
    const double velocity
        = 2.0 / (gamma + 1.0) * (-right_sound + 0.5 * (gamma - 1.0) * right.velocity + similarity);
    const double sound
        = 2.0 / (gamma + 1.0) * (right_sound - 0.5 * (gamma - 1.0) * (right.velocity - similarity));
    return {
        right.density * std::pow(sound / right_sound, 2.0 / (gamma - 1.0)),
        velocity,
        right.pressure * std::pow(sound / right_sound, 2.0 * gamma / (gamma - 1.0)),
    };
}

void validate_sod(const std::string& path,
                  double time,
                  double discontinuity,
                  double gamma,
                  double density_l1_tolerance,
                  double position_cell_tolerance)
{
    if (!(time > 0.0) || !(gamma > 1.0) || !(position_cell_tolerance > 0.0)) {
        throw std::invalid_argument("Sod validation parameters are invalid");
    }
    require_tolerance(density_l1_tolerance);
    const SodState left {1.0, 0.0, 1.0};
    const SodState right {0.125, 0.0, 0.1};
    const auto [star_pressure, star_velocity] = sod_star(left, right, gamma);
    const double right_sound = std::sqrt(gamma * right.pressure / right.density);
    const double shock_speed = right.velocity
        + right_sound
            * std::sqrt((gamma + 1.0) * star_pressure / (2.0 * gamma * right.pressure)
                        + (gamma - 1.0) / (2.0 * gamma));
    const double exact_contact = discontinuity + star_velocity * time;
    const double exact_shock = discontinuity + shock_speed * time;

    std::map<double, std::pair<double, std::size_t>> density_by_x;
    double density_l1 = 0.0;
    double velocity_l1 = 0.0;
    double pressure_l1 = 0.0;
    double minimum_density = std::numeric_limits<double>::infinity();
    double minimum_pressure = std::numeric_limits<double>::infinity();
    std::size_t cells = 0;
    const auto file = read_fields(path);
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        const auto& density = require_field(zone, "Density");
        const auto& velocity = require_field(zone, "VelocityX");
        const auto& pressure = require_field(zone, "Pressure");
        for (std::size_t cell = 0; cell < zone.cell_centers.size(); ++cell) {
            const auto exact = sample_sod((zone.cell_centers[cell][0] - discontinuity) / time,
                                          left,
                                          right,
                                          gamma,
                                          star_pressure,
                                          star_velocity);
            density_l1 += std::abs(density[cell] - exact.density);
            velocity_l1 += std::abs(velocity[cell] - exact.velocity);
            pressure_l1 += std::abs(pressure[cell] - exact.pressure);
            minimum_density = std::min(minimum_density, density[cell]);
            minimum_pressure = std::min(minimum_pressure, pressure[cell]);
            auto& average = density_by_x[zone.cell_centers[cell][0]];
            average.first += density[cell];
            ++average.second;
            ++cells;
        }
    }
    if (cells == 0 || density_by_x.size() < 4) {
        throw std::runtime_error("Sod output does not contain a valid section");
    }
    density_l1 /= static_cast<double>(cells);
    velocity_l1 /= static_cast<double>(cells);
    pressure_l1 /= static_cast<double>(cells);
    if (minimum_density <= 0.0 || minimum_pressure <= 0.0) {
        throw std::runtime_error("Sod output contains a non-positive state");
    }

    std::vector<std::pair<double, double>> section;
    for (const auto& [x, sum_count] : density_by_x) {
        section.push_back({x, sum_count.first / static_cast<double>(sum_count.second)});
    }
    const double split = 0.5 * (exact_contact + exact_shock);
    double contact_gradient = -1.0;
    double shock_gradient = -1.0;
    double observed_contact = 0.0;
    double observed_shock = 0.0;
    double minimum_spacing = std::numeric_limits<double>::infinity();
    for (std::size_t index = 1; index < section.size(); ++index) {
        const double midpoint = 0.5 * (section[index - 1].first + section[index].first);
        const double spacing = section[index].first - section[index - 1].first;
        const double gradient = std::abs(section[index].second - section[index - 1].second);
        minimum_spacing = std::min(minimum_spacing, spacing);
        if (midpoint > discontinuity && midpoint < split && gradient > contact_gradient) {
            contact_gradient = gradient;
            observed_contact = midpoint;
        }
        if (midpoint >= split && gradient > shock_gradient) {
            shock_gradient = gradient;
            observed_shock = midpoint;
        }
    }
    const double position_error_cells = std::max(std::abs(observed_contact - exact_contact),
                                                 std::abs(observed_shock - exact_shock))
        / minimum_spacing;
    if (density_l1 > density_l1_tolerance || position_error_cells > position_cell_tolerance) {
        throw std::runtime_error("Sod error exceeds the density or wave-position tolerance");
    }
    std::cout << std::setprecision(17) << "check=sod cells=" << cells << " rho_l1=" << density_l1
              << " u_l1=" << velocity_l1 << " p_l1=" << pressure_l1
              << " contact=" << observed_contact << " exact_contact=" << exact_contact
              << " shock=" << observed_shock << " exact_shock=" << exact_shock
              << " position_error_cells=" << position_error_cells << " min_rho=" << minimum_density
              << " min_p=" << minimum_pressure << '\n';
}

void validate_diagonal_symmetry(const std::string& path, double l1_tolerance)
{
    require_tolerance(l1_tolerance);
    struct Sample {
        double density = 0.0;
        double velocity_x = 0.0;
        double velocity_y = 0.0;
        double pressure = 0.0;
    };
    using Key = std::pair<long long, long long>;
    std::map<Key, Sample> samples;
    const auto key = [](double x, double y) {
        constexpr double scale = 1.0e12;
        return Key {std::llround(x * scale), std::llround(y * scale)};
    };
    const auto file = read_fields(path);
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        const auto& density = require_field(zone, "Density");
        const auto& velocity_x = require_field(zone, "VelocityX");
        const auto& velocity_y = require_field(zone, "VelocityY");
        const auto& pressure = require_field(zone, "Pressure");
        for (std::size_t cell = 0; cell < zone.cell_centers.size(); ++cell) {
            if (!samples
                     .emplace(
                         key(zone.cell_centers[cell][0], zone.cell_centers[cell][1]),
                         Sample {density[cell], velocity_x[cell], velocity_y[cell], pressure[cell]})
                     .second) {
                throw std::runtime_error("diagonal symmetry found duplicate cell centers");
            }
        }
    }
    double l1 = 0.0;
    double linf = 0.0;
    for (const auto& [location, sample] : samples) {
        const auto iterator = samples.find({location.second, location.first});
        if (iterator == samples.end()) {
            throw std::runtime_error("diagonal symmetry cannot find reflected cell");
        }
        const auto& reflected = iterator->second;
        const double error = std::max({
            std::abs(sample.density - reflected.density),
            std::abs(sample.pressure - reflected.pressure),
            std::abs(sample.velocity_x - reflected.velocity_y),
            std::abs(sample.velocity_y - reflected.velocity_x),
        });
        l1 += error;
        linf = std::max(linf, error);
    }
    l1 /= static_cast<double>(samples.size());
    if (l1 > l1_tolerance) {
        throw std::runtime_error("diagonal symmetry L1 error exceeds tolerance");
    }
    std::cout << std::setprecision(17) << "check=diagonal_symmetry cells=" << samples.size()
              << " l1=" << l1 << " linf=" << linf << " tolerance=" << l1_tolerance << '\n';
}

void validate_viscous_profile(const std::string& path,
                              const std::string& mode,
                              double reynolds,
                              double l2_tolerance,
                              double pressure_tolerance)
{
    if (mode != "couette" && mode != "conduction") {
        throw std::invalid_argument("viscous profile mode must be couette or conduction");
    }
    if (!(reynolds > 0.0)) {
        throw std::invalid_argument("viscous profile Reynolds number must be positive");
    }
    require_tolerance(l2_tolerance);
    require_tolerance(pressure_tolerance);
    const auto file = read_fields(path);
    double squared = 0.0;
    double maximum = 0.0;
    double pressure_maximum = 0.0;
    double density_maximum = 0.0;
    double sum_y = 0.0;
    double sum_q = 0.0;
    double sum_yy = 0.0;
    double sum_yq = 0.0;
    std::size_t cells = 0;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        const auto& density = require_field(zone, "Density");
        const auto& velocity_x = require_field(zone, "VelocityX");
        const auto& velocity_y = require_field(zone, "VelocityY");
        const auto& pressure = require_field(zone, "Pressure");
        const auto& temperature = require_field(zone, "Temperature");
        for (std::size_t cell = 0; cell < zone.cell_centers.size(); ++cell) {
            const double y = zone.cell_centers[cell][1];
            const double value = mode == "couette" ? velocity_x[cell] : temperature[cell];
            const double exact = mode == "couette" ? y : 1.0 + y;
            const double error = std::abs(value - exact);
            squared += error * error;
            maximum = std::max(maximum, error);
            pressure_maximum = std::max(pressure_maximum, std::abs(pressure[cell] - 1.0));
            if (mode == "conduction") {
                density_maximum
                    = std::max(density_maximum, std::abs(density[cell] - 1.0 / (1.0 + y)));
            }
            maximum = std::max(maximum, std::abs(velocity_y[cell]));
            if (mode == "conduction") {
                maximum = std::max(maximum, std::abs(velocity_x[cell]));
            }
            sum_y += y;
            sum_q += value;
            sum_yy += y * y;
            sum_yq += y * value;
            ++cells;
        }
    }
    if (cells == 0) throw std::runtime_error("viscous profile contains no cells");
    const double count = static_cast<double>(cells);
    const double denominator = count * sum_yy - sum_y * sum_y;
    if (!(denominator > 0.0)) {
        throw std::runtime_error("viscous profile has no wall-normal extent");
    }
    const double slope = (count * sum_yq - sum_y * sum_q) / denominator;
    const double l2 = std::sqrt(squared / count);
    if (l2 > l2_tolerance || maximum > 10.0 * l2_tolerance || pressure_maximum > pressure_tolerance
        || density_maximum > 10.0 * std::max(l2_tolerance, pressure_tolerance)) {
        throw std::runtime_error("viscous analytic profile exceeds tolerance");
    }
    std::cout << std::setprecision(17) << "check=viscous_profile mode=" << mode
              << " cells=" << cells << " l2=" << l2 << " linf=" << maximum
              << " pressure_linf=" << pressure_maximum << " density_linf=" << density_maximum
              << " slope=" << slope;
    if (mode == "couette") {
        std::cout << " wall_shear=" << slope / reynolds << " exact_wall_shear=" << 1.0 / reynolds
                  << " wall_shear_relative_error=" << std::abs(slope - 1.0);
    }
    std::cout << '\n';
}

void validate_poiseuille_profile(const std::string& path,
                                 double y0,
                                 double y1,
                                 double centerline_velocity,
                                 double velocity_l2_tolerance,
                                 double crossflow_tolerance,
                                 double pressure_span_tolerance,
                                 double homogeneity_tolerance)
{
    if (!(y1 > y0) || !std::isfinite(centerline_velocity)) {
        throw std::invalid_argument("Poiseuille profile bounds and velocity are invalid");
    }
    require_tolerance(velocity_l2_tolerance);
    require_tolerance(crossflow_tolerance);
    require_tolerance(pressure_span_tolerance);
    require_tolerance(homogeneity_tolerance);
    struct LayerRange {
        std::array<double, 6> minimum;
        std::array<double, 6> maximum;
        LayerRange()
        {
            minimum.fill(std::numeric_limits<double>::infinity());
            maximum.fill(-std::numeric_limits<double>::infinity());
        }
    };
    std::map<long long, LayerRange> layers;
    const auto file = read_fields(path);
    double weighted_squared_error = 0.0;
    double weighted_velocity = 0.0;
    double total_volume = 0.0;
    double maximum_velocity_error = 0.0;
    double maximum_crossflow = 0.0;
    double minimum_pressure = std::numeric_limits<double>::infinity();
    double maximum_pressure = -std::numeric_limits<double>::infinity();
    double minimum_density = std::numeric_limits<double>::infinity();
    double maximum_density = -std::numeric_limits<double>::infinity();
    double minimum_temperature = std::numeric_limits<double>::infinity();
    double maximum_temperature = -std::numeric_limits<double>::infinity();
    std::size_t cells = 0;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        if (zone.dimension != 3) {
            throw std::runtime_error("Poiseuille profile requires a three-dimensional field");
        }
        const auto& density = require_field(zone, "Density");
        const auto& velocity_x = require_field(zone, "VelocityX");
        const auto& velocity_y = require_field(zone, "VelocityY");
        const auto& velocity_z = require_field(zone, "VelocityZ");
        const auto& pressure = require_field(zone, "Pressure");
        const auto& temperature = require_field(zone, "Temperature");
        const auto& jacobian = require_field(zone, "Jacobian");
        for (std::size_t cell = 0; cell < zone.cell_centers.size(); ++cell) {
            const double eta = (zone.cell_centers[cell][1] - y0) / (y1 - y0);
            if (!(eta > 0.0 && eta < 1.0) || !(jacobian[cell] > 0.0)) {
                throw std::runtime_error("Poiseuille cell center or Jacobian is invalid");
            }
            const double exact_velocity = 4.0 * centerline_velocity * eta * (1.0 - eta);
            const double error = velocity_x[cell] - exact_velocity;
            weighted_squared_error += error * error * jacobian[cell];
            weighted_velocity += velocity_x[cell] * jacobian[cell];
            total_volume += jacobian[cell];
            maximum_velocity_error = std::max(maximum_velocity_error, std::abs(error));
            maximum_crossflow = std::max({
                maximum_crossflow,
                std::abs(velocity_y[cell]),
                std::abs(velocity_z[cell]),
            });
            minimum_pressure = std::min(minimum_pressure, pressure[cell]);
            maximum_pressure = std::max(maximum_pressure, pressure[cell]);
            minimum_density = std::min(minimum_density, density[cell]);
            maximum_density = std::max(maximum_density, density[cell]);
            minimum_temperature = std::min(minimum_temperature, temperature[cell]);
            maximum_temperature = std::max(maximum_temperature, temperature[cell]);
            auto& layer = layers[std::llround(zone.cell_centers[cell][1] * 1.0e12)];
            const std::array<double, 6> values {{
                density[cell],
                velocity_x[cell],
                velocity_y[cell],
                velocity_z[cell],
                pressure[cell],
                temperature[cell],
            }};
            for (std::size_t component = 0; component < values.size(); ++component) {
                if (!std::isfinite(values[component])) {
                    throw std::runtime_error("Poiseuille field contains a non-finite value");
                }
                layer.minimum[component] = std::min(layer.minimum[component], values[component]);
                layer.maximum[component] = std::max(layer.maximum[component], values[component]);
            }
            ++cells;
        }
    }
    if (cells == 0 || !(total_volume > 0.0)) {
        throw std::runtime_error("Poiseuille profile contains no cells");
    }
    double maximum_homogeneity_span = 0.0;
    for (const auto& [y, layer] : layers) {
        static_cast<void>(y);
        for (std::size_t component = 0; component < layer.minimum.size(); ++component) {
            maximum_homogeneity_span = std::max(
                maximum_homogeneity_span, layer.maximum[component] - layer.minimum[component]);
        }
    }
    const double velocity_l2 = std::sqrt(weighted_squared_error / total_volume);
    const double bulk_velocity = weighted_velocity / total_volume;
    const double exact_bulk_velocity = (2.0 / 3.0) * centerline_velocity;
    const double pressure_span = maximum_pressure - minimum_pressure;
    if (velocity_l2 > velocity_l2_tolerance || maximum_crossflow > crossflow_tolerance
        || pressure_span > pressure_span_tolerance
        || maximum_homogeneity_span > homogeneity_tolerance) {
        throw std::runtime_error("Poiseuille analytic profile exceeds tolerance");
    }
    std::cout << std::setprecision(17) << "check=poiseuille_profile cells=" << cells
              << " layers=" << layers.size() << " velocity_l2=" << velocity_l2
              << " velocity_linf=" << maximum_velocity_error << " bulk_velocity=" << bulk_velocity
              << " exact_bulk_velocity=" << exact_bulk_velocity
              << " crossflow_linf=" << maximum_crossflow << " pressure_span=" << pressure_span
              << " homogeneity_linf=" << maximum_homogeneity_span
              << " rho_range=" << minimum_density << ',' << maximum_density
              << " T_range=" << minimum_temperature << ',' << maximum_temperature << '\n';
}

void validate_uniform_source(const std::string& path,
                             double time,
                             const std::array<double, 5>& initial,
                             const std::array<double, 5>& source,
                             double tolerance)
{
    require_tolerance(tolerance);
    const std::array<const char*, 5> names {{
        "Density",
        "MomentumX",
        "MomentumY",
        "MomentumZ",
        "EnergyStagnationDensity",
    }};
    std::array<double, 5> exact {};
    for (std::size_t component = 0; component < exact.size(); ++component) {
        exact[component] = initial[component] + time * source[component];
    }
    const auto file = read_fields(path);
    double maximum = 0.0;
    std::size_t samples = 0;
    for (const auto& [zone_name, zone] : file) {
        static_cast<void>(zone_name);
        for (std::size_t component = 0; component < names.size(); ++component) {
            const auto& values = require_field(zone, names[component]);
            for (const double value : values) {
                maximum = std::max(maximum, std::abs(value - exact[component]));
                ++samples;
            }
        }
    }
    if (samples == 0 || maximum > tolerance) {
        throw std::runtime_error("uniform-source response exceeds tolerance");
    }
    std::cout << std::setprecision(17) << "check=uniform_source samples=" << samples
              << " max_abs=" << maximum << " tolerance=" << tolerance << '\n';
}

} // namespace

int main(int argc, char** argv)
{
    try {
        if (argc == 3 && std::string(argv[1]) == "finite") {
            validate_finite(argv[2]);
        } else if (argc == 4 && std::string(argv[1]) == "field-error") {
            report_field_errors(argv[2], argv[3]);
        } else if (argc == 4 && std::string(argv[1]) == "series-constant") {
            validate_constant_series(argv[2], parse_real(argv[3], "tolerance"));
        } else if (argc == 5 && std::string(argv[1]) == "compare") {
            compare_fields(argv[2], argv[3], parse_real(argv[4], "tolerance"));
        } else if (argc == 7 && std::string(argv[1]) == "compare-scaled-primitive") {
            compare_scaled_primitive_fields(argv[2],
                                            argv[3],
                                            parse_real(argv[4], "gamma"),
                                            parse_real(argv[5], "Mach"),
                                            parse_real(argv[6], "tolerance"));
        } else if (argc == 6 && std::string(argv[1]) == "compare-spatial") {
            compare_spatial_fields(argv[2],
                                   argv[3],
                                   parse_real(argv[4], "field tolerance"),
                                   parse_real(argv[5], "coordinate tolerance"));
        } else if (argc == 9 && std::string(argv[1]) == "uniform") {
            validate_uniform(argv[2],
                             {{"Density", parse_real(argv[3], "rho")},
                              {"VelocityX", parse_real(argv[4], "u")},
                              {"VelocityY", parse_real(argv[5], "v")},
                              {"VelocityZ", parse_real(argv[6], "w")},
                              {"Temperature", parse_real(argv[7], "temperature")}},
                             parse_real(argv[8], "tolerance"));
        } else if (argc == 13 && std::string(argv[1]) == "vortex") {
            validate_isentropic_vortex(argv[2],
                                       parse_real(argv[3], "time"),
                                       parse_real(argv[4], "length"),
                                       parse_real(argv[5], "x0"),
                                       parse_real(argv[6], "y0"),
                                       parse_real(argv[7], "beta"),
                                       parse_real(argv[8], "background_u"),
                                       parse_real(argv[9], "background_v"),
                                       parse_real(argv[10], "gamma"),
                                       parse_real(argv[11], "Mach"),
                                       parse_real(argv[12], "tolerance"));
        } else if (argc == 8 && std::string(argv[1]) == "sod") {
            validate_sod(argv[2],
                         parse_real(argv[3], "time"),
                         parse_real(argv[4], "discontinuity"),
                         parse_real(argv[5], "gamma"),
                         parse_real(argv[6], "density L1 tolerance"),
                         parse_real(argv[7], "position cell tolerance"));
        } else if (argc == 4 && std::string(argv[1]) == "diagonal-symmetry") {
            validate_diagonal_symmetry(argv[2], parse_real(argv[3], "L1 tolerance"));
        } else if (argc == 7 && std::string(argv[1]) == "viscous-profile") {
            validate_viscous_profile(argv[2],
                                     argv[3],
                                     parse_real(argv[4], "Reynolds number"),
                                     parse_real(argv[5], "L2 tolerance"),
                                     parse_real(argv[6], "pressure tolerance"));
        } else if (argc == 10 && std::string(argv[1]) == "poiseuille-profile") {
            validate_poiseuille_profile(argv[2],
                                        parse_real(argv[3], "lower wall coordinate"),
                                        parse_real(argv[4], "upper wall coordinate"),
                                        parse_real(argv[5], "centerline velocity"),
                                        parse_real(argv[6], "velocity L2 tolerance"),
                                        parse_real(argv[7], "crossflow tolerance"),
                                        parse_real(argv[8], "pressure span tolerance"),
                                        parse_real(argv[9], "homogeneity tolerance"));
        } else if (argc == 15 && std::string(argv[1]) == "uniform-source") {
            std::array<double, 5> initial {};
            std::array<double, 5> source {};
            for (std::size_t component = 0; component < 5; ++component) {
                initial[component] = parse_real(argv[4 + component], "initial conservative value");
                source[component] = parse_real(argv[9 + component], "uniform source value");
            }
            validate_uniform_source(argv[2],
                                    parse_real(argv[3], "time"),
                                    initial,
                                    source,
                                    parse_real(argv[14], "tolerance"));
        } else if (argc == 5 && std::string(argv[1]) == "tecplot-consistency") {
            validate_tecplot_consistency(argv[2], argv[3], parse_real(argv[4], "tolerance"));
        } else if (argc == 5 && std::string(argv[1]) == "nonzero") {
            validate_nonzero_field(argv[2], argv[3], parse_real(argv[4], "threshold"));
        } else if (argc == 4 && std::string(argv[1]) == "field-extrema") {
            report_field_extrema(argv[2], argv[3]);
        } else if (argc == 4 && std::string(argv[1]) == "flatplate-profile") {
            report_flatplate_profile(argv[2], parse_real(argv[3], "target x"));
        } else if (argc == 6 && std::string(argv[1]) == "flatplate-yplus") {
            report_flatplate_y_plus(argv[2],
                                    parse_real(argv[3], "target x"),
                                    parse_real(argv[4], "skin-friction coefficient"),
                                    parse_real(argv[5], "Reynolds number"));
        } else if (argc == 7 && std::string(argv[1]) == "derived") {
            validate_derived_fields(argv[2],
                                    parse_real(argv[3], "gamma"),
                                    parse_real(argv[4], "expected viscosity"),
                                    parse_real(argv[5], "expected Jacobian"),
                                    parse_real(argv[6], "tolerance"));
        } else {
            std::cerr << "usage:\n"
                         "  wcns_validate_release_case finite <field.cgns>\n"
                         "  wcns_validate_release_case field-error "
                         "<reference.cgns> <value.cgns>\n"
                         "  wcns_validate_release_case series-constant <series.txt> <tol>\n"
                         "  wcns_validate_release_case compare <lhs.cgns> <rhs.cgns> <tol>\n"
                         "  wcns_validate_release_case compare-scaled-primitive <lhs.cgns> "
                         "<rhs.cgns> <gamma> <Mach> <tol>\n"
                         "  wcns_validate_release_case compare-spatial <lhs.cgns> "
                         "<rhs.cgns> <field-tol> <coordinate-tol>\n"
                         "  wcns_validate_release_case uniform <field.cgns> "
                         "<rho> <u> <v> <w> <T> <tol>\n"
                         "  wcns_validate_release_case vortex <field.cgns> <time> "
                         "<length> <x0> <y0> <beta> <u0> <v0> <gamma> <Mach> "
                         "<L1-tol>\n"
                         "  wcns_validate_release_case sod <field.cgns> <time> "
                         "<x0> <gamma> <rho-L1-tol> <position-cell-tol>\n"
                         "  wcns_validate_release_case diagonal-symmetry "
                         "<field.cgns> <L1-tol>\n"
                         "  wcns_validate_release_case viscous-profile <field.cgns> "
                         "<couette|conduction> <Re> <L2-tol> <pressure-tol>\n"
                         "  wcns_validate_release_case poiseuille-profile <field.cgns> "
                         "<y0> <y1> <centerline-u> <velocity-L2-tol> "
                         "<crossflow-tol> <pressure-span-tol> <homogeneity-tol>\n"
                         "  wcns_validate_release_case uniform-source <field.cgns> "
                         "<time> <U0[5]> <S[5]> <tol>\n"
                         "  wcns_validate_release_case tecplot-consistency "
                         "<field.cgns> <field.dat> <tol>\n"
                         "  wcns_validate_release_case derived <field.cgns> "
                         "<gamma> <viscosity> <Jacobian> <tol>\n"
                         "  wcns_validate_release_case nonzero <field.cgns> "
                         "<field-name> <minimum-maximum-absolute-value>\n"
                         "  wcns_validate_release_case field-extrema <field.cgns> "
                         "<field-name>\n"
                         "  wcns_validate_release_case flatplate-profile <field.cgns> "
                         "<target-x>\n"
                         "  wcns_validate_release_case flatplate-yplus <field.cgns> "
                         "<target-x> <Cf> <Re>\n";
            return EXIT_FAILURE;
        }
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        std::cerr << "release validation failed: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
