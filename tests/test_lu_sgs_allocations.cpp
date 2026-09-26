#include <wcns/solver/implicit_time_integrator.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <new>
#include <stdexcept>

namespace {

bool probe_enabled = false;
std::size_t allocation_count = 0;

void* allocate(std::size_t bytes)
{
    const auto actual = bytes == 0 ? std::size_t {1} : bytes;
    if (void* memory = std::malloc(actual)) {
        if (probe_enabled) ++allocation_count;
        return memory;
    }
    throw std::bad_alloc();
}

template <class Function> std::size_t measure(Function&& function)
{
    allocation_count = 0;
    probe_enabled = true;
    try {
        function();
    } catch (...) {
        probe_enabled = false;
        throw;
    }
    probe_enabled = false;
    return allocation_count;
}

} // namespace

void* operator new(std::size_t bytes)
{
    return allocate(bytes);
}

void* operator new[](std::size_t bytes)
{
    return allocate(bytes);
}

void operator delete(void* memory) noexcept
{
    std::free(memory);
}

void operator delete[](void* memory) noexcept
{
    std::free(memory);
}

void operator delete(void* memory, std::size_t) noexcept
{
    std::free(memory);
}

void operator delete[](void* memory, std::size_t) noexcept
{
    std::free(memory);
}

int main()
{
    try {
        using namespace wcns;
        constexpr int variables = 2;
        const Extent3 extent {16, 8, 1};
        Field<Real> right_hand_side(extent, variables, 0, 0.0);
        Field<Real> diagonal_blocks(extent, variables * variables, 0, 0.0);
        Field<Real> scalar_diagonal(extent, variables, 0, 4.0);
        Field<Real> scalar_coupling(extent, 4, 0, 0.02);
        Field<Real> face_blocks(extent, 4 * variables * variables, 0, 0.0);
        for (int j = 0; j < extent.nj; ++j) {
            for (int i = 0; i < extent.ni; ++i) {
                right_hand_side(i, j, 0, 0) = 1.0 + 0.01 * i;
                right_hand_side(i, j, 0, 1) = 2.0 + 0.01 * j;
                diagonal_blocks(i, j, 0, 0) = 4.0;
                diagonal_blocks(i, j, 0, 1) = 0.1;
                diagonal_blocks(i, j, 0, 2) = 0.2;
                diagonal_blocks(i, j, 0, 3) = 3.5;
                for (int face = 0; face < 4; ++face) {
                    const int offset = face * variables * variables;
                    face_blocks(i, j, 0, offset) = 0.02;
                    face_blocks(i, j, 0, offset + 3) = 0.02;
                }
            }
        }
        LuSgsIterationConfig config;
        config.sweeps = 2;

        Real checksum = 0.0;
        const auto dense_allocations = measure([&] {
            const auto result = solve_block_lu_sgs(
                right_hand_side, diagonal_blocks, scalar_coupling, 2, config);
            checksum += result(7, 3, 0, 0) + result(7, 3, 0, 1);
        });
        const auto face_allocations = measure([&] {
            const auto result = solve_face_block_lu_sgs(
                right_hand_side, scalar_diagonal, face_blocks, 2, config);
            checksum += result(7, 3, 0, 0) + result(7, 3, 0, 1);
        });

        // Solver-owned fields and fixed scratch arrays may allocate at entry; cell and sweep
        // counts must not affect heap traffic.
        constexpr std::size_t maximum_fixed_allocations = 16;
        if (dense_allocations > maximum_fixed_allocations
            || face_allocations > maximum_fixed_allocations || !std::isfinite(checksum)) {
            throw std::runtime_error("LU-SGS allocation probe exceeded its fixed budget");
        }
        std::cout << "lu_sgs_allocation_probe"
                  << " cells=" << extent.size() << " sweeps=" << config.sweeps
                  << " dense_allocations=" << dense_allocations
                  << " face_allocations=" << face_allocations << '\n';
        return EXIT_SUCCESS;
    } catch (const std::exception& error) {
        probe_enabled = false;
        std::cerr << error.what() << '\n';
        return EXIT_FAILURE;
    }
}
