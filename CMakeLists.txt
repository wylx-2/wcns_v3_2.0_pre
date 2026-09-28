cmake_minimum_required(VERSION 3.20)

project(wcns VERSION 2.1.0 LANGUAGES C CXX)

set(WCNS_PROGRAM_VERSION "2.1" CACHE STRING "WCNS user-visible release label")

include(GNUInstallDirs)

execute_process(
    COMMAND git rev-parse --short=12 HEAD
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    OUTPUT_VARIABLE WCNS_GIT_COMMIT
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_QUIET
)
if(NOT WCNS_GIT_COMMIT AND EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/WCNS_SOURCE_REVISION")
    file(STRINGS
        "${CMAKE_CURRENT_SOURCE_DIR}/WCNS_SOURCE_REVISION"
        WCNS_GIT_COMMIT
        LIMIT_COUNT 1
    )
    string(STRIP "${WCNS_GIT_COMMIT}" WCNS_GIT_COMMIT)
endif()
if(NOT WCNS_GIT_COMMIT)
    set(WCNS_GIT_COMMIT unknown)
endif()

option(WCNS_BUILD_TESTS "Build WCNS unit tests" OFF)
option(WCNS_ENABLE_CGNS "Enable CGNS mesh I/O" ON)
option(WCNS_ENABLE_MPI "Enable MPI parallel execution" OFF)
option(WCNS_INSTALL_EXAMPLES "Install example and case templates" OFF)

add_library(wcns_core INTERFACE)
add_library(wcns::core ALIAS wcns_core)

target_include_directories(
    wcns_core
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)
target_compile_features(wcns_core INTERFACE cxx_std_20)

if(MSVC)
    target_compile_options(wcns_core INTERFACE /W4 /permissive-)
else()
    target_compile_options(wcns_core INTERFACE -Wall -Wextra -Wpedantic)
endif()

add_library(wcns_physics STATIC
    src/physics/double_mach_reflection.cpp
    src/physics/source_terms.cpp
    src/physics/thermodynamics.cpp
)
add_library(wcns::physics ALIAS wcns_physics)
target_link_libraries(wcns_physics PUBLIC wcns::core)

add_library(wcns_mesh STATIC
    src/mesh/algorithm_profile.cpp
    src/mesh/conservation_weights.cpp
    src/mesh/geometry_halo.cpp
    src/mesh/halo_exchange.cpp
    src/mesh/high_order_metrics.cpp
    src/mesh/linear_operators.cpp
    src/mesh/metrics.cpp
    src/mesh/structured_mesh.cpp
    src/mesh/topology.cpp
)
add_library(wcns::mesh ALIAS wcns_mesh)
target_link_libraries(wcns_mesh PUBLIC wcns::core wcns::physics)

set(WCNS_MPI_TARGET "")
if(WCNS_ENABLE_MPI)
    if(WIN32 AND MINGW)
        find_path(
            WCNS_INTEL_MPI_INCLUDE_DIR
            mpi.h
            HINTS
                "$ENV{I_MPI_ROOT}/include"
                "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/include"
        )
        find_file(
            WCNS_INTEL_MPI_LIBRARY
            impi.lib
            HINTS
                "$ENV{I_MPI_ROOT}/lib/release"
                "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/lib/release"
        )
        find_program(
            WCNS_MPIEXEC_EXECUTABLE
            NAMES mpiexec mpiexec.exe
            HINTS
                "$ENV{I_MPI_ROOT}/bin"
                "C:/Program Files (x86)/Intel/oneAPI/mpi/latest/bin"
        )
        if(NOT WCNS_INTEL_MPI_INCLUDE_DIR OR NOT WCNS_INTEL_MPI_LIBRARY)
            message(FATAL_ERROR "WCNS_ENABLE_MPI requires an Intel MPI installation")
        endif()
        add_library(wcns_intel_mpi UNKNOWN IMPORTED)
        set_target_properties(
            wcns_intel_mpi
            PROPERTIES
                IMPORTED_LOCATION "${WCNS_INTEL_MPI_LIBRARY}"
                INTERFACE_INCLUDE_DIRECTORIES "${WCNS_INTEL_MPI_INCLUDE_DIR}"
        )
        set(WCNS_MPI_TARGET wcns_intel_mpi)
    else()
        find_package(MPI REQUIRED COMPONENTS CXX)
        set(WCNS_MPI_TARGET MPI::MPI_CXX)
        set(WCNS_MPIEXEC_EXECUTABLE "${MPIEXEC_EXECUTABLE}")
    endif()
endif()

add_library(wcns_parallel STATIC
    src/parallel/block_distribution.cpp
    src/parallel/distributed_topology.cpp
    src/parallel/halo_exchanger.cpp
    src/parallel/mpi_runtime.cpp
)
add_library(wcns::parallel ALIAS wcns_parallel)
target_link_libraries(wcns_parallel PUBLIC wcns::mesh)
if(WCNS_ENABLE_MPI)
    target_link_libraries(wcns_parallel PUBLIC ${WCNS_MPI_TARGET})
    target_compile_definitions(wcns_parallel PUBLIC WCNS_HAS_MPI=1)
    if(MINGW)
        target_compile_definitions(wcns_parallel PUBLIC "__int64=long long")
    endif()
else()
    target_compile_definitions(wcns_parallel PUBLIC WCNS_HAS_MPI=0)
endif()

if(WCNS_ENABLE_CGNS)
    include(FetchContent)
    if(POLICY CMP0135)
        cmake_policy(SET CMP0135 NEW)
    endif()

    set(CGNS_ENABLE_HDF5 OFF CACHE BOOL "" FORCE)
    set(CGNS_ENABLE_64BIT OFF CACHE BOOL "" FORCE)
    set(CGNS_BUILD_SHARED OFF CACHE BOOL "" FORCE)
    set(CGNS_ENABLE_FORTRAN OFF CACHE BOOL "" FORCE)
    set(CGNS_ENABLE_TESTS OFF CACHE BOOL "" FORCE)
    set(CGNS_BUILD_TESTING OFF CACHE BOOL "" FORCE)

    FetchContent_Declare(
        wcns_cgns_source
        URL ${CMAKE_CURRENT_SOURCE_DIR}/third_party/cgns/CGNS-4.4.0.zip
        URL_HASH SHA256=72AD499936089102C1D77E3CDD78A623D61D1C6AD94C4DF3B23F36DB9DDFAD6A
    )
    FetchContent_MakeAvailable(wcns_cgns_source)

    add_library(wcns_cgns INTERFACE)
    add_library(wcns::cgns ALIAS wcns_cgns)
    target_link_libraries(wcns_cgns INTERFACE CGNS::cgns-static)
    if(MINGW)
        target_compile_definitions(wcns_cgns INTERFACE "__int64=long long")
    endif()

    add_library(wcns_io STATIC
        src/io/cgns_reader.cpp
    )
    add_library(wcns::io ALIAS wcns_io)
    target_link_libraries(wcns_io PUBLIC wcns::mesh PRIVATE wcns::cgns)
endif()

add_library(wcns_solver STATIC
    src/solver/boundary_conditions.cpp
    src/solver/euler.cpp
    src/solver/euler_solver.cpp
    src/solver/inviscid_flux.cpp
    src/solver/inviscid_wcns_solver.cpp
    src/solver/physical_boundary.cpp
    src/solver/riemann_solver.cpp
    src/solver/robustness.cpp
    src/solver/rans_two_equation_transport.cpp
    src/solver/sa_negative.cpp
    src/solver/sa_negative_transport.cpp
    src/solver/source_operator.cpp
    src/solver/spatial_operator.cpp
    src/solver/time_integrator.cpp
    src/solver/implicit_time_integrator.cpp
    src/solver/les_closure.cpp
    src/solver/les_filter.cpp
    src/solver/les_models.cpp
    src/solver/low_mach_preconditioner.cpp
    src/solver/two_equation_models.cpp
    src/solver/transport_model.cpp
    src/solver/turbulence_fields.cpp
    src/solver/turbulence_model.cpp
    src/solver/turbulence_transport.cpp
    src/solver/viscous_flux.cpp
    src/solver/viscous_gradient.cpp
    src/solver/viscous_halo.cpp
    src/solver/viscous_boundary.cpp
    src/solver/viscous_operator.cpp
    src/solver/viscous_wcns_solver.cpp
    src/solver/wall_distance.cpp
    src/solver/wcns_reconstruction.cpp
)
add_library(wcns::solver ALIAS wcns_solver)
target_link_libraries(wcns_solver PUBLIC wcns::parallel wcns::physics)

add_library(wcns_runtime STATIC
    src/runtime/boundary_output.cpp
    src/runtime/case_config.cpp
    src/runtime/flow_initializer.cpp
    src/runtime/output_manager.cpp
    src/runtime/quantity_registry.cpp
    src/runtime/simulation_driver.cpp
    src/runtime/stop_controller.cpp
    src/runtime/structured_partition.cpp
    src/runtime/weighted_statistics.cpp
)
add_library(wcns::runtime ALIAS wcns_runtime)
target_link_libraries(wcns_runtime PUBLIC wcns::solver)
target_compile_definitions(
    wcns_runtime
    PRIVATE
        WCNS_GIT_COMMIT="${WCNS_GIT_COMMIT}"
        WCNS_PROGRAM_VERSION="${WCNS_PROGRAM_VERSION}"
)

if(WCNS_ENABLE_CGNS)
    target_sources(wcns_runtime PRIVATE
        src/runtime/checkpoint.cpp
        src/runtime/field_output.cpp
    )
    target_link_libraries(wcns_runtime PRIVATE wcns::io wcns::cgns)
    add_executable(wcns_run src/app/wcns_run.cpp)
    target_link_libraries(wcns_run PRIVATE wcns::runtime wcns::io)

    add_executable(
        wcns_generate_release_cgns
        tools/generate_release_cgns.cpp
    )
    target_link_libraries(wcns_generate_release_cgns PRIVATE wcns::core wcns::cgns)

    add_executable(
        wcns_generate_couette_validation_cgns
        tools/generate_couette_validation_cgns.cpp
    )
    target_link_libraries(
        wcns_generate_couette_validation_cgns PRIVATE wcns::core wcns::cgns
    )

    add_executable(
        wcns_validate_release_case
        tools/validate_release_case.cpp
    )
    target_link_libraries(wcns_validate_release_case PRIVATE wcns::core wcns::cgns)

    add_executable(
        wcns_compare_metric_profiles
        tools/compare_metric_profiles.cpp
    )
    target_link_libraries(wcns_compare_metric_profiles PRIVATE wcns::io)

    add_executable(
        wcns_inspect_structured_mesh
        tools/inspect_structured_mesh.cpp
    )
    target_link_libraries(wcns_inspect_structured_mesh PRIVATE wcns::io)

    add_executable(
        wcns_convert_tmr_p2d_to_cgns
        tools/convert_tmr_p2d_to_cgns.cpp
    )
    target_link_libraries(wcns_convert_tmr_p2d_to_cgns PRIVATE wcns::core wcns::cgns)

    add_executable(
        wcns_convert_tmr_flatplate_p3d_to_cgns
        tools/convert_tmr_flatplate_p3d_to_cgns.cpp
    )
    target_link_libraries(
        wcns_convert_tmr_flatplate_p3d_to_cgns PRIVATE wcns::core wcns::cgns
    )

    add_executable(
        wcns_extract_rans_profile
        tools/extract_rans_profile.cpp
    )
    target_link_libraries(wcns_extract_rans_profile PRIVATE wcns::core wcns::cgns)
endif()

if(WCNS_BUILD_TESTS)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tests/CMakeLists.txt")
        message(FATAL_ERROR
            "WCNS_BUILD_TESTS requires the full development source tree; "
            "the compact release directory intentionally omits tests"
        )
    endif()
    include(CTest)
    add_subdirectory(tests)
endif()

if(WCNS_ENABLE_CGNS)
    install(
        TARGETS
            wcns_run
            wcns_generate_release_cgns
            wcns_generate_couette_validation_cgns
            wcns_validate_release_case
            wcns_compare_metric_profiles
            wcns_inspect_structured_mesh
            wcns_convert_tmr_p2d_to_cgns
            wcns_convert_tmr_flatplate_p3d_to_cgns
            wcns_extract_rans_profile
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif()
install(
    FILES
        README.md
        LICENSE.md
        THIRD_PARTY_NOTICES.md
        "算法补充.md"
    DESTINATION ${CMAKE_INSTALL_DATADIR}/wcns
)
install(
    DIRECTORY docs/
    DESTINATION ${CMAKE_INSTALL_DATADIR}/wcns/docs
    FILES_MATCHING PATTERN "*.md"
)
install(
    FILES examples/channel_retau180_from_t250_all_speed_roe.wcns
    DESTINATION ${CMAKE_INSTALL_DATADIR}/wcns/config-templates
)
if(WCNS_INSTALL_EXAMPLES)
    if(NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/examples"
       OR NOT EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/cases/config")
        message(FATAL_ERROR
            "WCNS_INSTALL_EXAMPLES requires the full development source tree; "
            "the compact release directory intentionally omits full examples and cases"
        )
    endif()
    install(
        DIRECTORY examples/
        DESTINATION ${CMAKE_INSTALL_DATADIR}/wcns/examples
    )
    install(
        DIRECTORY cases/config/
        DESTINATION ${CMAKE_INSTALL_DATADIR}/wcns/cases/config
        FILES_MATCHING PATTERN "*.wcns.in"
    )
endif()
