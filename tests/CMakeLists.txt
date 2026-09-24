add_executable(wcns_unit_tests
    test_main.cpp
    test_index.cpp
    test_array3d.cpp
    test_field.cpp
    test_geometry_operators.cpp
    test_geometry_halo.cpp
    test_high_order_metrics.cpp
    test_topology_field.cpp
    test_structured_block.cpp
    test_cgns_link.cpp
    test_conservation_weights.cpp
    test_metrics.cpp
    test_topology.cpp
    test_distribution.cpp
    test_euler.cpp
    test_thermodynamics.cpp
    test_wcns.cpp
    test_spatial_operator.cpp
    test_implicit_time_integrator.cpp
    test_source_terms.cpp
    test_physical_boundary.cpp
    test_double_mach_reflection.cpp
    test_inviscid_reconstruction.cpp
    test_inviscid_flux.cpp
    test_stage_j_sources.cpp
    test_viscous_flux.cpp
    test_viscous_gradient.cpp
    test_viscous_boundary.cpp
    test_robustness.cpp
    test_case_config.cpp
    test_simulation_driver.cpp
    test_stop_controller.cpp
    test_structured_partition.cpp
    test_flow_initializer.cpp
    test_output_manager.cpp
    test_boundary_output.cpp
    test_quantity_registry.cpp
    test_turbulence_infrastructure.cpp
    test_stage_z_les.cpp
)

target_link_libraries(
    wcns_unit_tests
    PRIVATE wcns::core wcns::mesh wcns::parallel wcns::solver wcns::runtime
)

add_executable(wcns_stage_t_allocation_probe test_stage_t_allocations.cpp)
target_link_libraries(
    wcns_stage_t_allocation_probe
    PRIVATE wcns::core wcns::mesh wcns::parallel wcns::solver
)
add_test(NAME wcns.stage_t.allocations COMMAND wcns_stage_t_allocation_probe)
if(WCNS_ENABLE_CGNS)
    target_link_libraries(wcns_unit_tests PRIVATE wcns::cgns)
    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    add_test(
        NAME wcns.validation_budget
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_validation_budget.py
            --checker ${CMAKE_CURRENT_SOURCE_DIR}/../tools/check_validation_budget.py
    )
    add_test(
        NAME wcns.case06_preparation
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_prepare_case06_naca0012.py
            --generator ${CMAKE_CURRENT_SOURCE_DIR}/../tools/prepare_case06_naca0012.py
            --manifest ${CMAKE_CURRENT_SOURCE_DIR}/../cases/validation/case06_naca0012/manifest.json
    )
    add_test(
        NAME wcns.rans_analyzer
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/test_rans_analyzer.py
            --analyzer ${CMAKE_CURRENT_SOURCE_DIR}/../tools/analyze_rans_run.py
    )

    add_executable(wcns_generate_test_cgns
        ${CMAKE_CURRENT_SOURCE_DIR}/../tools/generate_test_cgns.cpp
    )
    target_link_libraries(wcns_generate_test_cgns PRIVATE wcns::core wcns::cgns)

    add_test(
        NAME wcns.generate_test_cgns
        COMMAND wcns_generate_test_cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_3d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_invalid_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_3d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_one_sided_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_unknown_donor_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_partition_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_periodic_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_periodic_3d.cgns
    )
    set_tests_properties(wcns.generate_test_cgns PROPERTIES FIXTURES_SETUP cgns_files)

    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/freestream_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/freestream_case.wcns
        @ONLY
    )
    set(WCNS_IMPLICIT_CASE implicit-freestream-steady)
    set(WCNS_IMPLICIT_MODE steady)
    set(WCNS_IMPLICIT_RUN_LINES "# steady pseudo-time")
    set(WCNS_IMPLICIT_PRECONDITIONER none)
    set(WCNS_IMPLICIT_PRECONDITIONER_LINES "# no preconditioner cutoffs")
    set(WCNS_IMPLICIT_CHECKPOINT false)
    set(WCNS_IMPLICIT_RESTART_LINE "# fresh implicit run")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_freestream_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/implicit_freestream_steady.wcns
        @ONLY
    )
    set(WCNS_IMPLICIT_CASE implicit-freestream-dual-time)
    set(WCNS_IMPLICIT_MODE unsteady)
    set(WCNS_IMPLICIT_RUN_LINES "run.t_end = 0.002\ntime.physical.step = 0.001")
    set(WCNS_IMPLICIT_PRECONDITIONER weiss_smith)
    set(WCNS_IMPLICIT_PRECONDITIONER_LINES
        "preconditioner.mach_cutoff = 0.001\npreconditioner.viscous_cutoff = 1")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_freestream_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/implicit_freestream_dual_time.wcns
        @ONLY
    )
    set(WCNS_IMPLICIT_CASE implicit-checkpoint-source)
    set(WCNS_IMPLICIT_RUN_LINES "run.t_end = 0.001\ntime.physical.step = 0.001")
    set(WCNS_IMPLICIT_CHECKPOINT true)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_freestream_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/implicit_checkpoint_source.wcns
        @ONLY
    )
    set(WCNS_IMPLICIT_CASE implicit-checkpoint-restart)
    set(WCNS_IMPLICIT_RUN_LINES "run.t_end = 0.002\ntime.physical.step = 0.001")
    set(WCNS_IMPLICIT_CHECKPOINT false)
    set(WCNS_IMPLICIT_RESTART_LINE
        "restart.path = ${CMAKE_CURRENT_BINARY_DIR}/implicit-checkpoint-source/implicit-checkpoint-source.checkpoint.latest.cgns")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_freestream_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/implicit_checkpoint_restart.wcns
        @ONLY
    )
    foreach(les_model IN ITEMS
            smagorinsky
            scale_similarity
            mixed_smagorinsky_similarity
            dynamic_smagorinsky
            wale)
        set(WCNS_LES_MODEL ${les_model})
        if(les_model STREQUAL "smagorinsky")
            set(WCNS_LES_MODEL_LINES
                "les.smagorinsky.cs = 0.17\nles.smagorinsky.wall_damping = none")
        elseif(les_model STREQUAL "scale_similarity")
            set(WCNS_LES_MODEL_LINES
                "les.test_filter.ratio = 2\nles.similarity.cb = 1")
        elseif(les_model STREQUAL "mixed_smagorinsky_similarity")
            set(WCNS_LES_MODEL_LINES
                "les.test_filter.ratio = 2\nles.smagorinsky.cs = 0.17\nles.smagorinsky.wall_damping = none\nles.similarity.cb = 1")
        elseif(les_model STREQUAL "dynamic_smagorinsky")
            set(WCNS_LES_MODEL_LINES
                "les.test_filter.ratio = 2\nles.dynamic.average = local_box3_tensor\nles.dynamic.clipping = bounded\nles.dynamic.denominator_floor = 1e-20\nles.dynamic.coefficient_minimum = -0.05\nles.dynamic.coefficient_maximum = 0.09")
        else()
            set(WCNS_LES_MODEL_LINES "les.wale.cw = 0.325")
        endif()
        foreach(les_rank_label IN ITEMS serial mpi2 mpi4)
            set(WCNS_LES_RANK_LABEL ${les_rank_label})
            configure_file(
                ${CMAKE_CURRENT_SOURCE_DIR}/data/les_freestream_case.wcns.in
                ${CMAKE_CURRENT_BINARY_DIR}/les_${les_model}_${les_rank_label}.wcns
                @ONLY
            )
            configure_file(
                ${CMAKE_CURRENT_SOURCE_DIR}/data/les_implicit_convergence_case.wcns.in
                ${CMAKE_CURRENT_BINARY_DIR}/les_implicit_${les_model}_${les_rank_label}.wcns
                @ONLY
            )
        endforeach()
    endforeach()
    set(WCNS_LES_ZERO_CASE les-zero-none)
    set(WCNS_LES_ZERO_OUTPUT les-zero-none)
    set(WCNS_LES_ZERO_MODEL "turbulence.model = none")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/les_zero_equivalence_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/les_zero_none.wcns
        @ONLY
    )
    set(WCNS_LES_ZERO_CASE les-zero-smagorinsky)
    set(WCNS_LES_ZERO_OUTPUT les-zero-smagorinsky)
    set(WCNS_LES_ZERO_MODEL
        "turbulence.model = smagorinsky\nles.filter.type = box3_tensor\nles.filter.width_ratio = 1\nles.sgs_prandtl = 0.9\nles.smagorinsky.cs = 0\nles.smagorinsky.wall_damping = none")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/les_zero_equivalence_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/les_zero_smagorinsky.wcns
        @ONLY
    )
    set(WCNS_LES_STATS_CASE les-statistics-source)
    set(WCNS_LES_STATS_MAX_STEPS 1)
    set(WCNS_LES_STATS_END_TIME 0.001)
    set(WCNS_LES_STATS_OUTPUT les-statistics-source)
    set(WCNS_LES_STATS_CHECKPOINT true)
    set(WCNS_LES_STATS_RESTART "# fresh statistics run")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/les_statistics_restart_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/les_statistics_source.wcns
        @ONLY
    )
    foreach(les_stats_rank_label IN ITEMS serial mpi2 mpi4)
        set(WCNS_LES_STATS_CASE les-statistics-continuous-${les_stats_rank_label})
        set(WCNS_LES_STATS_MAX_STEPS 2)
        set(WCNS_LES_STATS_END_TIME 0.002)
        set(WCNS_LES_STATS_OUTPUT les-statistics-continuous-${les_stats_rank_label})
        set(WCNS_LES_STATS_CHECKPOINT false)
        set(WCNS_LES_STATS_RESTART "# fresh statistics run")
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/data/les_statistics_restart_case.wcns.in
            ${CMAKE_CURRENT_BINARY_DIR}/les_statistics_continuous_${les_stats_rank_label}.wcns
            @ONLY
        )
        set(WCNS_LES_STATS_CASE les-statistics-restart-${les_stats_rank_label})
        set(WCNS_LES_STATS_OUTPUT les-statistics-restart-${les_stats_rank_label})
        set(WCNS_LES_STATS_RESTART
            "restart.path = ${CMAKE_CURRENT_BINARY_DIR}/les-statistics-source/les-statistics-source.checkpoint.latest.cgns")
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/data/les_statistics_restart_case.wcns.in
            ${CMAKE_CURRENT_BINARY_DIR}/les_statistics_restart_${les_stats_rank_label}.wcns
            @ONLY
        )
    endforeach()
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_sa_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/implicit_sa_steady.wcns
        @ONLY
    )
    set(WCNS_TWO_EQUATION_IMPLICIT_TIME_LINES
        "time.integrator = lu_sgs\ntime.physical.scheme = bdf2\ntime.dual_time.max_iterations = 20\n# Restart smoke accepts the initialized physical layer; convergence quality is\n# covered by the dedicated RANS validation cases rather than this I/O fixture.\ntime.dual_time.absolute_tolerance = 1e6\ntime.dual_time.relative_tolerance = 1e-8\ntime.dual_time.cfl = 0.01\nlu_sgs.sweeps = 1\nlu_sgs.jacobian = scalar_spectral\nlu_sgs.relaxation = 0.5")
    set(WCNS_TWO_EQUATION_TIME_LINES ${WCNS_TWO_EQUATION_IMPLICIT_TIME_LINES})
    set(WCNS_TWO_EQUATION_MIN_STEPS 1)
    set(WCNS_TWO_EQUATION_CASE sst-steady)
    set(WCNS_TWO_EQUATION_OUTPUT_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/implicit-sst-steady)
    set(WCNS_TWO_EQUATION_MODEL k_omega_sst)
    set(WCNS_TWO_EQUATION_WALL resolved)
    set(WCNS_TWO_EQUATION_WALL_LINES "# no wall-function range")
    set(WCNS_TWO_EQUATION_MODE steady)
    set(WCNS_TWO_EQUATION_RUN_LINES "# steady pseudo-time")
    set(WCNS_TWO_EQUATION_CHECKPOINT false)
    set(WCNS_TWO_EQUATION_RESTART_LINE "# fresh two-equation run")
    set(WCNS_TWO_EQUATION_MAX_STEPS 1)
    set(WCNS_TWO_EQUATION_QUANTITIES
        "rho,u,p,k,omega,mu_t_over_mu,sst_f1,sst_f2,turbulence_production,turbulence_destruction,cross_diffusion,wall_distance,wall_y_plus")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_two_equation_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/implicit_sst_steady.wcns
        @ONLY
    )
    set(WCNS_TWO_EQUATION_CASE k-epsilon-steady)
    set(WCNS_TWO_EQUATION_OUTPUT_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/implicit-k-epsilon-steady)
    set(WCNS_TWO_EQUATION_MODEL k_epsilon)
    set(WCNS_TWO_EQUATION_WALL wall_function)
    set(WCNS_TWO_EQUATION_WALL_LINES
        "turbulence.experimental = true\nturbulence.wall_function.y_plus_min = 30\nturbulence.wall_function.y_plus_max = 300")
    set(WCNS_TWO_EQUATION_QUANTITIES
        "rho,u,p,k,epsilon,mu_t_over_mu,turbulence_production,turbulence_destruction,wall_distance,wall_y_plus")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_two_equation_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/implicit_k_epsilon_steady.wcns
        @ONLY
    )
    set(WCNS_TWO_EQUATION_TIME_LINES "time.integrator = ssprk3")
    set(WCNS_TWO_EQUATION_CHECKPOINT false)
    set(WCNS_TWO_EQUATION_RESTART_LINE "# fresh explicit two-equation run")
    foreach(two_equation_model IN ITEMS sst k_epsilon)
        if(two_equation_model STREQUAL "sst")
            set(WCNS_TWO_EQUATION_MODEL k_omega_sst)
            set(WCNS_TWO_EQUATION_WALL resolved)
            set(WCNS_TWO_EQUATION_WALL_LINES "# no wall-function range")
            set(WCNS_TWO_EQUATION_QUANTITIES
                "rho,u,p,k,omega,mu_t_over_mu,sst_f1,sst_f2,turbulence_production,turbulence_destruction,cross_diffusion,wall_distance,wall_y_plus")
        else()
            set(WCNS_TWO_EQUATION_MODEL k_epsilon)
            set(WCNS_TWO_EQUATION_WALL wall_function)
            set(WCNS_TWO_EQUATION_WALL_LINES
                "turbulence.experimental = true\nturbulence.wall_function.y_plus_min = 30\nturbulence.wall_function.y_plus_max = 300")
            set(WCNS_TWO_EQUATION_QUANTITIES
                "rho,u,p,k,epsilon,mu_t_over_mu,turbulence_production,turbulence_destruction,wall_distance,wall_y_plus")
        endif()
        set(WCNS_TWO_EQUATION_CASE ${two_equation_model}-explicit-steady)
        set(WCNS_TWO_EQUATION_OUTPUT_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}/implicit-${two_equation_model}-explicit-steady)
        set(WCNS_TWO_EQUATION_MODE steady)
        set(WCNS_TWO_EQUATION_RUN_LINES "# steady explicit pseudo-time")
        set(WCNS_TWO_EQUATION_MAX_STEPS 1)
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_two_equation_case.wcns.in
            ${CMAKE_CURRENT_BINARY_DIR}/explicit_${two_equation_model}_steady.wcns
            @ONLY
        )
    endforeach()
    set(WCNS_TWO_EQUATION_TIME_LINES ${WCNS_TWO_EQUATION_IMPLICIT_TIME_LINES})
    foreach(two_equation_model IN ITEMS sst k_epsilon)
        if(two_equation_model STREQUAL "sst")
            set(WCNS_TWO_EQUATION_MODEL k_omega_sst)
            set(WCNS_TWO_EQUATION_WALL resolved)
            set(WCNS_TWO_EQUATION_WALL_LINES "# no wall-function range")
            set(WCNS_TWO_EQUATION_QUANTITIES
                "rho,u,p,k,omega,mu_t_over_mu,sst_f1,sst_f2,turbulence_production,turbulence_destruction,cross_diffusion,wall_distance,wall_y_plus")
        else()
            set(WCNS_TWO_EQUATION_MODEL k_epsilon)
            set(WCNS_TWO_EQUATION_WALL wall_function)
            set(WCNS_TWO_EQUATION_WALL_LINES
                "turbulence.experimental = true\nturbulence.wall_function.y_plus_min = 30\nturbulence.wall_function.y_plus_max = 300")
            set(WCNS_TWO_EQUATION_QUANTITIES
                "rho,u,p,k,epsilon,mu_t_over_mu,turbulence_production,turbulence_destruction,wall_distance,wall_y_plus")
        endif()
        set(WCNS_TWO_EQUATION_CASE ${two_equation_model}-checkpoint-source)
        set(WCNS_TWO_EQUATION_OUTPUT_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}/implicit-${two_equation_model}-checkpoint-source)
        set(WCNS_TWO_EQUATION_MODE unsteady)
        set(WCNS_TWO_EQUATION_RUN_LINES "run.t_end = 0.001\ntime.physical.step = 0.001")
        set(WCNS_TWO_EQUATION_CHECKPOINT true)
        set(WCNS_TWO_EQUATION_RESTART_LINE "# fresh two-equation run")
        set(WCNS_TWO_EQUATION_MAX_STEPS 1)
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_two_equation_case.wcns.in
            ${CMAKE_CURRENT_BINARY_DIR}/implicit_${two_equation_model}_checkpoint_source.wcns
            @ONLY
        )
        set(WCNS_TWO_EQUATION_CASE ${two_equation_model}-checkpoint-restart)
        set(WCNS_TWO_EQUATION_OUTPUT_DIRECTORY
            ${CMAKE_CURRENT_BINARY_DIR}/implicit-${two_equation_model}-checkpoint-restart)
        set(WCNS_TWO_EQUATION_RUN_LINES "run.t_end = 0.002\ntime.physical.step = 0.001")
        set(WCNS_TWO_EQUATION_CHECKPOINT false)
        set(WCNS_TWO_EQUATION_MAX_STEPS 2)
        set(WCNS_TWO_EQUATION_RESTART_LINE
            "restart.path = ${CMAKE_CURRENT_BINARY_DIR}/implicit-${two_equation_model}-checkpoint-source/implicit-${two_equation_model}-checkpoint-source.checkpoint.latest.cgns")
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_two_equation_case.wcns.in
            ${CMAKE_CURRENT_BINARY_DIR}/implicit_${two_equation_model}_checkpoint_restart.wcns
            @ONLY
        )
    endforeach()
    set(WCNS_TWO_EQUATION_CASE wrong-model-checkpoint-restart)
    set(WCNS_TWO_EQUATION_OUTPUT_DIRECTORY
        ${CMAKE_CURRENT_BINARY_DIR}/implicit-wrong-model-checkpoint-restart)
    set(WCNS_TWO_EQUATION_MODEL k_epsilon)
    set(WCNS_TWO_EQUATION_WALL wall_function)
    set(WCNS_TWO_EQUATION_WALL_LINES
        "turbulence.experimental = true\nturbulence.wall_function.y_plus_min = 30\nturbulence.wall_function.y_plus_max = 300")
    set(WCNS_TWO_EQUATION_MODE unsteady)
    set(WCNS_TWO_EQUATION_RUN_LINES "run.t_end = 0.002\ntime.physical.step = 0.001")
    set(WCNS_TWO_EQUATION_CHECKPOINT false)
    set(WCNS_TWO_EQUATION_MAX_STEPS 2)
    set(WCNS_TWO_EQUATION_QUANTITIES
        "rho,u,p,k,epsilon,mu_t_over_mu,turbulence_production,turbulence_destruction,wall_distance,wall_y_plus")
    set(WCNS_TWO_EQUATION_RESTART_LINE
        "restart.path = ${CMAKE_CURRENT_BINARY_DIR}/implicit-sst-checkpoint-source/implicit-sst-checkpoint-source.checkpoint.latest.cgns")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_two_equation_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/implicit_wrong_model_checkpoint_restart.wcns
        @ONLY
    )
    set(WCNS_TWO_EQUATION_MODE steady)
    set(WCNS_TWO_EQUATION_RUN_LINES "# rank-consistency pseudo-time")
    set(WCNS_TWO_EQUATION_CHECKPOINT false)
    set(WCNS_TWO_EQUATION_RESTART_LINE "# fresh rank-consistency run")
    # The second step transports the wall-generated nonuniform model field
    # across runtime partition interfaces; a single uniform first step cannot
    # expose a missing gradient/diffusion-coefficient halo.
    set(WCNS_TWO_EQUATION_MAX_STEPS 2)
    set(WCNS_TWO_EQUATION_MIN_STEPS 2)
    foreach(two_equation_model IN ITEMS sst k_epsilon)
        if(two_equation_model STREQUAL "sst")
            set(WCNS_TWO_EQUATION_MODEL k_omega_sst)
            set(WCNS_TWO_EQUATION_WALL resolved)
            set(WCNS_TWO_EQUATION_WALL_LINES "# no wall-function range")
            set(WCNS_TWO_EQUATION_QUANTITIES
                "rho,u,p,k,omega,mu_t_over_mu,sst_f1,sst_f2,turbulence_production,turbulence_destruction,cross_diffusion,wall_distance,wall_y_plus")
        else()
            set(WCNS_TWO_EQUATION_MODEL k_epsilon)
            set(WCNS_TWO_EQUATION_WALL wall_function)
            set(WCNS_TWO_EQUATION_WALL_LINES
                "turbulence.experimental = true\nturbulence.wall_function.y_plus_min = 30\nturbulence.wall_function.y_plus_max = 300")
            set(WCNS_TWO_EQUATION_QUANTITIES
                "rho,u,p,k,epsilon,mu_t_over_mu,turbulence_production,turbulence_destruction,wall_distance,wall_y_plus")
        endif()
        set(WCNS_TWO_EQUATION_CASE ${two_equation_model}-rank-consistency)
        foreach(rank_label IN ITEMS serial mpi2 mpi4)
            set(WCNS_TWO_EQUATION_OUTPUT_DIRECTORY
                ${CMAKE_CURRENT_BINARY_DIR}/two-equation-${two_equation_model}-${rank_label})
            configure_file(
                ${CMAKE_CURRENT_SOURCE_DIR}/data/implicit_two_equation_case.wcns.in
                ${CMAKE_CURRENT_BINARY_DIR}/implicit_${two_equation_model}_${rank_label}.wcns
                @ONLY
            )
        endforeach()
    endforeach()
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/boundary_partition_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/boundary_partition_case.wcns
        @ONLY
    )
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/cylinder_o_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/cylinder_o_case.wcns
        @ONLY
    )
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/case07_boundary_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/case07_boundary_case.wcns
        @ONLY
    )
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/double_mach_reflection_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/double_mach_reflection_case.wcns
        @ONLY
    )
    set(WCNS_XZ_CASE_NAME xz-plane-statistics-smoke)
    set(WCNS_XZ_OUTPUT_DIRECTORY xz-plane-statistics-smoke)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/xz_plane_statistics_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/xz_plane_statistics_case.wcns
        @ONLY
    )
    foreach(les_boundary_rank IN ITEMS serial mpi2 mpi4)
        set(WCNS_LES_BOUNDARY_RANK ${les_boundary_rank})
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/data/les_boundary_case.wcns.in
            ${CMAKE_CURRENT_BINARY_DIR}/les_boundary_${les_boundary_rank}.wcns
            @ONLY
        )
    endforeach()
    if(WCNS_ENABLE_MPI)
        set(WCNS_XZ_CASE_NAME xz-plane-statistics-smoke-mpi)
        set(WCNS_XZ_OUTPUT_DIRECTORY xz-plane-statistics-smoke-mpi)
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/data/xz_plane_statistics_case.wcns.in
            ${CMAKE_CURRENT_BINARY_DIR}/xz_plane_statistics_mpi_case.wcns
            @ONLY
        )
    endif()
    set(WCNS_PERIODIC_CASE periodic-translation-2d)
    set(WCNS_PERIODIC_MESH generated_periodic_2d.cgns)
    set(WCNS_PERIODIC_PARTITION auto_split)
    set(WCNS_PERIODIC_U 0.2)
    set(WCNS_PERIODIC_W 0.0)
    set(WCNS_PERIODIC_BOUNDARY slip_wall)
    set(WCNS_PERIODIC_OUTPUT periodic-translation-2d)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/periodic_freestream_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/periodic_translation_2d.wcns
        @ONLY
    )
    set(WCNS_PERIODIC_CASE periodic-rotation-3d)
    set(WCNS_PERIODIC_MESH generated_periodic_3d.cgns)
    set(WCNS_PERIODIC_PARTITION zones_only)
    set(WCNS_PERIODIC_U 0.0)
    set(WCNS_PERIODIC_W 0.2)
    set(WCNS_PERIODIC_BOUNDARY farfield)
    set(WCNS_PERIODIC_OUTPUT periodic-rotation-3d)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/periodic_freestream_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/periodic_rotation_3d.wcns
        @ONLY
    )
    set(WCNS_VISCOUS_CASE couette-analytic)
    set(WCNS_VISCOUS_INITIAL couette)
    set(WCNS_VISCOUS_LOWER_U 0.0)
    set(WCNS_VISCOUS_UPPER_U 0.2)
    set(WCNS_VISCOUS_LOWER_T 1.0)
    set(WCNS_VISCOUS_UPPER_T 1.0)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/viscous_analytic_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/couette_analytic.wcns
        @ONLY
    )
    set(WCNS_VISCOUS_CASE linear-conduction-analytic)
    set(WCNS_VISCOUS_INITIAL linear_conduction)
    set(WCNS_VISCOUS_LOWER_U 0.0)
    set(WCNS_VISCOUS_UPPER_U 0.0)
    set(WCNS_VISCOUS_LOWER_T 1.0)
    set(WCNS_VISCOUS_UPPER_T 2.0)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/viscous_analytic_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/linear_conduction_analytic.wcns
        @ONLY
    )
    set(WCNS_TRANSPORT_CASE transport-sutherland-serial)
    set(WCNS_TRANSPORT_OUTPUT transport-sutherland-serial)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/transport_runtime_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/transport_sutherland_serial.wcns
        @ONLY
    )
    if(WCNS_ENABLE_MPI)
        set(WCNS_TRANSPORT_CASE transport-sutherland-mpi)
        set(WCNS_TRANSPORT_OUTPUT transport-sutherland-mpi)
        configure_file(
            ${CMAKE_CURRENT_SOURCE_DIR}/data/transport_runtime_case.wcns.in
            ${CMAKE_CURRENT_BINARY_DIR}/transport_sutherland_mpi.wcns
            @ONLY
        )
    endif()
    set(WCNS_OUTPUT_SUBDIR output-validation-serial)
    set(WCNS_HISTORY_FORMAT txt)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/freestream_output_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/freestream_output_serial.wcns
        @ONLY
    )
    set(WCNS_OUTPUT_SUBDIR output-validation-mpi)
    set(WCNS_HISTORY_FORMAT tecplot)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/freestream_output_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/freestream_output_mpi.wcns
        @ONLY
    )
    set(WCNS_CHECKPOINT_OUTPUT checkpoint-validation)
    set(WCNS_CHECKPOINT_GAMMA 1.4)
    set(WCNS_TRANSPORT_LINES "# default v1.0 transport migration")
    set(WCNS_CHECKPOINT_ENABLED true)
    set(WCNS_RESTART_FIELD false)
    set(WCNS_RESTART_LINE "# fresh start")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/checkpoint_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/checkpoint_source.wcns
        @ONLY
    )
    set(WCNS_UNSTEADY_CASE unsteady-freestream)
    set(WCNS_UNSTEADY_OUTPUT unsteady-validation)
    set(WCNS_UNSTEADY_ALLOW_EXISTING true)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/unsteady_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/unsteady_case.wcns
        @ONLY
    )
    set(WCNS_UNSTEADY_CASE existing-output-rejection)
    set(WCNS_UNSTEADY_OUTPUT existing-output-rejection)
    set(WCNS_UNSTEADY_ALLOW_EXISTING false)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/unsteady_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/existing_output_case.wcns
        @ONLY
    )
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/maximum_steps_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/maximum_steps_case.wcns
        @ONLY
    )
    set(WCNS_CHECKPOINT_ENABLED false)
    set(WCNS_RESTART_FIELD true)
    set(WCNS_RESTART_LINE
        "restart.path = ${CMAKE_CURRENT_BINARY_DIR}/checkpoint-validation/checkpoint-freestream.checkpoint.latest.cgns")
    set(WCNS_CHECKPOINT_GAMMA 1.4)
    set(WCNS_CHECKPOINT_OUTPUT restart-validation-serial)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/checkpoint_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/restart_serial.wcns
        @ONLY
    )
    set(WCNS_CHECKPOINT_OUTPUT restart-validation-mpi)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/checkpoint_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/restart_mpi.wcns
        @ONLY
    )
    set(WCNS_CHECKPOINT_GAMMA 1.3)
    set(WCNS_CHECKPOINT_OUTPUT restart-invalid)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/checkpoint_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/restart_invalid.wcns
        @ONLY
    )
    set(WCNS_CHECKPOINT_GAMMA 1.4)
    set(WCNS_CHECKPOINT_OUTPUT restart-transport-invalid)
    set(WCNS_TRANSPORT_LINES
        "transport.model = constant\ntransport.constant.viscosity_ratio = 2.0")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/checkpoint_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/restart_transport_invalid.wcns
        @ONLY
    )
    set(WCNS_TRANSPORT_LINES "# default v1.0 transport migration")
    set(WCNS_CONTINUITY_OUTPUT restart-continuity-source)
    set(WCNS_CONTINUITY_CHECKPOINT true)
    set(WCNS_CONTINUITY_RESTART "# fresh start")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/restart_continuity_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/restart_continuity_source.wcns
        @ONLY
    )
    set(WCNS_CONTINUITY_CHECKPOINT false)
    set(WCNS_CONTINUITY_RESTART
        "restart.path = ${CMAKE_CURRENT_BINARY_DIR}/restart-continuity-source/restart-continuity.checkpoint.step00000001.time3p000000000eM04.cgns")
    set(WCNS_CONTINUITY_OUTPUT restart-continuity-serial)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/restart_continuity_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/restart_continuity_serial.wcns
        @ONLY
    )
    set(WCNS_CONTINUITY_OUTPUT restart-continuity-mpi)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/restart_continuity_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/restart_continuity_mpi.wcns
        @ONLY
    )
    set(WCNS_SA_CONTINUITY_OUTPUT sa-restart-continuity-source)
    set(WCNS_SA_CONTINUITY_CHECKPOINT true)
    set(WCNS_SA_CONTINUITY_RESTART "# fresh start")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/sa_restart_continuity_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/sa_restart_continuity_source.wcns
        @ONLY
    )
    set(WCNS_SA_CONTINUITY_OUTPUT sa-restart-continuity-serial)
    set(WCNS_SA_CONTINUITY_CHECKPOINT false)
    set(WCNS_SA_CONTINUITY_RESTART
        "restart.path = ${CMAKE_CURRENT_BINARY_DIR}/sa-restart-continuity-source/sa-restart-continuity.checkpoint.step00000001.time1p000000000eM04.cgns")
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/sa_restart_continuity_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/sa_restart_continuity_serial.wcns
        @ONLY
    )
    set(WCNS_SA_CONTINUITY_OUTPUT sa-restart-continuity-mpi)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/data/sa_restart_continuity_case.wcns.in
        ${CMAKE_CURRENT_BINARY_DIR}/sa_restart_continuity_mpi.wcns
        @ONLY
    )
    add_test(
        NAME wcns.run.freestream.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/freestream_case.wcns
    )
    foreach(implicit_case IN ITEMS steady dual_time)
        add_test(
            NAME wcns.run.implicit_freestream.${implicit_case}.serial
            COMMAND wcns_run
                --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_freestream_${implicit_case}.wcns
        )
        set_tests_properties(
            wcns.run.implicit_freestream.${implicit_case}.serial
            PROPERTIES FIXTURES_REQUIRED cgns_files
        )
    endforeach()
    foreach(les_model IN ITEMS
            smagorinsky
            scale_similarity
            mixed_smagorinsky_similarity
            dynamic_smagorinsky
            wale)
        add_test(
            NAME wcns.run.les_${les_model}.serial
            COMMAND wcns_run
                --config ${CMAKE_CURRENT_BINARY_DIR}/les_${les_model}_serial.wcns
        )
        set_tests_properties(
            wcns.run.les_${les_model}.serial
            PROPERTIES
                FIXTURES_REQUIRED cgns_files
                FIXTURES_SETUP les_${les_model}_serial
        )
        add_test(
            NAME wcns.run.les_implicit_${les_model}.serial
            COMMAND wcns_run
                --config ${CMAKE_CURRENT_BINARY_DIR}/les_implicit_${les_model}_serial.wcns
        )
        set_tests_properties(
            wcns.run.les_implicit_${les_model}.serial
            PROPERTIES FIXTURES_REQUIRED cgns_files
        )
    endforeach()
    add_test(
        NAME wcns.run.les_zero.none.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/les_zero_none.wcns
    )
    set_tests_properties(
        wcns.run.les_zero.none.serial
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP les_zero_none
    )
    add_test(
        NAME wcns.run.les_zero.smagorinsky.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/les_zero_smagorinsky.wcns
    )
    set_tests_properties(
        wcns.run.les_zero.smagorinsky.serial
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP les_zero_smagorinsky
    )
    add_test(
        NAME wcns.check.les_zero_equivalence.serial
        COMMAND wcns_cgns_field_compare
            ${CMAKE_CURRENT_BINARY_DIR}/les-zero-none/les-zero-none.field.step00000001.time1p000000000eM05.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/les-zero-smagorinsky/les-zero-smagorinsky.field.step00000001.time1p000000000eM05.cgns
            0
    )
    set_tests_properties(
        wcns.check.les_zero_equivalence.serial
        PROPERTIES FIXTURES_REQUIRED "les_zero_none;les_zero_smagorinsky"
    )
    add_test(
        NAME wcns.run.les_statistics.source.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/les_statistics_source.wcns
    )
    set_tests_properties(
        wcns.run.les_statistics.source.serial
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP les_statistics_source
    )
    add_test(
        NAME wcns.run.les_statistics.continuous.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/les_statistics_continuous_serial.wcns
    )
    set_tests_properties(
        wcns.run.les_statistics.continuous.serial
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP les_statistics_continuous_serial
    )
    add_test(
        NAME wcns.run.les_statistics.restart.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/les_statistics_restart_serial.wcns
    )
    set_tests_properties(
        wcns.run.les_statistics.restart.serial
        PROPERTIES
            FIXTURES_REQUIRED les_statistics_source
            FIXTURES_SETUP les_statistics_restart_serial
    )
    add_test(
        NAME wcns.check.les_statistics.restart.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_les_time_statistics.py
            ${CMAKE_CURRENT_BINARY_DIR}/les-statistics-continuous-serial/les-statistics-continuous-serial.time_statistics.r1.txt
            ${CMAKE_CURRENT_BINARY_DIR}/les-statistics-restart-serial/les-statistics-restart-serial.time_statistics.r1.txt
    )
    set_tests_properties(
        wcns.check.les_statistics.restart.serial
        PROPERTIES FIXTURES_REQUIRED "les_statistics_continuous_serial;les_statistics_restart_serial"
    )
    add_test(
        NAME wcns.run.implicit_sa.steady.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_sa_steady.wcns
    )
    set_tests_properties(
        wcns.run.implicit_sa.steady.serial
        PROPERTIES FIXTURES_REQUIRED cgns_files
    )
    foreach(two_equation_case IN ITEMS sst k_epsilon)
        add_test(
            NAME wcns.run.explicit_${two_equation_case}.steady.serial
            COMMAND wcns_run
                --config ${CMAKE_CURRENT_BINARY_DIR}/explicit_${two_equation_case}_steady.wcns
        )
        set_tests_properties(
            wcns.run.explicit_${two_equation_case}.steady.serial
            PROPERTIES FIXTURES_REQUIRED cgns_files
        )
        add_test(
            NAME wcns.run.implicit_${two_equation_case}.steady.serial
            COMMAND wcns_run
                --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_${two_equation_case}_steady.wcns
        )
        set_tests_properties(
            wcns.run.implicit_${two_equation_case}.steady.serial
            PROPERTIES
                FIXTURES_REQUIRED cgns_files
                FIXTURES_SETUP implicit_${two_equation_case}_steady_files
        )
        if(two_equation_case STREQUAL "sst")
            set(WCNS_TWO_EQUATION_SECOND_FIELD omega)
            set(WCNS_TWO_EQUATION_STEADY_STEM implicit-sst-steady)
        else()
            set(WCNS_TWO_EQUATION_SECOND_FIELD epsilon)
            set(WCNS_TWO_EQUATION_STEADY_STEM implicit-k-epsilon-steady)
        endif()
        add_test(
            NAME wcns.${two_equation_case}_history_floor_diagnostics.serial
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/check_two_equation_history.py
                ${CMAKE_CURRENT_BINARY_DIR}/${WCNS_TWO_EQUATION_STEADY_STEM}/${WCNS_TWO_EQUATION_STEADY_STEM}.history.r1.txt
                ${CMAKE_CURRENT_BINARY_DIR}/${WCNS_TWO_EQUATION_STEADY_STEM}/${WCNS_TWO_EQUATION_STEADY_STEM}.manifest.r1.txt
                k
                ${WCNS_TWO_EQUATION_SECOND_FIELD}
        )
        set_tests_properties(
            wcns.${two_equation_case}_history_floor_diagnostics.serial
            PROPERTIES FIXTURES_REQUIRED implicit_${two_equation_case}_steady_files
        )
        add_test(
            NAME wcns.run.implicit_${two_equation_case}.checkpoint_source.serial
            COMMAND wcns_run
                --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_${two_equation_case}_checkpoint_source.wcns
        )
        set_tests_properties(
            wcns.run.implicit_${two_equation_case}.checkpoint_source.serial
            PROPERTIES
                FIXTURES_REQUIRED cgns_files
                FIXTURES_SETUP implicit_${two_equation_case}_checkpoint
                RESOURCE_LOCK wcns_implicit_${two_equation_case}_checkpoint
        )
        add_test(
            NAME wcns.run.implicit_${two_equation_case}.checkpoint_restart.serial
            COMMAND wcns_run
                --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_${two_equation_case}_checkpoint_restart.wcns
        )
        set_tests_properties(
            wcns.run.implicit_${two_equation_case}.checkpoint_restart.serial
            PROPERTIES
                FIXTURES_REQUIRED implicit_${two_equation_case}_checkpoint
                RESOURCE_LOCK wcns_implicit_${two_equation_case}_checkpoint
        )
    endforeach()
    add_test(
        NAME wcns.run.implicit_wrong_model_checkpoint_restart.serial
        COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:wcns_run>
            -DCONFIG=${CMAKE_CURRENT_BINARY_DIR}/implicit_wrong_model_checkpoint_restart.wcns
            -DEXPECTED_EXIT=1
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_expect_exit.cmake
    )
    set_tests_properties(
        wcns.run.implicit_wrong_model_checkpoint_restart.serial
        PROPERTIES FIXTURES_REQUIRED implicit_sst_checkpoint
    )
    foreach(two_equation_case IN ITEMS sst k_epsilon)
        add_test(
            NAME wcns.run.implicit_${two_equation_case}.rank_consistency.serial
            COMMAND wcns_run
                --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_${two_equation_case}_serial.wcns
        )
        set_tests_properties(
            wcns.run.implicit_${two_equation_case}.rank_consistency.serial
            PROPERTIES
                FIXTURES_REQUIRED cgns_files
                FIXTURES_SETUP two_equation_${two_equation_case}_rank_serial
        )
    endforeach()
    add_test(
        NAME wcns.run.implicit_checkpoint.source.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_checkpoint_source.wcns
    )
    set_tests_properties(
        wcns.run.implicit_checkpoint.source.serial
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP implicit_checkpoint
    )
    add_test(
        NAME wcns.run.implicit_checkpoint.restart.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_checkpoint_restart.wcns
    )
    set_tests_properties(
        wcns.run.implicit_checkpoint.restart.serial
        PROPERTIES FIXTURES_REQUIRED implicit_checkpoint
    )
    add_test(
        NAME wcns.run.transport.sutherland.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_transport_runtime.py
            --run $<TARGET_FILE:wcns_run>
            --config ${CMAKE_CURRENT_BINARY_DIR}/transport_sutherland_serial.wcns
            --output ${CMAKE_CURRENT_BINARY_DIR}/transport-sutherland-serial
    )
    set_tests_properties(
        wcns.run.transport.sutherland.serial
        PROPERTIES FIXTURES_REQUIRED cgns_files
    )
    add_test(
        NAME wcns.run.double_mach_reflection.smoke.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/double_mach_reflection_case.wcns
    )
    set_tests_properties(
        wcns.run.double_mach_reflection.smoke.serial
        PROPERTIES FIXTURES_REQUIRED cgns_files
    )
    add_test(
        NAME wcns.generate_xz_plane_statistics_mesh
        COMMAND wcns_generate_release_cgns
            ${CMAKE_CURRENT_BINARY_DIR}/xz_plane_statistics_mesh.cgns
            3 8 8 8 1 0 false
    )
    set_tests_properties(
        wcns.generate_xz_plane_statistics_mesh
        PROPERTIES FIXTURES_SETUP xz_plane_statistics_mesh
    )
    add_test(
        NAME wcns.generate_cylinder_o_mesh
        COMMAND wcns_generate_release_cgns cylinder-o
            ${CMAKE_CURRENT_BINARY_DIR}/cylinder_o_smoke.cgns
            32 20 4 1.0 8.0 2.5
    )
    set_tests_properties(
        wcns.generate_cylinder_o_mesh
        PROPERTIES FIXTURES_SETUP cylinder_o_mesh
    )
    add_test(
        NAME wcns.run.cylinder_o.dry_run.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/cylinder_o_case.wcns
            --dry-run
    )
    set_tests_properties(
        wcns.run.cylinder_o.dry_run.serial
        PROPERTIES FIXTURES_REQUIRED cylinder_o_mesh
    )
    add_test(
        NAME wcns.run.boundary_partition.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/boundary_partition_case.wcns
    )
    set_tests_properties(
        wcns.run.boundary_partition.serial
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP boundary_partition_1
            RESOURCE_LOCK boundary_partition_output
    )
    add_test(
        NAME wcns.run.cylinder_o.boundary.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/cylinder_o_case.wcns
    )
    set_tests_properties(
        wcns.run.cylinder_o.boundary.serial
        PROPERTIES
            FIXTURES_REQUIRED cylinder_o_mesh
            FIXTURES_SETUP cylinder_boundary_output
    )
    add_test(
        NAME wcns.check.cylinder_o.boundary.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_boundary_output.py
            ${CMAKE_CURRENT_BINARY_DIR}/cylinder-o-grid-smoke
    )
    set_tests_properties(
        wcns.check.cylinder_o.boundary.serial
        PROPERTIES FIXTURES_REQUIRED cylinder_boundary_output
    )
    add_test(
        NAME wcns.run.case07.boundary.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/case07_boundary_case.wcns
    )
    set_tests_properties(
        wcns.run.case07.boundary.serial
        PROPERTIES FIXTURES_SETUP case07_boundary_output
    )
    add_test(
        NAME wcns.check.case07.boundary.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_boundary_output.py
            ${CMAKE_CURRENT_BINARY_DIR}/case07-boundary-smoke
    )
    set_tests_properties(
        wcns.check.case07.boundary.serial
        PROPERTIES FIXTURES_REQUIRED case07_boundary_output
    )
    add_test(
        NAME wcns.run.xz_plane_statistics.smoke.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/xz_plane_statistics_case.wcns
    )
    set_tests_properties(
        wcns.run.xz_plane_statistics.smoke.serial
        PROPERTIES
            FIXTURES_REQUIRED xz_plane_statistics_mesh
            FIXTURES_SETUP xz_plane_statistics_output
    )
    add_test(
        NAME wcns.check.xz_plane_statistics.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_xz_statistics.py
            ${CMAKE_CURRENT_BINARY_DIR}/xz-plane-statistics-smoke/xz-plane-statistics-smoke.statistics.r1.txt
    )
    set_tests_properties(
        wcns.check.xz_plane_statistics.serial
        PROPERTIES FIXTURES_REQUIRED xz_plane_statistics_output
    )
    add_test(
        NAME wcns.check.boundary_3d.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_boundary_output.py
            ${CMAKE_CURRENT_BINARY_DIR}/xz-plane-statistics-smoke
            1
    )
    set_tests_properties(
        wcns.check.boundary_3d.serial
        PROPERTIES FIXTURES_REQUIRED xz_plane_statistics_output
    )
    add_test(
        NAME wcns.run.les_boundary.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/les_boundary_serial.wcns
    )
    set_tests_properties(
        wcns.run.les_boundary.serial
        PROPERTIES
            FIXTURES_REQUIRED xz_plane_statistics_mesh
            FIXTURES_SETUP les_boundary_serial
            RESOURCE_LOCK les_boundary_output
    )
    add_test(
        NAME wcns.check.les_boundary.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_les_boundary.py
            ${CMAKE_CURRENT_BINARY_DIR}/les-boundary-equivalence
    )
    set_tests_properties(
        wcns.check.les_boundary.serial
        PROPERTIES FIXTURES_REQUIRED les_boundary_serial
    )
    set_tests_properties(
        wcns.run.freestream.serial
        PROPERTIES FIXTURES_REQUIRED cgns_files
    )
    foreach(periodic_case IN ITEMS translation_2d rotation_3d)
        add_test(
            NAME wcns.run.periodic_${periodic_case}.serial
            COMMAND wcns_run
                --config ${CMAKE_CURRENT_BINARY_DIR}/periodic_${periodic_case}.wcns
        )
        set_tests_properties(
            wcns.run.periodic_${periodic_case}.serial
            PROPERTIES FIXTURES_REQUIRED cgns_files
        )
    endforeach()
    foreach(viscous_case IN ITEMS couette_analytic linear_conduction_analytic)
        add_test(
            NAME wcns.run.${viscous_case}.serial
            COMMAND wcns_run
                --config ${CMAKE_CURRENT_BINARY_DIR}/${viscous_case}.wcns
        )
        set_tests_properties(
            wcns.run.${viscous_case}.serial
            PROPERTIES
                FIXTURES_REQUIRED cgns_files
                FIXTURES_SETUP ${viscous_case}_boundary_serial
        )
    endforeach()
    add_test(
        NAME wcns.check.couette_boundary.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_viscous_boundary_output.py
            couette
            ${CMAKE_CURRENT_BINARY_DIR}/couette-analytic
            1
    )
    set_tests_properties(
        wcns.check.couette_boundary.serial
        PROPERTIES FIXTURES_REQUIRED couette_analytic_boundary_serial
    )
    add_test(
        NAME wcns.check.conduction_boundary.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_viscous_boundary_output.py
            conduction
            ${CMAKE_CURRENT_BINARY_DIR}/linear-conduction-analytic
            1
    )
    set_tests_properties(
        wcns.check.conduction_boundary.serial
        PROPERTIES FIXTURES_REQUIRED linear_conduction_analytic_boundary_serial
    )
    add_test(
        NAME wcns.release_matrix.smoke.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-matrix-smoke-serial
            --ranks 1
    )
    add_test(
        NAME wcns.release_matrix.smoke.periodic
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-matrix-smoke-periodic
            --ranks 1
            --cells-i 16
            --cells-j 8
            --zones-i 2
            --periodic-x
    )
    add_test(
        NAME wcns.release_matrix.smoke.warped3d
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-matrix-smoke-warped3d
            --ranks 1
            --dimension 3
            --cells-i 16
            --cells-j 8
            --cells-k 8
            --warp 0.05
    )
    add_test(
        NAME wcns.release.vortex.smoke.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_vortex_matrix.py
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/isentropic_vortex.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-vortex-smoke-serial
            --ranks 1
            --resolutions 16
            --end-time 0.01
            --cfl 0.1
            --finest-l1 1e-3
    )
    add_test(
        NAME wcns.release.sod.smoke.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_shock_matrix.py
            --case sod
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/sod.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-sod-smoke-serial
            --ranks 1
            --cells-i 80
            --cells-j 8
            --end-time 0.2
    )
    add_test(
        NAME wcns.release.quadrant.smoke.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_shock_matrix.py
            --case quadrant
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/quadrant_riemann.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-quadrant-smoke-serial
            --ranks 1
            --cells-i 64
            --cells-j 64
            --end-time 0.05
    )
    add_test(
        NAME wcns.release.viscous.conduction.smoke.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_viscous_matrix.py
            --case conduction
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/viscous_channel.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-viscous-conduction-smoke-serial
            --ranks 1
            --cells-i 8
            --cells-j 12
            --min-cells 4
    )
    add_test(
        NAME wcns.release.source.uniform.smoke.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_source_matrix.py
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/uniform_source.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-source-uniform-smoke-serial
            --ranks 1
            --resolution 12
            --end-time 0.02
    )
    add_test(
        NAME wcns.release.source.manufactured_2d.smoke.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_manufactured_matrix.py
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/manufactured_source.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-source-manufactured-2d-smoke-serial
            --ranks 1
            --dimension 2
            --cells-i 12
            --cells-j 12
            --end-time 0.02
    )
    add_test(
        NAME wcns.release.output_restart.smoke.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_output_restart_matrix.py
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/output_restart.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-output-restart-smoke-serial
            --restart-ranks 1
            --resolution 12
    )
    add_test(
        NAME wcns.release.failures.smoke.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_failure_matrix.py
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/failure_base.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-failures-smoke-serial
            --ranks 1
    )
    set_tests_properties(
        wcns.release.failures.smoke.serial
        PROPERTIES TIMEOUT 120
    )
    foreach(profile IN ITEMS phenglei_wcns scmm6_wcns)
        add_test(
            NAME wcns.release.freestream.${profile}.2d.serial
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-o2a-${profile}-2d-serial
                --ranks 1
                --cells-i 32
                --cells-j 16
                --steps 10
                --profile ${profile}
                --case-prefix o2a-${profile}-2d
        )
        add_test(
            NAME wcns.release.freestream.${profile}.3d.serial
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-o2a-${profile}-3d-serial
                --ranks 1
                --dimension 3
                --cells-i 16
                --cells-j 16
                --cells-k 8
                --warp 0.05
                --steps 10
                --profile ${profile}
                --case-prefix o2a-${profile}-3d
        )
    endforeach()
    set_tests_properties(
        wcns.release.freestream.phenglei_wcns.2d.serial
        PROPERTIES FIXTURES_SETUP release_o2a_single_zone
    )
    add_test(
        NAME wcns.release.freestream.native_multizone.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
            --run $<TARGET_FILE:wcns_run>
            --generator $<TARGET_FILE:wcns_generate_release_cgns>
            --validator $<TARGET_FILE:wcns_validate_release_case>
            --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
            --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-o2a-native-multizone-serial
            --reference-directory ${CMAKE_CURRENT_BINARY_DIR}/release-o2a-phenglei_wcns-2d-serial/output-r1
            --ranks 1
            --cells-i 32
            --cells-j 16
            --zones-i 2
            --steps 10
            --case-prefix o2a-native-multizone
    )
    set_tests_properties(
        wcns.release.freestream.native_multizone.serial
        PROPERTIES FIXTURES_REQUIRED release_o2a_single_zone
    )
    add_test(
        NAME wcns.run.output.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/freestream_output_serial.wcns
    )
    set_tests_properties(
        wcns.run.output.serial
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP output_serial_files
            RESOURCE_LOCK wcns_output_validation
    )

    add_executable(wcns_field_output_tests test_field_output.cpp)
    target_link_libraries(wcns_field_output_tests PRIVATE wcns::core wcns::cgns)
    add_test(
        NAME wcns.output_files.serial
        COMMAND wcns_field_output_tests
            ${CMAKE_CURRENT_BINARY_DIR}/output-validation-serial
            1
            txt
    )
    set_tests_properties(
        wcns.output_files.serial
        PROPERTIES FIXTURES_REQUIRED output_serial_files
    )
    add_executable(wcns_checkpoint_file_tests test_checkpoint_file.cpp)
    target_link_libraries(wcns_checkpoint_file_tests PRIVATE wcns::core wcns::cgns)
    add_executable(wcns_cgns_field_compare test_cgns_field_compare.cpp)
    target_link_libraries(wcns_cgns_field_compare PRIVATE wcns::core wcns::cgns)
    add_test(
        NAME wcns.run.checkpoint_source
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/checkpoint_source.wcns
    )
    set_tests_properties(
        wcns.run.checkpoint_source
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP checkpoint_file
            RESOURCE_LOCK wcns_checkpoint_validation
    )
    add_test(
        NAME wcns.checkpoint_file
        COMMAND wcns_checkpoint_file_tests
            ${CMAKE_CURRENT_BINARY_DIR}/checkpoint-validation/checkpoint-freestream.checkpoint.latest.cgns
            --strip-default-transport
    )
    set_tests_properties(
        wcns.checkpoint_file
        PROPERTIES
            FIXTURES_REQUIRED checkpoint_file
            FIXTURES_SETUP checkpoint_legacy_transport
    )
    add_test(
        NAME wcns.run.restart.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/restart_serial.wcns
    )
    set_tests_properties(
        wcns.run.restart.serial
        PROPERTIES
            FIXTURES_REQUIRED checkpoint_legacy_transport
            RESOURCE_LOCK wcns_checkpoint_validation
    )
    add_test(
        NAME wcns.run.unsteady_exact_time
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/unsteady_case.wcns
    )
    set_tests_properties(
        wcns.run.unsteady_exact_time
        PROPERTIES FIXTURES_REQUIRED cgns_files
    )
    add_test(
        NAME wcns.run.maximum_steps_exit
        COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:wcns_run>
            -DCONFIG=${CMAKE_CURRENT_BINARY_DIR}/maximum_steps_case.wcns
            -DEXPECTED_EXIT=2
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_expect_exit.cmake
    )
    set_tests_properties(
        wcns.run.maximum_steps_exit
        PROPERTIES FIXTURES_REQUIRED cgns_files
    )
    add_test(
        NAME wcns.run.incompatible_restart
        COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:wcns_run>
            -DCONFIG=${CMAKE_CURRENT_BINARY_DIR}/restart_invalid.wcns
            -DEXPECTED_EXIT=1
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_expect_exit.cmake
    )
    set_tests_properties(
        wcns.run.incompatible_restart
        PROPERTIES FIXTURES_REQUIRED checkpoint_legacy_transport
    )
    add_test(
        NAME wcns.run.incompatible_transport_restart
        COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:wcns_run>
            -DCONFIG=${CMAKE_CURRENT_BINARY_DIR}/restart_transport_invalid.wcns
            -DEXPECTED_EXIT=1
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_expect_exit.cmake
    )
    set_tests_properties(
        wcns.run.incompatible_transport_restart
        PROPERTIES FIXTURES_REQUIRED checkpoint_legacy_transport
    )
    add_test(
        NAME wcns.run.existing_output_rejection
        COMMAND ${CMAKE_COMMAND}
            -DPROGRAM=$<TARGET_FILE:wcns_run>
            -DCONFIG=${CMAKE_CURRENT_BINARY_DIR}/existing_output_case.wcns
            -DOUTPUT_DIRECTORY=${CMAKE_CURRENT_BINARY_DIR}/existing-output-rejection
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_existing_output_rejection.cmake
    )
    set_tests_properties(
        wcns.run.existing_output_rejection
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            RESOURCE_LOCK wcns_existing_output_rejection
    )
    add_test(
        NAME wcns.run.restart_continuity_source
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/restart_continuity_source.wcns
    )
    set_tests_properties(
        wcns.run.restart_continuity_source
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP restart_continuity_source_files
            RESOURCE_LOCK wcns_restart_continuity
    )
    add_test(
        NAME wcns.run.restart_continuity.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/restart_continuity_serial.wcns
    )
    set_tests_properties(
        wcns.run.restart_continuity.serial
        PROPERTIES
            FIXTURES_REQUIRED restart_continuity_source_files
            FIXTURES_SETUP restart_continuity_serial_files
            RESOURCE_LOCK wcns_restart_continuity
    )
    set(WCNS_CONTINUITY_FINAL
        restart-continuity.field.step00000002.time6p000000000eM04.cgns)
    add_test(
        NAME wcns.restart_continuity.serial
        COMMAND wcns_cgns_field_compare
            ${CMAKE_CURRENT_BINARY_DIR}/restart-continuity-source/${WCNS_CONTINUITY_FINAL}
            ${CMAKE_CURRENT_BINARY_DIR}/restart-continuity-serial/${WCNS_CONTINUITY_FINAL}
            1e-13
    )
    set_tests_properties(
        wcns.restart_continuity.serial
        PROPERTIES FIXTURES_REQUIRED restart_continuity_serial_files
    )
    add_test(
        NAME wcns.run.sa_restart_continuity_source
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/sa_restart_continuity_source.wcns
    )
    set_tests_properties(
        wcns.run.sa_restart_continuity_source
        PROPERTIES
            FIXTURES_REQUIRED cgns_files
            FIXTURES_SETUP sa_restart_continuity_source_files
            RESOURCE_LOCK wcns_sa_restart_continuity
    )
    add_test(
        NAME wcns.sa_history_schema.unsteady.serial
        COMMAND ${Python3_EXECUTABLE}
            ${CMAKE_CURRENT_SOURCE_DIR}/check_sa_history.py
            ${CMAKE_CURRENT_BINARY_DIR}/sa-restart-continuity-source/sa-restart-continuity.history.r1.txt
    )
    set_tests_properties(
        wcns.sa_history_schema.unsteady.serial
        PROPERTIES FIXTURES_REQUIRED sa_restart_continuity_source_files
    )
    add_test(
        NAME wcns.run.sa_restart_continuity.serial
        COMMAND wcns_run
            --config ${CMAKE_CURRENT_BINARY_DIR}/sa_restart_continuity_serial.wcns
    )
    set_tests_properties(
        wcns.run.sa_restart_continuity.serial
        PROPERTIES
            FIXTURES_REQUIRED sa_restart_continuity_source_files
            FIXTURES_SETUP sa_restart_continuity_serial_files
            RESOURCE_LOCK wcns_sa_restart_continuity
    )
    set(WCNS_SA_CONTINUITY_FINAL
        sa-restart-continuity.field.step00000002.time2p000000000eM04.cgns)
    add_test(
        NAME wcns.sa_restart_continuity.serial
        COMMAND wcns_cgns_field_compare
            ${CMAKE_CURRENT_BINARY_DIR}/sa-restart-continuity-source/${WCNS_SA_CONTINUITY_FINAL}
            ${CMAKE_CURRENT_BINARY_DIR}/sa-restart-continuity-serial/${WCNS_SA_CONTINUITY_FINAL}
            1e-13
            nu_tilde
    )
    set_tests_properties(
        wcns.sa_restart_continuity.serial
        PROPERTIES FIXTURES_REQUIRED sa_restart_continuity_serial_files
    )

    add_executable(wcns_cgns_reader_tests
        test_cgns_reader.cpp
    )
    target_link_libraries(
        wcns_cgns_reader_tests
        PRIVATE wcns::io wcns::mesh wcns::parallel
    )
    add_test(
        NAME wcns.cgns_reader
        COMMAND wcns_cgns_reader_tests
            ${CMAKE_CURRENT_BINARY_DIR}/generated_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_3d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_invalid_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_3d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_one_sided_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_unknown_donor_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_periodic_2d.cgns
            ${CMAKE_CURRENT_BINARY_DIR}/generated_periodic_3d.cgns
    )
    set_tests_properties(wcns.cgns_reader PROPERTIES FIXTURES_REQUIRED cgns_files)

    add_executable(wcns_halo_exchange_tests test_halo_exchange.cpp)
    target_link_libraries(
        wcns_halo_exchange_tests
        PRIVATE wcns::io wcns::parallel wcns::solver
    )
    add_test(
        NAME wcns.halo_exchange.serial
        COMMAND wcns_halo_exchange_tests
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
    )
    set_tests_properties(wcns.halo_exchange.serial PROPERTIES FIXTURES_REQUIRED cgns_files)

    add_executable(wcns_euler_solver_tests test_euler_solver.cpp)
    target_link_libraries(
        wcns_euler_solver_tests
        PRIVATE wcns::io wcns::parallel wcns::solver
    )
    add_test(
        NAME wcns.euler_solver.serial
        COMMAND wcns_euler_solver_tests
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
    )
    set_tests_properties(wcns.euler_solver.serial PROPERTIES FIXTURES_REQUIRED cgns_files)

    add_executable(wcns_geometry_plan_tests test_geometry_plan_mpi.cpp)
    target_link_libraries(
        wcns_geometry_plan_tests
        PRIVATE wcns::io wcns::parallel wcns::mesh wcns::solver
    )
    add_test(
        NAME wcns.geometry_plan.serial
        COMMAND wcns_geometry_plan_tests
            ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
    )
    set_tests_properties(wcns.geometry_plan.serial PROPERTIES FIXTURES_REQUIRED cgns_files)

    if(WCNS_ENABLE_MPI)
        add_test(
            NAME wcns.run.transport.sutherland.2
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/check_transport_runtime.py
                --run $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/transport_sutherland_mpi.wcns
                --output ${CMAKE_CURRENT_BINARY_DIR}/transport-sutherland-mpi
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 2
        )
        set_tests_properties(
            wcns.run.transport.sutherland.2
            PROPERTIES FIXTURES_REQUIRED cgns_files
        )
        add_test(
            NAME wcns.run.boundary_partition.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/boundary_partition_case.wcns
        )
        set_tests_properties(
            wcns.run.boundary_partition.2
            PROPERTIES
                FIXTURES_REQUIRED boundary_partition_1
                FIXTURES_SETUP boundary_partition_2
                RESOURCE_LOCK boundary_partition_output
        )
        add_test(
            NAME wcns.run.boundary_partition.4
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 4
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/boundary_partition_case.wcns
        )
        set_tests_properties(
            wcns.run.boundary_partition.4
            PROPERTIES
                FIXTURES_REQUIRED boundary_partition_2
                FIXTURES_SETUP boundary_partition_4
                RESOURCE_LOCK boundary_partition_output
        )
        add_test(
            NAME wcns.check.boundary_partition.mpi_equivalence
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/check_boundary_rank_equivalence.py
                ${CMAKE_CURRENT_BINARY_DIR}/boundary-partition-equivalence
        )
        set_tests_properties(
            wcns.check.boundary_partition.mpi_equivalence
            PROPERTIES FIXTURES_REQUIRED boundary_partition_4
        )
        add_test(
            NAME wcns.run.cylinder_o.boundary.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/cylinder_o_case.wcns
        )
        set_tests_properties(
            wcns.run.cylinder_o.boundary.2
            PROPERTIES
                FIXTURES_REQUIRED cylinder_boundary_output
                FIXTURES_SETUP cylinder_boundary_output_2
                RESOURCE_LOCK cylinder_boundary_output
        )
        add_test(
            NAME wcns.run.cylinder_o.boundary.4
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 4
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/cylinder_o_case.wcns
        )
        set_tests_properties(
            wcns.run.cylinder_o.boundary.4
            PROPERTIES
                FIXTURES_REQUIRED cylinder_boundary_output_2
                FIXTURES_SETUP cylinder_boundary_output_4
                RESOURCE_LOCK cylinder_boundary_output
        )
        add_test(
            NAME wcns.check.cylinder_o.boundary.mpi_equivalence
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/check_boundary_rank_equivalence.py
                ${CMAKE_CURRENT_BINARY_DIR}/cylinder-o-grid-smoke
        )
        set_tests_properties(
            wcns.check.cylinder_o.boundary.mpi_equivalence
            PROPERTIES FIXTURES_REQUIRED cylinder_boundary_output_4
        )
        add_test(
            NAME wcns.run.xz_plane_statistics.smoke.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/xz_plane_statistics_mpi_case.wcns
        )
        set_tests_properties(
            wcns.run.xz_plane_statistics.smoke.2
            PROPERTIES
                FIXTURES_REQUIRED xz_plane_statistics_mesh
                FIXTURES_SETUP xz_plane_statistics_output_mpi
        )
        add_test(
            NAME wcns.check.xz_plane_statistics.2
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/check_xz_statistics.py
                ${CMAKE_CURRENT_BINARY_DIR}/xz-plane-statistics-smoke-mpi/xz-plane-statistics-smoke-mpi.statistics.r2.txt
        )
        set_tests_properties(
            wcns.check.xz_plane_statistics.2
            PROPERTIES FIXTURES_REQUIRED xz_plane_statistics_output_mpi
        )
        add_test(
            NAME wcns.check.boundary_3d.2
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/check_boundary_output.py
                ${CMAKE_CURRENT_BINARY_DIR}/xz-plane-statistics-smoke-mpi
                2
        )
        set_tests_properties(
            wcns.check.boundary_3d.2
            PROPERTIES FIXTURES_REQUIRED xz_plane_statistics_output_mpi
        )
        add_test(
            NAME wcns.run.les_boundary.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/les_boundary_mpi2.wcns
        )
        set_tests_properties(
            wcns.run.les_boundary.2
            PROPERTIES
                FIXTURES_REQUIRED les_boundary_serial
                FIXTURES_SETUP les_boundary_2
                RESOURCE_LOCK les_boundary_output
        )
        add_test(
            NAME wcns.run.les_boundary.4
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 4
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/les_boundary_mpi4.wcns
        )
        set_tests_properties(
            wcns.run.les_boundary.4
            PROPERTIES
                FIXTURES_REQUIRED les_boundary_2
                FIXTURES_SETUP les_boundary_4
                RESOURCE_LOCK les_boundary_output
        )
        add_test(
            NAME wcns.check.les_boundary.mpi_equivalence
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/check_boundary_rank_equivalence.py
                ${CMAKE_CURRENT_BINARY_DIR}/les-boundary-equivalence
        )
        set_tests_properties(
            wcns.check.les_boundary.mpi_equivalence
            PROPERTIES FIXTURES_REQUIRED les_boundary_4
        )
        add_test(
            NAME wcns.run.freestream.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/freestream_case.wcns
        )
        set_tests_properties(
            wcns.run.freestream.2
            PROPERTIES FIXTURES_REQUIRED cgns_files
        )
        foreach(rank_count IN ITEMS 2 4)
            foreach(implicit_case IN ITEMS steady dual_time)
                add_test(
                    NAME wcns.run.implicit_freestream.${implicit_case}.${rank_count}
                    COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                        $<TARGET_FILE:wcns_run>
                        --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_freestream_${implicit_case}.wcns
                )
                set_tests_properties(
                    wcns.run.implicit_freestream.${implicit_case}.${rank_count}
                    PROPERTIES FIXTURES_REQUIRED cgns_files
                )
            endforeach()
            foreach(les_model IN ITEMS
                    smagorinsky
                    scale_similarity
                    mixed_smagorinsky_similarity
                    dynamic_smagorinsky
                    wale)
                add_test(
                    NAME wcns.run.les_${les_model}.${rank_count}
                    COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                        $<TARGET_FILE:wcns_run>
                        --config ${CMAKE_CURRENT_BINARY_DIR}/les_${les_model}_mpi${rank_count}.wcns
                )
                set_tests_properties(
                    wcns.run.les_${les_model}.${rank_count}
                    PROPERTIES
                        FIXTURES_REQUIRED cgns_files
                        FIXTURES_SETUP les_${les_model}_${rank_count}
                )
                add_test(
                    NAME wcns.check.les_${les_model}.${rank_count}
                    COMMAND wcns_cgns_field_compare
                        ${CMAKE_CURRENT_BINARY_DIR}/les-${les_model}-rserial/les-${les_model}.field.step00000001.time1p000000000eM03.cgns
                        ${CMAKE_CURRENT_BINARY_DIR}/les-${les_model}-rmpi${rank_count}/les-${les_model}.field.step00000001.time1p000000000eM03.cgns
                        2e-12
                        mu_sgs
                        mu_sgs_over_mu
                        sgs_energy_transfer
                        sgs_stress_xx
                        sgs_stress_yy
                        sgs_stress_zz
                        sgs_stress_xy
                        sgs_stress_xz
                        sgs_stress_yz
                        les_dynamic_coefficient
                        les_filter_width
                        les_grid_anisotropy
                )
                set_tests_properties(
                    wcns.check.les_${les_model}.${rank_count}
                    PROPERTIES
                        FIXTURES_REQUIRED "les_${les_model}_serial;les_${les_model}_${rank_count}"
                )
                add_test(
                    NAME wcns.run.les_implicit_${les_model}.${rank_count}
                    COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                        $<TARGET_FILE:wcns_run>
                        --config ${CMAKE_CURRENT_BINARY_DIR}/les_implicit_${les_model}_mpi${rank_count}.wcns
                )
                set_tests_properties(
                    wcns.run.les_implicit_${les_model}.${rank_count}
                    PROPERTIES FIXTURES_REQUIRED cgns_files
                )
            endforeach()
            add_test(
                NAME wcns.run.les_statistics.continuous.${rank_count}
                COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                    $<TARGET_FILE:wcns_run>
                    --config ${CMAKE_CURRENT_BINARY_DIR}/les_statistics_continuous_mpi${rank_count}.wcns
            )
            set_tests_properties(
                wcns.run.les_statistics.continuous.${rank_count}
                PROPERTIES
                    FIXTURES_REQUIRED cgns_files
                    FIXTURES_SETUP les_statistics_continuous_${rank_count}
            )
            add_test(
                NAME wcns.run.les_statistics.restart.${rank_count}
                COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                    $<TARGET_FILE:wcns_run>
                    --config ${CMAKE_CURRENT_BINARY_DIR}/les_statistics_restart_mpi${rank_count}.wcns
            )
            set_tests_properties(
                wcns.run.les_statistics.restart.${rank_count}
                PROPERTIES
                    FIXTURES_REQUIRED les_statistics_source
                    FIXTURES_SETUP les_statistics_restart_${rank_count}
            )
            add_test(
                NAME wcns.check.les_statistics.restart.${rank_count}
                COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/check_les_time_statistics.py
                    ${CMAKE_CURRENT_BINARY_DIR}/les-statistics-continuous-serial/les-statistics-continuous-serial.time_statistics.r1.txt
                    ${CMAKE_CURRENT_BINARY_DIR}/les-statistics-continuous-mpi${rank_count}/les-statistics-continuous-mpi${rank_count}.time_statistics.r${rank_count}.txt
                    ${CMAKE_CURRENT_BINARY_DIR}/les-statistics-restart-mpi${rank_count}/les-statistics-restart-mpi${rank_count}.time_statistics.r${rank_count}.txt
            )
            set_tests_properties(
                wcns.check.les_statistics.restart.${rank_count}
                PROPERTIES FIXTURES_REQUIRED "les_statistics_continuous_serial;les_statistics_continuous_${rank_count};les_statistics_restart_${rank_count}"
            )
            add_test(
                NAME wcns.run.implicit_sa.steady.${rank_count}
                COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                    $<TARGET_FILE:wcns_run>
                    --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_sa_steady.wcns
            )
            set_tests_properties(
                wcns.run.implicit_sa.steady.${rank_count}
                PROPERTIES FIXTURES_REQUIRED cgns_files
            )
            foreach(two_equation_case IN ITEMS sst k_epsilon)
                add_test(
                    NAME wcns.run.explicit_${two_equation_case}.steady.${rank_count}
                    COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                        $<TARGET_FILE:wcns_run>
                        --config ${CMAKE_CURRENT_BINARY_DIR}/explicit_${two_equation_case}_steady.wcns
                )
                set_tests_properties(
                    wcns.run.explicit_${two_equation_case}.steady.${rank_count}
                    PROPERTIES FIXTURES_REQUIRED cgns_files
                )
                add_test(
                    NAME wcns.run.implicit_${two_equation_case}.steady.${rank_count}
                    COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                        $<TARGET_FILE:wcns_run>
                        --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_${two_equation_case}_steady.wcns
                )
                set_tests_properties(
                    wcns.run.implicit_${two_equation_case}.steady.${rank_count}
                    PROPERTIES FIXTURES_REQUIRED cgns_files
                )
                add_test(
                    NAME wcns.run.implicit_${two_equation_case}.checkpoint_restart.${rank_count}
                    COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                        $<TARGET_FILE:wcns_run>
                        --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_${two_equation_case}_checkpoint_restart.wcns
                )
                set_tests_properties(
                    wcns.run.implicit_${two_equation_case}.checkpoint_restart.${rank_count}
                    PROPERTIES
                        FIXTURES_REQUIRED implicit_${two_equation_case}_checkpoint
                        RESOURCE_LOCK wcns_implicit_${two_equation_case}_checkpoint
                )
                add_test(
                    NAME wcns.run.implicit_${two_equation_case}.rank_consistency.${rank_count}
                    COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                        $<TARGET_FILE:wcns_run>
                        --config ${CMAKE_CURRENT_BINARY_DIR}/implicit_${two_equation_case}_mpi${rank_count}.wcns
                )
                set_tests_properties(
                    wcns.run.implicit_${two_equation_case}.rank_consistency.${rank_count}
                    PROPERTIES
                        FIXTURES_REQUIRED cgns_files
                        FIXTURES_SETUP two_equation_${two_equation_case}_rank_${rank_count}
                )
                if(two_equation_case STREQUAL "sst")
                    set(WCNS_TWO_EQUATION_COMPARISON_FIELD omega)
                else()
                    set(WCNS_TWO_EQUATION_COMPARISON_FIELD epsilon)
                endif()
                set(WCNS_TWO_EQUATION_RANK_FIELD
                    implicit-${two_equation_case}-rank-consistency.field.step00000002.time2p000000000eP00.cgns)
                add_test(
                    NAME wcns.implicit_${two_equation_case}.rank_consistency.${rank_count}
                    COMMAND wcns_cgns_field_compare
                        ${CMAKE_CURRENT_BINARY_DIR}/two-equation-${two_equation_case}-serial/${WCNS_TWO_EQUATION_RANK_FIELD}
                        ${CMAKE_CURRENT_BINARY_DIR}/two-equation-${two_equation_case}-mpi${rank_count}/${WCNS_TWO_EQUATION_RANK_FIELD}
                        1e-9
                        k
                        ${WCNS_TWO_EQUATION_COMPARISON_FIELD}
                )
                set_tests_properties(
                    wcns.implicit_${two_equation_case}.rank_consistency.${rank_count}
                    PROPERTIES FIXTURES_REQUIRED
                        "two_equation_${two_equation_case}_rank_serial;two_equation_${two_equation_case}_rank_${rank_count}"
                )
            endforeach()
        endforeach()
        foreach(periodic_case IN ITEMS translation_2d rotation_3d)
            add_test(
                NAME wcns.run.periodic_${periodic_case}.2
                COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                    $<TARGET_FILE:wcns_run>
                    --config ${CMAKE_CURRENT_BINARY_DIR}/periodic_${periodic_case}.wcns
            )
            set_tests_properties(
                wcns.run.periodic_${periodic_case}.2
                PROPERTIES FIXTURES_REQUIRED cgns_files
            )
        endforeach()
        add_test(
            NAME wcns.run.periodic_translation_2d.4
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 4
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/periodic_translation_2d.wcns
        )
        set_tests_properties(
            wcns.run.periodic_translation_2d.4
            PROPERTIES FIXTURES_REQUIRED cgns_files
        )
        foreach(viscous_case IN ITEMS couette_analytic linear_conduction_analytic)
            add_test(
                NAME wcns.run.${viscous_case}.2
                COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                    $<TARGET_FILE:wcns_run>
                    --config ${CMAKE_CURRENT_BINARY_DIR}/${viscous_case}.wcns
            )
            set_tests_properties(
                wcns.run.${viscous_case}.2
                PROPERTIES
                    FIXTURES_REQUIRED cgns_files
                    FIXTURES_SETUP ${viscous_case}_boundary_2
            )
        endforeach()
        add_test(
            NAME wcns.check.couette_boundary.2
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/check_viscous_boundary_output.py
                couette
                ${CMAKE_CURRENT_BINARY_DIR}/couette-analytic
                2
        )
        set_tests_properties(
            wcns.check.couette_boundary.2
            PROPERTIES FIXTURES_REQUIRED couette_analytic_boundary_2
        )
        add_test(
            NAME wcns.check.conduction_boundary.2
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/check_viscous_boundary_output.py
                conduction
                ${CMAKE_CURRENT_BINARY_DIR}/linear-conduction-analytic
                2
        )
        set_tests_properties(
            wcns.check.conduction_boundary.2
            PROPERTIES FIXTURES_REQUIRED linear_conduction_analytic_boundary_2
        )
        add_test(
            NAME wcns.release_matrix.smoke.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-matrix-smoke-mpi
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 1,2
        )
        add_test(
            NAME wcns.release.vortex.smoke.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_vortex_matrix.py
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/isentropic_vortex.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-vortex-smoke-mpi
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 1,2
                --resolutions 16
                --end-time 0.01
                --cfl 0.1
                --finest-l1 1e-3
        )
        add_test(
            NAME wcns.release.sod.smoke.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_shock_matrix.py
                --case sod
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/sod.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-sod-smoke-mpi
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 1,2
                --cells-i 80
                --cells-j 8
                --zones-i 2
                --end-time 0.2
        )
        add_test(
            NAME wcns.release.quadrant.smoke.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_shock_matrix.py
                --case quadrant
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/quadrant_riemann.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-quadrant-smoke-mpi
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 1,2
                --cells-i 64
                --cells-j 64
                --zones-i 2
                --end-time 0.05
        )
        add_test(
            NAME wcns.release.viscous.conduction.smoke.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_viscous_matrix.py
                --case conduction
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/viscous_channel.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-viscous-conduction-smoke-mpi
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 1,2
                --cells-i 8
                --cells-j 12
                --min-cells 4
        )
        add_test(
            NAME wcns.release.source.uniform.smoke.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_source_matrix.py
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/uniform_source.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-source-uniform-smoke-mpi
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 1,2
                --resolution 12
                --end-time 0.02
        )
        add_test(
            NAME wcns.release.source.manufactured_3d.smoke.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_manufactured_matrix.py
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/manufactured_source.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-source-manufactured-3d-smoke-mpi
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 1,2
                --dimension 3
                --cells-i 8
                --cells-j 8
                --cells-k 8
                --warp 0.05
                --end-time 0.02
        )
        add_test(
            NAME wcns.release.output_restart.smoke.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_output_restart_matrix.py
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/output_restart.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-output-restart-smoke-mpi
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --restart-ranks 1,2
                --resolution 12
        )
        add_test(
            NAME wcns.release.failures.smoke.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_failure_matrix.py
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/failure_base.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-failures-smoke-mpi
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 1,2
        )
        set_tests_properties(
            wcns.release.failures.smoke.mpi
            PROPERTIES TIMEOUT 180
        )
        foreach(profile IN ITEMS phenglei_wcns scmm6_wcns)
            add_test(
                NAME wcns.release.freestream.${profile}.2d.mpi
                COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
                    --run $<TARGET_FILE:wcns_run>
                    --generator $<TARGET_FILE:wcns_generate_release_cgns>
                    --validator $<TARGET_FILE:wcns_validate_release_case>
                    --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
                    --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-o2a-${profile}-2d-mpi
                    --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                    --ranks 1,2,4
                    --cells-i 32
                    --cells-j 16
                    --steps 10
                    --profile ${profile}
                    --case-prefix o2a-${profile}-2d-mpi
            )
            add_test(
                NAME wcns.release.freestream.${profile}.3d.mpi
                COMMAND ${Python3_EXECUTABLE}
                    ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
                    --run $<TARGET_FILE:wcns_run>
                    --generator $<TARGET_FILE:wcns_generate_release_cgns>
                    --validator $<TARGET_FILE:wcns_validate_release_case>
                    --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
                    --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-o2a-${profile}-3d-mpi
                    --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                    --ranks 1,2
                    --dimension 3
                    --cells-i 16
                    --cells-j 16
                    --cells-k 8
                    --warp 0.05
                    --steps 10
                    --profile ${profile}
                    --case-prefix o2a-${profile}-3d-mpi
            )
        endforeach()
        add_test(
            NAME wcns.release.freestream.native_multizone.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-o2a-native-multizone-mpi
                --reference-directory ${CMAKE_CURRENT_BINARY_DIR}/release-o2a-phenglei_wcns-2d-serial/output-r1
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 1,2
                --cells-i 32
                --cells-j 16
                --zones-i 2
                --steps 10
                --case-prefix o2a-native-multizone-mpi
        )
        set_tests_properties(
            wcns.release.freestream.native_multizone.mpi
            PROPERTIES FIXTURES_REQUIRED release_o2a_single_zone
        )
        add_test(
            NAME wcns.release.freestream.phenglei_wcns.3d_k_split.mpi
            COMMAND ${Python3_EXECUTABLE}
                ${CMAKE_CURRENT_SOURCE_DIR}/../tools/run_release_matrix.py
                --run $<TARGET_FILE:wcns_run>
                --generator $<TARGET_FILE:wcns_generate_release_cgns>
                --validator $<TARGET_FILE:wcns_validate_release_case>
                --template ${CMAKE_CURRENT_SOURCE_DIR}/../cases/config/freestream.wcns.in
                --work-dir ${CMAKE_CURRENT_BINARY_DIR}/release-o2a-phenglei-3d-k-split
                --mpi-exec ${WCNS_MPIEXEC_EXECUTABLE}
                --ranks 1,2
                --dimension 3
                --cells-i 8
                --cells-j 8
                --cells-k 16
                --warp 0.05
                --steps 10
                --profile phenglei_wcns
                --case-prefix o2a-phenglei-3d-k-split
        )
        add_test(
            NAME wcns.run.output.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/freestream_output_mpi.wcns
        )
        set_tests_properties(
            wcns.run.output.2
            PROPERTIES
                FIXTURES_REQUIRED cgns_files
                FIXTURES_SETUP output_mpi_files
                RESOURCE_LOCK wcns_output_validation
        )
        add_test(
            NAME wcns.output_files.2
        COMMAND wcns_field_output_tests
                ${CMAKE_CURRENT_BINARY_DIR}/output-validation-mpi
                2
                tecplot
        )
        set_tests_properties(
            wcns.output_files.2
            PROPERTIES FIXTURES_REQUIRED output_mpi_files
        )
        add_test(
            NAME wcns.run.restart.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/restart_mpi.wcns
        )
        set_tests_properties(
            wcns.run.restart.2
            PROPERTIES
                FIXTURES_REQUIRED checkpoint_legacy_transport
                RESOURCE_LOCK wcns_checkpoint_validation
        )
        add_test(
            NAME wcns.run.restart_continuity.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/restart_continuity_mpi.wcns
        )
        set_tests_properties(
            wcns.run.restart_continuity.2
            PROPERTIES
                FIXTURES_REQUIRED restart_continuity_source_files
                FIXTURES_SETUP restart_continuity_mpi_files
                RESOURCE_LOCK wcns_restart_continuity
        )
        add_test(
            NAME wcns.restart_continuity.2
            COMMAND wcns_cgns_field_compare
                ${CMAKE_CURRENT_BINARY_DIR}/restart-continuity-source/${WCNS_CONTINUITY_FINAL}
                ${CMAKE_CURRENT_BINARY_DIR}/restart-continuity-mpi/${WCNS_CONTINUITY_FINAL}
                1e-12
        )
        set_tests_properties(
            wcns.restart_continuity.2
            PROPERTIES FIXTURES_REQUIRED restart_continuity_mpi_files
        )
        add_test(
            NAME wcns.run.sa_restart_continuity.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_run>
                --config ${CMAKE_CURRENT_BINARY_DIR}/sa_restart_continuity_mpi.wcns
        )
        set_tests_properties(
            wcns.run.sa_restart_continuity.2
            PROPERTIES
                FIXTURES_REQUIRED sa_restart_continuity_source_files
                FIXTURES_SETUP sa_restart_continuity_mpi_files
                RESOURCE_LOCK wcns_sa_restart_continuity
        )
        add_test(
            NAME wcns.sa_restart_continuity.2
            COMMAND wcns_cgns_field_compare
                ${CMAKE_CURRENT_BINARY_DIR}/sa-restart-continuity-source/${WCNS_SA_CONTINUITY_FINAL}
                ${CMAKE_CURRENT_BINARY_DIR}/sa-restart-continuity-mpi/${WCNS_SA_CONTINUITY_FINAL}
                1e-12
                nu_tilde
        )
        set_tests_properties(
            wcns.sa_restart_continuity.2
            PROPERTIES FIXTURES_REQUIRED sa_restart_continuity_mpi_files
        )
        foreach(rank_count IN ITEMS 4 8)
            add_test(
                NAME wcns.run.freestream.${rank_count}
                COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n ${rank_count}
                    $<TARGET_FILE:wcns_run>
                    --config ${CMAKE_CURRENT_BINARY_DIR}/freestream_case.wcns
            )
            set_tests_properties(
                wcns.run.freestream.${rank_count}
                PROPERTIES FIXTURES_REQUIRED cgns_files
            )
        endforeach()
        add_test(
            NAME wcns.halo_exchange.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_halo_exchange_tests>
                ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
        )
        set_tests_properties(wcns.halo_exchange.2 PROPERTIES FIXTURES_REQUIRED cgns_files)
        add_test(
            NAME wcns.halo_exchange.4
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 4
                $<TARGET_FILE:wcns_halo_exchange_tests>
                ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
        )
        set_tests_properties(wcns.halo_exchange.4 PROPERTIES FIXTURES_REQUIRED cgns_files)

        add_test(
            NAME wcns.euler_solver.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_euler_solver_tests>
                ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
        )
        set_tests_properties(wcns.euler_solver.2 PROPERTIES FIXTURES_REQUIRED cgns_files)

        add_test(
            NAME wcns.geometry_plan.2
            COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
                $<TARGET_FILE:wcns_geometry_plan_tests>
                ${CMAKE_CURRENT_BINARY_DIR}/generated_multi_2d.cgns
        )
        set_tests_properties(wcns.geometry_plan.2 PROPERTIES FIXTURES_REQUIRED cgns_files)
    endif()
endif()

add_test(NAME wcns.unit COMMAND wcns_unit_tests)

add_executable(wcns_robustness_tests test_robustness_mpi.cpp)
target_link_libraries(
    wcns_robustness_tests
    PRIVATE wcns::parallel wcns::solver wcns::mesh
)
add_test(NAME wcns.robustness.serial COMMAND wcns_robustness_tests)

add_executable(wcns_mpi_runtime_tests test_mpi_runtime.cpp)
target_link_libraries(wcns_mpi_runtime_tests PRIVATE wcns::parallel wcns::solver)
add_test(NAME wcns.mpi_runtime.serial COMMAND wcns_mpi_runtime_tests)

add_executable(wcns_stage_j_solver_tests test_inviscid_wcns_solver_mpi.cpp)
target_link_libraries(
    wcns_stage_j_solver_tests
    PRIVATE wcns::parallel wcns::solver wcns::mesh
)
add_test(NAME wcns.stage_j_solver.serial COMMAND wcns_stage_j_solver_tests)

add_executable(wcns_stage_k_solver_tests test_viscous_wcns_solver_mpi.cpp)
target_link_libraries(
    wcns_stage_k_solver_tests
    PRIVATE wcns::parallel wcns::solver wcns::mesh
)
add_test(NAME wcns.stage_k_solver.serial COMMAND wcns_stage_k_solver_tests)

if(WCNS_ENABLE_MPI)
    if(NOT WCNS_MPIEXEC_EXECUTABLE)
        message(FATAL_ERROR "MPI tests require WCNS_MPIEXEC_EXECUTABLE")
    endif()
    add_test(
        NAME wcns.mpi_runtime.2
        COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2 $<TARGET_FILE:wcns_mpi_runtime_tests>
    )
    add_test(
        NAME wcns.mpi_runtime.4
        COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 4 $<TARGET_FILE:wcns_mpi_runtime_tests>
    )
    add_test(
        NAME wcns.stage_j_solver.2
        COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
            $<TARGET_FILE:wcns_stage_j_solver_tests>
    )
    add_test(
        NAME wcns.stage_k_solver.2
        COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
            $<TARGET_FILE:wcns_stage_k_solver_tests>
    )
    add_test(
        NAME wcns.robustness.2
        COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 2
            $<TARGET_FILE:wcns_robustness_tests>
    )
    add_test(
        NAME wcns.robustness.4
        COMMAND ${WCNS_MPIEXEC_EXECUTABLE} -n 4
            $<TARGET_FILE:wcns_robustness_tests>
    )
endif()
