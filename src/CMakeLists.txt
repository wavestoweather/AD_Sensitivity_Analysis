set(TYPES_SRCS
    microphysics/constants.cpp
    types/table_t.cpp
    microphysics/physical_parameterizations.cpp
    types/particle_model_constants_t.cpp
    types/gamma_table_t.cpp
    types/global_args_t.cpp
    types/input_parameters_t.cpp
    types/model_constants_t.cpp
    types/output_handle_t.cpp
    types/param_t.cpp
    types/segment_t.cpp
    types/checkpoint_t.cpp
    types/nc_parameters_t.cpp
    types/task_scheduler_t.cpp
    misc/general.cpp
    )

if( TARGET )
    if( TARGET MATCHES "simulation" )
        add_executable(trajectories microphysics/trajectories.cpp)
        target_link_libraries(trajectories types)
    endif()

    if( TARGET MATCHES "load_test" )
        add_executable(load_test scratch/load_test.cpp)
        target_link_libraries(load_test types)
    endif()

    if( TARGET MATCHES "netcdf_test" )
        add_executable(netcdf_test scratch/netcdf_test.cpp)
        target_link_libraries(netcdf_test types)
    endif()

    if( TARGET MATCHES "scan" )
        add_executable(scan scratch/scan.cpp)
        target_link_libraries(scan types)
    endif()
else()
    message("~~~ No target using -DTARGET specified. Using default 'simulation'")
    message("~~~ Possible targets: simulation load_test netcdf_test scan")
    add_executable(trajectories ${TYPES_SRCS} microphysics/trajectories.cpp )
    set_property(TARGET trajectories PROPERTY CXX_STANDARD 14)
    set_property(TARGET trajectories PROPERTY CXX_STANDARD_REQUIRED ON)

    if(MPI_CXX_LIBRARIES)
        target_link_libraries(trajectories ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES})
    else()
        target_link_libraries(trajectories ${GSL_LIBRARIES} ${NETCDF_LIBRARIES})
    endif()
endif()



