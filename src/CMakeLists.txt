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
    types/netcdf_reader_t.cpp
    types/task_scheduler_t.cpp
    misc/general.cpp
    )

if( TARGET )
    if( "${TARGET}" MATCHES "simulation" )
        add_executable(trajectories ${TYPES_SRCS} microphysics/trajectories.cpp )
        set_property(TARGET trajectories PROPERTY CXX_STANDARD 14)
        set_property(TARGET trajectories PROPERTY CXX_STANDARD_REQUIRED ON)

        if(MPI_CXX_LIBRARIES)
            target_link_libraries(trajectories ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES})
        else()
            target_link_libraries(trajectories ${GSL_LIBRARIES} ${NETCDF_LIBRARIES})
        endif()
    endif()
    if( "${TARGET}" MATCHES "timing" )
        add_executable(timing ${TYPES_SRCS} misc/timing.cpp)
        set_property(TARGET timing PROPERTY CXX_STANDARD 14)
        set_property(TARGET timing PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(timing ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES})
    endif()
    if( "${TARGET}" MATCHES "scan" )
        add_executable(scan ${TYPES_SRCS} misc/scan.cpp)
        set_property(TARGET scan PROPERTY CXX_STANDARD 14)
        set_property(TARGET scan PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(timing ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES})
    endif()
else()
    message("${BoldYellow}~~~ No target using -DTARGET specified. Using default 'simulation'")
    message("~~~ Possible targets: simulation timing scan${ColourReset}")
    add_executable(trajectories ${TYPES_SRCS} microphysics/trajectories.cpp )
    set_property(TARGET trajectories PROPERTY CXX_STANDARD 14)
    set_property(TARGET trajectories PROPERTY CXX_STANDARD_REQUIRED ON)

    if(MPI_CXX_LIBRARIES)
        target_link_libraries(trajectories ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES})
    else()
        target_link_libraries(trajectories ${GSL_LIBRARIES} ${NETCDF_LIBRARIES})
    endif()
endif()



