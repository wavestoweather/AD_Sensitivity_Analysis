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
set(REGRID_SRCS
    microphysics/constants.cpp
    types/netcdf_simulation_reader_t.cpp
    )

if( TARGET )
    if( "${TARGET}" MATCHES "simulation" )
        add_executable(trajectories ${TYPES_SRCS} microphysics/trajectories.cpp )
        set_property(TARGET trajectories PROPERTY CXX_STANDARD 17)
        set_property(TARGET trajectories PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(trajectories ${targets})
    endif()
    if( "${TARGET}" MATCHES "regrid" )
        add_executable(regrid ${REGRID_SRCS} misc/regrid.cpp)
        set_property(TARGET regrid PROPERTY CXX_STANDARD 17)
        set_property(TARGET regrid PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(regrid ${targets})
    endif()
    if( "${TARGET}" MATCHES "loom" )
        add_executable(loom ${REGRID_SRCS} misc/loom.cpp)
        set_property(TARGET loom PROPERTY CXX_STANDARD 17)
        set_property(TARGET loom PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(loom ${targets})
    endif()
    if( "${TARGET}" MATCHES "timing" )
        add_executable(timing ${TYPES_SRCS} misc/timing.cpp)
        set_property(TARGET timing PROPERTY CXX_STANDARD 17)
        set_property(TARGET timing PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(timing ${targets})
    endif()
    if( "${TARGET}" MATCHES "python_interface" )
        set(TYPES_SRCS
                microphysics/constants.cpp
                types/table_t.cpp
                microphysics/physical_parameterizations.cpp
                types/particle_model_constants_t.cpp
                types/gamma_table_t.cpp
                types/global_args_t.cpp
                types/input_parameters_t.cpp
                types/model_constants_t.cpp
                types/param_t.cpp
                types/segment_t.cpp
                types/checkpoint_t.cpp
                types/task_scheduler_t.cpp
                misc/general.cpp
                )
        set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}")
        add_library(python_interface MODULE ${TYPES_SRCS} interface/physics_t.cpp)
        set_property(TARGET python_interface PROPERTY CXX_STANDARD 17)
        set_property(TARGET python_interface PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(python_interface ${targets})
    endif()
else()
    message("${BoldYellow}~~~ No target using -DTARGET specified. Using default 'simulation'")
    message("~~~ Possible targets: simulation regrid loom timing python_interface${ColourReset}")
    add_executable(trajectories ${TYPES_SRCS} microphysics/trajectories.cpp )
    set_property(TARGET trajectories PROPERTY CXX_STANDARD 17)
    set_property(TARGET trajectories PROPERTY CXX_STANDARD_REQUIRED ON)
    target_link_libraries(trajectories ${targets})
endif()



