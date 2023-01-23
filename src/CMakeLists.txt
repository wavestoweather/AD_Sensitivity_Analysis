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
    if( "${TARGET}" MATCHES "timing" )
        add_executable(timing ${TYPES_SRCS} misc/timing.cpp)
        set_property(TARGET timing PROPERTY CXX_STANDARD 17)
        set_property(TARGET timing PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(timing ${targets})
    endif()
    if( "${TARGET}" MATCHES "scan" )
        add_executable(scan ${TYPES_SRCS} misc/scan.cpp)
        set_property(TARGET scan PROPERTY CXX_STANDARD 17)
        set_property(TARGET scan PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(scan ${targets})
    endif()
    if( "${TARGET}" MATCHES "scratch" )
        add_executable(test_tape ${TYPES_SRCS} scratch/test_tape.cpp)
        set_property(TARGET test_tape PROPERTY CXX_STANDARD 17)
        set_property(TARGET test_tape PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(test_tap ${targets})
    endif()
    if( "${TARGET}" MATCHES "reader_test" )
        add_executable(reader_test ${TYPES_SRCS} tests/reader_test.cpp)
        set_property(TARGET reader_test PROPERTY CXX_STANDARD 17)
        set_property(TARGET reader_test PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(reader_test ${targets})
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
    message("~~~ Possible targets: simulation timing scan${ColourReset}")
    add_executable(trajectories ${TYPES_SRCS} microphysics/trajectories.cpp )
    set_property(TARGET trajectories PROPERTY CXX_STANDARD 17)
    set_property(TARGET trajectories PROPERTY CXX_STANDARD_REQUIRED ON)
    target_link_libraries(trajectories ${targets})
endif()



