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
        set_property(TARGET trajectories PROPERTY CXX_STANDARD 17)
        set_property(TARGET trajectories PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(trajectories ${targets} nlohmann_json::nlohmann_json)
        # if(MPI_CXX_LIBRARIES)
        #     target_link_libraries(trajectories ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${CURL_LIBRARIES} ${HDF5_LIBRARIES} ${PNETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${Zip_LIBRARIES} ${AEC_LIBRARIES} ${CMAKE_DL_LIBS}) # CMAKE_DL_LIBS
        # else()
        #     target_link_libraries(trajectories ${GSL_LIBRARIES} ${CURL_LIBRARIES} ${HDF5_LIBRARIES} ${PNETCDF_LIBRARIES} ${NETCDF_LIBRARIES} ${ZLIB_LIBRARIES} ${Zip_LIBRARIES} ${AEC_LIBRARIES} ${CMAKE_DL_LIBS})
        # endif()
    endif()
    if( "${TARGET}" MATCHES "timing" )
        add_executable(timing ${TYPES_SRCS} misc/timing.cpp)
        set_property(TARGET timing PROPERTY CXX_STANDARD 17)
        set_property(TARGET timing PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(timing ${targets} nlohmann_json::nlohmann_json)
        # target_link_libraries(timing ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${CURL_LIBRARIES} ${HDF5_LIBRARIES} ${PNETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${Zip_LIBRARIES} ${AEC_LIBRARIES} ${CMAKE_DL_LIBS})
    endif()
    if( "${TARGET}" MATCHES "scan" )
        add_executable(scan ${TYPES_SRCS} misc/scan.cpp)
        set_property(TARGET scan PROPERTY CXX_STANDARD 17)
        set_property(TARGET scan PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(scan ${targets} nlohmann_json::nlohmann_json)
        # target_link_libraries(timing ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${CURL_LIBRARIES} ${HDF5_LIBRARIES} ${PNETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${Zip_LIBRARIES} ${AEC_LIBRARIES} ${CMAKE_DL_LIBS})
    endif()
    if( "${TARGET}" MATCHES "scratch" )
        add_executable(test_tape ${TYPES_SRCS} scratch/test_tape.cpp)
        set_property(TARGET test_tape PROPERTY CXX_STANDARD 17)
        set_property(TARGET test_tape PROPERTY CXX_STANDARD_REQUIRED ON)
        target_link_libraries(test_tap ${targets} nlohmann_json::nlohmann_json)
        # target_link_libraries(test_tape ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${CURL_LIBRARIES} ${HDF5_LIBRARIES} ${PNETCDF_LIBRARIES} ${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${Zip_LIBRARIES} ${AEC_LIBRARIES} ${CMAKE_DL_LIBS})
    endif()
else()
    message("${BoldYellow}~~~ No target using -DTARGET specified. Using default 'simulation'")
    message("~~~ Possible targets: simulation timing scan${ColourReset}")
    add_executable(trajectories ${TYPES_SRCS} microphysics/trajectories.cpp )
    set_property(TARGET trajectories PROPERTY CXX_STANDARD 17)
    set_property(TARGET trajectories PROPERTY CXX_STANDARD_REQUIRED ON)
    target_link_libraries(trajectories ${targets} nlohmann_json::nlohmann_json)
    # if(MPI_CXX_LIBRARIES)
    #     target_link_libraries(trajectories ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${CURL_LIBRARIES} ${HDF5_LIBRARIES}  ${MPI_CXX_LIBRARIES} ${ZLIB_LIBRARIES} ${Zip_LIBRARIES} ${AEC_LIBRARIES} ${CMAKE_DL_LIBS})
    # else()
    #     target_link_libraries(trajectories ${GSL_LIBRARIES} ${NETCDF_LIBRARIES} ${CURL_LIBRARIES} ${HDF5_LIBRARIES} ${PNETCDF_LIBRARIES} ${ZLIB_LIBRARIES} ${Zip_LIBRARIES} ${AEC_LIBRARIES} ${CMAKE_DL_LIBS})
    # endif()
endif()



