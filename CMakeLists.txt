cmake_minimum_required(VERSION 3.7.2)
if("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
    message(FATAL_ERROR "
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
~~~~                     Do *not* build in place.                      ~~~~
~~~~       Run cmake /path/to/source from a different directory        ~~~~
~~~~ Remove the newly created CMakeCache.txt in your source directory! ~~~~
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
")
endif("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

project(trajectories LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set all tools needed from subdirectory tools for searching libraries
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/tools)
include(BoostUtils)

# set convenient variables
set(AD_SIM_HOME ${CMAKE_SOURCE_DIR})

# set cmake flags
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

# default system packages
option(SYSTEM_PACKAGES "Use tools provided by the operating system" ON)

# Get system info
# Might be useful later
set(OSTYPE ${CMAKE_SYSTEM_NAME})
boost_report_value(OSTYPE)

set(OSVERSION ${CMAKE_SYSTEM_VERSION})
boost_report_value(OSVERSION)

set(ARCH ${CMAKE_SYSTEM_PROCESSOR})
boost_report_value(ARCH)

set(BUILDNAME "${OSTYPE}-${OSVERSION}/${ARCH}/${COMPILER_ID_TAG}" CACHE INTERNAL "buildname")
boost_report_value(BUILDNAME)

set(TOOLSET "${COMPILER_ID_TAG}/${ARCH}/${CMAKE_BUILD_TYPE}" CACHE INTERNAL "toolset")
boost_report_value(TOOLSET)

execute_process(COMMAND hostname
    COMMAND tr -d \\n
    OUTPUT_VARIABLE HOSTNAME)
boost_report_value(HOSTNAME)
set(SITE ${HOSTNAME})

# report cmake path and version
boost_report_pretty("CMake path" CMAKE_COMMAND)
if(NOT CMAKE_VERSION)
  set(CMAKE_VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}.${CMAKE_PATCH_VERSION})
endif(NOT CMAKE_VERSION)
math(EXPR CMAKE_VERSION_INT "${CMAKE_MAJOR_VERSION} * 10000 + ${CMAKE_MINOR_VERSION} * 100 + ${CMAKE_PATCH_VERSION}")
boost_report_pretty("CMake version" CMAKE_VERSION)

# Check gcc version and report a warning in case this had not been tested
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
    if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "6.3.0")
        message(WARNING "
    ~~~ GCC version below 6.3.0 has not been tested.
    ~~~ Fingers crossed this works.")
        # deprecation_warning(20160905 "Unsupported gcc version.")
    endif()
else()
    message(WARNING "~~~ You are not using GNU gcc. This has not been tested")
endif()

message(STATUS "Setting flags")

# Set flags
set(CXX_WARNING_FLAGS "-Wall")

# Check for build types such as -DCMAKE_BUILD_TYPE=debug
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE MATCHES "^None"))
    set(CMAKE_BUILD_TYPE "release")
endif()
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING "Choose the type of build, options are: None debug release" FORCE)

set(CMAKE_CXX_FLAGS_SIM "-DRK4ICE -DAUTUMN -DFLUX -DMET3D -DSB_CONV -DSB_SHAPE")
if( NPROCS MATCHES "^[0-9]+$")
    message("~~~ Using NPROCS=${NPROCS} (used with GNU parallel and not MPI!)")
else()
    message(WARNING "
    ~~~ A target number of processes has not been set.
    ~~~ Defaults to NPROCS=4. This is only important when
    ~~~ using GNU parallel instead of MPI.")
    set(NPROCS "4")
endif()

# Set flags used in all build types
set(CMAKE_CXX_FLAGS "-pipe ${CXX_STANDARD} ${CXX_WARNING_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_SIM} -DNPROCS=${NPROCS} -lgsl -lgslcblas -lm -DCODI_UseForcedInlines -fargument-noalias-global -ftree-loop-vectorize -lnetcdf_c++4 -lnetcdf")
string(REGEX REPLACE "[ ]+" " " CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING
    "Flags used by the compiler during all build types")

# Add release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O${RELOPTLEVEL} -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING
  "Flags used by compiler during release builds")

# Add debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -Og")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING
  "Flags used by compiler during debug builds")

message(STATUS "Searching for necessary libraries")

include_directories (${AD_SIM_HOME})

message("~~~ MPI")
find_package(MPI)

message("~~~ NetCDF")
set(NETCDF_CXX "Yes")
find_package(NetCDF)
if( NetCDF_found )
    include_directories (${NETCDF_INCLUDE_DIR})
else()
    if( NETCDF_INCLUDE_DIR )
        message(STATUS "Defaulting to NetCDF Include directory = ${NETCDF_INCLUDE_DIR}")
    else()
        message(FATAL_ERROR
"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CMake was not able to find your NetCDF installation (maybe it wasn't build using CMake).
Please define NETCDF_INCLUDE_DIR pointing to your installation.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
")
    endif()
endif()

message("~~~ CoDiPack")
find_package(CoDiPack)
if( CODIPACK_FOUND )
    include_directories (${CODIPACK_INCLUDE_DIRS})
else()
    message(FATAL_ERROR
"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CMake was not able to find your CoDiPack installation.
Please define CODIPACK_INCLUDEDIR pointing to 'codi.hpp'.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
")
endif()

message("~~~ GSL")
find_package(GSL)
include_directories (${GSL_INCLUDE_DIR})

# include the source directory for headers
include_directories(${AD_SIM_HOME})

# link_directories(${AD_SIM_HOME})
# include_directories(${AD_SIM_HOME}/include)
# include_directories(${AD_SIM_HOME}/include/types)
# include_directories(${AD_SIM_HOME}/include/microphysics)

message(STATUS "Adding executable and linking")

add_subdirectory(src)





# include_directories("${PROJECT_SOURCE_DIR}")