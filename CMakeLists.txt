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

project(trajectories LANGUAGES CXX C)

if(NOT WIN32)
    string(ASCII 27 Esc)
    set(ColourReset "${Esc}[m")
    set(ColourBold  "${Esc}[1m")
    set(Red         "${Esc}[31m")
    set(Green       "${Esc}[32m")
    set(Yellow      "${Esc}[33m")
    set(Blue        "${Esc}[34m")
    set(Magenta     "${Esc}[35m")
    set(Cyan        "${Esc}[36m")
    set(White       "${Esc}[37m")
    set(BoldRed     "${Esc}[1;31m")
    set(BoldGreen   "${Esc}[1;32m")
    set(BoldYellow  "${Esc}[1;33m")
    set(BoldBlue    "${Esc}[1;34m")
    set(BoldMagenta "${Esc}[1;35m")
    set(BoldCyan    "${Esc}[1;36m")
    set(BoldWhite   "${Esc}[1;37m")
endif()

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
    message(WARNING "~~~ You are not using a GNU compiler. This has not been tested")
endif()



# Check for build types such as -DBUILD_TYPE=debug
if ((NOT CMAKE_BUILD_TYPE) OR (CMAKE_BUILD_TYPE MATCHES "^None"))
    set(CMAKE_BUILD_TYPE "release")
endif()
boost_report_pretty("CMake build" CMAKE_BUILD_TYPE)
message(STATUS "Setting flags")

set(CXX_WARNING_FLAGS "-Wall -Wextra -Wpedantic")

if( NOT RELOPTLEVEL )
    set(RELOPTLEVEL 3)
endif()
set(CMAKE_BUILD_TYPE "${CMAKE_BUILD_TYPE}" CACHE STRING "Choose the type of build, options are: None debug release" FORCE)

set(CMAKE_CXX_FLAGS_SIM "-DRK4ICE -DAUTUMN -DFLUX -DMET3D -DSB_CONV -DSB_SHAPE")

# Set flags used in all build types
set(CMAKE_CXX_FLAGS "-pipe ${CXX_STANDARD} ${CXX_WARNING_FLAGS} ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_SIM} -DNPROCS=${NPROCS} -lm -DCODI_UseForcedInlines -fargument-noalias-global -ftree-loop-vectorize")
string(REGEX REPLACE "[ ]+" " " CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})

# Add release flags
set(CMAKE_CXX_FLAGS_RELEASE "-O${RELOPTLEVEL} -march=native")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}" CACHE STRING
  "Flags used by compiler during release builds")

# Add debug flags
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}" CACHE STRING
  "Flags used by compiler during debug builds")

if (SILENT_MODE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSILENT_MODE")
  message("${BoldYellow}~~~ Using silent mode. Suppress prints where possible.${ColourReset}")
endif()

if (COMPRESS_OUTPUT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DCOMPRESS_OUTPUT")
    message("${BoldYellow}~~~ Using compression (level 9) on all output columns.")
    message("~~~ Make sure you build NetCDF with HDF5 v1.10.2 or above.${ColourReset}")
endif()

# Different parameterization for collaboration with B8
if(B_EIGHT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DB_EIGHT")
    message("${BoldYellow}~~~ Using different CCN activation for ${ColourReset}")
    message("${BoldYellow}~~~ collaboration with B8: ${ColourReset}")
    message("${BoldYellow}~~~ Role of uncertainty in ice microphysical processes in warm conveyor belts ${ColourReset}")
endif()

if (TRUSTED_DATA)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRUSTED_DATA")
    message("${BoldYellow}~~~ Trusting the dataset.${ColourReset}")
    message("${BoldYellow}~~~ This way the program just stops simulations for every trajectory ${ColourReset}")
    message("${BoldYellow}~~~ with more than 10 consecutive NaNs for every input column ${ColourReset}")
    message("${BoldYellow}~~~ except longitude, latitude and time. ${ColourReset}")
    message("${BoldYellow}~~~ Without this option, the program assumes the dataset is broken ${ColourReset}")
    message("${BoldYellow}~~~ and terminates directly, giving information which ${ColourReset}")
    message("${BoldYellow}~~~ trajectory and time index starts with consecutive NaNs. ${ColourReset}")
endif()

# Add tracing flags
if (DEBUG_SEG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG_SEG")
    message("${BoldYellow}~~~ Output for debugging seg faults${ColourReset}")
endif()
if (TRACE_QV)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_QV")
    message("${BoldYellow}~~~ Tracing QV${ColourReset}")
endif()
if (TRACE_QC)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_QC")
    message("${BoldYellow}~~~ Tracing QC${ColourReset}")
endif()
if (TRACE_QR)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_QR")
    message("${BoldYellow}~~~ Tracing QR${ColourReset}")
endif()
if (TRACE_QI)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_QI")
    message("${BoldYellow}~~~ Tracing QI${ColourReset}")
endif()
if (TRACE_QS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_QS")
    message("${BoldYellow}~~~ Tracing QS${ColourReset}")
endif()
if (TRACE_QG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_QG")
    message("${BoldYellow}~~~ Tracing QG${ColourReset}")
endif()
if (TRACE_QH)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_QH")
    message("${BoldYellow}~~~ Tracing QH${ColourReset}")
endif()
if (TRACE_ENV)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_ENV")
    message("${BoldYellow}~~~ Tracing environment variables${ColourReset}")
endif()
if (TRACE_TIME)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_TIME")
    message("${BoldYellow}~~~ Tracing time step${ColourReset}")
endif()
if (TRACE_SAT)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_SAT")
    message("${BoldYellow}~~~ Tracing saturation adjustment${ColourReset}")
endif()
if (TRACE_COMM)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_COMM")
    message("${BoldYellow}~~~ Tracing MPI communication${ColourReset}")
endif()
if (TRACE_WRITE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_WRITE")
    message("${BoldYellow}~~~ Tracing writing output${ColourReset}")
endif()
if (DEVELOP)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEVELOP")
    message("${BoldYellow}~~~ Printing output for development${ColourReset}")
endif()
if (TRACE_READER)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DTRACE_READER")
    message("${BoldYellow}~~~ Printing output for netcdf_reader_t${ColourReset}")
endif()
if (LINE_INFO)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
    message("${BoldYellow}~~~ Compiling with -g.${ColourReset}")
endif()

message(STATUS "Searching for necessary libraries")

include_directories (${AD_SIM_HOME})

message("~~~ MPI")
find_package(MPI)
if( MPI_CXX_FOUND )
    include_directories(${MPI_CXX_INCLUDE_PATH})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_MPI")
    set(targets "${MPI_CXX_LIBRARIES}")
else()
    message(FATAL_ERROR
"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CMake was not able to find your MPI installation.
Please define MPI_HOME pointing to your installation. If this does not
work, take a look at https://cmake.org/cmake/help/latest/module/FindMPI.html.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
")
endif()

message("~~~ BOOST")
find_package(Boost 1.62 REQUIRED)

message("~~~ NLohmann/JSON")
find_package(nlohmann_json 3.9.1 REQUIRED)
# message(STATUS "${nlohmann_json::nlohmann_json}")

# set(targets "${targets};${nlohmann_json_LIBRARIES}")

message("~~~ zlib")
message(STATUS "Make sure this is the same zlib library that is used by HDF5")
find_package(ZLIB)
# Hack if for some reason it is not found in the usual directories
if( NOT ${ZLIB_FOUND} )
    message(WARNING "Did not find zlib. Looking at different directories now.")
    find_package(ZLIB PATHS ${HDF5_DIR})
    if( "${ZLIB_INCLUDE_DIR}" STREQUAL "" )
        message(ERROR "Did not find zlib!")
    endif()
    message(STATUS "Found zlib include: ${ZLIB_INCLUDE_DIR}")
    include_directories (${ZLIB_INCLUDE_DIR})
    if( "${ZLIB_LIBRARIES}" STREQUAL "" )
        # Hack if for some reason this is still empty
        string(REPLACE "include" "lib" ZLIB_LIBRARIES ${ZLIB_INCLUDE_DIR})
        set(ZLIB_LIBRARIES "${ZLIB_LIBRARIES}/libaec.a;${ZLIB_LIBRARIES}/libsz.a;${ZLIB_LIBRARIES}/libz.a")
    endif()
endif()
set(targets "${targets};${ZLIB_LIBRARIES}")

message("~~~ Zip")
find_library(Zip_LIBRARIES NAMES zip)
if( "${Zip_LIBRARIES}" STREQUAL "Zip_LIBRARIES-NOTFOUND" )
    message(STATUS "Zip not found. Make sure that HDF5 does not need Zip.")
    set(Zip_LIBRARIES "")
    set(AEC_LIBRARIES "")
else()
    message(STATUS "Found Zip: ${Zip_LIBRARIES}")
    message("~~~ AEC")
    find_library(AEC_LIBRARIES NAMES aec)
    message(STATUS "Found AEC: ${AEC_LIBRARIES}")
    set(targets "${targets};${Zip_LIBRARIES};${AEC_LIBRARIES}")
endif()

message("~~~ CURL")
find_package(CURL REQUIRED)
set(targets "${targets};${CURL_LIBRARIES}")

message("~~~ PnetCDF")
find_package(PnetCDF)
include_directories (${PNETCDF_INCLUDE_DIR})
if( PnetCDF_FOUND )
    set(targets "${targets};${PNETCDF_LIBRARIES}")
endif()

set(CMAKE_FIND_DEBUG_MODE TRUE)
message("~~~ NetCDF")
set(NETCDF "Yes")
find_package(NetCDF)
if( NetCDF_FOUND )
    include_directories (${NETCDF_INCLUDE_DIR})
    # message(STATUS "Found NetCDF include: ${NETCDF_INCLUDE_DIR}")
    set(targets "${targets};${NETCDF_LIBRARIES}")
else()
    message(FATAL_ERROR
"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
CMake was not able to find your NetCDF installation (maybe it wasn't build using CMake).
Please define NETCDF_DIR pointing to your installation as environment variable.
Make sure that -lnetcdf can be used by your compiler.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
")
endif()
set(CMAKE_FIND_DEBUG_MODE FALSE)

message("~~~ HDF5")
# I usually prefer the Config mode for CMake to look for HDF5 but for our supercomputer
# I have to use the module mode, whereas I want to use Config mode for my
# local machine and my docker container.
# set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)
set(HDF5_HINT ${HDF5_DIR})
find_package(HDF5)
# Hack if for some reason it is not found in the usual directories
if (NOT ${HDF5_FOUND} )
    message(WARNING "Could not find HDF5. Looking at other directories.")
    find_package(HDF5 REQUIRED PATHS ${HDF5_DIR})
    include_directories (${HDF5_INCLUDE_DIR})
    message(STATUS "HDF5 include dir: ${HDF5_INCLUDE_DIR}")
    if( "${HDF5_LIBRARIES}" STREQUAL "" )
        # Hack if for some reason this is still empty
        string(REPLACE "include" "lib" HDF5_LIBRARIES ${HDF5_INCLUDE_DIR})
        set(HDF5_LIBRARIES "${HDF5_LIBRARIES}/libhdf5_hl.a;${HDF5_LIBRARIES}/libhdf5.a")
    endif()
else()
    # Add high level library that are hard to find with CMake for some reason
    if( "${HDF5_HL_LIBRARIES}" STREQUAL "" )
        list(GET HDF5_INCLUDE_DIR 0 HEAD)
        string(REPLACE "include" "lib" HDF5_HL_LIBRARIES ${HEAD})
        set(HDF5_HL_LIBRARIES "${HDF5_HL_LIBRARIES}/libhdf5_hl.a;${HDF5_HL_LIBRARIES}/libhdf5.a")
    endif()
    set(targets "${targets};${HDF5_HL_LIBRARIES}")
endif()
message(STATUS "HDF5 Version: ${HDF5_VERSION}")
set(targets "${targets};${HDF5_LIBRARIES}")

message("~~~ CoDiPack")
find_package(CoDiPack)
if( CODIPACK_FOUND )
    include_directories (${CODIPACK_INCLUDE_DIR})
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
find_package(GSL REQUIRED)
include_directories (${GSL_INCLUDE_DIR})
set(targets "${targets};${GSL_LIBRARIES}")

include_directories(${AD_SIM_HOME})

message(STATUS "Adding executable and linking")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}" CACHE STRING
    "Flags used by the compiler during all build types")

set(targets "${targets};${CMAKE_DL_LIBS}")
message(STATUS ${targets})
add_subdirectory(src)

