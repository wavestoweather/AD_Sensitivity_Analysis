# - Find CoDiPack

if (CODIPACK_INCLUDES)
    set(CODIPACK_FIND_QUIETLY TRUE)
endif()

find_package(PkgConfig)
pkg_check_modules(PC_CODIPACK QUIET CoDiPack)
set(CODIPACK_DEFINITIONS ${PC_CODIPACK_CFLAGS_OTHER})

find_path(CODIPACK_INCLUDE_DIR codi.hpp
          HINTS ${CODIPACK_INCLUDEDIR} ${PC_CODIPACK_INCLUDE_DIRS}
          PATH_SUFFIXES include )
include(FindPackageHandleStandardArgs)

# handle the QUIETLY and REQUIRED arguments and set CODIPACK_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(CoDiPack  DEFAULT_MSG
                                CODIPACK_INCLUDE_DIR)

mark_as_advanced(CODIPACK_INCLUDE_DIR CODIPACK_LIBRARY )

set(CODIPACK_INCLUDE_DIRS ${CODIPACK_INCLUDE_DIR} )
