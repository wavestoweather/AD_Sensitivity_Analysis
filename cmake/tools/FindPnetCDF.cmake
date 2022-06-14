# - Find PnetCDF
# Find the native PnetCDF includes and library
#
#  PNETCDF_INCLUDES    - where to find netcdf.h, etc
#  PNETCDF_LIBRARIES   - Link these libraries when using NetCDF
#  PNETCDF_FOUND       - True if NetCDF found including required interfaces (see below)
#
# Your package can require certain interfaces to be FOUND by setting these
#
#  NETCDF_CXX         - require the C++ interface and link the C++ library
#  NETCDF_F77         - require the F77 interface and link the fortran library
#  NETCDF_F90         - require the F90 interface and link the fortran library
#
# The following are not for general use and are included in
# NETCDF_LIBRARIES if the corresponding option above is set.
#
#  NETCDF_LIBRARIES_C    - Just the C interface
#  NETCDF_LIBRARIES_CXX  - C++ interface, if available
#  NETCDF_LIBRARIES_F77  - Fortran 77 interface, if available
#  NETCDF_LIBRARIES_F90  - Fortran 90 interface, if available
#
# Normal usage would be:
#  set (NETCDF_F90 "YES")
#  find_package (NetCDF REQUIRED)
#  target_link_libraries (uses_f90_interface ${NETCDF_LIBRARIES})
#  target_link_libraries (only_uses_c_interface ${NETCDF_LIBRARIES_C})

if (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)
  # Already in cache, be silent
  set (PNETCDF_FIND_QUIETLY TRUE)
endif (PNETCDF_INCLUDES AND PNETCDF_LIBRARIES)

find_path (PNETCDF_INCLUDES pnetcdf.h
  HINTS PNETCDF_DIR ENV PNETCDF_DIR)

find_library (PNETCDF_LIBRARIES_C       NAMES pnetcdf)
mark_as_advanced(PNETCDF_LIBRARIES_C)

set (PNetCDF_has_interfaces "YES") # will be set to NO if we're missing any interfaces
set (PNetCDF_libs "${PNETCDF_LIBRARIES_C}")

get_filename_component (PNetCDF_lib_dirs "${PNETCDF_LIBRARIES_C}" PATH)

macro (PNetCDF_check_interface lang header libs)
  if (PNETCDF_${lang})
    find_path (PNETCDF_INCLUDES_${lang} NAMES ${header}
      HINTS "${PNETCDF_INCLUDES}" NO_DEFAULT_PATH)
    find_library (PNETCDF_LIBRARIES_${lang} NAMES ${libs}
      HINTS "${PNetCDF_lib_dirs}" NO_DEFAULT_PATH)
    mark_as_advanced (PNETCDF_INCLUDES_${lang} PNETCDF_LIBRARIES_${lang})
    if (PNETCDF_INCLUDES_${lang} AND PNETCDF_LIBRARIES_${lang})
      list (INSERT PNetCDF_libs 0 ${PNETCDF_LIBRARIES_${lang}}) # prepend so that -lnetcdf is last
    else (PNETCDF_INCLUDES_${lang} AND PNETCDF_LIBRARIES_${lang})
      set (PNetCDF_has_interfaces "NO")
      message (STATUS "Failed to find PnetCDF interface for ${lang}")
    endif (PNETCDF_INCLUDES_${lang} AND PNETCDF_LIBRARIES_${lang})
  endif (PNETCDF_${lang})
endmacro (PNetCDF_check_interface)

PNetCDF_check_interface (CXX pnetcdfcpp.h pnetcdf_c++)
PNetCDF_check_interface (F77 pnetcdf.inc  pnetcdff)
PNetCDF_check_interface (F90 pnetcdf.mod  pnetcdff)

set (PNETCDF_LIBRARIES "${PNetCDF_libs}" CACHE STRING "All PnetCDF libraries required for interface level")

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (PNetCDF DEFAULT_MSG PNETCDF_LIBRARIES PNETCDF_INCLUDES PNetCDF_has_interfaces)

mark_as_advanced (PNETCDF_LIBRARIES PNETCDF_INCLUDES)
