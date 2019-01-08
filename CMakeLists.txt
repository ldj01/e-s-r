cmake_minimum_required(VERSION 3.6)

project(espa-surface-reflectance C Fortran)

# Set the package version.
include(package_version.cmake)

include(FindPkgConfig)

if (NOT XML_INCLUDES)
    find_package(LibXml2)
    if (LIBXML2_FOUND)
        set (XML_INCLUDES ${LIBXML2_INCLUDE_DIR})
        set (XML_LIBS ${LIBXML2_LIBRARIES})
    endif ()
endif ()
set (XML_INCLUDES ${XML_INCLUDES} CACHE PATH "XML2 header file path" FORCE)
set (XML_LIBS ${XML_LIBS} CACHE STRING "XML2 libraries" FORCE)

# Check whether any HDFEOS flags or environment variables are set. If so,
# then use the HDFEOS library (including its GCTP library).  Otherwise,
# use the GCTP2 library.
if (HDFEOS_LIBS)
    set (HDFEOS_LIBS "-L${HDFEOS_LIBS} -lhdfeos")
    set (USE_HDFEOS TRUE)
elseif (DEFINED ENV{HDFEOS_LIB})
    set (HDFEOS_LIBS -L$ENV{HDFEOS_LIB} -lhdfeos)
    set (USE_HDFEOS TRUE)
elseif (USE_HDFEOS)
    pkg_check_modules(HDFEOS hdf-eos4)
    if (HDFEOS_FOUND)
        set (HDFEOS_LIBS ${HDFEOS_LDFLAGS})
    endif ()
endif ()
if (USE_HDFEOS)
    message(STATUS "Using HDFEOS library.")
    set (HDFEOS_LIBS ${HDFEOS_LIBS} CACHE STRING "HDFEOS libraries" FORCE)
    if (HDFEOS_GCTPLIB)
        set (GCTP_LIBS -L${HDFEOS_GCTPLIB} -lGctp -lm)
    elseif (DEFINED ENV{HDFEOS_GCTPLIB})
        set (GCTP_LIBS -L$ENV{HDFEOS_GCTPLIB} -lGctp -lm)
    endif ()
else ()
    message(STATUS "Using GCTP2 package.")
    if (NOT GCTP_LIBS)
        set (GCTP_LIBS "-L$ENV{GCTP2LIB} -lgctp -lm")
    endif ()
endif ()
set (GCTP_LIBS ${GCTP_LIBS} CACHE STRING "GCTP2 libraries" FORCE)

# Set a flag identifying whether to use the HDFEOS library.
# If not, the GCTP2 library will be used (default).
set (USE_HDFEOS ${USE_HDFEOS} CACHE BOOL "Use the HDFEOS library for GCTP?"
     FORCE)

if (NOT SZIP_INCLUDES)
    set (SZIP_INCLUDES $ENV{SZIPINC})
endif ()
if (NOT SZIP_LIBS)
    set (SZIP_LIBS "-L$ENV{SZIPLIB} -lsz")
endif ()
set (SZIP_INCLUDES ${SZIP_INCLUDES} CACHE PATH "SZip header file path" FORCE)
set (SZIP_LIBS ${SZIP_LIBS} CACHE STRING "SZIP libraries" FORCE)

if (NOT HDF4_INCLUDES)
    set (HDF4_INCLUDES $ENV{HDFINC})
endif ()
if (NOT HDF4_LIBS)
    set (HDF4_LIBS "-L$ENV{HDFLIB} -lmfhdf -ldf -ljpeg -lz -lm ${SZIP_LIBS}")
endif ()
if (NOT HDF4_LIBRARY_DIRS)
    set (HDF4_LIBRARY_DIRS $ENV{HDFLIB})
endif ()
set (HDF4_INCLUDES ${HDF4_INCLUDES} CACHE PATH "HDF 4 header file path" FORCE)
set (HDF4_LIBS ${HDF4_LIBS} CACHE STRING "HDF 4 libraries" FORCE)
set (HDF4_LIBRARY_DIRS ${HDF4_LIBRARY_DIRS} CACHE STRING "HDF 4 library path"
     FORCE)

if (NOT HDF5_INCLUDES)
    set (HDF5_INCLUDES $ENV{HDF5INC})
endif ()
if (NOT HDF5_LIBS)
    set (HDF5_LIBS "-L$ENV{HDF5LIB} -lhdf5_hl -lhdf5 -ldl")
endif ()
if (NOT HDF5_LIBRARY_DIRS)
    set (HDF5_LIBRARY_DIRS $ENV{HDF5LIB})
endif ()
set (HDF5_INCLUDES ${HDF5_INCLUDES} CACHE PATH "HDF 5 header file path" FORCE)
set (HDF5_LIBS ${HDF5_LIBS} CACHE STRING "HDF 5 libraries" FORCE)
set (HDF5_LIBRARY_DIRS ${HDF5_LIBRARY_DIRS} CACHE STRING "HDF 5 library path"
     FORCE)

if (NOT NETCDF_INCLUDES)
    if (DEFINED ENV{NCDF4INC})
        set (NETCDF_INCLUDES $ENV{NCDF4INC})
        set (NETCDF_LIBS "-L$ENV{NCDF4LIB} -lnetcdf")
        set (NETCDF_LIBRARY_DIRS $ENV{NCDF4LIB})
    else ()
        pkg_check_modules(NETCDF netcdf)
        if (NETCDF_FOUND)
            set (NETCDF_INCLUDES ${NETCDF_INCLUDE_DIRS})
            set (NETCDF_LIBS "${NETCDF_LDFLAGS}")
        endif ()
    endif ()
elseif (NOT NETCDF_LIBS)
    # If NETCDF_INCLUDES is set, it's assumed that NETCDF_LIBS is also set.
    message(FATAL_ERROR "NETCDF_LIBS must be set if NETCDF_INCLUDES is set.")
endif ()
set (NETCDF_INCLUDES ${NETCDF_INCLUDES} CACHE PATH "NetCDF header file path"
     FORCE)
set (NETCDF_LIBS ${NETCDF_LIBS} CACHE STRING "NetCDF libraries" FORCE)
set (NETCDF_LIBRARY_DIRS ${NETCDF_LIBRARY_DIRS} CACHE STRING
     "NETCDF library path" FORCE)

# Set the paths for the ESPA headers and libraries.
if (NOT ESPA_INCLUDES)
    if (DEFINED ENV{ESPAINC})
        set (ESPA_INCLUDES $ENV{ESPAINC})
    else ()
        set (ESPA_INCLUDES ${CMAKE_INSTALL_PREFIX}/include)
    endif ()
endif ()
if (NOT ESPA_LIBS)
    if (DEFINED ENV{ESPALIB})
        set (ESPA_LIBS $ENV{ESPALIB})
    else ()
        set (ESPA_LIBS ${CMAKE_INSTALL_PREFIX}/lib)
    endif ()
endif ()
set (ESPA_INCLUDES ${ESPA_INCLUDES} CACHE PATH "ESPA header file path" FORCE)
set (ESPA_LIBS ${ESPA_LIBS} CACHE PATH "ESPA library path" FORCE)

# Control application threading.
if (DEFINED ENV{ENABLE_THREADING})
    add_definitions(-fopenmp)
    set (OPENMP_LIBS "-lgomp")
endif ()

set (CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

enable_testing()

add_subdirectory(scripts)
add_subdirectory(lasrc)
add_subdirectory(ledaps)

add_custom_target(lasrc-aux
                  COMMAND "${CMAKE_COMMAND}"
                          -DCMAKE_INSTALL_COMPONENT=lasrc-aux
                          -P "${CMAKE_BINARY_DIR}/cmake_install.cmake")

add_custom_target(ledaps-aux
                  COMMAND "${CMAKE_COMMAND}"
                          -DCMAKE_INSTALL_COMPONENT=ledaps-aux
                          -P "${CMAKE_BINARY_DIR}/cmake_install.cmake")
