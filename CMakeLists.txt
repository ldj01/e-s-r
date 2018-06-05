cmake_minimum_required(VERSION 3.6)

project(espa-surface-reflectance C Fortran)
set (${PROJECT_NAME}_VERSION_MAJOR 1)
set (${PROJECT_NAME}_VERSION_MINOR 0)
set (${PROJECT_NAME}_VERSION_PATCH 0)
set (VERSION ${${PROJECT_NAME}_VERSION_MAJOR}.${${PROJECT_NAME}_VERSION_MINOR}.${${PROJECT_NAME}_VERSION_PATCH})

set (PACKAGE ${PROJECT_NAME})

set(CPACK_PACKAGE_VERSION ${VERSION})
set(CPACK_GENERATOR "TGZ")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_IGNORE_FILES "/build/;/.git/;/.gz/;~$;${CPACK_SOURCE_IGNORE_FILES}")
include(CPack)

if (XML2INC)
    set (XML_INCLUDES ${XML2INC})
    set (XML_LIBS ${XML2LIB})
else ()
    find_package(LibXml2)
    if (LIBXML2_FOUND)
        set (XML_INCLUDES ${LIBXML2_INCLUDE_DIR})
        set (XML_LIBS ${LIBXML2_LIBRARIES})
    else ()
        set (XML_INCLUDES $ENV{XML2INC})
        set (XML_LIBS -L$ENV{XML2LIB} -lxml2)
    endif ()
endif ()
set (XML_INCLUDES ${XML_INCLUDES} CACHE PATH "XML2 header file path" FORCE)
set (XML_LIBS ${XML_LIBS} CACHE STRING "XML2 libraries" FORCE)

if (HDFEOS_INC)
    set (HDFEOS_INCLUDES ${HDFEOS_INC})
    set (HDFEOS_LIBS ${HDFEOS_LIB})
else ()
    pkg_check_modules(HDFEOS hdf-eos4)
    if (HDFEOS_FOUND)
        set (HDFEOS_INCLUDES ${HDFEOS_INCLUDE_DIRS})
        set (HDFEOS_LIBS ${HDFEOS_LDFLAGS})
    else ()
        set (HDFEOS_INCLUDES $ENV{HDFEOS_INC})
        set (HDFEOS_LIBS -L$ENV{HDFEOS_LIB} -lhdfeos)
    endif ()
endif ()
set (HDFEOS_INCLUDES ${HDFEOS_INCLUDES} CACHE PATH "HDFEOS header file path"
     FORCE)
set (HDFEOS_LIBS ${HDFEOS_LIBS} CACHE STRING "HDFEOS libraries" FORCE)

if (NOT GCTP_INCLUDES)
    set (GCTP_INCLUDES $ENV{HDFEOS_GCTPINC})
endif ()
if (NOT GCTP_LIBS)
    set (GCTP_LIBS "-L$ENV{HDFEOS_GCTPLIB} -lGctp -lm")
endif ()
set (GCTP_INCLUDES ${GCTP_INCLUDES} CACHE PATH "GCTP header file path" FORCE)
set (GCTP_LIBS ${GCTP_LIBS} CACHE STRING "GCTP libraries" FORCE)

if (HDFINC)
    set (HDF4_INCLUDES ${HDFINC})
else ()
    set (HDF4_INCLUDES $ENV{HDFINC})
endif ()
if (HDFLIB)
    set (HDF4_LIBS "-L${HDFLIB} -lmfhdf -ldf -ljpeg -lz -lm")
else ()
    set (HDF4_LIBS "-L$ENV{HDFLIB} -lmfhdf -ldf -ljpeg -lz -lm")
endif ()
set (HDF4_INCLUDES ${HDF4_INCLUDES} CACHE PATH "HDF 4 header file path" FORCE)
set (HDF4_LIBS ${HDF4_LIBS} CACHE STRING "HDF 4 libraries" FORCE)

if (NCDF4INC)
    set (NETCDF_INCLUDES ${NCDF4INC})
    set (NETCDF_LIBS "-L${NCDF4LIB} -lnetcdf")
else ()
    pkg_check_modules(NETCDF netcdf)
    if (NETCDF_FOUND)
        set(NETCDF_INCLUDES ${NETCDF_INCLUDE_DIRS})
        set(NETCDF_LIBS ${NETCDF_LDFLAGS})
    else ()
        set (NETCDF_INCLUDES $ENV{NCDF4INC})
        set (NETCDF_LIBS "-L$ENV{NCDF4LIB} -lnetcdf")
    endif ()
endif ()
set (NETCDF_INCLUDES ${NETCDF_INCLUDES} CACHE PATH "NetCDF header file path"
     FORCE)
set (NETCDF_LIBS ${NETCDF_LIBS} CACHE STRING "NetCDF libraries" FORCE)

set (ESPA_INCLUDES $ENV{ESPAINC} CACHE PATH "ESPA header file path" FORCE)
set (ESPA_LIBS $ENV{ESPALIB} CACHE PATH "ESPA library path" FORCE)

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
