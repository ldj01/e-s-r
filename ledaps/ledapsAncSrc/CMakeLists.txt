include_directories(${NETCDF_INCLUDES})
include_directories(${HDF4_INCLUDES})

link_directories(${ESPA_LIBS})

add_executable(convert_ozone convert_ozone.c)
target_link_libraries(convert_ozone ${XML_LIBS} ${HDF4_LIBS} m)

add_executable(ncep_repackage ncep_reanalysis_surface_repackage.c)
target_link_libraries(ncep_repackage ${XML_LIBS} ${HDF4_LIBS} ${NETCDF_LIBS} m)

install(TARGETS convert_ozone ncep_repackage DESTINATION bin
        COMPONENT ledaps-aux EXCLUDE_FROM_ALL)
install(PROGRAMS updatencep.py updatetoms.py DESTINATION bin
        COMPONENT ledaps-aux EXCLUDE_FROM_ALL)
