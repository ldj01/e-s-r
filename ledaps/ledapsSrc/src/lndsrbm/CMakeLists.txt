include_directories(${ESPA_INCLUDES})
include_directories(${XML_INCLUDES})
include_directories(${HDFINC})
include_directories(${HDFEOS_INC})
include_directories(${ZLIBINC})
include_directories(${HDFEOS_GCTPINC})

link_directories(${ESPA_LIBS})

add_executable(lndsrbm lndsrbm.c)
target_link_libraries(lndsrbm -lespa_raw_binary -lespa_common
                      -lespa_format_conversion -lgctp3
                      ${XML_LIBS} ${HDF4_LIBS} ${HDFEOS_LIBS} m)

add_executable(comptemp comptemp.c)
target_link_libraries(comptemp -lespa_raw_binary -lespa_common
                      -lespa_format_conversion -lgctp3
                      ${XML_LIBS} ${HDF4_LIBS} ${HDFEOS_LIBS} m)

add_executable(dump_meta dump_meta.c)
target_link_libraries(dump_meta -lespa_raw_binary -lespa_common
                      -lespa_format_conversion -lgctp3
                      ${XML_LIBS} ${HDF4_LIBS} ${HDFEOS_LIBS} m)

add_executable(xy2geo LS_geoloc_driver.c LS_geoloc.c)
target_link_libraries(xy2geo -lespa_raw_binary -lespa_common
                      -lespa_format_conversion
                      ${XML_LIBS} ${HDF4_LIBS} ${HDFEOS_GCTPLIBS} m)

add_executable(geo2xy LS_geoloc_driver.c LS_geoloc.c)
target_link_libraries(geo2xy -lespa_raw_binary -lespa_common
                      -lespa_format_conversion
                      ${XML_LIBS} ${HDF4_LIBS} ${HDFEOS_GCTPLIBS} m)

add_executable(SDSreader3.0 SDSreader3.0.c)
target_link_libraries(SDSreader3.0 -lespa_raw_binary -lespa_common
                      -lespa_format_conversion -lgctp3
                      ${XML_LIBS} ${HDF4_LIBS} ${HDFEOS_LIBS} m)

install(TARGETS lndsrbm comptemp dump_meta xy2geo SDSreader3.0 DESTINATION bin)
install(PROGRAMS lndsrbm.py DESTINATION bin)
