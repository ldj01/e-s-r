include_directories(${ESPA_INCLUDES})
include_directories(${XML_INCLUDES})
include_directories(${HDFINC})

link_directories(${ESPA_LIBS})

file(GLOB SOURCES "*.c")

add_executable(lndcal ${SOURCES})
target_link_libraries(lndcal -lespa_raw_binary -lespa_common
                      -lespa_format_conversion
                      ${XML_LIBS} ${HDF4_LIBS} m)
install(TARGETS lndcal DESTINATION bin)
