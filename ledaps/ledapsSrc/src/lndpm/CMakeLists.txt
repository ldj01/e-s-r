include_directories(${ESPA_INCLUDES})
include_directories(${XML_INCLUDES})
include_directories(${HDFINC})

link_directories(${ESPA_LIBS})

add_executable(lndpm lndpm.c)
target_link_libraries(lndpm -lespa_raw_binary -lespa_common
                      -lespa_format_conversion
                      ${XML_LIBS} ${HDF4_LIBS} m)
install(TARGETS lndpm DESTINATION bin)

add_subdirectory(test)
