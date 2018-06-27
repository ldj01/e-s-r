include_directories(${ESPA_INCLUDES})
include_directories(${XML_INCLUDES})
include_directories(${HDF4_INCLUDES})

link_directories(${ESPA_LIBS})

add_executable(combine_l8_aux_data get_args.c combine_l8_aux_data.c)
target_link_libraries(combine_l8_aux_data -lespa_raw_binary -lespa_common
                      -lespa_format_conversion
                      ${XML_LIBS} ${HDF4_LIBS} m)

install(TARGETS combine_l8_aux_data DESTINATION bin COMPONENT lasrc-aux
        EXCLUDE_FROM_ALL)
