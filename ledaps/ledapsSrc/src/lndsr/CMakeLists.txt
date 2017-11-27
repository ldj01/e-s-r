include_directories(../lndpm)
include_directories(${ESPA_INCLUDES})
include_directories(${XML_INCLUDES})
include_directories(${HDFINC})
include_directories(${HDFEOS_INC})
include_directories(${ZLIBINC})

link_directories(${ESPA_LIBS})

add_executable(lndsr ar.c clouds.c date.c error.c grib.c input.c lndsr.c
               lut.c myhdf.c mystring.c output.c param.c prwv_input.c
               read_grib_tools.c sixs_runs.c sr.c CHAND.f CSALBR.f)
target_link_libraries(lndsr -lespa_raw_binary -lespa_common
                      -lespa_format_conversion -lgctp3
                      ${XML_LIBS} ${HDF4_LIBS} ${HDFEOS_LIBS} m)
install(TARGETS lndsr DESTINATION bin)
