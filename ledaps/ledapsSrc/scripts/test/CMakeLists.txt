# The following tests use input and reference data from the
# $ESPA_UNIT_TEST_DATA_DIR repository.

add_test(do_ledaps_l1_angle_bands ${CMAKE_CURRENT_SOURCE_DIR}/do_ledaps.sh 1
                                  ${CMAKE_CURRENT_SOURCE_DIR}/..)
add_test(do_ledaps_no_l1_angle_bands ${CMAKE_CURRENT_SOURCE_DIR}/do_ledaps.sh 0
                                     ${CMAKE_CURRENT_SOURCE_DIR}/..)
set_property(TEST do_ledaps_l1_angle_bands do_ledaps_no_l1_angle_bands
             PROPERTY ENVIRONMENT
             PATH=${CMAKE_BINARY_DIR}/ledaps/ledapsSrc/src/lndsr:${CMAKE_BINARY_DIR}/ledaps/ledapsSrc/src/lndcal:${CMAKE_BINARY_DIR}/ledaps/ledapsSrc/src/lndpm:$ENV{PATH})
