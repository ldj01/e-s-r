# The following tests use input and reference data from the
# $ESPA_UNIT_TEST_DATA_DIR repository.

add_test(do_lasrc_l1_angle_bands ${CMAKE_CURRENT_SOURCE_DIR}/do_lasrc.sh 1
                                 ${CMAKE_CURRENT_SOURCE_DIR}/..)
add_test(do_lasrc_no_l1_angle_bands ${CMAKE_CURRENT_SOURCE_DIR}/do_lasrc.sh 0
                                    ${CMAKE_CURRENT_SOURCE_DIR}/..)
set_property(TEST do_lasrc_l1_angle_bands do_lasrc_no_l1_angle_bands
             PROPERTY ENVIRONMENT
             PATH=${CMAKE_BINARY_DIR}/lasrc/c_version/src:$ENV{PATH})
