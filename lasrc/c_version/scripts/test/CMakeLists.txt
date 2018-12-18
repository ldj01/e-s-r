# The following tests use input and reference data from the
# $ESPA_UNIT_TEST_DATA_DIR repository.

add_test(do_lasrc ${CMAKE_CURRENT_SOURCE_DIR}/do_lasrc.sh
                  ${CMAKE_CURRENT_SOURCE_DIR}/..)
set_property(TEST do_lasrc PROPERTY ENVIRONMENT
             PATH=${CMAKE_BINARY_DIR}/lasrc/c_version/src:$ENV{PATH})
