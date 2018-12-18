# The following tests use input and reference data from the
# $ESPA_UNIT_TEST_DATA_DIR repository.

add_test(do_ledaps ${CMAKE_CURRENT_SOURCE_DIR}/do_ledaps.sh
                   ${CMAKE_CURRENT_SOURCE_DIR}/..)
set_property(TEST do_ledaps PROPERTY ENVIRONMENT
             PATH=${CMAKE_BINARY_DIR}/ledaps/ledapsSrc/src/lndsr:${CMAKE_BINARY_DIR}/ledaps/ledapsSrc/src/lndcal:${CMAKE_BINARY_DIR}/ledaps/ledapsSrc/src/lndpm:$ENV{PATH})
