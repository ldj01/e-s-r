# The following tests use input and reference data from the
# $LEVEL2_UNIT_TEST_DATA repository.

add_test(lndsr ${CMAKE_CURRENT_SOURCE_DIR}/lndsr.sh
               ${CMAKE_CURRENT_BINARY_DIR}/..
               ${CMAKE_CURRENT_SOURCE_DIR}/../../lndpm/test/lndpm_ref)

set_property(TEST lndsr PROPERTY ENVIRONMENT
             PATH=${CMAKE_CURRENT_BINARY_DIR}/../../6sV-1.0B:$ENV{PATH})
