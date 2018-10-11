# The following tests use input and reference data from the
# $LEVEL2_UNIT_TEST_DATA repository.

add_test(lndcal ${CMAKE_CURRENT_SOURCE_DIR}/lndcal.sh
               ${CMAKE_CURRENT_BINARY_DIR}/..
               ${CMAKE_CURRENT_SOURCE_DIR}/../../lndpm/test/lndpm_ref)
