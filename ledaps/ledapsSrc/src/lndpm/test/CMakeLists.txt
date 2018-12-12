# The following tests use input and reference data from the
# $ESPA_UNIT_TEST_DATA_DIR repository.

add_test(lndpm ${CMAKE_CURRENT_SOURCE_DIR}/lndpm.sh
               ${CMAKE_CURRENT_BINARY_DIR}/..
               ${CMAKE_CURRENT_SOURCE_DIR}/lndpm_ref)
