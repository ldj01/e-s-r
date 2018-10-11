# The following tests use input and reference data from the
# $LEVEL2_UNIT_TEST_DATA repository.

add_test(lasrc ${CMAKE_CURRENT_SOURCE_DIR}/lasrc.sh
               ${CMAKE_CURRENT_BINARY_DIR}/..)
