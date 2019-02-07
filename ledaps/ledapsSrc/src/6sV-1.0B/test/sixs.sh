#!/bin/bash
#
# Test 6S.

if [ "$#" -ne 2 ]; then
    echo "Usage:  $0 <path_to_sixsV1.0B> <input_data_path>"
    exit 1
fi
bin_dir=$1
input_dir=$2

$bin_dir/sixsV1.0B < ${input_dir}/sixs_input > sixs_output
if [ $? -ne 0 ]; then
    echo "Error: sixsV1.0B processing failed."
    exit 1
fi

# Compare the results.
status=0
diff ${input_dir}/sixs_output.ref sixs_output
if [ $? -ne 0 ]; then
    echo "sixs_output differs from reference version."
    status=1
fi

if [ $status -ne 0 ]; then
    echo "Test differences found."
    exit 1
fi

rm sixs_output

echo "Test completed successfully."
exit 0
