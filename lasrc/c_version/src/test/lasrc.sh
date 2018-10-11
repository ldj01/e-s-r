#!/bin/bash
#
# Test lasrc.

base_scene="LC08_L1TP_040036_20141219_20160822_01_T1"
aux_file="L8ANC2014353.hdf_fused"

if [ "$#" -ne 1 ]; then
    echo "Usage:  $0 <path_to_lasrc>"
    exit 1
fi
bin_dir=$1

data_files=(${LEVEL2_UNIT_TEST_DATA}/espa-surface-reflectance/lasrc_ref/*)
input_dir=$LEVEL2_UNIT_TEST_DATA/espa-surface-reflectance/input/l8

mkdir -p lasrc && cd lasrc

cp $input_dir/*.img .
cp $input_dir/*.hdr .

sed -e s%LEVEL2_UNIT_TEST_DATA%${LEVEL2_UNIT_TEST_DATA}% \
    ${input_dir}/${base_scene}.xml > ${base_scene}.xml
$bin_dir/lasrc --xml=${base_scene}.xml --aux=${aux_file}
if [ $? -ne 0 ]; then
    echo "Error: lasrc processing failed."
    exit 1
fi

# Compare the results.
status=0
for i in "${data_files[@]}"; do
    base_name=`basename $i`

    echo "Comparing $base_name..."

    # If the file is an ASCII header file, use diff.  Otherwise, assume it's
    # a binary file, and use cmp to dump the octal differences.
    ext="${i##*.}"
    if [ "$ext" = "hdr" ]; then
        diff $i $base_name
        if [ $? -ne 0 ]; then
            echo "${base_name} differs from reference version."
            status=1
        fi
    else
        lines=`grep lines ${i%.img}.hdr | awk '{print $3}'`
        samples=`grep samples ${i%.img}.hdr | awk '{print $3}'`
        dtype=`grep "data type" ${i%.img}.hdr | awk '{print $4}'`
        if [ "$dtype" = "1" ]; then
            depth=8
        else
            depth=16
        fi
        numdiffs=`compare -metric AE -depth $depth -size ${samples}x${lines} \
                      gray:$i gray:${base_name} null: 2>&1`
        if [ $numdiffs != "0" ]; then
            echo "${base_name} differs from reference version in $numdiffs " \
                 "pixels."
            status=1
        fi
    fi
done

if [ $status -ne 0 ]; then
    echo "Test differences found."
    exit 1
fi

cd ..
rm -r lasrc

echo "Test completed successfully."
exit 0
