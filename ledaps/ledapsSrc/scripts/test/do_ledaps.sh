#!/bin/bash
#
# Test the ledaps Python script.

base_scene="LE07_L1TP_043028_20020419_20180206_01_T1"
pfile="lndsr.${base_scene}.txt"

if [ "$#" -ne 1 ]; then
    echo "Usage:  $0 <path_to_ledaps>"
    exit 1
fi
bin_dir=$1

data_files=(${ESPA_UNIT_TEST_DATA_DIR}/espa-surface-reflectance/lndsr_ref/*)
input_dir=$ESPA_UNIT_TEST_DATA_DIR/espa-surface-reflectance/input_l7

rm -rf ledaps
mkdir ledaps && cd ledaps

ln -s $input_dir/*.img .
ln -s $input_dir/*.hdr .
cp $input_dir/${base_scene}.xml .
chmod u+w ${base_scene}.xml

$bin_dir/do_ledaps.py --xml=${base_scene}.xml --use_l1_angle_bands
if [ $? -ne 0 ]; then
    echo "Error: ledaps processing failed."
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
rm -r ledaps

echo "Test completed successfully."
exit 0
