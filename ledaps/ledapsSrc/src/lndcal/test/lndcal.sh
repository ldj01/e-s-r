#!/bin/bash
#
# Test lndcal.

base_scene="LE07_L1TP_043028_20020419_20180206_01_T1"
pfile="lndcal.${base_scene}.txt"

if [ "$#" -ne 2 ]; then
    echo "Usage:  $0 <path_to_lndcal> <pfile_path>"
    exit 1
fi
bin_dir=$1
pfile_dir=$2

data_files=(${LEVEL2_UNIT_TEST_DATA}/espa-surface-reflectance/lndcal_ref/*)
input_dir=$LEVEL2_UNIT_TEST_DATA/espa-surface-reflectance/input

mkdir -p lndcal && cd lndcal

cp $input_dir/* .

sed -e s%LEVEL2_UNIT_TEST_DATA%${input_dir}% \
    -e s%LEDAPS_AUX_DIR%${LEDAPS_AUX_DIR}% $pfile_dir/$pfile > pfile.local
$bin_dir/lndcal --pfile pfile.local
if [ $? -ne 0 ]; then
    echo "Error: lndcal processing failed."
    exit 1
fi

status=0
for i in "${data_files[@]}"; do
    diff $i `basename $i`
    if [ $? -ne 0 ]; then
        echo `basename $i` "differs from reference version."
        status=1
    fi
done

if [ $status -ne 0 ]; then
    echo "Test differences found."
    exit 1
fi

cd ..
rm -r lndcal

echo "Test completed successfully."
exit 0
