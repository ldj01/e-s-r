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
input_dir=$LEVEL2_UNIT_TEST_DATA/espa-surface-reflectance/input_l7

rm -rf lndcal
mkdir lndcal && cd lndcal

cp $input_dir/* .
chmod u+w *

sed -e s%LEDAPS_AUX_DIR%${LEDAPS_AUX_DIR}% $pfile_dir/$pfile > pfile.local
$bin_dir/lndcal --pfile pfile.local
if [ $? -ne 0 ]; then
    echo "Error: lndcal processing failed."
    exit 1
fi

status=0
for i in "${data_files[@]}"; do
    base_i=`basename $i`

    # For the XML file, ignore the records that vary from one run to the next.
    if [ "$base_i" = "${base_scene}.xml" ]; then
        sed -e 's%<production_date>.*<%<production_date><%' $i > tmp1.xml
        sed -e 's%<production_date>.*<%<production_date><%' $base_i > tmp2.xml
        diff tmp1.xml tmp2.xml
    else
        diff $i $base_i
    fi
    if [ $? -ne 0 ]; then
        echo $base_i "differs from reference version."
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
