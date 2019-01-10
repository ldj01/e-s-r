#!/bin/bash
#
# Test lndsr.

base_scene="LE07_L1TP_043028_20020419_20180206_01_T1"
pfile="lndsr.${base_scene}.txt"

if [ "$#" -ne 2 ]; then
    echo "Usage:  $0 <path_to_lndsr> <pfile_path>"
    exit 1
fi
bin_dir=$1
pfile_dir=$2

data_files=(${ESPA_UNIT_TEST_DATA_DIR}/espa-surface-reflectance/lndsr_ref/*)
input_dir=$ESPA_UNIT_TEST_DATA_DIR/espa-surface-reflectance

rm -rf lndsr
mkdir lndsr && cd lndsr

ln -s $input_dir/input_l7/*.hdr .
ln -s $input_dir/input_l7/*.img .
ln -s $input_dir/lndcal_ref/*.hdr .
ln -s $input_dir/lndcal_ref/*.img .
cp $input_dir/lndcal_ref/*.xml .
chmod u+w *.xml

sed -e s%LEDAPS_AUX_DIR%${LEDAPS_AUX_DIR}% $pfile_dir/$pfile > pfile.local
$bin_dir/lndsr --pfile pfile.local
if [ $? -ne 0 ]; then
    echo "Error: lndsr processing failed."
    exit 1
fi

# Compare the results.
status=0
for i in "${data_files[@]}"; do
    base_name=`basename $i`

    echo "Comparing $base_name..."

    # For the XML file, ignore the records that vary from one run to the next.
    # If the file is an ASCII header file, use diff.  Otherwise, assume it's
    # a binary file, and use cmp to dump the octal differences.
    ext="${i##*.}"
    if [ "$ext" = "xml" ]; then
        sed -e 's%<production_date>.*<%<production_date><%' $i > tmp1.xml
        sed -e 's%<production_date>.*<%<production_date><%' $base_name > \
                                                                      tmp2.xml
        diff tmp1.xml tmp2.xml
        if [ $? -ne 0 ]; then
            echo "${base_name} differs from reference version."
            status=1
        fi
    elif [ "$ext" = "hdr" ]; then
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
        if [ "x$numdiffs" != "x0" ]; then
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
rm -r lndsr

echo "Test completed successfully."
exit 0
