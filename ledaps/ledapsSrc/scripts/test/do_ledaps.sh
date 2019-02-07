#!/bin/bash
#
# Test the ledaps Python script.
#
# Arguments:
#   L1 angle band flag:   0 = don't use L1 angle bands (i.e., create angle
#                             bands)
#                         1 = use L1 angle bands
#  path to ledaps script

base_scene="LE07_L1TP_043028_20020419_20180206_01_T1"
pfile="lndsr.${base_scene}.txt"

if [ "$#" -ne 2 ]; then
    echo "Usage:  $0 <l1_angle_band_flag> <path_to_ledaps>"
    exit 1
fi
angle_band_flag=$1
bin_dir=$2

data_files=(${ESPA_UNIT_TEST_DATA_DIR}/espa-surface-reflectance/lndsr_ref/*)
input_dir=$ESPA_UNIT_TEST_DATA_DIR/espa-surface-reflectance/input_l7

rm -rf ledaps_$angle_band_flag
mkdir ledaps_$angle_band_flag && cd ledaps_$angle_band_flag

ln -s $input_dir/*.img .
ln -s $input_dir/*.hdr .
cp $input_dir/${base_scene}.xml .
chmod u+w ${base_scene}.xml

# If the test creates angle bands, remove the links to the L1 angle bands
# in the input data dir.  Also remove L1 angle bands from the input XML file.
if [ $1 -eq 0 ]; then
    angle_band_opt=""
    ln -s $input_dir/${base_scene}_ANG.txt .
    rm *solar* *sensor*
    sed_remove_angles='/product=.angle_bands/,/<.band>/ d'
else
    angle_band_opt="--use_l1_angle_bands"
    sed_remove_angles=''
fi

sed -e "$sed_remove_angles" \
    ${input_dir}/${base_scene}.xml > ${base_scene}.xml
$bin_dir/do_ledaps.py --xml=${base_scene}.xml $angle_band_opt
if [ $? -ne 0 ]; then
    echo "Error: ledaps processing failed."
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
        sed -e 's%<production_date>.*<%<production_date><%' \
            -e "$sed_remove_angles" \
            $i > tmp1.xml
        sed -e 's%<production_date>.*<%<production_date><%' \
            -e "$sed_remove_angles" \
            $base_name > tmp2.xml
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
rm -r ledaps_$angle_band_flag

echo "Test completed successfully."
exit 0
