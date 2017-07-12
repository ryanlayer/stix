#!/bin/bash

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest

STOP_ON_FAIL=0

STIX="../../bin/stix"
if [ ! -e "$STIX" ]; then
    echo "stix not found at $STIX"
    exit 1
fi

SQLITE3=`which sqlite3`

GIGGLE=`which giggle`
if [ -z "$GIGGLE" ]; then
    if [ -e "../../../giggle/bin/giggle" ];then
        GIGGLE="../../../giggle/bin/giggle"
    else
        echo "giggle not found in path or at ../../../giggle/bin/giggle"
        exit 1
    fi
fi

# Make the index
run make_giggle_index \
    $GIGGLE index \
        -i "../data/four_alt_sort/*gz" \
        -o ../data/four_alt_sort_b \
        -f -s
assert_exit_code 0
assert_in_stderr "Indexed 87491 intervals."

# Make PED db
run make_ped_db_no_col \
    $STIX \
        -i ../data/four_alt_sort_b \
        -p ../data/four.ped \
        -d ../data/four.ped.db 
assert_exit_code 1

run make_ped_db \
    $STIX \
        -i ../data/four_alt_sort_b \
        -p ../data/four.ped \
        -d ../data/four.ped.db \
        -c 6
assert_exit_code 0
if [ -n $SQLITE3 ];then
    assert_equal 2 $( $SQLITE3 ../data/four.ped.db "SELECT COUNT(*) FROM ped WHERE Population=='CEU'" )
    assert_equal 2 $( $SQLITE3 ../data/four.ped.db "SELECT COUNT(*) FROM ped WHERE Sex==1" )
fi

run query_four \
    $STIX \
        -i ../data/four_alt_sort_b \
        -d ../data/four.ped.db \
        -s 500 \
        -f ../data/1kg.four.13.14.vcf.gz 
assert_exit_code 0
assert_equal 3068 $( cat $STDOUT_FILE | grep -c -v "^#" )
assert_equal 3068 $( cat $STDOUT_FILE | grep -v "^#" | grep -c STIX_QUANT_DEPTHS )
assert_equal 3068 $( cat $STDOUT_FILE | grep -v "^#" | grep -c STIX_QUANTS )
assert_equal 3068 $( cat $STDOUT_FILE | grep -v "^#" | grep -c STIX_ZERO )
assert_equal 3068 $( cat $STDOUT_FILE | grep -v "^#" | grep -c STIX_ONE )
assert_equal 12 $( cat $STDOUT_FILE | grep -v "^#" | grep -c "STIX_ZERO=0;" )
assert_equal 19 $( cat $STDOUT_FILE | grep -v "^#" | grep -c "STIX_ZERO=1;" )
assert_equal 20 $( cat $STDOUT_FILE | grep -v "^#" | grep -c "STIX_ZERO=2;" )
assert_equal 70 $( cat $STDOUT_FILE | grep -v "^#" | grep -c "STIX_ZERO=3;" )
assert_equal 2947 $( cat $STDOUT_FILE | grep -v "^#" | grep -c "STIX_ZERO=4;" )
assert_equal 0 $( cat $STDOUT_FILE | grep -v "^#" | grep -c "STIX_SAMPLE_DEPTH" )

run query_four \
    $STIX \
        -i ../data/four_alt_sort_b \
        -d ../data/four.ped.db \
        -s 500 \
        -f ../data/1kg.four.13.14.vcf.gz \
        -v Alt_file
assert_exit_code 0
assert_equal 3068 $( cat $STDOUT_FILE | grep -c -v "^#" )
assert_equal 3068 $( cat $STDOUT_FILE | grep -v "^#" | grep -c STIX_QUANT_DEPTHS )
assert_equal 3068 $( cat $STDOUT_FILE | grep -v "^#" | grep -c STIX_QUANTS )
assert_equal 3068 $( cat $STDOUT_FILE | grep -v "^#" | grep -c STIX_ZERO )
assert_equal 3068 $( cat $STDOUT_FILE | grep -v "^#" | grep -c STIX_ONE )
assert_equal $( cat $STDOUT_FILE | grep -v "^#" | grep -c "STIX_ZERO=4;" ) \
             $( cat $STDOUT_FILE | grep -v "^#" | grep -c -v "STIX_SAMPLE_DEPTH" )
assert_equal $( cat $STDOUT_FILE | grep -v "^#" | grep -c -v "STIX_ZERO=4;" ) \
             $( cat $STDOUT_FILE | grep -v "^#" | grep -c "STIX_SAMPLE_DEPTH" )


$GIGGLE index \
    -i "../data/low_alt/*gz" \
    -o ../data/low_alt_b \
    -f -s \
> /dev/null 2> /dev/null

$STIX \
    -i ../data/low_alt_b \
    -p ../data/low.ped \
    -d ../data/low.ped.db \
    -c 2 \
> /dev/null 2> /dev/null

for S in NA12878 NA12889 NA12890; do
    run query_low_cover_del_$S \
        $STIX \
            -i ../data/low_alt_b \
            -d ../data/low.ped.db \
            -t DEL \
            -s 500 \
            -l 19:12694867-12694867 \
            -r 19:12698924-12698924 \
            -F "Sample=\"$S\""
    assert_exit_code 0
    BT=$( ((bedtools intersect \
                -b <(echo -e "19\t12694367\t12694867") \
                -a ../data/low_alt/$S.bed.gz) \
           | cut -f5-7 \
           | bedtools intersect \
                -a stdin \
                -b <(echo -e "19\t12698924\t12699424")
          ) | wc -l ) 
    assert_equal $BT $( cat $STDOUT_FILE | tail -n 1 | awk '{print $4}' )
done

for S in NA12878 NA12889 NA12890; do
    run query_low_cover_inv_$S \
        $STIX \
            -i ../data/low_alt_b \
            -d ../data/low.ped.db \
            -t INV \
            -s 500 \
            -l 12:12544868-12544868 \
            -r 12:12546613-12546613 \
            -F "Sample=\"$S\""
    assert_exit_code 0
    BT=$( ((bedtools intersect \
                -b <(echo -e "12\t12544368\t12545368") \
                -a ../data/low_alt/$S.bed.gz) \
           | cut -f5-7 \
           | bedtools intersect \
                -a stdin \
                -b <(echo -e "12\t12546113\t12547113")
          ) | wc -l ) 
    assert_equal $BT $( cat $STDOUT_FILE | tail -n 1 | awk '{print $4}' )
done
