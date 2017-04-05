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


