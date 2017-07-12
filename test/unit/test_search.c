#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>

#include "unity.h"

#include <src/giggle_index.h>
#include <src/util.h>
#include <src/lists.h>

#include <htslib/hfile.h>
#include <htslib/vcf.h>

#include "ped.h"
#include "search.h"

void setUp(void) { }
void tearDown(void) { }

//13 37414451 37414451    13 37414779 37414779
//13 39934236 39934236    13 39935552 39935552
//13 51069347 51069347    13 51075082 51075082
//14 83852035 83852035    14 83852358 83852358
//13 67950542 67950542    13 67950594 67950594

//{{{void test_stix_region_to_breakpoint(void)
void test_stix_region_to_breakpoint(void)
{
    char *region = strdup("1:1-2");
    struct stix_breakpoint *bp = stix_region_to_breakpoint(region);
    TEST_ASSERT_EQUAL(0, strcmp("1", bp->chrm));
    TEST_ASSERT_EQUAL(1, bp->start);
    TEST_ASSERT_EQUAL(2, bp->end);
    TEST_ASSERT_EQUAL(0, bp->strand);

    free(region);
    free(bp->chrm);
    free(bp);

    region = strdup("X:10-20");
    bp = stix_region_to_breakpoint(region);
    TEST_ASSERT_EQUAL(0, strcmp("X", bp->chrm));
    TEST_ASSERT_EQUAL(10, bp->start);
    TEST_ASSERT_EQUAL(20, bp->end);
    TEST_ASSERT_EQUAL(0, bp->strand);

    free(region);
    free(bp->chrm);
    free(bp);
}
//}}}

void test_stix_check_inv(void)
{
    struct stix_breakpoint q_l,
                           q_r,
                           i_l,
                           i_r;
    uint32_t evidence_type;
    uint32_t slop = 10;



    /*
     *       i_l   i_r
     *        v     v
     *     q_l| | | |q_r
     * ------^ABCDEFG$-----
     *     +.......+       
     *       +....+
     *        -........-
     *          -....-
     *     +========-
     *        -========+
     *
     *     +..-   +..- 
     *       +..-   +..-      
     *     +==+     +==+
     * ------^GFEDCBA$-----
     */



    q_l.chrm = "1";
    q_l.start = 100;
    q_l.end = 101;
    q_l.strand = -1;

    q_r.chrm = "1";
    q_r.start = 500;
    q_r.end = 501;
    q_r.strand = 1;


    // Paired end
    //PASS
    /*                100|             |500
     *           q_l           q_r
     *             s|----|e      s|----|e
     *       ------------ABCDEFGHIJKLMNO--------------
     *                 +..............+
     *       ------------ONMLKJIHGFEDCBA--------------
     *                 +..-
     */
    i_l.chrm = "1"; i_l.start = 98; i_l.end = 99; i_l.strand = 1;
    i_r.chrm = "1"; i_r.start = 499; i_r.end = 500; i_r.strand = 1;
    TEST_ASSERT_EQUAL(1, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(1, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));
    //FAIL
    //        ------------ABCDEFGHIJKLMNO--------------
    //                  -..............-
    i_l.chrm = "1"; i_l.start = 98; i_l.end = 99; i_l.strand = -1;
    i_r.chrm = "1"; i_r.start = 499; i_r.end = 500; i_r.strand = -1;
    TEST_ASSERT_EQUAL(0, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));
    //        ------------ABCDEFGHIJKLMNO--------------
    //                  +...+
    i_l.chrm = "1"; i_l.start = 98; i_l.end = 99; i_l.strand = 1;
    i_r.chrm = "1"; i_r.start = 299; i_r.end = 300; i_r.strand = 1;
    TEST_ASSERT_EQUAL(0, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));
    //        ------------ABCDEFGHIJKLMNO--------------
    //            +....................+
    i_l.chrm = "1"; i_l.start = 9; i_l.end = 10; i_l.strand = 1;
    i_r.chrm = "1"; i_r.start = 499; i_r.end = 500; i_r.strand = 1;
    TEST_ASSERT_EQUAL(0, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));
    //        ------------ABCDEFGHIJKLMNO--------------
    //                  +....................+
    i_l.chrm = "1"; i_l.start = 98; i_l.end = 99; i_l.strand = 1;
    i_r.chrm = "1"; i_r.start = 599; i_r.end = 600; i_r.strand = 1;
    TEST_ASSERT_EQUAL(0, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));

    //        ------------ABCDEFGHIJKLMNO--------------
    //            +..........................+
    i_l.chrm = "1"; i_l.start = 9; i_l.end = 10; i_l.strand = 1;
    i_r.chrm = "1"; i_r.start = 599; i_r.end = 600; i_r.strand = 1;
    TEST_ASSERT_EQUAL(0, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));

    //PASS
    /*                100|             |500
     *                q_l           q_r
     *                  s|----|e      s|----|e
     *       ------------ABCDEFGHIJKLMNO--------------
     *                    -..............-
     *       ------------ONMLKJIHGFEDCBA--------------
     *                                +..-
     */
    i_l.chrm = "1"; i_l.start = 101; i_l.end = 102; i_l.strand = -1;
    i_r.chrm = "1"; i_r.start = 501; i_r.end = 502; i_r.strand = -1;
    TEST_ASSERT_EQUAL(1, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(1, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));
    //FAIL
    //        ------------ABCDEFGHIJKLMNO--------------
    //                     +..............+
    i_l.chrm = "1"; i_l.start = 101; i_l.end = 102; i_l.strand = 1;
    i_r.chrm = "1"; i_r.start = 501; i_r.end = 502; i_r.strand = 1;
    TEST_ASSERT_EQUAL(0, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));
    //        ------------ABCDEFGHIJKLMNO--------------
    //                                -...-
    i_l.chrm = "1"; i_l.start = 301; i_l.end = 302; i_l.strand = -1;
    i_r.chrm = "1"; i_r.start = 601; i_r.end = 602; i_r.strand = -1;
    TEST_ASSERT_EQUAL(0, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));
    //        ------------ABCDEFGHIJKLMNO--------------
    //                 -..................-
    i_l.chrm = "1"; i_l.start = 81; i_l.end = 82; i_l.strand = -1;
    i_r.chrm = "1"; i_r.start = 501; i_r.end = 502; i_r.strand = -1;
    TEST_ASSERT_EQUAL(0, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));
    //        ------------ABCDEFGHIJKLMNO--------------
    //                     -.....................-
    i_l.chrm = "1"; i_l.start = 101; i_l.end = 102; i_l.strand = -1;
    i_r.chrm = "1"; i_r.start = 601; i_r.end = 602; i_r.strand = -1;
    TEST_ASSERT_EQUAL(0, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));
    //        ------------ABCDEFGHIJKLMNO--------------
    //                 -.........................-
    i_l.chrm = "1"; i_l.start = 81; i_l.end = 82; i_l.strand = -1;
    i_r.chrm = "1"; i_r.start = 601; i_r.end = 602; i_r.strand = -1;
    TEST_ASSERT_EQUAL(0, stix_check_inv(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, INV));

    // Split read
}

//{{{void test_stix_check_dup(void)
void test_stix_check_dup(void)
{
    struct stix_breakpoint q_l,
                           q_r,
                           i_l,
                           i_r;
    uint32_t evidence_type;
    uint32_t slop = 10;



    /*
     *       i_l   i_r
     *        v     v
     *     q_l|-| |-|q_r
     * ------^ABCDEFG$-----
     *          -...+
     *        -...+
     *         +==+
     *
     *            +..- 
     *              +..- 
     *             +==+
     * ------^ABCDEFGABCDEFG$-----
     */

    i_l.chrm = "1";
    i_l.start = 100;
    i_l.end = 101;
    i_l.strand = -1;

    i_r.chrm = "1";
    i_r.start = 500;
    i_r.end = 501;
    i_r.strand = 1;

    q_l.chrm = "1";
    q_l.start = 100;
    q_l.end = 101;
    q_l.strand = -1;

    q_r.chrm = "1";
    q_r.start = 500;
    q_r.end = 501;
    q_r.strand = 1;


    // Paired end
    //PASS
    TEST_ASSERT_EQUAL(1, stix_check_dup(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(1, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DUP));
    //FAIL
    i_l.strand = 1;
    i_r.strand = -1;
    TEST_ASSERT_EQUAL(0, stix_check_dup(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DUP));
    i_l.strand = 1;
    i_r.strand = 1;
    TEST_ASSERT_EQUAL(0, stix_check_dup(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DUP));
    i_l.strand = 1;
    i_r.strand = 1;
    TEST_ASSERT_EQUAL(0, stix_check_dup(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DUP));

    i_l.start = 10;
    i_l.end = 11;
    i_l.strand = -1;
    i_r.start = 500;
    i_r.end = 501;
    i_r.strand = 1;
    TEST_ASSERT_EQUAL(0, stix_check_dup(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DUP));
    i_l.start = 100;
    i_l.end = 101;
    i_l.strand = -1;
    i_r.start = 5000;
    i_r.end = 5001;
    i_r.strand = 1;
    TEST_ASSERT_EQUAL(0, stix_check_dup(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DUP));
    i_l.start = 100;
    i_l.end = 101;
    i_l.strand = -1;
    i_r.start = 300;
    i_r.end = 301;
    i_r.strand = 1;
    TEST_ASSERT_EQUAL(0, stix_check_dup(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DUP));
    i_l.start = 200;
    i_l.end = 201;
    i_l.strand = -1;
    i_r.start = 500;
    i_r.end = 501;
    i_r.strand = 1;


    //split read
    TEST_ASSERT_EQUAL(0, stix_check_dup(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DUP));
    i_l.start = 100;
    i_l.end = 101;
    i_l.strand = 1;
    i_r.strand = 1;
    TEST_ASSERT_EQUAL(1, stix_check_dup(&q_l, &q_r, &i_l, &i_r, 1, slop));
    TEST_ASSERT_EQUAL(1, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 1, slop, DUP));
    i_l.strand = -1;
    i_r.strand = -1;
    TEST_ASSERT_EQUAL(1, stix_check_dup(&q_l, &q_r, &i_l, &i_r, 1, slop));
    TEST_ASSERT_EQUAL(1, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 1, slop, DUP));

}
//}}}

//{{{void test_stix_check_del(void)
void test_stix_check_del(void)
{
    struct stix_breakpoint q_l,
                           q_r,
                           i_l,
                           i_r;
    uint32_t evidence_type;
    uint32_t slop = 100;


    /*   
     *      i_l     i_r
     *       v       v
     * q_l | |       | | q_r
     * ------^ABCDEFG$-----
     *       +.........- 
     *     +.........- 
     *
     *       +==- 
     *     +==- 
     * ------^$-------
     */

    q_l.chrm = "1";
    q_l.start = 100;
    q_l.end = 110;
    q_l.strand = 1;

    q_r.chrm = "1";
    q_r.start = 500;
    q_r.end = 510;
    q_r.strand = -1;

    i_l.chrm = "1";
    i_l.start = 100;
    i_l.end = 110;
    i_l.strand = 1;

    i_r.chrm = "1";
    i_r.start = 500;
    i_r.end = 510;
    i_r.strand = -1;

    TEST_ASSERT_EQUAL(1, stix_check_del(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(1, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DEL));
    TEST_ASSERT_EQUAL(0, stix_check_del(&q_l, &q_r, &i_l, &i_r, 1, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 1, slop, DEL));

    i_r.start = 1500;
    i_r.end = 1510;
    TEST_ASSERT_EQUAL(0, stix_check_del(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DEL));
    TEST_ASSERT_EQUAL(0, stix_check_del(&q_l, &q_r, &i_l, &i_r, 1, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 1, slop, DEL));

    i_r.start = 50;
    i_r.end = 51;
    TEST_ASSERT_EQUAL(0, stix_check_del(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DEL));
    TEST_ASSERT_EQUAL(0, stix_check_del(&q_l, &q_r, &i_l, &i_r, 1, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 1, slop, DEL));

    i_l.chrm = "1";
    i_r.chrm = "2";
    TEST_ASSERT_EQUAL(0, stix_check_del(&q_l, &q_r, &i_l, &i_r, 0, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 0, slop, DEL));
    TEST_ASSERT_EQUAL(0, stix_check_del(&q_l, &q_r, &i_l, &i_r, 1, slop));
    TEST_ASSERT_EQUAL(0, stix_check_sv(&q_l, &q_r, &i_l, &i_r, 1, slop, DEL));
}
//}}}

//{{{ void test_stix_parse_result(void)
void test_stix_parse_result(void)
{

    struct stix_breakpoint *left = NULL, *right = NULL;

    char *result = NULL;
    asprintf(&result, "1\t2\t3\t1\t4\t5\t6\t-1\t0");
    uint32_t evidence_type;

    uint32_t ret = stix_parse_result(result, &left, &right, &evidence_type);

    TEST_ASSERT_EQUAL(0, strcmp("1", left->chrm));
    TEST_ASSERT_EQUAL(2, left->start);
    TEST_ASSERT_EQUAL(3, left->end);
    TEST_ASSERT_EQUAL(1, left->strand);

    TEST_ASSERT_EQUAL(0, strcmp("4", right->chrm));
    TEST_ASSERT_EQUAL(5, right->start);
    TEST_ASSERT_EQUAL(6, right->end);
    TEST_ASSERT_EQUAL(-1, right->strand);

    TEST_ASSERT_EQUAL(0, evidence_type);
}
//}}}

//{{{ void test_stix_run_giggle_query(void)
void test_stix_run_giggle_query(void)
{
    char *input_path_name = "../data/four_alt_sort/*gz";
    char *output_path_name = "tmp_test_stix_run_giggle_query";

    uint64_t indexed_intervals = giggle_bulk_insert(input_path_name,
                                                    output_path_name,
                                                    1);

    struct stix_breakpoint left = {"13", 37414451, 37414452, 0}, 
                           right = {"13", 37414779, 37414780, 0}; 

    struct giggle_index *gi = NULL;
    struct uint_pair *sample_alt_depths = NULL;

    //13 37414451 37414451    13 37414779 37414779
    uint32_t num_sample_alt_depths = 
            stix_run_giggle_query(&gi,
                                  output_path_name,
                                  DEL,
                                  &left,
                                  &right,
                                  500,
                                  NULL,
                                  0,
                                  &sample_alt_depths);
    /*
      echo -e "13\t37413951\t37414452" \
      | bedtools intersect \
           -a ../data/four_alt_sort/HG00672.13.14.bed.gz \
           -b stdin \
      | wc -l
      5
     
      echo -e "13\t37413951\t37414452" \
      | bedtools intersect \
           -a ../data/four_alt_sort/HG00674.13.14.bed.gz \
           -b stdin \
      | wc -l
      5
     
      echo -e "13\t37413951\t37414452" \
      | bedtools intersect \
           -a ../data/four_alt_sort/NA12812.13.14.bed.gz \
           -b stdin \
      | wc -l
      1
     
      echo -e "13\t37413951\t37414452" \
      | bedtools intersect \
           -a ../data/four_alt_sort/NA12878.13.14.bed.gz \
           -b stdin \
      | wc -l
      6
     */

    TEST_ASSERT_TRUE(gi != NULL);
    TEST_ASSERT_EQUAL(4, num_sample_alt_depths);

    uint32_t A_sample_alt_depths[4] = {5,5,1,6};

    uint32_t i;
    for (i = 0; i < num_sample_alt_depths; ++i) {
        TEST_ASSERT_EQUAL(A_sample_alt_depths[i], sample_alt_depths[i].first);
        TEST_ASSERT_EQUAL(0, sample_alt_depths[i].second);
    }

    /*
      echo -e "13\t39929236\t39934236" \
      | bedtools intersect \
           -a ../data/four_alt_sort/HG00672.13.14.bed.gz \
           -b stdin \
      | awk '$5==13' \
      | wc -l
      5
     
      echo -e "13\t39929236\t39934236" \
      | bedtools intersect \
           -a ../data/four_alt_sort/HG00674.13.14.bed.gz \
           -b stdin \
      | awk '$5==13' \
      | wc -l
      12      
     
      echo -e "13\t39929236\t39934236" \
      | bedtools intersect \
           -a ../data/four_alt_sort/NA12812.13.14.bed.gz \
           -b stdin \
      | awk '$5==13' \
      | wc -l
      4
     
      echo -e "13\t39929236\t39934236" \
      | bedtools intersect \
           -a ../data/four_alt_sort/NA12878.13.14.bed.gz \
           -b stdin \
      | awk '$5==13' \
      | wc -l
      0
     */
    left.start = 39934236;
    left.end = 39934236;
    right.start = 39935552;
    right.end = 39935552;
    num_sample_alt_depths = stix_run_giggle_query(&gi,
                                                  output_path_name,
                                                  DEL,
                                                  &left,
                                                  &right,
                                                  500,
                                                  NULL,
                                                  0,
                                                  &sample_alt_depths);
    TEST_ASSERT_EQUAL(4, num_sample_alt_depths);
    uint32_t A_sample_alt_depths_1[4] = {5,12,4,0};
    for (i = 0; i < num_sample_alt_depths; ++i) {
        TEST_ASSERT_EQUAL(A_sample_alt_depths_1[i], sample_alt_depths[i].first);
        TEST_ASSERT_EQUAL(0, sample_alt_depths[i].second);
    }

    free(sample_alt_depths);
    sample_alt_depths = NULL;

    uint32_t test_sub_group[2] = {1,3};
    num_sample_alt_depths = stix_run_giggle_query(&gi,
                                                  output_path_name,
                                                  DEL,
                                                  &left,
                                                  &right,
                                                  500,
                                                  test_sub_group,
                                                  2,
                                                  &sample_alt_depths);
    TEST_ASSERT_EQUAL(2, num_sample_alt_depths);
    uint32_t A_sample_alt_depths_2[2] = {12,0};
    for (i = 0; i < num_sample_alt_depths; ++i) {
        TEST_ASSERT_EQUAL(A_sample_alt_depths_2[i], sample_alt_depths[i].first);
        TEST_ASSERT_EQUAL(0, sample_alt_depths[i].second);
    }

    giggle_index_destroy(&gi);

    rmrf(output_path_name);
}
//}}}

//{{{void test_ped_get_uniq_col_groups_stix_run_giggle_query(void)
void test_ped_get_uniq_col_groups_stix_run_giggle_query(void)
{

    char *input_path_name = "../data/four_alt_sort/*gz";
    char *output_path_name = "tmp_test_ped_create_db";

    uint64_t indexed_intervals = giggle_bulk_insert(input_path_name,
                                                    output_path_name,
                                                    1);

    char *ped_file_name = "../data/four.ped";
    char *ped_file_name_db = "four.ped.db";

    uint32_t num_rows = ped_create_db(ped_file_name,
                                      ped_file_name_db,
                                      output_path_name,
                                      6);



    char *cols_2[2] = {"Sex", "Population"};
    char ***uniq_vals;
    uint32_t num_uniq_vals;
    uint32_t **uniq_groups_ids;
    uint32_t *uniq_groups_sizes;
    sqlite3 *db = NULL;
    num_uniq_vals = ped_get_uniq_col_groups(ped_file_name_db,
                                            &db,
                                            cols_2,
                                            2,
                                            NULL,
                                            &uniq_vals,
                                            &uniq_groups_ids,
                                            &uniq_groups_sizes);
    /* 
     sqlite3 four.ped.db \
        "SELECT Sex,Population FROM ped GROUP BY Sex,Population"
     1|CEU
     1|CHS
     2|CEU
     2|CHS
     */
    TEST_ASSERT_EQUAL(0, strcmp("1",   uniq_vals[0][0]));
    TEST_ASSERT_EQUAL(0, strcmp("CEU", uniq_vals[0][1]));
    TEST_ASSERT_EQUAL(0, strcmp("1",   uniq_vals[1][0]));
    TEST_ASSERT_EQUAL(0, strcmp("CHS", uniq_vals[1][1]));
    TEST_ASSERT_EQUAL(0, strcmp("2",   uniq_vals[2][0]));
    TEST_ASSERT_EQUAL(0, strcmp("CEU", uniq_vals[2][1]));
    TEST_ASSERT_EQUAL(0, strcmp("2",   uniq_vals[3][0]));
    TEST_ASSERT_EQUAL(0, strcmp("CHS", uniq_vals[3][1]));

    /*
     sqlite3 four.ped.db \
        "SELECT Giggle_File_Id,Alt_File FROM ped WHERE Sex=='1' AND Population=='CEU'"
     2|NA12812.13.14.bed.gz
     sqlite3 four.ped.db \
        "SELECT Giggle_File_Id,Alt_File FROM ped WHERE Sex=='1' AND Population=='CHS'"
     1|HG00674.13.14.bed.gz
     sqlite3 four.ped.db \
        "SELECT Giggle_File_Id,Alt_File FROM ped WHERE Sex=='2' AND Population=='CEU'"
     3|NA12878.13.14.bed.gz
     sqlite3 four.ped.db \
        "SELECT Giggle_File_Id,Alt_File FROM ped WHERE Sex=='2' AND Population=='CHS'"
     0|HG00672.13.14.bed.gz
     */
    TEST_ASSERT_EQUAL(1, uniq_groups_sizes[0]);
    TEST_ASSERT_EQUAL(2, uniq_groups_ids[0][0]);

    TEST_ASSERT_EQUAL(1, uniq_groups_sizes[1]);
    TEST_ASSERT_EQUAL(1, uniq_groups_ids[1][0]);

    TEST_ASSERT_EQUAL(1, uniq_groups_sizes[2]);
    TEST_ASSERT_EQUAL(3, uniq_groups_ids[2][0]);

    TEST_ASSERT_EQUAL(1, uniq_groups_sizes[3]);
    TEST_ASSERT_EQUAL(0, uniq_groups_ids[3][0]);


    uint32_t *union_group_ids;
    uint32_t num_union_group_ids = 
            ped_union_groups(num_uniq_vals,
                             uniq_groups_ids,
                             uniq_groups_sizes,
                             &union_group_ids);

    TEST_ASSERT_EQUAL(4, num_union_group_ids);
    TEST_ASSERT_EQUAL(2, union_group_ids[0]);
    TEST_ASSERT_EQUAL(1, union_group_ids[1]);
    TEST_ASSERT_EQUAL(3, union_group_ids[2]);
    TEST_ASSERT_EQUAL(0, union_group_ids[3]);

    struct giggle_index *gi = NULL;
    struct uint_pair *sample_alt_depths = NULL;
    struct stix_breakpoint left = {"13", 37414451, 37414452, 0}, 
                           right = {"13", 37414779, 37414780, 0}; 

    uint32_t num_sample_alt_depths =
            stix_run_giggle_query(&gi,
                                  output_path_name,
                                  DEL,
                                  &left,
                                  &right,
                                  500,
                                  union_group_ids,
                                  num_union_group_ids,
                                  &sample_alt_depths);

    /*
      echo -e "13\t37413951\t37414452" \
      | bedtools intersect \
           -a ../data/four_alt_sort/NA12812.13.14.bed.gz \
           -b stdin \
      | wc -l
      1

      echo -e "13\t37413951\t37414452" \
      | bedtools intersect \
           -a ../data/four_alt_sort/HG00674.13.14.bed.gz \
           -b stdin \
      | wc -l
      5
   
      echo -e "13\t37413951\t37414452" \
      | bedtools intersect \
           -a ../data/four_alt_sort/NA12878.13.14.bed.gz \
           -b stdin \
      | wc -l
      6

      echo -e "13\t37413951\t37414452" \
      | bedtools intersect \
           -a ../data/four_alt_sort/HG00672.13.14.bed.gz \
           -b stdin \
      | wc -l
      5
     */
    uint32_t A_sample_alt_depths[4] = {1,5,6,5};

    uint32_t i;
    for (i = 0; i < num_sample_alt_depths; ++i) {
        TEST_ASSERT_EQUAL(A_sample_alt_depths[i], sample_alt_depths[i].first);
        TEST_ASSERT_EQUAL(0, sample_alt_depths[i].second);
    }

    sqlite3_close(db);
    giggle_index_destroy(&gi);
    free(union_group_ids);
    remove(ped_file_name_db);
    rmrf(output_path_name);
}
//}}}

//{{{void test_stix_get_uniq(void)
void test_stix_get_uniq(void)
{
    uint32_t s1[15] = {5, 5, 5, 5, 5, 1, 2, 2, 3, 3, 3, 4, 4, 4, 6};
    uint32_t A_u1[6] = {1, 2, 3, 4, 5, 6};
    uint32_t *u1;
    uint32_t num_uniq = stix_get_uniq(s1, 15, &u1);

    TEST_ASSERT_EQUAL(6, num_uniq);
    uint32_t i;
    for (i = 0; i < num_uniq; ++i)
        TEST_ASSERT_EQUAL(A_u1[i], u1[i]);

    free(u1);

    uint32_t num_rand = 10000;
    srand(1);
    uint32_t *r = (uint32_t  *)malloc(num_rand*sizeof(uint32_t));
    for (i = 0; i < num_rand; ++i) {
        r[i] = rand();
    }

    num_uniq = stix_get_uniq(r, num_rand, &u1);
    TEST_ASSERT_EQUAL(num_rand, num_uniq);

    free(u1);

    uint32_t uniqs = 0;
    for (i = 0; i < num_rand; ) {
        uint32_t j;
        uint32_t r_i = rand();
        uniqs += 1;
        for (j = 0; j < rand()%5 + 1; ++j) {
            r[i] = r_i;
            i+=1;
            if (i == num_rand)
                break;
        }
    }

    num_uniq = stix_get_uniq(r, num_rand, &u1);
    TEST_ASSERT_EQUAL(uniqs, num_uniq);
}
//}}}

//{{{ void test_stix_get_quartile_counts(void)
void test_stix_get_quartile_counts(void)
{
    
    uint32_t full[11] = {6, 7, 15, 36, 39, 40, 41, 42, 43, 47, 49};
    uint32_t Q1, Q2, Q3; 
    int32_t counts[4];


    uint32_t ret = stix_get_quartile_counts(full,
                                            11,
                                            &Q1, 
                                            &Q2, 
                                            &Q3, 
                                            counts);

    TEST_ASSERT_EQUAL(15, Q1); 
    TEST_ASSERT_EQUAL(40, Q2);
    TEST_ASSERT_EQUAL(43, Q3);
    TEST_ASSERT_EQUAL(2, counts[0]);
    TEST_ASSERT_EQUAL(3, counts[1]);
    TEST_ASSERT_EQUAL(3, counts[2]);
    TEST_ASSERT_EQUAL(3, counts[3]);
    TEST_ASSERT_EQUAL(11, counts[0]+counts[1]+counts[2]+counts[3]);
}
//}}}

//{{{void test_stix_get_vcf_breakpoints(void)
void test_stix_get_vcf_breakpoints(void)
{
    char *vcf_file_name =
        "../data/ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.13.vcf.gz";

    struct stix_breakpoint *left = NULL, *right = NULL;

    htsFile *fp = hts_open(vcf_file_name, "r");
    if ( !fp )
        err(EX_DATAERR, "Could not read file: %s", vcf_file_name);

    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr)
        errx(EX_DATAERR, "Header not found: %s\n", vcf_file_name);

    bcf1_t *line = bcf_init1();

    left = (struct stix_breakpoint *)
            malloc(sizeof(struct stix_breakpoint));
    left->chrm = NULL;

    right = (struct stix_breakpoint *)
            malloc(sizeof(struct stix_breakpoint));
    right->chrm = NULL;
    
    //{{{uint32_t A_pos[] = {
    uint32_t A_pos[1366] = {
        19137951, 19201274, 19327680, 19332586, 19343199, 19409662, 19431519,
        19517515, 19538091, 19568038, 19657637, 19748081, 19750981, 19751677,
        19769450, 19868503, 19916784, 20078039, 20193717, 20198191, 20249620,
        20342595, 20404660, 20432571, 20459958, 20466889, 20521072, 20556376,
        20650106, 20690063, 20703947, 20725019, 20825502, 20901598, 20904283,
        20974950, 20984878, 21122508, 21124682, 21131732, 21192264, 21253021,
        21455655, 21478092, 21486094, 21494548, 21521082, 21523861, 21531959,
        21676154, 21680348, 21684132, 21699973, 21729270, 21742331, 21742530,
        22052127, 22084932, 22199021, 22281752, 22335544, 22418080, 22418824,
        22422759, 22423929, 22526340, 22557059, 22574676, 22649437, 22772926,
        22815546, 22929101, 23124062, 23185045, 23285354, 23321660, 23325587,
        23386576, 23415573, 23427512, 23559430, 23632424, 23759334, 23791326,
        23817415, 23825664, 24029787, 24103778, 24224347, 24292906, 24353464,
        24426300, 24448439, 24508855, 24533332, 24593476, 24594438, 24631697,
        24646109, 24813663, 24842640, 24848738, 24910253, 24915267, 24930061,
        25000966, 25006832, 25024743, 25033844, 25157128, 25224222, 25229177,
        25248005, 25380706, 25452842, 25569778, 25570291, 25604232, 25608061,
        25633732, 25675323, 25928898, 26024762, 26040258, 26068178, 26083979,
        26176913, 26202676, 26290786, 26449350, 26770657, 26771638, 26814426,
        26833612, 27050020, 27058865, 27084595, 27146638, 27419654, 27466389,
        27789413, 27799101, 28069522, 28268151, 28531819, 28616137, 28654503,
        28663492, 28664302, 28669631, 28693868, 29133938, 29293532, 29343210,
        29479325, 29549817, 29602204, 29634964, 29750820, 29874696, 29915331,
        30127542, 30305730, 30427112, 30436507, 30617053, 30620964, 30662842,
        30745625, 30781681, 30955694, 30962247, 30985676, 31043601, 31196962,
        31238278, 31418158, 31510657, 31529016, 31564262, 31580834, 31586529,
        31630207, 31652399, 31703504, 31943658, 32026546, 32151258, 32152169,
        32155694, 32261571, 32325271, 32473786, 32475166, 32477801, 32478934,
        32532611, 32571771, 32671727, 32748396, 33130540, 33138320, 33138833,
        33177068, 33207476, 33390940, 33396504, 33444606, 33449762, 33575754,
        33932508, 34132580, 34132787, 34135729, 34192594, 34445402, 34641611,
        34824085, 34955842, 35024653, 35079029, 35138379, 35175772, 35355315,
        35459861, 35558633, 35564214, 35629445, 35759866, 35917527, 35989695,
        36004497, 36018573, 36031960, 36266865, 36470881, 36621733, 36629832,
        36798633, 36800966, 36806770, 36818914, 36853343, 36879180, 36992488,
        37245844, 37414451, 37434425, 37474506, 37497329, 37866506, 37871695,
        38048420, 38054350, 38070390, 38071941, 38282648, 38284604, 38388850,
        38411352, 38413103, 38425303, 38469767, 38476656, 38557766, 38683184,
        38701628, 39057317, 39108149, 39119964, 39154253, 39379047, 39384022,
        39528210, 39747518, 39803228, 39934236, 39949745, 39949861, 40058283,
        40058532, 40118735, 40193625, 40432610, 40454957, 40615339, 40739300,
        40856674, 40987879, 41076764, 41083165, 41103707, 41126533, 41235679,
        41270467, 41336593, 41342478, 41343494, 41645337, 41848920, 41850132,
        41862020, 41964015, 41990987, 41995088, 42003505, 42103909, 42393271,
        42518974, 42590782, 42642364, 43053246, 43203341, 43300174, 43324654,
        43505477, 43537316, 43617959, 43632181, 43649646, 43730050, 43801077,
        43820808, 43855567, 44157706, 44164101, 44514711, 44569426, 44680114,
        44817784, 44828572, 44920443, 45042429, 45043557, 45089908, 45120082,
        45163226, 45208067, 45278784, 45373416, 45396007, 45407777, 45407920,
        45597154, 45755053, 45793595, 45867926, 45926499, 45942983, 45951543,
        45965244, 45987925, 46115892, 46187271, 46252087, 46272789, 46292265,
        46329517, 46336590, 46336979, 46357796, 46605221, 46620528, 46628565,
        46634503, 46694688, 46725129, 46773603, 46791846, 46803688, 46851547,
        46853391, 46897499, 46914555, 46959924, 47040677, 47116976, 47184714,
        47281287, 47450372, 47499171, 47517876, 47754499, 47813923, 47836225,
        47941349, 47979453, 48007007, 48015694, 48026449, 48131313, 48271262,
        48376163, 48432866, 48446585, 48491687, 48664569, 48778665, 48796677,
        48881747, 49117576, 49129111, 49141939, 49153099, 49202718, 49308956,
        49533344, 49830324, 49888567, 49951399, 49975026, 50040827, 50195239,
        50383240, 50420122, 50474431, 50474496, 50485366, 50513107, 50579077,
        50737594, 50816917, 50954424, 50976406, 51001531, 51031205, 51069347,
        51264793, 51444745, 51452890, 51512934, 51634556, 52426982, 52936213,
        53035198, 53037516, 53141039, 53159492, 53175570, 53337761, 53357909,
        53404126, 53491628, 53499301, 53509531, 53522252, 53555507, 53618359,
        53908426, 54106954, 54234765, 54239148, 54254732, 54357929, 54423928,
        54691869, 54726311, 54786024, 54899481, 54963912, 55168991, 55253911,
        55262417, 55622332, 55628383, 55652760, 55660490, 55725083, 55776646,
        55872248, 55898609, 55916655, 55941166, 55953583, 55970790, 55972109,
        56151096, 56158689, 56241930, 56244426, 56267802, 56366307, 56449028,
        56704080, 56824061, 56851591, 56887566, 56913507, 57034959, 57201021,
        57234396, 57251122, 57309233, 57372476, 57395823, 57405047, 57445368,
        57550612, 57617322, 57656897, 57706322, 57715850, 57755853, 57768513,
        57787093, 57812347, 57905549, 57908967, 57944699, 57971829, 57977661,
        58047341, 58049571, 58053922, 58125856, 58141970, 58170351, 58230122,
        58364827, 58369964, 58386042, 58429647, 58450881, 58501696, 58503079,
        58589036, 58709299, 58713272, 58724187, 58769855, 58846616, 58910776,
        58946874, 59021863, 59173493, 59227220, 59346064, 59371397, 59372231,
        59468049, 59557464, 59706151, 59732387, 59793967, 59845024, 59856817,
        59893265, 59920657, 60195594, 60441789, 60441949, 60881444, 61210105,
        61211424, 61234615, 61470068, 61505666, 61526951, 61551693, 61552987,
        61706364, 61711881, 61955908, 61972560, 62012122, 62015482, 62145338,
        62188977, 62199159, 62410977, 62528684, 62654462, 62701484, 62778367,
        62807244, 62901975, 62913811, 63102170, 63221783, 63232491, 63283529,
        63312341, 63316943, 63377012, 63391261, 63527431, 63590763, 63654332,
        63669367, 63692671, 63732368, 63734729, 63813301, 63861850, 63887835,
        63945862, 63973697, 64009797, 64092323, 64224634, 64259306, 64311708,
        64313714, 64412737, 64413115, 64415383, 64663178, 64701555, 64708313,
        64747544, 64809684, 64816746, 64880332, 64918629, 64960282, 65065553,
        65185508, 65259027, 65269221, 65342437, 65352157, 65383658, 65443709,
        65525844, 65553795, 65611730, 65627743, 65637801, 65672871, 65751249,
        65760188, 65862105, 66061578, 66192707, 66322948, 66454066, 66583084,
        66715451, 66827310, 66891604, 66893461, 66963964, 66985485, 67098389,
        67133809, 67140826, 67174780, 67219662, 67225113, 67239841, 67278789,
        67336014, 67354887, 67455650, 67512194, 67774476, 67844490, 67861804,
        67947821, 67950542, 67979923, 68046883, 68090914, 68100521, 68135372,
        68148293, 68162195, 68280917, 68361106, 68375480, 68397289, 68412126,
        68512083, 68538974, 68671779, 68834333, 68842989, 68876787, 68971851,
        68974475, 68977877, 69061733, 69113265, 69199224, 69244692, 69385826,
        69390582, 69402666, 69446730, 69451996, 69452083, 69487563, 69699225,
        69737856, 69762469, 69762610, 69837680, 69936959, 69945225, 69992080,
        70054738, 70079642, 70125481, 70196108, 70401965, 70559014, 70643810,
        70646032, 70735724, 70747968, 70781285, 70782022, 70806748, 70808717,
        70832742, 70834619, 70923321, 70931613, 71009862, 71010008, 71069243,
        71069516, 71123947, 71169278, 71187931, 71195662, 71197723, 71249688,
        71274807, 71281394, 71382622, 71416116, 71446636, 71677047, 71786943,
        71841988, 71865440, 71876312, 71912648, 71917480, 71924371, 71994738,
        72026449, 72353291, 72353387, 72386637, 72453942, 72477264, 72478323,
        72522618, 72538150, 72708034, 72774393, 72843031, 72886169, 72893591,
        73014489, 73299589, 73418866, 73721912, 73847728, 73948540, 74005683,
        74088212, 74357216, 74479408, 74714807, 74956939, 74958044, 74987886,
        74987942, 75126449, 75436056, 75437624, 75439503, 75501355, 75603674,
        75625199, 75629554, 75678825, 75809234, 75855728, 75870962, 76086980,
        76190024, 76363073, 76386568, 76468948, 76475002, 76475380, 76563144,
        76589467, 76608863, 76656282, 76657932, 76680448, 76771729, 76825484,
        76868478, 76892786, 76963996, 76993254, 76993592, 77001712, 77065606,
        77116818, 77965821, 78116716, 78171272, 78234266, 78387406, 78391438,
        78455964, 78503659, 78507580, 78556446, 78764971, 78806448, 78852265,
        78936926, 78966861, 78997598, 79203003, 79245919, 79463205, 79487066,
        79547696, 79633826, 79757405, 79882783, 79998433, 80454912, 80463449,
        80595659, 80636620, 80651619, 80681809, 80750200, 80769215, 80773097,
        80882593, 80962476, 80998713, 81027703, 81032178, 81035751, 81037160,
        81280265, 81322565, 81403339, 81486103, 81738623, 81765330, 81805787,
        81835856, 81857908, 81867481, 81971382, 82083860, 82169385, 82208815,
        82264637, 82293556, 82355641, 82366876, 82462079, 82510664, 82526320,
        82537560, 82585179, 82692429, 82693737, 82756456, 82759262, 82901862,
        83049480, 83065129, 83131132, 83165934, 83209741, 83393905, 83417949,
        83439655, 83557604, 83639070, 83649455, 83709206, 83739497, 83788210,
        83801227, 83813363, 83829668, 83841962, 83890389, 83909763, 83921689,
        83946149, 84035973, 84036242, 84050412, 84101668, 84102127, 84127645,
        84176355, 84184868, 84209439, 84258707, 84282816, 84294973, 84326054,
        84466620, 84482284, 84488058, 84528154, 84532224, 84543906, 84564298,
        84590695, 84628249, 84628527, 84667990, 84902320, 85045318, 85045669,
        85048214, 85054244, 85055132, 85056636, 85105644, 85182827, 85185373,
        85198016, 85223246, 85270032, 85275060, 85363388, 85398548, 85414626,
        85472210, 85547745, 85662155, 85712619, 85800252, 85845094, 85874699,
        85929515, 85936698, 86056785, 86264810, 86385463, 86507908, 86535639,
        86579415, 86635010, 86652256, 86652511, 86659202, 86673791, 86674318,
        86727177, 86975520, 86985094, 87171110, 87213429, 87250630, 87325078,
        87326003, 87369477, 87378200, 87383743, 87396748, 87428591, 87451419,
        87546118, 87604541, 87641809, 87683395, 87780600, 87826217, 87831226,
        87891387, 87943918, 87956012, 87973360, 87987546, 88002167, 88022747,
        88277139, 88398093, 88469260, 88518614, 88632584, 88662182, 88956615,
        88963898, 89015296, 89149180, 89180086, 89235681, 89306848, 89410134,
        89421966, 89738894, 89846868, 89926170, 90002158, 90011334, 90205547,
        90205641, 90210860, 90248499, 90282175, 90348126, 90353472, 90427195,
        90536709, 90642915, 90682364, 90685691, 90726943, 90771710, 90780145,
        90785641, 90862850, 90910915, 91024120, 91046169, 91104395, 91211719,
        91306922, 91316300, 91328908, 91475871, 91518995, 91565667, 91588665,
        91609391, 91623829, 91717987, 91730945, 91754060, 91901930, 92111709,
        92179524, 92386317, 92419045, 92493008, 92520182, 92532464, 92550570,
        92641017, 92691606, 92711494, 92810438, 92857572, 92916846, 93000101,
        93007985, 93009206, 93056849, 93075342, 93126579, 93244533, 93462398,
        93541170, 93703457, 93809360, 93825443, 93883567, 93893263, 94006505,
        94011781, 94040798, 94059533, 94087904, 94090465, 94115322, 94142615,
        94229965, 94255922, 94267869, 94392427, 94512191, 94700666, 94769505,
        94829588, 94868594, 95003326, 95070439, 95196983, 95290106, 95438745,
        95568488, 95639513, 95716495, 95776503, 95781155, 95837866, 95887798,
        95893200, 95969028, 95979815, 95981746, 95984375, 96012524, 96020121,
        96061433, 96215514, 96417941, 96788630, 96800745, 96859122, 96869270,
        96919749, 97003590, 97041104, 97050602, 97085587, 97259860, 97271197,
        97304546, 97338253, 97345755, 97347834, 97490423, 97533065, 97607763,
        97627382, 97640600, 97895054, 97941665, 98047204, 98103788, 98147612,
        98161206, 98255649, 98409007, 98409905, 98463178, 98529935, 98600665,
        98896542, 98898311, 99067183, 99094072, 99128485, 99151177, 99235813,
        99254021, 99256683, 99279600, 99300107, 99313344, 99545457, 99631023,
        99686987, 99789002, 99940467, 100006835, 100081323, 100112931,
        100225836, 100493743, 100573019, 100605170, 100704290, 100734298,
        100737563, 100755889, 100775567, 100778540, 100821968, 100854103,
        100899294, 101073711, 101205430, 101209219, 101396190, 101406453,
        101414621, 101456032, 101468183, 101608612, 101625773, 101768518,
        101894135, 101915002, 101990710, 102038992, 102121830, 102176440,
        102282155, 102421141, 102440059, 102562811, 102635422, 102649543,
        102658009, 102858878, 102928662, 102937997, 102952659, 103097005,
        103128050, 103357637, 103500932, 103516186, 103599027, 103687334,
        103711721, 103711779, 103712336, 103875493, 104015378, 104026134,
        104034356, 104041403, 104091776, 104166006, 104276663, 104365381,
        104365533, 104527993, 104606440, 104672818, 104780402, 104807531,
        104883059, 104983248, 104998186, 105097588, 105107567, 105131439,
        105195941, 105235750, 105350008, 105350325, 105428174, 105528926,
        105753209, 105793861, 105809814, 105953655, 106109511, 106191651,
        106225759, 106226828, 106350892, 106364012, 106449327, 106475334,
        106660671, 106706335, 106707033, 106957625, 107012250, 107033444,
        107086233, 107162937, 107206205, 107427685, 107437827, 107453721,
        107510021, 107520787, 107700908, 107711721, 107835832, 108115980,
        108191306, 108298387, 108444027, 108451626, 108560567, 108603779,
        108773920, 108849020, 108886626, 108891563, 108910536, 108947031,
        108999898, 109054242, 109093036, 109267189, 109361752, 109419916,
        109476778, 109508305, 109648933, 109677382, 109679682, 109679772,
        109679876, 109782622, 109794696, 109883290, 110197589, 110221791,
        110228812, 110273747, 110321555, 110359533, 110359650, 110381571,
        110451516, 110461026, 110473028, 110732474, 110741110, 110798244,
        110969682, 111125559, 111136027, 111410060, 111432999, 111587190,
        111601037, 111628194, 111703587, 111716033, 111748379, 111769191,
        112044488, 112154522, 112321839, 112340112, 112590743, 112624416,
        112648970, 112690980, 112828690, 113263252, 113281665, 113324987,
        113419734, 113454412, 113481313, 113492442, 113539571, 113541911,
        113559755, 113616375, 113640435, 113663233, 113687428, 113694664,
        113702488, 113938754, 113964701, 113991046, 114043692, 114050112,
        114086267, 114130467, 114135930, 114156692, 114159895, 114214674,
        114505558, 114764511, 114823113, 114890473, 114905898, 114990593,
        115061617};

        uint32_t A_end[1366] = {
            19147431, 19204213, 19347002, 19339623, 19361275, 19411281,
            19445686, 19518397, 19599694, 19620460, 19659609, 19748859,
            19752968, 19771378, 19770950, 19872000, 19918067, 20080506,
            20197300, 20199621, 20250754, 20343472, 20408192, 20434228,
            20471203, 20472852, 20521213, 20556918, 20651780, 20705752,
            20704360, 20733843, 20828834, 20904064, 20904900, 20978060,
            20997163, 21122651, 21124756, 21135063, 21193175, 21257972,
            21458145, 21480764, 21487043, 21495275, 21524308, 21525395,
            21535657, 21678161, 21682064, 21685339, 21701660, 21732128,
            21746496, 21750536, 22055126, 22084993, 22201914, 22288721,
            22335600, 22421895, 22419917, 22422811, 22424895, 22527765,
            22557393, 22574770, 22653734, 22780009, 22815666, 22933278,
            23128203, 23185603, 23296101, 23327919, 23327223, 23389359,
            23489906, 23448310, 23560743, 23637293, 23766578, 23795752,
            23819792, 23827300, 24039104, 24114771, 24224478, 24385284,
            24355712, 24432158, 24448580, 24511711, 24533382, 24596223,
            24595942, 24633601, 24653774, 24814238, 24843745, 24852319,
            24922174, 24915319, 24930598, 25061346, 25006918, 25025385,
            25039750, 25158505, 25228334, 25241170, 25253051, 25393008,
            25453272, 25570956, 25580925, 25606280, 25609202, 25633961,
            25683998, 25942013, 26026113, 26040788, 26068982, 26084994,
            26189344, 26202899, 26292348, 26449911, 26773395, 26773256,
            26816083, 26840078, 27052181, 27067602, 27085350, 27148462,
            27421007, 27468345, 27791302, 27799171, 28081191, 28269060,
            28532703, 28618137, 28674467, 28664483, 28665497, 28670811,
            28697260, 29143887, 29303329, 29348287, 29484569, 29553087,
            29602276, 29636254, 29751681, 29874763, 29920085, 30135892,
            30317686, 30448001, 30436932, 30619790, 30621152, 30663758,
            30746220, 30784622, 30958042, 30964839, 30986932, 31043825,
            31197776, 31238380, 31418481, 31514667, 31543717, 31567922,
            31583646, 31590439, 31630698, 31652449, 31705692, 31945999,
            32026714, 32154395, 32227017, 32206685, 32261657, 32328562,
            32479944, 32494801, 32494960, 32488181, 32538777, 32573002,
            32674730, 32749222, 33138124, 33138393, 33140390, 33179127,
            33211783, 33439544, 33450209, 33444805, 33456520, 33576795,
            33934866, 34135900, 34138895, 34144821, 34193335, 34454014,
            34644363, 34824292, 34955923, 35024868, 35080394, 35138499,
            35179099, 35357657, 35460209, 35564882, 35569182, 35629587,
            35760592, 35922262, 35997598, 36041742, 36043604, 36032040,
            36270937, 36473054, 36625756, 36641142, 36802177, 36811320,
            36806861, 36818975, 36855456, 36881626, 36992614, 37246886,
            37414779, 37435801, 37476700, 37500736, 37866655, 37875553,
            38049311, 38055832, 38120730, 38085578, 38284565, 38286872,
            38404417, 38412764, 38416966, 38433149, 38469832, 38478360,
            38561548, 38684405, 38703540, 39060204, 39114322, 39180969,
            39154704, 39383428, 39394772, 39529994, 39751624, 39803397,
            39935552, 39949939, 39951214, 40068688, 40070642, 40125605,
            40194618, 40439476, 40461017, 40616052, 40739459, 40856730,
            40989064, 41081348, 41085049, 41104506, 41126585, 41235762,
            41279393, 41336684, 41342557, 41343553, 41647343, 41850520,
            41861223, 41873573, 41965318, 41992386, 42024798, 42004440,
            42105936, 42523109, 42524008, 42592414, 42647731, 43053893,
            43228591, 43302630, 43326939, 43516976, 43548586, 43621659,
            43633647, 43652895, 43730127, 43801149, 43820941, 43858487,
            44159624, 44168478, 44514763, 44569871, 44685688, 44819684,
            44835854, 44920505, 45043106, 45043642, 45099486, 45125268,
            45163284, 45213265, 45284656, 45373507, 45397218, 45427206,
            45428601, 45597483, 45759601, 45795321, 45868868, 45935191,
            45943788, 45954145, 45979586, 46005904, 46129916, 46187609,
            46253595, 46274088, 46297418, 46341834, 46451573, 46352505,
            46456676, 46607516, 46623011, 46685327, 46635417, 46748231,
            46793698, 46777299, 46792496, 46821404, 46980336, 46929851,
            46897729, 46916431, 46962313, 47041662, 47117463, 47184841,
            47282235, 47463283, 47500040, 47518025, 47755817, 47824386,
            47840640, 47944863, 47979881, 48011376, 48015798, 48134521,
            48131365, 48281633, 48384315, 48433070, 48446763, 48496089,
            48665060, 48785796, 48798173, 48882557, 49119314, 49129845,
            49142174, 49153175, 49203531, 49312042, 49536617, 49831663,
            49891284, 49953862, 49975083, 50044445, 50197793, 50392241,
            50421068, 50480139, 50477448, 50487736, 50514880, 50581748,
            50739569, 50816996, 50954925, 50978224, 51001881, 51063025,
            51075082, 51267763, 51446358, 51452950, 51515917, 51635745,
            52427049, 52951379, 53061757, 53038864, 53152228, 53178037,
            53210481, 53339298, 53363315, 53404196, 53491705, 53499400,
            53518226, 53523322, 53555562, 53633714, 53908494, 54115930,
            54240253, 54241743, 54261705, 54361006, 54425457, 54699860,
            54733012, 54822016, 54901731, 54970778, 55257495, 55253999,
            55265358, 55625067, 55635720, 55659574, 55660837, 55725826,
            55780553, 55872323, 56042341, 55964052, 55941469, 55959650,
            56035083, 56051938, 56151219, 56158770, 56248550, 56244487,
            56267877, 56366503, 56546071, 56749798, 56825196, 56852464,
            56890242, 56918798, 57044669, 57235683, 57247704, 57256313,
            57312717, 57372704, 57398112, 57406285, 57449180, 57553516,
            57617418, 57672785, 57714754, 57826089, 57802908, 57786936,
            57788378, 57890273, 57944016, 57911599, 57951793, 57979979,
            57978570, 58048239, 58157848, 58089414, 58125918, 58150652,
            58173064, 58231764, 58372204, 58373821, 58386447, 58433077,
            58454565, 58523057, 58516904, 58599899, 58765842, 58713346,
            58727299, 58770708, 58853729, 58910877, 58960672, 59022937,
            59177665, 59235014, 59355899, 59372028, 59374139, 59471629,
            59557948, 59706507, 59737393, 59800498, 59845121, 59856890,
            59894382, 59921655, 60195650, 60493573, 60474174, 60881693,
            61211553, 61248645, 61237883, 61470781, 61507199, 61527710,
            61551873, 61553884, 61717024, 61718622, 61956094, 61974100,
            62013500, 62015638, 62155814, 62197361, 62205915, 62414474,
            62531608, 62657669, 62721066, 62782467, 62889932, 62902078,
            62913864, 63102883, 63224295, 63233673, 63293178, 63313044,
            63387121, 63464700, 63425179, 63598538, 63591541, 63664409,
            63673901, 63697762, 63733420, 63736420, 63814240, 63868388,
            63891711, 63949810, 63987879, 64011258, 64528935, 64237601,
            64375358, 64313266, 64315765, 64413759, 64415233, 64416522,
            64697274, 64726845, 64711947, 64769967, 64809769, 64821658,
            64894220, 64924986, 64960349, 65067408, 65195765, 65271829,
            65283174, 65344433, 65353710, 65383724, 65444103, 65530611,
            65562858, 65616438, 65631242, 65638216, 65675236, 65753153,
            65760241, 65862166, 66062510, 66207455, 66325576, 66456256,
            66584492, 66720796, 66923397, 66901672, 66899951, 66965862,
            66986365, 67139235, 67134595, 67170019, 67177425, 67222592,
            67318828, 67253565, 67326541, 67341407, 67396340, 67470239,
            67527771, 67787397, 67846242, 67861974, 67953705, 67950594,
            67982305, 68099045, 68094373, 68101466, 68136555, 68152516,
            68172613, 68284746, 68374480, 68377269, 68397946, 68412187,
            68533309, 68554375, 68691502, 68836370, 68845819, 68877506,
            68973919, 69002942, 68978201, 69066139, 69125270, 69202971,
            69268758, 69393681, 69391667, 69402948, 69449959, 69465016,
            69452333, 69488357, 69699710, 69741991, 69810380, 69792803,
            69837730, 69941920, 69946272, 69992192, 70071044, 70079696,
            70128221, 70201994, 70406793, 70569774, 70651855, 70654767,
            70775452, 70750716, 70796505, 70816532, 70823665, 70826589,
            70857091, 70839574, 71123472, 70934619, 71010044, 71010061,
            71088624, 71092705, 71124554, 71172632, 71189517, 71213416,
            71200001, 71255603, 71274905, 71296302, 71385148, 71431669,
            71447111, 71677219, 71792190, 71846729, 71865510, 71876362,
            71922001, 71921984, 71924424, 72003419, 72035702, 72353439,
            72353445, 72386687, 72454025, 72480580, 72482473, 72525153,
            72542066, 72709143, 72774482, 72846923, 72886529, 72898446,
            73024610, 73300302, 73418920, 73724475, 73857221, 73950436,
            74006858, 74090683, 74360269, 74479952, 74721469, 75008211,
            74958257, 74987977, 74988127, 75166259, 75439445, 75439001,
            75444035, 75501574, 75603826, 75625277, 75634206, 75678893,
            75810995, 75901135, 75877145, 76087573, 76195761, 76368108,
            76386955, 76469011, 76551443, 76528975, 76566177, 76592004,
            76612760, 76658442, 76659581, 76680512, 76775561, 76827288,
            76868542, 76892840, 76964276, 76999419, 76993644, 77002287,
            77065719, 77117169, 77973099, 78120570, 78172191, 78235339,
            78388016, 78391500, 78459660, 78509602, 78513329, 78560100,
            78767838, 78808705, 78855970, 78937549, 78970957, 78998715,
            79208540, 79246160, 79463327, 79496405, 79551302, 79633881,
            79757461, 79884162, 79998631, 80460159, 80463570, 80595902,
            80636690, 80655639, 80681863, 80750580, 80770098, 80789382,
            80884906, 80963067, 81002583, 81030521, 81032253, 81052868,
            81045957, 81281194, 81322649, 81436344, 81486206, 81746251,
            81765389, 81814929, 81962887, 81858072, 81873856, 81972927,
            82085867, 82169564, 82211859, 82265408, 82293644, 82355696,
            82395869, 82462420, 82510741, 82590676, 82587653, 82590057,
            82748550, 82694976, 82756683, 82761532, 82909246, 83049762,
            83068547, 83132207, 83172206, 83214144, 83394710, 83418010,
            83439881, 83564027, 83641765, 83649574, 83714159, 83739670,
            83792579, 83801437, 83817064, 83829778, 83875904, 83893718,
            83912581, 83928218, 83954563, 84037400, 84036376, 84051492,
            84162902, 84158661, 84130121, 84262016, 84186833, 84264073,
            84261338, 84286556, 84296390, 84328536, 84663208, 84487785,
            84531371, 84536864, 84555458, 84570142, 84673008, 84637784,
            84702446, 84699144, 84689637, 84905651, 85045549, 85057064,
            85049644, 85058927, 85064837, 85057018, 85106498, 85193761,
            85185429, 85198120, 85223355, 85270084, 85282949, 85370903,
            85403264, 85419741, 85477896, 85555600, 85670196, 85713914,
            85803606, 85848855, 85883917, 85929693, 85937632, 86059489,
            86269513, 86388291, 86515057, 86548599, 86580272, 86637189,
            86652380, 86655549, 86659417, 86676244, 86675121, 86752505,
            86976665, 86994802, 87172740, 87215061, 87320900, 87353759,
            87345985, 87370513, 87436002, 87390588, 87434385, 87437688,
            87453551, 87552671, 87605422, 87642390, 87685446, 87780653,
            87831617, 87838132, 87910938, 87963792, 88098472, 88033399,
            87987653, 88003710, 88022806, 88277221, 88399089, 88472105,
            88519972, 88638758, 88683247, 88958850, 88964906, 89021661,
            89156315, 89186533, 89240133, 89328332, 89414766, 89422930,
            89739118, 89850213, 89926260, 90018535, 90013523, 90243749,
            90270699, 90213131, 90249594, 90337083, 90504437, 90355395,
            90460217, 90549572, 90643023, 90684527, 90688461, 90726998,
            90773616, 90784356, 90785702, 90864805, 90914148, 91027630,
            91046227, 91108627, 91220379, 91320133, 91316366, 91330237,
            91510784, 91519741, 91593166, 91703008, 91614248, 91724663,
            91725620, 91734356, 91754179, 91904843, 92111830, 92197506,
            92392377, 92420049, 92499229, 92520236, 92534600, 92553031,
            92641371, 92734932, 92715382, 92825376, 92860157, 92955761,
            93000176, 93008852, 93018904, 93057421, 93082163, 93139233,
            93245963, 93471647, 93557015, 93707759, 93810687, 93825647,
            93883679, 93893987, 94047653, 94103934, 94045542, 94149294,
            94092721, 94091194, 94149294, 94183775, 94230577, 94267989,
            94269404, 94392525, 94521781, 94702458, 94770323, 94832055,
            94868675, 95004679, 95070601, 95201279, 95294778, 95438795,
            95575775, 95644952, 95716554, 95790088, 95783044, 95838521,
            95892358, 95894047, 95971179, 95980662, 96023675, 95994432,
            96015479, 96024581, 96084863, 96215575, 96420106, 96793262,
            96801552, 96860558, 96881494, 96925919, 97014736, 97043086,
            97051641, 97086775, 97265609, 97271264, 97331903, 97340081,
            97347899, 97347922, 97491016, 97537120, 97618323, 97632717,
            97642010, 97900703, 97944188, 98073067, 98108407, 98155405,
            98177292, 98262984, 98443552, 98410408, 98463604, 98533077,
            98602049, 98898599, 98904285, 99067544, 99094844, 99130481,
            99151248, 99240420, 99258251, 99256778, 99280491, 99301396,
            99316179, 99548740, 99631758, 99688674, 99794272, 99943204,
            100007265, 100082459, 100113023, 100226405, 100495966, 100574184,
            100605224, 100704379, 100740238, 100738820, 100759019, 100775641,
            100783344, 100826073, 100855528, 100901146, 101085419, 101208768,
            101211529, 101397458, 101417256, 101443442, 101486973, 101510586,
            101608662, 101625871, 101768616, 101896428, 101915895, 101992403,
            102044138, 102133956, 102178791, 102283051, 102459014, 102443730,
            102565982, 102635496, 102649653, 102669629, 102863870, 102928730,
            102944172, 102952751, 103097217, 103129601, 103359123, 103501963,
            103517049, 103603544, 103688922, 103712463, 103722284, 103722538,
            103875570, 104015573, 104040753, 104055382, 104042724, 104092258,
            104174333, 104279157, 104384569, 104397435, 104534248, 104607695,
            104683644, 104785363, 104810017, 104883804, 104990005, 104999305,
            105106630, 105118455, 105132809, 105195995, 105238751, 105565248,
            105425266, 105523082, 105555638, 105762061, 105794060, 105815568,
            105956497, 106113308, 106194716, 106309632, 106307937, 106356335,
            106369179, 106449921, 106476386, 106667718, 106718214, 106721103,
            106989086, 107014129, 107038693, 107088714, 107163124, 107206256,
            107427737, 107444348, 107465040, 107518160, 107525463, 107701450,
            107713652, 107837480, 108120035, 108191834, 108303236, 108448581,
            108451712, 108561472, 108603844, 108774809, 108849196, 108886818,
            108908350, 108910598, 108951803, 109005966, 109057822, 109093114,
            109267251, 109362442, 109421222, 109477147, 109512882, 109658619,
            109677459, 109689974, 109686517, 109695151, 109787941, 109796068,
            109884074, 110262756, 110222462, 110229127, 110281933, 110329887,
            110383698, 110385387, 110381653, 110465092, 110469229, 110481807,
            110761471, 110763477, 110798307, 110969740, 111126509, 111136077,
            111428336, 111436305, 111589000, 111601098, 111630984, 111703667,
            111716101, 111748432, 111776187, 112046438, 112155606, 112321962,
            112340180, 112598333, 112624555, 112649034, 112698684, 112837911,
            113267813, 113282325, 113326195, 113441090, 113456736, 113481644,
            113492529, 113539647, 113541977, 113562282, 113616435, 113640557,
            113663296, 113692863, 113695196, 113703258, 113939729, 113979268,
            113991116, 114048904, 114054638, 114094988, 114130543, 114137087,
            114171876, 114176793, 114214806, 114506637, 114764591, 114826687,
            114891658, 114906027, 114991881, 115063113};
            
        int32_t A_cipos_0[1366] = {
            -45, -9, -1000, -1000, -1000, -50, -500, -204, 0, -500, -3, -1,
            -142, -22, 0, -7, -61, -267, -2, -28, -337, -393, -50, -47, -182,
            -1000, 0, -32, -1, -32, 0, -199, -12, -30, -19, -339, -150, 0, 0,
            -50, -50, -30, -82, -33, -19, -4, -110, -48, -4, -50, 0, -48, -58,
            -150, -107, -181, -1000, 0, -1000, -270, 0, -222, -170, 0, -56,
            -301, 0, 0, -328, -19, 0, -342, -244, -12, -500, -148, -443, -27,
            0, -500, -408, -361, -39, -50, -108, -202, -26, -19, 0, -8, -78,
            -50, 0, -9, 0, -201, -25, -325, -259, -335, -2, -12, -1000, 0, 0,
            -50, 0, -6, -30, -2, -50, -1000, -19, -1000, -250, -35, -60, -150,
            -45, 0, -27, -1, -35, 0, -3, -13, -10, 0, -133, -12, -4, -3, -59,
            -3, 0, -500, -170, -31, -42, -279, -42, 0, -6, -1, -28, -297,
            -1000, -182, -88, -244, -13, -26, -1000, -44, -500, -500, 0, -13,
            -1, 0, -1, -35, -150, -150, -240, -26, 0, -62, -110, -10, -262,
            -40, 0, 0, -2, 0, 0, -23, -7, -13, -34, -17, 0, 0, -193, -208, 0,
            -383, 0, -500, 0, -89, -14, -47, 0, -28, -6, -58, -1, 0, -1000, 0,
            -76, -27, -26, -500, 0, 0, -500, -268, -21, -500, -8, 0, -314, -14,
            -267, 0, 0, 0, -92, 0, -14, -240, 0, 0, -13, 0, -15, -17, -500,
            -1000, 0, 0, -87, -26, -6, -500, -174, -47, 0, 0, -17, -38, 0,
            -292, 0, -24, -40, -265, 0, -253, -387, -33, -289, 0, -57, -23,
            -19, -9, -5, -17, 0, -24, -26, -154, -124, 0, -34, -9, 0, -126,
            -165, -11, 0, 0, -330, 0, -1, -150, 0, -3, -159, -78, -263, -50, 0,
            0, -40, -14, -283, -19, 0, 0, -225, 0, 0, 0, -25, -50, -11, -47,
            -2, -35, -150, -2, -286, -27, -96, -10, -150, -1, -1000, -6, -50,
            -79, -1000, 0, -164, -38, 0, 0, 0, -10, -13, -500, 0, -3, -31, -20,
            -5, 0, 0, 0, -50, -78, 0, -500, -4, 0, -2, -17, 0, 0, -267, -99,
            -22, -1000, -40, -149, -1000, -119, -17, 0, -5, -82, -193, -29,
            -17, 0, 0, -32, -31, -500, -50, -500, -500, -74, -87, -500, 0,
            -500, 0, -14, -12, 0, -196, 0, -8, -500, -10, 0, -295, -2, -125,
            -14, 0, -150, 0, -500, 0, -500, -68, 0, 0, -16, -87, -342, -98, -6,
            -25, 0, 0, 0, -198, -63, -432, -39, -1, 0, 0, -32, -1000, -36, -54,
            -62, -150, -58, -1, -150, -42, 0, 0, -415, 0, -17, 0, -58, -9, 0,
            -1, -84, 0, -50, -500, -2, -500, -150, -500, -15, -291, 0, 0, 0,
            -128, -218, 0, -1000, 0, -7, 0, -220, -33, -27, -1, -170, -192,
            -49, -4, -2, -500, 0, -7, -249, -50, -59, 0, -122, -234, 0, -31, 0,
            0, -85, -41, 0, 0, 0, -1000, 0, 0, 0, -3, -500, -234, -362, -150,
            -29, -1, -500, -500, -46, -49, 0, -19, -7, -50, -70, 0, -13, -349,
            -269, -8, -1000, 0, -1000, -500, -12, -261, -114, -397, -88, 0,
            -500, 0, -26, -15, -108, -1000, -172, -1, -179, 0, 0, -500, -1000,
            -150, 0, -28, 0, -500, 0, -255, -41, -11, -232, -291, 0, -34, -76,
            0, 0, -11, -4, 0, 0, -13, -31, 0, -10, 0, 0, -12, 0, -173, -26,
            -198, -14, 0, -351, -17, -217, 0, -14, -50, 0, -1000, -18, -30,
            -50, -3, -1000, -7, -50, -16, 0, 0, -186, -54, -6, -3, -50, 0,
            -500, -1000, 0, -2, -20, -315, -50, -17, -8, -19, -20, -25, -324,
            -31, -23, -1000, -1, -1000, -312, -158, -277, -12, -317, -500,
            -500, -50, -500, 0, -8, -60, -8, 0, -14, -500, -21, -500, 0, -341,
            0, 0, -50, -353, -44, -45, -537, -23, -3, 0, 0, -12, -392, -10, -3,
            0, -1000, -21, 0, -408, -210, -33, -3, -3, -164, 0, -228, -35, -13,
            -1000, -1, -500, -11, -51, -500, -120, 0, -10, 0, -22, -11, -43,
            -1, -232, -1000, -3, -25, -1000, -40, -5, 0, -21, -389, -500, -58,
            -229, -162, -301, -66, 0, -18, -16, -251, -302, -17, -6, 0, -252,
            -500, 0, -247, 0, -106, 0, -1000, 0, -98, -262, 0, -500, 0, -219,
            -13, 0, -500, -1000, -23, -2, -1, -500, -1000, -265, 0, -4, -232,
            -32, -10, 0, 0, -1000, 0, -359, -200, -36, 0, -31, -500, 0, -27,
            -50, -7, -13, 0, -18, -40, 0, 0, -500, -7, 0, -337, -134, 0, 0, 0,
            0, 0, -177, -165, -50, -93, 0, -432, 0, -276, -1000, -304, 0, -382,
            -16, -111, -23, -218, -254, -242, -22, -12, 0, 0, 0, -20, -199, -1,
            -134, 0, 0, 0, -160, 0, -24, -500, -4, -94, -500, 0, -4, 0, 0,
            -500, -8, -2, 0, -4, -5, 0, -500, -11, 0, 0, 0, -23, 0, 0, 0, -5,
            -229, -46, 0, -318, -23, 0, -275, -264, -104, -231, 0, 0, -142, 0,
            -257, 0, -368, 0, 0, -31, -26, 0, 0, -22, 0, -418, 0, 0, 0, -2, 0,
            0, 0, -500, -83, 0, -233, -218, 0, -150, -8, -290, 0, -150, 0, -25,
            0, -244, -500, 0, -5, -199, -31, 0, -51, 0, 0, 0, -500, 0, 0, 0,
            -500, -56, -500, 0, 0, -169, -500, 0, -202, -22, 0, -278, -30, 0,
            0, -2, -101, 0, -26, 0, -299, 0, -21, 0, -24, -11, -65, -235, -500,
            -25, 0, -38, 0, -331, -3, -20, -17, 0, -6, -47, -210, -13, 0,
            -1000, -500, -500, 0, -500, -51, -500, 0, -1000, -295, -50, 0, -43,
            -1, -28, -58, 0, -3, -1000, 0, 0, 0, 0, -315, -216, -1000, -153,
            -12, -500, -5, -19, -1000, -50, -3, 0, 0, -30, -6, -22, -65, -1000,
            0, -235, 0, -17, 0, 0, -150, -9, -310, -135, -21, -6, -31, 0, -60,
            -132, 0, -500, -500, -18, 0, -16, -27, 0, -186, 0, -30, -216,
            -1000, -18, 0, -500, 0, -20, 0, 0, -30, -11, -69, -1000, -30, -50,
            -62, -6, -157, -219, -500, -1000, -50, -395, 0, -5, 0, -15, -1,
            -500, 0, -166, -20, 0, 0, -19, -500, -150, 0, -35, -84, 0, -5, -37,
            0, -443, -5, -7, 0, -26, -500, -500, 0, -96, -16, 0, -500, 0,
            -1000, -6, -18, -34, 0, -2, 0, -11, -354, -21, -149, 0, -50, -102,
            -334, -9, -24, -500, -209, -2, 0, -198, -309, -348, -8, -65, -602,
            -268, -26, -142, -10, 0, 0, 0, 0, -500, -2, 0, -141, 0, 0, -29,
            -49, -25, -10, 0, -2, -5, -2, -55, 0, -259, 0, 0, -17, 0, -2, -31,
            0, -1000, -183, -5, -43, -232, -21, -50, 0, -150, -104, -1000,
            -1000, 0, -275, -13, -28, -101, -10, -11, -192, -7, -3, -49, -17,
            0, -8, -52, -17, 0, -279, -9, -9, -99, -248, -15, -14, -500, -199,
            -38, -500, -77, -16, -276, 0, -362, -14, -92, -500, -34, -20, -203,
            0, -246, -1000, 0, -15, -5, -233, -13, -18, -157, -31, -16, 0, -52,
            0, -29, -43, -207, 0, 0, -500, -98, -89, 0, -8, -22, -308, -118,
            -51, -371, -19, -9, -102, -500, -500, 0, 0, 0, 0, -449, 0, -15,
            -500, -17, -8, -177, -22, -8, -197, 0, 0, -323, -3, 0, -16, 0, 0,
            -10, -3, -179, 0, -3, -220, 0, -500, 0, 0, 0, -1000, -1000, -19,
            -250, -9, -318, -500, 0, -161, -332, -16, -4, 0, -5, -253, -7,
            -1000, -9, -39, 0, -4, -30, 0, 0, 0, -20, 0, -175, -18, -24, -500,
            0, -15, -211, -19, 0, -181, -73, -500, 0, -33, -294, -27, -28, 0,
            0, 0, -22, -11, -244, -13, -1, -341, -1, -7, -3, -500, -305, 0,
            -476, 0, -216, 0, 0, -150, 0, -323, -150, -47, 0, 0, 0, -28, 0, -1,
            -43, 0, -268, -177, 0, -169, -205, -96, -11, -1, 0, -46, -500, 0,
            -9, 0, -1000, -12, -139, 0, -150, 0, 0, -100, 0, -1000, -11, -3, 0,
            -242, 0, 0, 0, -11, -225, -2, 0, 0, -55, 0, 0, -121, -1000, -150,
            -50, -30, -500, -231, -10, 0, 0, 0, -8, 0, 0, 0, -35, -66, -50,
            -332, -500, 0, -500, -190, -150, 0, -269, -1000, -500, 0, 0, 0,
            -166, -7, 0, 0, -35};

        int32_t A_cipos_1[1366] = {
            45, 9, 500, 500, 500, 50, 0, 73, 0, 0, 3, 1, 0, 23, 1, 7, 61, 0, 3,
            29, 0, 0, 50, 48, 76, 500, 0, 32, 1, 33, 0, 0, 12, 31, 20, 0, 150,
            0, 0, 50, 50, 31, 50, 34, 19, 4, 50, 48, 4, 50, 1, 48, 58, 150, 0,
            0, 500, 0, 500, 0, 0, 277, 0, 0, 57, 0, 0, 0, 0, 20, 0, 0, 0, 13,
            0, 0, 0, 27, 0, 0, 0, 0, 0, 50, 0, 0, 26, 20, 0, 8, 50, 51, 0, 10,
            0, 0, 25, 0, 0, 0, 3, 12, 500, 0, 0, 50, 0, 7, 31, 2, 50, 500, 19,
            500, 0, 0, 60, 0, 45, 0, 27, 1, 36, 0, 3, 13, 10, 0, 0, 13, 5, 3,
            59, 4, 1, 0, 0, 32, 43, 0, 42, 0, 7, 2, 29, 0, 500, 134, 50, 0, 13,
            26, 500, 45, 0, 0, 0, 13, 1, 0, 2, 36, 150, 150, 0, 26, 0, 50, 50,
            11, 0, 41, 1, 0, 2, 0, 0, 24, 7, 14, 34, 18, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 50, 15, 48, 0, 28, 6, 59, 2, 0, 500, 0, 0, 27, 27, 0, 0, 0,
            0, 0, 22, 0, 8, 0, 0, 14, 0, 0, 0, 0, 0, 0, 14, 0, 0, 1, 14, 0, 15,
            18, 0, 500, 0, 0, 0, 26, 7, 0, 0, 48, 0, 0, 17, 39, 0, 0, 0, 25,
            41, 0, 0, 0, 0, 33, 0, 0, 57, 23, 19, 10, 6, 0, 0, 25, 27, 107, 0,
            0, 34, 9, 0, 0, 0, 12, 0, 0, 50, 0, 1, 150, 0, 4, 0, 0, 0, 50, 0,
            0, 41, 15, 0, 19, 0, 0, 0, 0, 0, 0, 25, 50, 11, 47, 3, 36, 150, 3,
            0, 27, 0, 11, 150, 1, 500, 6, 50, 50, 500, 1, 0, 39, 0, 0, 0, 11,
            13, 0, 0, 3, 31, 21, 6, 0, 0, 0, 50, 50, 0, 0, 4, 0, 2, 17, 0, 0,
            0, 0, 22, 500, 40, 50, 500, 0, 17, 0, 6, 50, 0, 29, 18, 0, 0, 32,
            32, 0, 50, 0, 0, 0, 50, 0, 0, 0, 0, 14, 13, 0, 0, 0, 9, 0, 10, 0,
            0, 2, 0, 14, 1, 150, 0, 0, 0, 0, 50, 0, 0, 16, 0, 0, 50, 6, 25, 0,
            0, 0, 0, 0, 0, 39, 1, 0, 0, 32, 500, 36, 54, 62, 0, 58, 2, 150, 42,
            0, 0, 0, 0, 17, 0, 0, 9, 0, 2, 50, 0, 50, 0, 2, 0, 150, 0, 15, 0,
            0, 0, 0, 0, 0, 0, 500, 0, 7, 0, 0, 33, 28, 2, 0, 0, 49, 4, 2, 0, 0,
            7, 0, 50, 0, 0, 0, 0, 0, 31, 0, 1, 0, 41, 0, 0, 0, 500, 0, 0, 0, 4,
            0, 0, 0, 0, 29, 1, 0, 0, 46, 50, 0, 19, 8, 50, 50, 0, 13, 0, 0, 9,
            500, 1, 500, 0, 12, 0, 97, 0, 50, 0, 0, 0, 26, 15, 0, 500, 0, 1, 0,
            1, 0, 0, 500, 150, 0, 28, 0, 0, 0, 0, 41, 12, 0, 0, 0, 34, 50, 0,
            0, 12, 4, 0, 0, 13, 31, 0, 10, 0, 0, 13, 0, 0, 26, 0, 14, 0, 0, 17,
            50, 0, 14, 50, 0, 500, 18, 31, 50, 3, 500, 8, 50, 16, 0, 0, 0, 55,
            7, 4, 50, 0, 0, 500, 0, 3, 21, 0, 50, 17, 8, 19, 21, 26, 0, 31, 24,
            500, 2, 500, 0, 0, 0, 13, 0, 0, 0, 50, 0, 0, 9, 0, 8, 0, 15, 0, 21,
            0, 0, 0, 0, 0, 50, 0, 45, 46, 0, 24, 4, 0, 0, 13, 0, 10, 4, 1, 500,
            22, 0, 0, 0, 34, 3, 3, 0, 0, 0, 35, 14, 500, 1, 0, 11, 51, 0, 0, 0,
            11, 0, 23, 11, 44, 1, 0, 500, 3, 26, 500, 0, 5, 0, 21, 0, 0, 50, 0,
            0, 0, 50, 0, 18, 17, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 0, 500, 0,
            50, 0, 0, 0, 0, 0, 13, 0, 0, 500, 23, 2, 1, 0, 500, 0, 0, 5, 0, 33,
            11, 0, 0, 500, 0, 0, 0, 37, 0, 0, 0, 0, 28, 50, 8, 14, 0, 18, 41,
            0, 0, 0, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 50, 68, 0, 0, 0, 0, 500,
            0, 0, 0, 16, 50, 23, 0, 0, 0, 23, 12, 0, 0, 0, 21, 0, 1, 0, 0, 0,
            0, 0, 0, 25, 0, 4, 50, 0, 1, 5, 0, 0, 0, 9, 2, 1, 5, 5, 0, 0, 12,
            0, 0, 0, 23, 0, 0, 0, 5, 0, 46, 1, 0, 23, 0, 0, 0, 50, 0, 1, 1, 0,
            0, 0, 0, 0, 0, 0, 32, 26, 0, 0, 23, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 150, 9, 0, 0, 150, 0, 26, 0, 141, 0, 0, 6, 0, 32, 0,
            52, 0, 0, 0, 0, 0, 0, 0, 0, 57, 0, 1, 0, 0, 0, 0, 0, 22, 0, 0, 30,
            0, 0, 3, 0, 0, 26, 0, 0, 0, 22, 0, 25, 11, 66, 0, 0, 25, 0, 39, 0,
            0, 3, 21, 18, 0, 6, 0, 0, 14, 0, 500, 0, 0, 0, 0, 52, 0, 0, 500, 0,
            50, 0, 44, 1, 29, 59, 0, 3, 500, 0, 0, 0, 0, 0, 0, 500, 0, 12, 0,
            6, 20, 500, 50, 3, 0, 0, 30, 7, 23, 66, 500, 1, 0, 0, 18, 0, 1,
            150, 9, 0, 0, 21, 7, 31, 0, 61, 0, 0, 0, 0, 19, 1, 17, 27, 0, 0, 0,
            31, 0, 500, 18, 0, 0, 0, 21, 0, 0, 30, 11, 50, 500, 30, 50, 0, 7,
            0, 0, 0, 500, 50, 0, 0, 5, 0, 16, 2, 0, 0, 0, 20, 0, 0, 19, 0, 150,
            0, 36, 50, 0, 6, 37, 0, 0, 5, 7, 0, 26, 0, 0, 0, 50, 17, 0, 0, 0,
            500, 6, 18, 34, 0, 2, 0, 11, 0, 21, 0, 0, 50, 0, 0, 10, 24, 0, 0,
            3, 0, 0, 0, 0, 8, 50, 0, 0, 26, 0, 10, 0, 0, 0, 0, 0, 3, 0, 0, 0,
            0, 30, 50, 26, 10, 0, 3, 5, 2, 55, 0, 0, 0, 1, 17, 0, 3, 31, 0,
            500, 0, 6, 43, 0, 21, 55, 0, 150, 0, 500, 500, 0, 0, 13, 28, 0, 10,
            11, 0, 8, 3, 49, 17, 0, 8, 0, 18, 0, 0, 9, 9, 0, 0, 16, 15, 0, 0,
            39, 0, 50, 16, 0, 0, 0, 15, 0, 0, 35, 20, 0, 0, 0, 500, 0, 16, 6,
            0, 14, 19, 0, 31, 17, 0, 52, 0, 29, 43, 0, 0, 0, 0, 50, 50, 0, 9,
            23, 0, 50, 51, 0, 20, 9, 0, 0, 0, 0, 0, 0, 0, 0, 1, 15, 0, 17, 8,
            0, 22, 9, 0, 0, 0, 0, 3, 0, 17, 0, 0, 10, 4, 0, 0, 3, 0, 0, 0, 0,
            0, 0, 500, 500, 20, 50, 9, 0, 0, 0, 0, 0, 16, 5, 1, 5, 0, 8, 500,
            10, 40, 0, 4, 31, 0, 0, 0, 21, 0, 0, 18, 24, 0, 0, 16, 0, 20, 0, 0,
            0, 0, 0, 34, 0, 27, 29, 0, 0, 0, 22, 12, 0, 0, 2, 0, 1, 7, 4, 0, 0,
            0, 0, 0, 0, 0, 0, 150, 0, 0, 150, 47, 0, 0, 0, 28, 0, 2, 44, 0, 0,
            0, 0, 0, 0, 50, 12, 2, 0, 47, 0, 0, 10, 0, 500, 13, 0, 0, 150, 0,
            0, 50, 0, 500, 12, 3, 0, 0, 0, 0, 0, 12, 0, 3, 0, 0, 55, 0, 0, 0,
            500, 150, 50, 31, 0, 0, 10, 0, 0, 0, 9, 0, 0, 0, 35, 50, 50, 0, 0,
            0, 0, 0, 150, 0, 0, 500, 0, 0, 0, 0, 0, 7, 0, 0, 36};

        int32_t A_ciend_0[1366] = {
            -45, -9, -500, -500, -500, -50, 0, -50, 0, 0, -3, -1, 0, -22, 0,
            -7, -61, 0, -2, -28, 0, 0, -50, -47, -50, -500, 0, -32, -1, -32, 0,
            0, -12, -30, -19, 0, -150, 0, 0, -50, -50, -30, -50, -33, -19, -4,
            -50, -48, -4, -50, 0, -48, -58, -150, 0, 0, -500, 0, -500, 0, 0,
            -80, 0, 0, -56, 0, 0, 0, 0, -19, 0, 0, 0, -12, 0, 0, 0, -27, 0, 0,
            0, 0, 0, -50, 0, 0, -26, -19, 0, -8, -50, -50, 0, -9, 0, 0, -25, 0,
            0, 0, -2, -12, -500, 0, 0, -50, 0, -6, -30, -2, -50, -500, -19,
            -500, 0, 0, -60, 0, -45, 0, -27, -1, -35, 0, -3, -13, -10, 0, 0,
            -12, -4, -3, -59, -3, 0, 0, 0, -31, -42, 0, -42, 0, -6, -1, -28, 0,
            -500, -67, -50, 0, -13, -26, -500, -44, 0, 0, 0, -13, -1, 0, -1,
            -35, -150, -150, 0, -26, 0, -50, -50, -10, 0, -40, 0, 0, -2, 0, 0,
            -23, -7, -13, -34, -17, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50, -14, -47,
            0, -28, -6, -58, -1, 0, -500, 0, 0, -27, -26, 0, 0, 0, 0, 0, -21,
            0, -8, 0, 0, -14, 0, 0, 0, 0, 0, 0, -14, 0, 0, 0, -13, 0, -15, -17,
            0, -500, 0, 0, 0, -26, -6, 0, 0, -47, 0, 0, -17, -38, 0, 0, 0, -24,
            -40, 0, 0, 0, 0, -33, 0, 0, -57, -23, -19, -9, -5, 0, 0, -24, -26,
            -264, 0, 0, -34, -9, 0, 0, 0, -11, 0, 0, -123, 0, -1, -150, 0, -3,
            0, 0, 0, -50, 0, 0, -40, -14, 0, -19, 0, 0, 0, 0, 0, 0, -25, -50,
            -11, -47, -2, -35, -150, -2, 0, -27, 0, -10, -150, -1, -500, -6,
            -50, -50, -500, 0, 0, -38, 0, 0, 0, -10, -13, 0, 0, -3, -31, -20,
            -5, 0, 0, 0, -50, -50, 0, 0, -4, 0, -2, -17, 0, 0, 0, 0, -22, -500,
            -40, -50, -500, 0, -17, 0, -5, -50, 0, -29, -17, 0, 0, -32, -31, 0,
            -67, 0, 0, 0, -50, 0, 0, 0, 0, -14, -12, 0, 0, 0, -8, 0, -10, 0, 0,
            -2, 0, -14, 0, -150, 0, 0, 0, 0, -50, 0, 0, -16, 0, 0, -50, -6,
            -25, 0, 0, 0, 0, 0, 0, -39, -1, 0, 0, -32, -500, -36, -54, -62, 0,
            -58, -1, -150, -42, 0, 0, 0, 0, -17, 0, 0, -9, 0, -1, -50, 0, -50,
            0, -2, 0, -150, 0, -15, 0, 0, 0, 0, 0, 0, 0, -500, 0, -7, 0, 0,
            -33, -27, -1, 0, 0, -49, -4, -2, 0, 0, -7, 0, -50, 0, 0, 0, 0, 0,
            -31, 0, 0, 0, -41, 0, 0, 0, -500, 0, 0, 0, -3, 0, 0, 0, 0, -29, -1,
            0, 0, -46, -49, 0, -19, -7, -50, -50, 0, -13, 0, 0, -8, -500, 0,
            -500, 0, -12, 0, -94, 0, -50, 0, 0, 0, -26, -15, 0, -500, 0, -1, 0,
            0, 0, 0, -500, -150, 0, -28, 0, 0, 0, 0, -41, -11, 0, 0, 0, -34,
            -50, 0, 0, -11, -4, 0, 0, -13, -31, 0, -10, 0, 0, -12, 0, 0, -26,
            0, -14, 0, 0, -17, -71, 0, -14, -50, 0, -500, -18, -30, -50, -3,
            -500, -7, -50, -16, 0, 0, 0, -54, -6, -3, -50, 0, 0, -500, 0, -2,
            -20, 0, -50, -17, -8, -19, -20, -25, 0, -31, -23, -500, -1, -500,
            0, 0, 0, -12, 0, 0, 0, -50, 0, 0, -8, 0, -8, 0, -14, 0, -21, 0, 0,
            0, 0, 0, -96, 0, -44, -45, 0, -23, -3, 0, 0, -12, 0, -10, -3, 0,
            -500, -21, 0, 0, 0, -33, -3, -3, 0, 0, 0, -35, -13, -500, -1, 0,
            -11, -51, 0, 0, 0, -10, 0, -22, -11, -43, -1, 0, -500, -3, -25,
            -500, 0, -5, 0, -21, 0, 0, -50, 0, 0, 0, -50, 0, -18, -16, 0, 0, 0,
            -6, 0, 0, 0, 0, 0, 0, 0, 0, -500, 0, -50, 0, 0, 0, 0, 0, -13, 0, 0,
            -500, -23, -2, -1, 0, -500, 0, 0, -4, 0, -32, -10, 0, 0, -500, 0,
            0, 0, -36, 0, 0, 0, 0, -27, -50, -7, -13, 0, -18, -40, 0, 0, 0, -7,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -50, -50, 0, 0, 0, 0, -500, 0, 0, 0,
            -16, -50, -23, 0, 0, 0, -22, -12, 0, 0, 0, -20, 0, -1, 0, 0, 0, 0,
            0, 0, -24, 0, -4, -50, 0, 0, -4, 0, 0, 0, -8, -2, 0, -4, -5, 0, 0,
            -11, 0, 0, 0, -23, 0, 0, 0, -5, 0, -46, 0, 0, -23, 0, 0, 0, -50, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, -31, -26, 0, 0, -22, 0, 0, 0, 0, 0, -2,
            0, 0, 0, 0, 0, 0, 0, 0, 0, -150, -8, 0, 0, -150, 0, -25, 0, -90, 0,
            0, -5, 0, -31, 0, -51, 0, 0, 0, 0, 0, 0, 0, 0, -56, 0, 0, 0, 0, 0,
            0, 0, -22, 0, 0, -30, 0, 0, -2, 0, 0, -26, 0, 0, 0, -21, 0, -24,
            -11, -65, 0, 0, -25, 0, -38, 0, 0, -3, -20, -17, 0, -6, 0, 0, -13,
            0, -500, 0, 0, 0, 0, -51, 0, 0, -500, 0, -50, 0, -43, -1, -28, -58,
            0, -3, -500, 0, 0, 0, 0, 0, 0, -500, 0, -12, 0, -5, -19, -500, -50,
            -3, 0, 0, -30, -6, -22, -65, -500, 0, 0, 0, -17, 0, 0, -150, -9, 0,
            0, -21, -6, -31, 0, -60, 0, 0, 0, 0, -18, 0, -16, -27, 0, 0, 0,
            -30, 0, -500, -18, 0, 0, 0, -20, 0, 0, -30, -11, -50, -500, -30,
            -50, 0, -6, 0, 0, 0, -500, -50, 0, 0, -5, 0, -15, -1, 0, 0, 0, -20,
            0, 0, -19, 0, -150, 0, -35, -50, 0, -5, -37, 0, 0, -5, -7, 0, -26,
            0, 0, 0, -50, -16, 0, 0, 0, -500, -6, -18, -34, 0, -2, 0, -11, 0,
            -21, 0, 0, -50, 0, 0, -9, -24, 0, 0, -2, 0, 0, 0, 0, -8, -50, 0, 0,
            -26, 0, -10, 0, 0, 0, 0, 0, -2, 0, 0, 0, 0, -29, -49, -25, -10, 0,
            -2, -5, -2, -55, 0, 0, 0, 0, -17, 0, -2, -31, 0, -500, 0, -5, -43,
            0, -21, -50, 0, -150, 0, -500, -500, 0, 0, -13, -28, 0, -10, -11,
            0, -7, -3, -49, -17, 0, -8, 0, -17, 0, 0, -9, -9, 0, 0, -15, -14,
            0, 0, -38, 0, -50, -16, 0, 0, 0, -14, 0, 0, -34, -20, 0, 0, 0,
            -500, 0, -15, -5, 0, -13, -18, 0, -31, -16, 0, -52, 0, -29, -43, 0,
            0, 0, 0, -50, -199, 0, -8, -22, 0, -50, -51, 0, -19, -9, 0, 0, 0,
            0, 0, 0, 0, 0, 0, -15, 0, -17, -8, 0, -22, -8, 0, 0, 0, 0, -3, 0,
            -16, 0, 0, -10, -3, 0, 0, -3, 0, 0, 0, 0, 0, 0, -500, -500, -19,
            -50, -9, 0, 0, 0, 0, 0, -16, -4, 0, -5, 0, -7, -500, -9, -39, 0,
            -4, -30, 0, 0, 0, -20, 0, 0, -18, -24, 0, 0, -15, 0, -19, 0, 0, 0,
            0, 0, -33, 0, -27, -28, 0, 0, 0, -22, -11, 0, 0, -1, 0, -1, -7, -3,
            0, 0, 0, 0, 0, 0, 0, 0, -150, 0, 0, -150, -47, 0, 0, 0, -28, 0, -1,
            -43, 0, 0, 0, 0, 0, 0, -50, -11, -1, 0, -46, 0, 0, -9, 0, -500,
            -12, 0, 0, -150, 0, 0, -50, 0, -500, -11, -3, 0, 0, 0, 0, 0, -11,
            0, -2, 0, 0, -55, 0, 0, 0, -500, -150, -50, -30, 0, 0, -10, 0, 0,
            0, -8, 0, 0, 0, -35, -50, -50, 0, 0, 0, 0, 0, -150, 0, 0, -500, 0,
            0, 0, 0, 0, -7, 0, 0, -35};

        int32_t A_ciend_1[1366] = {
            45, 9, 1000, 1000, 1000, 50, 500, 204, 0, 500, 3, 1, 131, 23, 1, 7,
            61, 238, 3, 29, 261, 393, 50, 48, 182, 1000, 0, 32, 1, 33, 0, 201,
            12, 31, 20, 81, 150, 0, 0, 50, 50, 31, 82, 34, 19, 4, 110, 48, 4,
            50, 1, 48, 58, 150, 83, 45, 1000, 0, 1000, 236, 0, 222, 184, 0, 57,
            162, 0, 0, 268, 20, 0, 343, 283, 13, 500, 154, 336, 27, 0, 500,
            223, 311, 149, 50, 111, 198, 26, 20, 0, 8, 78, 51, 0, 10, 0, 80,
            25, 361, 240, 354, 3, 12, 1000, 0, 0, 50, 0, 7, 31, 2, 50, 1000,
            19, 1000, 305, 57, 60, 236, 45, 0, 27, 1, 36, 0, 3, 13, 10, 0, 135,
            13, 5, 3, 59, 4, 1, 500, 120, 32, 43, 284, 42, 0, 7, 2, 29, 268,
            1000, 182, 88, 285, 13, 26, 1000, 45, 500, 500, 0, 13, 1, 0, 2, 36,
            150, 150, 195, 26, 0, 62, 110, 11, 335, 41, 1, 0, 2, 0, 0, 24, 7,
            14, 34, 18, 0, 0, 230, 145, 0, 261, 0, 500, 0, 89, 15, 48, 0, 28,
            6, 59, 2, 0, 1000, 0, 354, 27, 27, 500, 0, 0, 500, 266, 22, 500, 8,
            0, 299, 14, 268, 0, 0, 0, 94, 0, 14, 227, 0, 1, 14, 0, 15, 18, 500,
            1000, 0, 0, 251, 26, 7, 500, 166, 48, 0, 0, 17, 39, 0, 248, 0, 25,
            41, 271, 0, 240, 403, 33, 274, 0, 57, 23, 19, 10, 6, 56, 0, 25, 27,
            154, 146, 0, 34, 9, 0, 328, 54, 12, 0, 0, 330, 0, 1, 150, 0, 4,
            173, 48, 267, 50, 0, 0, 41, 15, 270, 19, 0, 0, 210, 0, 0, 0, 25,
            50, 11, 47, 3, 36, 150, 3, 295, 27, 105, 11, 150, 1, 1000, 6, 50,
            79, 1000, 1, 111, 39, 0, 0, 0, 11, 13, 500, 0, 3, 31, 21, 6, 0, 0,
            0, 50, 78, 0, 500, 4, 0, 2, 17, 0, 0, 336, 97, 22, 1000, 40, 149,
            1000, 46, 17, 0, 6, 82, 194, 29, 18, 0, 0, 32, 32, 500, 50, 500,
            500, 138, 87, 500, 0, 500, 0, 14, 13, 0, 213, 0, 9, 500, 10, 0,
            301, 2, 156, 14, 1, 150, 0, 500, 0, 500, 68, 0, 0, 16, 237, 324,
            98, 6, 25, 0, 0, 0, 205, 76, 311, 39, 1, 0, 0, 32, 1000, 36, 54,
            62, 290, 58, 2, 150, 42, 0, 0, 424, 0, 17, 0, 54, 9, 0, 2, 84, 0,
            50, 500, 2, 500, 150, 500, 15, 257, 0, 0, 0, 118, 232, 0, 1000, 0,
            7, 0, 211, 33, 28, 2, 159, 198, 49, 4, 2, 500, 0, 7, 248, 50, 8, 0,
            157, 315, 0, 31, 0, 1, 94, 41, 0, 0, 0, 1000, 0, 0, 0, 4, 500, 360,
            414, 163, 29, 1, 500, 500, 46, 50, 0, 19, 8, 50, 70, 0, 13, 373,
            253, 9, 1000, 1, 1000, 500, 12, 211, 114, 426, 88, 0, 500, 0, 26,
            15, 106, 1000, 302, 1, 284, 1, 0, 500, 1000, 150, 0, 28, 0, 500, 0,
            241, 41, 12, 242, 253, 0, 34, 76, 0, 0, 12, 4, 0, 0, 13, 31, 0, 10,
            0, 0, 13, 0, 196, 26, 187, 14, 0, 400, 17, 217, 0, 14, 50, 0, 1000,
            18, 31, 50, 3, 1000, 8, 50, 16, 0, 0, 144, 55, 7, 4, 50, 0, 500,
            1000, 0, 3, 21, 247, 50, 17, 8, 19, 21, 26, 280, 31, 24, 1000, 2,
            1000, 409, 318, 292, 13, 303, 500, 500, 50, 500, 0, 9, 83, 8, 0,
            15, 500, 21, 500, 0, 350, 0, 0, 50, 385, 45, 46, 485, 24, 4, 0, 0,
            13, 259, 10, 4, 1, 1000, 22, 0, 255, 205, 34, 3, 3, 159, 0, 221,
            35, 14, 1000, 1, 500, 11, 51, 500, 123, 0, 11, 0, 23, 11, 44, 1,
            275, 1000, 3, 26, 1000, 52, 5, 0, 21, 257, 500, 58, 253, 147, 207,
            66, 0, 18, 17, 319, 297, 71, 6, 0, 262, 500, 0, 238, 0, 160, 0,
            1000, 0, 98, 261, 0, 500, 0, 282, 13, 0, 500, 1000, 23, 2, 1, 500,
            1000, 246, 0, 5, 319, 33, 11, 0, 0, 1000, 0, 413, 211, 37, 0, 23,
            500, 0, 28, 50, 8, 14, 0, 18, 41, 0, 0, 500, 7, 0, 304, 135, 0, 0,
            0, 0, 0, 233, 198, 50, 93, 0, 448, 0, 333, 1000, 367, 0, 375, 16,
            111, 23, 212, 343, 299, 23, 12, 0, 0, 0, 21, 192, 1, 117, 0, 0, 0,
            207, 0, 25, 500, 4, 94, 500, 1, 5, 0, 0, 500, 9, 2, 1, 5, 5, 0,
            500, 12, 0, 0, 0, 23, 0, 0, 0, 5, 229, 46, 1, 401, 23, 0, 369, 263,
            104, 257, 1, 1, 143, 0, 240, 0, 320, 0, 0, 32, 26, 0, 0, 23, 0,
            243, 0, 0, 0, 2, 0, 0, 0, 500, 106, 0, 241, 219, 0, 150, 9, 427, 0,
            150, 0, 26, 0, 244, 500, 0, 6, 166, 32, 0, 52, 0, 0, 0, 500, 0, 0,
            0, 500, 57, 500, 1, 0, 186, 500, 0, 215, 22, 0, 280, 30, 0, 0, 3,
            263, 0, 26, 0, 299, 0, 22, 0, 25, 11, 66, 226, 500, 25, 0, 39, 0,
            528, 3, 21, 18, 0, 6, 79, 220, 14, 0, 1000, 500, 500, 0, 500, 52,
            500, 0, 1000, 304, 50, 0, 44, 1, 29, 59, 0, 3, 1000, 0, 0, 0, 0,
            250, 188, 1000, 176, 12, 500, 6, 20, 1000, 50, 3, 0, 0, 30, 7, 23,
            66, 1000, 1, 213, 0, 18, 0, 1, 150, 9, 218, 114, 21, 7, 31, 0, 61,
            126, 0, 500, 500, 19, 1, 17, 27, 0, 167, 0, 31, 212, 1000, 18, 0,
            500, 0, 21, 0, 0, 30, 11, 69, 1000, 30, 50, 92, 7, 154, 218, 500,
            1000, 50, 457, 0, 5, 0, 16, 2, 500, 0, 286, 20, 0, 0, 19, 500, 150,
            0, 36, 84, 0, 6, 37, 0, 503, 5, 7, 0, 26, 500, 500, 0, 96, 17, 0,
            500, 0, 1000, 6, 18, 34, 0, 2, 0, 11, 326, 21, 239, 0, 50, 126,
            281, 10, 24, 500, 193, 3, 0, 200, 316, 424, 8, 65, 319, 269, 26,
            172, 10, 0, 0, 0, 0, 500, 3, 0, 128, 0, 0, 30, 50, 26, 10, 0, 3, 5,
            2, 55, 0, 220, 0, 1, 17, 0, 3, 31, 0, 1000, 234, 6, 43, 218, 21,
            50, 0, 150, 111, 1000, 1000, 0, 220, 13, 28, 137, 10, 11, 274, 8,
            3, 49, 17, 0, 8, 39, 18, 0, 257, 9, 9, 108, 230, 16, 15, 500, 70,
            39, 500, 77, 16, 305, 0, 394, 15, 95, 500, 35, 20, 24, 0, 235,
            1000, 0, 16, 6, 302, 14, 19, 203, 31, 17, 0, 52, 0, 29, 43, 215, 0,
            0, 500, 98, 89, 0, 9, 23, 302, 118, 51, 236, 20, 9, 85, 500, 500,
            0, 0, 0, 0, 546, 1, 15, 500, 17, 8, 198, 22, 9, 199, 0, 0, 317, 3,
            0, 17, 0, 0, 10, 4, 265, 0, 3, 203, 0, 500, 0, 0, 0, 1000, 1000,
            20, 250, 9, 336, 500, 0, 201, 309, 16, 5, 1, 5, 259, 8, 1000, 10,
            40, 0, 4, 31, 0, 0, 0, 21, 0, 172, 18, 24, 500, 0, 16, 193, 20, 0,
            139, 58, 500, 0, 34, 312, 27, 29, 0, 0, 0, 22, 12, 274, 38, 2, 381,
            1, 7, 4, 500, 281, 0, 464, 0, 209, 0, 0, 150, 0, 364, 150, 47, 0,
            0, 0, 28, 0, 2, 44, 0, 272, 184, 0, 165, 179, 96, 12, 2, 0, 47,
            500, 0, 10, 0, 1000, 13, 114, 0, 150, 0, 0, 100, 0, 1000, 12, 3, 0,
            223, 0, 0, 0, 12, 212, 3, 0, 0, 55, 0, 0, 268, 1000, 150, 50, 31,
            500, 255, 10, 0, 0, 0, 9, 0, 0, 0, 35, 66, 50, 315, 500, 0, 500,
            137, 150, 0, 177, 1000, 500, 0, 0, 0, 84, 7, 0, 0, 36};
//}}}

    uint32_t num_sv_dels = 0;
    uint32_t i = 0;
    enum stix_sv_type sv_type;
    while (bcf_read(fp, hdr, line) == 0) {
        if (stix_get_vcf_breakpoints(fp,
                                     hdr,
                                     line,
                                     left,
                                     right,
                                     &sv_type) == 0) {
            if (sv_type == DEL) {
                num_sv_dels+=1;
                TEST_ASSERT_EQUAL(0, strcmp("13", left->chrm));
                TEST_ASSERT_EQUAL(0, strcmp("13", right->chrm));
                TEST_ASSERT_EQUAL(A_pos[i] - 1 + A_cipos_0[i], left->start);
                TEST_ASSERT_EQUAL(A_pos[i] - 1 + A_cipos_1[i], left->end);
                TEST_ASSERT_EQUAL(A_end[i]  + A_ciend_0[i], right->start);
                TEST_ASSERT_EQUAL(A_end[i]  + A_ciend_1[i], right->end);
                i+=1;
            }
        }
    }

    /*
     * bcftools view -H -i 'SVTYPE="DEL"'
     * ALL.wgs.integrated_sv_map_v2.20130502.svs.genotypes.13.vcf.gz | wc -l
     * 1366
     */
    TEST_ASSERT_EQUAL(1366, num_sv_dels);
}
//}}}
