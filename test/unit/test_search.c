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

//{{{void test_stix_check_del(void)
void test_stix_check_del(void)
{
    struct stix_breakpoint q_l,
                           q_r,
                           i_l,
                           i_r;
    uint32_t evidence_type;
    uint32_t slop = 100;

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
                                      5);



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
