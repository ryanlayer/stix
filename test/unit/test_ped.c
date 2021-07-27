#ifndef _GNU_SOURCE
#define _GNU_SOURCE 1
#endif
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <inttypes.h>
#include <string.h>
#include <sqlite3.h>

#include "unity.h"

#include <src/giggle_index.h>
#include <src/util.h>

#include "ped.h"

void setUp(void) { }
void tearDown(void) { }

//{{{void test_is_int(void)
void test_is_int(void)
{

    uint32_t i;
    for (i = 0; i < 10000; ++i) {
        char *rand_s = NULL;
        int rand_i = rand();
        int ret = asprintf(&rand_s, "%d", rand_i);
        int v;
        TEST_ASSERT_EQUAL(1, is_int(rand_s, &v));
        TEST_ASSERT_EQUAL(rand_i, v);
        free(rand_s);
    }

    for (i = 0; i < 10000; ++i) {
        char *rand_s = NULL;
        int rand_i = rand();
        int ret = asprintf(&rand_s, "%dX", rand_i);
        int v;
        TEST_ASSERT_EQUAL(0, is_int(rand_s, &v));
        free(rand_s);
    }

    for (i = 0; i < 10000; ++i) {
        char *rand_s = NULL;
        int rand_i = rand();
        int ret = asprintf(&rand_s, "X%d", rand_i);
        int v;
        TEST_ASSERT_EQUAL(0, is_int(rand_s, &v));
        free(rand_s);
    }

    for (i = 0; i < 10000; ++i) {
        char *rand_s = NULL;
        int rand_i = rand();
        int ret = asprintf(&rand_s, "X%dX", rand_i);
        int v;
        TEST_ASSERT_EQUAL(0, is_int(rand_s, &v));
        free(rand_s);
    }

    int v;
    TEST_ASSERT_EQUAL(0, is_int("134a", &v));
    TEST_ASSERT_EQUAL(0, is_int("134a", &v));

    TEST_ASSERT_EQUAL(0, is_int("134.3", &v));

    TEST_ASSERT_EQUAL(0, is_int(" 13 ", &v));

    // overflow
    TEST_ASSERT_EQUAL(0, is_int("928374238742387429340234", &v));
    TEST_ASSERT_EQUAL(0, is_int("2147483648", &v));


    TEST_ASSERT_EQUAL(1, is_int("13", &v));
    TEST_ASSERT_EQUAL(13, v);

    TEST_ASSERT_EQUAL(1, is_int(" 14", &v));
    TEST_ASSERT_EQUAL(14, v);

    TEST_ASSERT_EQUAL(1, is_int("2147483647", &v));
    TEST_ASSERT_EQUAL(INT_MAX, v);
}
//}}}

//{{{ void test_check_field_name(void)
void test_check_field_name(void)
{
    // cannot start with a number
    TEST_ASSERT_EQUAL(0, check_field_name("123A"));

    // Test the invlid ranges
    int i;
    char test_str[3] = "A B";
    for (i = 32; i < 48; ++i) {
        test_str[1] = i;
        TEST_ASSERT_EQUAL(1, check_field_name(test_str));
    }

    for (i = 58; i < 65; ++i) {
        test_str[1] = i;
        TEST_ASSERT_EQUAL(1, check_field_name(test_str));
    }

    for (i = 91; i < 96; ++i) {
        if (i == '_')
            continue;
        test_str[1] = i;
        TEST_ASSERT_EQUAL(1, check_field_name(test_str));
    }

    for (i = 173; i < 128; ++i) {
        test_str[1] = i;
        TEST_ASSERT_EQUAL(1, check_field_name(test_str));
    }

    TEST_ASSERT_EQUAL(-1, check_field_name("A"));
    TEST_ASSERT_EQUAL(-1, check_field_name("A234"));
    TEST_ASSERT_EQUAL(-1, check_field_name("A_234"));

    char nums[10] = {'0','1','2','3','4','5','6','7','8','9'};

    char *starts_with_num = strdup("Xtest");
    for (i = 0; i < 10; ++i) {
        starts_with_num[0] = nums[i];
        TEST_ASSERT_EQUAL(0, check_field_name(starts_with_num));
    }

    char bad_chars[11] = {'[',']','^','`',':',';','<','=','>','?','@'};
    char *has_bad_char = strdup("testXtest");
    for (i = 0; i < 10; ++i) {
        has_bad_char[4] = bad_chars[i];
        TEST_ASSERT_EQUAL(4, check_field_name(has_bad_char));
    }

    char *has_good_char = strdup("testXtest");
    for (i = 0; i < 10; ++i) {
        has_good_char[4] = nums[i];
        TEST_ASSERT_EQUAL(-1, check_field_name(has_good_char));
    }

    has_good_char[4] = '_';
    TEST_ASSERT_EQUAL(-1, check_field_name(has_good_char));


    free(starts_with_num);
    free(has_bad_char);
    free(has_good_char);
}
//}}}

//{{{ void test_ped_get_column_names(void)
void test_ped_get_column_names(void)
{
    char *ped_file_name = "../data/four.ped";
    char **ped_field_names = NULL;
    int *ped_field_is_int = NULL;
    int num_ped_fields = ped_get_column_names_types(ped_file_name,
                                                    &ped_field_names,
                                                    &ped_field_is_int);

    char *A_ped_field_names[5] = { "Sample",
                                   "Sex",
                                   "Population",
                                   "Super_Population",
                                   "Alt_File" };
    int A_ped_field_is_int[5] = { 0, 1, 0, 0, 0};

    TEST_ASSERT_EQUAL(5, num_ped_fields);

    uint32_t i;
    for (i = 0; i < num_ped_fields; ++i) {
        TEST_ASSERT_EQUAL(0,
                          strcmp(A_ped_field_names[i],
                                 ped_field_names[i]));
        TEST_ASSERT_EQUAL(A_ped_field_is_int[i], ped_field_is_int[i]);
    }

    for (i = 0; i < num_ped_fields; ++i) 
        free(ped_field_names[i]);

    free(ped_field_names);
    free(ped_field_is_int);
}
///}}}

//{{{ void test_ped_create_db(void)
void test_ped_create_db(void)
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

    TEST_ASSERT_EQUAL(4, num_rows);
    remove(ped_file_name_db);
    rmrf(output_path_name);
}
//}}}

//{{{void test_ped_get_matching_sample_ids(void)
void test_ped_get_matching_sample_ids(void)
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

    char *select_query_1 = "Sex == 1";
    uint32_t *sample_ids = NULL;
    uint32_t num_samples = ped_get_matching_sample_ids(ped_file_name_db,
                                                       select_query_1,
                                                       &sample_ids);
    TEST_ASSERT_EQUAL(2, num_samples);
    TEST_ASSERT_EQUAL(1, sample_ids[0]);
    TEST_ASSERT_EQUAL(2, sample_ids[1]);

    free(sample_ids);

    char *select_query_2 = "(Sex == 2) AND (Population == 'CEU')";
    num_samples = ped_get_matching_sample_ids(ped_file_name_db,
                                              select_query_2,
                                              &sample_ids);
    TEST_ASSERT_EQUAL(1, num_samples);
    TEST_ASSERT_EQUAL(3, sample_ids[0]);

    free(sample_ids);

    char *select_query_3 = "(Sex == 3)";
    num_samples = ped_get_matching_sample_ids(ped_file_name_db,
                                              select_query_3,
                                              &sample_ids);
    TEST_ASSERT_EQUAL(0, num_samples);

    remove(ped_file_name_db);
    rmrf(output_path_name);
}    
//}}}

//{{{void test_ped_get_uniq_col_groups
void test_ped_get_uniq_col_groups(void)
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
 
    TEST_ASSERT_EQUAL(4, num_uniq_vals);
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
        "SELECT Giggle_File_Id FROM ped WHERE Sex=='1' AND Population=='CEU'"
     2
     sqlite3 four.ped.db \
        "SELECT Giggle_File_Id FROM ped WHERE Sex=='1' AND Population=='CHS'"
     1
     sqlite3 four.ped.db \
        "SELECT Giggle_File_Id FROM ped WHERE Sex=='2' AND Population=='CEU'"
     3
     sqlite3 four.ped.db \
        "SELECT Giggle_File_Id FROM ped WHERE Sex=='2' AND Population=='CHS'"
     0
     */
    TEST_ASSERT_EQUAL(1, uniq_groups_sizes[0]);
    TEST_ASSERT_EQUAL(2, uniq_groups_ids[0][0]);

    TEST_ASSERT_EQUAL(1, uniq_groups_sizes[1]);
    TEST_ASSERT_EQUAL(1, uniq_groups_ids[1][0]);

    TEST_ASSERT_EQUAL(1, uniq_groups_sizes[2]);
    TEST_ASSERT_EQUAL(3, uniq_groups_ids[2][0]);

    TEST_ASSERT_EQUAL(1, uniq_groups_sizes[3]);
    TEST_ASSERT_EQUAL(0, uniq_groups_ids[3][0]);

    /* 
     sqlite3 four.ped.db "SELECT Population FROM ped GROUP BY Population"
     CEU
     CHS
     */

    char *cols_1[1] = {"Population"};
    num_uniq_vals = ped_get_uniq_col_groups(ped_file_name_db,
                                            &db,
                                            cols_1,
                                            1,
                                            NULL,
                                            &uniq_vals,
                                            &uniq_groups_ids,
                                            &uniq_groups_sizes);

    TEST_ASSERT_EQUAL(2, num_uniq_vals);
    TEST_ASSERT_EQUAL(0, strcmp("CEU",   uniq_vals[0][0]));
    TEST_ASSERT_EQUAL(0, strcmp("CHS",   uniq_vals[1][0]));

    /*
     SELECT Giggle_File_Id FROM ped WHERE Population=='CEU'
     2
     3
     SELECT Giggle_File_Id FROM ped WHERE Population=='CHS'
     0
     1
     */
    TEST_ASSERT_EQUAL(2, uniq_groups_sizes[0]);
    TEST_ASSERT_EQUAL(2, uniq_groups_ids[0][0]);
    TEST_ASSERT_EQUAL(3, uniq_groups_ids[0][1]);

    TEST_ASSERT_EQUAL(2, uniq_groups_sizes[1]);
    TEST_ASSERT_EQUAL(0, uniq_groups_ids[1][0]);
    TEST_ASSERT_EQUAL(1, uniq_groups_ids[1][1]);

    /*
     sqlite3 four.ped.db \
        "SELECT Giggle_File_Id FROM ped WHERE Population=='CEU' AND Sex==1"
     2
     sqlite3 four.ped.db \
        "SELECT Giggle_File_Id FROM ped WHERE Population=='CHS' AND Sex==1"
     1
     */
     
    char *select_q = "Sex == 1";
    num_uniq_vals = ped_get_uniq_col_groups(ped_file_name_db,
                                            &db,
                                            cols_1,
                                            1,
                                            select_q,
                                            &uniq_vals,
                                            &uniq_groups_ids,
                                            &uniq_groups_sizes);

    TEST_ASSERT_EQUAL(2, num_uniq_vals);
    TEST_ASSERT_EQUAL(0, strcmp("CEU", uniq_vals[0][0]));
    TEST_ASSERT_EQUAL(1, uniq_groups_sizes[0]);
    TEST_ASSERT_EQUAL(2, uniq_groups_ids[0][0]);

    TEST_ASSERT_EQUAL(0, strcmp("CHS", uniq_vals[1][0]));
    TEST_ASSERT_EQUAL(1, uniq_groups_sizes[1]);
    TEST_ASSERT_EQUAL(1, uniq_groups_ids[1][0]);

    remove(ped_file_name_db);
    rmrf(output_path_name);
    sqlite3_close(db);
}    
//}}}

//{{{void test_ped_union_groups(void)
void test_ped_union_groups(void)
{

    uint32_t g0_sizes[3] = {2, 1, 4};
    uint32_t *g0[3];
    g0[0] = (uint32_t *)malloc(g0_sizes[0]*sizeof(uint32_t));
    g0[1] = (uint32_t *)malloc(g0_sizes[1]*sizeof(uint32_t));
    g0[2] = (uint32_t *)malloc(g0_sizes[2]*sizeof(uint32_t));

    g0[0][0] = 10;
    g0[0][1] = 11;

    g0[1][0] = 20;

    g0[2][0] = 30;
    g0[2][1] = 31;
    g0[2][2] = 32;
    g0[2][3] = 33;

    uint32_t *g0_union;
    uint32_t g0_total_size = ped_union_groups(3, g0, g0_sizes, &g0_union);

    TEST_ASSERT_EQUAL(7, g0_total_size);

    uint32_t A_g0_union[7] = { 10, 11, 20, 30, 31, 32, 33 };
    uint32_t i;
    for (i = 0; i < g0_total_size; ++i)
        TEST_ASSERT_EQUAL(A_g0_union[i], g0_union[i]);

    free(g0[0]);
    free(g0[1]);
    free(g0[2]);
}
//}}}

//{{{void test_ped_get_matching_sample_ids(void)
void test_ped_get_cols_info_by_id(void)
{
    char *input_path_name = "../data/four_alt_sort/*gz";
    char *output_path_name = "tmp_ped_get_cols_info_by_id";

    uint64_t indexed_intervals = giggle_bulk_insert(input_path_name,
                                                    output_path_name,
                                                    1);

    char *ped_file_name = "../data/four.ped";
    char *ped_file_name_db = "four.ped.db";

    uint32_t num_rows = ped_create_db(ped_file_name,
                                      ped_file_name_db,
                                      output_path_name,
                                      5);

    char *cols[2] = {"Sex", "Population"};
    char **col_vals = NULL, **col_names = NULL;
    sqlite3 *db = NULL;
    uint32_t num_col_vals = ped_get_cols_info_by_id(ped_file_name_db,
                                                    &db,
                                                    cols,
                                                    2,
                                                    1,
                                                    &col_vals,
                                                    &col_names);
    /*
     * sqlite3 ../data/four.ped.db "select *from ped"
     * 2|NA12812|1|CEU|CEU Utah Residents (CEPH) with Northern and Western European Ancestry|EUR|NA12812.13.14.bed.gz
     * 0|HG00672|2|CHS|Southern Han Chinese|EAS|HG00672.13.14.bed.gz
     * 3|NA12878|2|CEU|CEU Utah Residents (CEPH) with Northern and Western European Ancestry|EUR|NA12878.13.14.bed.gz
     * 1|HG00674|1|CHS|Southern Han Chinese|EAS|HG00674.13.14.bed.gz
     */

    TEST_ASSERT_EQUAL(2, num_col_vals);
    TEST_ASSERT_EQUAL(0, strcmp("1", col_vals[0]));
    TEST_ASSERT_EQUAL(0, strcmp("CHS", col_vals[1]));

    TEST_ASSERT_EQUAL(0, strcmp("Sex", col_names[0]));
    TEST_ASSERT_EQUAL(0, strcmp("Population", col_names[1]));

    free(col_vals[0]);
    free(col_vals[1]);
    free(col_vals);

    free(col_names[0]);
    free(col_names[1]);
    free(col_names);

    char *cols_1[3] = {"Sex", "Population", "Super_Population"};
    num_col_vals = ped_get_cols_info_by_id(ped_file_name_db,
                                           &db,
                                           cols_1,
                                           3,
                                           3,
                                           &col_vals,
                                           &col_names);
    TEST_ASSERT_EQUAL(3, num_col_vals);
    TEST_ASSERT_EQUAL(0, strcmp("2", col_vals[0]));
    TEST_ASSERT_EQUAL(0, strcmp("CEU", col_vals[1]));
    TEST_ASSERT_EQUAL(0, strcmp("EUR", col_vals[2]));

    TEST_ASSERT_EQUAL(0, strcmp("Sex", col_names[0]));
    TEST_ASSERT_EQUAL(0, strcmp("Population", col_names[1]));
    TEST_ASSERT_EQUAL(0, strcmp("Super_Population", col_names[2]));

    free(col_vals[0]);
    free(col_vals[1]);
    free(col_vals[2]);
    free(col_vals);

    free(col_names[0]);
    free(col_names[1]);
    free(col_names[2]);
    free(col_names);

    num_col_vals = ped_get_cols_info_by_id(ped_file_name_db,
                                           &db,
                                           NULL,
                                           0,
                                           1,
                                           &col_vals,
                                           &col_names);
    TEST_ASSERT_EQUAL(6, num_col_vals);
    TEST_ASSERT_EQUAL(0, strcmp("Giggle_File_Id", col_names[0]));
    TEST_ASSERT_EQUAL(0, strcmp("Sample", col_names[1]));
    TEST_ASSERT_EQUAL(0, strcmp("Sex", col_names[2]));
    TEST_ASSERT_EQUAL(0, strcmp("Population", col_names[3]));
    TEST_ASSERT_EQUAL(0, strcmp("Super_Population", col_names[4]));
    TEST_ASSERT_EQUAL(0, strcmp("Alt_File", col_names[5]));

    sqlite3_close(db);
    rmrf(output_path_name);
    remove(ped_file_name_db);
}    
//}}}
