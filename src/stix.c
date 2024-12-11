#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <sysexits.h>
#include <err.h>

#include <sqlite/sqlite3.h>

#include <htslib/hfile.h>
#include <htslib/vcf.h>

#include "ped.h"
#include "search.h"
#include "sharding_utils.h"

char *stix_sv_type_strings[5] = {"DEL", "DUP", "INS", "INV", "BND"};

uint32_t parse_aggregate_csv(char *aggregate,
                             char ***agg_cols);

void print_results(struct giggle_index *gi,
                   char *ped_db_file_name,
                   struct stix_breakpoint *left,
                   struct stix_breakpoint *right,
                   uint32_t *sample_ids,
                   char *filter,
                   struct uint_pair *sample_alt_depths,
                   uint32_t num_samples,
                   char **agg_cols,
                   uint32_t num_agg_cols,
                   uint32_t json_out,
                   uint32_t summary_only,
                   uint32_t depths_only);

void update_vcf_header(bcf_hdr_t *hdr,
                       uint32_t v_is_set,
                       char *sample_column);

//{{{int help(int exit_code)
int help(int exit_code)
{
    fprintf(stderr,
            "STIX: Structural Variant Index\n"
            "Version: 2.0.0\n"
            "usage:   stix <options>\n"
            "         options:\n"
            "             -i  index directory\n"
            "             -s  slop\n"
            "             -P  padding base piars for query insertion(default 50)\n"
            "             -p  PED file\n"
            "             -B  Sharding file (giggle_index \\t STIX_index)\n"
            "             -Q  Batch queries in table format (left \\t right \\t len \\t SVtype \\t ID)\n"
            "             -T  Threshold to define STIX_ONE (number of supporting reads required to define a hit) (default: 1)\n"
            "             -c  Alt file column (default 1)\n"
            "             -d  PED database file\n"
            "             -r  right SV region\n"
            "             -l  left SV region\n"
            "             -f  VCF file\n"
            "             -a  List of columns to aggregate over\n"
            "             -F  Filter samples by PED field query\n"
            "             -j  JSON output\n"
            "             -t  SV type (DEL,INS,INV,DUP,BND) \n"
            "             -L  Length of Insertion(use it when -t INS) \n"
            "             -R  Relative Erorr Threshold to compare the length of query INS and targeted INS (0.0-1.0) (default:0.15) \n"
            "             -v  Add sample depth to VCF file\n"
            "             -S  Give only summary\n"
            "             -D  Give sample depth array\n"
            "             -V  Verbose mode(print debug information, will increase output file size greatly)\n");
    return exit_code;
}
//}}}

//{{{void print_results(struct giggle_index *gi,
void print_results(struct giggle_index *gi,
                   char *ped_db_file_name,
                   struct stix_breakpoint *left,
                   struct stix_breakpoint *right,
                   uint32_t *sample_ids,
                   char *filter,
                   struct uint_pair *sample_alt_depths,
                   uint32_t num_samples,
                   char **agg_cols,
                   uint32_t num_agg_cols,
                   uint32_t json_out,
                   uint32_t summary_only,
                   uint32_t depths_only)
{
    uint32_t i; 
    uint32_t num_col_vals;
    int32_t zero_count, one_count;
    uint32_t Q1, Q2, Q3, min, max;
    int32_t counts[4];
    uint32_t ret = stix_get_summary(sample_alt_depths,
                                    sample_ids,
                                    num_samples,
                                    &zero_count,
                                    &one_count,
                                    &Q1,
                                    &Q2,
                                    &Q3,
                                    &min,
                                    &max,
                                    counts);

    if (json_out == 1)
        printf("{ \"results\": {\n");

    if (json_out == 1)
    {
        printf("\"summary\": [\n");

        printf("{"
               "\"name\":\"Total\", "
               "\"zero_count\":\"%d\", "
               "\"one_count\":\"%d\","
               "\"quantiles\":[\"%u\",\"%u\",\"%u\"], "
               "\"counts\":[\"%u\",\"%u\",\"%u\",\"%u\"]",
               zero_count,
               one_count,
               Q1, Q2, Q3,
               counts[0], counts[1], counts[2], counts[3]);
    }
    else
    {
        printf("Total\t0:1\t%d:%d\t%u:%u:%u\t%d:%d:%d:%d",
               zero_count, one_count,
               Q1, Q2, Q3,
               counts[0], counts[1], counts[2], counts[3]);
    }

    if (depths_only == 1)
    {
        uint32_t *sample_depths = NULL;
        ret = stix_get_sample_depths(sample_alt_depths,
                                     NULL,
                                     num_samples,
                                     &sample_depths);
        if (json_out == 0)
        {
            for (i = 0; i < num_samples; ++i)
                printf("\t%d", sample_depths[i]);
        }
        else
        {
            printf(",\"depths\":[");
            for (i = 0; i < num_samples; ++i)
            {
                if (i != 0)
                    printf(",");
                printf("%d", sample_depths[i]);
            }
            printf("]");
        }

        free(sample_depths);
    }

    if (json_out == 1)
        printf("}");
    printf("\n");

    sqlite3 *db = NULL;

    if (num_agg_cols > 0)
    {
        char ***uniq_vals;
        uint32_t num_uniq_vals;
        uint32_t **uniq_groups_ids;
        uint32_t *uniq_groups_sizes;
        num_uniq_vals = ped_get_uniq_col_groups(ped_db_file_name,
                                                &db,
                                                agg_cols,
                                                num_agg_cols,
                                                filter,
                                                &uniq_vals,
                                                &uniq_groups_ids,
                                                &uniq_groups_sizes);

        for (i = 0; i < num_uniq_vals; ++i)
        {
            uint32_t j;
            char *group_name_tmp, *group_name;
            ret = asprintf(&group_name, "%s", uniq_vals[i][0]);
            for (j = 1; j < num_agg_cols; ++j)
            {
                ret = asprintf(&group_name_tmp,
                               "%s,%s",
                               group_name,
                               uniq_vals[i][j]);
                free(group_name);
                group_name = group_name_tmp;
            }

            ret = stix_get_summary(sample_alt_depths,
                                   uniq_groups_ids[i],
                                   uniq_groups_sizes[i],
                                   &zero_count,
                                   &one_count,
                                   &Q1,
                                   &Q2,
                                   &Q3,
                                   &min,
                                   &max,
                                   counts);

            if (json_out)
            {
                printf(",{"
                       "\"name\":\"%s\", "
                       "\"zero_count\":\"%d\", "
                       "\"one_count\":\"%d\","
                       "\"quantiles\":[\"%u\",\"%u\",\"%u\"], "
                       "\"counts\":[\"%u\",\"%u\",\"%u\",\"%u\"]",
                       group_name,
                       zero_count,
                       one_count,
                       Q1, Q2, Q3,
                       counts[0], counts[1], counts[2], counts[3]);
            }
            else
            {
                printf("%s\t0:1\t%d:%d\t%u:%u:%u\t%d:%d:%d:%d",
                       group_name,
                       zero_count, one_count,
                       Q1, Q2, Q3,
                       counts[0], counts[1], counts[2], counts[3]);
            }

            if (depths_only == 1)
            {
                uint32_t *sample_depths = NULL;
                ret = stix_get_sample_depths(sample_alt_depths,
                                             uniq_groups_ids[i],
                                             uniq_groups_sizes[i],
                                             &sample_depths);
                int k;
                if (json_out == 0)
                {
                    for (k = 0; k < uniq_groups_sizes[i]; ++k)
                        printf("\t%d", sample_depths[k]);
                }
                else
                {
                    printf(",\"depths\":[");
                    for (k = 0; k < uniq_groups_sizes[i]; ++k)
                    {
                        if (k != 0)
                            printf(",");
                        printf("%d", sample_depths[k]);
                    }
                    printf("]");
                }
                free(sample_depths);
            }

            if (json_out == 1)
                printf("}");
            printf("\n");

            free(group_name);
        }
    }

    if (summary_only == 1)
    {
        if (json_out == 1)
            printf("]}}");
        return;
    }

    if (depths_only == 0)
    {
        if (json_out)
            printf("],\n\"samples\": [\n");

        char **col_vals = NULL, **col_names = NULL;
        // printf("num_samples:%d\n",num_samples);
        for (i = 0; i < num_samples; ++i)
        {
            uint32_t idx;
            if (sample_ids != NULL)
                idx = sample_ids[i];
            else
                idx = i;

            // printf("ped_get_cols_info_by_id:%d \n",idx);
            num_col_vals = ped_get_cols_info_by_id(ped_db_file_name,
                                                   &db,
                                                   NULL,
                                                   0,
                                                   idx,
                                                   &col_vals,
                                                   &col_names);

            uint32_t j;
            if (json_out == 1)
            {
                if (i > 0)
                    printf(",");

                printf("{");
                for (j = 0; j < num_col_vals; ++j)
                {
                    if (j > 0)
                        printf(",");
                    printf("\"%s\":\"%s\"", col_names[j], col_vals[j]);
                }
                printf(",\"Pairend\":\"%u\",\"Split\":\"%u\"}\n",
                       sample_alt_depths[i].first,
                       sample_alt_depths[i].second);
            }
            else
            {
                if (i == 0)
                {
                    for (j = 0; j < num_col_vals; ++j)
                        printf("%s\t", col_names[j]);
                    printf("Pairend\tSplit\n");
                }

                for (j = 0; j < num_col_vals; ++j)
                    printf("%s\t", col_vals[j]);

                printf("%u\t%u\n",
                       sample_alt_depths[i].first,
                       sample_alt_depths[i].second);
            }
        }
    }

    if (json_out == 1)
        printf("]}}");

    sqlite3_close(db);
}
//}}}

//{{{void print_results_shard(struct giggle_index *gi,
void print_results_shard(
    // struct giggle_index *gi_all[],
    //    char *ped_db_file_name,
    Shard shard_array[],
    struct stix_breakpoint *left,
    struct stix_breakpoint *right,
    uint32_t *sample_ids_all[],
    char *filter,
    struct uint_pair *sample_alt_depths_all[],
    uint32_t num_samples_all[],
    uint32_t size_num_samples_all,
    char **agg_cols,
    uint32_t num_agg_cols,
    uint32_t json_out,
    uint32_t summary_only,
    uint32_t depths_only)
{
    uint32_t i;

    uint32_t num_col_vals;
    int32_t zero_count, one_count;
    uint32_t Q1, Q2, Q3, min, max;
    int32_t counts[4];
    // uint32_t size_num_samples_all =  sizeof(num_samples_all) / sizeof(num_samples_all[0]);

    // if(V_is_set ==1)
    //     fprintf(stderr,"<<<<<<<<<<num_samples_all[0]:%d\tnum_samples_all[1]:%d,",num_samples_all[0],num_samples_all[1]);
    uint32_t ret = stix_get_summary_shard(sample_alt_depths_all,
                                          sample_ids_all,
                                          num_samples_all,
                                          size_num_samples_all,
                                          &zero_count,
                                          &one_count,
                                          &Q1,
                                          &Q2,
                                          &Q3,
                                          &min,
                                          &max,
                                          counts);

    if (json_out == 1)
        printf("{ \"results\": {\n");

    if (json_out == 1)
    {
        printf("\"summary\": [\n");

        printf("{"
               "\"name\":\"Total\", "
               "\"zero_count\":\"%d\", "
               "\"one_count\":\"%d\","
               "\"quantiles\":[\"%u\",\"%u\",\"%u\"], "
               "\"counts\":[\"%u\",\"%u\",\"%u\",\"%u\"]",
               zero_count,
               one_count,
               Q1, Q2, Q3,
               counts[0], counts[1], counts[2], counts[3]);
    }
    else
    {
        printf("Total\t0:1\t%d:%d\t%u:%u:%u\t%d:%d:%d:%d",
               zero_count, one_count,
               Q1, Q2, Q3,
               counts[0], counts[1], counts[2], counts[3]);
    }

    if (depths_only == 1)
    {

        for (size_t num_sample_idx = 0; num_sample_idx < size_num_samples_all; num_sample_idx++)
        {
            uint32_t *sample_depths = NULL;
            struct uint_pair *sample_alt_depths = sample_alt_depths_all[num_sample_idx];
            ret = stix_get_sample_depths(sample_alt_depths, // todo! -> sharding
                                         NULL,
                                         num_samples_all[num_sample_idx],
                                         &sample_depths);
            if (json_out == 0)
            {
                for (i = 0; i < num_samples_all[num_sample_idx]; ++i)
                    printf("\t%d", sample_depths[i]);
            }
            else
            {
                printf(",\"depths_shard%ld\":[", num_sample_idx);
                for (i = 0; i < num_samples_all[num_sample_idx]; ++i)
                {
                    if (i != 0)
                        printf(",");
                    printf("%d", sample_depths[i]);
                }
                printf("]");
            }

            free(sample_depths);
        }
    }

    if (json_out == 1)
        printf("}");
    printf("\n");

    // sqlite3 *db = NULL;
    /*---------------------*/
    if (num_agg_cols > 0)
    {
        // printf(">>> %d",num_agg_cols);
        for (size_t num_sample_idx = 0; num_sample_idx < size_num_samples_all; num_sample_idx++)
        {
            char ***uniq_vals;
            uint32_t num_uniq_vals;
            uint32_t **uniq_groups_ids;
            uint32_t *uniq_groups_sizes;
            sqlite3 *db = NULL;
            num_uniq_vals = ped_get_uniq_col_groups(shard_array[num_sample_idx].stixdb_path,
                                                    &db,
                                                    agg_cols,
                                                    num_agg_cols,
                                                    filter,
                                                    &uniq_vals,
                                                    &uniq_groups_ids,
                                                    &uniq_groups_sizes);

            for (i = 0; i < num_uniq_vals; ++i)
            {
                uint32_t j;
                char *group_name_tmp, *group_name;

                ret = asprintf(&group_name, "%s", uniq_vals[i][0]);
                for (j = 1; j < num_agg_cols; ++j)
                {
                    ret = asprintf(&group_name_tmp,
                                   "%s,%s",
                                   group_name,
                                   uniq_vals[i][j]);
                    free(group_name);
                    group_name = group_name_tmp;
                }
                // fprintf(stderr,">>>>>>>>>>>>>>>>>>>%s",group_name);
                ret = stix_get_summary(sample_alt_depths_all[num_sample_idx],
                                       uniq_groups_ids[i],
                                       uniq_groups_sizes[i],
                                       &zero_count,
                                       &one_count,
                                       &Q1,
                                       &Q2,
                                       &Q3,
                                       &min,
                                       &max,
                                       counts);

                if (json_out)
                {
                    printf(",{"
                           "\"name\":\"%s\", "
                           "\"zero_count\":\"%d\", "
                           "\"one_count\":\"%d\","
                           "\"quantiles\":[\"%u\",\"%u\",\"%u\"], "
                           "\"counts\":[\"%u\",\"%u\",\"%u\",\"%u\"]",
                           group_name,
                           zero_count,
                           one_count,
                           Q1, Q2, Q3,
                           counts[0], counts[1], counts[2], counts[3]);
                }
                else
                {
                    printf("%s\t0:1\t%d:%d\t%u:%u:%u\t%d:%d:%d:%d",
                           group_name,
                           zero_count, one_count,
                           Q1, Q2, Q3,
                           counts[0], counts[1], counts[2], counts[3]);
                }

                if (depths_only == 1)
                {
                    uint32_t *sample_depths = NULL;
                    ret = stix_get_sample_depths(sample_alt_depths_all[num_sample_idx],
                                                 uniq_groups_ids[i],
                                                 uniq_groups_sizes[i],
                                                 &sample_depths);
                    int k;
                    if (json_out == 0)
                    {
                        for (k = 0; k < uniq_groups_sizes[i]; ++k)
                            printf("\t%d", sample_depths[k]);
                    }
                    else
                    {
                        printf(",\"depths\":[");
                        for (k = 0; k < uniq_groups_sizes[i]; ++k)
                        {
                            if (k != 0)
                                printf(",");
                            printf("%d", sample_depths[k]);
                        }
                        printf("]");
                    }
                    free(sample_depths);
                }

                if (json_out == 1)
                    printf("}");
                printf("\n");

                free(group_name);
            }
        }
    }
    /*---------------------*/

    if (summary_only == 1)
    {
        if (json_out == 1)
            printf("]}}");
        return;
    }

    if (depths_only == 0)
    {
        int current_line_no = 0;
        for (size_t num_sample_idx = 0; num_sample_idx < size_num_samples_all; num_sample_idx++)
        {

            if (json_out)
                printf("],\n\"samples_shard%ld\": [\n", num_sample_idx);

            sqlite3 *db = NULL;
            // printf("num_samples:%d\n",num_samples);
            uint32_t num_samples = num_samples_all[num_sample_idx];
            for (i = 0; i < num_samples; ++i)
            {
                uint32_t idx;
                idx = i;
                char **col_vals = NULL, **col_names = NULL;

                // printf("ped_get_cols_info_by_id:%d,num_sample_idx:%ld \n",idx,num_sample_idx);
                // printf("shard_array[num_sample_idx].stixdb_path:%s\n",shard_array[num_sample_idx].stixdb_path);
                num_col_vals = ped_get_cols_info_by_id(shard_array[num_sample_idx].stixdb_path, // todo! -> sharding
                                                       &db,
                                                       NULL,
                                                       0,
                                                       idx,
                                                       &col_vals,
                                                       &col_names);

                uint32_t j;
                if (json_out == 1)
                {
                    if (i > 0)
                        printf(",");

                    printf("{");
                    for (j = 0; j < num_col_vals; ++j)
                    {
                        if (j > 0)
                            printf(",");
                        printf("\"%s\":\"%s\"", col_names[j], col_vals[j]);
                    }
                    printf(",\"Pairend\":\"%u\",\"Split\":\"%u\"}\n",
                           sample_alt_depths_all[num_sample_idx][i].first,
                           sample_alt_depths_all[num_sample_idx][i].second);
                    // fflush(stdout);
                }
                else
                {
                    if (num_sample_idx == 0 && i == 0) // only print header once
                    {
                        for (j = 0; j < num_col_vals; ++j)
                            printf("%s\t", col_names[j]);
                        printf("Pairend\tSplit\n");
                    }

                    // printf("%d\t",current_line_no);
                    // current_line_no ++;

                    for (j = 0; j < num_col_vals; ++j)
                        printf("%s\t", col_vals[j]);

                    printf("%u\t%u\n",
                           sample_alt_depths_all[num_sample_idx][i].first,
                           sample_alt_depths_all[num_sample_idx][i].second);
                    // fflush(stdout);
                }
            }
            sqlite3_close(db);
            fflush(stdout);
        }
    }

    if (json_out == 1)
    {
        printf("]}}");
        fflush(stdout);
    }

    // sqlite3_close(db);
}
//}}}

//{{{uint32_t parse_aggregate_csv(char *aggregate,
uint32_t parse_aggregate_csv(char *aggregate,
                             char ***agg_cols)
{
    uint32_t i;
    uint32_t num_agg_cols = 1;
    for (i = 0; i < strlen(aggregate); ++i)
    {
        if (aggregate[i] == ',')
            num_agg_cols += 1;
    }

    *agg_cols = (char **)malloc(num_agg_cols * sizeof(char *));

    char *p = strtok(aggregate, ",");
    uint32_t col_i = 0;
    while (p != NULL)
    {
        (*agg_cols)[col_i] = p;
        col_i += 1;
        p = strtok(NULL, ",");
    }

    return num_agg_cols;
}
//}}}

//{{{void update_vcf_header(bcf_hdr_t *hdr)
void update_vcf_header(bcf_hdr_t *hdr,
                       uint32_t v_is_set,
                       char *sample_column)
{
    char *STIX_QUANTS_LINE =
        "##INFO=<ID=STIX_QUANTS,Number=3,Type=Integer,"
        "Description=\"STIX quantile values.\">";
    if (bcf_hdr_append(hdr, STIX_QUANTS_LINE) != 0)
        errx(EX_DATAERR, "Error updating header\n");

    char *STIX_QUANT_DEPTHS_LINE =
        "##INFO=<ID=STIX_QUANT_DEPTHS,Number=4,Type=Integer,"
        "Description=\"STIX quantile depths.\">";
    if (bcf_hdr_append(hdr, STIX_QUANT_DEPTHS_LINE) != 0)
        errx(EX_DATAERR, "Error updating header\n");

    char *STIX_ZERO_LINE =
        "##INFO=<ID=STIX_ZERO,Number=1,Type=Integer,"
        "Description=\"STIX samples with zero evidence.\">";
    if (bcf_hdr_append(hdr, STIX_ZERO_LINE) != 0)
        errx(EX_DATAERR, "Error updating header\n");

    char *STIX_ONE_LINE =
        "##INFO=<ID=STIX_ONE,Number=1,Type=Integer,"
        "Description=\"STIX samples with one piece of evidence.\">";
    if (bcf_hdr_append(hdr, STIX_ONE_LINE) != 0)
        errx(EX_DATAERR, "Error updating header\n");

    char *MIN_SUP_READS_LINE =
        "##INFO=<ID=STIX_ONE_MIN_READS,Number=1,Type=Integer,"
        "Description=\"STIX How many supporting reads define a hit,Used for calucation of STIX_ONE and STIX_ZERO.\">";
    if (bcf_hdr_append(hdr, MIN_SUP_READS_LINE) != 0)
        errx(EX_DATAERR, "Error updating header\n");

    if (v_is_set == 1)
    {
        char *STIX_SAMPLE_DEPTH_LINE =
            "##INFO=<ID=STIX_SAMPLE_DEPTH,Number=.,Type=String,"
            "Description=\"STIX sample-level depth information.\">";
        if (bcf_hdr_append(hdr, STIX_SAMPLE_DEPTH_LINE) != 0)
            errx(EX_DATAERR, "Error updating header\n");
    }
}
//}}}

int V_is_set = 0; // init global Verbose mode
int L_is_set = 0; // init global Length of Insertions
int R_is_set = 0;
int T_is_set = 0;
uint32_t length_of_insertion = 0;
float ovpct_threshold = 0.15;
uint32_t min_supporting_reads = 1;

//{{{int main(int argc, char **argv)
int main(int argc, char **argv)
{
    int c;
    int i_is_set = 0,
        s_is_set = 0,
        p_is_set = 0,
        c_is_set = 0,
        d_is_set = 0,
        r_is_set = 0,
        l_is_set = 0,
        f_is_set = 0,
        a_is_set = 0,
        F_is_set = 0,
        j_is_set = 0,
        t_is_set = 0,
        v_is_set = 0,
        P_is_set = 0,
        Q_is_set = 0,
        B_is_set = 0;

    char *index_dir_name = NULL;
    char *ped_file_name = NULL;
    char *ped_db_file_name = NULL;
    char *vcf_file_name = NULL;
    char *r_region = NULL;
    char *l_region = NULL;
    char *aggregate = NULL;
    char *filter = NULL;
    char *sv_type = NULL;
    char *sample_column = NULL;
    char *sharding_file_name = NULL;
    char *table_query_file_name = NULL;

    uint32_t slop = 0;
    uint32_t ins_padding = 50;
    uint32_t col_id = 1;
    uint32_t summary_only = 0;
    uint32_t depths_only = 0;

    while ((c = getopt(argc, argv, "i:P:s:T:B:Q:p:c:d:r:l:f:a:F:jt:v:SDVL:R:")) != -1)
    {
        switch (c)
        {
        case 'i':
            i_is_set = 1;
            index_dir_name = optarg;
            break;
        case 's':
            s_is_set = 1;
            slop = atoi(optarg);
            break;
        case 'B':
            B_is_set = 1;
            sharding_file_name = optarg;
            break;
        case 'Q':
            Q_is_set = 1;
            table_query_file_name = optarg;
            break;
        case 'p':
            p_is_set = 1;
            ped_file_name = optarg;
            break;
        case 'P':
            P_is_set = 1;
            ins_padding = atoi(optarg);
            break;
        case 'T':
            T_is_set = 1;
            min_supporting_reads = atoi(optarg);
            break;
        case 'c':
            c_is_set = 1;
            col_id = atoi(optarg);
            break;
        case 'd':
            d_is_set = 1;
            ped_db_file_name = optarg;
            break;
        case 'r':
            r_is_set = 1;
            r_region = optarg;
            break;
        case 'l':
            l_is_set = 1;
            l_region = optarg;
            break;
        case 'f':
            f_is_set = 1;
            vcf_file_name = optarg;
            break;
        case 'a':
            a_is_set = 1;
            aggregate = optarg;
            break;
        case 'F':
            F_is_set = 1;
            filter = optarg;
            break;
        case 'j':
            j_is_set = 1;
            break;
        case 't':
            t_is_set = 1;
            sv_type = optarg;
            break;
        case 'L':
            L_is_set = 1;
            length_of_insertion = atoi(optarg);
            break;
        case 'R':
            R_is_set = 1;
            ovpct_threshold = atof(optarg);
            break;
        case 'v':
            v_is_set = 1;
            sample_column = optarg;
            break;
        case 'S':
            summary_only = 1;
            break;
        case 'D':
            depths_only = 1;
            break;
        case 'V':
            V_is_set = 1;
            break;
        case 'h':
            return help(EX_OK);
        case '?':
            if ((optopt == 'i') ||
                (optopt == 's') ||
                (optopt == 'B') ||
                (optopt == 'Q') ||
                (optopt == 'T') ||
                (optopt == 'p') ||
                (optopt == 'P') ||
                (optopt == 'd') ||
                (optopt == 'r') ||
                (optopt == 'l') ||
                (optopt == 'f') ||
                (optopt == 'a') ||
                (optopt == 'F') ||
                (optopt == 't') ||
                (optopt == 'L') ||
                (optopt == 'R') ||
                (optopt == 'v'))
                fprintf(stderr,
                        "Option -%c requires an argument.\n",
                        optopt);
            else if (isprint(optopt))
                fprintf(stderr,
                        "Unknown option `-%c'.\n",
                        optopt);
            else
                fprintf(stderr,
                        "Unknown option character `\\x%x'.\n",
                        optopt);

            return help(EX_USAGE);
        default:
            return help(EX_OK);
        }
    }

    char **agg_cols = NULL;
    uint32_t num_agg_cols = 0;
    if (a_is_set == 1)
        num_agg_cols = parse_aggregate_csv(aggregate, &agg_cols);

    /*create index*/
    if ((i_is_set == 1) &&
        (p_is_set == 1) &&
        (d_is_set == 1))
    {

        uint32_t num_rows = ped_create_db(ped_file_name,
                                          ped_db_file_name,
                                          index_dir_name,
                                          col_id);
        fprintf(stderr,
                "PED database %s with %u rows created from %s.\n",
                ped_db_file_name,
                num_rows,
                ped_file_name);
        return EX_OK;
    }

    /*search*/
    if (B_is_set == 0)
    {                          // single db mode
        if ((i_is_set == 1) && // giggle index
            (d_is_set == 1) && // ped db
            (s_is_set == 1))
        { // slop

            uint32_t *sample_ids = NULL;
            uint32_t num_samples = 0;
            if (F_is_set == 1)
                num_samples = ped_get_matching_sample_ids(ped_db_file_name,
                                                          filter,
                                                          &sample_ids);

            struct giggle_index *gi = NULL;

            struct uint_pair *sample_alt_depths = NULL;
            struct stix_breakpoint *left = NULL, *right = NULL;

            /*Single query mode*/
            if ((l_is_set == 1) && // left interval
                (r_is_set == 1) && // right interval
                (t_is_set == 1))
            { // sv type

                left = stix_region_to_breakpoint(l_region);
                right = stix_region_to_breakpoint(r_region);

                enum stix_sv_type query_type = DEL;

                if (strcmp(sv_type, "DUP") == 0)
                    query_type = DUP;
                else if (strcmp(sv_type, "INV") == 0)
                    query_type = INV;
                else if (strcmp(sv_type, "BND") == 0)
                    query_type = BND;
                else if (strcmp(sv_type, "INS") == 0)
                    query_type = INS;

                uint32_t num_sample_alt_depths =
                    stix_run_giggle_query(&gi,
                                          index_dir_name,
                                          query_type,
                                          left,
                                          right,
                                          slop,
                                          ins_padding,
                                          sample_ids,
                                          num_samples,
                                          &sample_alt_depths,
                                          0);
                // DEBUG
                //  printf("num_sample_alt_depths:%d\n",num_sample_alt_depths);
                print_results(gi,
                              ped_db_file_name,
                              left,
                              right,
                              sample_ids,
                              filter,
                              sample_alt_depths,
                              num_sample_alt_depths,
                              agg_cols,
                              num_agg_cols,
                              j_is_set,
                              summary_only,
                              depths_only);

                free(left->chrm);
                free(left);
                free(right->chrm);
                free(right);
            }
            else if (f_is_set == 1) /*VCF input mode*/
            {

                htsFile *fp = hts_open(vcf_file_name, "r");
                if (!fp)
                    err(EX_DATAERR, "Could not read file: %s", vcf_file_name);

                bcf_hdr_t *hdr = bcf_hdr_read(fp);
                if (!hdr)
                    errx(EX_DATAERR, "Header not found: %s\n", vcf_file_name);

                update_vcf_header(hdr, v_is_set, sample_column);

                htsFile *out_f = hts_open("-", "w");
                bcf_hdr_write(out_f, hdr);

                bcf1_t *line = bcf_init1();

                left = (struct stix_breakpoint *)
                    malloc(sizeof(struct stix_breakpoint));
                left->chrm = NULL;

                right = (struct stix_breakpoint *)
                    malloc(sizeof(struct stix_breakpoint));
                right->chrm = NULL;

                enum stix_sv_type type;

                sqlite3 *db = NULL;
                char **cols = NULL;
                uint32_t num_cols = 0;

                if (v_is_set == 1)
                {
                    cols = (char **)malloc(sizeof(char *));
                    cols[0] = sample_column;
                    num_cols = 1;
                }

                while (bcf_read(fp, hdr, line) == 0)
                {
                    if (stix_get_vcf_breakpoints(fp,
                                                 hdr,
                                                 line,
                                                 left,
                                                 right,
                                                 &type) == 0)
                    {

                        uint32_t num_sample_alt_depths =
                            stix_run_giggle_query(&gi,
                                                  index_dir_name,
                                                  type,
                                                  left,
                                                  right,
                                                  slop,
                                                  ins_padding,
                                                  sample_ids,
                                                  num_samples,
                                                  &sample_alt_depths,
                                                  0);

                        uint32_t i;
                        uint32_t num_col_vals;
                        int32_t zero_count, one_count;
                        uint32_t Q1, Q2, Q3, min, max;
                        int32_t counts[4];

                        uint32_t ret = stix_get_summary(sample_alt_depths,
                                                        sample_ids,
                                                        num_sample_alt_depths,
                                                        &zero_count,
                                                        &one_count,
                                                        &Q1,
                                                        &Q2,
                                                        &Q3,
                                                        &min,
                                                        &max,
                                                        counts);

                        if (v_is_set == 1)
                        {
                            char *stix_sample_depth_string = NULL,
                                 *tmp_string = NULL;

                            for (i = 0; i < num_sample_alt_depths; ++i)
                            {
                                if (sample_alt_depths[i].first +
                                        sample_alt_depths[i].second >
                                    0)
                                {

                                    char **col_vals = NULL, **col_names = NULL;
                                    num_col_vals = ped_get_cols_info_by_id(
                                        ped_db_file_name,
                                        &db,
                                        cols,
                                        num_cols,
                                        i,
                                        &col_vals,
                                        &col_names);
                                    if (stix_sample_depth_string != NULL)
                                    {
                                        ret = asprintf(&tmp_string,
                                                       "%s,%s|%u",
                                                       stix_sample_depth_string,
                                                       col_vals[0],
                                                       sample_alt_depths[i].first +
                                                           sample_alt_depths[i].second);
                                    }
                                    else
                                    {
                                        ret = asprintf(&tmp_string,
                                                       "%s|%u",
                                                       col_vals[0],
                                                       sample_alt_depths[i].first +
                                                           sample_alt_depths[i].second);
                                    }
                                    free(stix_sample_depth_string);
                                    stix_sample_depth_string = tmp_string;

                                    free(col_vals[0]);
                                    free(col_vals);
                                    free(col_names[0]);
                                    free(col_names);
                                }
                            }

                            if (stix_sample_depth_string != NULL)
                            {
                                ret = bcf_update_info_string(
                                    hdr,
                                    line,
                                    "STIX_SAMPLE_DEPTH",
                                    stix_sample_depth_string);
                                if (ret != 0)
                                    errx(EX_DATAERR,
                                         "Error adding STIX_SAMPLE_DEPTH to "
                                         "info field.\n");
                            }
                        }

                        ret = bcf_update_info_int32(hdr,
                                                    line,
                                                    "STIX_ZERO",
                                                    &zero_count,
                                                    1);
                        if (ret != 0)
                            errx(EX_DATAERR,
                                 "Error adding STIX_ZERO to info field.\n");

                        ret = bcf_update_info_int32(hdr,
                                                    line,
                                                    "STIX_ONE",
                                                    &one_count,
                                                    1);
                        if (ret != 0)
                            errx(EX_DATAERR,
                                 "Error adding STIX_ONE to info field.\n");

                        uint32_t quants[3] = {Q1, Q2, Q3};
                        ret = bcf_update_info_int32(hdr,
                                                    line,
                                                    "STIX_QUANTS",
                                                    quants,
                                                    3);
                        if (ret != 0)
                            errx(EX_DATAERR,
                                 "Error adding STIX_QUANTS to info field.\n");

                        ret = bcf_update_info_int32(hdr,
                                                    line,
                                                    "STIX_QUANT_DEPTHS",
                                                    counts,
                                                    4);
                        if (ret != 0)
                            errx(EX_DATAERR,
                                 "Error adding STIX_QUANT_DEPTHS to info field.\n");

                        bcf_write(out_f, hdr, line);

                        /*
                        fprintf(stdout,
                                "0:1\t%d:%d\t%u:%u:%u\t%d:%d:%d:%d\n",
                                zero_count, one_count,
                                Q1, Q2, Q3,
                                counts[0], counts[1], counts[2], counts[3]);
                        */

                        // free(sample_alt_depths);
                    }
                }

                free(left->chrm);
                free(left);
                free(right->chrm);
                free(right);

                bcf_destroy(line);
                bcf_hdr_destroy(hdr);
                hts_close(fp);
                hts_close(out_f);
            }

            if (sample_alt_depths != NULL)
                free(sample_alt_depths);

            if (gi != NULL)
            {
                giggle_index_destroy(&gi);
                cache.destroy();
            }

            return EX_OK;
        }
    }
    else
    { // sharding mode

        int sharding_arr_length = 0;
        Shard *sharding_arr = NULL;
        sharding_arr = read_shards_from_file(sharding_file_name, &sharding_arr_length);
        // fprintf(stderr, "Detected sharding file: %s\n%d shards loaded...\n", sharding_file_name, sharding_arr_length);
        // for (int i = 0; i < sharding_arr_length; i++)
        //     fprintf(stderr, "sharding_arr[i].stixdb_path:%s,%d\n",
        //             sharding_arr[i].stixdb_path,
        //             i);

        /*Sharded Single query mode*/
        if ((l_is_set == 1) && // left interval
            (r_is_set == 1) && // right interval
            (t_is_set == 1))
        { // sv type

            struct stix_breakpoint *left = NULL, *right = NULL;
            left = stix_region_to_breakpoint(l_region);
            right = stix_region_to_breakpoint(r_region);
            enum stix_sv_type query_type = DEL;

            if (strcmp(sv_type, "DUP") == 0)
                query_type = DUP;
            else if (strcmp(sv_type, "INV") == 0)
                query_type = INV;
            else if (strcmp(sv_type, "BND") == 0)
                query_type = BND;
            else if (strcmp(sv_type, "INS") == 0)
                query_type = INS;

            /*from xinchang
            setup array for all independent query.
            */
            uint32_t *sample_ids_all[sharding_arr_length];
            uint32_t num_samples_all[sharding_arr_length];
            uint32_t num_sample_alt_depths_all[sharding_arr_length];
            struct uint_pair *sample_alt_depths_all[sharding_arr_length];
            // struct giggle_index *gi_all[sharding_arr_length];

            for (int i = 0; i < sharding_arr_length; i++)
            {
                /*Individual search in one shard*/
                // fprintf(stderr,"i=%d,sharding_arr[i].stixdb_path=%s",i,sharding_arr[i].stixdb_path);
                char *ped_db_file_name_shard = sharding_arr[i].stixdb_path;
                char *index_dir_name_shard = sharding_arr[i].giggle_path;

                uint32_t num_samples = 0;
                uint32_t *sample_ids = NULL;
                // fprintf(stderr,"%d,index_dir_name_shard:%s\n",i,index_dir_name_shard);
                if (F_is_set == 1)
                    num_samples = ped_get_matching_sample_ids(ped_db_file_name_shard,
                                                              filter,
                                                              &sample_ids);
                struct giggle_index *gi = NULL;
                // struct giggle_index *gi = NULL;
                struct uint_pair *sample_alt_depths = NULL;

                uint32_t num_sample_alt_depths =
                    stix_run_giggle_query(&gi,
                                          index_dir_name_shard,
                                          query_type,
                                          left,
                                          right,
                                          slop,
                                          ins_padding,
                                          sample_ids,
                                          num_samples,
                                          &sample_alt_depths,
                                          i);
                /*
                from xinchang
                sample_alt_depths[0].first Pairend count for the first sample
                sample_alt_depths[0].second Split count for the first sample
                */

                sample_ids_all[i] = sample_ids;
                num_sample_alt_depths_all[i] = num_sample_alt_depths;
                sample_alt_depths_all[i] = sample_alt_depths;
                num_samples_all[i] = num_samples;
            }

            uint32_t size_num_sample_alt_depths_all = sizeof(num_sample_alt_depths_all) / sizeof(num_sample_alt_depths_all[0]);

            print_results_shard(
                // gi_all,
                sharding_arr,
                left,
                right,
                sample_ids_all,
                filter,
                sample_alt_depths_all,
                num_sample_alt_depths_all,
                size_num_sample_alt_depths_all,
                agg_cols,
                num_agg_cols,
                j_is_set,
                summary_only,
                depths_only);

            free(left->chrm);
            free(left);
            free(right->chrm);
            free(right);
        }

        else if (Q_is_set == 1)
        {

            int table_query_length = 0;
            TableQuery *table_query_arr =
                read_table_queries_from_file(table_query_file_name, &table_query_length);
            if (V_is_set)
                printf("%s,%s,%s,%s,%s\n",
                       table_query_arr[0].left_str,
                       table_query_arr[0].right_str,
                       table_query_arr[0].len,
                       table_query_arr[0].svtype,
                       table_query_arr[0].ID);

            for (int table_q_idx = 0; table_q_idx < table_query_length; table_q_idx++)
            {
                fprintf(stdout, ">>>%s\t%s\t%s\t%s\t%s\n", table_query_arr[table_q_idx].left_str,
                        table_query_arr[table_q_idx].right_str,
                        table_query_arr[table_q_idx].len,
                        table_query_arr[table_q_idx].svtype,
                        table_query_arr[table_q_idx].ID);
                fflush(stdout);

                struct stix_breakpoint *left = NULL, *right = NULL;
                left = stix_region_to_breakpoint(table_query_arr[table_q_idx].left_str);
                right = stix_region_to_breakpoint(table_query_arr[table_q_idx].right_str);
                enum stix_sv_type query_type = DEL;

                if (strcmp(table_query_arr[table_q_idx].svtype, "DUP") == 0)
                    query_type = DUP;
                else if (strcmp(table_query_arr[table_q_idx].svtype, "INV") == 0)
                    query_type = INV;
                else if (strcmp(table_query_arr[table_q_idx].svtype, "BND") == 0)
                    query_type = BND;
                else if (strcmp(table_query_arr[table_q_idx].svtype, "INS") == 0)
                    query_type = INS;

                /*from xinchang
                setup array for all independent query.
                */
                uint32_t *sample_ids_all[sharding_arr_length];
                uint32_t num_samples_all[sharding_arr_length];
                uint32_t num_sample_alt_depths_all[sharding_arr_length];
                struct uint_pair *sample_alt_depths_all[sharding_arr_length];
                // struct giggle_index *gi_all[sharding_arr_length];

                for (int i = 0; i < sharding_arr_length; i++)
                {
                    /*Individual search in one shard*/
                    // fprintf(stderr,"i=%d,sharding_arr[i].stixdb_path=%s",i,sharding_arr[i].stixdb_path);
                    char *ped_db_file_name_shard = sharding_arr[i].stixdb_path;
                    char *index_dir_name_shard = sharding_arr[i].giggle_path;

                    uint32_t num_samples = 0;
                    uint32_t *sample_ids = NULL;
                    // fprintf(stderr,"%d,index_dir_name_shard:%s\n",i,index_dir_name_shard);
                    if (F_is_set == 1)
                        num_samples = ped_get_matching_sample_ids(ped_db_file_name_shard,
                                                                  filter,
                                                                  &sample_ids);
                    struct giggle_index *gi = NULL;
                    // struct giggle_index *gi = NULL;
                    struct uint_pair *sample_alt_depths = NULL;

                    uint32_t num_sample_alt_depths =
                        stix_run_giggle_query(&gi,
                                              index_dir_name_shard,
                                              query_type,
                                              left,
                                              right,
                                              slop,
                                              ins_padding,
                                              sample_ids,
                                              num_samples,
                                              &sample_alt_depths,
                                              i);
                    sample_ids_all[i] = sample_ids;
                    num_sample_alt_depths_all[i] = num_sample_alt_depths;
                    sample_alt_depths_all[i] = sample_alt_depths;
                    num_samples_all[i] = num_samples;

                    if (gi != NULL)
                    {
                        giggle_index_destroy(&gi);
                        cache.destroy();
                    }
                }

                uint32_t size_num_sample_alt_depths_all = sizeof(num_sample_alt_depths_all) / sizeof(num_sample_alt_depths_all[0]);

                print_results_shard(
                    // gi_all,
                    sharding_arr,
                    left,
                    right,
                    sample_ids_all,
                    filter,
                    sample_alt_depths_all,
                    num_sample_alt_depths_all,
                    size_num_sample_alt_depths_all,
                    agg_cols,
                    num_agg_cols,
                    j_is_set,
                    summary_only,
                    depths_only);



                /*
                Free all corresponding variables
                */
               
                free(left->chrm);
                free(left);
                free(right->chrm);
                free(right);

                for (int i = 0; i < sharding_arr_length; i++)
                {

                    if (sample_alt_depths_all[i] != NULL)
                    {
                        free(sample_alt_depths_all[i]);
                    }
                    if (sample_ids_all[i] != NULL)
                    {
                        free(sample_ids_all[i]);
                    }
                }

                if (table_query_arr != NULL)
                {
                    for (int i=0; i < table_query_length; i++)
                    {
                        free(table_query_arr[i].left_str);
                        free(table_query_arr[i].right_str);
                        free(table_query_arr[i].len);
                        free(table_query_arr[i].svtype);
                        free(table_query_arr[i].ID);
                    }
                }
            }
        }
        else if (f_is_set == 1) /*VCF input mode*/
        {

            htsFile *fp = hts_open(vcf_file_name, "r");
            if (!fp)
                err(EX_DATAERR, "Could not read file: %s", vcf_file_name);

            bcf_hdr_t *hdr = bcf_hdr_read(fp);
            if (!hdr)
                errx(EX_DATAERR, "Header not found: %s\n", vcf_file_name);

            update_vcf_header(hdr, v_is_set, sample_column);

            htsFile *out_f = hts_open("-", "w");
            bcf_hdr_write(out_f, hdr);

            bcf1_t *line = bcf_init1();

            enum stix_sv_type query_type;

            /*----------*/
            // struct giggle_index *gi_all[sharding_arr_length];
            // for (int i = 0; i < sharding_arr_length; ++i)
            // {
            //     gi_all[i] = NULL;
            // }
            int Iter_number = 0;
            const int max_iter = 4;
            while (bcf_read(fp, hdr, line) == 0)
            {
                Iter_number++;
                struct stix_breakpoint *left = NULL, *right = NULL;
                left = (struct stix_breakpoint *)
                    malloc(sizeof(struct stix_breakpoint));
                left->chrm = NULL;

                right = (struct stix_breakpoint *)
                    malloc(sizeof(struct stix_breakpoint));
                right->chrm = NULL;

                if (stix_get_vcf_breakpoints(fp,
                                             hdr,
                                             line,
                                             left,
                                             right,
                                             &query_type) == 0)
                {

                    uint32_t *sample_ids_all[sharding_arr_length];
                    uint32_t num_samples_all[sharding_arr_length];
                    uint32_t num_sample_alt_depths_all[sharding_arr_length];
                    struct uint_pair *sample_alt_depths_all[sharding_arr_length];

                    for (int i = 0; i < sharding_arr_length; i++)
                    {

                        char *ped_db_file_name_shard = sharding_arr[i].stixdb_path;
                        char *index_dir_name_shard = sharding_arr[i].giggle_path;

                        uint32_t num_samples = 0;
                        uint32_t *sample_ids = NULL;

                        struct giggle_index *gi = NULL;
                        // fprintf(stderr,"@@@@@@@@@@@@@%p\n",gi_all[i]);
                        // if (gi_all[i] != NULL)
                        // {
                        //     gi = gi_all[i];
                        // }
                        // struct giggle_index *gi =NULL;
                        // struct giggle_index *gi = NULL;
                        struct uint_pair *sample_alt_depths = NULL;

                        // if (F_is_set == 1)
                        //     num_samples = ped_get_matching_sample_ids(ped_db_file_name_shard,
                        //                                               filter,
                        //                                               &sample_ids);
                        // fprintf(stderr,">>>>>>>>>>%p\n",gi);
                        uint32_t num_sample_alt_depths =
                            stix_run_giggle_query(&gi,
                                                  index_dir_name_shard,
                                                  query_type,
                                                  left,
                                                  right,
                                                  slop,
                                                  ins_padding,
                                                  sample_ids,
                                                  num_samples,
                                                  &sample_alt_depths,
                                                  i);
                        sample_ids_all[i] = sample_ids;
                        num_sample_alt_depths_all[i] = num_sample_alt_depths;
                        sample_alt_depths_all[i] = sample_alt_depths;
                        num_samples_all[i] = num_samples;
                        // if (gi_all[i] != NULL)
                        //     gi_all[i] = gi;

                        // if (gi_all[i] == NULL)
                        // {
                        //     // fprintf(stderr,"<<<<<<<<<%p,%p\n",gi,gi_all[i]);
                        //     gi_all[i] = gi;
                        // }

                        if (gi != NULL)
                        {
                            giggle_index_destroy(&gi);
                            cache.destroy();
                        }
                    }

                    /*from xinchang
                    不知道catch怎么来的，无法批量删除。所以每次调用都要销毁giggle index，damn!
                    */
                    // for (int i = 0; i < sharding_arr_length; ++i)
                    // {
                    //     giggle_index_destroy(&gi_all[i]);
                    //     cache.destroy();
                    // }

                    // uint32_t i;
                    // uint32_t num_col_vals;
                    int32_t zero_count, one_count;
                    uint32_t Q1, Q2, Q3, min, max;
                    int32_t counts[4];
                    // uint32_t size_num_samples_all = sizeof(num_samples_all) / sizeof(num_samples_all[0]);
                    uint32_t size_num_sample_alt_depths_all = sizeof(num_sample_alt_depths_all) / sizeof(num_sample_alt_depths_all[0]);
                    // if(V_is_set)
                    //     fprintf(stderr,">>>>>>>>>>num_samples_all[0]:%d\tnum_samples_all[1]:%d,",num_sample_alt_depths_all[0],num_sample_alt_depths_all[1]);
                    uint32_t ret = stix_get_summary_shard(sample_alt_depths_all,
                                                          sample_ids_all,
                                                          num_sample_alt_depths_all,
                                                          size_num_sample_alt_depths_all,
                                                          &zero_count,
                                                          &one_count,
                                                          &Q1,
                                                          &Q2,
                                                          &Q3,
                                                          &min,
                                                          &max,
                                                          counts);

                    /*----------------------------------------

                    if (v_is_set == 1)
                    {

                        // fprintf(stderr,">>>>>>>>>>>");
                        uint32_t num_cols = 0;
                        char **cols = NULL;
                        cols = (char **)malloc(sizeof(char *));
                        cols[0] = sample_column;
                        num_cols = 1;
                        char *stix_sample_depth_string = NULL;

                        for (int shard_idx = 0; shard_idx < sharding_arr_length; shard_idx++)
                        {


                            char *tmp_string = stix_sample_depth_string;
                            uint32_t num_sample_alt_depths = num_sample_alt_depths_all[shard_idx];
                            struct uint_pair *sample_alt_depths = sample_alt_depths_all[shard_idx];
                            // fprintf(stderr,"num_sample:%d\n",num_sample_alt_depths);
                            for (uint32_t i = 0; i < num_sample_alt_depths; ++i)
                            {

                                fprintf(stderr,"sample_alt_depths:%d,%d\n",sample_alt_depths[i].first,sample_alt_depths[i].second);
                                if (sample_alt_depths[i].first +
                                        sample_alt_depths[i].second >
                                    0)
                                {
                                    fprintf(stderr,"22ddcols:%s\n",*cols);
                                    fprintf(stderr,"sharding_arr[i].stixdb_path:%s\n",sharding_arr[shard_idx].stixdb_path);
                                    sqlite3 *db = NULL;
                                    char **col_vals = NULL, **col_names = NULL;
                                    uint32_t num_col_vals = ped_get_cols_info_by_id(
                                        sharding_arr[i].stixdb_path,
                                        &db,
                                        cols,
                                        num_cols,
                                        i,
                                        &col_vals,
                                        &col_names);
                                    fprintf(stderr,"stix_sample_depth_string:%s,col_vals:%s,col_names:%d\n",stix_sample_depth_string,col_vals[0],*col_names[0]);
                                    if (stix_sample_depth_string != NULL)
                                    {
                                        ret = asprintf(&tmp_string,
                                                       "%s,%s|%u",
                                                       stix_sample_depth_string,
                                                       col_vals[0],
                                                       sample_alt_depths[i].first +
                                                           sample_alt_depths[i].second);
                                    }
                                    else
                                    {
                                        ret = asprintf(&tmp_string,
                                                       "%s|%u",
                                                       col_vals[0],
                                                       sample_alt_depths[i].first +
                                                           sample_alt_depths[i].second);

                                    }

                                    stix_sample_depth_string = tmp_string;

                                }
                            }



                            if (stix_sample_depth_string != NULL)
                            {
                                ret = bcf_update_info_string(
                                    hdr,
                                    line,
                                    "STIX_SAMPLE_DEPTH",
                                    stix_sample_depth_string);
                                if (ret != 0)
                                    errx(EX_DATAERR,
                                         "Error adding STIX_SAMPLE_DEPTH to "
                                         "info field.\n");
                            }
                            // free(col_vals[0]);
                            // free(col_vals);
                            // free(col_names[0]);
                            // free(col_names);
                            // sqlite3_close(db);
                        }

                        free(stix_sample_depth_string);

                    }
                    ----------------------------------------
                    */

                    ret = bcf_update_info_int32(hdr,
                                                line,
                                                "STIX_ZERO",
                                                &zero_count,
                                                1);
                    if (ret != 0)
                        errx(EX_DATAERR,
                             "Error adding STIX_ZERO to info field.\n");

                    ret = bcf_update_info_int32(hdr,
                                                line,
                                                "STIX_ONE",
                                                &one_count,
                                                1);
                    if (ret != 0)
                        errx(EX_DATAERR,
                             "Error adding STIX_ONE to info field.\n");

                    ret = bcf_update_info_int32(hdr,
                                                line,
                                                "STIX_ONE_MIN_READS",
                                                &min_supporting_reads,
                                                1);
                    if (ret != 0)
                        errx(EX_DATAERR,
                             "Error adding STIX_ONE to info field.\n");

                    uint32_t quants[3] = {Q1, Q2, Q3};
                    ret = bcf_update_info_int32(hdr,
                                                line,
                                                "STIX_QUANTS",
                                                quants,
                                                3);
                    if (ret != 0)
                        errx(EX_DATAERR,
                             "Error adding STIX_QUANTS to info field.\n");

                    ret = bcf_update_info_int32(hdr,
                                                line,
                                                "STIX_QUANT_DEPTHS",
                                                counts,
                                                4);
                    if (ret != 0)
                        errx(EX_DATAERR,
                             "Error adding STIX_QUANT_DEPTHS to info field.\n");

                    bcf_write(out_f, hdr, line);

                    /*
                    fprintf(stdout,
                            "0:1\t%d:%d\t%u:%u:%u\t%d:%d:%d:%d\n",
                            zero_count, one_count,
                            Q1, Q2, Q3,
                            counts[0], counts[1], counts[2], counts[3]);
                    */
                    for (int i = 0; i < sharding_arr_length; i++)
                    {
                        if (sample_ids_all[i] != NULL)
                        {

                            free(sample_ids_all[i]); // 假设这里存储的是通过malloc分配的uint32_t*
                        }
                        if (sample_alt_depths_all[i] != NULL)
                        {
                            free(sample_alt_depths_all[i]); // 假设这里存储的是通过malloc分配的struct uint_pair*
                        }
                    }

                    //    if (Iter_number % max_iter == 0)
                    //    {
                    //        for (int i = 0; i < sharding_arr_length; ++i)
                    //        {
                    //            giggle_index_destroy(&gi_all[i]);
                    //            cache.destroy();
                    //            gi_all[i] =NULL;
                    //        }
                    //    }
                    //    else
                    //    {
                    //        fprintf(stderr, "new index..\n");
                    //    }
                }

                free(left->chrm);
                free(left);
                free(right->chrm);
                free(right);
            }

            bcf_destroy(line);
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            hts_close(out_f);

            /*----------*/
        }

        /*
        from xinchang
        destory gi at the end of all query in the shard mode
        */
        // todo!
        return EX_OK;
    }

    return help(EX_OK);
}
//}}}
