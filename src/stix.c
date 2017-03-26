#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <sysexits.h>
#include <err.h>
#include <sqlite3.h>

#include "ped.h"
#include "search.h"

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
                   uint32_t json_out);

//{{{int help(int exit_code)
int help(int exit_code)
{
    fprintf(stderr,
            "usage:   stix <options>\n"
            "         options:\n"
            "             -i  index directory\n"
            "             -s  slop\n"
            "             -p  PED file\n"
            "             -c  Alt file column (default 1)\n"
            "             -d  PED database file\n"
            "             -r  right SV region\n"
            "             -l  left SV region\n"
            "             -f  VCF file\n"
            "             -a  List of columns to aggregate over\n"
            "             -F  Filter samples by PED field query\n"
            "             -j  JSON output\n");
    return exit_code;
}
//}}}

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
        j_is_set = 0;

    char *index_dir_name = NULL;
    char *ped_file_name= NULL;
    char *ped_db_file_name= NULL;
    char *vcf_file_name= NULL;
    char *r_region = NULL;
    char *l_region = NULL;
    char *aggregate = NULL;
    char *filter = NULL;
    uint32_t slop = 0;
    uint32_t col_id = 1;

    while((c = getopt (argc, argv, "i:s:p:c:d:r:l:f:a:F:j")) != -1) {
        switch (c) {
            case 'i':
                i_is_set = 1;
                index_dir_name = optarg;
                break;
            case 's':
                s_is_set = 1;
                slop = atoi(optarg); 
                break;
            case 'p':
                p_is_set = 1;
                ped_file_name = optarg;
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
            case 'h':
                return help(EX_OK);
            case '?':
                if ( (optopt == 'i') ||
                     (optopt == 's') ||
                     (optopt == 'p') ||
                     (optopt == 'd') ||
                     (optopt == 'r') ||
                     (optopt == 'l') ||
                     (optopt == 'f') ||
                     (optopt == 'a') ||
                     (optopt == 'F') )
                    fprintf (stderr,
                             "Option -%c requires an argument.\n",
                             optopt);
                else if (isprint (optopt))
                    fprintf (stderr,
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

    if ((i_is_set == 1) && 
        (p_is_set == 1) && 
        (d_is_set == 1)) {

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

    if ((i_is_set == 1) &&  // giggle index
        (d_is_set == 1) &&  // ped db
        (s_is_set == 1) &&  // slop
        (l_is_set == 1) &&  // left interval
        (r_is_set == 1) ) { // right interval

        uint32_t *sample_ids = NULL;
        uint32_t num_samples = 0;
        if (F_is_set == 1) 
            num_samples = ped_get_matching_sample_ids(ped_db_file_name,
                                                      filter,
                                                      &sample_ids);

        struct giggle_index *gi = NULL;
        struct uint_pair *sample_alt_depths = NULL;

        struct stix_breakpoint *left = stix_region_to_breakpoint(l_region);
        struct stix_breakpoint *right = stix_region_to_breakpoint(r_region);

        uint32_t num_sample_alt_depths = 
                stix_run_giggle_query(&gi,
                                      index_dir_name,
                                      DEL,
                                      left,
                                      right,
                                      slop,
                                      sample_ids,
                                      num_samples,
                                      &sample_alt_depths);


        char **agg_cols = NULL;
        uint32_t num_agg_cols = 0;
        if (a_is_set == 1) 
            num_agg_cols = parse_aggregate_csv(aggregate, &agg_cols);

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
                      j_is_set);
      
        free(left->chrm);
        free(left);
        free(right->chrm);
        free(right);
        free(sample_alt_depths);
        giggle_index_destroy(&gi);
        cache.destroy();
        return EX_OK;
    }

    return help(EX_OK);
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
                   uint32_t json_out)
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

    if (json_out == 1) {
        printf("{ \"results\": {\n");
    }


    if (json_out == 1) {
        printf("\"summary\": [\n");

        printf("{"
               "\"name\":\"Total\", "
               "\"zero_count\":\"%d\", "
               "\"one_count\":\"%d\","
               "\"quantiles\":[\"%u\",\"%u\",\"%u\"], "
               "\"counts\":[\"%u\",\"%u\",\"%u\",\"%u\"]"
               "}\n",
               zero_count,
               one_count,
               Q1, Q2, Q3,
               counts[0], counts[1], counts[2], counts[3]);
    } else {
        printf("Total\t0:1\t%d:%d\t%u:%u:%u\t%d:%d:%d:%d\n",
               zero_count, one_count, 
               Q1, Q2, Q3,
               counts[0], counts[1], counts[2], counts[3]);
    }

    sqlite3 *db = NULL;
         
    if (num_agg_cols > 0) {
        char ***uniq_vals;
        uint32_t num_uniq_vals;
        uint32_t **uniq_groups_ids;
        uint32_t *uniq_groups_sizes;
        fprintf(stderr, "filter:%p\n", filter);
        num_uniq_vals = ped_get_uniq_col_groups(ped_db_file_name,
                                                &db,
                                                agg_cols,
                                                num_agg_cols,
                                                filter,
                                                &uniq_vals,
                                                &uniq_groups_ids,
                                                &uniq_groups_sizes);

        for (i = 0; i < num_uniq_vals; ++i) {
            uint32_t j;
            char *group_name_tmp, *group_name;
            ret = asprintf(&group_name, "%s", uniq_vals[i][0]);
            for (j = 1; j < num_agg_cols; ++j) {
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

            if (json_out) {
                printf(",{"
                       "\"name\":\"%s\", "
                       "\"zero_count\":\"%d\", "
                       "\"one_count\":\"%d\","
                       "\"quantiles\":[\"%u\",\"%u\",\"%u\"], "
                       "\"counts\":[\"%u\",\"%u\",\"%u\",\"%u\"]"
                       "}\n",
                       group_name,
                       zero_count,
                       one_count,
                       Q1, Q2, Q3,
                       counts[0], counts[1], counts[2], counts[3]);
            } else  {
                printf("%s\t0:1\t%d:%d\t%u:%u:%u\t%d:%d:%d:%d\n",
                       group_name,
                       zero_count, one_count, 
                       Q1, Q2, Q3,
                       counts[0], counts[1], counts[2], counts[3]);
            }

            free(group_name);
        }
    }

    if (json_out)
        printf("],\n\"samples\": [\n");

         

    char **col_vals = NULL, **col_names = NULL;
    for ( i = 0; i < num_samples; ++i ) {
        uint32_t idx;
        if (sample_ids != NULL)
            idx = sample_ids[i];
        else
            idx = i;

        num_col_vals = ped_get_cols_info_by_id(ped_db_file_name,
                                               &db,
                                               NULL,
                                               0,
                                               idx,
                                               &col_vals,
                                               &col_names);


        uint32_t j;
        if (json_out == 1) {
            if (i > 0)
                printf(",");

            printf("{");
            for ( j = 0; j < num_col_vals; ++j ) {
                if (j > 0)
                    printf(",");
                printf("\"%s\":\"%s\"", col_names[j], col_vals[j]);
            }
            printf(",\"Pairend\":\"%u\",\"Split\":\"%u\"}\n",
                    sample_alt_depths[i].first,
                    sample_alt_depths[i].second);

        } else {
            if (i == 0) {
                for ( j = 0; j < num_col_vals; ++j)
                    printf("%s\t", col_names[j]);
                printf("Pairend\tSplit\n");
            }

            for ( j = 0; j < num_col_vals; ++j )
                printf("%s\t", col_vals[j]);

            printf("%u\t%u\n",
                    sample_alt_depths[i].first,
                    sample_alt_depths[i].second);
        }
    }
    if (json_out == 1) 
        printf("]}}");

    sqlite3_close(db);
}
//}}}

//{{{uint32_t parse_aggregate_csv(char *aggregate,
uint32_t parse_aggregate_csv(char *aggregate,
                             char ***agg_cols)
{
    uint32_t i;
    uint32_t num_agg_cols = 1;
    for (i = 0; i < strlen(aggregate); ++i) {
        if (aggregate[i] == ',')
            num_agg_cols += 1;
    }

    *agg_cols = (char **)malloc(num_agg_cols*sizeof(char*));

    char *p = strtok(aggregate, ",");
    uint32_t col_i = 0;
    while (p != NULL) {
        (*agg_cols)[col_i] = p;
        col_i+=1;
        p = strtok(NULL, ",");
    }

    return num_agg_cols;
}
//}}}
