#define _GNU_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <sysexits.h>
#include <err.h>
#include <sqlite3.h>

#include <htslib/hfile.h>
#include <htslib/vcf.h>

#include "ped.h"
#include "search.h"

char *stix_sv_type_strings[5] = { "DEL", "DUP", "INS", "INV", "BND" };

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
            "usage:   stix <options>\n"
            "         options:\n"
            "             -i  index directory\n"
            "             -s  slop\n"
            "             -P  padding base piars for query insertion(default 5)\n"
            "             -p  PED file\n"
            "             -c  Alt file column (default 1)\n"
            "             -d  PED database file\n"
            "             -r  right SV region\n"
            "             -l  left SV region\n"
            "             -f  VCF file\n"
            "             -a  List of columns to aggregate over\n"
            "             -F  Filter samples by PED field query\n"
            "             -j  JSON output\n"
            "             -t  SV type \n"
            "             -v  Add sample depth to VCF file\n"
            "             -S  Give only summary\n"
            "             -D  Give sample depth array\n");
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

    if (json_out == 1) {
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
    } else {
        printf("Total\t0:1\t%d:%d\t%u:%u:%u\t%d:%d:%d:%d",
               zero_count, one_count, 
               Q1, Q2, Q3,
               counts[0], counts[1], counts[2], counts[3]);
    }

    if (depths_only == 1) {
        uint32_t *sample_depths = NULL;
        ret = stix_get_sample_depths(sample_alt_depths,
                                     NULL,
                                     num_samples,
                                     &sample_depths);
        if (json_out == 0) {
            for ( i = 0; i < num_samples; ++i) 
                printf("\t%d", sample_depths[i]);
        } else {
            printf(",\"depths\":[");
            for ( i = 0; i < num_samples; ++i) {
                if (i != 0 )
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
         
    if (num_agg_cols > 0) {
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
                       "\"counts\":[\"%u\",\"%u\",\"%u\",\"%u\"]",
                       group_name,
                       zero_count,
                       one_count,
                       Q1, Q2, Q3,
                       counts[0], counts[1], counts[2], counts[3]);
            } else  {
                printf("%s\t0:1\t%d:%d\t%u:%u:%u\t%d:%d:%d:%d",
                       group_name,
                       zero_count, one_count, 
                       Q1, Q2, Q3,
                       counts[0], counts[1], counts[2], counts[3]);
            }

            if (depths_only == 1) {
                uint32_t *sample_depths = NULL;
                ret = stix_get_sample_depths(sample_alt_depths,
                                             uniq_groups_ids[i],
                                             uniq_groups_sizes[i],
                                             &sample_depths);
                int k;
                if (json_out == 0) {
                    for ( k = 0; k < uniq_groups_sizes[i]; ++k)
                        printf("\t%d", sample_depths[k]);
                } else {
                    printf(",\"depths\":[");
                    for ( k = 0; k < uniq_groups_sizes[i]; ++k) {
                        if (k != 0 )
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

    if (summary_only == 1) {
        if (json_out == 1) 
            printf("]}}");
        return;
    }

   
    if (depths_only == 0) {
        if (json_out)
            printf("],\n\"samples\": [\n");

        char **col_vals = NULL, **col_names = NULL;
        // printf("num_samples:%d\n",num_samples);
        for ( i = 0; i < num_samples; ++i ) {
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

    if (v_is_set == 1) {
        char *STIX_SAMPLE_DEPTH_LINE = 
            "##INFO=<ID=STIX_SAMPLE_DEPTH,Number=.,Type=String,"
            "Description=\"STIX sample-level depth information.\">";
        if (bcf_hdr_append(hdr, STIX_SAMPLE_DEPTH_LINE) != 0)
            errx(EX_DATAERR, "Error updating header\n");
    }
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
        j_is_set = 0,
        t_is_set = 0,
        v_is_set = 0,
        P_is_set = 0;

    char *index_dir_name = NULL;
    char *ped_file_name= NULL;
    char *ped_db_file_name= NULL;
    char *vcf_file_name= NULL;
    char *r_region = NULL;
    char *l_region = NULL;
    char *aggregate = NULL;
    char *filter = NULL;
    char *sv_type = NULL;
    char *sample_column = NULL;
    uint32_t slop = 0;
    uint32_t ins_padding = 5;
    uint32_t col_id = 1;
    uint32_t summary_only = 0;
    uint32_t depths_only = 0;

    while((c = getopt (argc, argv, "i:P:s:p:c:d:r:l:f:a:F:jt:v:SD")) != -1) {
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
            case 'P':
                P_is_set = 1;
                ins_padding = atoi(optarg);
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
            case 'h':
                return help(EX_OK);
            case '?':
                if ( (optopt == 'i') ||
                     (optopt == 's') ||
                     (optopt == 'p') ||
                     (optopt == 'P') ||
                     (optopt == 'd') ||
                     (optopt == 'r') ||
                     (optopt == 'l') ||
                     (optopt == 'f') ||
                     (optopt == 'a') ||
                     (optopt == 'F') ||
                     (optopt == 't') ||
                     (optopt == 'v') )
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

    char **agg_cols = NULL;
    uint32_t num_agg_cols = 0;
    if (a_is_set == 1) 
        num_agg_cols = parse_aggregate_csv(aggregate, &agg_cols);

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
        (s_is_set == 1) ) { // slop

        uint32_t *sample_ids = NULL;
        uint32_t num_samples = 0;
        if (F_is_set == 1) 
            num_samples = ped_get_matching_sample_ids(ped_db_file_name,
                                                      filter,
                                                      &sample_ids);

        struct giggle_index *gi = NULL;

        struct uint_pair *sample_alt_depths = NULL;
        struct stix_breakpoint *left = NULL, *right = NULL;

        if ((l_is_set == 1) &&  // left interval
            (r_is_set == 1) &&  // right interval
            (t_is_set == 1) ) { // sv type

            left = stix_region_to_breakpoint(l_region);
            right = stix_region_to_breakpoint(r_region);

            enum stix_sv_type query_type = DEL;

            if (strcmp(sv_type,"DUP") == 0)
                query_type = DUP;
            else if (strcmp(sv_type,"INV") == 0)
                query_type = INV;
            else if (strcmp(sv_type,"BND") == 0)
                query_type = BND;
            else if (strcmp(sv_type,"INS") == 0)
                query_type = INS;
        
            //DEBUG
            // 这里的num_sample_alt_depths是3，所以导致print_results出错
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
                                          &sample_alt_depths);
            //DEBUG
            // printf("num_sample_alt_depths:%d\n",num_sample_alt_depths);
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
        } else if (f_is_set == 1) {

            htsFile *fp = hts_open(vcf_file_name, "r");
            if ( !fp )
                err(EX_DATAERR, "Could not read file: %s", vcf_file_name);

            bcf_hdr_t *hdr = bcf_hdr_read(fp);
            if (!hdr)
                errx(EX_DATAERR, "Header not found: %s\n", vcf_file_name);

            update_vcf_header(hdr, v_is_set, sample_column);

            htsFile *out_f = hts_open("-","w");
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

            if (v_is_set == 1) {
                cols = (char **) malloc(sizeof(char *));
                cols[0] = sample_column;
                num_cols = 1;
            }

            while (bcf_read(fp, hdr, line) == 0) {
                if (stix_get_vcf_breakpoints(fp,
                                             hdr,
                                             line,
                                             left,
                                             right,
                                             &type) == 0) {

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
                                                  &sample_alt_depths);

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
                
                    if (v_is_set == 1) {
                        char *stix_sample_depth_string = NULL,
                             *tmp_string = NULL;

                        for (i = 0; i < num_sample_alt_depths; ++i) {
                            if (sample_alt_depths[i].first + 
                                sample_alt_depths[i].second > 0) {
                                
                                char **col_vals = NULL, **col_names = NULL;
                                num_col_vals = ped_get_cols_info_by_id(
                                        ped_db_file_name,
                                        &db,
                                        cols,
                                        num_cols,
                                        i,
                                        &col_vals,
                                        &col_names);
                                if (stix_sample_depth_string != NULL) {
                                    ret = asprintf(&tmp_string,
                                                   "%s,%s|%u",
                                                   stix_sample_depth_string,
                                                   col_vals[0],
                                                   sample_alt_depths[i].first + 
                                                   sample_alt_depths[i].second);
                                } else {
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

                        if (stix_sample_depth_string != NULL) {
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

                    uint32_t quants[3] = {Q1,Q2,Q3};
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

                    //free(sample_alt_depths);
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
        if (gi != NULL) {
            giggle_index_destroy(&gi);
            cache.destroy();
        }

        return EX_OK;
    }

    return help(EX_OK);
}
//}}}
