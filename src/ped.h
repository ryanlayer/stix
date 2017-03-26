#ifndef __PED_H__
#define __PED_H__

#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <sqlite3.h>
#include <sys/stat.h>

struct uint32_t_str_pair {
    uint32_t uint;
    char *str;
};

int uint32_t_str_pair_cmp(const void *a, const void *b);

struct uint32_t_ll_node {
    uint32_t v;
    struct uint32_t_ll_node *next;
};

struct uint32_t_ll {
    struct uint32_t_ll_node *head, *tail;
    uint32_t len;
};

struct str_ll_node {
    char **str;
    uint32_t num_strs;
    struct str_ll_node *next;
};

struct str_ll {
    struct str_ll_node *head, *tail;
    uint32_t len;
};

struct str_col_ll_node {
    char **str;
    char **col;
    uint32_t num_strs;
    struct str_col_ll_node *next;
};

struct str_col_ll {
    struct str_col_ll_node *head, *tail;
    uint32_t len;
};


static int uint32_t_ll_callback(void *ll_p,
                                int argc,
                                char **argv,
                                char **col_name);

static int str_ll_callback(void *ll_p,
                           int argc,
                           char **argv,
                           char **col_name);

static int str_col_ll_callback(void *ll_p,
                               int argc,
                               char **argv,
                               char **col_name);

int check_field_name(char *field_name);

int is_int(char *s, int *v);

int ped_get_column_names_types(char *ped_file_name,
                               char ***ped_field_names,
                               int **ped_field_is_int);

uint32_t ped_create_db(char *ped_file_name,
                       char *ped_db_file_name,
                       char *giggle_index_dir,
                       uint32_t file_name_col);

uint32_t ped_get_matching_sample_ids(char *ped_file_name_db,
                                     char *select_query,
                                     uint32_t **sample_ids);

uint32_t ped_get_uniq_col_groups(char *ped_file_name_db,
                                 sqlite3 **db,
                                 char **cols,
                                 uint32_t num_cols,
                                 char *select_query,
                                 char ****uniq_col_values,
                                 uint32_t ***uniq_groups_ids,
                                 uint32_t **uniq_groups_sizes);

uint32_t ped_union_groups(uint32_t num_groups,
                          uint32_t **uniq_groups_ids,
                          uint32_t *uniq_groups_sizes,
                          uint32_t **union_group_ids);

uint32_t ped_get_cols_info_by_id(char *ped_file_name_db,
                                 sqlite3 **db,
                                 char **cols,
                                 uint32_t num_cols,
                                 uint32_t sample_id,
                                 char ***col_vals,
                                 char ***col_names);
#endif
