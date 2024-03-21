#define _GNU_SOURCE
#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <sqlite3.h>
#include <sys/stat.h>
#include <htslib/vcf.h>
#include <htslib/faidx.h>
#include <err.h>
#include <sysexits.h>
#include <inttypes.h>
#include <errno.h>
#include <errno.h>
#include <libgen.h>

#include <src/giggle_index.h>

#include "ped.h"

//{{{ static int uint32_t_ll_callback(void *ll_p,
static int uint32_t_ll_callback(void *ll_p,
                    int argc,
                    char **argv,
                    char **col_name)
{
    struct uint32_t_ll *ll = (struct uint32_t_ll *)ll_p;

    struct uint32_t_ll_node *new_node = (struct uint32_t_ll_node *)
        malloc(sizeof(struct uint32_t_ll_node));
    if (!new_node )
        err(EX_OSERR, "malloc error");
    new_node->v = atoi(argv[0]);
    new_node->next = NULL;

    if (ll->head == NULL)
        ll->head = new_node;
    else
        ll->tail->next = new_node;

    ll->tail = new_node;

    ll->len = ll->len + 1;
    return 0;
}
//}}}

//{{{ static int str_ll_callback(void *ll_p,
static int str_ll_callback(void *ll_p,
                           int argc,
                           char **argv,
                           char **col_name)
{
    struct str_ll *ll = (struct str_ll *)ll_p;

    struct str_ll_node *new_node = (struct str_ll_node *)
        malloc(sizeof(struct str_ll_node));
    if (!new_node )
        err(EX_OSERR, "ERROR str_ll_callback(): malloc error");

    new_node->num_strs = argc;
    new_node->str = (char **)malloc(argc * sizeof(char *));

    uint32_t i;
    for (i = 0; i < new_node->num_strs; ++i) 
        new_node->str[i] = strdup(argv[i]);

    if (ll->head == NULL)
        ll->head = new_node;
    else
        ll->tail->next = new_node;

    ll->tail = new_node;

    ll->len = ll->len + 1;

    return 0;
}
//}}}

//{{{ static int str_col_ll_callback(void *ll_p,
static int str_col_ll_callback(void *ll_p,
                               int argc,
                               char **argv,
                               char **col_name)
{
    struct str_col_ll *ll = (struct str_col_ll *)ll_p;

    struct str_col_ll_node *new_node = (struct str_col_ll_node *)
        malloc(sizeof(struct str_col_ll_node));
    if (!new_node )
        err(EX_OSERR, "ERROR str_col_ll_callback(): malloc error");

    new_node->num_strs = argc;
    new_node->str = (char **)malloc(argc * sizeof(char *));
    new_node->col = (char **)malloc(argc * sizeof(char *));

    uint32_t i;
    for (i = 0; i < new_node->num_strs; ++i)  {
        new_node->str[i] = strdup(argv[i]);
        new_node->col[i] = strdup(col_name[i]);
    }

    if (ll->head == NULL)
        ll->head = new_node;
    else
        ll->tail->next = new_node;

    ll->tail = new_node;

    ll->len = ll->len + 1;

    return 0;
}
//}}}

//{{{ int check_field_name(char *field_name)
int check_field_name(char *field_name)
{
    // The first character cannot be a numer

    if ((field_name[0] >= '0') && (field_name[0] <= '9'))
        return 0;

    int i;

    for (i = 0; i < strlen(field_name); ++i) {
        if ( (field_name[i] < '0') ||
            ((field_name[i] >= ':') && (field_name[i] <= '@')) ||
            ((field_name[i] >= '[') && (field_name[i] <= '`') &&
                (field_name[i] != '_')) ||
             (field_name[i] > 'z') )
            return i;
    }

    return -1;
}
//}}}

//{{{ int is_int(char *s, int *v)
//base on http://rus.har.mn/blog/2014-05-19/strtol-error-checking/
// 1: is an int
// 0: is text
int is_int(char *s, int *v)
{
    errno = 0;
    char *endptr;
    long val = strtol(s, &endptr, 10);
    if ( ((errno != 0 ) ||(*endptr != '\0')) || (val>INT_MAX))
        return 0;
    else {
        *v = (int) val;
        return 1;
    }
}
//}}}

//{{{int uint32_t_str_pair_cmp(void *a, void *b)
int uint32_t_str_pair_cmp(const void *a, const void *b)
{
    struct uint32_t_str_pair *_a = (struct uint32_t_str_pair *)a;
    struct uint32_t_str_pair *_b = (struct uint32_t_str_pair *)b;
    return strcmp(_a->str, _b->str);
}
//}}}

//{{{int ped_get_column_names_types(char *ped_file_name,
int ped_get_column_names_types(char *ped_file_name,
                               char ***ped_field_names,
                               int **ped_field_is_int)
{
    uint32_t i, j;
    int num_ped_fields = 0;

    // Figure out which fields will be in the DB
    if (ped_file_name != NULL) {
        FILE *ped_f = fopen(ped_file_name, "r");
        if (!ped_f)
            err(EX_NOINPUT,
                "ERROR ped_get_column_names_types(): "
                "Cannot open file '%s'",
                ped_file_name);

        char *line = NULL, *tmp_line;
        size_t len = 0;

        ssize_t read = getline(&line, &len, ped_f);
        if (read == -1) {
            if (feof(ped_f))
                errx(EX_NOINPUT,
                     "ERROR ped_get_column_names_types(): "
                     "Error reading file '%s': End of file",
                     ped_file_name);
            err(EX_NOINPUT,
                "ERROR ped_get_column_names_types(): "
                "Error reading file '%s'", ped_file_name);
        }

        // Scan the first line to get the names of the fields
        if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = '\0';

        if (line[0] == '#'){
            line++;
        }

        char *word;
        tmp_line = (char *) malloc((strlen(line)+1) * sizeof(char));
        if (!tmp_line )
            err(EX_OSERR,
                "ERROR ped_get_column_names_types(): "
                "malloc error");

        strcpy(tmp_line, line);

        word = strtok(tmp_line, "\t");
        while (word != NULL) {
            num_ped_fields += 1;
            word = strtok(NULL, "\t");
        }

        if (num_ped_fields == 0) {
            errx(EX_NOINPUT,
                 "ERROR ped_get_column_names_types(): "
                 "Empty PED file '%s'.", ped_file_name);
        } else {
            *ped_field_names = (char **) 
                    malloc(num_ped_fields * sizeof(char *));
            if (*ped_field_names == NULL)
                err(EX_OSERR,
                    "ERROR ped_get_column_names_types(): "
                    "malloc error");
            strcpy(tmp_line, line);

            // Set field names
            char *p = strtok(tmp_line, "\t");
            (*ped_field_names)[0] = strdup(p);
            for (i = 1; i < num_ped_fields; ++i) {
                p = strtok(NULL, "\t");
                (*ped_field_names)[i] = strdup(p);
            }

            // Convert " " to "_"
            for (i = 0; i < num_ped_fields; ++i) {
                for (j = 0; j < strlen( (*ped_field_names)[i]); ++j) {
                    if ( (*ped_field_names)[i][j] == ' ')
                         (*ped_field_names)[i][j] = '_';
                }
            }

            // Check for problems with field names
            for (i = 0; i < num_ped_fields; ++i) {
                int r = check_field_name( (*ped_field_names) [i]);
                if (r >= 0) {
                    errx(EX_NOINPUT, 
                         "ERROR ped_get_column_names_types(): "
                         "Invalid character '%c' in field name '%s' from file "
                         "'%s'",
                         (*ped_field_names)[i][r],
                         (*ped_field_names)[i],
                         ped_file_name);
                }
            }

            // Set field types
            *ped_field_is_int = (int *) malloc(num_ped_fields * sizeof(int));
            if (*ped_field_is_int == NULL)
                err(EX_OSERR,
                    "ERROR ped_get_column_names_types(): "
                    "malloc error");

            for (i = 0; i < num_ped_fields; ++i)
                (*ped_field_is_int)[i] = 1;

            uint32_t line_no = 2;

            while ( (read = getline(&line, &len, ped_f)) != -1) {
                if (line[strlen(line) - 1] == '\n')
                    line[strlen(line) - 1] = '\0';

                uint32_t j;
                word = strtok(line, "\t");

                for (i = 0; i < num_ped_fields; ++i) {
                    if (word == NULL) {
                        errx(EX_NOINPUT,
                             "ERROR ped_get_column_names_types(): "
                             "Missing field in file '%s' on line %u.\n",
                             ped_file_name,
                             line_no);
                    }

                    int v;
                    (*ped_field_is_int)[i] &= is_int(word, &v);
                    word = strtok(NULL, "\t");
                }

                // unparsed data on this line
                if (word != NULL)
                    errx(EX_NOINPUT,
                         "ERROR ped_get_column_names_types(): "
                         "Extra field in file '%s' on line %u.\n",
                         ped_file_name,
                         line_no);

                line_no += 1;
            }

            // check to see if no data is read
            if (line_no == 2) 
                errx(EX_NOINPUT, 
                     "ERROR ped_get_column_names_types(): "
                     "No data in PED file '%s'",
                     ped_file_name);

            if (ferror(ped_f) != 0 )
                err(EX_NOINPUT,
                    "ERROR ped_get_column_names_types(): "
                    "Error reading file '%s'",
                    ped_file_name);
        }
        fclose(ped_f);
        free(line);
        free(tmp_line);
    }

    return num_ped_fields;
}
//}}}

//{{{ int ped_create_db(char *ped_db_file_name,
uint32_t ped_create_db(char *ped_file_name,
                       char *ped_db_file_name,
                       char *giggle_index_dir,
                       uint32_t file_name_col)
{
    // Get the indexed files from the giggle index so we can match the sample
    // file from the index with the ped file
    char **names = NULL;
    uint32_t *num_intervals = NULL;
    double *mean_interval_sizes = NULL;
    uint32_t num_files = giggle_get_indexed_files(giggle_index_dir,
                                                  &names,
                                                  &num_intervals,
                                                  &mean_interval_sizes);

    struct uint32_t_str_pair *giggle_names_order =
            (struct uint32_t_str_pair *)
            malloc(num_files * sizeof(struct uint32_t_str_pair));
    uint32_t i;
    for (i = 0; i < num_files; ++i) {
        giggle_names_order[i].uint = i;
        // strip out the path info from names
        giggle_names_order[i].str = strdup(basename(names[i]));
    }

    // ### quick sort, the last parameter is the sort function.
    qsort(giggle_names_order, 
          num_files,
          sizeof(struct uint32_t_str_pair),
          uint32_t_str_pair_cmp);


    // Get ped file headers
    char **ped_field_names = NULL;
    int *ped_field_is_int = NULL;

    int num_ped_fields = ped_get_column_names_types(ped_file_name,
                                                    &ped_field_names,
                                                    &ped_field_is_int);

    char *q_insert_postfix = NULL, *q_insert_postfix_tmp = NULL;
    int ret = asprintf(&q_insert_postfix, "?1");
    if (ret == -1) errx(1, "ERROR ped_create_db(): asprintf error");

    char *q_base_insert = NULL, *q_base_insert_tmp = NULL;
    ret = asprintf(&q_base_insert, "INSERT INTO ped(Giggle_File_id");
    if (ret == -1) errx(1, "ERROR ped_create_db(): asprintf error");

    char *q_create_table = NULL, *q_create_table_tmp = NULL;
    ret = asprintf(&q_create_table,
                   "CREATE TABLE ped(Giggle_File_Id INTEGER");
    if (ret == -1) errx(1, "ERROR ped_create_db(): asprintf error");

    // Add fields from PED to the table
    for (i = 0; i < num_ped_fields; ++i) {
        if (ped_field_is_int[i] == 1) {
            ret = asprintf(&q_create_table_tmp,
                           "%s, %s INTEGER",
                           q_create_table,
                           ped_field_names[i]);
            if (ret == -1) errx(1, "ERROR ped_create_db(): asprintf error");

        } else {
            ret = asprintf(&q_create_table_tmp,
                           "%s, %s TEXT",
                           q_create_table,
                           ped_field_names[i]);
            if (ret == -1) errx(1, "ERROR ped_create_db(): asprintf error");
        }
        free(q_create_table);
        q_create_table = q_create_table_tmp;


        ret = asprintf(&q_base_insert_tmp,
                       "%s, %s",
                       q_base_insert,
                       ped_field_names[i]);
        if (ret == -1) errx(1, "ERROR ped_create_db(): asprintf error");
        free(q_base_insert);
        q_base_insert = q_base_insert_tmp;

        ret = asprintf(&q_insert_postfix_tmp,
                       "%s, ?%u",
                       q_insert_postfix,
                       i+2);
        if (ret == -1) errx(1, "ERROR ped_create_db(): asprintf error");
        free(q_insert_postfix);
        q_insert_postfix = q_insert_postfix_tmp;
    }

    ret = asprintf(&q_create_table_tmp, "%s);", q_create_table);
    if (ret == -1) errx(1, "ERROR ped_create_db(): asprintf error");
    free(q_create_table);
    q_create_table = q_create_table_tmp;

    ret = asprintf(&q_base_insert_tmp,
                   "%s) VALUES (%s);",
                   q_base_insert,
                   q_insert_postfix);
    if (ret == -1) errx(1, "ERROR ped_create_db(): asprintf error");
    free(q_base_insert);
    q_base_insert = q_base_insert_tmp;
    free(q_insert_postfix);

    //fprintf(stderr, "%s\n", q_create_table);
    //fprintf(stderr, "%s\n", q_base_insert);


    // removed DB if it is there
    struct stat buffer;
    ret = stat(ped_db_file_name, &buffer);
    if (ret == 0)
        remove(ped_db_file_name);

    // create the table
    sqlite3 *db;
    char *err_msg = NULL;
    ret = sqlite3_open(ped_db_file_name, &db);
    if( ret != SQLITE_OK )
        err(EX_SOFTWARE,
            "ERROR ped_create_db(): SQL error '%s' for database '%s'",
            err_msg,
            ped_db_file_name);

    //sqlite3_exec(db, "BEGIN TRANSACTION;", NULL, NULL, NULL);

    ret = sqlite3_exec(db, q_create_table, NULL, 0, &err_msg);
    if( ret != SQLITE_OK )
        err(EX_SOFTWARE,
            "ERROR ped_create_db(): SQL error '%s' in query '%s'",
            err_msg,
            q_create_table);


    // Get rows from ped file
    FILE *ped_f = fopen(ped_file_name, "r");
    if (!ped_f)
        err(EX_NOINPUT,
            "ERROR ped_create_db(): "
            "Cannot open file '%s'",
            ped_file_name);

    char *line = NULL;
    size_t len = 0;
    uint32_t line_no = 1;

    // get rid of header
    ssize_t read = getline(&line, &len, ped_f);
    line_no += 1;
    if (read == -1) {
        if (feof(ped_f))
            errx(EX_NOINPUT,
                 "ERROR ped_create_db(): "
                 "Error reading file '%s': End of file",
                 ped_file_name);
        err(EX_NOINPUT,
            "ERROR ped_create_db(): "
            "Error reading file '%s'", ped_file_name);
    }


    uint32_t giggle_file_id = 0;
    sqlite3_stmt *base_insert_stmt = NULL;

    uint32_t inserts = 0;
    while ( (read = getline(&line, &len, ped_f)) != -1) {
        line_no += 1;
        if (line[strlen(line) - 1] == '\n')
            line[strlen(line) - 1] = '\0';

        uint32_t col_i;
        char *word = strtok(line, "\t");

        ret = sqlite3_prepare_v2(db,
                                 q_base_insert,
                                 strlen(q_base_insert),
                                 &base_insert_stmt,
                                 NULL);

        for (col_i = 0; col_i < num_ped_fields; ++col_i) {
            if (word == NULL) {
                errx(EX_NOINPUT,
                     "ERROR ped_create_db(): "
                     "Missing field in file '%s' on line %u.\n",
                     ped_file_name,
                     line_no);
            }

            if (ret != SQLITE_OK) 
                errx(EX_SOFTWARE,
                     "ERROR ped_create_db(): "
                     "Can't prepare insert statment %s (%i): %s",
                     q_base_insert, ret, sqlite3_errmsg(db));


            if (ped_field_is_int[col_i] == 1) {
                ret = sqlite3_bind_int(base_insert_stmt,
                                       col_i + 2,
                                       atoi(word));
                if (ret != SQLITE_OK) 
                    errx(EX_SOFTWARE,
                         "ERROR ped_create_db(): "
                         "Error binding value in insert (%i): %s",
                         ret, sqlite3_errmsg(db));
            } else {
                ret = sqlite3_bind_text(base_insert_stmt,
                                        col_i + 2,
                                        word,
                                        strlen(word),
                                        NULL);
                if (ret != SQLITE_OK) 
                    errx(EX_SOFTWARE,
                         "ERROR ped_create_db(): "
                        "Error binding value in insert (%i): %s",
                        ret, sqlite3_errmsg(db));
            }

            if (col_i+1 == file_name_col) {
                struct uint32_t_str_pair tmp_pair;
                tmp_pair.str = word;

                //struct uint32_t_str_pair *   =

                struct uint32_t_str_pair * giggle_file_i  =
                    bsearch(&tmp_pair,
                            giggle_names_order, 
                            num_files,
                            sizeof(struct uint32_t_str_pair),
                            uint32_t_str_pair_cmp);

                if (giggle_file_i == NULL)
                    errx(1,
                         "ERROR ped_create_db(): "
                         "PED file %s not found in giggle index.\n",
                         word);

                ret = sqlite3_bind_int(base_insert_stmt,
                                       1,
                                       giggle_file_i->uint);
                if (ret != SQLITE_OK) 
                    errx(EX_SOFTWARE,
                         "ERROR ped_create_db(): "
                         "Error binding value in insert (%i): %s",
                         ret, sqlite3_errmsg(db));

            }

            word = strtok(NULL, "\t");

        }
        if (word != NULL)
            errx(EX_NOINPUT,
                 "ERROR ped_create_db(): "
                 "Extra field in file '%s' on line %u.\n",
                 ped_file_name,
                 line_no);

        ret = sqlite3_step(base_insert_stmt);
        if (ret != SQLITE_DONE)
            errx(EX_SOFTWARE,
                 "ERROR ped_create_db(): "
                "Error on insert(%i): %s",
                ret, sqlite3_errmsg(db));

        ret = sqlite3_finalize(base_insert_stmt);
        if (ret != SQLITE_OK)
            errx(EX_SOFTWARE,
                 "ERROR ped_create_db(): "
                "Error finalizing insert(%i): %s",
                ret, sqlite3_errmsg(db));
        inserts += 1;
    }
    if (line_no == 2) 
        errx(EX_NOINPUT, 
             "ERROR ped_create_db(): "
             "No data in PED file '%s'",
             ped_file_name);

    if (ferror(ped_f) != 0 )
        err(EX_NOINPUT,
            "ERROR ped_create_db(): "
            "Error reading file '%s'",
            ped_file_name);

    fclose(ped_f);
    free(line);

    sqlite3_close(db);

    return inserts;
}
//}}}

//{{{ uint32_t ped_get_matching_sample_ids(char *ped_file_name_db,
uint32_t ped_get_matching_sample_ids(char *ped_file_name_db,
                                     char *select_query,
                                     uint32_t **sample_ids)
{
    sqlite3 *db;
    char *err_msg = NULL;

    int ret = sqlite3_open(ped_file_name_db, &db);
    if( ret != SQLITE_OK )
        err(EX_NOINPUT,
            "ERROR ped_get_matching_sample_ids(): "
            "SQL error '%s' for database '%s'",
            err_msg, ped_file_name_db);

    char *test_q;

    ret = asprintf(&test_q,
                   "SELECT Giggle_File_Id FROM ped WHERE %s "
                   "ORDER BY Giggle_File_Id;",
                   select_query);

    if (ret == -1)
        err(EX_OSERR, "ERROR ped_get_matching_sample_ids(): asprintf error");

    struct uint32_t_ll ll;
    ll.head = NULL;
    ll.tail = NULL;
    ll.len = 0;

    ret = sqlite3_exec(db, test_q, uint32_t_ll_callback, &ll, &err_msg);
    if( ret != SQLITE_OK )
        err(EX_SOFTWARE,
            "ERROR ped_get_matching_sample_ids(): "
            "SQL error '%s' in query '%s'", err_msg, test_q);

    if (ll.len > 0) {
        *sample_ids = (uint32_t *) malloc(ll.len * sizeof(uint32_t));
        if (*sample_ids == NULL)
            err(EX_OSERR, "ERROR ped_get_matching_sample_ids(): malloc error");

        struct uint32_t_ll_node *tmp, *curr = ll.head;
        uint32_t i;
        for (i = 0; i < ll.len; ++i) {
            (*sample_ids)[i] = curr->v;
            tmp = curr->next;
            free(curr);
            curr = tmp;
        }
    }
    sqlite3_close(db);
    return ll.len;
}
//}}}

//{{{ uint32_t ped_get_uniq_col_groups(char *ped_file_name_db,
uint32_t ped_get_uniq_col_groups(char *ped_file_name_db,
                                 sqlite3 **db,
                                 char **cols,
                                 uint32_t num_cols,
                                 char *select_query,
                                 char ****uniq_col_values,
                                 uint32_t ***uniq_groups_ids,
                                 uint32_t **uniq_groups_sizes)
{
    char *err_msg = NULL;
    int ret;

    if (*db == NULL) {
        ret = sqlite3_open(ped_file_name_db, db);
        if( ret != SQLITE_OK )
            err(EX_NOINPUT,
                "ERROR ped_get_uniq_col_groups(): "
                "SQL error '%s' for database '%s'",
                err_msg, ped_file_name_db);
    }

    char *col_list = NULL;
    ret = asprintf(&col_list, "%s", cols[0]);
    if (ret == -1)
        err(EX_OSERR, "ERROR ped_get_uniq_col_groups(): asprintf error");


    // agregate cols into a comma sep list
    uint32_t i;
    for (i = 1; i < num_cols; ++i) {
        char *col_list_tmp;
        ret = asprintf(&col_list_tmp, "%s,%s", col_list,cols[i]);
        if (ret == -1)
            err(EX_OSERR, "ERROR ped_get_uniq_col_groups(): asprintf error");

        free(col_list);
        col_list = col_list_tmp;
    }

    // get col uniq values
    char *q_select = NULL;
    if (select_query == NULL )
        ret = asprintf(&q_select,
                       "SELECT %s FROM ped GROUP BY %s",
                       col_list, 
                       col_list);
    else
        ret = asprintf(&q_select,
                       "SELECT %s FROM ped WHERE %s GROUP BY %s",
                       col_list, 
                       select_query,
                       col_list);

    if (ret == -1)
        err(EX_OSERR, "ERROR ped_get_uniq_col_groups(): asprintf error");
    
    struct str_ll s_ll;
    s_ll.head = NULL;
    s_ll.tail = NULL;
    s_ll.len = 0;

    ret = sqlite3_exec(*db, q_select, str_ll_callback, &s_ll, &err_msg);
    if( ret != SQLITE_OK )
        err(EX_SOFTWARE,
            "ERROR ped_get_uniq_col_groups(): "
            "SQL error '%s' in query '%s'", err_msg, q_select);

    uint32_t num_uniq_groups = s_ll.len;

    char ***c = (char ***)malloc(num_uniq_groups * sizeof(char **));
    struct str_ll_node *s_ll_curr = s_ll.head, *s_ll_tmp;

    for (i = 0; i < num_uniq_groups; ++i) {
        c[i] = (char **)malloc(s_ll_curr->num_strs * sizeof(char *));
        uint32_t j;
        for (j = 0; j < s_ll_curr->num_strs; ++j) {
            c[i][j] = strdup(s_ll_curr->str[j]);
            free(s_ll_curr->str[j]);
        }

        free(s_ll_curr->str);

        s_ll_tmp = s_ll_curr;
        s_ll_curr = s_ll_curr->next;
        free(s_ll_tmp);
    }

    *uniq_col_values = c;

    free(q_select);

    // collect the ids of individuals that fit into each group
    *uniq_groups_ids = (uint32_t **)malloc(num_uniq_groups*sizeof(uint32_t *));
    *uniq_groups_sizes = (uint32_t *)malloc(num_uniq_groups*sizeof(uint32_t));

    for (i = 0; i < num_uniq_groups; ++i) {
        // Concat the list of uniq values to one string
        char *concat = NULL, *concat_tmp = NULL;
        ret = asprintf(&concat,"%s=='%s'", cols[0], c[i][0]);
        uint32_t j;
        for (j = 1; j < num_cols;++j) {
            ret = asprintf(&concat_tmp,
                           "%s AND %s=='%s'",
                           concat,
                           cols[j],
                           c[i][j]);

            if (ret == -1)
                err(EX_OSERR,
                    "ERROR ped_get_uniq_col_groups(): asprintf error");
            free(concat);
            concat = concat_tmp;
        }

        if (select_query == NULL)
            ret = asprintf(&q_select,
                    "SELECT Giggle_File_Id FROM ped WHERE %s",
                    concat);
        else
            ret = asprintf(&q_select,
                    "SELECT Giggle_File_Id FROM ped WHERE (%s) AND (%s)",
                    select_query,
                    concat);


        //fprintf(stderr, "%s\n", q_select);

        struct uint32_t_ll u_ll;
        u_ll.head = NULL;
        u_ll.tail = NULL;
        u_ll.len = 0;

        ret = sqlite3_exec(*db,
                           q_select,
                           uint32_t_ll_callback,
                           &u_ll,
                           &err_msg);
        if( ret != SQLITE_OK )
            err(EX_SOFTWARE,
                "ERROR ped_get_uniq_col_groups(): "
                "SQL error '%s' in query '%s'",
                err_msg,
                q_select);

        (*uniq_groups_sizes)[i] = u_ll.len;
        (*uniq_groups_ids)[i] =
                (uint32_t *)malloc((*uniq_groups_sizes)[i]*sizeof(uint32_t));
        struct uint32_t_ll_node *u_ll_curr = u_ll.head;
        for (j = 0; j < (*uniq_groups_sizes)[i]; ++j) {
            (*uniq_groups_ids)[i][j] = u_ll_curr->v;
            struct uint32_t_ll_node *u_ll_tmp = u_ll_curr->next;
            free(u_ll_curr);
            u_ll_curr = u_ll_tmp;
        }

        free(q_select);
    }

    return  num_uniq_groups;
}
//}}}

//{{{uint32_t ped_union_groups(uint32_t num_groups,
uint32_t ped_union_groups(uint32_t num_groups,
                          uint32_t **uniq_groups_ids,
                          uint32_t *uniq_groups_sizes,
                          uint32_t **union_group_ids)
{
    uint32_t i, num_union_group_ids = 0;
    for (i = 0; i < num_groups; ++i)
        num_union_group_ids += uniq_groups_sizes[i];

    *union_group_ids = (uint32_t *)malloc(num_union_group_ids*sizeof(uint32_t));

    uint32_t j, k;
    i = 0;
    for (j = 0; j < num_groups; ++j) {
        for (k = 0; k < uniq_groups_sizes[j]; ++k) {
            (*union_group_ids)[i] = uniq_groups_ids[j][k];
            i+=1;
        }
    }

    return num_union_group_ids;
}
//}}}

//{{{ uint32_t ped_get_cols_info_by_id(char *ped_file_name_db,
uint32_t ped_get_cols_info_by_id(char *ped_file_name_db,
                                 sqlite3 **db,
                                 char **cols,
                                 uint32_t num_cols,
                                 uint32_t sample_id,
                                 char ***col_vals,
                                 char ***col_names)
{
    char *err_msg = NULL;
    int ret;

    if (*db == NULL) {
        ret = sqlite3_open(ped_file_name_db, db);
        if( ret != SQLITE_OK )
            err(EX_NOINPUT,
                "ERROR ped_get_cols_info_by_id(): "
                "SQL error '%s' for database '%s'",
                err_msg, ped_file_name_db);
    }


    char *all_cols = NULL;
    if (num_cols > 0) {
        char *all_cols_tmp; 
        all_cols = strdup(cols[0]);
        uint32_t i;
        for(i = 1; i < num_cols; ++i) {
            ret = asprintf(&all_cols_tmp,
                           "%s,%s",
                           all_cols,
                           cols[i]);
            if (ret == -1)
                err(EX_OSERR,
                    "ERROR ped_get_cols_info_by_id(): "
                    "asprintf error");
            free(all_cols);
            all_cols = all_cols_tmp;
        }
    } else {
        all_cols = strdup("*");
    }

    char *q_select;
    ret = asprintf(&q_select,
                   "SELECT %s FROM ped WHERE "
                   "Giggle_File_Id == %u",
                   all_cols,
                   sample_id);
    if (ret == -1)
        err(EX_OSERR, 
            "ERROR ped_get_cols_info_by_id(): "
            "asprintf error");

    struct str_col_ll s_ll;
    s_ll.head = NULL;
    s_ll.tail = NULL;
    s_ll.len = 0;

    ret = sqlite3_exec(*db, q_select, str_col_ll_callback, &s_ll, &err_msg);
    if( ret != SQLITE_OK )
        err(EX_SOFTWARE,"SQL error '%s' in query '%s'", err_msg, q_select);

    if (s_ll.len != 1)
        errx(1,
             "ERROR ped_get_cols_info_by_id(): "
             "Expected one row, %u returned. %s\n",
             s_ll.len,
             q_select);

    *col_vals = (char **)malloc(s_ll.head->num_strs * sizeof(char *));
    *col_names = (char **)malloc(s_ll.head->num_strs * sizeof(char *));

    uint32_t i;
    for (i = 0; i < s_ll.head->num_strs; ++i) {
        (*col_vals)[i] = strdup(s_ll.head->str[i]);
        (*col_names)[i] = strdup(s_ll.head->col[i]);
        free(s_ll.head->str[i]);
        free(s_ll.head->col[i]);
    }

    free(s_ll.head->str);
    free(s_ll.head->col);
    uint32_t len = s_ll.head->num_strs;
    free(s_ll.head);

    return len;
}
//}}}
