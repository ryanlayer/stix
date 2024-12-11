#include "sharding_utils.h"
#include <ctype.h>
#include <err.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <sysexits.h>

Shard *read_shards_from_file(const char *filename, int *shard_len) {
    Shard *shards_array = (Shard *)malloc(sizeof(Shard)); // array of Shards

    *shard_len = 0;
    if (filename != NULL) {

        FILE *fhandler = fopen(filename, "r");
        if (!fhandler)
            err(EX_NOINPUT,
                "ERROR read_shards_from_file(): "
                "Cannot open file '%s'",
                filename);

        char *line = NULL, *tmp_line;
        size_t len = 0;
        ssize_t read;
        int num_ped_fields = 2;
        char *word;
        while ((read = getline(&line, &len, fhandler)) != -1) {
            (*shard_len)++;
            // Shard *shards_array_temp = shards_array;
            if (*shard_len >= 2) {
                Shard *shards_array_temp =
                    realloc(shards_array, (*shard_len) * sizeof(Shard));

                if (shards_array_temp == NULL) {
                    // fprintf(stderr, "Memory reallocation failed\n");
                    free(shards_array); // 在失败前确保释放原始内存
                    errx(
                        EX_TEMPFAIL,
                        "ERROR read_shards_from_file(): "
                        "Expand shards_array failed on line %u.\n",
                        *shard_len
                    );
                }
                shards_array = shards_array_temp;
            }

            if (line[strlen(line) - 1] == '\n')
                line[strlen(line) - 1] = '\0';
            if (line[0] == '#') {
                continue;
            }

            uint32_t j;
            word = strtok(line, "\t");
            Shard shard;
            for (uint32_t i = 0; i < num_ped_fields; ++i) {
                // fprintf(stderr, "word:%s,i=%d\n", word, i);
                if (word == NULL) {
                    errx(
                        EX_NOINPUT,
                        "ERROR read_shards_from_file(): "
                        "Missing field in file '%s' on line %u.\n",
                        filename,
                        *shard_len
                    );
                }
                // fprintf(stderr, "start_switch\n");
                // fprintf(stderr,"%s,%s,%s,%d\n",shard.giggle_path,shard.stixdb_path,shard.bed_path,*shard_len);
                switch (i) {
                case 0:
                    shard.giggle_path = malloc(strlen(word) + 1);
                    if (shard.giggle_path) {
                        strcpy(shard.giggle_path, word);
                    } else {
                        for (int i = 0; i < *shard_len - 1; ++i) {
                            free(shards_array[i].giggle_path);
                            free(shards_array[i].stixdb_path);
                        }
                        free(shards_array);
                        free(line);
                        fclose(fhandler);
                        errx(EX_TEMPFAIL, "Memory allocation failed for giggle_path");
                    }
                    break;
                case 1:
                    shard.stixdb_path = malloc(strlen(word) + 1);
                    if (shard.stixdb_path) {

                        strcpy(shard.stixdb_path, word);
                    } else {
                        for (int i = 0; i < *shard_len - 1; ++i) {
                            free(shards_array[i].giggle_path);
                            free(shards_array[i].stixdb_path);
                        }
                        free(shards_array);
                        free(line);
                        fclose(fhandler);
                        errx(EX_TEMPFAIL, "Memory allocation failed for giggle_path");
                    }
                    break;
                // case 2:
                //     shard.bed_path = malloc(strlen(word) + 1);
                //     if (shard.bed_path) strcpy(shard.bed_path, word);
                //     break;
                default:
                    errx(
                        EX_NOINPUT,
                        "ERROR read_shards_from_file(): "
                        "Ignore extra field in file '%s' on line %u.\n",
                        filename,
                        *shard_len
                    );
                    break;
                }
                word = strtok(NULL, "\t");
            }
            // fprintf(stderr,">>%s,%s,%s,%d\n",shard.giggle_path,shard.stixdb_path,shard.bed_path,*shard_len);
            shards_array[*shard_len - 1] = shard;
            // free(shards_array);
            // shards_array = shards_array_temp;
        }
        // fhandler
        // shard_len = line_no;
        free(line);
        fclose(fhandler);
    }

    // fprintf(stderr,">>>>%s,%s,%s,%d\n",
    //     shards_array[0].giggle_path, shards_array[0].stixdb_path,
    //     shards_array[0].bed_path,*shard_len);
    return shards_array;
}

/*fromxinchang
read batch queries in table format
*/
TableQuery *read_table_queries_from_file(const char *filename, int *table_len) {
    TableQuery *tablequery_array = (TableQuery *)malloc(sizeof(TableQuery));

    *table_len = 0;
    if (filename != NULL) {

        FILE *fhandler = fopen(filename, "r");
        if (!fhandler)
            err(EX_NOINPUT,
                "ERROR read_table_queries_from_file(): "
                "Cannot open file '%s'",
                filename);

        char *line = NULL, *tmp_line;
        size_t len = 0;
        ssize_t read;
        int num_ped_fields = 5;
        char *word;
        while ((read = getline(&line, &len, fhandler)) != -1) {
            (*table_len)++;
            if (*table_len >= 2) {
                TableQuery *tablequery_array_temp =
                    realloc(tablequery_array, (*table_len) * sizeof(TableQuery));

                if (tablequery_array_temp == NULL) {
                    // fprintf(stderr, "Memory reallocation failed\n");
                    free(tablequery_array); // 在失败前确保释放原始内存
                    errx(
                        EX_TEMPFAIL,
                        "ERROR read_table_queries_from_file(): "
                        "Expand tablequery_array failed on line %u.\n",
                        *table_len
                    );
                }
                tablequery_array = tablequery_array_temp;
            }

            if (line[strlen(line) - 1] == '\n')
                line[strlen(line) - 1] = '\0';
            if (line[0] == '#') {

                for (int i = 0; i < *table_len; i++) {
                    free(tablequery_array[i].left_str);
                    free(tablequery_array[i].right_str);
                    free(tablequery_array[i].len);
                    free(tablequery_array[i].svtype);
                    free(tablequery_array[i].ID);
                }
                free(tablequery_array);
                continue;
            }

            uint32_t j;
            word = strtok(line, "\t");
            TableQuery tablequery;
            for (uint32_t i = 0; i < num_ped_fields; ++i) {
                // fprintf(stderr, "word:%s,i=%d\n", word, i);
                if (word == NULL) {
                    errx(
                        EX_NOINPUT,
                        "ERROR read_table_queries_from_file(): "
                        "Missing field in file '%s' on line %u.\n",
                        filename,
                        *table_len
                    );
                }
                // fprintf(stderr, "start_switch\n");
                // fprintf(stderr,"%s,%s,%s,%d\n",shard.giggle_path,shard.stixdb_path,shard.bed_path,*shard_len);
                switch (i) {
                case 0:
                    tablequery.left_str = malloc(strlen(word) + 1);
                    if (tablequery.left_str) {

                        strcpy(tablequery.left_str, word);
                    } else {
                        if (!tablequery.left_str) {

                            free(tablequery.right_str);
                            free(tablequery.len);
                            free(tablequery.svtype);
                            free(tablequery.ID);
                            errx(
                                EX_TEMPFAIL,
                                "ERROR: Memory allocation failed for left_str."
                            );
                        }
                    }
                    break;
                case 1:
                    tablequery.right_str = malloc(strlen(word) + 1);
                    if (tablequery.right_str) {
                        strcpy(tablequery.right_str, word);
                    } else {
                        if (!tablequery.right_str) {

                            free(tablequery.left_str);
                            free(tablequery.len);
                            free(tablequery.svtype);
                            free(tablequery.ID);
                            errx(
                                EX_TEMPFAIL,
                                "ERROR: Memory allocation failed for left_str."
                            );
                        }
                    }
                    break;
                case 2:
                    tablequery.len = malloc(strlen(word) + 1);
                    if (tablequery.len) {
                        strcpy(tablequery.len, word);
                    } else {
                        if (!tablequery.len) {

                            free(tablequery.left_str);
                            free(tablequery.right_str);
                            free(tablequery.svtype);
                            free(tablequery.ID);
                            errx(
                                EX_TEMPFAIL,
                                "ERROR: Memory allocation failed for left_str."
                            );
                        }
                    }
                    break;
                case 3:
                    tablequery.svtype = malloc(strlen(word) + 1);
                    if (tablequery.svtype) {
                        strcpy(tablequery.svtype, word);
                    } else {
                        if (!tablequery.svtype) {

                            free(tablequery.left_str);
                            free(tablequery.right_str);
                            free(tablequery.len);
                            free(tablequery.ID);
                            errx(
                                EX_TEMPFAIL,
                                "ERROR: Memory allocation failed for left_str."
                            );
                        }
                    }
                    break;
                case 4:
                    tablequery.ID = malloc(strlen(word) + 1);
                    if (tablequery.ID) {
                        strcpy(tablequery.ID, word);
                    } else {
                        if (!tablequery.ID) {

                            free(tablequery.left_str);
                            free(tablequery.right_str);
                            free(tablequery.len);
                            free(tablequery.svtype);
                            errx(
                                EX_TEMPFAIL,
                                "ERROR: Memory allocation failed for left_str."
                            );
                        }
                    }
                    break;
                default:
                    errx(
                        EX_NOINPUT,
                        "ERROR read_shards_from_file(): "
                        "Ignore extra field in file '%s' on line %u.\n",
                        filename,
                        *table_len
                    );
                    break;
                }
                word = strtok(NULL, "\t");
            }
            tablequery_array[*table_len - 1] = tablequery;
        }
        fclose(fhandler);
    }
    return tablequery_array;
}
