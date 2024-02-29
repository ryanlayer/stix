#include <stdio.h>
#include <strings.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <sysexits.h>
#include <err.h>
#include <stdint.h>
#include "sharding_utils.h"

Shard *read_shards_from_file(const char *filename)
{
    Shard *shards_array = (Shard *)malloc(sizeof(Shard)); // array of Shards

    if (filename != NULL)
    {

        FILE *fhandler = fopen(filename, "r");
        if (!fhandler)
            err(EX_NOINPUT,
                "ERROR read_shards_from_file(): "
                "Cannot open file '%s'",
                filename);

        char *line = NULL, *tmp_line;
        size_t len = 0;
        // ssize_t read = getline(&line, &len, fhandler);
        ssize_t read;
        // if (read == -1) {
        //     if (feof(fhandler))
        //         errx(EX_NOINPUT,
        //              "ERROR read_shards_from_file(): "
        //              "Error reading file '%s': End of file",
        //              filename);
        //     err(EX_NOINPUT,
        //         "ERROR read_shards_from_file(): "
        //         "Error reading file '%s'", filename);
        // }

        int num_ped_fields = 3;
        int line_no = 0;
        char *word;
        while ((read = getline(&line, &len, fhandler)) != -1)
        {

            line_no++;
            Shard *shards_array_temp = shards_array_temp;
            if (line_no >= 2)
            {
                Shard *shards_array_temp = realloc(shards_array, line_no * sizeof(Shard));

                if (shards_array_temp == NULL)
                {
                    // fprintf(stderr, "Memory reallocation failed\n");
                    free(shards_array); // 在失败前确保释放原始内存
                    errx(EX_TEMPFAIL,
                         "ERROR read_shards_from_file(): "
                         "Expand shards_array failed on line %u.\n",
                         line_no);
                }
            }

            if (line[strlen(line) - 1] == '\n')
                line[strlen(line) - 1] = '\0';
            if (line[0] == '#')
            {
                continue;
            }

            uint32_t j;
            word = strtok(line, "\t");
            Shard shard;
            for (uint32_t i = 0; i < num_ped_fields; ++i)
            {
                // fprintf(stderr, "word:%s,i=%d\n", word, i);
                if (word == NULL)
                {
                    errx(EX_NOINPUT,
                         "ERROR read_shards_from_file(): "
                         "Missing field in file '%s' on line %u.\n",
                         filename,
                         line_no);
                }
                fprintf(stderr, "start_switch");
                switch (i)
                {
                case 0:
                    shard.giggle_path = word;
                    break;
                case 1:
                    shard.stixdb_path = word;
                    break;
                case 2:
                    shard.bed_path = word;
                    break;
                default:
                    break;
                }
                word = strtok(NULL, "\t");
            }
            shards_array_temp[line_no - 1] = shard;
            // free(shards_array);
            shards_array = shards_array_temp;
        }
    }

    return shards_array;
}