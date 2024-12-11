#ifndef __SHARDING_H__
#define __SHARDING_H__

// struct shard {
//     char* giggle_path,
//     char* stixdb_path,
//     char* bed_path
// }

#define MAX_NAME_LENGTH 512

typedef struct shard {
    char *giggle_path;
    char *stixdb_path;
    // char* bed_path;

} Shard;

Shard *read_shards_from_file(const char *filename, int *shard_len);

typedef struct tablequery {
    char *left_str;
    char *right_str;
    char *len;
    char *svtype;
    char *ID;
    // char* bed_path;

} TableQuery;

TableQuery *read_table_queries_from_file(const char *filename, int *table_len);

#endif
