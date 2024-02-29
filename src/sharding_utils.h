#ifndef __SHARDING_H__
#define __SHARDING_H__

// struct shard {
//     char* giggle_path,
//     char* stixdb_path,
//     char* bed_path    
// }

#define MAX_NAME_LENGTH 512

typedef struct shard {
    char* giggle_path;
    char* stixdb_path;
    char* bed_path;

} Shard;

Shard* read_shards_from_file(const char *filename);


// uint32_t stix_run_giggle_query_sharding(Shard* sharding_array,
//                                enum stix_sv_type sv_type,
//                                struct stix_breakpoint *q_left_bp,
//                                struct stix_breakpoint *q_right_bp,
//                                uint32_t slop,
//                                uint32_t ins_padding,
//                                uint32_t *sample_ids,
//                                uint32_t num_samples,
//                                struct uint_pair **sample_alt_depths);


#endif