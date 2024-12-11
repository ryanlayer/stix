#ifndef __SEARCH_H__
#define __SEARCH_H__

#include <giggle/src/giggle_index.h>
#include <giggle/src/lists.h>
#include <stdint.h>

#include <htslib/hfile.h>
#include <htslib/vcf.h>

uint32_t safe_sub(uint32_t a, uint32_t b);

char *stix_sv_type_strings[5];

enum stix_sv_type { DEL, DUP, INS, INV, BND };

struct stix_breakpoint {
    char *chrm;
    uint32_t start, end;
    int32_t strand;
};

uint32_t stix_parse_result(
    char *result,
    struct stix_breakpoint **left,
    struct stix_breakpoint **right,
    uint32_t *evidence_type
);

uint32_t stix_run_giggle_query(
    struct giggle_index **gi,
    char *giggle_index_dir,
    enum stix_sv_type sv_type,
    struct stix_breakpoint *q_left_bp,
    struct stix_breakpoint *q_right_bp,
    uint32_t slop,
    uint32_t ins_padding,
    uint32_t *sample_ids,
    uint32_t num_samples,
    struct uint_pair **sample_alt_depths,
    int mode
);

uint32_t stix_check_sv(
    struct stix_breakpoint *q_left_bp,
    struct stix_breakpoint *q_right_bp,
    struct stix_breakpoint *in_left_bp,
    struct stix_breakpoint *in_right_bp,
    uint32_t evidence_type,
    uint32_t slop,
    uint32_t ins_padding,
    enum stix_sv_type sv_type
);

uint32_t stix_check_inv(
    struct stix_breakpoint *q_left_bp,
    struct stix_breakpoint *q_right_bp,
    struct stix_breakpoint *in_left_bp,
    struct stix_breakpoint *in_right_bp,
    uint32_t evidence_type,
    uint32_t slop
);

uint32_t stix_check_del(
    struct stix_breakpoint *q_left_bp,
    struct stix_breakpoint *q_right_bp,
    struct stix_breakpoint *in_left_bp,
    struct stix_breakpoint *in_right_bp,
    uint32_t evidence_type,
    uint32_t slop
);

uint32_t stix_check_ins(
    struct stix_breakpoint *q_left_bp,
    struct stix_breakpoint *q_right_bp,
    struct stix_breakpoint *in_left_bp,
    struct stix_breakpoint *in_right_bp,
    uint32_t evidence_type,
    uint32_t slop,
    uint32_t ins_padding
);

uint32_t stix_check_bnd(
    struct stix_breakpoint *q_left_bp,
    struct stix_breakpoint *q_right_bp,
    struct stix_breakpoint *in_left_bp,
    struct stix_breakpoint *in_right_bp,
    uint32_t evidence_type,
    uint32_t slop
);

uint32_t stix_check_dup(
    struct stix_breakpoint *q_left_bp,
    struct stix_breakpoint *q_right_bp,
    struct stix_breakpoint *in_left_bp,
    struct stix_breakpoint *in_right_bp,
    uint32_t evidence_type,
    uint32_t slop
);

struct stix_breakpoint *stix_region_to_breakpoint(char *region);

uint32_t stix_get_uniq(uint32_t *full, uint32_t num_full, uint32_t **uniq);

int32_t stix_bsearch_seq(
    uint32_t key,
    uint32_t *D,
    uint32_t D_size,
    int32_t lo,
    int32_t hi
);

uint32_t stix_get_quartile_counts(
    uint32_t *full,
    uint32_t num_full,
    uint32_t *Q1,
    uint32_t *Q2,
    uint32_t *Q3,
    int32_t *counts
);

uint32_t stix_get_summary(
    struct uint_pair *sample_alt_depths,
    uint32_t *sample_ids,
    uint32_t num_samples,
    int32_t *zero_count,
    int32_t *one_count,
    uint32_t *Q1,
    uint32_t *Q2,
    uint32_t *Q3,
    uint32_t *min,
    uint32_t *max,
    int32_t *counts
);

uint32_t stix_get_summary_shard(
    struct uint_pair *sample_alt_depths_all[],
    uint32_t *sample_ids_all[],
    uint32_t num_samples_all[],
    uint32_t num_samples_size,
    int32_t *zero_count,
    int32_t *one_count,
    uint32_t *Q1,
    uint32_t *Q2,
    uint32_t *Q3,
    uint32_t *min,
    uint32_t *max,
    int32_t *counts
);

uint32_t stix_get_vcf_breakpoints(
    htsFile *fp,
    bcf_hdr_t *hdr,
    bcf1_t *line,
    struct stix_breakpoint *left,
    struct stix_breakpoint *right,
    enum stix_sv_type *sv_type
);

uint32_t stix_get_sample_depths(
    struct uint_pair *sample_alt_depths,
    uint32_t *sample_ids,
    uint32_t num_samples,
    uint32_t **sample_depths
);

#endif
