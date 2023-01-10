
#ifndef BC_STATS_H
#define BC_STATS_H

#include "inttypes.h"
#include "htslib/khash.h"
#include "atac_data.h"
#include "rna_data.h"
#include "str_util.h"

// hash table of flags to see if index is present
KHASH_INIT(kh_cnt, int32_t, char, 0, kh_int_hash_func, kh_int_hash_equal);

typedef struct bc_stats_t {
    char *bc;
    uint32_t counts;
    uint32_t atac_counts;
    uint32_t rna_counts;
    khash_t(kh_cnt) *genes;
    khash_t(kh_cnt) *vars;
    khash_t(kh_cnt) *peaks;
    uint32_t n_gene;
    uint32_t n_var;
    uint32_t n_peak;
    float frip;
} bc_stats_t;

int bc_stats_init(bc_stats_t *bcc);
bc_stats_t *bc_stats_alloc();
void bc_stats_free(bc_stats_t *bcc);
void bc_stats_dstry(bc_stats_t *bcc);

#endif // BC_STATS_H

