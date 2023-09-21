
#ifndef BC_STATS_H
#define BC_STATS_H

#include "inttypes.h"
#include "htslib/khash.h"
#include "atac_data.h"
#include "rna_data.h"
#include "str_util.h"

// hash table of flags to see if index is present
// key is index, no value since hash is not a map/dict.
KHASH_INIT(kh_cnt, int32_t, char, 0, kh_int_hash_func, kh_int_hash_equal);

typedef struct bc_stats_t {
    char *bc;
    uint32_t counts; // total number of UMIs + fragments
    uint32_t atac_counts; // total number of atac fragments
    uint32_t rna_counts; // total number of rna UMIs

    uint32_t n_atac_vars; // total number of atac variants detected
    uint32_t n_rna_vars; // total number of rna variants detected

    uint32_t n_gene; // total number of rna genes detected
    uint32_t n_peak; // total number of atac peaks detected

    float frip; // fraction of reads in peaks
    float frig; // fraction of reads in genes

    float rna_mt; // percent of rna umis from mitochondria
    float atac_mt; // percent of atac fragments from mitochondria

    // used to store how many genes, peaks, and variants have been detected.
    khash_t(kh_cnt) *genes;
    khash_t(kh_cnt) *atac_vars;
    khash_t(kh_cnt) *rna_vars;
    khash_t(kh_cnt) *peaks;
} bc_stats_t;

int bc_stats_init(bc_stats_t *bcc);
bc_stats_t *bc_stats_alloc();
void bc_stats_free(bc_stats_t *bcc);
void bc_stats_dstry(bc_stats_t *bcc);

#endif // BC_STATS_H

