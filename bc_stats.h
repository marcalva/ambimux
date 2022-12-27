
#ifndef BC_STATS_H
#define BC_STATS_H

#include "inttypes.h"
#include "htslib/khash.h"
#include "atac_data.h"
#include "rna_data.h"
#include "str_util.h"

// hash table of flags to see if index is present
KHASH_INIT(kh_cnt, int32_t, char, 0, kh_int_hash_func, kh_int_hash_equal);

typedef struct bc_counts {
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
} bc_counts;

int bc_counts_init(bc_counts *bcc);
bc_counts *bc_counts_alloc();
void bc_counts_free(bc_counts *bcc);
void bc_counts_dstry(bc_counts *bcc);
/* Create a copy of bc_counts
 *
 * @param bcc Pointer to bc_counts object to copy
 * @param n Number of elements in bcc.
 * @return Array of bc_counts object. Must be freed by caller.
 */
bc_counts *bc_counts_copy(const bc_counts *bcc, size_t n);

KHASH_INIT(kh_bcs, char *, bc_counts *, 1, kh_str_hash_func, kh_str_hash_equal);

typedef struct {
    khash_t(kh_bcs) *kh_bcc;
    str_map *barcodes;
    bc_counts **bcc; // internal array for sorting.
    size_t n_bcc;
} bcs_stats;

bcs_stats *bc_stats_alloc();
void bcs_stats_free(bcs_stats *bcss);
int bcs_stats_fill(bcs_stats *bcss, bam_rna_t *br, bam_atac_t *ba);
int bcs_stats_sort(bcs_stats *bcss, int ascend);
uint32_t bcs_stats_min(bcs_stats *bcss, uint32_t top_n, int *ret);

void bcs_stats_print_sorted(bcs_stats *bcss, int64_t max_n);

static inline
int cmp_bc_counts(const void *bc1, const void *bc2){
    bc_counts *bcc1 = *(bc_counts **)bc1;
    bc_counts *bcc2 = *(bc_counts **)bc2;

    if (bcc1->counts < bcc2->counts) return(-1);
    if (bcc1->counts > bcc2->counts) return(1);
    return(0);
}

static inline
int cmp_bc_counts_rev(const void *bc1, const void *bc2){
    bc_counts *bcc1 = *(bc_counts **)bc1;
    bc_counts *bcc2 = *(bc_counts **)bc2;

    if (bcc1->counts > bcc2->counts) return(-1);
    if (bcc1->counts < bcc2->counts) return(1);
    return(0);
}

#endif // BC_STATS_H

