
#ifndef BAM_DAT_H
#define BAM_DAT_H

#include <stdlib.h>
#include "rna_data.h"
#include "atac_data.h"
#include "variants.h"
#include "region.h"
#include "bc_stats.h"

/*! @typedef
 * @abstract Store barcode data
 */
typedef struct {
    // RNA
    khash_t(khrmn) *rna_mols; // indexed by UMI

    // ATAC
    khash_t(khap) *atac_pairs; // indexed by query name hash
    khash_t(khaf) *atac_frags; // indexed by location

    // barcode stats/counts
    bc_stats_t *bc_stats;

} bc_data_t;

// hash table for bc_data_t
// key is barcode string, value is pointer to bc_data_t
KHASH_INIT(kh_bc_dat, char *, bc_data_t *, 1, kh_str_hash_func, kh_str_hash_equal);

bc_data_t *bc_data_init();
void bc_data_dstry(bc_data_t *bcdat);

/* Add RNA read to bc_data_t
 * 
 * Adds the read @p r to @f rna_dups.
 * A copy is created and added.
 *
 * @return 0 on success, -1 on error.
 */
int bc_data_add_rna_read(bc_data_t *bcdat, const rna_read1_t *r, const char *name);

/* De-duplicate RNA reads
 *
 * De-duplicates the reads in @f rna_dups to form @f rna_mols.
 *
 * @return 0 on success, -1 on error.
 */
int bc_data_rna_dedup(bc_data_t *bcdat);

/* Call RNA variants in bc_data_t struct
 */
int bc_data_rna_var_call(bc_data_t *bcdat, g_var_t *gv, contig_map *cmap, 
        uint8_t min_qual);

/* Add ATAC read to bc_data_t
 *
 * @return 0 on success, -1 on error.
 */
int bc_data_add_atac_read(bc_data_t *bcdat, const atac_read1_t *ar, qshort qname);

/* De-duplicate the PCR duplicates and form frags.
 *
 * Deduplicate the reads in each frag to form a deduplicated frag.
 * This will remove destroy the underlying atac_dups_t struct and free 
 * its associated memory, leaving dups equal to null. See 'atac_frag_dedup' 
 * for more detail.
 *
 * The duplicated read pairs are present in the frags object. If a 
 * read pair is chimeric, do nothing with the read pair. Otherwise, 
 * deduplicate it.
 *
 * @return 0 on success, -1 on error.
 */
int bc_data_dedup_atac(bc_data_t *bcdat);

/* Call ATAC variants in bc_data
 */
int bc_data_atac_var_call(bc_data_t *bcdat, g_var_t *gv, contig_map *cmap, 
        uint8_t min_qual);

/* Call ATAC peaks in bc_data
 */
int bc_data_atac_peak_call(bc_data_t *bcdat, iregs_t *pks, contig_map *cmap);

/*******************************************************************************
 * bam_data_t
 ******************************************************************************/

/*! @typedef
 * @abstract Structure to hold the atac and rna pileup
 */
typedef struct {
    khash_t(kh_bc_dat) *bc_data;

    // flags
    uint8_t has_rna;
    uint8_t has_atac;
    uint8_t has_stats;

    str_map *bcs;
} bam_data_t;

bam_data_t *bam_data_init();
void bam_data_dstry(bam_data_t *bam_dat);
int bam_data_atac_free_pairs(bam_data_t *bam_data);
int bam_data_rna_add_read(bam_data_t *bam_data, const char *bc, 
        const rna_read1_t *r, const char *name);
int bam_data_rna_dedup(bam_data_t *bam_data);
int bam_data_rna_var_call(bam_data_t *bam_data, g_var_t *gv, 
        contig_map *cmap, uint8_t min_qual);
int bam_data_atac_add_read(bam_data_t *bam_data, const char *bc, 
        const atac_read1_t *r, qshort qname);
int bam_data_atac_dedup(bam_data_t *bam_data);
int bam_data_atac_var_call(bam_data_t *bam_data, g_var_t *gv, 
        contig_map *cmap, uint8_t min_qual);
int bam_data_atac_peak_call(bam_data_t *bam_data, iregs_t *pks, 
        contig_map *cmap);

int bam_data_fill_bcs(bam_data_t *b, str_map *bcs);
int bam_data_fill_stats(bam_data_t *bam_data);

typedef struct {
    char *bc;
    uint32_t count;
} bc_count_t;

static inline
int cmp_bc_count_rev(const void *bc1, const void *bc2){
    bc_count_t *bcc1 = (bc_count_t *)bc1;
    bc_count_t *bcc2 = (bc_count_t *)bc2;

    if (bcc1->count > bcc2->count) return(-1);
    if (bcc1->count < bcc2->count) return(1);
    return(0);
}

uint32_t bam_data_count_of_n(bam_data_t *bam_data, int top_n, int *ret);

#endif // BAM_DAT_H
