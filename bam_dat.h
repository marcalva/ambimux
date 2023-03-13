
#ifndef BAM_DAT_H
#define BAM_DAT_H

#include <stdlib.h>
#include <pthread.h>
#include "rna_data.h"
#include "atac_data.h"
#include "variants.h"
#include "region.h"
#include "bc_stats.h"
#include "htslib/sam.h"

/*******************************************************************************
 * bc_data_t
 ******************************************************************************/

/*! @typedef
 * @abstract Store barcode data
 * @field rna_mols hash table of RNA molecules, indexed by UMI char*
 * @field atac_pairs hash table of ATAC read pairs, indexed by qshort (uint64)
 * @field atac_frags hash table of ATAC fragments, indexed by genomic location
 * @field bc_stats array of barcode stats
 */
typedef struct {
    // RNA
    khash_t(khrmn) *rna_mols; // indexed by UMI

    // ATAC
    khash_t(khap) *atac_pairs; // indexed by query name hash
    khash_t(khaf) *atac_frags; // indexed by location

    // list implementations for traversal
    // these are just placeholders to avoid traversing the hash table, 
    // do not free the pointers in the nodes themselves
    ml_t(mlar) mols_l;
    ml_t(mlaf) frags_l;

    // barcode stats/counts
    bc_stats_t *bc_stats;

    pthread_mutex_t bc_lock;
} bc_data_t;

// hash table for bc_data_t
// key is barcode string, value is pointer to bc_data_t
KHASH_INIT(kh_bc_dat, char *, bc_data_t *, 1, kh_str_hash_func, kh_str_hash_equal);

/* Initialize empty bc_data_t struct.
 * Hash tables are initialized.
 * Return null on error.
 */
bc_data_t *bc_data_init();

/* Destroy and free all memory of bcdat
 */
void bc_data_dstry(bc_data_t *bcdat);
/* Free pairs in bcdat struct
 */
void bc_data_free_atac_pairs(bc_data_t *bcdat);

int bc_data_fill_list(bc_data_t *bcdat);

/*******************************************************************************
 * BC RNA
 ******************************************************************************/

/* Add RNA read rna_read1_t to bc_data_t.
 * See bam_data_rna_add_read.
 * @return 0 on success, -1 on error.
 */
int bc_data_rna_add_read(bc_data_t *bcdat, const rna_read1_t *r, umishort umih);

/* De-duplicates the reads in @f rna_mols
 * See bam_data_rna_dedup.
 * @return 0 on success, -1 on error.
 */
int bc_data_rna_dedup(bc_data_t *bcdat);

/* Call RNA variants in bc_data_t struct
 * See bam_data_rna_var_call.
 * @return -1 on error, or the number of variants called on success
 */
int bc_data_rna_var_call(bc_data_t *bcdat, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual);

/*******************************************************************************
 * BC ATAC
 ******************************************************************************/

/* Add ATAC read to bc_data_t
 * See bam_data_atac_add_read.
 * @return 0 on success, -1 on error.
 */
int bc_data_atac_add_read(bc_data_t *bcdat, const atac_read1_t *ar, qshort qname);

/* De-duplicate the PCR duplicates and form frags.
 * See bam_data_atac_dedup.
 * @return 0 on success, -1 on error.
 */
int bc_data_atac_dedup(bc_data_t *bcdat);

/* Call ATAC variants in bc_data
 * see bam_data_atac_var_call.
 * @return -1 on error, or the number of variants called on success
 */
int bc_data_atac_var_call(bc_data_t *bcdat, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual);

/* Call ATAC peaks in bc_data
 * See bam_data_atac_peak_call
 * @return -1 on error, or the number of peaks called on success
 */
int bc_data_atac_peak_call(bc_data_t *bcdat, iregs_t *pks, str_map *cmap);

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

    pthread_mutex_t bam_lock;
} bam_data_t;

// Initialize empty bam_data object
// Return NULL on error, or pointer to the allocated object
bam_data_t *bam_data_init();

// Free all memory in and including @p bam_dat
void bam_data_dstry(bam_data_t *bam_dat);

// Free the ATAC read pairs in bam_data
int bam_data_atac_free_pairs(bam_data_t *bam_data);

/* Set barcodes in bam data.
 *
 * If @f b->bcs is present, they are removed.
 * If @p bcs is given, these barcodes are first added to b->bcs.
 * Then, barcodes in @f b->bc_data are added if not found in 
 * @p bcs.
 *
 * @return -1 on error, 0 on success.
 */
int bam_data_fill_bcs(bam_data_t *b, str_map *bcs);

/*******************************************************************************
 * BAM RNA
 ******************************************************************************/

/* Add RNA read rna_read1_t to bam_data.
 * 
 * Adds the read @p r to the @f bc_data and @f rna_mols structs.
 * The strings @p bc and @p name are duplicated.
 * The read object @p r is duplicated.
 *
 * @param bam_data pointer to bam_data_t object to add RNA read to.
 * @param bc string of the barcode ID for the RNA read.
 * @param r Pointer to rna_read1_t object.
 * @param name name of the UMI or read used for deduplicating.
 *
 * @return 0 on success, -1 on error.
 */
int bam_data_rna_add_read(bam_data_t *bam_data, const char *bc, 
        const rna_read1_t *r, umishort umih);

/* De-duplicate RNA reads in bam_data_t object.
 * Loops through the fragments and collapses duplicate 
 * reads. See rna_read1_equal for definition of equality.
 *
 * @return 0 on success, -1 on error.
 */
int bam_data_rna_dedup(bam_data_t *bam_data);

/* Call variants in RNA bam_data_t object.
 * If no reads are present, do nothing.
 *
 * @param bam_data pointer to bam_data_t object.
 * @param gv pointer to g_var_t object.
 * @param cmap pointer to str_map object with chromosome name to indices.
 * @param min_qual minimum base phred quality in read to consider for 
 *  calling a variant allele.
 *
 * @return -1 on error, or the number of variants called.
 */
int bam_data_rna_var_call(bam_data_t *bam_data, g_var_t *gv, 
        str_map *cmap, uint8_t min_qual);

/*******************************************************************************
 * BAM ATAC
 ******************************************************************************/

/* Add ATAC read atac_read1_t to bam_data.
 * 
 * Adds the read @p r to the @f bc_data and @f atac_pairs structs.
 * The strings @p bc and @p qname are duplicated.
 * The read object @p r is duplicated.
 *
 * @param bam_data pointer to bam_data_t object to add ATAC read to.
 * @param bc string of the barcode ID for the ATAC read.
 * @param r Pointer to atac_read1_t object.
 * @param name name of the query read used to match pairs.
 *
 * @return 0 on success, -1 on error.
 */
int bam_data_atac_add_read(bam_data_t *bam_data, const char *bc, 
        const atac_read1_t *r, qshort qname);

/* De-duplicate the PCR duplicates in frags.
 * Loop through the frags and deduplicate the underlying reads.
 * If no barcodes or frags are present, do nothing.
 * See atac_rd_pair_equal for how equality is defined.
 *
 * @param bam_data pointer to bam_data_t object.
 *
 * @return 0 on success, -1 on error.
 */
int bam_data_atac_dedup(bam_data_t *bam_data);

/* Call variants in ATAC bam_data_t object.
 * If no fragments are present, do nothing.
 *
 * @param bam_data pointer to bam_data_t object.
 * @param gv pointer to g_var_t object.
 * @param cmap pointer to str_map object with chromosome name to indices.
 * @param min_qual minimum base phred quality in read to consider for 
 *  calling a variant allele.
 *
 * @return -1 on error, or the number of variants called.
 */
int bam_data_atac_var_call(bam_data_t *bam_data, g_var_t *gv, 
        str_map *cmap, uint8_t min_qual);

/* Call peaks in ATAC bam_data_t object.
 *
 * @param bam_data pointer to bam_data_t
 * @param pks Pointer to iregs_t object containing the peaks.
 * @param cmap pointer to str_map object that maps chromosome names
 *  to indices
 *
 * @return -1 on error, or the number of peaks called on success
 */
int bam_data_atac_peak_call(bam_data_t *bam_data, iregs_t *pks, 
        str_map *cmap);

int bam_data_fill_list(bam_data_t *bam_data);

/*******************************************************************************
 * bc_stats
 ******************************************************************************/

int bam_data_fill_stats(bam_data_t *bam_data, 
        const sam_hdr_t *rna_hdr, const sam_hdr_t *atac_hdr);

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
