
#ifndef ATAC_DATA_H
#define ATAC_DATA_H

/* 
 * Data structures and methods for scATAC-seq
 */

#include <stdlib.h>
#include "str_util.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "variants.h"
#include "counts.h"

typedef uint64_t qshort;

#define kh_qname_hash_func(i) (kh_int64_hash_func((i)))
#define kh_qname_hash_equal(a, b) ((a) == (b))

static inline qshort qname2qshort(const char *s){
    qshort qs = (qshort)(str_int64_hash(s));
    return(qs);
}

// collision check
// KHASH_INIT(khcol, qshort, uint8_t, 1, kh_qname_hash_func, kh_qname_hash_equal);

/*! @typedef
 * @abstract Single read
 *
 * @field qshort query name of read
 * @field loc read region.
 * @field bases A seq_bases_t object to store observed bases.
 * @field next pointer to next atac_read in a linked list, NULL otherwise.
 */
typedef struct atac_read1_t {
    g_region loc;
    seq_blist_t bases;
} atac_read1_t;

/*@ @typedef
 * @abstract Read pair
 *
 * @field qshort Query name of read pair.
 * @field loc Genomic region of atac read pair.
 * @field r1 Pointer to first atac read in pair.
 * @field r2 Pointer to secong atac read in pair.
 * @field s Number of reads (0, 1, or 2).
 */
typedef struct atac_rd_pair_t {
    atac_read1_t r1;
    atac_read1_t r2;
    uint8_t s;
} atac_rd_pair_t;

/* read pairs keyed by query name */
KHASH_INIT(khap, qshort, atac_rd_pair_t *, 1, kh_qname_hash_func, kh_qname_hash_equal);

typedef struct atac_dups_t {
    struct _atac_dup1_t {
        atac_rd_pair_t rd; // read pair
        size_t n_rd; // number of (duplicate) reads per pair.
    } *dups;
    size_t size, m;
} atac_dups_t;

/* duplicates keyed by region */
KHASH_INIT(khad, g_reg_pair, atac_dups_t *, 1, kh_reg_pair_hash, kh_reg_pair_equal);

/*! @typedef
 @abstract Fragment (deduplicated read pair)

 @field loc Genomic region.
 @field bases A seq_blist_t object to store observed bases.
 @field pks ATAC peaks the fragment overlaps.
 @field s Number of supporting reads for the fragment.

 @note
 */
typedef struct atac_frag_t {
    seq_blist_t bases;
    vacs_t vacs;
    iregn_t pks;
    size_t s;
} atac_frag_t;

KHASH_INIT(khaf, g_reg_pair, atac_frag_t *, 1, kh_reg_pair_hash, kh_reg_pair_equal);

/*! @typedef
 * @abstract Store a list of frags, duplicates, and reads for each barcode
 *
 * @field frags frags object. Accessed by region
 * @field dups Read duplicates. Accessed by region
 * @field reads Hanging reads from the pairs. Accessed by query name
 */
typedef struct bc_atac_t {
    khash_t(khaf) *frags;
    khash_t(khad) *dups;
    khash_t(khap) *pairs;
} bc_atac_t;

// hash table of key: barcode string val: bc_atac object
KHASH_INIT(khab, char *, bc_atac_t *, 1, kh_str_hash_func, kh_str_hash_equal);

/*! @typedef
 * @abstract Store reads and fragments for each barcode.
 *

 * @field bc_dat hash table of bc_atac objects keyed by barcode (char*)
 * @field barcodes str_map of barcodes
 */
typedef struct bam_atac_t {
    khash_t(khab) *bc_dat;
} bam_atac_t;

/*******************************************************************************
 * atac_read1_t
 ******************************************************************************/

/* Set members of atac read to 0.
 *
 * If NULL is passed instead of pointer, nothing happens.
 * Sets the query name to NULL, initializes the seq blist object, and 
 * sets next to NULL.
 *
 * @param ar Pointer to atac read object.
 */
void atac_read_set0(atac_read1_t *ar);

/* Initialize atac read
 *
 * Allocates and zeroes the members.
 *
 * @return Pointer to dynamically allocated object.
 * @note The object must be freed by caller.
 */
atac_read1_t *atac_read_init();

/* Free an atac_read object's members.
 * This does not destroy the object pointed to by @p ar, but frees 
 * the memory of its members.
 */
void atac_read_free(atac_read1_t *ar);

/* Free an atac_read object and its members.
 * The object pointed to by @p ar is freed.
 * 
 * If NULL is passes, nothing happens.
 */
void atac_read_dstry(atac_read1_t *ar);

/* Create a deep copy of an atac read object.
 *
 * The object returned must be freed by the caller.
 *
 * Note that the next field of the copy is set to NULL to avoid any 
 * mistakes.
 *
 * If @p r is NULL, returns NULL.
 *
 * The int pointed to by ret is modified after the function call. If an 
 * error occured, it is modified to -1. Otherwise, it is set to 0.
 *
 * @param r Pointer to atac read object.
 * @param ret Pointer to int to store return status.
 */
atac_read1_t *atac_read_cpy(const atac_read1_t *r, int *ret);

/* Compare if two atac reads are equal.
 *
 * Two atac reads are equal if they have the same genomic position, the same 
 * number of pos bases, and the same observed bases at these positions. The 
 * base quality is not considered in the comparison.
 *
 * If both arguments are NULL, 1 for equality is returned.
 * If only one argument is NULL, 0 for inequality is returned.
 * 
 * @param r1 Pointer to atac_read1_t.
 * @param r2 Pointer to atac_read1_t.
 * @return 1 if the reads are equal, 0 if they are not equal
 */
int atac_read_equal(atac_read1_t r1, atac_read1_t r2);

/* Add base to atac_read1_t object.
 *
 * The base argument is copied and added to the seq_blist_t bases field in @p r.
 *
 * @return 0 on success, -1 on error.
 */
int atac_read1_add_base(atac_read1_t *r, const seq_base_t *base);

/*******************************************************************************
 * atac_rd_pair_t
 ******************************************************************************/

/* Set 0 for read pair
 */
void atac_rd_pair_set0(atac_rd_pair_t *rp);

/* Initialize an atac_rd_pair_t object.
 *
 * @return Pointer to allocated object, NULL on error.
 */
atac_rd_pair_t *atac_rd_pair_init();

/* Destroy an atac_rd_pair.
 *
 * Free all the underlying memory and the object itself.
 */
void atac_rd_pair_free(atac_rd_pair_t *rp);
void atac_rd_pair_dstry(atac_rd_pair_t *rp);

/* Copy an atac_rd_pair_t object.
 *
 * If rp is NULL, return NULL.
 * 
 * Create a deep copy of the atac_rd_pair_t object pointed to by rp.
 *
 * The object must be freed by the caller.
 *
 * @return NULL, or newly allocated object.
 */
atac_rd_pair_t *atac_rd_pair_cpy(const atac_rd_pair_t *rp, int *ret);

/* Add an atac_read1_t to read pair.
 *
 * If no reads are present in the pair, add to @f r1. If r1 is 
 * present add the read such that r1 <= r2 relative to genomic region.
 *
 * Expects non-null arguments.
 * Expects the query name of @p rp is set and equal to that of @p ar.
 *
 * @param rp Pointer to atac_rd_pair_t object to add to.
 * @param ar Pointer to atac_read1_t object to add.
 * @return 0 on success, -1 on error.
 */
int atac_rd_pair_add_read(atac_rd_pair_t *rp, const atac_read1_t *ar);

/* Test if two atac read pairs are equal.
 *
 * Expects both arguments are not null.
 *
 * If the number of reads are not the same, 0 for inequality is returned.
 * Read 1 of @p rp1 and @p rp2 are tested for equality, then read 2 
 * of each argument are tested. If either are inequal, the 0 for 
 * inequality is returned.
 *
 * If the number of reads and each read is equal, then 1 is returned for 
 * equality.
 *
 * @return 1 for equality, 0 for inequality, -1  for error.
 */
int atac_rd_pair_equal(const atac_rd_pair_t *rp1, const atac_rd_pair_t *rp2);

/* Set the highest base quality in a read of the two reads 
 *
 * Expects @p rp and @p cmp have the same bases.
 * Compares the quality score at each base between the two reads. If the 
 * base quality is higher in @p cmp, then set the quality in @p rp 
 * to that of @p cmp.
 *
 * @param rp Pointer to atac_rd_pair_t to modify the qualities of.
 * @param cmp Pointer to atac_rd_pair_t to compare the qualities.
 * @return 0 on success, -1 on error.
 */
int atac_rd_pair_match_qual(atac_rd_pair_t *rp, const atac_rd_pair_t *cmp);

/*******************************************************************************
 * atac_dups_t
 ******************************************************************************/

/* Allocate and initialize and atac dups object.
 *
 * Memory is allocated and the members are initialized.
 * @return Allocated and initialized object, or NULL on error.
 */
atac_dups_t *atac_dups_init();

/* Free the memory an atac dups object and its memers.
 *
 * All read pairs are freed.
 */
void atac_dups_free(atac_dups_t *d);
void atac_dups_dstry(atac_dups_t *d);

/* Add a read pair to an atac dups object.
 *
 * The read pair is added to atac dups. If the read pair was 
 * found, increment @f rds_per_dup in @p d. If the read pair was not found, 
 * then copy the read pair, add to @f dups, and set the rds_per_dup to 1.
 * Sets the region from @p rp.
 *
 * Expects non-null arguments.
 * w
 * @param d A pointer to the atac_dups_t object to add the read pair to.
 * @param rp A pointer to the read pair object to add.
 * @return 0 on success, -1 on error or if memory allocation fails.
 */
int atac_dups_add_pair(atac_dups_t *d, const atac_rd_pair_t *rp);

/*******************************************************************************
 * atac_frag_t
 ******************************************************************************/

atac_frag_t *atac_frag_init();
void atac_frag_dstry(atac_frag_t *f);

/* Deduplicate an atac dups object and form a frag.
 *
 * Expects non-null @p d.
 * Expects at least 1 PCR duplicate in @p d with at least one supporting read.
 *
 * Forms a fragment from a PCR duplicate with the most supporting reads.
 * If two distinct PCR duplicates have the same number of supporting reads, 
 * return NULL to discard the PCR duplicate.
 *
 * Returns a fragment from the best PCR duplicate with the most supporting 
 * reads. The (paired) region is copied, and the bases are added to the frag 
 * from the best PCR duplicate read pair.
 */
atac_frag_t *atac_dups_dedup(const atac_dups_t *d, int *ret);

/* call variants from atac fragments
 *
 * Calls the seq_blist_call_var function.
 *
 * Expects non-null arguments
 * @param f Pointer to atac_frag_t object.
 * @param gv Pointer to GenomeVar object.
 * @param cmap Contig map between integer IDs and chromosome names.
 */
int atac_frag_var_call(atac_frag_t *f, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual);

/* Call peaks for a fragment.
 *
 * @param f The fragment to get overlapping peaks for.
 * @param reg The g_reg_pair of the fragment.
 * @param pks Pointer to the iregs_t which contains the peaks.
 * @param cmap Pointer to the contig map.
 * @return The number of overlapping peaks added, or -1 on error.
 */
int atac_frag_peak_call(atac_frag_t *f, g_reg_pair reg, iregs_t *pks, 
        contig_map *cmap);

/*******************************************************************************
 * bc_atac_t
 ******************************************************************************/

bc_atac_t *bc_atac_init();
void bc_atac_dstry_pairs(bc_atac_t *bca);
void bc_atac_dstry_dups(bc_atac_t *bca);
void bc_atac_dstry_frags(bc_atac_t *bca);
void bc_atac_dstry(bc_atac_t *bca);

/* Add read to bc atac object.
 *
 * Adds the read to the bca pairs field.
 *
 * Both arguments must not be NULL.
 * The query name of @p ar must be set.
 * The pairs field of @p bca must be set.
 *
 * @param bca Pointer to bc_atac_t object to add to.
 * @param ar Pointer to atac_read1_t object to add.
 * @return 0 on success, -1 on error.
 */
int bc_atac_add_read(bc_atac_t *bca, const atac_read1_t *ar, qshort qname);

/* Add a read pair to dups object in bc_atac
 *
 * @return 0 on success, -1 on error.
 */
int bc_atac_add_dup(bc_atac_t *bca, atac_rd_pair_t *rp);

/* Form duplicates from atac reads
 *
 * Expects non-null @p bca argument.
 * The read pairs in the bca cannot be NULL.
 * Only fully formed read pairs with two reads are added.
 * Duplicates are reads with equal read pair locations.
 * Note with PCR duplicates, there may be differing bases.
 * Assumes the read pairs are ordered.
 *
 * @param bca Pointer to bc_atac_t object to form duplicates for.
 * @return 0 on success, -1 on error.
 */
int bc_atac_form_dups(bc_atac_t *bca);

/* De-duplicate the PCR duplicates and form frags.
 *
 * Expects non-null @p bca argument and @f dups is non-null.
 * Calls atac_dups_dedup on the @f dups field and fills the 
 * @f frags field.
 * Assumes the read pairs are ordered.
 *
 * @param bca Pointer to bc_atac_t object to form duplicates for.
 * @return 0 on success, -1 on error.
 */
int bc_atac_dedup(bc_atac_t *bca);

/* Call variants from frags in a bc_atac object.
 *
 * Loop through frags and call variants from the observed bases.
 *
 * @param bca Pointer to bc_atac_t object.
 * @param gv Pointer to GenomeVar object.
 * @param cmap Pointer to contig_map object.
 * @return The number of variants called, or -1 on error.
 */
int bc_atac_var_call(bc_atac_t *bca, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual);

/* Call peaks in a barcode
 *
 * @return The number of overlapping peaks added, or -1 on error.
 */
int bc_atac_peak_call(bc_atac_t *bca, iregs_t *pks, contig_map *cmap);

/*******************************************************************************
 * bam_atac_t
 ******************************************************************************/

bam_atac_t *bam_atac_init();
void bam_atac_free_pairs(bam_atac_t *bam_a);
void bam_atac_free_dups(bam_atac_t *bam_a);
void bam_atac_dstry(bam_atac_t *bam_a);
int bam_atac_add_read(bam_atac_t *bam_a, const char *bc, const atac_read1_t *ar, 
        qshort qname);
int bam_atac_form_dups(bam_atac_t *bam_a);
int bam_atac_dedup(bam_atac_t *bam_a);

/* Call variants in a bam atac object.
 *
 * Loop through the barcodes and call variants in each fragment.
 * Expects non-null arguments.
 * The frags in the bam atac object must be formed from calling 
 * bam_atac_dedup.
 *
 * @param bama_a Pointer to bam_atac_t object.
 * @param gv Pointer to GenomeVar object.
 * @param cmap Pointer to contig_map object.
 * 
 * @return The number of variants called, or -1 on error.
 */
int bam_atac_var_call(bam_atac_t *bam_a, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual);

/* Call peaks in a bam object.
 *
 * @return The number of peaks called, or -1 on error.
 */
int bam_atac_peak_call(bam_atac_t *bam_a, iregs_t *pks, contig_map *cmap);

/*******************************************************************************
 * miscellaneous
 ******************************************************************************/

/* Count number of fragments in bam_atac
 */
void count_frags(bam_atac_t *bam_a, int *n_frag, int *n_reads, int *n_bc);

uint64_t get_n_buckets(bam_atac_t *bam_a);

void print_atac_read(atac_read1_t *ar);

void print_vac_bam(bam_atac_t *b, GenomeVar *gv);

void print_frag_dup(bam_atac_t *b);

void print_dups_n(bam_atac_t *b);

#endif // ATAC_DATA_H
