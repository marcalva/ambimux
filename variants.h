
#ifndef VARIANTS_H
#define VARIANTS_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/vcf.h"
#include "str_util.h"
#include "bins.h"
#include "g_list.h"
#include "kvec.h"

/* 
#define REF 0
#define ALT 1
#define OTHER 2
#define NA_ALLELE 3
*/
enum alleles {REF, ALT, OTHER, NA_ALLELE};
#define N_ALLELE 4

/************************************************************************
 * var_t
 ***********************************************************************/

typedef struct var_t { 
    bcf1_t *b;
    int32_t vix;
} var_t;

/* var_t vector type
 * mv_t(vcfr_vec)
 */
mv_declare(vcfr_vec, var_t);

// TODO: see if I can keep this list sorted by position instead of index
/* var_t list type
 * ml_t(vcfr_list)
 */
// return -1 if v1 < v2 in position, first compare rid then pos
// return 0 if the same
static inline int var_pos_cmp(var_t v1, var_t v2){
    if (v1.b->rid < v2.b->rid) return(-1);
    else if (v1.b->rid > v2.b->rid) return(1);

    if (v1.b->pos < v2.b->pos) return(-1);
    else if (v1.b->pos > v2.b->pos) return(1);

    return(0);
}
// #define var_t_cmp(p, q) ( ((q).vix < (p).vix) - ((p).vix < (q).vix) )
ml_declare(vcfr_list, var_t, var_pos_cmp);

#define vcfr_list_init(ll) ml_init(vcfr_list, ll)
#define vcfr_list_free(ll) ml_free(vcfr_list, ll)
#define vcfr_list_insert(ll, vcfr, skip_dup, dup_ok) \
    ml_insert(vcfr_list, ll, vcfr, skip_dup, dup_ok)

// KHASH_INIT(var, var_id_t*, var_t*, 1, _var_id_t_hash_func, _var_id_t_hash_equal);
KHASH_INIT(var, char*, var_t*, 1, kh_str_hash_func, kh_str_hash_equal);

/************************************************************************
 * chr_bins_t
 ***********************************************************************/

typedef struct chr_bins_t {
    ml_t(vcfr_list) bins[MAX_BIN];
} chr_bins_t;

/* chr_bins_t vector type
 * mv_t(cbin_vec)
 */
mv_declare(cbin_vec, chr_bins_t *);

/************************************************************************
 * g_var_t
 ***********************************************************************/

typedef struct {
    str_map *chrm_ix; // chromosome  
    mv_t(cbin_vec) chr_bins;

    mv_t(vcfr_vec) vix2var; // variant index to var_t object.

    bcf_hdr_t *vcf_hdr;
} g_var_t;

// Functions

/************************************************************************
 * bcf functions
 ***********************************************************************/

int load_vcf(const char *vcffn, const char *region, int region_set, 
        bcf_srs_t **sr, bcf_hdr_t **vcf_hdr);

int sub_vcf_samples(bcf_hdr_t **vcf_hdr, const char *samplefn);

/* Create a contig map from a bcf header.
 *
 * Stores contigs from bcf hdr in a contig map. The resulting indexes 
 * match the rids of the bcf header.
 *
 * @param hdr Pointer to bcf_hdr_t object.
 * @param cm Pointer to cm object.
 *
 * @return 0 on success, -1 on error.
 */
int bcf_hdr_chr_ix(const bcf_hdr_t *hdr, str_map *cm);

/* Return a variant ID for a BCF line.
 *
 * Since the ID column may be missing, a variant can be identified using 
 * CHR, POS, ID, REF, and ALT. The variant ID is tab-delimited.
 *
 * This returns a pointer to char that contains the NULL terminated string. 
 * The function allocates the memory, and it is the caller's job to free 
 * the memory.
 *
 * @param h pointer to VCF header object.
 * @param b pointer to bcf1_t object.
 * @param delim delimiter between fields
 *
 * @return pointer to char that contains the NULL terminated variant ID string.
 * variant ID is of the format CHRM\tPOS\tID\tREF\tALT. If there are multiple 
 * ALT alleles, they are separated by commas.
 */
char *var_id(const bcf_hdr_t *h, bcf1_t *b, char delim);

/* get SNP allele
 *
 * return whether the given base is REF, ALT, or OTHER in the VCF/BCF record.
 * @param b bcf record with ref and alt allele to test against
 * @param base base allele to test
 * @return REF if base == b->ref, ALT if base == b->alt, OTHER if neither
 */
uint8_t base_ref_alt(bcf1_t *b, char base);

int get_var_len(bcf1_t *b);

int get_bcf1_bin(bcf1_t *b);

/* Check if a SNP is bi-allelic.
 *
 * Must have 1 ref and 1 alt allele, and must be a SNP.
 * @return 0 if bi-allelic SNP, 1 otherwise.
 */
int is_biallelic_snp(bcf1_t *b);

/* get number of allele or samples missing
 *
 * All pointers must be non-null, otherwise undefined behaviour.
 * Only bi-allelic variants are supported.
 * Uses the GT format field, error if missing.
 *
 * @param vcf_hdr VCF header
 * @param b VCF line
 * @param nmiss updated number of missing alleles
 * @param n_allele updated number of alternate allele counts
 * @param total updated total number of alleles
 * @return -1 on error, 0 on success.
 */
int bcf1_t_miss_maf(bcf_hdr_t *vcf_hdr, bcf1_t *b, int *nmiss, int *n_allele, 
        int *ntotal);

/* check if fmt tag is valid for getting allele probabilities from a VCF line
 * Assumes the GT or GP tag is present in the header.
 *
 * to be valid, the variant must
 *  have the fmt data from the GT or GP tag
 *  bi-allelic
 *  fmt must be present in the header
 *  fmt must be present in the VCF line
 *  is monoploid or diploid
 *  has at least one non-missing sample
 *
 * @return 
 *   -2 if not bi-allelic
 *   -1 if fmt is missing
 *   0 if valid 
 */
int is_gt_valid(bcf_hdr_t *vcf_hdr, bcf1_t *b);
int is_gp_valid(bcf_hdr_t *vcf_hdr, bcf1_t *b);

/************************************************************************
 * var_t functions
 ***********************************************************************/

var_t *var_alloc();

// free underlying memory, including destroying the bcf1 struct
void var_free(var_t *var);

/************************************************************************
 * chr_bins_t functions
 ***********************************************************************/

// Allocate and initialize, return NULL on error
chr_bins_t *chr_bins_alloc();

// free memory underlying chr_bins. This does not free the bcf records.
// call g_var_free_bcfr first.
void chr_bins_free(chr_bins_t *bins);

/************************************************************************
 * g_var_t functions
 ***********************************************************************/

g_var_t *g_var_alloc();

/* add a chromosome name to g_var_t
 * adds the string to the index, and allocates a chr_bins_t object.
 * @return integer index of chromosome name, or -1 on error
 */
int g_var_add_chr(g_var_t *gv, const char *chr);

/* add variant to g_var from a bcf record.
 * return 0 on success, -1 on error.
 */
int g_var_add_var(g_var_t *gv, bcf1_t *b, const bcf_hdr_t *hdr);

/* Reads bi-allelic snps from vcf file.
 *
 * @param max_miss ignore variants with number missing alleles > max_miss.
 *   Set to negative value to ignore.
 * @param maf_cut ignore variants with maf <= maf_cut or maf >= (1-maf_cut).
 *   Set to negative value to ignore.
 */
g_var_t *g_var_read_vcf(bcf_srs_t *sr, bcf_hdr_t *vcf_hdr, int max_miss, double maf_cut);

/*
 * free the underlying memory associated with the bcf records.
 */
void g_var_free_bcfr(g_var_t *gv);

/* free memory associated with gv, but not object itself.
 * frees the bcf records
 * does not free the vcf header, must be freed separately.
 */
void g_var_free(g_var_t *gv);

/*
 * bcf records should be freed with g_var_free_bcfr before calling this.
 */
void g_var_dstry(g_var_t *gv);

/************************************************************************
 ***********************************************************************/

/* get allele prob. for bi-allelic SNPs 
 *
 * This returns the allele prob. of the alternate allele for all samples 
 * in the bcf line. Only bi-allelic SNPs are valid.
 * If missing, the dose is equal to -1 for that sample.
 * If the samples were subsetted with bcf_hdr_set_samples, 
 * then only those samples will be returned.
 * Samples are returned in the same order as given in the header
 *
 * 
 * @param vcf_hdr VCF header
 * @param b VCF line
 * @param extra add this many elements to the allocated array
 * @param tag_id the numeric ID of the genotype tag (GT or GP)
 * @return double array of dose, or NULL on failure.
 * 
 * The returned array must be freed by caller.
 */
float *bcf1_ap_gt(bcf_hdr_t *vcf_hdr, bcf1_t *b, int extra);
float *bcf1_ap_gp(bcf_hdr_t *vcf_hdr, bcf1_t *b, int extra);

/* Return dose in a two-dimensional array.
 *
 * @param gv pointer to g_var_t object.
 * @param vcf_hdr pointer to VCF header of the VCF lines in @p gv.
 * @param ids Array of the variant integer IDs to return dose for.
 *   The order of these variants will be preserved in @p dose.
 * @param ni length of @p ids array.
 * @param field one of "GT" or "GP".
 * @return float matrix, NULL on error.
 */
float **ap_array_gt(g_var_t *gv, bcf_hdr_t *vcf_hdr, int32_t *ids, int ni, char *field);

/* Get the variants that overlap region [beg, end)
 *
 * Variants are added to the list ml_t(vcfr_list) in @p vars. 
 * Assumes that this points to a valid list, undefined behaviour 
 * otherwise.
 *
 * Returns an error if gv or ref is null.
 *
 * The region is specified as 0-based [beg, end).
 * Will return an error if end <= beg, or if beg or end are less than 0.
 *
 *
 * @param gv g_var_t object to retrieve variants from.
 * @param ref Reference sequence name (chromosome) in character array.
 * @param beg 0-based start position of region.
 * @param end 0-based position the base pair after the end. Open end position.
 * @param vars pointer to var_t pointer.
 *
 * @return Number of overlapping variants found.
 *  -1 if the reference is not found in g_var_t. -2 on error
 *
 * @note The containers in vars must be freed, but not the actual 
 * contents.
 */
int g_var_get_region_vars(g_var_t *gv, const char* ref, int32_t beg, 
        int32_t end, ml_t(vcfr_list) *vars);

int n_snp(g_var_t *gv, int *n_snp);

/* return the variant at index ix
 *
 * @return NULL if ix is invalid, pointer to var_t if successful.
 */
var_t *gv_vari(g_var_t *gv, int32_t ix);

#endif // VARIANTS_H

