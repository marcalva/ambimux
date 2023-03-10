
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

/* 
#define REF 0
#define ALT 1
#define OTHER 2
#define NA_ALLELE 3
*/
enum alleles {REF, ALT, OTHER, NA_ALLELE};
#define N_ALLELE 4

// #define MAX_BIN (((1<<18)-1)/7)

/************************************************************************
 * VCF file
 ***********************************************************************/

int load_vcf(const char *vcffn, const char *region, int region_set, 
        bcf_srs_t **sr, bcf_hdr_t **vcf_hdr);

int sub_vcf_samples(bcf_hdr_t **vcf_hdr, const char *samplefn);

/************************************************************************
 * Var
 ***********************************************************************/

typedef struct Var { 
    bcf1_t *b;
    int32_t vix;
    struct Var *next;
} Var;

// KHASH_INIT(var, var_id_t*, Var*, 1, _var_id_t_hash_func, _var_id_t_hash_equal);
KHASH_INIT(var, char*, Var*, 1, kh_str_hash_func, kh_str_hash_equal);

typedef struct chrmVar {
    Var *bins[MAX_BIN];
    uint16_t vars_n[MAX_BIN];
} ChrmVar;

typedef struct {
    str_map *chrm_ix; // chromosome index
    ChrmVar **chrms; // array to array of ChrmVar
    int chrms_m; // max number of elements
    Var **ix2var; // variant index to Var object
    int32_t n_v, n_e, n_a; // num. variants, num. elements, num. allocated
    bcf_hdr_t *vcf_hdr;
} GenomeVar;

// Functions

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

int get_var_bin(bcf1_t *b);

Var *init_var();

GenomeVar *init_genomevar();

/* initialize chrmVar object
 * return NULL if memory wasn't allocated
 */ 
ChrmVar *init_ChrmVar();

int add_var(GenomeVar *gv, bcf1_t *b, const bcf_hdr_t *hdr);

/*
 * @param max_miss ignore variatns with number missing alleles > max_miss.
 *   Set to negative value to ignore.
 * @param maf_cut ignore variants with maf <= maf_cut or maf >= (1-maf_cut).
 *   Set to negative value to ignore.
 */
GenomeVar *vcf2gv(bcf_srs_t *sr, bcf_hdr_t *vcf_hdr, int max_miss, double maf_cut);

void destroy_gv(GenomeVar *gv);

// destroy bcf1_t records to preserve memory
void free_bcf1_t(GenomeVar *gv);

/* Check if a SNP is bi-allelic.
 *
 * Must have 1 ref and 1 alt allele, and must be a SNP.
 * @return 0 if bi-allelic SNP, 1 otherwise.
 */
int is_biallelic_snp(bcf1_t *b);

/* get number of allele or samples missing
 *
 * @param vcf_hdr VCF header
 * @param b VCF line
 * @param nmiss updated number of missing alleles
 * @param na updated number of alternate allele counts
 * @param total updated total number of alleles
 * @return -1 on error, 0 on success.
 */
int bcf1_t_miss_maf(bcf_hdr_t *vcf_hdr, bcf1_t *b, int *nmiss, int *na, int *ntotal);

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
 * @param gv pointer to GenomeVar object.
 * @param vcf_hdr pointer to VCF header of the VCF lines in @p gv.
 * @param ids Array of the variant integer IDs to return dose for.
 *   The order of these variants will be preserved in @p dose.
 * @param ni length of @p ids array.
 * @param field one of "GT" or "GP".
 * @return float matrix, NULL on error.
 */
float **ap_array_gt(GenomeVar *gv, bcf_hdr_t *vcf_hdr, int32_t *ids, int ni, char *field);

/* Get the variants that overlap region [beg, end)
 *
 * If @p *vars is NULL, then the pointer pointed to by vars is replaced 
 * with a linked list of Var object. If @p vars already points to a linked 
 * list of Var objects, then the overlapping variants are appened to the 
 * list. There is no dummy head node.
 *
 * @param gv GenomeVar object to retrieve variants from.
 * @param ref Reference sequence name (chromosome) in character array.
 * @param beg 0-based start position of region.
 * @param end 0-based position the base pair after the end. Open end position.
 * @param vars pointer to Var pointer.
 *
 * @return Number of overlapping variants found.
 *  -1 if the reference is not found in GenomeVar. -2 on error
 *
 * @note The containers in vars must be freed, but not the actual 
 * contents.
 */
int region_vars(GenomeVar *gv, const char* ref, hts_pos_t beg, 
        hts_pos_t end, Var **vars);

/* Get variants that overlap given region [beg, end).
 * Use region_vars function above instead.
 *
 * If the chromosome @p ref is not present, do nothing and return 0.
 * Returns the overlapping in the array of the argument @p vars.
 * This function reallocates the memory if necessary.
 *
 * @param gv GenomeVar object to retrieve variants from.
 * @param ref Reference sequence name (chromosome) in character array.
 * @param beg 0-based start position of region.
 * @param end 0-based position the base pair after the end. Open end position.
 * @param vars pointer to array of Var pointers.
 * @param vars_m Pointer to integer that contains the length of the vars array.
 *
 * @return Number of overlapping variants in @p vars.
 *  -1 if the reference is not found in GenomeVar. -2 on error
 *
 * The @p vars array must be an allocated array, or NULL. The @p vars_len 
 * gives the length of the @p vars array. The function can increase the size of 
 * @p vars by calling realloc, and will change @p vars_len to reflect the udpated 
 * size.
 */
int vars_from_region(GenomeVar *gv, const char* ref, hts_pos_t beg, 
        hts_pos_t end, Var ***vars, int *vars_m);

int n_snp(GenomeVar *gv, int *n_snp);

/* return the variant at index ix
 *
 * @return NULL if ix is invalid, pointer to Var if successful.
 */
Var *gv_vari(GenomeVar *gv, int32_t ix);

#endif // VARIANTS_H

