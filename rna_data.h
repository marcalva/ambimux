
#ifndef RNA_DATA_H
#define RNA_DATA_H

/* 
 * Data structures and methods for scRNA-seq
 */

#include <stdlib.h>
#include "str_util.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "kavl.h"
#include "variants.h"
#include "counts.h"
#include "gtf_anno.h"

// store UMI as hashed 64 bit integer
typedef uint64_t umishort;

#define kh_umi_hash_func(i) (kh_int64_hash_func((i)))
#define kh_umi_hash_equal(a, b) ((a) == (b))

static inline umishort umi2umishort(const char *s){
    umishort us = (umishort)(str_int64_hash(s));
    return(us);
}

/* @brief Structure to store data of a single RNA read
 */
typedef struct rna_read1_t {
    g_region loc; // location of the rna read
    ml_t(seq_base_l) bl; // observed bases
    ml_t(seq_gene_l) gl; // observed genes
    umishort umi;
    uint32_t n_reads; // number of supporting reads
} rna_read1_t;

/* Compare if two rna reads are equal.
 *
 * Compare rna reads @p r1 and @p r2.
 * Expects non-null arguments.
 * If both names are set, the strings are compared for equality.
 * Two reads are equal if 
 *  1) the same names (if they are both non-null, ignored otherwise)
 *  2) the same region,
 *  3) the same number of seq bases
 *  4) the same number of genes
 *  5) the same base positions
 *  6) each position has the same base.
 *  7) the same gene IDs
 *  8) the same splice for each gene ID
 *  9) the same qualities at each base (if @p cmp_qual != 0).
 *
 * @return 1 if equal, 0 if not equal, -1 on error.
 */
int rna_read1_cmp(const rna_read1_t *r1, const rna_read1_t *r2, int cmp_qual);
int rna_read1_equal(const rna_read1_t *r1, const rna_read1_t *r2, int cmp_qual);

// tree for rna_read1_t.
// Key is `rna_read1_t`, val is `uint32_t` for number of supporting seq reads.
// Sorted by key using `rna_read1_cmp`.
#define rna_read1_keycmp(key1, key2) (rna_read1_cmp(&(key1), &(key2), 0))
KBTREE_INIT(kb_rr, rna_read1_t, rna_read1_keycmp);

#define kh_rna_read1_hf(r) (kh_umi_hash_func((r).umi))
#define kh_rna_read1_he(r1, r2) ((r1).umi == (r2).umi)
KHASH_INIT(kh_rr, rna_read1_t, char, 0, kh_rna_read1_hf, kh_rna_read1_he);

mv_declare(mv_rr, rna_read1_t);

/*******************************************************************************
 * rna dups
 ******************************************************************************/

/* @brief Store duplicate RNA read
 *
 * Designed to store RNA reads (which may differ) with the same UMI.
 * The arrays dups and rds_per_dup have the same size and allocation.
 */
typedef struct rna_dups_t {
    umishort umi;
    mv_t(mv_rr) rds;
    // kbtree_t(kb_rr) *bt_rds;
} rna_dups_t;

#define kh_rna_dups_hf(r) (kh_umi_hash_func((r).umi))
#define kh_rna_dups_he(r1, r2) ((r1).umi == (r2).umi)
KHASH_INIT(kh_rd, rna_dups_t, char, 0, kh_rna_dups_hf, kh_rna_dups_he);

/*******************************************************************************
 * rna mol
 ******************************************************************************/

/* @brief Store de-duplicated RNA reads representing molecules.
 *
 */
typedef struct rna_mol_t {
    umishort umi;
    g_region loc; // location of the rna molecule
    ml_t(seq_base_l) bl; // observed bases
    ml_t(seq_vac_l) vl; // observed variant calls
    ml_t(seq_gene_l) gl; // observed genes
    uint32_t n_reads; // number of supporting reads
} rna_mol_t;

#define rna_mol_cmp(p, q) ( ( ((q).umi) < ((p).umi) ) - ( ((p).umi) < ((q).umi) ) )
// always insert at the beginning since we don't need ordering
//  makes insertions much much faster
#define always_lt(p, q) -1
ml_declare(ml_rm, rna_mol_t, always_lt)

/*******************************************************************************
 * rna read1
 ******************************************************************************/

void rna_read1_init(rna_read1_t *r);
rna_read1_t *rna_read1_alloc();

void rna_read1_free(rna_read1_t *r);
void rna_read1_dstry(rna_read1_t *r);

/* Copy rna_read1_t object.
 *
 * If @p r is null, return null successfully.
 * Allocate a new rna_read1 object and copy the contents.
 * The return status is held by @p ret. 0 on success, -1 on error.
 *
 * @return Pointer to rna_read1_t copy.
 */
rna_read1_t *rna_read1_cpy(const rna_read1_t *r, int *ret);

/* Add gene to rna_read1_t object.
 *
 * The gene argument is copied and added to the seq_glist_t genes field in @p r.
 *
 * @return 0 on success, -1 on error.
 */
int rna_read1_add_gene(rna_read1_t *r, const seq_gene_t *gene);
int rna_read1_add_base(rna_read1_t *r, seq_base_t base);

/* Set highest base quality in an rna read
 *
 * Compare base qualities between @p r and @p cmp and set the 
 * highest quality of the two in @p r.
 *
 * Expects r and cmp to be equal (have the same bases).
 *
 * @return 0 on success, -1 on error.
 */
int rna_read1_match_qual(rna_read1_t *r, const rna_read1_t *cmp);

/*******************************************************************************
 * rna dups
 ******************************************************************************/

int rna_dups_init(rna_dups_t *rd);
rna_dups_t *rna_dups_alloc();

void rna_dups_free(rna_dups_t *rd);
void rna_dups_dstry(rna_dups_t *rd);

/**
 * @brief Add a read to the rna_dups_t structure.
 *
 * This function adds a read to the rna_dups_t structure. If the read already
 * exists in the structure, the base qualities of the existing read are updated
 * to the highest values between the existing and new reads, and the count of
 * duplicates for that read is incremented. If the read does not exist, it is
 * added to the structure with a count of 1.
 *
 * If the input arguments are null, return -1 with an error message
 * If the dups or rds_per_dup arrays in the rd structure are
 * it returns -1 with an an error message.
 *
 * The function then iterates over the existing reads in the rd structure using
 * the rna_read1_cmp function to compare the input read with each existing read.
 * If rna_read1_cmp returns a negative value, indicating an error, the function
 * returns an error message. If a match is found (rna_read1_cmp returns 0), the
 * rna_read1_match_qual function is called to update the base qualities of the
 * existing read with the highest values between the existing and new reads. The
 * count of duplicates for that read is then incremented, and the function
 * returns 0 to indicate success.
 *
 * If no match is found, the function reallocates memory for the dups and
 * rds_per_dup arrays if necessary, using the realloc function. If the
 * reallocation fails, it returns an error message.
 *
 * The function then creates a copy of the input read using the rna_read1_cpy
 * function and adds it to the dups array. The count of duplicates for the new
 * read is set to 1 in the rds_per_dup array. The size of the rd structure is
 * incremented, and the function returns 0 to indicate success.
 *
 * If any error occurs during the execution of the function, it returns -1 to
 * indicate failure.
 *
 * @param rd Pointer to the rna_dups_t structure.
 * @param r Pointer to the rna_read1_t structure representing the read to add.
 *
 * @return 0 on success, -1 on error.
 */
int rna_dups_add_read(rna_dups_t *rd, const rna_read1_t *r);

/**
 * @brief Deduplicate RNA reads and return the best read based on the number of supporting reads.
 *
 * This function analyzes a set of duplicate RNA reads and selects the one with the highest 
 * number of supporting reads. If there is an ambiguity or error, it returns NULL and sets 
 * the return code appropriately.
 *
 * @param dups A pointer to the RNA duplicates structure (`rna_dups_t`) containing reads.
 * @param ret  A pointer to an integer where the function will store the return code:
 *             - `0` on success
 *             - non-zero on error
 * 
 * @return A pointer to an allocated `rna_mol_t` structure containing the best RNA read,
 *         or NULL in case of errors or ambiguity.
 */
rna_mol_t *rna_dups_dedup(rna_dups_t *dups, int *ret);

/*******************************************************************************
 * rna_mol_t
 ******************************************************************************/

int rna_mol_init(rna_mol_t *rmol);
rna_mol_t *rna_mol_alloc();

void rna_mol_free(rna_mol_t *rmol);
void rna_mol_dstry(rna_mol_t *rmol);

rna_mol_t *rna_mol_cpy(const rna_mol_t *m, int *ret);

/**
 * Call genetic variants in an RNA molecule.
 *
 * This function examines the base list of an RNA molecule to call genetic variants. 
 * It uses a provided variant calling algorithm and then frees the base list of the 
 * RNA molecule.
 *
 * @param m A pointer to the RNA molecule (`rna_mol_t`) to analyze.
 * @param gv A pointer to genetic variant information (`g_var_t`).
 * @param cmap A pointer to a chromosome map (`str_map`) used in the variant calling process.
 * @param min_qual The minimum quality threshold for calling a variant.
 *
 * @return An integer status code:
 *         - `0` on success
 *         - non-zero error code on failure
 */
int rna_mol_var_call(rna_mol_t *m, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual);

#endif // RNA_DATA_H
