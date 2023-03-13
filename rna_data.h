
#ifndef RNA_DATA_H
#define RNA_DATA_H

/* 
 * Data structures and methods for scRNA-seq
 */

#include <stdlib.h>
#include "str_util.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
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

/*! @typedef
 * @abstract Structure to store data of a single RNA read
 *
 * @field loc read region.
 * @field bases A seq_blist_t object to store observed bases.
 * @field genes A seq_glist_t object to store overlapping genes.
 * @field next pointer to next atac_read in a linked list, NULL otherwise.
 */
typedef struct rna_read1_t {
    g_region loc;
    ml_t(seq_base_l) bl;
    ml_t(seq_gene_l) gl;
} rna_read1_t;

/*! @typedef
 * @abstract Store duplicate RNA read
 *
 * @field dups Pointer to rna_read1_t array.
 * @field rds_per_dup Number of supporting reads for each read in dups.
 * @field size Number of reads in dups (length of valid elements in dups).
 * @field m Allocated size of dups and rds_per_dup.
 */
typedef struct rna_dups_t {
    rna_read1_t *dups;
    uint32_t *rds_per_dup;
    size_t size, m;
} rna_dups_t;

/*! @typedef
 * @abstract Store de-duplicated RNA reads representing molecules.
 *
 * @field dups Pointer to rna_read1_t array.
 * @field rds_per_dup Number of supporting reads for each read in dups.
 * @field size Number of reads in dups (length of valid elements in dups).
 * @field m Allocated size of dups and rds_per_dup.
 * @field n_reads Number of supporting reads
 */
typedef struct rna_mol_t {
    rna_dups_t *dups;
    g_region loc;
    ml_t(seq_base_l) bl;
    ml_t(seq_vac_l) vl;
    ml_t(seq_gene_l) gl;
    size_t n_reads;
    uint8_t _dedup;
} rna_mol_t;

// store reads by name
// KHASH_INIT(khrmn, char*, rna_mol_t *, 1, kh_str_hash_func, kh_str_hash_equal);

// store reads by UMI (hash)
KHASH_INIT(khrmn, umishort, rna_mol_t *, 1, kh_umi_hash_func, kh_umi_hash_equal);

#define mlar_lt(p, q) -1
ml_declare(mlar, rna_mol_t *, mlar_lt);

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

int rna_dups_init(rna_dups_t *rd);
rna_dups_t *rna_dups_alloc();

void rna_dups_free(rna_dups_t *rd);
void rna_dups_dstry(rna_dups_t *rd);

int rna_dups_add_read(rna_dups_t *rd, const rna_read1_t *r);

int rna_mol_init(rna_mol_t *rmol);
rna_mol_t *rna_mol_alloc();

void rna_mol_free(rna_mol_t *rmol);
void rna_mol_dstry(rna_mol_t *rmol);

rna_mol_t *rna_mol_cpy(const rna_mol_t *m, int *ret);

/* Add read to dups in rna_mol_t struct
 */
int rna_mol_add_read(rna_mol_t *mol, const rna_read1_t *r);

/* De-duplicate/get best RNA read and form rna_mol_t
 *
 * Expects non-null arguments.
 * If @p rd has no reads, return NULL successfully.
 */
int rna_mol_dedup(rna_mol_t *mol);

int rna_mol_var_call(rna_mol_t *m, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual);

/*
void count_umis(bam_rna_t *bam_r, int *n_umi, int *n_reads, int *n_bc);
void print_bam_rna(bam_rna_t *b, gene_anno_t *anno, bcf_hdr_t *vcf_hdr, 
        g_var_t *gv);
*/

#endif // RNA_DATA_H
