
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
    seq_blist_t bases;
    seq_glist_t genes;
} rna_read1_t;

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
int rna_read1_equal(const rna_read1_t *r1, const rna_read1_t *r2, int cmp_qual);

/* Add gene to rna_read1_t object.
 *
 * The gene argument is copied and added to the seq_glist_t genes field in @p r.
 *
 * @return 0 on success, -1 on error.
 */
int rna_read1_add_gene(rna_read1_t *r, const seq_gene_t *gene);
int rna_read1_add_base(rna_read1_t *r, const seq_base_t *base);

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
    size_t *rds_per_dup;
    size_t size, m;
} rna_dups_t;

int rna_dups_init(rna_dups_t *rd);
rna_dups_t *rna_dups_alloc();

void rna_dups_free(rna_dups_t *rd);
void rna_dups_dstry(rna_dups_t *rd);

int rna_dups_add_read(rna_dups_t *rd, const rna_read1_t *r);

// store reads by name
KHASH_INIT(khrdn, char *, rna_dups_t *, 1, kh_str_hash_func, kh_str_hash_equal);

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
    g_region loc;
    seq_blist_t bases;
    seq_glist_t genes;
    vacs_t vacs;
    size_t n_reads;
} rna_mol_t;

void rna_mol_init(rna_mol_t *rmol);
rna_mol_t *rna_mol_alloc();

void rna_mol_free(rna_mol_t *rmol);
void rna_mol_dstry(rna_mol_t *rmol);

rna_mol_t *rna_mol_cpy(const rna_mol_t *m, int *ret);

/* De-duplicate/get best RNA read and form rna_mol_t
 *
 * Expects non-null arguments.
 * If @p rd has no reads, return NULL successfully.
 */
rna_mol_t *rna_dups_dedup(rna_dups_t *rd, int *ret);

int rna_mol_var_call(rna_mol_t *m, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual);

// store reads by name
KHASH_INIT(khrmn, char*, rna_mol_t *, 1, kh_str_hash_func, kh_str_hash_equal);

/*! @typedef
 * @abstract Store barcode RNA data
 *
 * @field dups Pointer to hash table containing duplicates, such as reads 
 *  from the same UMI
 * @field mols Pointer to hash table containing de-duplicated RNA molecules.
 */
typedef struct bc_rna_t {
    khash_t(khrdn) *dups;
    khash_t(khrmn) *mols;
} bc_rna_t;

void bc_rna_init(bc_rna_t *br);
bc_rna_t *bc_rna_alloc();

void bc_rna_free_dups(bc_rna_t *br);
void bc_rna_free_mols(bc_rna_t *br);
void bc_rna_free(bc_rna_t *br);
void bc_rna_dstry(bc_rna_t *br);

/* Add read to bc_rna_t object.
 * 
 * Adds the read @p r to the dups field in @p br.
 * A copy is created and added.
 *
 * @return 0 on success, -1 on error.
 */
int bc_rna_add_read(bc_rna_t *br, const rna_read1_t *r, const char *name);

/* Add molecule to bc_rna_t object.
 *
 * Adds a copy of @p m to @p br. This copies @p name and places as the 
 * key for the allocated @p m.
 *
 * @param name UMI name of the molecule.
 * @return 0 on success, -1 on error.
 */
int bc_rna_add_mol(bc_rna_t *br, const rna_mol_t *m, const char *name);

/* De-duplicate RNA reads
 *
 * De-duplicates the reads in dups to form molecules, and adds them to 
 * the mols field.
 *
 * @return 0 on success, -1 on error.
 */
int bc_rna_dedup(bc_rna_t *br);

/* Call variants at sequenced bases
 *
 * @return Number of variants called, or -1 on error.
 */
int bc_rna_var_call(bc_rna_t *br, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual);

KHASH_INIT(khrbc, char *, bc_rna_t *, 1, kh_str_hash_func, kh_str_hash_equal);

/*! @typedef
 * @abstract Store bam RNA data
 *
 * @field bc_rna Pointer to bc_rna_t
 */
typedef struct bam_rna_t {
    khash_t(khrbc) *bc_rna;
} bam_rna_t;

void bam_rna_init(bam_rna_t *bam_r);
bam_rna_t *bam_rna_alloc();

void bam_rna_free_dups(bam_rna_t *bam_r);
void bam_rna_free(bam_rna_t *bam_r);
void bam_rna_dstry(bam_rna_t *bam_r);

int bam_rna_add_read(bam_rna_t *bam_r, const char *bc, const rna_read1_t *r, 
        const char *name);

int bam_rna_dedup(bam_rna_t *bam_r);

/* Call variants in bam rna file
 *
 * @return Number of variants called, or -1 on error.
 */
int bam_rna_var_call(bam_rna_t *bam_r, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual);

void count_umis(bam_rna_t *bam_r, int *n_umi, int *n_reads, int *n_bc);
void print_bam_rna(bam_rna_t *b, Annotation *anno, bcf_hdr_t *vcf_hdr, 
        GenomeVar *gv);

#endif // RNA_DATA_H
