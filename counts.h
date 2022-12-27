
#ifndef COUNTS_H
#define COUNTS_H

#include <stdlib.h>
#include "str_util.h"
#include "region.h"
#include "variants.h"
#include "htslib/hts.h"
#include "htslib/khash.h"

#define count_t uint16_t
#define n_bases 16 // 2^4

/*******************************************************************************
 * sequencing base
 ******************************************************************************/

/*! @typedef
 * @abstract sequencing base
 *
 * 
 * @field pos genomic position in g_pos struct object.
 * @field base A uint8_t of an observed nucleotide.
 *  This uses 4 bit encoding A:1, C:2, G:4, T:8, N:15.
 *  To go from 4-bit encoding to the 3-bit here, use seq_nt16_int function.
 *  To go from 4-bit encoding to char encoded nucleotide, use seq_nt16_str function.
 * @field qual The quality score of the observed base.
 * @field next Pointer to next pbc_t in a list
 */
typedef struct seq_base_t {
    g_pos pos;
    uint8_t base;
    uint8_t qual;
    struct seq_base_t *next;
} seq_base_t;

/*! @typedef
 * @abstract Store multiple pbc_t objects
 *
 * To get the next seq_base_t object, use p = p->next.
 *
 * @field head Pointer to the first seq_base_t object.
 * @field n The number of seq_base_t objects.
 */
typedef struct seq_blist_t {
    seq_base_t *head;
    size_t n;
} seq_blist_t;

/* Properly initialize a seq_base_t object.
 *
 * @param sbase Pointer o seq_base_t object to initialize.
 */
void init_seq_base(seq_base_t *sbase);
seq_base_t *alloc_seq_base();

/* Destroy a seq_base_t object
 *
 * @param b Pointer to seq_base_t object to destroy.
 * @return Pointer ot the next seq_base_t object in the list.
 */
seq_base_t *dstry_seq_base(seq_base_t *b);

/* Allocate and return a copy of a seq_base_t object.
 *
 * @param b Pointer to seq_base_t object.
 * @param ret Pointer to int to hold return status. 0 on success, -1 on error.
 * @return Pointer to copy of seq_base_t object.
 *  NULL on error.
 */
seq_base_t *copy_seq_base(const seq_base_t *b, int *ret);

/* Allocate and initialize a seq_blist_t object.
 *
 * @return Pointer to allocated seq_blist_t object.
 */
void init_seq_blist(seq_blist_t *s);
seq_blist_t *alloc_seq_blist();

/* Free the memory underlying a blist */
void free_seq_blist(seq_blist_t *s);

/* Check if two seq_blist objects are the same.
 *
 * @param l1 A seq_blist_t object.
 * @param l2 A seq_blist_t object.
 * @param qual_cmp Integer flag indicating whether to compare base quality
 *  in addition to observed base. To ignore qual, set to 0.
 * @return 1 if they are equal, 0 if not.
 */
int seq_blist_equal(seq_blist_t l1, seq_blist_t l2, int qual_cmp);

/* Add seq_base_t to a seq_blist_t object.
 *
 * First checks to see if position is present. If the position is present and 
 * if dup_ok != 0, there is no error. If dup_ok = 0, then 
 * will return -1 in error.
 *
 * If skip_dup = 1, then nothing is added if the position was already present. 
 * If skip_dup = 0, then the seq_base_t object is added anyways.
 *
 * If @p s or @p b are NULL, return -1 as error.
 *
 * The seq_base_t object is always copied if added.
 *
 * @param s Pointer to seq_blist_t to add @pb  to.
 * @param b Pointer to seq_base_t to add.
 * @param dup_ok Integer flag to indicate if duplicate position is OK. If 
 *  0, then return -1 on error if duplicate encountered.
 * @param skip_dup Integer flag to indicate whether to skip adding @p b when 
 * it is a duplicate position.
 *
 * @return -1 on error, 0 if the position was present and not added (when 
 *  skip_dup = 1), 1 if the position was present and added anyways, and 
 *  2 if the position was not present and added.
 * @note The seq_base_t object is always copied when added.
 */
int blist_add_base(seq_blist_t *s, const seq_base_t *b, int dup_ok, int skip_dup);

/* Copy a seq_blist_t object.
 *
 * @param s A pointer to a seq_blist_t object.
 * @param ret Pointer to int to hold return status. 0 on success, -1 on error.
 *
 * @return A pointer to the copied seq_blist_t object, or NULL on error.
 */
seq_blist_t *copy_seq_blist(const seq_blist_t *s, int *ret);

/* Set highest base quality in seq_blist
 *
 * Expects both arguments to be equal, possibly differing by base quality 
 * at the bases.
 *
 * @return 0 if success, -1 on error.
 */
int seq_blist_match_qual(seq_blist_t *bl, const seq_blist_t *cmp);

/*******************************************************************************
 * variant allele calls
 ******************************************************************************/

/*
 * 0 is reference allele
 * 1 is first alternate allele
 * 2 is second alternate allele
 * ...
 * 15 is N.
 */
#define MAX_ALLELE 16

/*! @typedef
 * @abstract variant allele call
 *
 * @field vix Variant integer ID. Map to chr name with contig_map.
 *  Missing values are encoded as 15 (0xf)
 * @field allele Index of variant allele. Ref is 0, first alt allele is 1, second is 2, ... 
 */
typedef struct vac_t {
    int32_t vix;
    uint8_t allele;
    uint8_t qual;
    struct vac_t *next;
} vac_t;

/* An unsorted list of vac_t objects.
 */
typedef struct vacs_t {
    vac_t *head;
    size_t n;
} vacs_t;

/* Initialize/allocate vac object.
 *
 * alloc returns dynamically allocated vac.
 * init sets the members to NA values.
 *
 * @return Pointer to vac object.
 */
void init_vac(vac_t *v);
vac_t *alloc_vac();
vac_t *dstry_vac(vac_t *v);

/* Initialize members of a vacs object.
 *
 * If @p vacs is null, do nothing.
 *
 * @param vacs Pointer to a vacs object.
 */
void init_vacs(vacs_t *vacs);

/* Free memory underlying a vacs_t object.
 *
 * Free the vac_t objects in a vacs_t object.
 * Does not free the vacs_t object itself.
 *
 * @param vacs Pointer to vacs_t object.
 */
void free_vacs(vacs_t *vacs);

/* Add a vac to vacs
 *
 * @param vacs Pointer to vacs object to add vac to.
 * @param v Pointer to a vac_t object to add
 *
 * @return 0 on success, -1 on error
 */
int vacs_add(vacs_t *vacs, vac_t *v);

/* Call variants from sequenced bases
 *
 * For each overlapping SNV, record the allele in a vac_t object 
 * as the integer of the allele (0 for ref, 1 for first alt, ...).
 * Then add the vac_t object to @p vacs.
 * If the vase in @p b is 'N', skip and return 0.
 * If the base quality is < min_qual, skip and return 0.
 */
int seq_base_call_var(seq_base_t *b, vacs_t *vacs, GenomeVar *gv, 
    contig_map *cmap, uint8_t min_qual);

int seq_blist_call_var(seq_blist_t *s, vacs_t *vacs, GenomeVar *gv, 
        contig_map *cmap, uint8_t min_qual);

/*******************************************************************************
 * sequencing gene
 ******************************************************************************/

/*! @typedef
 * @abstract Structure to store feature alignment
 *
 * @field gene_id Integer ID of the gene.
 * @field splice Integer code of splice.
 * @field next Pointer to next object in a list.
 */
typedef struct seq_gene_t {
   int32_t gene_id;
   uint8_t splice;
   struct seq_gene_t *next;
} seq_gene_t;

/* Initialize/allocate a seq_gene_t object.
 */
void seq_gene_init(seq_gene_t *g);
seq_gene_t *seq_gene_alloc();

seq_gene_t *seq_gene_dstry(seq_gene_t *g);

/* copy a seq_gene_t object.
 *
 * If NULL is passed, return NULL successfully.
 * The src object is deep copied, where a seq_gene_t object is allocated 
 * and the contents of src are copied to the allocated object.
 * The next member of the allocated copy is set to NULL.
 *
 * The return status is set in @p ret. The variable is modified, where 
 * 0 indicates success, -1 indiciates error (no memory).
 */
seq_gene_t *seq_gene_cpy(const seq_gene_t *src, int *ret);

typedef struct seq_glist_t {
    seq_gene_t *head;
    size_t n;
} seq_glist_t;

void seq_glist_init(seq_glist_t *gl);
seq_glist_t *seq_glist_alloc();

void seq_glist_free(seq_glist_t *gl);

seq_glist_t *seq_glist_cpy(const seq_glist_t *src, int *ret);

int seq_glist_add_gene(seq_glist_t *gl, const seq_gene_t *gene);

#endif // COUNTS_H
