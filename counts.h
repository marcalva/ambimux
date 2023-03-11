
#ifndef COUNTS_H
#define COUNTS_H

#include <stdlib.h>
#include "str_util.h"
#include "region.h"
#include "variants.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "g_list.h"

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
} seq_base_t;

static inline int base_cmp(seq_base_t b1, seq_base_t b2){
    return( poscmp(b1.pos, b2.pos) );
}

static inline int seq_base_cmp(seq_base_t b1, seq_base_t b2, int cmp_qual){
    int pc = poscmp(b1.pos, b2.pos);
    if (pc != 0) return(pc);
    int bc = b1.base - b2.base;
    if (bc != 0) return(bc);
    if (cmp_qual){
        int qc = b1.qual - b2.qual;
        if (qc != 0) return(qc);
    }
    return(0);
}

// return 1 if equal, 0 if not
static inline int seq_base_equal(seq_base_t b1, seq_base_t b2, int cmp_qual){
    int qe = cmp_qual ? b1.qual == b2.qual : 1;
    if ( (poscmp(b1.pos, b2.pos) == 0) && 
         b1.base == b2.base && 
         qe ){
        return(1);
    } else {
        return(0);
    }
}

/* initialize sbase.
 * Undefined behaviour if sbase is null.
 */
void seq_base_init(seq_base_t *sbase);

/* Properly initialize a seq_base_t object.
 *
 * @param sbase Pointer o seq_base_t object to initialize.
 */
void seq_base_init(seq_base_t *sbase);
seq_base_t *seq_base_alloc();

/* Destroy a seq_base_t object
 *
 * @param b Pointer to seq_base_t object to destroy.
 * @return Pointer ot the next seq_base_t object in the list.
 */
void seq_base_dstry(seq_base_t *b);

/* Allocate and return a copy of a seq_base_t object.
 *
 * @param b Pointer to seq_base_t object.
 * @param ret Pointer to int to hold return status. 0 on success, -1 on error.
 * @return Pointer to copy of seq_base_t object.
 *  NULL on error.
 */
seq_base_t *seq_base_dup(const seq_base_t *b, int *ret);

/*******************************************************************************
 * base list
 ******************************************************************************/

/* declare base_list type
 * list type: ml_t(base_list)
 * node type: ml_node_t(base_list)
 *
 */
ml_declare(base_list, seq_base_t, base_cmp);

// base_list functions
#define base_list_init(ll) ml_init(base_list, ll)
#define base_list_free(ll) ml_free(base_list, ll)

/* insert seq_base_t into base list object.
 * @param ll list of bases of type ml_t(base_list)
 * @param base base of type seq_base_t
 * @param skip_dup if base already found and skip_dup=1, then do nothing
 * @param dup_ok if dup_ok=0 and base already found, return error.
 */
#define base_list_insert(ll, base, skip_dup, dup_ok) \
    ml_insert(base_list, ll, base, skip_dup, dup_ok)

#define base_list_cpy(dest, src) ml_cpy(base_list, dest, src)

int base_list_cmp(ml_t(base_list) bl1, ml_t(base_list) bl2, int cmp_qual);
int base_list_equal(ml_t(base_list) bl1, ml_t(base_list) bl2, int cmp_qual);
int base_list_match_qual(ml_t(base_list) *bl, const ml_t(base_list) *cmp);

/*******************************************************************************
 * base_listist
 ******************************************************************************/


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
 * @field vix Variant integer ID.
 *  Missing values are encoded as 15 (0xf)
 * @field allele Index of variant allele. Ref is 0, first alt allele is 1, second is 2, ... 
 */
typedef struct seq_vac_t {
    int32_t vix;
    uint8_t allele;
    uint8_t qual;
} seq_vac_t;

/* Initialize/allocate vac object.
 *
 * alloc returns dynamically allocated vac.
 * init sets the members to NA values.
 *
 * @return Pointer to vac object.
 */
void seq_vac_init(seq_vac_t *v);
seq_vac_t *seq_vac_alloc();
void seq_vac_dstry(seq_vac_t *v);

/*******************************************************************************
 * vac list
 ******************************************************************************/

#define vac_cmp(v1, v2) ((v1).vix - (v2).vix)

/* declare vac_list type
 * list type: ml_t(vac_list)
 * node type: ml_node_t(vac_list)
 *
 */
ml_declare(vac_list, seq_vac_t, vac_cmp);

// base_list functions
#define vac_list_init(ll) ml_init(vac_list, ll)
#define vac_list_free(ll) ml_free(vac_list, ll)

/* insert seq_base_t into base list object.
 * @param ll list of bases of type ml_t(base_list)
 * @param base base of type seq_base_t
 * @param skip_dup if base already found and skip_dup=1, then do nothing
 * @param dup_ok if dup_ok=0 and base already found, return error.
 */
#define vac_list_insert(ll, vac, skip_dup, dup_ok) \
    ml_insert(vac_list, ll, vac, skip_dup, dup_ok)

#define vac_list_cpy(dest, src) ml_cpy(vac_list, dest, src)


/* Call variants from sequenced bases
 *
 * For each overlapping SNV, record the allele in a seq_vac_t object 
 * as the integer of the allele (0 for ref, 1 for first alt, ...).
 * Then add the seq_vac_t object to @p vacs.
 * If the vase in @p b is 'N', skip and return 0.
 * If the base quality is < min_qual, skip and return 0.
 */

int seq_base_call_var(seq_base_t b, ml_t(vac_list) *vl, g_var_t *gv, 
        str_map *cmap, uint8_t min_qual);

int vac_list_call_var(ml_t(base_list) *bl, ml_t(vac_list) *vl, g_var_t *gv, 
        str_map *cmap, uint8_t min_qual);

// int base_list_equal(ml_t(base_list) bl1, ml_t(base_list) bl2, int cmp_qual);

/*******************************************************************************
 ******************************************************************************/

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
} seq_gene_t;

/* Initialize/allocate a seq_gene_t object.
 */
void seq_gene_init(seq_gene_t *g);
seq_gene_t *seq_gene_alloc();

void seq_gene_dstry(seq_gene_t *g);

/* compare two genes
 * return -1 if g1 < g2, 0 if g1 == g2, 1 if g1 > g2
 */
int seq_gene_cmp(seq_gene_t g1, seq_gene_t g2);

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
seq_gene_t *seq_gene_dup(const seq_gene_t *src, int *ret);

/*******************************************************************************
 * gene list
 ******************************************************************************/

#define gene_cmp(g1, g2) ((g1).gene_id - (g2).gene_id)

/* declare gene_list type
 * list type: ml_t(gene_list)
 * node type: ml_node_t(gene_list)
 *
 */
ml_declare(gene_list, seq_gene_t, gene_cmp);

// base_list functions
#define gene_list_init(ll) ml_init(gene_list, ll)
#define gene_list_free(ll) ml_free(gene_list, ll)

/* insert seq_base_t into base list object.
 * @param ll list of bases of type ml_t(base_list)
 * @param base base of type seq_base_t
 * @param skip_dup if base already found and skip_dup=1, then do nothing
 * @param dup_ok if dup_ok=0 and base already found, return error.
 */
#define gene_list_insert(ll, gene, skip_dup, dup_ok) \
    ml_insert(gene_list, ll, gene, skip_dup, dup_ok)

#define gene_list_cpy(dest, src) ml_cpy(gene_list, dest, src)

int gene_list_cmp(ml_t(gene_list) gl1, ml_t(gene_list) gl2);
int gene_list_equal(ml_t(gene_list) gl1, ml_t(gene_list) gl2);

// int base_list_equal(ml_t(base_list) bl1, ml_t(base_list) bl2, int cmp_qual);

/*******************************************************************************
 ******************************************************************************/

#endif // COUNTS_H
