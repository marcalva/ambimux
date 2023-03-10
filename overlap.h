
#ifndef OVERLAP_H
#define OVERLAP_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/kstring.h"
#include "htslib/sam.h"
#include "gtf_anno.h"
#include "variants.h"
#include "counts.h"

/* Overlap a bam1_t read with features given Annotation Find the features that
 * the read fully intersects, i.e. is fully contained.  For each feature, give
 * whether the read is spliced, unspliced, or ambiguous.
 *
 * The features and splice status are returned as integer indices in an array.
 * The number of elements in the array are stored in @p n_feat, while the
 * allocated size of the array is stored in @p m_feat.
 *
 * @param h pointer to header for bam record @p b.  @param b pointer to bam1_t
 * object for overlap.  @param a pointer to Annotation object with features.
 * @param n_feat pointer to integer giving the number of features that overlap.
 * @param m_feat pointer to integer giving the allocated size of the updated
 * arrays @p feat and @p splice.  @param feat pointer to array of char arrays
 * from @p a containing variant IDs.  @param splice integer array containing
 * splice status for each overlapping feature in @p feat.  @return 0 if
 * success, -1 if fail
 */
int overlap_bam1_feats(const sam_hdr_t *h, bam1_t *b, const Annotation *a,
        uint8_t *n_feat, uint8_t *m_feat, char ***feat, uint8_t **splice);

/*
 * @return Number of genes overlapping, or -1 on error.
 */
int bam1_feat_overlap(const sam_hdr_t *h, bam1_t *b, const Annotation *a, 
        seq_glist_t **gl);

/* Calculate whether a read is spliced, unspliced, or ambiguous.
 * 
 * @param b pointer to bam record.
 * @param gene pointer to gene.
 * @param ret set to 0 on success, -1 on error after calling.
 *
 * @return One of SPLICE UNSPLICE AMBIG
 */
uint8_t bam1_spliced(bam1_t *b, Gene *g, int *ret);

/* Find variants that overlap an alignment
 *
 * This function returns Var objects in a linked list pointer to by 
 * the argument in @p vars.The list @p vars is updated to contain the 
 * overlapping Var objects.
 *
 * Proper usage would be 
 *  Var *var_list;
 *  bam1_vars_overlap(sam_hdr, bam1, gv, &var_list);
 *
 * @param h pointer to header for bam record @p b.
 * @param b pointer to bam1_t object for overlap.
 * @param gv pointer to GenomeVar object with variants for overlap.
 * @param var pointer to address pointing to first Var object in a linked 
 *  list.
 * @return -1 on error, 0 on success.
 */
int bam1_vars_overlap(const sam_hdr_t *h, bam1_t *b, GenomeVar *gv, 
        Var **vars);

/* Get sequenced bases at overlapping variant sites in alignment.
 *
 */
int bam1_seq_base(const sam_hdr_t *h, bam1_t *b, GenomeVar *gv, seq_blist_t **bl);

/* Return the number of bases that overlap between [a1,a2) and [b1,b2) features. 
 * a1 must be less than a2, and b1 must be less than b2.
 * If a_strand and b_strand are opposite strands ('+' or '-'), return 0 overlap.
 * Otherwise, the strand is ignored.
 * the a and b regions are half open, with a2 and b2 non-inclusive
 *
 * @return Number of base pairs that overlap, or -1 on error.
 * */
int64_t bp_overlap(int64_t a1, int64_t a2, char a_strand, int64_t b1, int64_t b2, char b_strand);

#endif // OVERLAP_H
