
#ifndef SAM_READ_H
#define SAM_READ_H

#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/hts.h"
#include "str_util.h"

/* Get start and end 0-based coordinates of alignment.
 *
 * Provides 0-based [start,end) position of reference-consuming alignment.
 * The variables pointer to by start and end are updated after calling
 * If any of the arguments passed to the function are NULL, return -1 
 * without updating the pointers.
 *
 * @param b Pointer to bam1_t object.
 * @param tid Pointer to int32_t object to update chromosome ID.
 * @param start Pointer to int32_t object to update start coordinate.
 * @param end Pointer to int32_t object to update end coordinate.
 * @return 0 on success, -1 on error.
 */
int get_rcoord_bam(const bam1_t *b, int32_t *tid, int32_t *start, int32_t *end, 
        int adj_soft_clip);

/* Check if bam read is unmapped.
 *
 * @param b Pointer to bam1_t object.
 * 
 * @return 1 if unmapped, 0 otherwise.
 */
static inline int bam1_unmapped(const bam1_t *b){
    int f = 0, f2 = 0;
    if ((b->core.pos + 1) == (bam_endpos(b)))
        f = 1;
    if (((b->core.flag)&(BAM_FUNMAP)) > 0)
        f2 = 1;
    if (f != f2){
        err_msg(-1, 0, "Conflicting flag and position unmapped results"
                " %u-%u %i\n", b->core.pos + 1, bam_endpos(b), f2);
    }
    return f;
}

/* Load BAM file, header, and index
 *
 * The samFile, sam_hdr_t, and hts_idx_t argument objects are 
 * modified after successful call. They must be freed.
 *
 * @param bamfn BAM file name string.
 * @param bam Double pointer to samFile object.
 * @param bam_hdr Double pointer to sam_hdr_t object.
 * @param bam_idx Double pointer to hts_idx_t object.
 *
 * @return 0 on success, -1 on error.
 */
int load_bam(const char *bamfn, samFile **bam, sam_hdr_t **bam_hdr, 
        hts_idx_t **bam_idx);

// Overlap tid
int ovrlp_tid(sam_hdr_t *sam_hdr, bcf_hdr_t *bcf_hdr, int **t1, int **t2);

// Return the base that overlaps the BAM record. If no overlap, return 'N'
// base_ref is the tid chromosome that corresponds to the BAM record
int get_ovrlp_base(const bam1_t *b, int32_t base_ref, int32_t base_pos, 
        char *base, uint8_t *qual);

/* Get base pair of aligned bam query given reference position.
 *
 * The bam1_t is an aligned read sequence. The parameters ref and pos give 
 * the coordinates of the reference position to seek. The ref gives the int code of the 
 * reference chromosome, and pos gives the 0-based coordinate position. 
 * This function gets the overlapping base from the query sequence in b that 
 * aligns to the given reference position.
 * The passed arguments base and qual, which are the base pair and quality, 
 * are updated after the call if the read aligns and overlaps with the given position.
 *
 * @param b Pointer to bam1_t object.
 * @param ref Integer of the reference chromosome, given by the sam header.
 * @param pos Integer of the base pair position in the reference chromosome.
 * @param base Pointer to uint8_t to store the base pair. Encoded in 4-bit, where 
 *  A=1, C=2, G=4, T=8, N=15.
 * @param qual Pointer to uint8_t that stores the quality score.
 *
 * @return -1 on error, 0 if not found, 1 if found.
 */
int bam1_site_base(const bam1_t *b, int32_t ref, int32_t pos, 
        uint8_t *base, uint8_t *qual);

char *get_tag(const bam1_t *b, char tag[2]);

int64_t get_itag(const bam1_t *b, char tag[2]);

void print_bam1_t(const bam1_t *b);

#endif // SAM_READ_H
