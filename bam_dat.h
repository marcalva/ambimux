
#ifndef BAM_DAT_H
#define BAM_DAT_H

#include <stdlib.h>
#include "rna_data.h"
#include "atac_data.h"

/*! @typedef Structure to hold the atac and rna pileup
 *
 * @field rna A pointer to a bam_rna_t object.
 * @field atac A pointer to a bam_atac_t object.
 * @field bcs A pointer to a str_map containing the barcodes.
 *  This list contains all valid barcodes, observed in the data 
 *  and given in a whitelist.
 *  The order is preserved, such that any output will follow this order.
 */
typedef struct {
    bam_rna_t *rna;
    bam_atac_t *atac;
    str_map *bcs;
} bam_data_t;

bam_data_t *bam_data_init();
void bam_data_dstry(bam_data_t *bam_dat);

int bam_data_fill_bcs(bam_data_t *b, str_map *bcs);


#endif // BAM_DAT_H
