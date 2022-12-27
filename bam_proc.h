
#ifndef BAM_PROC_H
#define BAM_PROC_H

#include "clopts.h"
#include "atac_data.h"
#include "rna_data.h"
#include "bam_dat.h"

/* These are the main pileup functions
 *
 */
int run_atac(obj_pars *objs, bam_atac_t **bam_atac);
int run_rna(obj_pars *objs, bam_rna_t **bam_rna);

// generate and write out gene and variant counts
int bam_count(bam_data_t *bam_dat, obj_pars *objs, char *filename);

#endif // BAM_PROC_H
