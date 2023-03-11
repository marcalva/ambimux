
#ifndef BAM_PROC_H
#define BAM_PROC_H

#include "clopts.h"
#include "atac_data.h"
#include "rna_data.h"
#include "bam_dat.h"

typedef struct {
    bam1_t **a;
    uint32_t n;
    obj_pars *objs;
    bam_data_t *bam_data;
    int ret;
} bam_thrd_args;

/* These are the main pileup functions
 *
 */
int proc_atac1(bam1_t *bam_r, obj_pars *objs, bam_data_t *bam_data);
int run_atac(obj_pars *objs, bam_data_t *bam_data);
int proc_rna1(bam1_t *bam_r, obj_pars *objs, bam_data_t *bam_data);
int run_rna(obj_pars *objs, bam_data_t *bam_data);

// generate and write out gene and variant counts
int bam_count(bam_data_t *bam_dat, obj_pars *objs, char *filename);

void *atac_thrd_fx(void *arg);
void *rna_thrd_fx(void *arg);
int rna_atac_thrd_call(bam1_t **a, uint32_t n, bam_data_t *bam_data, obj_pars *objs, 
        int rna_atac);
#endif // BAM_PROC_H
