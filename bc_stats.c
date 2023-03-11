
#include "bc_stats.h"
#include "str_util.h"
#include <stdlib.h>

int bc_stats_init(bc_stats_t *bcc){
    if (bcc == NULL)
        return err_msg(-1, 0, "bc_stats_init: bcc is NULL");

    bcc->bc = NULL;
    bcc->counts = 0;
    bcc->atac_counts = 0;
    bcc->rna_counts = 0;

    bcc->n_atac_vars = 0;
    bcc->n_rna_vars = 0;

    bcc->n_gene = 0;
    bcc->n_peak = 0;

    bcc->frip = -1;

    bcc->rna_mt = 0;
    bcc->atac_mt = 0;

    bcc->genes = kh_init(kh_cnt);
    bcc->atac_vars = kh_init(kh_cnt);
    bcc->rna_vars = kh_init(kh_cnt);
    bcc->peaks = kh_init(kh_cnt);

    if (bcc->genes == NULL || bcc->atac_vars == NULL 
        || bcc->rna_vars == NULL || bcc->peaks == NULL)
        return err_msg(-1, 0, "bc_stats_init: could not initialize hash tables");

    return(0);
}

bc_stats_t *bc_stats_alloc(){
    bc_stats_t *bcc = calloc(1, sizeof(bc_stats_t));
    if (bcc == NULL)
        err_msg(-1, 0, "bc_stats_alloc: %s", strerror(errno));

    if (bc_stats_init(bcc) < 0)
        return(NULL);

    return(bcc);
}

void bc_stats_free(bc_stats_t *bcc){
    if (bcc == NULL)
        return;

    kh_destroy(kh_cnt, bcc->genes);
    kh_destroy(kh_cnt, bcc->rna_vars);
    kh_destroy(kh_cnt, bcc->atac_vars);
    kh_destroy(kh_cnt, bcc->peaks);
}

void bc_stats_dstry(bc_stats_t *bcc){
    if (bcc == NULL)
        return;

    bc_stats_free(bcc);
    free(bcc);
}

