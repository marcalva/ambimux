
#include "bam_dat.h"

bam_data_t *bam_data_init(){
    bam_data_t *b = (bam_data_t *)calloc(1, sizeof(bam_data_t));
    if (b == NULL){
        err_msg(-1, 0, "bam_data_init: %s", strerror(errno));
        return(NULL);
    }

    b->rna = NULL;
    b->atac = NULL;
    b->bcs = init_str_map();
    if (b->bcs == NULL)
        return(NULL);

    return(b);
}

void bam_data_dstry(bam_data_t *bam_dat){
    if (bam_dat == NULL)
        return;

    if (bam_dat->rna)
        bam_rna_dstry(bam_dat->rna);

    if (bam_dat->atac)
        bam_atac_dstry(bam_dat->atac);

    destroy_str_map(bam_dat->bcs);
    free(bam_dat);
}

int bam_data_fill_bcs(bam_data_t *b, str_map *bcs){
    if (b == NULL)
        return err_msg(-1, 0, "bam_data_fill_bcs: arguments are NULL");

    if (b->bcs){
        destroy_str_map(b->bcs);
        b->bcs = init_str_map();
    }

    // fill in barcodes from bcs argument
    if (bcs != NULL){
        b->bcs = str_map_copy(bcs);
        return(0);
    }

    // fill in barcodes from data
    // barcode is only added if not present.
    int found;
    khint_t k;
    if (b->rna != NULL){
        for (k = kh_begin(b->rna->bc_rna); k != kh_end(b->rna->bc_rna); ++k){
            if (!kh_exist(b->rna->bc_rna, k)) continue;
            char *bc_key = kh_key(b->rna->bc_rna, k);
            if (bc_key == NULL) continue;
            if (add2str_map(b->bcs, bc_key, &found) < 0)
                return(-1);
        }
    }
    if (b->atac != NULL){
        for (k = kh_begin(b->atac->bc_dat); k != kh_end(b->atac->bc_dat); ++k){
            if (!kh_exist(b->atac->bc_dat, k)) continue;
            char *bc_key = kh_key(b->atac->bc_dat, k);
            if (bc_key == NULL) continue;
            if (add2str_map(b->bcs, bc_key, &found) < 0)
                return(-1);
        }
    }

    return(0);
}

