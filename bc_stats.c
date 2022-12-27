
#include "bc_stats.h"
#include "str_util.h"
#include <stdlib.h>

int bc_counts_init(bc_counts *bcc){
    if (bcc == NULL)
        return err_msg(-1, 0, "bc_counts_init: bcc is NULL");

    bcc->bc = NULL;
    bcc->counts = 0;
    bcc->atac_counts = 0;
    bcc->rna_counts = 0;
    bcc->genes = kh_init(kh_cnt);
    bcc->vars = kh_init(kh_cnt);
    bcc->peaks = kh_init(kh_cnt);
    if (bcc->genes == NULL || bcc->vars == NULL || bcc->peaks == NULL)
        return err_msg(-1, 0, "bc_counts_init: could not initialize hash tables");
    bcc->n_gene = 0;
    bcc->n_var = 0;
    bcc->n_peak = 0;
    bcc->frip = -1;

    return(0);
}

bc_counts *bc_counts_alloc(){
    bc_counts *bcc = calloc(1, sizeof(bc_counts));
    if (bcc == NULL)
        err_msg(-1, 0, "bc_counts_alloc: %s", strerror(errno));

    if (bc_counts_init(bcc) < 0)
        return(NULL);

    return(bcc);
}

void bc_counts_free(bc_counts *bcc){
    if (bcc == NULL)
        return;

    kh_destroy(kh_cnt, bcc->genes);
    kh_destroy(kh_cnt, bcc->vars);
    kh_destroy(kh_cnt, bcc->peaks);
}

void bc_counts_dstry(bc_counts *bcc){
    if (bcc == NULL)
        return;

    bc_counts_free(bcc);
    free(bcc);
}

bcs_stats *bc_stats_alloc(){
    bcs_stats *bcss = (bcs_stats *)calloc(1, sizeof(bcs_stats));
    if (bcss == NULL){
        err_msg(-1, 0, "bc_stats_alloc: %s", strerror(errno));
        return(NULL);
    }

    bcss->bcc = NULL;
    bcss->n_bcc = 0;
    bcss->kh_bcc = kh_init(kh_bcs);
    if(bcss->kh_bcc == NULL)
        return(NULL);
    bcss->barcodes = init_str_map();
    if(bcss->barcodes == NULL)
        return(NULL);

    return(bcss);
}

void bcs_stats_free(bcs_stats *bcss){
    if (bcss == NULL)
        return;

    khint_t k;
    for (k = kh_begin(bcss->kh_bcc); k != kh_end(bcss->kh_bcc); ++k){
        if (!kh_exist(bcss->kh_bcc, k)) continue;
        bc_counts *bcc = kh_val(bcss->kh_bcc, k);
        bc_counts_dstry(bcc);
    }
    kh_destroy(kh_bcs, bcss->kh_bcc);
    if (bcss->bcc) free(bcss->bcc);
    destroy_str_map(bcss->barcodes);
}

int bcs_stats_fill(bcs_stats *bcss, bam_rna_t *br, bam_atac_t *ba){
    if (br == NULL && ba == NULL)
        return err_msg(-1, 0, "bcs_stats_fill: arguments are NULL");

    khint_t k;
    if (br != NULL){
        for (k = kh_begin(br->bc_rna); k != kh_end(br->bc_rna); ++k){
            if (!kh_exist(br->bc_rna, k)) continue;
            char *bc_key = kh_key(br->bc_rna, k);
            bc_rna_t *bc_rna = kh_val(br->bc_rna, k);
            if (bc_key == NULL || bc_rna == NULL) continue;

            // check if barcode is present
            khint_t k_bc;
            int found;
            int bc_ix = add2str_map(bcss->barcodes, bc_key, &found);
            if (bc_ix < 0)
                return(-1);
            if (found == 0){
                char *bc_cpy = str_map_str(bcss->barcodes, bc_ix);
                int ret = 0;
                k_bc = kh_put(kh_bcs, bcss->kh_bcc, bc_cpy, &ret);
                if (ret < 0)
                    return err_msg(-1, 0, "bcs_stats_fill: failed to add to hash table");
                bc_counts *bctmp = bc_counts_alloc();
                if (bctmp == NULL)
                    return(-1);
                bctmp->bc = bc_cpy;
                kh_val(bcss->kh_bcc, k_bc) = bctmp;
            } else {
                k_bc = kh_get(kh_bcs, bcss->kh_bcc, bc_key);
            }

            if (k_bc == kh_end(bcss->kh_bcc))
                return err_msg(-1, 0, "bcs_stats_fill: bc not found, there is a bug");

            bc_counts *bcc = kh_val(bcss->kh_bcc, k_bc);

            khint_t k_mol;
            for (k_mol = kh_begin(bc_rna->mols); k_mol != kh_end(bc_rna->mols); ++k_mol){
                if (!kh_exist(bc_rna->mols, k_mol)) continue;
                rna_mol_t *m = kh_val(bc_rna->mols, k_mol);
                if (m == NULL) continue;
                seq_gene_t *gene = m->genes.head;
                for (gene = m->genes.head; gene != NULL; gene = gene->next){
                    int ret;
                    kh_put(kh_cnt, bcc->genes, gene->gene_id, &ret);
                    if (ret < 0)
                        return err_msg(-1, 0, "bcs_stats_fill: failed to add to gene ID");
                }
                vac_t *v = m->vacs.head;
                for (v = m->vacs.head; v != NULL; v = v->next){
                    int ret;
                    kh_put(kh_cnt, bcc->vars, v->vix, &ret);
                    if (ret < 0)
                        return err_msg(-1, 0, "bcs_stats_fill: failed to add to var ID");
                }
                ++bcc->counts;
                ++bcc->rna_counts;
            }
        }
    }
    if (ba != NULL){
        for (k = kh_begin(ba->bc_dat); k != kh_end(ba->bc_dat); ++k){
            if (!kh_exist(ba->bc_dat, k)) continue;
            char *bc_key = kh_key(ba->bc_dat, k);
            bc_atac_t *bc_atac = kh_val(ba->bc_dat, k);
            if (bc_key == NULL || bc_atac == NULL) continue;
            
            // check if barcode is present
            khint_t k_bc;
            int found;
            int bc_ix = add2str_map(bcss->barcodes, bc_key, &found);
            if (bc_ix < 0)
                return(-1);
            if (found == 0){
                char *bc_cpy = str_map_str(bcss->barcodes, bc_ix);
                int ret = 0;
                k_bc = kh_put(kh_bcs, bcss->kh_bcc, bc_cpy, &ret);
                if (ret < 0)
                    return err_msg(-1, 0, "bcs_stats_fill: failed to add to hash table");
                bc_counts *bctmp = bc_counts_alloc();
                if (bctmp == NULL)
                    return(-1);
                bctmp->bc = bc_cpy;
                kh_val(bcss->kh_bcc, k_bc) = bctmp;
            } else {
                k_bc = kh_get(kh_bcs, bcss->kh_bcc, bc_key);
            }

            if (k_bc == kh_end(bcss->kh_bcc))
                return err_msg(-1, 0, "bcs_stats_fill: bc not found, there is a bug");

            bc_counts *bcc = kh_val(bcss->kh_bcc, k_bc);

            uint32_t in_pk = 0;

            khint_t k_f;
            for (k_f = kh_begin(bc_atac->frags); k_f != kh_end(bc_atac->frags); ++k_f){
                if (!kh_exist(bc_atac->frags, k_f)) continue;
                atac_frag_t *f = kh_val(bc_atac->frags, k_f);
                if (f == NULL) continue;
                int p_i;
                for (p_i = 0; p_i < f->pks.n; ++p_i){
                    int ret;
                    kh_put(kh_cnt, bcc->peaks, f->pks.ix[p_i], &ret);
                    if (ret < 0)
                        return err_msg(-1, 0, "bcs_stats_fill: failed to add to peak ID");
                }
                if (f->pks.n) ++in_pk;
                vac_t *v = f->vacs.head;
                for (v = f->vacs.head; v != NULL; v = v->next){
                    int ret;
                    kh_put(kh_cnt, bcc->vars, v->vix, &ret);
                    if (ret < 0)
                        return err_msg(-1, 0, "bcs_stats_fill: failed to add to var ID");
                }
                ++bcc->counts;
                ++bcc->atac_counts;
            }
            if (bcc->atac_counts)
                bcc->frip = (float)in_pk / (float)bcc->atac_counts;
            else
                bcc->frip = -1;
        }
    }

    // fill n_gene and n_var
    for (k = kh_begin(bcss->kh_bcc); k != kh_end(bcss->kh_bcc); ++k){
        if (!kh_exist(bcss->kh_bcc, k)) continue;
        bc_counts *bcc = kh_val(bcss->kh_bcc, k);
        bcc->n_gene = kh_size(bcc->genes);
        bcc->n_var = kh_size(bcc->vars);
        bcc->n_peak = kh_size(bcc->peaks);
    }

    return(0);
}

int bcs_stats_fill_bcc(bcs_stats *bcss){
    if (bcss == NULL)
        return err_msg(-1, 0, "bcs_stats_fill_bcc: arg is NULL");

    size_t n_bcs = kh_size(bcss->kh_bcc);
    bcss->bcc = realloc(bcss->bcc, n_bcs * sizeof(bc_counts *));
    if (bcss->bcc == NULL)
        return err_msg(-1, 0, "bcs_stats_fill_bcc: %s", strerror(errno));

    size_t bc_ix = 0;
    khint_t k;
    for (k = kh_begin(bcss->kh_bcc); k != kh_end(bcss->kh_bcc); ++k){
        if (!kh_exist(bcss->kh_bcc, k)) continue;
        bcss->bcc[bc_ix++] = kh_val(bcss->kh_bcc, k);
    }
    if (bc_ix != n_bcs)
        return err_msg(-1, 0, "bcs_stats_fill_bcc: bc_ix != n_bcs");

    bcss->n_bcc = n_bcs;

    return(0);
}

int bcs_stats_sort(bcs_stats *bcss, int ascend){
    if (bcss == NULL)
        return err_msg(-1, 0, "bcs_stats_sort: arguments are NULL");

    if (bcss->bcc == NULL && bcs_stats_fill_bcc(bcss) < 0) return(-1);
    size_t n_bcs = bcss->n_bcc;
    

    if (ascend){
        qsort(bcss->bcc, n_bcs, sizeof(bc_counts *), cmp_bc_counts);
    } else {
        qsort(bcss->bcc, n_bcs, sizeof(bc_counts *), cmp_bc_counts_rev);
    }

    return(0);
}


uint32_t bcs_stats_min(bcs_stats *bcss, uint32_t top_n, int *ret){
    *ret = 0;
    if (bcss == NULL){
        *ret = err_msg(-1, 0, "bcs_stats_min: bcss is NULL");
        return(0);
    }

    // fill bcc array
    if (bcs_stats_fill_bcc(bcss) < 0){
        *ret = -1;
        return(0);
    }

    // sort bcc array
    if (bcs_stats_sort(bcss, 0) < 0){
        *ret = -1;
        return(0);
    }

    if (top_n >= (uint32_t)bcss->n_bcc){
        *ret = err_msg(-1, 0, "bcs_stats_min: top_n=%zu >= n_bcc=%zu", 
                top_n, bcss->n_bcc);
        return(0);
    }

    uint32_t min_c = bcss->bcc[top_n]->counts;

    // free bcc array
    free(bcss->bcc);
    bcss->bcc = NULL;
    bcss->n_bcc = 0;

    // get the count of the top_n ranked barcode.
    return(min_c);
}

void bcs_stats_print_sorted(bcs_stats *bcss, int64_t max_n){
    if (bcss == NULL) return;

    size_t n_bcs = bcss->n_bcc;

    max_n = max_n > n_bcs ? n_bcs : max_n;

    fprintf(stdout, "BC counts atac_counts rna_counts n_gene n_var\n");
    int i;
    for (i = 0; i < max_n; ++i){
        fprintf(stdout, "%s: ", bcss->bcc[i]->bc);
        fprintf(stdout, " %u", bcss->bcc[i]->counts);
        fprintf(stdout, " %u", bcss->bcc[i]->atac_counts);
        fprintf(stdout, " %u", bcss->bcc[i]->rna_counts);
        fprintf(stdout, " %u", bcss->bcc[i]->n_gene);
        fprintf(stdout, " %u\n", bcss->bcc[i]->n_var);
    }
}
