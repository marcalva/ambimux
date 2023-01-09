
#include "bam_dat.h"

bc_data_t *bc_data_init(){
    bc_data_t *bcdat = calloc(1, sizeof(bc_data_t));
    if (bcdat == NULL){
        err_msg(-1, 0, "bc_data_init: %s", strerror(errno));
        return(NULL);
    }
    bcdat->rna_mols = kh_init(khrmn);
    bcdat->atac_pairs = kh_init(khap);
    bcdat->atac_frags = kh_init(khaf);
    bcdat->bc_stats = NULL;
    return(bcdat);
}

void bc_data_dstry(bc_data_t *bcdat){
    if (bcdat == NULL) return;

    if (bcdat->rna_mols != NULL){
        kh_destroy(khrmn, bcdat->rna_mols);
    }
    if (bcdat->atac_pairs != NULL){
        kh_destroy(khap, bcdat->atac_pairs);
    }
    if (bcdat->atac_frags != NULL){
        kh_destroy(khaf, bcdat->atac_frags);
    }
    bc_counts_dstry(bcdat->bc_stats);

    free(bcdat);
}

void bc_data_free_atac_pairs(bc_data_t *bcdat){
    if (bcdat == NULL) return;
    if (bcdat->atac_pairs == NULL) return;
    khash_t(khap) *pairs = bcdat->atac_pairs;
    khint_t k;
    for (k = kh_begin(pairs); k != kh_end(pairs); ++k){
        if (!kh_exist(pairs, k)) continue;
        atac_rd_pair_t *rp = kh_val(pairs, k);
        atac_rd_pair_dstry(rp);
    }
    kh_destroy(khap, pairs);
    bcdat->atac_pairs = NULL;
}

/*******************************************************************************
 * RNA
 ******************************************************************************/

int bc_data_add_rna_read(bc_data_t *bcdat, const rna_read1_t *r, const char *name){
    // check input
    if (bcdat == NULL || r == NULL)
        return err_msg(-1, 0, "bc_data_add_rna_read: 'bcdat' or 'r' is null");
    if (name == NULL)
        return err_msg(-1, 0, "bc_data_add_rna_read: 'name' is null");

    khash_t(khrmn) *mols = bcdat->rna_mols;
    khint_t k = kh_get(khrmn, mols, (char *)name); // add read by name (UMI barcode)

    // if read isn't present, allocate and add
    int kret = 0;
    if (k == kh_end(mols)){
        char *name_cpy = strdup(name);
        if (name_cpy == NULL)
            return err_msg(-1 , 0, "bc_data_add_rna_read: %s", strerror(errno));

        k = kh_put(khrmn, mols, name_cpy, &kret);
        if (kret < 0)
            return err_msg(-1, 0, "bc_data_add_rna_read: failed to add name to mols hash table");

        if ( (kh_val(mols, k) = rna_mol_alloc()) == NULL )
            return(-1);
    }

    rna_mol_t *mol = kh_val(mols, k);
    if (rna_mol_add_read(mol, r) < 0) return(-1);

    return 0;
}

int bc_data_rna_dedup(bc_data_t *bcdat){
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_rna_dedup: 'bcdat' is null");

    khash_t(khrmn) *mols = bcdat->rna_mols;

    khint_t kd;
    for (kd = kh_begin(mols); kd != kh_end(mols); ++kd){
        if (!kh_exist(mols, kd)) continue;

        rna_mol_t *mol = kh_val(mols, kd);
        if (mol == NULL)
            return err_msg(-1, 0, "bc_data_rna_dedup: a molecule is null");

        char *name = kh_key(mols, kd);
        if (name == NULL)
            return err_msg(-1, 0, "bc_data_rna_dedup: read name is null");

        // get best RNA molecule by majority
        if (rna_mol_dedup(mol) < 0) return(-1);
    }
    return(0);
}

int bc_data_rna_var_call(bc_data_t *bcdat, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual){
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_rna_var_call: br is NULL");

    if (bcdat->rna_mols == NULL)
        return err_msg(-1, 0, "bc_data_rna_var_call: mols is NULL");

    khash_t(khrmn) *mols = bcdat->rna_mols;
    int n_add = 0;
    khint_t km;
    for (km = kh_begin(mols); km != kh_end(mols); ++km){
        if (!kh_exist(mols, km)) continue;

        rna_mol_t *mol = kh_val(mols, km);
        if (mol == NULL)
            return err_msg(-1, 0, "bc_data_rna_var_call: a molecule is NULL");

        int a = rna_mol_var_call(mol, gv, cmap, min_qual);
        if (a < 0)
            return(-1);
        n_add += a;
    }
    return(n_add);
}

/*******************************************************************************
 * ATAC
 ******************************************************************************/

int bc_data_add_atac_read(bc_data_t *bcdat, const atac_read1_t *ar, qshort qname){
    if (bcdat == NULL || ar == NULL)
        return err_msg(-1, 0, "bc_add_atac_read: 'bcdat' or 'ar' is null");

    khash_t(khap) *pairs = bcdat->atac_pairs;
    khash_t(khaf) *frags = bcdat->atac_frags;

    int ret;
    khint_t kp;
    atac_rd_pair_t *rp;
    kp = kh_get(khap, pairs, qname);

    // if read pair isn't present, allocate and add read pair t
    if (kp == kh_end(pairs)){
        kp = kh_put(khap, pairs, qname, &ret);
        if (ret < 0)
            return err_msg(-1, 0, "bc_add_atac_read: failed to add qname to hash");
        
        if ( (kh_val(pairs, kp) = atac_rd_pair_init()) == NULL ) return(-1);
    }
    rp = kh_val(pairs, kp);

    // add read (copy) to read pair
    if ( atac_rd_pair_add_read(rp, ar) < 0 )
        return err_msg(-1, 0, "bc_add_atac_read: failed to add read to pair");

    // If a read pair was formed, add to duplicates and destroy read.
    if (rp->s == 2){
        if (khaf_add_dup(frags, rp) < 0)
            return(-1);
        atac_rd_pair_dstry(rp);
        kh_del(khap, pairs, kp);
        if (kh_size(pairs) > 4 && 
            (kh_n_buckets(pairs) > (10 * kh_size(pairs))) )
            kh_resize(khap, pairs, kh_size(pairs));
    }
    return(0);
}

int bc_data_dedup_atac(bc_data_t *bcdat){
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_dedup_atac: 'bcdat' is null");

    if (bcdat->atac_frags == NULL)
        return err_msg(-1, 0, "bc_dedup_atac: atac_frags is null");

    khint_t kf;
    khash_t(khaf) *frags = bcdat->atac_frags;

    // for each dup, deduplicate into frag and add to frags.
    for (kf = kh_begin(frags); kf != kh_end(frags); ++kf){
        if (!kh_exist(frags, kf)) continue;

        g_reg_pair reg = kh_key(frags, kf);

        // skip the duplicate if chimeric
        if (reg.r1.rid != reg.r2.rid){
            continue;
        }

        atac_frag_t *frag = kh_val(frags, kf);
        if (atac_frag_dedup(frag) < 0)
            return(-1);
    }

    return 0;
}

int bc_data_atac_var_call(bc_data_t *bcdat, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual){
    if (bcdat == NULL || gv == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_atac_var_call: bcdat is null");

    if (bcdat->atac_frags == NULL)
        return err_msg(-1, 0, "bc_data_atac_var_call: bcdat frags is null");

    khash_t(khaf) *frags = bcdat->atac_frags;
    int n_add = 0;
    khint_t kf;
    for (kf = kh_begin(frags); kf != kh_end(frags); ++kf){
        if (!kh_exist(frags, kf)) continue;

        atac_frag_t *f = kh_val(frags, kf);
        if (!f)
            return err_msg(-1, 0, "bc_data_atac_var_call: frag is NULL");

        int a = atac_frag_var_call(f, gv, cmap, min_qual);
        if (a < 0)
            return(-1);
        n_add += a;
    }
    return n_add;
}

int bc_data_atac_peak_call(bc_data_t *bcdat, iregs_t *pks, contig_map *cmap){
    if (bcdat == NULL || pks == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_atac_peak_call: arguments are null");

    if (bcdat->atac_frags == NULL)
        return err_msg(-1, 0, "bc_data_atac_peak_call: bcdat frags is null");

    khash_t(khaf) *frags = bcdat->atac_frags;
    int n_add = 0;
    khint_t kf;
    for (kf = kh_begin(frags); kf != kh_end(frags); ++kf){
        if (!kh_exist(frags, kf)) continue;

        atac_frag_t *f = kh_val(frags, kf);
        if (!f)
            return err_msg(-1, 0, "bc_data_atac_peak_call: frag is NULL");
        g_reg_pair reg = kh_key(frags, kf);

        int np;
        if ( (np = atac_frag_peak_call(f, reg, pks, cmap)) < 0)
            return(-1);
        n_add += np;
    }
    return(n_add);
}

/*******************************************************************************
 * bam_data_t
 ******************************************************************************/

bam_data_t *bam_data_init(){
    bam_data_t *b = (bam_data_t *)calloc(1, sizeof(bam_data_t));
    if (b == NULL){
        err_msg(-1, 0, "bam_data_init: %s", strerror(errno));
        return(NULL);
    }

    b->bc_data = kh_init(kh_bc_dat);
    b->has_rna = 0;
    b->has_atac = 0;

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

int bam_data_atac_free_pairs(bam_data_t *bam_data){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_free_pairs: 'bam_data' is null");
    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t kbc;
    for (kbc = kh_begin(bcs_hash); kbc != kh_end(bcs_hash); ++kbc){
        if (!kh_exist(bcs_hash, kbc)) continue;
        bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
        bc_data_free_atac_pairs(bc_dat);
    }
    return(0);
}


int bam_data_rna_add_read(bam_data_t *bam_data, const char *bc, 
        const rna_read1_t *r, const char *name){
    if (bam_data == NULL || r == NULL)
        return err_msg(-1, 0, "bam_rna_add_read: arguments are NULL");

    bam_data->has_rna = 1;

    int ret;
    khint_t kbc;
    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    kbc = kh_get(kh_bc_dat, bcs_hash, (char *)bc);
    // if barcode isn't present, add barcode and bc_data
    if (kbc == kh_end(bcs_hash)){
        char *bc_cpy = strdup(bc);
        kbc = kh_put(kh_bc_dat, bcs_hash, bc_cpy, &ret);
        if (ret < 0)
            return err_msg(-1, 0, "bam_rna_add_read: failed to add read to bcs");

        if ( (kh_val(bcs_hash, kbc) = bc_data_init()) == NULL )
            return(-1);
    }
    bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
    if (bc_data_add_rna_read(bc_dat, r, name) < 0)
        return err_msg(-1, 0, "bam_rna_add_read: failed to add read to bam");
    return(0);
}

int bam_data_rna_dedup(bam_data_t *bam_data){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_rna_dedup: 'bam_data' is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_rna_dedup: barcode is null");

        if (bc_data_rna_dedup(bc_dat) < 0)
            return -1;
    }
    return 0;
}

int bam_data_rna_var_call(bam_data_t *bam_data, GenomeVar *gv, 
        contig_map *cmap, uint8_t min_qual){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_rna_var_call: 'bam_data' is null");
    if (gv == NULL)
        return err_msg(-1, 0, "bam_data_rna_peak_call: 'gv' is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_rna_var_call: barcode is null");

        if (bc_data_rna_var_call(bc_dat, gv, cmap, min_qual) < 0)
            return -1;
    }
    return 0;
}

int bam_data_atac_add_read(bam_data_t *bam_data, const char *bc, 
        const atac_read1_t *r, qshort qname){
    if (bam_data == NULL || r == NULL)
        return err_msg(-1, 0, "bam_data_atac_add_read: 'bam_data' or 'r' is null");

    bam_data->has_atac = 1;

    int ret;
    khint_t kbc;
    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    kbc = kh_get(kh_bc_dat, bcs_hash, (char *)bc);
    // if barcode isn't present, add barcode and bc_data
    if (kbc == kh_end(bcs_hash)){
        char *bc_cpy = strdup(bc);
        kbc = kh_put(kh_bc_dat, bcs_hash, bc_cpy, &ret);
        if (ret < 0)
            return err_msg(-1, 0, "bam_data_atac_add_read: failed to add read to bcs");

        if ( (kh_val(bcs_hash, kbc) = bc_data_init()) == NULL )
            return(-1);
    }

    bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
    if (bc_data_add_atac_read(bc_dat, r, qname) < 0)
        return err_msg(-1, 0, "bam_data_atac_add_read: failed to add read to bam");
    return(0);
}

int bam_data_atac_dedup(bam_data_t *bam_data){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_dedup: 'bam_data' is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t kbc;
    for (kbc = kh_begin(bcs_hash); kbc != kh_end(bcs_hash); ++kbc){
        if (!kh_exist(bcs_hash, kbc)) continue;
        bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_atac_dedup: barcode is null");

        if (bc_data_dedup_atac(bc_dat) < 0)
            return -1;
    }
    return(0);
}

int bam_data_atac_var_call(bam_data_t *bam_data, GenomeVar *gv, 
        contig_map *cmap, uint8_t min_qual){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_var_call: 'bam_data' is null");
    if (gv == NULL)
        return err_msg(-1, 0, "bam_data_atac_peak_call: 'gv' is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_atac_var_call: barcode is null");

        if (bc_data_atac_var_call(bc_dat, gv, cmap, min_qual) < 0)
            return -1;
    }
    return 0;
}

int bam_data_atac_peak_call(bam_data_t *bam_data, iregs_t *pks, 
        contig_map *cmap){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_peak_call: 'bam_data' is null");
    if (pks == NULL)
        return err_msg(-1, 0, "bam_data_atac_peak_call: 'gv' is null");
    if (cmap == NULL)
        return err_msg(-1, 0, "bam_data_atac_peak_call: 'cmap' is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_atac_peak_call: barcode is null");

        if (bc_data_atac_peak_call(bc_dat, pks, cmap) < 0)
            return -1;
    }

    return 0;
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
    for (k = kh_begin(b->bc_data); k != kh_end(b->bc_data); ++k){
        if (!kh_exist(b->bc_data, k)) continue;
        char *bc_key = kh_key(b->bc_data, k);
        if (bc_key == NULL) continue;
        if (add2str_map(b->bcs, bc_key, &found) < 0)
            return(-1);
    }
    return(0);
}

/*******************************************************************************
 * bc_stats
 ******************************************************************************/

int bam_data_fill_stats(bam_data_t *bam_data){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_fill_stats: 'bam_data' is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_fill_stats: barcode is null");

        bc_counts *bcc = bc_counts_alloc();

        if (bam_data->has_rna){
            khash_t(khrmn) *mols = bc_dat->rna_mols;
            khint_t k_mol;
            for (k_mol = kh_begin(mols); k_mol != kh_end(mols); ++k_mol){
                if (!kh_exist(mols, k_mol)) continue;
                rna_mol_t *mol = kh_val(mols, k_mol);
                if (mol == NULL) continue;
                seq_gene_t *gene = mol->genes.head;
                for (gene = mol->genes.head; gene != NULL; gene = gene->next){
                    int ret;
                    kh_put(kh_cnt, bcc->genes, gene->gene_id, &ret);
                    if (ret < 0)
                        return err_msg(-1, 0, "bcs_stats_fill: failed to add to gene ID");
                }
                vac_t *v = mol->vacs.head;
                for (v = mol->vacs.head; v != NULL; v = v->next){
                    int ret;
                    kh_put(kh_cnt, bcc->vars, v->vix, &ret);
                    if (ret < 0)
                        return err_msg(-1, 0, "bcs_stats_fill: failed to add to var ID");
                }
                ++bcc->counts;
                ++bcc->rna_counts;
            }
        }
        if (bam_data->has_atac){
            uint32_t in_pk = 0;
            khash_t(khaf) *frags = bc_dat->atac_frags;
            khint_t k_f;
            for (k_f = kh_begin(frags); k_f != kh_end(frags); ++k_f){
                if (!kh_exist(frags, k_f)) continue;
                atac_frag_t *frag = kh_val(frags, k_f);
                if (frag == NULL) continue;
                int p_i;
                for (p_i = 0; p_i < frag->pks.n; ++p_i){
                    int ret;
                    kh_put(kh_cnt, bcc->peaks, frag->pks.ix[p_i], &ret);
                    if (ret < 0)
                        return err_msg(-1, 0, "bcs_stats_fill: failed to add to peak ID");
                }
                if (frag->pks.n) ++in_pk;
                vac_t *v = frag->vacs.head;
                for (v = frag->vacs.head; v != NULL; v = v->next){
                    int ret;
                    kh_put(kh_cnt, bcc->vars, v->vix, &ret);
                    if (ret < 0)
                        return err_msg(-1, 0, "bcs_stats_fill: failed to add to var ID");
                }
                ++bcc->counts;
                ++bcc->atac_counts;
                if (bcc->atac_counts)
                    bcc->frip = (float)in_pk / (float)bcc->atac_counts;
                else
                    bcc->frip = -1;
            }
        }

        // fill n_gene and n_var
        bcc->n_gene = kh_size(bcc->genes);
        bcc->n_var = kh_size(bcc->vars);
        bcc->n_peak = kh_size(bcc->peaks);
        bc_dat->bc_stats = bcc;
    }
    return 0;
}

uint32_t bam_data_count_of_n(bam_data_t *bam_data, int top_n, int *ret){
    *ret = 0;

    khint_t k;
    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    int n_bcs = kh_size(bcs_hash);
    bc_count_t *bc_array = (bc_count_t *)calloc(n_bcs, sizeof(bc_count_t));

    int bci = 0;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        char *bc_key = kh_key(bcs_hash, k);
        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_fill_stats: barcode is null");

        bc_array[bci].bc = bc_key;
        bc_array[bci].count = bc_dat->bc_stats->counts;
        ++bci;
    }
    qsort(bc_array, n_bcs, sizeof(bc_counts), cmp_bc_count_rev);

    if (top_n > n_bcs){
        *ret = err_msg(-1, 0, "bam_data_count_of_n: top_n%zu >= n_bcs=%zu", 
                top_n, n_bcs);
        return(0);
    }
    uint32_t n_c = bc_array[top_n].count;
    free(bc_array);
    return(n_c);
}

