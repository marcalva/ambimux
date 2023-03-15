
#include "r_count.h"
#include "str_util.h"
#include "gtf_anno.h"

/*******************************************************************************
 * barcode counts
 ******************************************************************************/

bc_counts_t *bc_counts_init(){
    bc_counts_t *bcc = calloc(1, sizeof(bc_counts_t));
    if (bcc == NULL){
        err_msg(-1, 0, "bc_counts_init: %s", strerror(errno));
        return(NULL);
    }

    bcc->atac_ac = kb_init(kh_cnode, KB_DEFAULT_SIZE);
    bcc->atac_pc = kb_init(kh_cnode, KB_DEFAULT_SIZE);
    bcc->rna_ac = kb_init(kh_cnode, KB_DEFAULT_SIZE);
    bcc->rna_gc = kb_init(kh_cnode, KB_DEFAULT_SIZE);
    if (bcc->atac_ac == NULL || bcc->atac_pc == NULL 
        || bcc->rna_ac == NULL || bcc->rna_gc == NULL){
        err_msg(-1, 0, "bc_counts_init: failed to init kbtree"); 
        return(NULL);
    }

    return(bcc);
}

void bc_counts_dstry(bc_counts_t *bcc){
    if (bcc == NULL) return;

    kb_destroy(kh_cnode, bcc->atac_ac);
    kb_destroy(kh_cnode, bcc->atac_pc);
    kb_destroy(kh_cnode, bcc->rna_ac);
    kb_destroy(kh_cnode, bcc->rna_gc);
    free(bcc);
}

/*******************************************************************************
 * BAM counts
 ******************************************************************************/

bam_counts_t *bam_counts_init(){
    bam_counts_t *agc = calloc(1, sizeof(bam_counts_t));
    if (agc == NULL){
        err_msg(-1, 0, "bc_ag_counts_init: %s", strerror(errno));
        return(NULL);
    }

    agc->gene_ix = NULL;
    agc->bc_ix = NULL;
    agc->pks = NULL;

    agc->bc_counts = kh_init(kh_bc_cnode);

    agc->has_atac_ac = 0;
    agc->has_atac_pc = 0;
    agc->has_rna_ac = 0;
    agc->has_rna_gc = 0;

    agc->rna_gcs_nz = 0;
    agc->rna_acs_nz = 0;
    agc->atac_acs_nz = 0;
    agc->atac_pcs_nz = 0;

    return(agc);
}

void bam_counts_dstry(bam_counts_t *agc){
    if (agc == NULL) return;

    if (agc->gene_ix) destroy_str_map(agc->gene_ix);
    if (agc->bc_ix) destroy_str_map(agc->bc_ix);

    if (agc->bc_counts){
        khint_t k_bc;
        for (k_bc = kh_begin(agc->bc_counts); k_bc != kh_end(agc->bc_counts); ++k_bc){
            if (!kh_exist(agc->bc_counts, k_bc)) continue;
            bc_counts_t *bcc = kh_val(agc->bc_counts, k_bc);
            bc_counts_dstry(bcc);
        }
        kh_destroy(kh_bc_cnode, agc->bc_counts);
    }

    free(agc);
}

int bam_counts_add_gv(bam_counts_t *agc, g_var_t *gv){
    if (agc == NULL || gv == NULL)
        return err_msg(-1, 0, "bam_counts_add_gv_map: argument is null");

    agc->gv = gv;
    return(0);
}

int bam_counts_add_gene_map(bam_counts_t *agc, str_map *sm){
    if (agc == NULL || sm == NULL)
        return err_msg(-1, 0, "bam_counts_add_gene_map: argument is null");

    destroy_str_map(agc->gene_ix);
    agc->gene_ix = str_map_copy(sm);
    if (agc->gene_ix == NULL) return -1;
    return 0;
}

int bam_counts_add_bc_map(bam_counts_t *agc, str_map *sm){
    if (agc == NULL || sm == NULL)
        return err_msg(-1, 0, "bam_counts_add_bc_map: argument is null");

    agc->bc_ix = str_map_copy(sm);
    if (agc->bc_ix == NULL) return -1;
    return 0;
}

int bam_counts_add_peaks(bam_counts_t *agc, iregs_t *pks){
    if (agc == NULL || pks == NULL)
        return err_msg(-1, 0, "bam_counts_add_peaks: argument is null");

    agc->pks = pks;
    return 0;
}

/*******************************************************************************
 * barcode counts
 ******************************************************************************/

int bam_counts_count(bam_counts_t *agc, bam_data_t *bam_data){
    if (agc == NULL || bam_data == NULL)
        return err_msg(-1, 0, "bam_counts_count: argument is null");

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_counts_count: 'bc_data' is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k_bc;
    for (k_bc = kh_begin(bcs_hash); k_bc != kh_end(bcs_hash); ++k_bc){
        if (!kh_exist(bcs_hash, k_bc)) continue;
        char *bc_key = kh_key(bcs_hash, k_bc);
        bc_data_t *bc_data = kh_val(bcs_hash, k_bc);

        int bci, found;
        if ( (bci = add2str_map(agc->bc_ix, bc_key, &found)) < 0)
            return -1;
        bc_key = str_map_str(agc->bc_ix, bci);

        // add barcode to kbtree
        int ret;
        khint_t k_bt = kh_put(kh_bc_cnode, agc->bc_counts, bc_key, &ret);
        if (ret == -1){
            return err_msg(-1, 0, "bam_counts_count: failed to add barcode %s to "
                    "hash table", bc_key);
        } else if (ret == 0){
            return err_msg(-1, 0, "bam_counts_count: barcode %s found twice\n", bc_key);
        }
        bc_counts_t *bcc = bc_counts_init();
        if (bcc == NULL) return(-1);

        if (bam_data->has_rna){
            khash_t(khrmn) *mols = bc_data->rna_mols;
            khint_t k_m;
            for (k_m = kh_begin(mols); k_m != kh_end(mols); ++k_m){
                if (!kh_exist(mols, k_m)) continue;
                rna_mol_t *mol = kh_val(mols, k_m);
                if (mol == NULL) continue;
                ml_node_t(seq_gene_l) *gn;
                if (ml_size(&mol->gl) != 1) continue; // discard multigene UMIs
                for (gn = ml_begin(&mol->gl); gn; gn = ml_node_next(gn)){
                    seq_gene_t gene = ml_node_val(gn);
                    int32_t fix = gene.gene_id;
                    uint8_t spl = gene.splice;
                    assert(spl <= N_SPL);
                    cnt_node_t *p, t;
                    memset(&t, 0, sizeof(cnt_node_t));
                    t.ix = (int)fix;
                    t.counts[spl] = 1;
                    p = kb_getp(kh_cnode, bcc->rna_gc, &t);
                    if (!p){
                        kb_putp(kh_cnode, bcc->rna_gc, &t);
                        ++agc->rna_gcs_nz;
                    } else {
                        p->counts[spl] += 1;
                    }
                    agc->has_rna_gc = 1;
                }
                ml_node_t(seq_vac_l) *vn;
                for (vn = ml_begin(&mol->vl); vn; vn = ml_node_next(vn)){
                    seq_vac_t vac = ml_node_val(vn);
                    int32_t vix = vac.vix;
                    uint8_t allele = vac.allele;
                    if (allele >= MAX_ALLELE)
                        return err_msg(-1, 0, "bam_counts_count: "
                                "allele %u > %i", allele, MAX_ALLELE);
                    cnt_node_t *p, t;
                    memset(&t, 0, sizeof(cnt_node_t));
                    t.ix = vix;
                    t.counts[allele] = 1;
                    p = kb_getp(kh_cnode, bcc->rna_ac, &t);
                    if (!p){
                        kb_putp(kh_cnode, bcc->rna_ac, &t);
                        ++agc->rna_acs_nz;
                    }
                    else {
                        p->counts[allele] += 1;
                    }
                    agc->has_rna_ac = 1;
                }
            }
        }
        if (bam_data->has_atac){
            khash_t(khaf) *frags = bc_data->atac_frags;
            khint_t k_m;
            for (k_m = kh_begin(frags); k_m != kh_end(frags); ++k_m){
                if (!kh_exist(frags, k_m)) continue;
                atac_frag_t *frag = kh_val(frags, k_m);
                if (frag == NULL){
                    printf("frag is null in bam_counts_count, is this ok?\n");
                    continue;
                }
                size_t ip;
                for (ip = 0; ip < mv_size(&frag->pks); ++ip){
                    cnt_node_t *p, t;
                    memset(&t, 0, sizeof(cnt_node_t));
                    t.ix = mv_i(&frag->pks, ip);
                    t.counts[0] = 1;
                    p = kb_getp(kh_cnode, bcc->atac_pc, &t);
                    if (p == NULL){
                        kb_putp(kh_cnode, bcc->atac_pc, &t);
                        ++agc->atac_pcs_nz;
                    } else {
                        p->counts[0] += 1;
                    }
                    agc->has_atac_pc = 1;
                }
                ml_node_t(seq_vac_l) *vn;
                for (vn = ml_begin(&frag->vl); vn; vn = ml_node_next(vn)){
                    seq_vac_t vac = ml_node_val(vn);
                    int32_t vix = vac.vix;
                    uint8_t allele = vac.allele;
                    if (allele >= MAX_ALLELE)
                        return err_msg(-1, 0, "bam_atac_ac_count: "
                                "allele %u > %i", allele, MAX_ALLELE);

                    cnt_node_t *p, t;
                    memset(&t, 0, sizeof(cnt_node_t));
                    t.ix = vix;
                    t.counts[allele] = 1;
                    p = kb_getp(kh_cnode, bcc->atac_ac, &t);
                    if (!p){
                        kb_putp(kh_cnode, bcc->atac_ac, &t);
                        ++agc->atac_acs_nz;
                    }
                    else{
                        p->counts[allele] += 1;
                    }
                    agc->has_atac_ac = 1;
                }
            }
        }
        kh_val(agc->bc_counts, k_bt) = bcc;
    }

    return(0);
}

/*******************************************************************************
 * write out counts
 ******************************************************************************/

int write_num(BGZF *fp, int n, char *c, char **intstrp, size_t *intstrp_m){
    int il, ret;
    if ( (il = int2strp(n, intstrp, intstrp_m)) < 0 ) return(-1);
    ret = bgzf_write(fp, *intstrp, il);
    if (ret < 0)
        return err_msg(-1, 0, "write_num: failed to write %i", n);
    ret = bgzf_write(fp, c, 1);
    if (ret < 0)
        return err_msg(-1, 0, "write_num: failed to write %s", c);
    return(0);
}


int bam_counts_write_atac_ac(bam_counts_t *agc, char *fn){
    char delim[] = " ";

    if (agc == NULL)
        return(0);

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "bam_counts_write_atac_ac: "
                "failed to create output directory for %s", fn);

    BGZF *fp[3];
    char *ofn[3];
    int ret;

    // create names and open files
    char mtx_fn[3][20] = {".atac.ac.ref.mtx.gz", ".atac.ac.alt.mtx.gz",".atac.ac.oth.mtx.gz"};
    int i;
    for (i = 0; i < 3; ++i){
        ofn[i] = strcat2((const char*)fn, (const char*)mtx_fn[i]);
        if (ofn[i] == NULL) return -1;
    }
    for (i = 0; i < 3; ++i){
        fp[i] = bgzf_open(ofn[i], "wg1");
        if (fp[i] == NULL){
            err_msg(-1, 0, "bam_counts_write_atac_ac: failed to open file %s", ofn[i]);
            return -1;
        }
    }

    size_t len;
    char *strs[3];

    // var bc non-zero lengths
    int n_vars = (int)mv_size(&agc->gv->vix2var);
    int n_bcs = (int)agc->bc_ix->n;
    int n_nz = (int)agc->atac_acs_nz;
    strs[0] = int2str(n_vars, &len);
    strs[1] = int2str(n_bcs, &len);
    strs[2] = int2str(n_nz, &len);

    // write header
    for (i = 0; i < 3; ++i){
        if (mtx_write_hdr(fp[i], strs[0], strs[1], strs[2], delim) < 0){
            for (i = 0; i < 3; ++i) bgzf_close(fp[i]);
            for (i = 0; i < 3; ++i) free(ofn[i]);
            return err_msg(-1, 0, "bam_counts_write_atac_ac: failed to write to file %s", fn);
        }
    }
    for (i = 0; i < 3; ++i) free(strs[i]);

    // write counts
    int il; // intstrp string length
    size_t intstrp_m = 1; // intstrp allocated size
    char *intstrp = malloc(sizeof(char) * intstrp_m);
    int k, bci = 1;
    for (k = 0; k < agc->bc_ix->n; ++k){
        char *bc_key = str_map_str(agc->bc_ix, k);
        if (bc_key == NULL) continue;

        khint_t k_bc = kh_get(kh_bc_cnode, agc->bc_counts, bc_key);
        if (k_bc == kh_end(agc->bc_counts)){
            bci++;
            continue;
        }
        bc_counts_t *bc_counts = kh_val(agc->bc_counts, k_bc);
        kbtree_t(kh_cnode) *bt = bc_counts->atac_ac;

        kbitr_t itr;
        kb_itr_first(kh_cnode, bt, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kh_cnode, bt, &itr)){
            cnt_node_t *n = &kb_itr_key(cnt_node_t, &itr);
            int vix = n->ix + 1; // convert index to 1-based for mtx file
            int bc_ix = k + 1;

            // write variant index
            if ((il = int2strp(vix, &intstrp, &intstrp_m)) < 0) return -1;
            for (i = 0; i < 3; ++i){
                ret = bgzf_write(fp[i], intstrp, il);
                ret = bgzf_write(fp[i], delim, 1);
            }

            // write barcode index
            if ((il = int2strp(bc_ix, &intstrp, &intstrp_m)) < 0) return -1;
            for (i = 0; i < 3; ++i){
                ret = bgzf_write(fp[i], intstrp, il);
                ret = bgzf_write(fp[i], delim, 1);
            }

            // write ref allele counts
            if ((il = int2strp((int)n->counts[0], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[0], intstrp, il);

            // write alt allele counts
            if ((il = int2strp((int)n->counts[1], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[1], intstrp, il);

            // write all other counts
            uint32_t oth_counts = 0;
            for (i = 2; i < MAX_ALLELE; ++i)
                oth_counts += n->counts[i];

            if ((il = int2strp((int)oth_counts, &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[2], intstrp, il);
            
            for (i = 0; i < 3; ++i)
                ret = bgzf_write(fp[i], "\n", 1);

            if (ret < 0)
                return err_msg(-1, 0, "bam_counts_write_atac_ac: failed to write to file %s", fn);
        }
        bci++;
    }
    if (k != (bci - 1)){
        fprintf(stdout, "k=%i bci=%i\n", k, bci);
    }
    for (i = 0; i < 3; ++i) bgzf_close(fp[i]);
    for (i = 0; i < 3; ++i) free(ofn[i]);

    free(intstrp);

    return(0);
}

int bam_counts_write_atac_pc(bam_counts_t *agc, char *fn){
    char delim[] = " ";
    if (agc == NULL) return(0);
    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "bam_counts_write_atac_pc: "
                "failed to create output directory for %s", fn);

    BGZF *fp;
    char *ofn;
    int i;

    // create names and open files
    char mtx_fn[20] = ".atac.peaks.mtx.gz";
    ofn = strcat2((const char*)fn, (const char*)mtx_fn);
    if (ofn == NULL) return(-1);
    fp = bgzf_open(ofn, "wg1");
    if (fp == NULL)
        return err_msg(-1, 0, "bam_counts_write_atac_pc: failed to open file %s", ofn);

    size_t len;
    char *strs[3];

    strs[0] = int2str((int)agc->pks->n, &len);
    strs[1] = int2str((int)agc->bc_ix->n, &len);
    strs[2] = int2str((int)agc->atac_pcs_nz, &len);

    // write out header
    if (mtx_write_hdr(fp, strs[0], strs[1], strs[2], delim) < 0){
        bgzf_close(fp);
        free(ofn);
        return err_msg(-1, 0, "bam_counts_write_atac_pc: failed to write to file %s", fn);
    }
    for (i = 0; i < 3; ++i) free(strs[i]);

    size_t intstrp_m = 1;
    char *intstrp = malloc(sizeof(char) * intstrp_m);
    int k, bci = 1;
    for (k = 0; k < agc->bc_ix->n; ++k){
        char *bc_key = str_map_str(agc->bc_ix, k);
        if (bc_key == NULL) continue;

        khint_t k_bc = kh_get(kh_bc_cnode, agc->bc_counts, bc_key);
        if (k_bc == kh_end(agc->bc_counts)){
            bci++;
            continue;
        }
        bc_counts_t *bc_counts = kh_val(agc->bc_counts, k_bc);
        kbtree_t(kh_cnode) *bt = bc_counts->atac_pc;

        kbitr_t itr;
        kb_itr_first(kh_cnode, bt, &itr);
        for (; kb_itr_valid(&itr); kb_itr_next(kh_cnode, bt, &itr)){
            cnt_node_t *n = &kb_itr_key(cnt_node_t, &itr);

            // write peak index
            if (write_num(fp, n->ix+1, delim, &intstrp, &intstrp_m) < 0) return(-1);

            // write barcode index
            if (write_num(fp, k+1, delim, &intstrp, &intstrp_m) < 0) return(-1);

            // write peak counts
            if (write_num(fp, (int)n->counts[0], "\n", &intstrp, &intstrp_m) < 0) return(-1);
        }
        ++bci;
    }
    if (k != (bci - 1))
        fprintf(stdout, "k=%i bci=%i\n", k, bci);

    bgzf_close(fp);
    free(ofn);

    free(intstrp);

    return(0);
}

int bam_counts_write_rna_ac(bam_counts_t *agc, char *fn){
    char delim[] = " ";

    if (agc == NULL)
        return(0);

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "bam_counts_write_rna_ac: "
                "failed to create output directory for %s", fn);

    BGZF *fp[3];
    char *ofn[3];
    int ret;

    // create names and open files
    char mtx_fn[3][20] = {".rna.ac.ref.mtx.gz", ".rna.ac.alt.mtx.gz",".rna.ac.oth.mtx.gz"};
    int i;
    for (i = 0; i < 3; ++i){
        ofn[i] = strcat2((const char*)fn, (const char*)mtx_fn[i]);
        if (ofn[i] == NULL) return -1;
    }
    for (i = 0; i < 3; ++i){
        fp[i] = bgzf_open(ofn[i], "wg1");
        if (fp[i] == NULL){
            err_msg(-1, 0, "bam_counts_write_rna_ac: failed to open file %s", ofn[i]);
            return -1;
        }
    }

    size_t len;
    char *strs[3];

    // var bc non-zero lengths
    int n_vars = (int)mv_size(&agc->gv->vix2var);
    int n_bcs = (int)agc->bc_ix->n;
    int n_nz = (int)agc->rna_acs_nz;
    strs[0] = int2str(n_vars, &len);
    strs[1] = int2str(n_bcs, &len);
    strs[2] = int2str(n_nz, &len);

    // write header
    for (i = 0; i < 3; ++i){
        if (mtx_write_hdr(fp[i], strs[0], strs[1], strs[2], delim) < 0){
            for (i = 0; i < 3; ++i) bgzf_close(fp[i]);
            for (i = 0; i < 3; ++i) free(ofn[i]);
            return err_msg(-1, 0, "bam_counts_write_rna_ac: failed to write to file %s", fn);
        }
    }
    for (i = 0; i < 3; ++i) free(strs[i]);

    // write counts
    int il; // intstrp string length
    size_t intstrp_m = 1; // intstrp allocated size
    char *intstrp = malloc(sizeof(char) * intstrp_m);
    int k, bci = 1;
    for (k = 0; k < agc->bc_ix->n; ++k){
        char *bc_key = str_map_str(agc->bc_ix, k);
        if (bc_key == NULL) continue;

        khint_t k_bc = kh_get(kh_bc_cnode, agc->bc_counts, bc_key);
        if (k_bc == kh_end(agc->bc_counts)){
            bci++;
            continue;
        }
        bc_counts_t *bc_counts = kh_val(agc->bc_counts, k_bc);
        kbtree_t(kh_cnode) *bt = bc_counts->rna_ac;

        kbitr_t itr;
        kb_itr_first(kh_cnode, bt, &itr); 

        // loop over each variant
        for (; kb_itr_valid(&itr); kb_itr_next(kh_cnode, bt, &itr)){
            cnt_node_t *n = &kb_itr_key(cnt_node_t, &itr);
            int vix = n->ix + 1; // convert index to 1-based for mtx file
            int bc_ix = k + 1;

            // write variant index
            if ((il = int2strp(vix, &intstrp, &intstrp_m)) < 0) return -1;
            for (i = 0; i < 3; ++i){
                ret = bgzf_write(fp[i], intstrp, il);
                ret = bgzf_write(fp[i], delim, 1);
            }

            // write barcode index
            if ((il = int2strp(bc_ix, &intstrp, &intstrp_m)) < 0) return -1;
            for (i = 0; i < 3; ++i){
                ret = bgzf_write(fp[i], intstrp, il);
                ret = bgzf_write(fp[i], delim, 1);
            }

            // write ref allele counts
            if ((il = int2strp((int)n->counts[0], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[0], intstrp, il);

            // write alt allele counts
            if ((il = int2strp((int)n->counts[1], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[1], intstrp, il);

            // write all other counts
            uint32_t oth_counts = 0;
            for (i = 2; i < MAX_ALLELE; ++i)
                oth_counts += n->counts[i];

            if ((il = int2strp((int)oth_counts, &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[2], intstrp, il);
            
            for (i = 0; i < 3; ++i)
                ret = bgzf_write(fp[i], "\n", 1);

            if (ret < 0)
                return err_msg(-1, 0, "bam_counts_write_rna_ac: failed to write to file %s", fn);
        }
        bci++;
    }
    if (k != (bci - 1)){
        fprintf(stdout, "k=%i bci=%i\n", k, bci);
    }
    for (i = 0; i < 3; ++i) bgzf_close(fp[i]);
    for (i = 0; i < 3; ++i) free(ofn[i]);

    free(intstrp);

    return(0);
}

int bam_counts_write_rna_gc(bam_counts_t *agc, char *fn){
    char delim[] = " ";

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "bam_counts_write: "
                "failed to create output directory for %s", fn);

    BGZF *fp[3];
    char *ofn[3];
    int ret;

    // create names and open files
    char mtx_fn[3][15] = {".gc.spl.mtx.gz", ".gc.uns.mtx.gz", ".gc.amb.mtx.gz"};
    int i;
    for (i = 0; i < 3; ++i){
        ofn[i] = strcat2((const char*)fn, (const char*)mtx_fn[i]);
        if (ofn[i] == NULL) return -1;
    }
    for (i = 0; i < 3; ++i){
        fp[i] = bgzf_open(ofn[i], "wg1");
        if (fp[i] == NULL){
            err_msg(-1, 0, "bc_gc_write: failed to open file %s", ofn[i]);
            return -1;
        }
    }

    size_t len;
    char *strs[3];

    // var bc non-zero lengths
    strs[0] = int2str((int)agc->gene_ix->n, &len);
    strs[1] = int2str((int)agc->bc_ix->n, &len);
    strs[2] = int2str((int)agc->rna_gcs_nz, &len);

    // write header
    for (i = 0; i < 3; ++i){
        if (mtx_write_hdr(fp[i], strs[0], strs[1], strs[2], delim) < 0){
            for (i = 0; i < 3; ++i) bgzf_close(fp[i]);
            for (i = 0; i < 3; ++i) free(ofn[i]);
            return err_msg(-1, 0, "bam_counts_write: failed to write to file %s", fn);
        }
    }
    for (i = 0; i < 3; ++i) free(strs[i]);

    // write counts
    int il; // intstrp string length
    size_t intstrp_m = 1; // intstrp allocated size
    char *intstrp = malloc(sizeof(char) * intstrp_m);
    int k, bci = 1;
    for (k = 0; k < agc->bc_ix->n; ++k){
        char *bc_key = str_map_str(agc->bc_ix, k);
        if (bc_key == NULL) continue;

        khint_t k_bc = kh_get(kh_bc_cnode, agc->bc_counts, bc_key);
        if (k_bc == kh_end(agc->bc_counts)){
            bci++;
            continue;
        }
        bc_counts_t *bc_counts = kh_val(agc->bc_counts, k_bc);
        kbtree_t(kh_cnode) *bt = bc_counts->rna_gc;

        kbitr_t itr;
        kb_itr_first(kh_cnode, bt, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kh_cnode, bt, &itr)){
            cnt_node_t *n = &kb_itr_key(cnt_node_t, &itr);
            int fix = n->ix + 1; // convert index to 1-based for mtx file
            int bc_ix = k + 1;

            for (i = 0; i < 3; ++i){
                // write gene index
                if ((il = int2strp(fix, &intstrp, &intstrp_m)) < 0) return -1;
                ret = bgzf_write(fp[i], intstrp, il);

                ret = bgzf_write(fp[i], delim, 1);

                // write barcode index
                if ((il = int2strp(bc_ix, &intstrp, &intstrp_m)) < 0) return -1;
                ret = bgzf_write(fp[i], intstrp, il);

                ret = bgzf_write(fp[i], delim, 1);
            }

            if ((il = int2strp((int)n->counts[SPLICE], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[0], intstrp, il);

            if ((il = int2strp((int)n->counts[UNSPLICE], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[1], intstrp, il);

            if ((il = int2strp((int)n->counts[AMBIG], &intstrp, &intstrp_m)) < 0) return -1;
            ret = bgzf_write(fp[2], intstrp, il);
            
            for (i = 0; i < 3; ++i)
                ret = bgzf_write(fp[i], "\n", 1);

            if (ret < 0)
                return err_msg(-1, 0, "bam_counts_write: failed to write to file %s", fn);
        }
        bci++;
    }
    if (k != (bci - 1)){
        fprintf(stdout, "k=%i bci=%i\n", k, bci);
    }
    for (i = 0; i < 3; ++i) bgzf_close(fp[i]);
    for (i = 0; i < 3; ++i) free(ofn[i]);

    free(intstrp);

    return 0;
}

int bam_counts_write(bam_counts_t *agc, gene_anno_t *anno, g_var_t *gv, char *fn){
    if (agc == NULL)
        return err_msg(-1, 0, "bam_counts_write: agc is NULL");

    if (agc->has_rna_ac){
        if (bam_counts_write_rna_gc(agc, fn) < 0)
            return -1;
    }
    if (agc->has_rna_ac){
        if (bam_counts_write_rna_ac(agc, fn) < 0)
            return -1;
    }
    if (agc->has_atac_ac){
        if (bam_counts_write_atac_ac(agc, fn) < 0)
            return -1;
    }
    if (agc->has_atac_pc){
        if (bam_counts_write_atac_pc(agc, fn) < 0)
            return(-1);
    }

    BGZF *fp;
    char *ofn;
    // write genes 
    if (agc->has_rna_gc){
        char gene_fn[] = ".gene.txt.gz";
        ofn = strcat2((const char*)fn, (const char*)gene_fn);
        if (ofn == NULL) return -1;
        fp = bgzf_open(ofn, "wg1");
        if (fp == NULL)
            return err_msg(-1, 0, "bam_counts_write: failed to open file %s", ofn);
        if (write_gene_data(fp, anno, agc->gene_ix) < 0)
            return -1;
        bgzf_close(fp);
        // int write_gene_data(BGZF *fp, gene_anno_t *anno, str_map *gene_ix){
        //     if (write_str_map(agc->gene_ix, ofn) < 0)
        free(ofn);
    }

    // write variants
    char var_fn[] = ".var.txt.gz";
    ofn = strcat2((const char*)fn, (const char*)var_fn);
    if (ofn == NULL) return -1;
    fp = bgzf_open(ofn, "wg1");
    if (fp == NULL)
        return err_msg(-1, 0, "bam_counts_write: failed to open file %s", ofn);
    int i;
    int32_t ne = mv_size(&gv->vix2var);
    for (i = 0; i < ne; ++i){
        var_t *var = gv_vari(gv, i);
        if (var == NULL) continue;
        char *vid = var_id(agc->gv->vcf_hdr, var->b, '\t');
        if (bgzf_write(fp, vid, strlen(vid)) < 0)
            return err_msg(-1, 0, "bam_counts_write: failed to write variants to %s", ofn);
        if (bgzf_write(fp, "\n", 1) < 0)
            return err_msg(-1, 0, "bam_counts_write: failed to write variants to %s", ofn);
        free(vid);
    }
    bgzf_close(fp);
    free(ofn);

    // write peaks
    if (agc->has_atac_pc){
        char peak_fn[] = ".peaks.txt.gz";
        ofn = strcat2((const char *)fn, (const char *)peak_fn);
        if (ofn == NULL) return(-1);
        fp = bgzf_open(ofn, "wg1");
        if (fp == NULL)
            return err_msg(-1, 0, "bam_counts_write: failed to open file %s", ofn);
        if (iregs_write(agc->pks, fp) < 0) return(-1);
        bgzf_close(fp);
        free(ofn);
    }
    
    // write barcodes
    char bc_fn[] = ".barcodes.txt.gz";
    ofn = strcat2((const char*)fn, (const char*)bc_fn);
    if (ofn == NULL) return -1;
    if (write_str_map(agc->bc_ix, ofn) < 0)
        return -1;
    free(ofn);
    return 0;
}

int mtx_write_hdr(BGZF *fp, char *len1, char *len2, char *len3, char *delim){
    char mtx_hdr[] = "%%MatrixMarket matrix coordinate integer general\n";
    int ret;
    ret = bgzf_write(fp, mtx_hdr, strlen(mtx_hdr));
    ret = bgzf_write(fp, "%\n", 2);

    ret = bgzf_write(fp, len1, strlen(len1));
    ret = bgzf_write(fp, delim, 1);
    ret = bgzf_write(fp, len2, strlen(len2));
    ret = bgzf_write(fp, delim, 1);
    ret = bgzf_write(fp, len3, strlen(len3));
    ret = bgzf_write(fp, "\n", 1);
    if (ret < 0)
        return err_msg(-1, 0, "mtx_write_hdr: failed to write header");
    return(0);
}

