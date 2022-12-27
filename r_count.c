
#include "r_count.h"
#include "str_util.h"
#include "gtf_anno.h"

/*******************************************************************************
 * allele counts
 ******************************************************************************/

ac_node *init_ac_node(){
    ac_node *n = calloc(1, sizeof(ac_node));

    if (n == NULL){
        err_msg(-1, 0, "init_ac_node: %s", strerror(errno));
        return NULL;
    }

    ac_node_set_zero(n);
    return n;
}

ac_node *destroy_ac_node(ac_node *n){
    ac_node *next = n->next;
    free(n);
    return next;
}

void ac_node_set_zero(ac_node *n){
    n->ix = -1;
    int i;
    for (i = 0; i < MAX_ALLELE; ++i) n->counts[i] = 0;
    n->next = NULL;
}

/*******************************************************************************
 * gene counts
 ******************************************************************************/

gc_node *init_gc_node(){
    gc_node *n = calloc(1, sizeof(gc_node));

    if (n == NULL){
        err_msg(-1, 0, "init_gc_node: %s", strerror(errno));
        return NULL;
    }

    gc_node_set_zero(n);

    return n;
}

gc_node *destroy_gc_node(gc_node *n){
    gc_node *next = n->next;
    free(n);
    return next;
}

void gc_node_set_zero(gc_node *n){
    n->ix = -1;
    int i;
    for (i = 0; i < N_SPL; ++i) n->counts[i] = 0;
    n->next = NULL;
}

/*******************************************************************************
 * barcode counts
 ******************************************************************************/

bam_ag_t *bam_ag_init(){
    bam_ag_t *agv = calloc(1, sizeof(bam_ag_t));
    if (agv == NULL){
        err_msg(-1, 0, "bc_ag_counts_init: %s", strerror(errno));
        return(NULL);
    }

    agv->gene_ix = NULL;
    agv->bc_ix = NULL;
    agv->pks = NULL;

    agv->bt_atac_ac = kh_init(kh_bt);
    agv->bt_atac_pc = kh_init(kh_bt);
    agv->bt_rna_ac = kh_init(kh_bt);
    agv->bt_rna_gc = kh_init(kh_bt);

    agv->rna_gcs_nz = 0;
    agv->rna_acs_nz = 0;
    agv->atac_acs_nz = 0;

    return(agv);
}

void bam_ag_dstry(bam_ag_t *agc){
    if (agc == NULL) return;

    if (agc->gene_ix) destroy_str_map(agc->gene_ix);
    if (agc->bc_ix) destroy_str_map(agc->bc_ix);

    if (agc->bt_atac_ac){
        khint_t k_bc;
        for (k_bc = kh_begin(agc->bt_atac_ac); k_bc != kh_end(agc->bt_atac_ac); ++k_bc){
            if (!kh_exist(agc->bt_atac_ac, k_bc)) continue;
            kbtree_t(str) *bt = kh_val(agc->bt_atac_ac, k_bc).b;
            kb_destroy(str, bt);
        }
        kh_destroy(kh_bt, agc->bt_atac_ac);
    }
    if (agc->bt_rna_ac){
        khint_t k_bc;
        for (k_bc = kh_begin(agc->bt_rna_ac); k_bc != kh_end(agc->bt_rna_ac); ++k_bc){
            if (!kh_exist(agc->bt_rna_ac, k_bc)) continue;
            kbtree_t(str) *bt = kh_val(agc->bt_rna_ac, k_bc).b;
            kb_destroy(str, bt);
        }
        kh_destroy(kh_bt, agc->bt_rna_ac);
    }
    if (agc->bt_rna_gc){
        khint_t k_bc;
        for (k_bc = kh_begin(agc->bt_rna_gc); k_bc != kh_end(agc->bt_rna_gc); ++k_bc){
            if (!kh_exist(agc->bt_rna_gc, k_bc)) continue;
            kbtree_t(str) *ct = kh_val(agc->bt_rna_gc, k_bc).b;
            kb_destroy(str, ct);
        }
        kh_destroy(kh_bt, agc->bt_rna_gc);
    }
    if (agc->bt_atac_pc){
        khint_t k_bc;
        for (k_bc = kh_begin(agc->bt_atac_pc); k_bc != kh_end(agc->bt_atac_pc); ++k_bc){
            if (!kh_exist(agc->bt_atac_pc, k_bc)) continue;
            kbtree_t(str) *bt = kh_val(agc->bt_atac_pc, k_bc).b;
            kb_destroy(str, bt);
        }
        kh_destroy(kh_bt, agc->bt_atac_pc);
    }
    free(agc);
}

int bam_ag_add_gv(bam_ag_t *agc, GenomeVar *gv){
    if (agc == NULL || gv == NULL)
        return err_msg(-1, 0, "bam_ag_add_gv_map: arguments are NULL");

    agc->gv = gv;
    return(0);
}

int bam_ag_add_gene_map(bam_ag_t *agc, str_map *sm){
    if (agc == NULL || sm == NULL)
        return err_msg(-1, 0, "bam_ag_add_gene_map: arguments are NULL");

    destroy_str_map(agc->gene_ix);
    agc->gene_ix = str_map_copy(sm);
    if (agc->gene_ix == NULL) return -1;
    return 0;
}

int bam_ag_add_bc_map(bam_ag_t *agc, str_map *sm){
    if (agc == NULL || sm == NULL)
        return err_msg(-1, 0, "bam_ag_add_bc_map: arguments are NULL");

    agc->bc_ix = str_map_copy(sm);
    if (agc->bc_ix == NULL) return -1;
    return 0;
}

int bam_ag_add_peaks(bam_ag_t *agc, iregs_t *pks){
    if (agc == NULL || pks == NULL)
        return err_msg(-1, 0, "bam_ag_add_bc_map: arguments are NULL");

    agc->pks = pks;
    return 0;
}

/*******************************************************************************
 * barcode counts
 ******************************************************************************/

int bam_rna_gc_count(bam_ag_t *agc, bam_data_t *bam_dat){

    bam_rna_t *br = bam_dat->rna;
    if (br == NULL)
        return 0;
    khint_t k_bc;
    for (k_bc = kh_begin(br->bc_rna); k_bc != kh_end(br->bc_rna); ++k_bc){
        if (!kh_exist(br->bc_rna, k_bc)) continue;
        char *bc_key = kh_key(br->bc_rna, k_bc);
        bc_rna_t *bc_rna = kh_val(br->bc_rna, k_bc);

        int bci, found;
        if ( (bci = add2str_map(agc->bc_ix, bc_key, &found)) < 0)
            return -1;
        bc_key = str_map_str(agc->bc_ix, bci);

        // add barcode to hash table
        int ret;
        khint_t k_gcs = kh_put(kh_bt, agc->bt_rna_gc, bc_key, &ret);
        if (ret == -1){
            return err_msg(-1, 0, "bam_rna_gc_count: failed to add barcode %s to "
                    "hash table", bc_key);
        } else if (ret == 0){
            return err_msg(-1, 0, "bam_rna_gc_count: barcode %s found twice\n", bc_key);
        }
        kbtree_t(str) *bt = kb_init(str, KB_DEFAULT_SIZE);
        if (bt == NULL)
            return err_msg(-1, 0, "bam_rna_gc_count: failed to init kbtree %s", 
                    strerror(errno));

        // add counts
        khint_t k_m;
        for (k_m = kh_begin(bc_rna->mols); k_m != kh_end(bc_rna->mols); ++k_m){
            if (!kh_exist(bc_rna->mols, k_m)) continue;
            rna_mol_t *mol = kh_val(bc_rna->mols, k_m);
            if (mol == NULL) continue;
            seq_gene_t *gene = mol->genes.head;
            size_t n_feat = mol->genes.n;
            if (n_feat > 1) continue; // only count UMIs overlapping one gene
            for (; gene != NULL; gene = gene->next){
                int32_t fix = gene->gene_id;
                uint8_t spl = gene->splice;
                if (spl > N_SPL)
                    return err_msg(-1, 0, "bam_rna_gc_count: "
                            "spl %u > %i", spl, N_SPL);
                cnt_node_t *p, t = {0};
                t.ix = (int)fix;
                t.counts[spl] = 1;
                p = kb_getp(str, bt, &t);
                if (!p){
                    kb_putp(str, bt, &t);
                    ++agc->rna_gcs_nz;
                }
                else ++p->counts[spl];
            }
        }
        kh_val(agc->bt_rna_gc, k_gcs).b = bt;
    }

    return 0;
}

int bam_rna_ac_count(bam_ag_t *agc, bam_data_t *bam_dat){

    bam_rna_t *ba = bam_dat->rna;
    if (ba == NULL)
        return 0;
    khint_t k_bc;
    for (k_bc = kh_begin(ba->bc_rna); k_bc != kh_end(ba->bc_rna); ++k_bc){
        if (!kh_exist(ba->bc_rna, k_bc)) continue;
        char *bc_key = kh_key(ba->bc_rna, k_bc);
        bc_rna_t *bc_rna = kh_val(ba->bc_rna, k_bc);

        int bci, found;
        if ( (bci = add2str_map(agc->bc_ix, bc_key, &found)) < 0)
            return -1;
        bc_key = str_map_str(agc->bc_ix, bci);

        // add barcode to hash table
        int ret;
        khint_t k_acs;

        // init tree
        k_acs = kh_put(kh_bt, agc->bt_rna_ac, bc_key, &ret);
        if (ret == -1){
            return err_msg(-1, 0, "bam_rna_ac_count: failed to add barcode %s to "
                    "hash table", bc_key);
        } else if (ret == 0){
            return err_msg(-1, 0, "bam_rna_ac_count: barcode %s found twice\n", bc_key);
        }
        kbtree_t(str) *bt = kb_init(str, KB_DEFAULT_SIZE);

        // add counts
        khint_t k_m;
        for (k_m = kh_begin(bc_rna->mols); k_m != kh_end(bc_rna->mols); ++k_m){
            if (!kh_exist(bc_rna->mols, k_m)) continue;
            rna_mol_t *mol = kh_val(bc_rna->mols, k_m);
            if (mol == NULL) continue;
            vac_t *vac = mol->vacs.head;
            for (; vac != NULL; vac = vac->next){
                int32_t vix = vac->vix;
                uint8_t allele = vac->allele;
                if (allele >= MAX_ALLELE)
                    return err_msg(-1, 0, "bam_rna_ac_count: "
                            "allele %u > %i", allele, MAX_ALLELE);
                cnt_node_t *p, t = {0};
                t.ix = vix;
                t.counts[allele] = 1;
                p = kb_getp(str, bt, &t);
                if (!p){
                    kb_putp(str, bt, &t);
                    ++agc->rna_acs_nz;
                }
                else ++p->counts[allele];
            }
        }
        kh_val(agc->bt_rna_ac, k_acs).b = bt;
    }

    return 0;
}

int bam_atac_ac_count(bam_ag_t *agc, bam_data_t *bam_dat){
    if (bam_dat == NULL)
        return err_msg(-1, 0, "bam_atac_ac_count: bam_dat or atac is null ");

    bam_atac_t *ba = bam_dat->atac;
    if (ba == NULL)
        return 0;
    if (ba->bc_dat == NULL)
        return err_msg(-1, 0, "bam_atac_ac_count: bc_dat is null ");
    if (agc->bc_ix == NULL)
        return err_msg(-1, 0, "bam_atac_ac_count: agc->bc_ix is null ");
    if (agc->bt_atac_ac == NULL)
        return err_msg(-1, 0, "bam_atac_ac_count: agc->bt_atac_ac is null ");

    uint32_t bc_num = 0;
    khint_t k_bc;
    for (k_bc = kh_begin(ba->bc_dat); k_bc != kh_end(ba->bc_dat); ++k_bc){
        if (!kh_exist(ba->bc_dat, k_bc)) continue;
        char *bc_key = kh_key(ba->bc_dat, k_bc);
        bc_atac_t *bc_atac = kh_val(ba->bc_dat, k_bc);
        if (bc_atac == NULL)
            return err_msg(-1, 0, "bam_atac_ac_count: bc_atac is null");

        if (bc_atac->frags == NULL)
            return err_msg(-1, 0, "bam_atac_ac_count: bc_atac->frags is null");

        bc_num++;
        int bci, found;
        if ( (bci = add2str_map(agc->bc_ix, bc_key, &found)) < 0)
            return -1;
        bc_key = str_map_str(agc->bc_ix, bci);

        // add barcode to a->gcs hash table
        khint_t k_acs;
        int ret;

        // init tree
        k_acs = kh_put(kh_bt, agc->bt_atac_ac, bc_key, &ret);
        if (ret == -1){
            return err_msg(-1, 0, "bam_atac_ac_count: failed to add barcode %s to "
                    "hash table", bc_key);
        } else if (ret == 0){
            return err_msg(-1, 0, "bam_atac_ac_count: barcode %s found twice\n", bc_key);
        }
        kbtree_t(str) *bt = kb_init(str, KB_DEFAULT_SIZE);

        // add counts
        uint32_t n_var = 0;
        khint_t k_m;
        for (k_m = kh_begin(bc_atac->frags); k_m != kh_end(bc_atac->frags); ++k_m){
            if (!kh_exist(bc_atac->frags, k_m)) continue;
            atac_frag_t *frag = kh_val(bc_atac->frags, k_m);
            if (frag == NULL) continue;
            vac_t *vac = frag->vacs.head;
            for (; vac != NULL; vac = vac->next){
                int32_t vix = vac->vix;
                uint8_t allele = vac->allele;
                if (allele >= MAX_ALLELE)
                    return err_msg(-1, 0, "bam_atac_ac_count: "
                            "allele %u > %i", allele, MAX_ALLELE);

                cnt_node_t *p, t = {0};
                t.ix = vix;
                t.counts[allele] = 1;
                p = kb_getp(str, bt, &t);
                if (!p){
                    kb_putp(str, bt, &t);
                    ++agc->atac_acs_nz;
                }
                else{
                    ++p->counts[allele];
                }
                ++n_var;
            }
        }
        kh_val(agc->bt_atac_ac, k_acs).b = bt;
    }

    return 0;
}

int bam_atac_pc_count(bam_ag_t *agc, bam_data_t *bam_dat){
    if (bam_dat == NULL)
        return err_msg(-1, 0, "bam_atac_pc_count: bam_dat or atac is null ");

    bam_atac_t *ba = bam_dat->atac;
    if (ba == NULL)
        return 0;
    if (ba->bc_dat == NULL)
        return err_msg(-1, 0, "bam_atac_pc_count: bc_dat is null ");
    if (agc->bc_ix == NULL)
        return err_msg(-1, 0, "bam_atac_pc_count: agc->bc_ix is null ");
    if (agc->bt_atac_pc == NULL)
        return err_msg(-1, 0, "bam_atac_pc_count: agc->by_atac_pc is null ");

    uint32_t bc_num = 0;
    khint_t k_bc;
    for (k_bc = kh_begin(ba->bc_dat); k_bc != kh_end(ba->bc_dat); ++k_bc){
        if (!kh_exist(ba->bc_dat, k_bc)) continue;
        char *bc_key       = kh_key(ba->bc_dat, k_bc);
        bc_atac_t *bc_atac = kh_val(ba->bc_dat, k_bc);
        if (bc_atac == NULL)
            return err_msg(-1, 0, "bam_atac_pc_count: bc_atac is null");

        ++bc_num;
        int bci, found;
        if ( (bci = add2str_map(agc->bc_ix, bc_key, &found)) < 0)
            return -1;
        bc_key = str_map_str(agc->bc_ix, bci);

        // add barcode to a->gcs hash table
        khint_t k_c;
        int ret;

        // init tree
        k_c = kh_put(kh_bt, agc->bt_atac_pc, bc_key, &ret);
        if (ret == -1){
            return err_msg(-1, 0, "bam_atac_pc_count: failed to add barcode %s to "
                    "hash table", bc_key);
        } else if (ret == 0){
            return err_msg(-1, 0, "bam_atac_pc_count: barcode %s found twice\n", bc_key);
        }
        kbtree_t(str) *bt = kb_init(str, KB_DEFAULT_SIZE);

        // loop over each fragment and add reg counts
        khint_t k_m;
        for (k_m = kh_begin(bc_atac->frags); k_m != kh_end(bc_atac->frags); ++k_m){
            if (!kh_exist(bc_atac->frags, k_m)) continue;
            atac_frag_t *f = kh_val(bc_atac->frags, k_m);
            if (f == NULL){
                printf("frag is null in bam_atac_pc_count, is this ok?\n");
                continue;
            }
            uint32_t c_add = 1;
            iregn_t pk = f->pks;
            int i;
            for (i = 0; i < pk.n; ++i){
                cnt_node_t *p, t = {0};
                t.ix = pk.ix[i];
                t.counts[0] = c_add;
                p = kb_getp(str, bt, &t);
                if (p == NULL){
                    kb_putp(str, bt, &t);
                    ++agc->atac_pcs_nz;
                } else {
                    p->counts[0] += c_add;
                }
            }
        }
        kh_val(agc->bt_atac_pc, k_c).b = bt;
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


int bam_ag_write_atac_ac(bam_ag_t *agc, char *fn){
    char delim[] = " ";

    if (agc == NULL)
        return(0);

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "bam_ag_write_atac_ac: "
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
            err_msg(-1, 0, "bam_ag_write_atac_ac: failed to open file %s", ofn[i]);
            return -1;
        }
    }

    size_t len;
    char *strs[3];

    // var bc non-zero lengths
    strs[0] = int2str((int)agc->gv->n_v, &len);
    strs[1] = int2str((int)agc->bc_ix->n, &len);
    strs[2] = int2str((int)agc->atac_acs_nz, &len);

    // write header
    for (i = 0; i < 3; ++i){
        if (mtx_write_hdr(fp[i], strs[0], strs[1], strs[2], delim) < 0){
            for (i = 0; i < 3; ++i) bgzf_close(fp[i]);
            for (i = 0; i < 3; ++i) free(ofn[i]);
            return err_msg(-1, 0, "bam_ag_write_atac_ac: failed to write to file %s", fn);
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

        khint_t k_bc = kh_get(kh_bt, agc->bt_atac_ac, bc_key);
        if (k_bc == kh_end(agc->bt_atac_ac)){
            bci++;
            continue;
        }
        kbtree_t(str) *bt = kh_val(agc->bt_atac_ac, k_bc).b;

        kbitr_t itr;
        kb_itr_first(str, bt, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(str, bt, &itr)){
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
                return err_msg(-1, 0, "bam_ag_write_atac_ac: failed to write to file %s", fn);
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

int bam_ag_write_atac_pc(bam_ag_t *agc, char *fn){
    char delim[] = " ";
    if (agc == NULL) return(0);
    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "bam_ag_write_atac_pc: "
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
        return err_msg(-1, 0, "bam_ag_write_atac_pc: failed to open file %s", ofn);

    size_t len;
    char *strs[3];

    strs[0] = int2str((int)agc->pks->n, &len);
    strs[1] = int2str((int)agc->bc_ix->n, &len);
    strs[2] = int2str((int)agc->atac_pcs_nz, &len);

    // write out header
    if (mtx_write_hdr(fp, strs[0], strs[1], strs[2], delim) < 0){
        bgzf_close(fp);
        free(ofn);
        return err_msg(-1, 0, "bam_ag_write_atac_pc: failed to write to file %s", fn);
    }
    for (i = 0; i < 3; ++i) free(strs[i]);

    size_t intstrp_m = 1;
    char *intstrp = malloc(sizeof(char) * intstrp_m);
    int k, bci = 1;
    for (k = 0; k < agc->bc_ix->n; ++k){
        char *bc_key = str_map_str(agc->bc_ix, k);
        if (bc_key == NULL) continue;

        khint_t k_bc = kh_get(kh_bt, agc->bt_atac_pc, bc_key);
        // if barcode not found in peak counts
        if (k_bc == kh_end(agc->bt_atac_pc)){
            ++bci;
            continue;
        }
        kbtree_t(str) *bt = kh_val(agc->bt_atac_pc, k_bc).b;

        kbitr_t itr;
        kb_itr_first(str, bt, &itr);
        for (; kb_itr_valid(&itr); kb_itr_next(str, bt, &itr)){
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

int bam_ag_write_rna_ac(bam_ag_t *agc, Annotation *anno, char *fn){
    char delim[] = " ";

    if (agc == NULL)
        return(0);

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "bam_ag_write_rna_ac: "
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
            err_msg(-1, 0, "bam_ag_write_rna_ac: failed to open file %s", ofn[i]);
            return -1;
        }
    }

    size_t len;
    char *strs[3];

    // var bc non-zero lengths
    strs[0] = int2str((int)agc->gv->n_v, &len);
    strs[1] = int2str((int)agc->bc_ix->n, &len);
    strs[2] = int2str((int)agc->rna_acs_nz, &len);

    // write header
    for (i = 0; i < 3; ++i){
        if (mtx_write_hdr(fp[i], strs[0], strs[1], strs[2], delim) < 0){
            for (i = 0; i < 3; ++i) bgzf_close(fp[i]);
            for (i = 0; i < 3; ++i) free(ofn[i]);
            return err_msg(-1, 0, "bam_ag_write_rna_ac: failed to write to file %s", fn);
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

        khint_t k_bc = kh_get(kh_bt, agc->bt_rna_ac, bc_key);
        if (k_bc == kh_end(agc->bt_rna_ac)){
            bci++;
            continue;
        }
        kbtree_t(str) *bt = kh_val(agc->bt_rna_ac, k_bc).b;

        kbitr_t itr;
        kb_itr_first(str, bt, &itr); 

        // loop over each variant
        for (; kb_itr_valid(&itr); kb_itr_next(str, bt, &itr)){
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
                return err_msg(-1, 0, "bam_ag_write_rna_ac: failed to write to file %s", fn);
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

int bam_ag_write_rna_gc(bam_ag_t *agc, Annotation *anno, char *fn){
    char delim[] = " ";

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "bam_ag_write: "
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
            return err_msg(-1, 0, "bam_ag_write: failed to write to file %s", fn);
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

        khint_t k_bc = kh_get(kh_bt, agc->bt_rna_gc, bc_key);
        if (k_bc == kh_end(agc->bt_rna_gc)){
            bci++;
            // there can be whitelist barcodes not added in data, 
            // so don't throw an error.
            continue;
        }
        kbtree_t(str) *bt = kh_val(agc->bt_rna_gc, k_bc).b;

        kbitr_t itr;
        kb_itr_first(str, bt, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(str, bt, &itr)){
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
                return err_msg(-1, 0, "bam_ag_write: failed to write to file %s", fn);
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

int bam_ag_write(bam_ag_t *agc, Annotation *anno, GenomeVar *gv, char *fn){
    if (agc == NULL)
        return err_msg(-1, 0, "bam_ag_write: agc is NULL");

    if (kh_size(agc->bt_rna_gc) > 0){
        if (bam_ag_write_rna_gc(agc, anno, fn) < 0)
            return -1;
    }
    if (kh_size(agc->bt_rna_ac) > 0){
        if (bam_ag_write_rna_ac(agc, anno, fn) < 0)
            return -1;
    }
    if (kh_size(agc->bt_atac_ac) > 0){
        if (bam_ag_write_atac_ac(agc, fn) < 0)
            return -1;
    }
    if (kh_size(agc->bt_atac_pc) > 0){
        if (bam_ag_write_atac_pc(agc, fn) < 0)
            return(-1);
    }

    BGZF *fp;
    char *ofn;
    // write genes 
    if (kh_size(agc->bt_rna_gc) > 0){
        char gene_fn[] = ".gene.txt.gz";
        ofn = strcat2((const char*)fn, (const char*)gene_fn);
        if (ofn == NULL) return -1;
        fp = bgzf_open(ofn, "wg1");
        if (fp == NULL)
            return err_msg(-1, 0, "bam_ag_write: failed to open file %s", ofn);
        if (write_gene_data(fp, anno, agc->gene_ix) < 0)
            return -1;
        bgzf_close(fp);
        // int write_gene_data(BGZF *fp, Annotation *anno, str_map *gene_ix){
        //     if (write_str_map(agc->gene_ix, ofn, ' ', '\n') < 0)
        free(ofn);
    }

    // write variants
    char var_fn[] = ".var.txt.gz";
    ofn = strcat2((const char*)fn, (const char*)var_fn);
    if (ofn == NULL) return -1;
    fp = bgzf_open(ofn, "wg1");
    if (fp == NULL)
        return err_msg(-1, 0, "bam_ag_write: failed to open file %s", ofn);
    int i;
    int32_t ne = agc->gv->n_e;
    for (i = 0; i < ne; ++i){
        Var *var = gv_vari(gv, i);
        if (var == NULL) continue;
        char *vid = var_id(agc->gv->vcf_hdr, var->b, '\t');
        if (bgzf_write(fp, vid, strlen(vid)) < 0)
            return err_msg(-1, 0, "bam_ag_write: failed to write variants to %s", ofn);
        free(vid);
    }
    bgzf_close(fp);
    free(ofn);

    // write peaks
    if (kh_size(agc->bt_atac_pc) > 0){
        char peak_fn[] = ".peaks.txt.gz";
        ofn = strcat2((const char *)fn, (const char *)peak_fn);
        if (ofn == NULL) return(-1);
        BGZF *fp = bgzf_open(ofn, "wg1");
        if (fp == NULL)
            return err_msg(-1, 0, "bam_ag_write: failed to open file %s", ofn);
        if (iregs_write(agc->pks, fp) < 0) return(-1);
        bgzf_close(fp);
        free(ofn);
    }
    
    // write barcodes
    char bc_fn[] = ".barcodes.txt.gz";
    ofn = strcat2((const char*)fn, (const char*)bc_fn);
    if (ofn == NULL) return -1;
    if (write_str_map(agc->bc_ix, ofn, ' ', '\n') < 0)
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

