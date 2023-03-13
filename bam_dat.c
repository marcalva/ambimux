
#include "bam_dat.h"
#include "sam_read.h"
#include "counts.h"

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

    if (bcdat->rna_mols == NULL || 
        bcdat->atac_pairs == NULL || 
        bcdat->atac_frags == NULL){
        err_msg(-1, 0, "bc_data_init: failed to initialize hash table");
        return(NULL);
    }

    ml_init(mlaf, &bcdat->frags_l);
    ml_init(mlar, &bcdat->mols_l);

    int err;
    if ((err = pthread_mutex_init(&bcdat->bc_lock, NULL)) != 0){
        err_msg(-1, 0, "bc_data_init: failed to initialize mutex %i", err);
        return(NULL);
    }

    return(bcdat);
}

void bc_data_dstry(bc_data_t *bcdat){
    if (bcdat == NULL) return;

    khint_t k;

    if (bcdat->rna_mols != NULL){
        for (k = kh_begin(bcdat->rna_mols); k != kh_end(bcdat->rna_mols); ++k){
            if (!kh_exist(bcdat->rna_mols, k)) continue;
            umishort umih = kh_key(bcdat->rna_mols, k);
            rna_mol_t *mol = kh_val(bcdat->rna_mols, k);
            rna_mol_dstry(mol);
        }
        kh_destroy(khrmn, bcdat->rna_mols);
    }
    if (bcdat->atac_pairs != NULL){
        for (k = kh_begin(bcdat->atac_pairs); k != kh_end(bcdat->atac_pairs); ++k){
            if (!kh_exist(bcdat->atac_pairs, k)) continue;
            atac_rd_pair_t *rp = kh_val(bcdat->atac_pairs, k);
            atac_rd_pair_dstry(rp);
        }
        kh_destroy(khap, bcdat->atac_pairs);
    }
    if (bcdat->atac_frags != NULL){
        for (k = kh_begin(bcdat->atac_frags); k != kh_end(bcdat->atac_frags); ++k){
            if (!kh_exist(bcdat->atac_frags, k)) continue;
            atac_frag_t *frag = kh_val(bcdat->atac_frags, k);
            atac_frag_dstry(frag);
        }
        kh_destroy(khaf, bcdat->atac_frags);
    }
    bc_stats_dstry(bcdat->bc_stats);

    ml_free(mlaf, &bcdat->frags_l);
    ml_free(mlar, &bcdat->mols_l);

    pthread_mutex_destroy(&bcdat->bc_lock);

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

int bc_data_fill_list(bc_data_t *bcdat){
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_fill_list: argument is null");

    khash_t(khrmn) *mol_hash = bcdat->rna_mols;
    khash_t(khaf) *frag_hash = bcdat->atac_frags;

    khint_t k_bc;
    for (k_bc = kh_begin(mol_hash); k_bc != kh_end(mol_hash); ++k_bc){
        if (!kh_exist(mol_hash, k_bc)) continue;
        rna_mol_t *mol = kh_val(mol_hash, k_bc);

        if (ml_insert(mlar, &bcdat->mols_l, mol, 1, 1) < 0)
            return(-1);
    }

    for (k_bc = kh_begin(frag_hash); k_bc != kh_end(frag_hash); ++k_bc){
        if (!kh_exist(frag_hash, k_bc)) continue;
        atac_frag_t *frag = kh_val(frag_hash, k_bc);

        if (ml_insert(mlaf, &bcdat->frags_l, frag, 1, 1) < 0)
            return(-1);
    }
    return(1);
}

/*******************************************************************************
 * RNA
 ******************************************************************************/

int bc_data_rna_add_read(bc_data_t *bcdat, const rna_read1_t *r, umishort umih){
    // check input
    if (bcdat == NULL || r == NULL)
        return err_msg(-1, 0, "bc_data_rna_add_read: argument is null");

    if (bcdat->rna_mols == NULL)
        return err_msg(-1, 0, "bc_data_rna_add_read: bdat->rna_mols is null");

    khash_t(khrmn) *mols = bcdat->rna_mols;
    khint_t k = kh_get(khrmn, mols, umih); // add read by UMI hash

    // if mol isn't present, allocate and add
    int kret = 0;
    if (k == kh_end(mols)){
        k = kh_put(khrmn, mols, umih, &kret);
        if (kret < 0)
            return err_msg(-1, 0, "bc_data_rna_add_read: failed to add UMI key "
                    "to mols hash table");

        if ( (kh_val(mols, k) = rna_mol_alloc()) == NULL )
            return(-1);
    }

    rna_mol_t *mol = kh_val(mols, k);
    if (rna_mol_add_read(mol, r) < 0) return(-1);

    return 0;
}

int bc_data_rna_dedup(bc_data_t *bcdat){
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_rna_dedup: bcdat is null");

    if (bcdat->rna_mols == NULL)
        return err_msg(-1, 0, "bc_data_rna_dedup: bcdat->rna_mols is null");

    khash_t(khrmn) *mols = bcdat->rna_mols;

    khint_t kd;
    for (kd = kh_begin(mols); kd != kh_end(mols); ++kd){
        if (!kh_exist(mols, kd)) continue;

        rna_mol_t *mol = kh_val(mols, kd);
        assert(mol != NULL);

        // get best RNA molecule by majority
        if (rna_mol_dedup(mol) < 0) return(-1);
    }
    return(0);
}

int bc_data_rna_var_call(bc_data_t *bcdat, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (bcdat == NULL || gv == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_rna_var_call: argument is null");

    if (bcdat->rna_mols == NULL)
        return err_msg(-1, 0, "bc_data_rna_var_call: bcdat->rna_mols is null");

    khash_t(khrmn) *mols = bcdat->rna_mols;
    int n_add = 0;
    khint_t km;
    for (km = kh_begin(mols); km != kh_end(mols); ++km){
        if (!kh_exist(mols, km)) continue;

        rna_mol_t *mol = kh_val(mols, km);
        assert(mol != NULL);

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

int bc_data_atac_add_read(bc_data_t *bcdat, const atac_read1_t *ar, qshort qname){
    if (bcdat == NULL || ar == NULL)
        return err_msg(-1, 0, "bc_data_atac_add_read: argument is null");

    khash_t(khap) *pairs = bcdat->atac_pairs;
    khash_t(khaf) *frags = bcdat->atac_frags;

    if (pairs == NULL || frags == NULL)
        return err_msg(-1, 0, "bc_data_atac_add_read: 'pairs' or 'frags' "
                "is null");

    int ret;
    khint_t kp;
    kp = kh_get(khap, pairs, qname);

    // if read pair isn't present, allocate and add read pair t
    if (kp == kh_end(pairs)){
        kp = kh_put(khap, pairs, qname, &ret);
        if (ret < 0)
            return err_msg(-1, 0, "bc_data_atac_add_read: failed to add "
                    "'qname' to hash");
        
        if ( (kh_val(pairs, kp) = atac_rd_pair_init()) == NULL ) return(-1);
    }
    atac_rd_pair_t *rp = kh_val(pairs, kp);

    // add read (copy) to read pair
    if ( atac_rd_pair_add_read(rp, ar) < 0 )
        return err_msg(-1, 0, "bc_data_atac_add_read: failed to add read to pair");

    // If a read pair was formed, add to duplicates and destroy read.
    if (rp->s == 2){
        if (khaf_add_dup(frags, rp, 1) < 0)
            return(-1);
        atac_rd_pair_dstry(rp);
        kh_del(khap, pairs, kp);
    }
    return(0);
}

int bc_data_atac_dedup(bc_data_t *bcdat){
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_atac_dedup: bcdat is null");

    if (bcdat->atac_frags == NULL)
        return err_msg(-1, 0, "bc_data_atac_dedup: bcdat->atac_frags is null");

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
        assert(frag != NULL);

        if (atac_frag_dedup(frag) < 0)
            return(-1);
    }
    return 0;
}

int bc_data_atac_var_call(bc_data_t *bcdat, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (bcdat == NULL || gv == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_atac_var_call: argument is null");

    if (bcdat->atac_frags == NULL)
        return err_msg(-1, 0, "bc_data_atac_var_call: bcdat->atac_frags is null");

    khash_t(khaf) *frags = bcdat->atac_frags;
    int n_add = 0;
    khint_t kf;
    for (kf = kh_begin(frags); kf != kh_end(frags); ++kf){
        if (!kh_exist(frags, kf)) continue;

        atac_frag_t *frag = kh_val(frags, kf);
        if (frag == NULL)
            return err_msg(-1, 0, "bc_data_atac_var_call: frag is NULL");

        int a = atac_frag_var_call(frag, gv, cmap, min_qual);
        if (a < 0){
            g_reg_pair greg = kh_key(frags, kf);
            print_g_region(stderr, greg.r1);
            print_g_region(stderr, greg.r2);
            return(-1);
        }
        n_add += a;
    }
    return n_add;
}

int bc_data_atac_peak_call(bc_data_t *bcdat, iregs_t *pks, str_map *cmap){
    if (bcdat == NULL || pks == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_atac_peak_call: argument is null");

    if (bcdat->atac_frags == NULL)
        return err_msg(-1, 0, "bc_data_atac_peak_call: bcdat->atac_frags is null");

    khash_t(khaf) *frags = bcdat->atac_frags;
    int n_add = 0;
    khint_t kf;
    for (kf = kh_begin(frags); kf != kh_end(frags); ++kf){
        if (!kh_exist(frags, kf)) continue;

        atac_frag_t *frag = kh_val(frags, kf);
        assert(frag != NULL);
        g_reg_pair reg = kh_key(frags, kf);

        int np;
        if ( (np = atac_frag_peak_call(frag, reg, pks, cmap)) < 0)
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
    b->has_stats = 0;

    b->bcs = init_str_map();
    if (b->bcs == NULL)
        return(NULL);

    int err;
    if ((err = pthread_mutex_init(&b->bam_lock, NULL)) != 0){
        err_msg(-1, 0, "bam_data_init: failed to initialize mutex errno=%i", err);
        return(NULL);
    }

    return(b);
}

void bam_data_dstry(bam_data_t *bam_data){
    if (bam_data == NULL)
        return;

    if (bam_data->bc_data != NULL){
        khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
        khint_t kbc;
        for (kbc = kh_begin(bcs_hash); kbc != kh_end(bcs_hash); ++kbc){
            if (!kh_exist(bcs_hash, kbc)) continue;
            char *bc_key = kh_key(bcs_hash, kbc);
            free(bc_key);
            bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
            bc_data_dstry(bc_dat);
        }
        kh_destroy(kh_bc_dat, bam_data->bc_data);
    }

    destroy_str_map(bam_data->bcs);

    pthread_mutex_destroy(&bam_data->bam_lock);

    free(bam_data);
}

int bam_data_atac_free_pairs(bam_data_t *bam_data){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_free_pairs: bam_data is null");
    if (bam_data->bc_data == NULL)
        return(0);
    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t kbc;
    for (kbc = kh_begin(bcs_hash); kbc != kh_end(bcs_hash); ++kbc){
        if (!kh_exist(bcs_hash, kbc)) continue;
        bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
        bc_data_free_atac_pairs(bc_dat);
    }
    return(0);
}

int bam_data_fill_bcs(bam_data_t *bam_data, str_map *bcs){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_fill_bcs: argument is NULL");

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_fill_bcs: bam_data->bc_data is null");

    if (bam_data->bcs){
        destroy_str_map(bam_data->bcs);
        bam_data->bcs = init_str_map();
    }

    // fill in barcodes from bcs argument if present
    if (bcs != NULL){
        destroy_str_map(bam_data->bcs);
        bam_data->bcs = str_map_copy(bcs);
        return(0);
    }

    // fill in barcodes from data
    // barcode is only added if not present.
    int found;
    khint_t k;
    for (k = kh_begin(bam_data->bc_data); k != kh_end(bam_data->bc_data); ++k){
        if (!kh_exist(bam_data->bc_data, k)) continue;
        char *bc_key = kh_key(bam_data->bc_data, k);
        if (bc_key == NULL) continue;
        if (add2str_map(bam_data->bcs, bc_key, &found) < 0)
            return(-1);
    }
    return(0);
}

/*******************************************************************************
 * BAM RNA
 ******************************************************************************/

int bam_data_rna_add_read(bam_data_t *bam_data, const char *bc, 
        const rna_read1_t *r, umishort umih){
    if (bam_data == NULL || bc == NULL || r == NULL)
        return err_msg(-1, 0, "bam_data_rna_add_read: argument is null");

    bam_data->has_rna = 1;

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_rna_dedup: bam_data->bc_data is null");

    if (pthread_mutex_lock(&bam_data->bam_lock) != 0)
        return err_msg(-1, 0, "bam_data_rna_add_read: failed mutex lock");

    int ret;
    khint_t kbc;
    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    kbc = kh_get(kh_bc_dat, bcs_hash, (char *)bc);
    // if barcode isn't present, add barcode and bc_data
    if (kbc == kh_end(bcs_hash)){
        char *bc_cpy = strdup(bc);
        if (bc_cpy == NULL){
            pthread_mutex_unlock(&bam_data->bam_lock);
            return err_msg(-1, 0, "bam_data_rna_add_read: %s", strerror(errno));
        }

        kbc = kh_put(kh_bc_dat, bcs_hash, bc_cpy, &ret);
        if (ret < 0){
            pthread_mutex_unlock(&bam_data->bam_lock);
            return err_msg(-1, 0, "bam_data_rna_add_read: failed to add read to bcs");
        }

        if ( (kh_val(bcs_hash, kbc) = bc_data_init()) == NULL ){
            pthread_mutex_unlock(&bam_data->bam_lock);
            return(-1);
        }
    }
    bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
    if (pthread_mutex_unlock(&bam_data->bam_lock) != 0)
        return err_msg(-1, 0, "bam_data_rna_add_read: failed mutex unlock");

    if (pthread_mutex_lock(&bc_dat->bc_lock) != 0)
        return err_msg(-1, 0, "bam_data_rna_add_read: failed mutex lock");
    ret = bc_data_rna_add_read(bc_dat, r, umih);
    if (pthread_mutex_unlock(&bc_dat->bc_lock) != 0)
        return err_msg(-1, 0, "bam_data_rna_add_read: failed mutex unlock");
    if (ret < 0)
        return err_msg(-1, 0, "bam_data_rna_add_read: failed to add read to bam");
    return(0);
}

int bam_data_rna_dedup(bam_data_t *bam_data){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_rna_dedup: bam_data is null");

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_rna_dedup: bam_data->bc_data is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        assert(bc_dat != NULL);

        if (bc_data_rna_dedup(bc_dat) < 0)
            return -1;
    }
    return 0;
}

int bam_data_rna_var_call(bam_data_t *bam_data, g_var_t *gv, 
        str_map *cmap, uint8_t min_qual){
    if (bam_data == NULL || gv == NULL || cmap == NULL)
        return err_msg(-1, 0, "bam_data_rna_var_call: argument is null");

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_rna_var_call: bam_data->bc_data is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        assert(bc_dat != NULL);

        if (bc_data_rna_var_call(bc_dat, gv, cmap, min_qual) < 0)
            return -1;
    }
    return 0;
}

/*******************************************************************************
 * BAM ATAC
 ******************************************************************************/

int bam_data_atac_add_read(bam_data_t *bam_data, const char *bc, 
        const atac_read1_t *r, qshort qname){
    if (bam_data == NULL || bc == NULL || r == NULL)
        return err_msg(-1, 0, "bam_data_atac_add_read: argument is null");

    bam_data->has_atac = 1;

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_add_read: 'bc_data' is null");

    if (pthread_mutex_lock(&bam_data->bam_lock) != 0)
        return err_msg(-1, 0, "bam_data_atac_add_read: failed bam mutex lock");

    int ret;
    khint_t kbc;
    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    kbc = kh_get(kh_bc_dat, bcs_hash, (char *)bc);
    // if barcode isn't present, add barcode and bc_data
    if (kbc == kh_end(bcs_hash)){
        char *bc_cpy = strdup(bc);
        if (bc_cpy == NULL)
            return err_msg(-1, 0, "bam_data_atac_add_read: %s", strerror(errno));

        kbc = kh_put(kh_bc_dat, bcs_hash, bc_cpy, &ret);
        if (ret < 0)
            return err_msg(-1, 0, "bam_data_atac_add_read: failed to add read to bcs");

        if ( (kh_val(bcs_hash, kbc) = bc_data_init()) == NULL )
            return(-1);
    }
    if (pthread_mutex_unlock(&bam_data->bam_lock) != 0)
        return err_msg(-1, 0, "bam_data_atac_add_read: failed bam mutex unlock");

    bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
    if (pthread_mutex_lock(&bc_dat->bc_lock) != 0)
        return err_msg(-1, 0, "bam_data_atac_add_read: failed bc mutex lock");
    ret = bc_data_atac_add_read(bc_dat, r, qname);
    if (pthread_mutex_unlock(&bc_dat->bc_lock) != 0)
        return err_msg(-1, 0, "bam_data_atac_add_read: failed bc mutex unlock");
    if (ret < 0)
        return err_msg(-1, 0, "bam_data_atac_add_read: failed to add read to bam");
    return(0);
}

int bam_data_atac_dedup(bam_data_t *bam_data){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_dedup: bam_data is null");

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_dedup: bam_data->bc_data is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t kbc;
    for (kbc = kh_begin(bcs_hash); kbc != kh_end(bcs_hash); ++kbc){
        if (!kh_exist(bcs_hash, kbc)) continue;
        bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
        assert(bc_dat != NULL);

        if (bc_data_atac_dedup(bc_dat) < 0)
            return -1;
    }
    return(0);
}

int bam_data_atac_var_call(bam_data_t *bam_data, g_var_t *gv, 
        str_map *cmap, uint8_t min_qual){
    if (bam_data == NULL || gv == NULL || cmap == NULL)
        return err_msg(-1, 0, "bam_data_atac_var_call: argument is null");

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_var_call: bam_data->bc_data is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        assert(bc_dat != NULL);

        if (bc_data_atac_var_call(bc_dat, gv, cmap, min_qual) < 0)
            return -1;
    }
    return 0;
}

int bam_data_atac_peak_call(bam_data_t *bam_data, iregs_t *pks, 
        str_map *cmap){
    if (bam_data == NULL || pks == NULL || cmap == NULL)
        return err_msg(-1, 0, "bam_data_atac_peak_call: argument is null");

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_peak_call: bam_data->bc_data is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        assert(bc_dat != NULL);

        if (bc_data_atac_peak_call(bc_dat, pks, cmap) < 0)
            return -1;
    }

    return 0;
}

int bam_data_fill_list(bam_data_t *bam_data){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_fill_list: argument is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        assert(bc_dat != NULL);

        if (bc_data_fill_list(bc_dat) < 0)
            return(-1);
    }

    return(0);
}

/*******************************************************************************
 * bc_stats
 ******************************************************************************/

int bc_data_rna_fill_stats(bc_stats_t *bcc, const sam_hdr_t *rna_hdr, bc_data_t *bc_dat){
    if (bcc == NULL || rna_hdr == NULL || bc_dat == NULL)
        return err_msg(-1, 0, "bc_data_rna_fill: arguments are null");

    if (bc_dat->rna_mols == NULL)
        return err_msg(-1, 0, "bc_data_rna_fill: rna_mols is null");

    khash_t(khrmn) *mols = bc_dat->rna_mols;

    float n_mt = 0;
    khint_t k_mol;
    for (k_mol = kh_begin(mols); k_mol != kh_end(mols); ++k_mol){
        if (!kh_exist(mols, k_mol)) continue;
        rna_mol_t *mol = kh_val(mols, k_mol);
        if (mol == NULL) continue;

        // count genes
        ml_node_t(seq_gene_l) *gn;
        for (gn = ml_begin(&mol->gl); gn; gn = ml_node_next(gn)){
            seq_gene_t gene = ml_node_val(gn);
            int ret;
            kh_put(kh_cnt, bcc->genes, gene.gene_id, &ret);
            if (ret < 0)
                return err_msg(-1, 0, "bc_data_rna_fill_stats: failed to add to gene ID");
        }

        // count variant calls
        ml_node_t(seq_vac_l) *vn;
        for (vn = ml_begin(&mol->vl); vn; vn = ml_node_next(vn)){
            int ret;
            seq_vac_t v = ml_node_val(vn);
            kh_put(kh_cnt, bcc->rna_vars, v.vix, &ret);
            if (ret < 0)
                return err_msg(-1, 0, "bc_data_rna_fill_stats: failed to add to var ID");
        }

        ++bcc->counts;
        ++bcc->rna_counts;

        // count MT
        const char *chr = sam_hdr_tid2name(rna_hdr, mol->loc.rid);
        if (chr == NULL)
            return err_msg(-1, 0, "bc_data_rna_fill_stats: could not get rna chromosome name for %i", mol->loc.rid);
        int is_mt = 0;
        if ( (is_mt = chr_is_mt(chr)) < 0 ) return(-1);
        if (is_mt) n_mt = n_mt + 1.0;
    }
    if (bcc->rna_counts) bcc->rna_mt = n_mt / (float)bcc->rna_counts;
    else bcc->rna_mt = 0;

    bcc->n_gene = kh_size(bcc->genes);
    bcc->n_rna_vars = kh_size(bcc->rna_vars);

    return(0);
}

int bc_data_atac_fill_stats(bc_stats_t *bcc, const sam_hdr_t *atac_hdr, bc_data_t *bc_dat){
    if (bcc == NULL || atac_hdr == NULL || bc_dat == NULL)
        return err_msg(-1, 0, "bc_data_atac_fill: arguments are null");

    if (bc_dat->atac_frags == NULL)
        return err_msg(-1, 0, "bc_data_atac_fill: atac_frags is null");

    uint32_t in_pk = 0;
    khash_t(khaf) *frags = bc_dat->atac_frags;
    float n_mt = 0;
    khint_t k_f;
    for (k_f = kh_begin(frags); k_f != kh_end(frags); ++k_f){
        if (!kh_exist(frags, k_f)) continue;
        atac_frag_t *frag = kh_val(frags, k_f);
        g_reg_pair reg = kh_key(frags, k_f);
        if (frag == NULL) continue;

        // peak counts
        size_t p_i;
        for (p_i = 0; p_i < mv_size(&frag->pks); ++p_i){
            int ret;
            kh_put(kh_cnt, bcc->peaks, mv_i(&frag->pks, p_i), &ret);
            if (ret < 0)
                return err_msg(-1, 0, "bc_data_atac_fill_stats: failed to add to peak ID");
        }
        if (mv_size(&frag->pks)) ++in_pk;

        // count variant calls
        ml_node_t(seq_vac_l) *vn;
        for (vn = ml_begin(&frag->vl); vn; vn = ml_node_next(vn)){
            seq_vac_t v = ml_node_val(vn);
            int ret;
            kh_put(kh_cnt, bcc->atac_vars, v.vix, &ret);
            if (ret < 0)
                return err_msg(-1, 0, "bc_data_atac_fill_stats: failed to add to var ID");
        }

        ++bcc->counts;
        ++bcc->atac_counts;

        // count MT
        const char *chr = sam_hdr_tid2name(atac_hdr, reg.r1.rid);
        if (chr == NULL)
            return err_msg(-1, 0, "bc_data_atac_fill_stats: could not get atac chromosome name for %i", reg.r1.rid);
        int is_mt = 0;
        if ( (is_mt = chr_is_mt(chr)) < 0 ) return(-1);
        if (is_mt) n_mt = n_mt + 1.0;
    }
    if (bcc->atac_counts) bcc->atac_mt = n_mt / (float)bcc->atac_counts;
    else bcc->atac_mt = 0;

    if (bcc->atac_counts)
        bcc->frip = (float)in_pk / (float)bcc->atac_counts;
    else
        bcc->frip = 0;

    bcc->n_peak = kh_size(bcc->peaks);
    bcc->n_atac_vars = kh_size(bcc->atac_vars);

    return(0);
}

int bam_data_fill_stats(bam_data_t *bam_data, 
        const sam_hdr_t *rna_hdr, const sam_hdr_t *atac_hdr){
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_fill_stats: 'bam_data' is null");

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_fill_stats: 'bc_data' is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_fill_stats: barcode is null");

        bc_stats_t *bcc = bc_stats_alloc();

        if (bam_data->has_rna){
            if (bc_data_rna_fill_stats(bcc, rna_hdr, bc_dat) < 0) return(-1);
        }
        if (bam_data->has_atac){
            if (bc_data_atac_fill_stats(bcc, atac_hdr, bc_dat) < 0) return(-1);
        }
        bc_dat->bc_stats = bcc;
    }
    bam_data->has_stats = 1;
    return 0;
}

uint32_t bam_data_count_of_n(bam_data_t *bam_data, int top_n, int *ret){
    *ret = 0;
    if (bam_data == NULL){
        *ret = err_msg(-1, 0, "bam_data_count_of_n: 'bam_data' is null");
        return(0);
    }
    if (bam_data->bc_data == NULL){
        *ret = err_msg(-1, 0, "bam_data_count_of_n: 'bc_data' is null");
        return(0);
    }

    khint_t k;
    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    int n_bcs = kh_size(bcs_hash);
    bc_count_t *bc_array = (bc_count_t *)calloc(n_bcs, sizeof(bc_count_t));
    if (bc_array == NULL){
        *ret = err_msg(-1, 0, "bam_data_count_of_n: %s", strerror(errno));
        return(0);
    }

    if (top_n < 0 || top_n > n_bcs){
        *ret = err_msg(-1, 0, "bam_data_count_of_n: top_n=%zu must be between ", 
                ">= 0 and <= %i", top_n, n_bcs);
        return(0);
    }

    int bci = 0;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        char *bc_key = kh_key(bcs_hash, k);
        if (bc_key == NULL){
            *ret = err_msg(-1, 0, "bam_data_count_of_n: 'bc_key' is null");
            return(0);
        }
        bc_data_t *bc_dat = kh_val(bcs_hash, k);
        if (bc_dat == NULL){
            *ret = err_msg(-1, 0, "bam_data_count_of_n: 'bc_dat' is null");
            return(0);
        }
        if (bc_dat->bc_stats == NULL){
            *ret = err_msg(-1, 0, "bam_data_count_of_n: bc 'bc_stats' is null");
            return(0);
        }

        bc_array[bci].bc = bc_key;
        bc_array[bci].count = bc_dat->bc_stats->counts;
        ++bci;
    }
    if (bci != n_bcs){
        *ret = err_msg(-1, 0, "bam_data_count_of_n: bci=%i != n_bcs=%i", bci, n_bcs);
        return(0);
    }
    qsort(bc_array, n_bcs, sizeof(bc_count_t), cmp_bc_count_rev);

    uint32_t n_c = bc_array[top_n].count;
    free(bc_array);
    return(n_c);
}

