
#include "bam_dat.h"
#include "atac_data.h"
#include "sam_read.h"
#include "counts.h"

bc_data_t *bc_data_init(){
    bc_data_t *bcdat = calloc(1, sizeof(bc_data_t));
    if (bcdat == NULL){
        err_msg(-1, 0, "bc_data_init: %s", strerror(errno));
        return(NULL);
    }
    if (mt_init(rna_dups, &bcdat->rna_dupls) < 0) {
        free(bcdat);
        return NULL;
    }
    if (mt_init(rna_mols, &bcdat->rna_mlcls) < 0) {
        free(bcdat);
        return NULL;
    }
    if (rna_dups_bag_init(&bcdat->rna_dups) < 0) {
        free(bcdat);
        return NULL;
    }
    if (rna_mlc_bag_init(&bcdat->rna_mlcs) < 0) {
        free(bcdat);
        return NULL;
    }
    if (mt_init(atac_pair, &bcdat->atac_pairs) < 0) {
        free(bcdat);
        return NULL;
    }
    if (atac_pair_bag_init(&bcdat->atac_prs) < 0) {
        free(bcdat);
        return NULL;
    }
    if (mt_init(atac_dups, &bcdat->atac_dupls) < 0) {
        free(bcdat);
        return NULL;
    }
    if (atac_dup_bag_init(&bcdat->atac_dups) < 0) {
        free(bcdat);
        return NULL;
    }
    if (mt_init(atac_frags, &bcdat->atac_frags) < 0) {
        free(bcdat);
        return NULL;
    }
    if (atac_frag_bag_init(&bcdat->atac_frgs) < 0) {
        free(bcdat);
        return NULL;
    }
    bcdat->bc_stats = NULL;

    int err;
    if ((err = pthread_mutex_init(&bcdat->bc_lock, NULL)) != 0){
        err_msg(-1, 0, "bc_data_init: failed to initialize mutex %i", err);
        return(NULL);
    }

    return(bcdat);
}

void bc_data_free_reads(bc_data_t *bcdat){
    if (bcdat == NULL) return;

    mt_free(rna_dups, &bcdat->rna_dupls);
    mt_free(rna_mols, &bcdat->rna_mlcls);
    rna_dups_bag_free(&bcdat->rna_dups);
    rna_mlc_bag_free(&bcdat->rna_mlcs);
    mt_free(atac_pair, &bcdat->atac_pairs);
    mt_free(atac_dups, &bcdat->atac_dupls);
    mt_free(atac_frags, &bcdat->atac_frags);
    atac_pair_bag_free(&bcdat->atac_prs);
    atac_dup_bag_free(&bcdat->atac_dups);
    atac_frag_bag_free(&bcdat->atac_frgs);
}

void bc_data_dstry(bc_data_t *bcdat){
    if (bcdat == NULL) return;

    bc_data_free_reads(bcdat);

    bc_stats_dstry(bcdat->bc_stats);

    pthread_mutex_destroy(&bcdat->bc_lock);

    free(bcdat);
}

// TODO: check that everything is freed appropriately
void bc_data_free_atac_pairs(bc_data_t *bcdat){
    if (bcdat == NULL) return;
    mt_free(atac_pair, &bcdat->atac_pairs);
    atac_pair_bag_free(&bcdat->atac_prs);
}

/*******************************************************************************
 * RNA
 ******************************************************************************/

int bc_data_rna_add_read(bc_data_t *bcdat, const rna_read1_t *r, umishort umih){
    // check input
    if (bcdat == NULL || r == NULL)
        return err_msg(-1, 0, "bc_data_rna_add_read: argument is null");

    rna_dups_t rna_dup;
    if (rna_dups_init(&rna_dup) < 0)
        return err_msg(-1, 0, "bc_data_rna_add_read: failed to initialize rna_dups");
    mt_node_t(rna_dups) *node = NULL;
    node = mt_add(rna_dups, &bcdat->rna_dupls, umih, rna_dup);
    if (node == NULL)
        return err_msg(-1, 0, "bc_data_rna_add_read: failed to add rna_dups to bcdat");
    rna_dups_t *p = &node->value;
    if (rna_dups_add_read(p, r) < 0)
        return err_msg(-1, 0, "bc_data_rna_add_read: failed to add read to rna_dups");

    if (rna_dups_bag_add_read(&bcdat->rna_dups, r, umih) < 0)
        return -1;

    return 0;
}

int bc_data_rna_dedup(bc_data_t *bcdat){
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_rna_dedup: bcdat is null");

    mt_itr_t(rna_dups) itrt;
    if (mt_itr_first(rna_dups, &bcdat->rna_dupls, &itrt) < 0)
        return err_msg(-1, 0, "bc_data_rna_dedup: failed to initialize iterator");

    for (; mt_itr_valid(&itrt); mt_itr_next(rna_dups, &itrt)) {
        umishort *umih = mt_itr_key(rna_dups, &itrt);
        rna_dups_t *dupls = mt_itr_val(rna_dups, &itrt);
        int ret;

        rna_mol_t *mol = rna_dups_dedup(dupls, &ret);
        if (ret < 0)
            return err_msg(-1, 0, "bc_data_rna_dedup: failed to deduplicate reads");
        if (mol == NULL) // if no reads, do nothing
            return 0;

        rna_dups_free(dupls);

        mt_node_t(rna_mols) *node = NULL;
        node = mt_add(rna_mols, &bcdat->rna_mlcls, *umih, *mol);
        if (node == NULL)
            return err_msg(-1, 0, "bc_data_rna_dedup: failed to add rna_mol to bcdat");

        free(mol);
    }
    
    rna_dups_bag_itr itr, *itrp = &itr;
    rna_dups_bag_itr_first(itrp, &bcdat->rna_dups);
    while (rna_dups_bag_itr_alive(itrp)) {
        umishort *umih = rna_dups_bag_itr_key(itrp);
        rna_dups_t *dups = rna_dups_bag_itr_val(itrp);
        int ret;

        rna_mol_t *mol = rna_dups_dedup(dups, &ret);
        if (ret < 0)
            return -1;
        rna_dups_free(dups);

        ret = rna_mlc_bag_add(&bcdat->rna_mlcs, mol, *umih);
        free(mol);

        rna_dups_bag_itr_next(itrp);
    }

    return(0);
}

int bc_data_rna_var_call(bc_data_t *bcdat, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (bcdat == NULL || gv == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_rna_var_call: argument is null");

    int n_add = 0;
    mt_itr_t(rna_mols) itrt;
    if (mt_itr_first(rna_mols, &bcdat->rna_mlcls, &itrt) < 0)
        return err_msg(-1, 0, "bc_data_rna_var_call: failed to initialize iterator");

    for (; mt_itr_valid(&itrt); mt_itr_next(rna_mols, &itrt)) {
        rna_mol_t *mol = mt_itr_val(rna_mols, &itrt);
        umishort *umih = mt_itr_key(rna_mols, &itrt);
        int ret = rna_mol_var_call(mol, gv, cmap, min_qual);
        if (ret < 0)
            return err_msg(-1, 0, "bc_data_rna_var_call: failed to call variant");
        n_add += ret;
    }

    n_add = 0;
    rna_mlc_bag_itr itr, *itrp = &itr;
    rna_mlc_bag_itr_first(itrp, &bcdat->rna_mlcs);
    while (rna_mlc_bag_itr_alive(itrp)) {
        rna_mol_t *mol = rna_mlc_bag_itr_val(itrp);
        int ret = rna_mol_var_call(mol, gv, cmap, min_qual);
        if (ret < 0)
            return -1;
        n_add += ret;

        rna_mlc_bag_itr_next(itrp);
    }

    return(n_add);
}

/*******************************************************************************
 * ATAC
 ******************************************************************************/

int bc_data_atac_add_read(bc_data_t *bcdat, const atac_read1_t *ar, qshort qname){
    if (bcdat == NULL || ar == NULL)
        return err_msg(-1, 0, "bc_data_atac_add_read: argument is null");

    atac_rd_pair_t atac_pair;
    atac_rd_pair_set0(&atac_pair);
    mt_node_t(atac_pair) *node = NULL;
    node = mt_add(atac_pair, &bcdat->atac_pairs, qname, atac_pair);
    int ret = atac_rd_pair_add_read(&node->value, ar);
    if (ret < 0)
        return err_msg(-1, 0, "bc_data_atac_add_read: failed to add read to atac pair");

    if (atac_pair_bag_add_read(&bcdat->atac_prs, ar, qname) < 0)
        return err_msg(-1, 0, "bc_data_atac_add_read: failed to add read to atac pair");

    return(0);
}

int bc_data_atac_get_dups(bc_data_t *bcdat) {
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_atac_get_dups: argument is null");

    mt_itr_t(atac_pair) itrt;
    if (mt_itr_first(atac_pair, &bcdat->atac_pairs, &itrt) < 0)
        return err_msg(-1, 0, "bc_data_atac_get_dups: failed to initialize iterator");

    for (; mt_itr_valid(&itrt); mt_itr_next(atac_pair, &itrt)) {
        atac_rd_pair_t *ap = mt_itr_val(atac_pair, &itrt);
        g_reg_pair reg_pair = get_reg_pair(ap->r1.loc, ap->r2.loc); // new key for dups
        mt_node_t(atac_dups) *node = mt_find(atac_dups, &bcdat->atac_dupls, reg_pair);
        if (node == NULL) {
            atac_dups_t dupls;
            atac_dups_init(&dupls);
            node = mt_add(atac_dups, &bcdat->atac_dupls, reg_pair, dupls);
        }
        if (atac_dups_add_pair(&node->value, ap) < 0)
            return err_msg(-1, 0, "bc_data_atac_get_dups: failed to add read to dups");
    }
    
    atac_pair_bag_itr itr, *itrp = &itr;
    atac_pair_bag_itr_first(itrp, &bcdat->atac_prs);
    for (; atac_pair_bag_itr_alive(itrp); atac_pair_bag_itr_next(itrp)) {
        atac_rd_pair_t *ap = atac_pair_bag_itr_val(itrp);
        assert(ap);
        if (atac_dup_bag_add_read(&bcdat->atac_dups, ap, 1) < 0)
            return -1;
    }
    // free atac pairs bag
    atac_pair_bag_free(&bcdat->atac_prs);
    return 0;
}

int bc_data_atac_dedup(bc_data_t *bcdat){
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_atac_dedup: bcdat is null");

    mt_itr_t(atac_dups) itrt;
    if (mt_itr_first(atac_dups, &bcdat->atac_dupls, &itrt) < 0)
        return err_msg(-1, 0, "bc_data_atac_get_dups: failed to initialize iterator");

    for (; mt_itr_valid(&itrt); mt_itr_next(atac_dups, &itrt)) {
        g_reg_pair *reg_pair = mt_itr_key(atac_dups, &itrt);
        atac_dups_t *dups = mt_itr_val(atac_dups, &itrt);
        int ret;

        atac_frag_t *frag = atac_dups_dedup(dups, &ret);
        if (ret < 0)
            return -1;

        // free space after deduplication
        atac_dups_free(dups);

        // add to frags
        mt_node_t(atac_frags) *node;
        node = mt_add(atac_frags, &bcdat->atac_frags, *reg_pair, *frag);
        free(frag); // free since we copied to new node
    }

    atac_dup_bag_itr itr, *itrp = &itr;
    atac_dup_bag_itr_first(itrp, &bcdat->atac_dups);
    for (; atac_dup_bag_itr_alive(itrp); atac_dup_bag_itr_next(itrp)) {
        g_reg_pair *reg = atac_dup_bag_itr_key(itrp);
        atac_dups_t *dups = atac_dup_bag_itr_val(itrp);
        int ret;

        atac_frag_t *frag = atac_dups_dedup(dups, &ret);
        if (ret < 0)
            return -1;
        atac_dups_free(dups);

        if (frag == NULL)
            continue;

        ret = atac_frag_bag_add(&bcdat->atac_frgs, frag, *reg);
        if (ret < 0)
            return -1;
        free(frag);

    }
    atac_dup_bag_free(&bcdat->atac_dups);

    return 0;
}

int bc_data_atac_var_call(bc_data_t *bcdat, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (bcdat == NULL || gv == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_atac_var_call: argument is null");

    int n_add = 0;
    mt_itr_t(atac_frags) itrt;
    if (mt_itr_first(atac_frags, &bcdat->atac_frags, &itrt) < 0)
        return err_msg(-1, 0, "bc_data_atac_var_call: failed to initialize iterator");

    for (; mt_itr_valid(&itrt); mt_itr_next(atac_frags, &itrt)) {
        atac_frag_t *frag = mt_itr_val(atac_frags, &itrt);
        g_reg_pair *reg = mt_itr_key(atac_frags, &itrt);
        int ret = atac_frag_var_call(frag, gv, cmap, min_qual);
        if (ret < 0) {
            print_g_region(stderr, reg->r1);
            print_g_region(stderr, reg->r2);
            return err_msg(-1, 0, "bc_data_atac_var_call: failed to call variant");
        }
        n_add += ret;
    }

    n_add = 0;
    atac_frag_bag_itr itr, *itrp = &itr;
    atac_frag_bag_itr_first(itrp, &bcdat->atac_frgs);
    for (; atac_frag_bag_itr_alive(itrp); atac_frag_bag_itr_next(itrp)) {
        g_reg_pair *reg = atac_frag_bag_itr_key(itrp);
        atac_frag_t *frag = atac_frag_bag_itr_val(itrp);

        int a = atac_frag_var_call(frag, gv, cmap, min_qual);
        if (a < 0){
            print_g_region(stderr, reg->r1);
            print_g_region(stderr, reg->r2);
            return(-1);
        }
        n_add += a;
    }

    return n_add;
}

int bc_data_atac_peak_call(bc_data_t *bcdat, iregs_t *pks, str_map *cmap){
    if (bcdat == NULL || pks == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_atac_peak_call: argument is null");

    int n_add = 0;

    mt_itr_t(atac_frags) itrt;
    if (mt_itr_first(atac_frags, &bcdat->atac_frags, &itrt) < 0)
        return err_msg(-1, 0, "bc_data_atac_var_call: failed to initialize iterator");

    for (; mt_itr_valid(&itrt); mt_itr_next(atac_frags, &itrt)) {
        atac_frag_t *frag = mt_itr_val(atac_frags, &itrt);
        g_reg_pair *reg = mt_itr_key(atac_frags, &itrt);
        int np;
        if ((np = atac_frag_peak_call(frag, *reg, pks, cmap)) < 0)
            return err_msg(-1, 0, "bc_data_atac_peak_call: failed to call peaks");
        n_add += np;
    }

    n_add = 0;
    atac_frag_bag_itr itr, *itrp = &itr;
    atac_frag_bag_itr_first(itrp, &bcdat->atac_frgs);
    for (; atac_frag_bag_itr_alive(itrp); atac_frag_bag_itr_next(itrp)) {
        g_reg_pair *reg = atac_frag_bag_itr_key(itrp);
        atac_frag_t *frag = atac_frag_bag_itr_val(itrp);

        int np;
        if ( (np = atac_frag_peak_call(frag, *reg, pks, cmap)) < 0)
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

void bam_data_free_bcs(bam_data_t *bam_data){
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
        bam_data->bc_data = NULL;
    }
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
        return err_msg(-1, 0, "bam_data_rna_add_read: bam_data->bc_data is null");

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

        // free dups struct to save space
        mt_free(rna_dups, &bc_dat->rna_dupls);
        rna_dups_bag_free(&bc_dat->rna_dups);
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

int bam_data_atac_get_dups(bam_data_t *bam_data) {
    if (bam_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_get_dups: bam_data is null");

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_atac_get_dups: bam_data->bc_data is null");

    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    khint_t kbc;
    for (kbc = kh_begin(bcs_hash); kbc != kh_end(bcs_hash); ++kbc){
        if (!kh_exist(bcs_hash, kbc)) continue;
        bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
        assert(bc_dat != NULL);

        if (bc_data_atac_get_dups(bc_dat) < 0)
            return -1;
    }
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

/*******************************************************************************
 * bc_stats
 ******************************************************************************/

int bc_data_rna_fill_stats(bc_stats_t *bcc, const sam_hdr_t *rna_hdr, bc_data_t *bc_dat){
    if (bcc == NULL || rna_hdr == NULL || bc_dat == NULL)
        return err_msg(-1, 0, "bc_data_rna_fill: arguments are null");

    float n_mt = 0;
    uint32_t in_gene = 0;

    mt_itr_t(rna_mols) itrt;
    if (mt_itr_first(rna_mols, &bc_dat->rna_mlcls, &itrt) < 0)
        return err_msg(-1, 0, "bc_data_rna_fill_stats: failed to initialize iterator");

    for (; mt_itr_valid(&itrt); mt_itr_next(rna_mols, &itrt)) {
        umishort *umih = mt_itr_key(rna_mols, &itrt);
        rna_mol_t *mol = mt_itr_val(rna_mols, &itrt);

        // count genes
        ml_node_t(seq_gene_l) *gn;
        for (gn = ml_begin(&mol->gl); gn; gn = ml_node_next(gn)) {
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
            return err_msg(-1, 0, "bc_data_rna_fill_stats: could not get rna chromosome "
                           "name for %i", mol->loc.rid);
        int is_mt = 0;
        if ( (is_mt = chr_is_mt(chr)) < 0 ) return(-1);
        if (is_mt) n_mt = n_mt + 1.0;
    }

    if (bcc->rna_counts) bcc->rna_mt = n_mt / (float)bcc->rna_counts;
    else bcc->rna_mt = 0;

    if (bcc->rna_counts)
        bcc->frig = (float)in_gene / (float)bcc->rna_counts;
    else
        bcc->frig = 0;

    bcc->n_gene = kh_size(bcc->genes);
    bcc->n_rna_vars = kh_size(bcc->rna_vars);

    return(0);
}

int bc_data_atac_fill_stats(bc_stats_t *bcc, const sam_hdr_t *atac_hdr, bc_data_t *bc_dat){
    if (bcc == NULL || atac_hdr == NULL || bc_dat == NULL)
        return err_msg(-1, 0, "bc_data_atac_fill: arguments are null");

    float n_mt = 0;
    uint32_t in_pk = 0;

    mt_itr_t(atac_frags) itrt;
    if (mt_itr_first(atac_frags, &bc_dat->atac_frags, &itrt) < 0)
        return err_msg(-1, 0, "bc_data_atac_fill_stats: failed to initialize iterator");

    for (; mt_itr_valid(&itrt); mt_itr_next(atac_frags, &itrt)) {
        atac_frag_t *frag = mt_itr_val(atac_frags, &itrt);
        g_reg_pair reg = *mt_itr_key(atac_frags, &itrt);

        // peaks present
        size_t p_i;
        for (p_i = 0; p_i < mv_size(&frag->pks); ++p_i){
            int ret;
            int32_t pk_ix = mv_i(&frag->pks, p_i);
            kh_put(kh_cnt, bcc->peaks, pk_ix, &ret);
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

