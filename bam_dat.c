
#include "bam_dat.h"
#include "atac_data.h"
#include "kbtree.h"
#include "rna_data.h"
#include "sam_read.h"
#include "counts.h"
#include "str_util.h"

bc_data_t *bc_data_init(){
    bc_data_t *bcdat = calloc(1, sizeof(bc_data_t));
    if (bcdat == NULL){
        err_msg(-1, 0, "bc_data_init: %s", strerror(errno));
        return(NULL);
    }
    bcdat->rna_dups = kh_init(kh_rd);
    if (bcdat->rna_dups == NULL){
        free(bcdat);
        err_msg(-1, 0, "bc_data_init: %s", strerror(errno));
        return(NULL);
    }
    bcdat->rna_mols = ml_alloc(ml_rm);
    if (bcdat->rna_mols == NULL){
        free(bcdat);
        err_msg(-1, 0, "bc_data_init: %s", strerror(errno));
        return(NULL);
    }
    bcdat->atac_prs = kh_init(kh_arp);
    if (bcdat->atac_prs == NULL){
        free(bcdat);
        err_msg(-1, 0, "bc_data_init: %s", strerror(errno));
        return(NULL);
    }
    bcdat->atac_dups = kh_init(kh_ad);
    if (bcdat->atac_dups == NULL){
        free(bcdat);
        err_msg(-1, 0, "bc_data_init: %s", strerror(errno));
        return(NULL);
    }
    bcdat->atac_frgs = ml_alloc(ml_af);
    if (bcdat->atac_frgs == NULL){
        free(bcdat);
        err_msg(-1, 0, "bc_data_init: %s", strerror(errno));
        return(NULL);
    }
    bcdat->bc_stats = NULL;

    int err;
    if ((err = pthread_mutex_init(&bcdat->bc_lock, NULL)) != 0){
        err_msg(-1, 0, "bc_data_init: failed to initialize mutex %i", err);
        return(NULL);
    }

    return(bcdat);
}

void bc_data_free_rna_dupls(bc_data_t *bcdat) {
    if (bcdat == NULL) return;
    if (bcdat->rna_dups != NULL) {
        khash_t(kh_rd) *dups = bcdat->rna_dups;
        khint_t k;
        for (k = kh_begin(dups); k != kh_end(dups); ++k) {
            if (!kh_exist(dups, k)) continue;
            rna_dups_t *rd = &kh_key(dups, k);
            rna_dups_free(rd);
        }
        kh_destroy(kh_rd, bcdat->rna_dups);
        bcdat->rna_dups = NULL;
    }
}

void bc_data_free_rna_mlcls(bc_data_t *bcdat) {
    if (bcdat == NULL) return;
    if (bcdat->rna_mols != NULL) {
        ml_node_t(ml_rm) *node = ml_begin(bcdat->rna_mols);
        for (; node != NULL; node = ml_node_next(node)) {
            rna_mol_t *mol = &ml_node_val(node);
            rna_mol_free(mol);
        }
        ml_dstry(ml_rm, bcdat->rna_mols);
        bcdat->rna_mols = NULL;
    }
}

void bc_data_free_atac_pairs(bc_data_t *bcdat){
    if (bcdat == NULL) return;
    if (bcdat->atac_prs != NULL) {
        khash_t(kh_arp) *pairs = bcdat->atac_prs;
        khint_t pair;
        for (pair = kh_begin(pairs); pair != kh_end(pairs); ++pair) {
            if (!kh_exist(pairs, pair)) continue;
            atac_rd_pair_t *p = &kh_key(pairs, pair);
            atac_rd_pair_free(p);
        }
        kh_destroy(kh_arp, bcdat->atac_prs);
        bcdat->atac_prs = NULL;
    }
}

void bc_data_free_atac_dupls(bc_data_t *bcdat){
    if (bcdat == NULL) return;
    if (bcdat->atac_dups != NULL) {
        khash_t(kh_ad) *dups = bcdat->atac_dups;
        khint_t dup;
        for (dup = kh_begin(dups); dup != kh_end(dups); ++dup) {
            if (!kh_exist(dups, dup)) continue;
            atac_dups_t *ad = &kh_key(dups, dup);
            atac_dups_dstry(ad);
        }
        kh_destroy(kh_ad, bcdat->atac_dups);
        bcdat->atac_dups = NULL;
    }
}

void bc_data_free_atac_frags(bc_data_t *bcdat) {
    if (bcdat == NULL) return;
    if (bcdat->atac_frgs != NULL) {
        ml_node_t(ml_af) *node = ml_begin(bcdat->atac_frgs);
        for (; node != NULL; node = ml_node_next(node)) {
            atac_frag_t *af = &ml_node_val(node);
            atac_frag_free(af);
        }
        ml_dstry(ml_af, bcdat->atac_frgs);
        bcdat->atac_frgs = NULL;
    }
}

void bc_data_free_reads(bc_data_t *bcdat){
    if (bcdat == NULL) return;

    bc_data_free_rna_dupls(bcdat);

    bc_data_free_rna_mlcls(bcdat);
    
    bc_data_free_atac_pairs(bcdat);

    bc_data_free_atac_dupls(bcdat);

    bc_data_free_atac_frags(bcdat);
}

void bc_data_dstry(bc_data_t *bcdat){
    if (bcdat == NULL) return;

    bc_data_free_reads(bcdat);

    bc_stats_dstry(bcdat->bc_stats);

    pthread_mutex_destroy(&bcdat->bc_lock);

    free(bcdat);
}

/*******************************************************************************
 * RNA
 ******************************************************************************/

int bc_data_rna_add_read(bc_data_t *bcdat, const rna_read1_t *r, umishort umih){
    // check input
    if (bcdat == NULL || r == NULL)
        return err_msg(-1, 0, "bc_data_rna_add_read: argument is null");

    khash_t(kh_rd) *dups = bcdat->rna_dups;
    if (dups == NULL)
        return err_msg(-1, 0, "bc_data_rna_add_read: bcdat->rna_dups is null");

    // check if read is uninitialized
    if (r->loc.rid == -1)
        return err_msg(-1, 0, "bc_data_rna_add_read: read region is null");

    // first, add an empty rna_dup_t to rna_dups
    rna_dups_t rna_dup;
    if (rna_dups_init(&rna_dup) < 0)
        return err_msg(-1, 0, "bc_data_rna_add_read: failed to initialize rna_dups");
    rna_dup.umi = umih;

    int ret;
    khint_t k = kh_get(kh_rd, dups, rna_dup);
    if (k == kh_end(dups)) { // if not found
        k = kh_put(kh_rd, dups, rna_dup, &ret);
        if (k == kh_end(dups) || ret < 0) {
            rna_dups_free(&rna_dup);
            return err_msg(-1, 0, "bc_data_rna_add_read: failed to add read to rna_dups");
        }
    } else {
        rna_dups_free(&rna_dup);
    }
    rna_dups_t *p = &kh_key(dups, k);
    if (p == NULL) {
        rna_dups_free(&rna_dup);
        return err_msg(-1, 0, "bc_data_rna_add_read: failed to retrieve rna_dups");
    }

    if (rna_dups_add_read(p, r) < 0)
        return err_msg(-1, 0, "bc_data_rna_add_read: failed to add read to rna_dups");

    return 0;
}

int bc_data_rna_dedup(bc_data_t *bcdat){
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_rna_dedup: bcdat is null");

    khash_t(kh_rd) *dups = bcdat->rna_dups;
    if (dups == NULL)
        return err_msg(-1, 0, "bc_data_rna_dedup: bcdat->rna_dups is null");

    if (kh_size(dups) == 0)
        return 0;

    khint_t k;
    for (k = kh_begin(dups); k != kh_end(dups); ++k) {
        if (!kh_exist(dups, k)) continue;
        rna_dups_t *dupls = &kh_key(dups, k);
        if (dupls == NULL)
            return err_msg(-1, 0, "bc_data_rna_dedup: read is null");

        int ret;
        rna_mol_t *mol = rna_dups_dedup(dupls, &ret);
        rna_dups_free(dupls);
        if (ret < 0)
            return err_msg(-1, 0, "bc_data_rna_dedup: failed to deduplicate reads");

        if (mol == NULL)
            continue;

        ret = ml_insert(ml_rm, bcdat->rna_mols, *mol, 0, 0);
        if (ret < 0) {
            rna_mol_dstry(mol);
            return err_msg(-1, 0, "bc_data_rna_dedup: failed to add rna mol to bcdat");
        }
        free(mol);
    }

    kh_destroy(kh_rd, bcdat->rna_dups);
    bcdat->rna_dups = NULL;

    return 0;
}

int bc_data_rna_var_call(bc_data_t *bcdat, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (bcdat == NULL || gv == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_rna_var_call: argument is null");

    int n_add = 0;

    ml_t(ml_rm) *mols = bcdat->rna_mols;
    if (mols == NULL)
        return err_msg(-1, 0, "bc_data_rna_var_call: bcdat->rna_mols is null");

    ml_node_t(ml_rm) *node = ml_begin(mols);
    for (; node != NULL; node = ml_node_next(node)) {
        rna_mol_t *mol = &ml_node_val(node);
        if (mol == NULL)
            return err_msg(-1, 0, "bc_data_rna_var_call: mol is null");
        int ret = rna_mol_var_call(mol, gv, cmap, min_qual);
        if (ret < 0)
            return err_msg(-1, 0, "bc_data_rna_var_call: failed to call variant");
        n_add += ret;
    }

    return(n_add);
}

/*******************************************************************************
 * ATAC
 ******************************************************************************/

int bc_data_atac_add_read(bc_data_t *bcdat, const atac_read1_t *ar, qshort qname){
    if (bcdat == NULL || ar == NULL)
        return err_msg(-1, 0, "bc_data_atac_add_read: argument is null");

    // add atac pair struct to bcdat->atac_prs
    atac_rd_pair_t atac_pair;
    atac_rd_pair_set0(&atac_pair);
    atac_pair.n_reads = 1;
    atac_pair.qname = qname;

    khash_t(kh_arp) *pairs = bcdat->atac_prs;
    if (pairs == NULL)
        return err_msg(-1, 0, "bc_data_atac_add_read: bcdat->atac_prs is null");

    khint_t kp;
    kp = kh_get(kh_arp, pairs, atac_pair);

    // if pair isn't present, add it
    int ret = 0;
    if (kp == kh_end(pairs)){
        kp = kh_put(kh_arp, pairs, atac_pair, &ret);
        if (ret < 0)
            return err_msg(-1, 0, "bc_data_atac_add_read: failed to add read to bcdat");
    }
    atac_rd_pair_t *p = &kh_key(pairs, kp);
    if (p == NULL)
        return err_msg(-1, 0, "bc_data_atac_add_read: unexpected error occurred retrieving atac pair");

    // add read1 to atac pair struct
    if (atac_rd_pair_add_read(p, ar) < 0)
        return err_msg(-1, 0, "bc_data_atac_add_read: failed to add read to atac pair");

    // if we form a complete read pair, add to dups
    if (p->s == 2) {
        if (bc_data_atac_add_pair(bcdat, p) < 0)
            return err_msg(-1, 0, "bc_data_atac_add_read: failed to add pair to dups");
        atac_rd_pair_free(p);
        // p->qname = qname; // important for keeping hash valid
        kh_del(kh_arp, pairs, kp);
    }

    return(0);
}

int bc_data_atac_add_pair(bc_data_t *bcdat, const atac_rd_pair_t *ap) {
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_atac_add_pair: argument is null");
    if (ap == NULL)
        return err_msg(-1, 0, "bc_data_atac_add_pair: read pair is null");
    if (ap->s == 0)
        return err_msg(-1, 0, "bc_data_atac_add_pair: read pair is empty");
    if (ap->s != 2)
        return 0;

    int is_chim = atac_rd_pair_is_chimeric(ap);
    if (is_chim < 0)
        return err_msg(-1, 0, "bc_data_atac_add_pair: failed to check if pair is chimeric");
    if (is_chim)
        return 0;

    atac_dups_t dupls;
    atac_dups_init(&dupls);
    dupls.reg = get_reg_pair(ap->r1.loc, ap->r2.loc); // new key for dups

    khash_t(kh_ad) *dups = bcdat->atac_dups;
    if (dups == NULL)
        return err_msg(-1, 0, "bc_data_atac_add_pair: bcdat->atac_dups is null");

    khint_t kp;
    kp = kh_get(kh_ad, dups, dupls);

    // if duplicate not found, add it
    int ret = 0;
    if (kp == kh_end(dups)){
        kp = kh_put(kh_ad, dups, dupls, &ret);
        if (ret < 0)
            return err_msg(-1, 0, "bc_data_atac_add_pair: failed to add dups to bcdat");
    }
    atac_dups_t *ptr = &kh_key(dups, kp);
    if (ptr == NULL)
        return err_msg(-1, 0, "bc_data_atac_add_pair: failed to retrieve dups");

    if (atac_dups_add_pair(ptr, ap) < 0)
        return err_msg(-1, 0, "bc_data_atac_add_pair: failed to add read to dups");

    return 1;
}

int bc_data_atac_dedup(bc_data_t *bcdat){
    if (bcdat == NULL)
        return err_msg(-1, 0, "bc_data_atac_dedup: bcdat is null");

    if (bcdat->atac_dups == NULL)
        return err_msg(-1, 0, "bc_data_atac_dedup: bcdat->atac_dups is null");

    if (kh_size(bcdat->atac_dups) == 0)
        return 0;
    
    khash_t(kh_ad) *bc_dups = bcdat->atac_dups;

    khint_t kp;
    for (kp = kh_begin(bc_dups); kp != kh_end(bc_dups); ++kp) {
        if (!kh_exist(bc_dups, kp)) continue;
        atac_dups_t *dups = &kh_key(bc_dups, kp);
        if (dups == NULL)
            return err_msg(-1, 0, "bc_data_atac_dedup: read pair is null");

        int ret;
        atac_frag_t *frag = atac_dups_dedup(dups, &ret);
        atac_dups_free(dups); // free space after deduplication
        if (ret < 0)
            return err_msg(-1, 0, "bc_data_atac_dedup: failed to deduplicate reads");

        if (frag == NULL)
            continue;

        // add to frags
        ret = ml_insert(ml_af, bcdat->atac_frgs, *frag, 0, 0);
        if (ret < 0) {
            atac_frag_dstry(frag);
            return err_msg(-1, 0, "bc_data_atac_dedup: failed to add frag to bcdat");
        }
        free(frag);
    }

    // free hash
    kh_destroy(kh_ad, bc_dups);
    bcdat->atac_dups = NULL;

    return 0;
}

int bc_data_atac_var_call(bc_data_t *bcdat, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (bcdat == NULL || gv == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_atac_var_call: argument is null");

    ml_t(ml_af) *frags = bcdat->atac_frgs;
    if (frags == NULL)
        return err_msg(-1, 0, "bc_data_atac_var_call: bcdat->atac_frgs is null");

    int n_add = 0;

    ml_node_t(ml_af) *node = ml_begin(frags);
    for (; node != NULL; node = ml_node_next(node)) {
        atac_frag_t *frag = &ml_node_val(node);
        if (frag == NULL)
            return err_msg(-1, 0, "bc_data_atac_var_call: frag is null");
        int ret = atac_frag_var_call(frag, gv, cmap, min_qual);
        if (ret < 0)
            return err_msg(-1, 0, "bc_data_atac_var_call: failed to call variant");
        n_add += ret;
    }

    return n_add;
}

int bc_data_atac_peak_call(bc_data_t *bcdat, iregs_t *pks, str_map *cmap){
    if (bcdat == NULL || pks == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_data_atac_peak_call: argument is null");

    ml_t(ml_af) *frags = bcdat->atac_frgs;
    if (frags == NULL)
        return err_msg(-1, 0, "bc_data_atac_peak_call: bcdat->atac_frgs is null");

    int n_add = 0;

    ml_node_t(ml_af) *node = ml_begin(frags);
    for (; node != NULL; node = ml_node_next(node)) {
        atac_frag_t *frag = &ml_node_val(node);
        g_reg_pair reg = frag->reg;
        if (frag == NULL)
            return err_msg(-1, 0, "bc_data_atac_peak_call: frag is null");
        int np = atac_frag_peak_call(frag, reg, pks, cmap);
        if (np < 0)
            return err_msg(-1, 0, "bc_data_atac_peak_call: failed to call peaks");
        n_add += np;
    }

    return n_add;
}

/*******************************************************************************
 * bam_data_t
 ******************************************************************************/

bam_data_t *bam_data_init(){
    bam_data_t *b = (bam_data_t *)calloc(1, sizeof(bam_data_t));
    if (b == NULL){
        err_msg(-1, 0, "bam_data_init: %s", strerror(errno));
        return NULL;
    }

    b->bc_data = kh_init(kh_bc_dat);
    if (b->bc_data == NULL) {
        free(b);
        err_msg(-1, 0, "bam_data_init: %s", strerror(errno));
        return NULL;
    }
    b->has_rna = 0;
    b->has_atac = 0;
    b->has_stats = 0;

    b->bcs = NULL;
    b->bcs = init_str_map();
    if (b->bcs == NULL) {
        kh_destroy(kh_bc_dat, b->bc_data);
        free(b);
        err_msg(-1, 0, "bam_data_init: %s", strerror(errno));
        return NULL;
    }

    int err;
    if ((err = pthread_mutex_init(&b->bam_lock, NULL)) != 0){
        err_msg(-1, 0, "bam_data_init: failed to initialize mutex errno=%i", err);
        kh_destroy(kh_bc_dat, b->bc_data);
        destroy_str_map(b->bcs);
        free(b);
        return NULL;
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

    bam_data_free_bcs(bam_data);

    destroy_str_map(bam_data->bcs);

    pthread_mutex_destroy(&bam_data->bam_lock);

    free(bam_data);
}

void bam_data_atac_free_pairs(bam_data_t *bam_data){
    if (bam_data == NULL)
        return;
    if (bam_data->bc_data == NULL)
        return;
    khash_t(kh_bc_dat) *bcs_hash = bam_data->bc_data;
    if (bcs_hash == NULL)
        return;
    khint_t kbc;
    for (kbc = kh_begin(bcs_hash); kbc != kh_end(bcs_hash); ++kbc){
        if (!kh_exist(bcs_hash, kbc)) continue;
        bc_data_t *bc_dat = kh_val(bcs_hash, kbc);
        bc_data_free_atac_pairs(bc_dat);
    }
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
        return 0;
    }

    // fill in barcodes from data
    // barcode is only added if not present.
    int found;
    khint_t k;
    for (k = kh_begin(bam_data->bc_data); k != kh_end(bam_data->bc_data); ++k){
        if (!kh_exist(bam_data->bc_data, k)) continue;
        char *bc_key = kh_key(bam_data->bc_data, k);
        if (bc_key == NULL) continue;
        if (add2str_map(bam_data->bcs, bc_key, &found) < 0) {
            destroy_str_map(bam_data->bcs);
            bam_data->bcs = NULL;
            return err_msg(-1, 0, "bam_data_fill_bcs: failed to add barcode");
        }
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

    if (bam_data->bc_data == NULL)
        return err_msg(-1, 0, "bam_data_rna_add_read: bam_data->bc_data is null");

    if (strlen(bc) < 5)
        fprintf(stderr, "strange BC=%s\n", bc);

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
            free(bc_cpy);
            pthread_mutex_unlock(&bam_data->bam_lock);
            return err_msg(-1, 0, "bam_data_rna_add_read: failed to add read to bcs");
        }

        if ( (kh_val(bcs_hash, kbc) = bc_data_init()) == NULL ){
            kh_del(kh_bc_dat, bcs_hash, kbc);
            free(bc_cpy);
            pthread_mutex_unlock(&bam_data->bam_lock);
            return err_msg(-1, 0, "bam_data_rna_add_read: failed to add barcode to bcs");
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

    bam_data->has_rna = 1;

    return 0;
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
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_rna_dedup: barcode data bc_dat is null");

        if (bc_data_rna_dedup(bc_dat) < 0)
            return err_msg(-1, 0, "bam_data_rna_dedup: failed to deduplicate reads");
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
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_rna_var_call: barcode data bc_dat is null");

        if (bc_data_rna_var_call(bc_dat, gv, cmap, min_qual) < 0)
            return err_msg(-1, 0, "bam_data_rna_var_call: failed to call variants");
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
        if (bc_cpy == NULL) {
            pthread_mutex_unlock(&bam_data->bam_lock);
            return err_msg(-1, 0, "bam_data_atac_add_read: %s", strerror(errno));
        }

        kbc = kh_put(kh_bc_dat, bcs_hash, bc_cpy, &ret);
        if (ret < 0) {
            free(bc_cpy);
            pthread_mutex_unlock(&bam_data->bam_lock);
            return err_msg(-1, 0, "bam_data_atac_add_read: failed to add barcode to bcs");
        }

        if ( (kh_val(bcs_hash, kbc) = bc_data_init()) == NULL ) {
            kh_del(kh_bc_dat, bcs_hash, kbc);  // Remove the entry from the hash table
            free(bc_cpy);
            pthread_mutex_unlock(&bam_data->bam_lock);
            return err_msg(-1, 0, "bam_data_atac_add_read: failed to add barcode to bcs");
        }
    }
    bc_data_t *bc_dat = kh_val(bcs_hash, kbc);

    if (pthread_mutex_unlock(&bam_data->bam_lock) != 0)
        return err_msg(-1, 0, "bam_data_atac_add_read: failed bam mutex unlock");

    if (pthread_mutex_lock(&bc_dat->bc_lock) != 0)
        return err_msg(-1, 0, "bam_data_atac_add_read: failed bc mutex lock");

    ret = bc_data_atac_add_read(bc_dat, r, qname);

    if (pthread_mutex_unlock(&bc_dat->bc_lock) != 0)
        return err_msg(-1, 0, "bam_data_atac_add_read: failed bc mutex unlock");

    if (ret < 0)
        return err_msg(-1, 0, "bam_data_atac_add_read: failed to add read to bam");

    bam_data->has_atac = 1;

    return 0;
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
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_atac_dedup: barcode data bc_dat is null");

        if (bc_data_atac_dedup(bc_dat) < 0)
            return err_msg(-1, 0, "bam_data_atac_dedup: failed to deduplicate reads");
    }

    return 0;
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
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_atac_var_call: barcode data bc_dat is null");

        if (bc_data_atac_var_call(bc_dat, gv, cmap, min_qual) < 0)
            return err_msg(-1, 0, "bam_data_atac_var_call: failed to call variants");
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
        if (bc_dat == NULL)
            return err_msg(-1, 0, "bam_data_atac_peak_call: barcode data bc_dat is null");

        if (bc_data_atac_peak_call(bc_dat, pks, cmap) < 0)
            return err_msg(-1, 0, "bam_data_atac_peak_call: failed to call peaks");
    }

    return 0;
}

/*******************************************************************************
 * bc_stats
 ******************************************************************************/

int bc_data_rna_fill_stats(bc_stats_t *bcc, const sam_hdr_t *rna_hdr, bc_data_t *bc_dat){
    if (bcc == NULL || rna_hdr == NULL || bc_dat == NULL)
        return err_msg(-1, 0, "bc_data_rna_fill: arguments are null");

    if (bc_dat->rna_mols == NULL)
        return err_msg(-1, 0, "bc_data_rna_fill_stats: bc_dat->rna_mols is null");

    float n_mt = 0;
    uint32_t in_gene = 0;

    ml_t(ml_rm) *mols = bc_dat->rna_mols;
    ml_node_t(ml_rm) *node = ml_begin(mols);
    for (; node; node = ml_node_next(node)) {
        rna_mol_t *mol = &ml_node_val(node);
        if (mol == NULL)
            return err_msg(-1, 0, "bc_data_rna_fill_stats: mol is null");

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

    ml_t(ml_af) *frags = bc_dat->atac_frgs;
    if (frags == NULL)
        return err_msg(-1, 0, "bc_data_atac_fill_stats: bc_dat->atac_frgs is null");

    float n_mt = 0;
    uint32_t in_pk = 0;

    ml_node_t(ml_af) *node = ml_begin(frags);
    for (; node; node = ml_node_next(node)) {
        atac_frag_t *frag = &ml_node_val(node);
        g_reg_pair reg = frag->reg;
        if (frag == NULL)
            return err_msg(-1, 0, "bc_data_atac_fill_stats: frag is null");

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
        ml_node_t(seq_vac_l) *vn = ml_begin(&frag->vl);
        for (; vn; vn = ml_node_next(vn)){
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
        if ( (is_mt = chr_is_mt(chr)) < 0 )
            return(-1);
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
            if (bc_data_rna_fill_stats(bcc, rna_hdr, bc_dat) < 0) {
                bc_stats_free(bcc);
                return err_msg(-1, 0, "bam_data_fill_stats: failed to fill rna stats");
            }
        }
        if (bam_data->has_atac){
            if (bc_data_atac_fill_stats(bcc, atac_hdr, bc_dat) < 0) {
                bc_stats_free(bcc);
                return err_msg(-1, 0, "bam_data_fill_stats: failed to fill atac stats");
            }
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
        free(bc_array);
        *ret = err_msg(-1, 0, "bam_data_count_of_n: top_n=%zu must be between ", 
                ">= 0 and <= %i", top_n, n_bcs);
        return(0);
    }

    int bci = 0;
    for (k = kh_begin(bcs_hash); k != kh_end(bcs_hash); ++k){
        if (!kh_exist(bcs_hash, k)) continue;

        char *bc_key = kh_key(bcs_hash, k);

        if (bc_key == NULL) {
            *ret = err_msg(-1, 0, "bam_data_count_of_n: 'bc_key' is null");
            free(bc_array);
            return 0;
        }

        bc_data_t *bc_dat = kh_val(bcs_hash, k);

        if (bc_dat == NULL) {
            *ret = err_msg(-1, 0, "bam_data_count_of_n: 'bc_dat' is null");
            free(bc_array);
            return 0;
        }

        if (bc_dat->bc_stats == NULL) {
            *ret = err_msg(-1, 0, "bam_data_count_of_n: bc 'bc_stats' is null");
            free(bc_array);
            return 0;
        }

        bc_array[bci].bc = bc_key;
        bc_array[bci].count = bc_dat->bc_stats->counts;
        ++bci;
    }
    if (bci != n_bcs){
        *ret = err_msg(-1, 0, "bam_data_count_of_n: bci=%i != n_bcs=%i", bci, n_bcs);
        free(bc_array);
        return(0);
    }
    qsort(bc_array, n_bcs, sizeof(bc_count_t), cmp_bc_count_rev);

    uint32_t n_c = bc_array[top_n].count;
    free(bc_array);
    return(n_c);
}
