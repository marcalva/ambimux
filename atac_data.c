
#include "atac_data.h"
#include <stdlib.h>
#include <string.h>
#include "str_util.h"
#include "variants.h"
#include "region.h"
#include "counts.h"
#include "overlap.h"
#include "htslib/hts.h"

/*******************************************************************************
 * atac_read1_t
 ******************************************************************************/

void atac_read_set0(atac_read1_t *ar){
    if (ar == NULL) return;

    init_g_region(&ar->loc);
    init_seq_blist(&ar->bases);
}

atac_read1_t *atac_read_init(){
    atac_read1_t *ar = (atac_read1_t *)calloc(1, sizeof(atac_read1_t));

    if (ar == NULL){
        err_msg(-1, 0, "atac_read_init: %s", strerror(errno));
        return NULL;
    }

    atac_read_set0(ar);

    return(ar);
}

void atac_read_free(atac_read1_t *ar){
    if (ar == NULL) return;

    free_seq_blist(&ar->bases);
}

void atac_read_dstry(atac_read1_t *ar){
    if (ar == NULL) return;
    atac_read_free(ar);
    free(ar);
}

atac_read1_t *atac_read_cpy(const atac_read1_t *r, int *ret){
    *ret = 0;
    if (r == NULL) return(NULL);

    atac_read1_t *cpy = atac_read_init();
    if (cpy == NULL){
        *ret = -1;
        return(NULL);
    }

    cpy->loc = r->loc;

    int cret;
    seq_blist_t *b_cpy = copy_seq_blist(&r->bases, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    cpy->bases = *b_cpy;
    free(b_cpy);

    return(cpy);
}

int atac_read_equal(atac_read1_t r1, atac_read1_t r2){

    /* compare region */
    int rcmp = regioncmp(r1.loc, r2.loc);
    if (rcmp != 0)
        return(0);

    /* compare bases, ignore qual */
    int bcmp = seq_blist_equal(r1.bases, r2.bases, 0);
    if (bcmp != 1)
        return(0);
    return(1);
}

int atac_read1_add_base(atac_read1_t *r, const seq_base_t *base){
    if (r == NULL || base == NULL)
        return err_msg(-1, 0, "atac_read1_add_base: arguments are NULL");

    if (blist_add_base(&r->bases, base, 1, 1) < 0) return(-1);

    return(0);
}

/*******************************************************************************
 * atac_rd_pair_t
 ******************************************************************************/

void atac_rd_pair_set0(atac_rd_pair_t *rp){
    if (rp == NULL) return;
    atac_read_set0(&rp->r1);
    atac_read_set0(&rp->r2);
    rp->s = 0;
}

atac_rd_pair_t *atac_rd_pair_init(){
    atac_rd_pair_t *rp = (atac_rd_pair_t *)calloc(1, sizeof(atac_rd_pair_t));
    if (rp == NULL){
        err_msg(-1, 0, "atac_rd_pair_init: %s", strerror(errno));
        return(NULL);
    }

    atac_rd_pair_set0(rp);
    
    return(rp);
}

void atac_rd_pair_free(atac_rd_pair_t *rp){
    if (rp == NULL) return;

    atac_read_free(&rp->r1);
    atac_read_free(&rp->r2);
}

void atac_rd_pair_dstry(atac_rd_pair_t *rp){
    if (rp == NULL) return;

    atac_rd_pair_free(rp);
    free(rp);
}

atac_rd_pair_t *atac_rd_pair_cpy(const atac_rd_pair_t *rp, int *ret){
    *ret = 0;
    if (rp == NULL) return(NULL);

    atac_rd_pair_t *cp = (atac_rd_pair_t *)calloc(1, sizeof(atac_rd_pair_t));
    if (cp == NULL){
        *ret = err_msg(-1, 0, "atac_rd_pair_cpy: %s", strerror(errno));
        return(NULL);
    }

    atac_read1_t *rtmp = NULL;

    int cret = 0;
    rtmp = atac_read_cpy(&rp->r1, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    cp->r1 = *rtmp;
    free(rtmp);
    rtmp = atac_read_cpy(&rp->r2, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    cp->r2 = *rtmp;
    free(rtmp);
    cp->s = rp->s;
    return(cp);
}

int atac_rd_pair_add_read(atac_rd_pair_t *rp, const atac_read1_t *ar){
    if (rp == NULL || ar == NULL)
        return err_msg(-1, 0, "atac_tpl_add_read: arguments must not be NULL");

    int cret = 0;
    atac_read1_t *rtmp = NULL;

    if (rp->s == 0){
        rtmp = atac_read_cpy(ar, &cret);
        rp->r1 = *rtmp;
        if (cret < 0) return(-1);
        free(rtmp);
    } else if (rp->s == 1){
        rtmp = atac_read_cpy(ar, &cret);
        rp->r2 = *rtmp;
        if (cret < 0) return(-1);
        free(rtmp);

        /* switch order if r2 comes before r1 */
        int rcmp = regioncmp(rp->r2.loc, rp->r1.loc);
        if (rcmp < 0){
            atac_read1_t tmp = rp->r2;
            rp->r2 = rp->r1;
            rp->r1 = tmp;
        }
    } else {
        return err_msg(-1, 0, "atac_rd_pair_add_read: trying to add read to full pair "
                "(query name was found three times)");
    }
    ++rp->s;
    return(0);
}

int atac_rd_pair_equal(const atac_rd_pair_t *rp1, const atac_rd_pair_t *rp2){
    if (rp1 == NULL || rp2 == NULL)
        return err_msg(-1, 0, "atac_rd_pair_equal: arguments cannot be NULL");

    if (rp1->s != rp2->s)
        return(0);

    int rcmp1 = atac_read_equal(rp1->r1, rp2->r1); // compare read 1
    if (rcmp1 != 1) return(0);

    int rcmp2 = atac_read_equal(rp1->r2, rp2->r2); // compare read 1
    if (rcmp2 != 1) return(0);

    return(1);
}

int atac_rd_pair_match_qual(atac_rd_pair_t *rp, const atac_rd_pair_t *cmp){
    if (rp == NULL || cmp == NULL)
        return err_msg(-1, 0, "atac_rd_pair_match_qual: arguments are NULL");
    if (rp->s != cmp->s)
        return err_msg(-1, 0, "atac_rd_pair_match_qual: "
                "number of reads don't match");
    if (rp->s == 0)
        return(0);

    if (seq_blist_match_qual(&rp->r1.bases, &cmp->r1.bases) < 0)
        return(-1);
    if (rp->s == 2){
        if (seq_blist_match_qual(&rp->r2.bases, &cmp->r2.bases) < 0)
            return(-1);
    }

    return(0);
}

/*******************************************************************************
 * atac_dups_t
 ******************************************************************************/

atac_dups_t *atac_dups_init(){
    atac_dups_t *d = (atac_dups_t *)calloc(1, sizeof(atac_dups_t));
    if (d == NULL){
        err_msg(-1, 0, "atac_dups_init: %s", strerror(errno));
        return(NULL);
    }
    d->size = 0;
    d->m = 1;
    d->dups = calloc(d->m, sizeof(struct _atac_dup1_t));
    if (d->dups == NULL){
        err_msg(-1, 0, "atac_dups_init: %s", strerror(errno));
        return(NULL);
    }
    return(d);
}

void atac_dups_free(atac_dups_t *d){
    if (d == NULL) return;
    int i;
    for (i = 0; i < d->size; ++i)
        atac_rd_pair_free(&d->dups[i].rd); // free underlying reads
    free(d->dups);
}

void atac_dups_dstry(atac_dups_t *d){
    if (d == NULL) return;
    atac_dups_free(d);
    free(d);
}

int atac_dups_add_pair(atac_dups_t *d, const atac_rd_pair_t *rp){
    if (d == NULL || rp == NULL)
        return err_msg(-1, 0, "atac_dups_add_pair: arguments must not be NULL");

    int i, ret = 0, tcmp = 0;
    for (i = 0; i < d->size; ++i){
        tcmp = atac_rd_pair_equal(&d->dups[i].rd, rp);
        if (tcmp < 0)
            return(-1);
        if (tcmp == 1){
            if (atac_rd_pair_match_qual(&d->dups[i].rd, rp) < 0)
                return(-1);
            ++d->dups[i].n_rd;;
            return(0);
        }
    }
    if (d->size >= d->m){
        d->m = d->size + 1;
        d->dups = realloc(d->dups, d->m * sizeof(struct _atac_dup1_t));
        if (d->dups == NULL)
            return err_msg(-1, 0, "atac_dups_add_pair: %s", strerror(errno));
    }

    // copy read pair
    atac_rd_pair_t *tmp = atac_rd_pair_cpy(rp, &ret);
    if (ret < 0)
        return(-1);

    d->dups[d->size].rd = *tmp;
    d->dups[d->size].n_rd = 1;
    free(tmp);

    ++d->size;

    return(0);
}

/*******************************************************************************
 * atac_frag_t
 ******************************************************************************/

atac_frag_t *atac_frag_init(){
    atac_frag_t *f = (atac_frag_t *)calloc(1, sizeof(atac_frag_t));
    if (f == NULL){
        err_msg(-1, 0, "atac_frag_init: %s", strerror(errno));
        return(NULL);
    }
    init_seq_blist(&f->bases);
    init_vacs(&f->vacs);
    f->pks.ix = NULL;
    f->pks.n = 0;
    f->pks.m = 0;
    f->s = 0;
    return(f);
}

void atac_frag_dstry(atac_frag_t *f){
    if (f == NULL) return;

    free_seq_blist(&f->bases);
    free_vacs(&f->vacs);
    free(f->pks.ix);
    free(f);
}

atac_frag_t *atac_dups_dedup(const atac_dups_t *d, int *ret){
    *ret = 0;
    if (d == NULL){
        err_msg(-1, 0, "atac_dups_dedup: d is NULL");
        *ret = -1;
        return(NULL);
    }

    // check for bugs
    if (d->size == 0){
        err_msg(-1, 0, "atac_dups_dedup: initialized atac_dups_t has no "
                "reads/duplicates, there is a bug");
        *ret = -1;
        return(NULL);
    }

    int ix_best;
    size_t max_c = 0; // store read count
    size_t rp_w_max = 0; // number of read pairs with max_c read counts

    int i;
    for (i = 0; i < d->size; ++i){
        // check for bugs
        if (d->dups[i].n_rd == 0){
            err_msg(-1, 0, "atac_dups_dedup: a dup read pair has 0 "
                    "supporting reads, there is a bug");
            *ret = -1;
            return(NULL);
        }

        if (d->dups[i].n_rd == max_c){
            ++rp_w_max;
        } else if (d->dups[i].n_rd > max_c){
            max_c = d->dups[i].n_rd;
            ix_best = i;
            rp_w_max = 1;
        }
    }

    if (max_c == 0 || rp_w_max == 0){
        err_msg(-1, 0, "atac_dups_dedup: no best duplicate found "
                ", there is a bug");
        *ret = -1;
        return(NULL);
    }

    /* initialize fragment */
    atac_frag_t *f = atac_frag_init();
    atac_rd_pair_t rpb = d->dups[ix_best].rd;

    /* add region */
    // f->loc = get_reg_pair(rpb.r1->loc, rpb.r2->loc);

    int bret;
    seq_base_t *base;

    // check for bugs
    if (rpb.s < 2){
        err_msg(-1, 0, "atac_dups_dedup: dup read pair has less than 2 reads"
                ", there is a bug");
        *ret = -1;
        return(NULL);
    }

    /* add bases from read 1 */
    base = rpb.r1.bases.head;
    while (base){
        bret = blist_add_base(&f->bases, base, 1, 1);
        if (bret < 0){
            *ret = -1;
            return(NULL);
        }
        base = base->next;
    }

    /* add bases from read 2 */
    base = rpb.r2.bases.head;
    while (base){
        bret = blist_add_base(&f->bases, base, 1, 1);
        if (bret < 0){
            *ret = -1;
            return(NULL);
        }
        base = base->next;
    }

    /* add number of supporting reads */
    f->s = d->dups[ix_best].n_rd;

    return(f);
}

int atac_frag_var_call(atac_frag_t *f, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual){
    if (f == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "atac_frag_var_call: arguments must not be NULL");
    }

    int ret = seq_blist_call_var(&f->bases, &f->vacs, gv, cmap, min_qual);

    // free the bases
    free_seq_blist(&f->bases);

    return(ret);
}

int atac_frag_peak_call(atac_frag_t *f, g_reg_pair reg, iregs_t *pks, contig_map *cmap){
    if (f == NULL)
        return err_msg(-1, 0, "atac_frag_peak_call: f is NULL");
    if (pks == NULL)
        return err_msg(-1, 0, "atac_frag_peak_call: pks is NULL");

    if (reg.r1.rid != reg.r2.rid)
        return err_msg(-1, 0, "atac_frag_peak_call: chromosomes in given region pair don't matc");
    int rid = (int)reg.r1.rid;
    const char *chr = cm_ix_to_chr(cmap, rid);

    hts_pos_t beg = reg.r1.start;
    hts_pos_t end = reg.r2.end - 1;
    if (end < beg)
        return err_msg(-1, 0, "atac_frag_peak_call: end %" PRIhts_pos " < beg %" PRIhts_pos, end, beg);

    iregn_t overlaps = {NULL,0,0};
    int ret = iregs_overlap(pks, chr, beg, end, &overlaps);
    if (ret < 0) return(-1);

    // check if peak overlaps by at least one base
    // if so, add to f->pks
    int i;
    for (i = 0; i < overlaps.n; ++i){
        int ix = overlaps.ix[i];
        if (ix >= pks->n)
            return err_msg(-1, 0, "atac_frag_peak_call: %i index >= num of peaks %i", ix, pks->n);
        g_region pk_reg = pks->reg[ix];
        int64_t ovrlp;
        // check cut site 1 (+4 based on 10X)
        ovrlp = bp_overlap(beg+4, beg+5, '.', pk_reg.start, pk_reg.end, '.');
        if (ovrlp > 0 && iregn_add_ix(&f->pks, overlaps.ix[i]) < 0)
            return(-1);

        // check cut site 2 (-5 based on 10X)
        ovrlp = bp_overlap(end-5, end-4, '.', pk_reg.start, pk_reg.end, '.');
        if (ovrlp > 0 && iregn_add_ix(&f->pks, overlaps.ix[i]) < 0)
            return(-1);
    }

    free(overlaps.ix);

    return(f->pks.n);
}

/*******************************************************************************
 * bc_atac_t
 ******************************************************************************/

bc_atac_t *bc_atac_init(){
    bc_atac_t *bca = (bc_atac_t *)calloc(1, sizeof(bc_atac_t));
    if (bca == NULL){
        err_msg(-1, 0, "bc_atac_init: %s", strerror(errno));
        return(NULL);
    }
    bca->frags = kh_init(khaf);
    bca->dups = kh_init(khad);
    bca->pairs = kh_init(khap);
    if ( bca->frags == NULL || bca->dups == NULL || bca->pairs == NULL ){
        err_msg(-1, 0, "bc_atac_init: %s", strerror(errno));
        return(NULL);
    }
    return(bca);
}

void bc_atac_dstry_pairs(bc_atac_t *bca){
    if (bca == NULL || bca->pairs == NULL) return;

    khint_t k;
    for (k = kh_begin(bca->pairs); k != kh_end(bca->pairs); ++k){
        if (!kh_exist(bca->pairs, k)) continue;
        atac_rd_pair_t *rp = kh_val(bca->pairs, k);
        atac_rd_pair_dstry(rp);
    }
    kh_destroy(khap, bca->pairs);
    bca->pairs = NULL;
}

void bc_atac_dstry_dups(bc_atac_t *bca){
    if (bca == NULL || bca->dups == NULL) return;

    khint_t k;
    for (k = kh_begin(bca->dups); k != kh_end(bca->dups); ++k){
        if (!kh_exist(bca->dups, k)) continue;
        atac_dups_t *d = kh_val(bca->dups, k);
        atac_dups_dstry(d);
    }
    kh_destroy(khad, bca->dups);
    bca->dups = NULL;
}

void bc_atac_dstry_frags(bc_atac_t *bca){
    if (bca == NULL || bca->frags == NULL) return;


    khint_t k;
    for (k = kh_begin(bca->frags); k != kh_end(bca->frags); ++k){
        if (!kh_exist(bca->frags, k)) continue;
        atac_frag_t *f = kh_val(bca->frags, k);
        atac_frag_dstry(f);
    }
    kh_destroy(khaf, bca->frags);
    bca->frags = NULL;
}

void bc_atac_dstry(bc_atac_t *bca){
    if (!bca) return;

    bc_atac_dstry_frags(bca);
    bc_atac_dstry_dups(bca);
    bc_atac_dstry_pairs(bca);

    free(bca);
}

int bc_atac_add_read(bc_atac_t *bca, const atac_read1_t *ar, qshort qname){
    if (bca == NULL || ar == NULL)
        return err_msg(-1, 0, "bc_atac_add_read: arguments cannot be NULL");

    if (bca->pairs == NULL)
        return err_msg(-1, 0, "bc_atac_add_read: pairs is NULL, likely a bug");

    khash_t(khap) *pairs = bca->pairs;

    int ret;
    khint_t kp;
    atac_rd_pair_t *rp;
    kp = kh_get(khap, pairs, qname);

    // if read pair isn't present, allocate and add read pair t
    if (kp == kh_end(pairs)){
        kp = kh_put(khap, pairs, qname, &ret);
        if (ret < 0)
            return err_msg(-1, 0, "bc_atac_add_reads: failed to add qname to khap");
        
        if ( (kh_val(pairs, kp) = atac_rd_pair_init()) == NULL ) return(-1);
    }
    rp = kh_val(pairs, kp);

    // add read (copy) to read pair
    if ( atac_rd_pair_add_read(rp, ar) < 0 )
        return err_msg(-1, 0, "bc_atac_add_reads: failed to add read to pair");

    // If a read pair was formed, add to duplicates and destroy read.
    if (rp->s == 2){
        if (bc_atac_add_dup(bca, rp) < 0)
            return(-1);
        atac_rd_pair_dstry(rp);
        kh_del(khap, bca->pairs, kp);
        if ( kh_size(bca->pairs) > 4 && 
                (kh_n_buckets(bca->pairs) > (10 * kh_size(bca->pairs))) )
            kh_resize(khap, bca->pairs, kh_size(bca->pairs));
    }

    return(0);
}

int bc_atac_add_dup(bc_atac_t *bca, atac_rd_pair_t *rp){
    if (bca == NULL)
        return err_msg(-1, 0, "bc_atac_add_dup: arguments must not be NULL");
    if (rp == NULL)
        return(0);
    int kret;
    khint_t kd;
    g_reg_pair reg_pair = get_reg_pair(rp->r1.loc, rp->r2.loc);
    kd = kh_get(khad, bca->dups, reg_pair);
    if (kd == kh_end(bca->dups)){
        // add region key
        kd = kh_put(khad, bca->dups, reg_pair, &kret);// TODO: high heap alloc
        if (kret < 0){
            return err_msg(-1 , 0, "bc_atac_add_dup: failed to add "
                    "duplicate reg key to hash table");
        }
        // add allocated dups // TODO: high heap alloc
        if ( (kh_val(bca->dups, kd) = atac_dups_init()) == NULL) return(-1);
    }
    // add read pair to dups
    atac_dups_t *dups = kh_val(bca->dups, kd);
    if (atac_dups_add_pair(dups, rp) < 0)
        return(-1);

    return(0);
}

int bc_atac_form_dups(bc_atac_t *bca){
    if (!bca)
        return err_msg(-1, 0, "bc_atac_form_dups: arguments must not be NULL");

    khint_t kp;
    // loop through read pairs in barcode
    for (kp = kh_begin(bca->pairs); kp != kh_end(bca->pairs); ++kp){
        if (!kh_exist(bca->pairs, kp)) continue;
        atac_rd_pair_t *rp = kh_val(bca->pairs, kp);
        if (rp == NULL) return err_msg(-1, 0, "bc_atac_form_dups: read pair is NULL");
        if (rp->s != 2) continue;
        
        if (bc_atac_add_dup(bca, rp) < 0)
            return(-1);

        // destroy read pair after forming dup
        atac_rd_pair_dstry(rp);
        kh_del(khap, bca->pairs, kp);
    }
    return(0);
}

int bc_atac_dedup(bc_atac_t *bca){
    if (!bca)
        return err_msg(-1, 0, "bc_atac_dedup: arguments must not be NULL");

    if (!bca->dups)
        return err_msg(-1, 0, "bc_atac_dedup: no duplicates found");

    int kret;
    khint_t kd, kf;
    for (kd = kh_begin(bca->dups); kd != kh_end(bca->dups); ++kd){
        if (!kh_exist(bca->dups, kd)) continue;

        atac_dups_t *d = kh_val(bca->dups, kd);
        if (!d) return err_msg(-1, 0, "bc_atac_dedup: duplicate is NULL");
        g_reg_pair reg = kh_key(bca->dups, kd);

        // skip the duplicate if chimeric
        if (reg.r1.rid != reg.r2.rid){
            atac_dups_dstry(d);
            kh_del(khad, bca->dups, kd);
            continue;
        }

        // add the region key for the frag
        kf = kh_put(khaf, bca->frags, reg, &kret);// TODO: high heap alloc
        // should not be present
        if (kret == 0)
            return err_msg(-1, 0, "bc_atac_dedup: region already found, "
                    "this is likely because of a bug");
        if (kret < 0)
            return err_msg(-1, 0, "bc_atac_dedup: could not add key to frag table");

        // form de-duplicated frag
        int dret = 0;
        atac_frag_t *f = atac_dups_dedup(d, &dret);// TODO: high heap alloc
        if (dret < 0) return(-1);

        // add to hash
        kh_val(bca->frags, kf) = f;

        // destroy dups after forming frag
        atac_dups_dstry(d);
        kh_del(khad, bca->dups, kd);
    }
    kh_resize(khad, bca->dups, kh_size(bca->dups));

    return(0);
}

int bc_atac_var_call(bc_atac_t *bca, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual){
    if (bca == NULL)
        return err_msg(-1, 0, "bc_atac_var_call: bca is null");

    if (bca->frags == NULL)
        return err_msg(-1, 0, "bc_atac_var_call: bca frags is null");

    int n_add = 0;
    khint_t kf;
    for (kf = kh_begin(bca->frags); kf != kh_end(bca->frags); ++kf){
        if (!kh_exist(bca->frags, kf)) continue;

        atac_frag_t *f = kh_val(bca->frags, kf);
        if (!f)
            return err_msg(-1, 0, "bc_atac_var_call: frag is NULL");

        int a = atac_frag_var_call(f, gv, cmap, min_qual);
        if (a < 0)
            return(-1);
        n_add += a;
    }
    return n_add;
}

int bc_atac_peak_call(bc_atac_t *bca, iregs_t *pks, contig_map *cmap){
    if (bca == NULL || pks == NULL || cmap == NULL)
        return err_msg(-1, 0, "bc_atac_peak_call: arguments are null");

    if (bca->frags == NULL)
        return err_msg(-1, 0, "bc_atac_peak_call: bca frags is null");

    int n_add = 0;
    khint_t kf;
    for (kf = kh_begin(bca->frags); kf != kh_end(bca->frags); ++kf){
        if (!kh_exist(bca->frags, kf)) continue;

        atac_frag_t *f = kh_val(bca->frags, kf);
        if (!f)
            return err_msg(-1, 0, "bc_atac_peak_call: frag is NULL");
        g_reg_pair reg = kh_key(bca->frags, kf);

        int np;
        if ( (np = atac_frag_peak_call(f, reg, pks, cmap)) < 0)
            return(-1);
        n_add += np;
    }
    return(n_add);
}

/*******************************************************************************
 * bam_atac_t
 ******************************************************************************/

bam_atac_t *bam_atac_init(){
    bam_atac_t *bam_a = (bam_atac_t *)calloc(1, sizeof(bam_atac_t));
    if (!bam_a){
        err_msg(-1, 0, "bam_atac_init: %s:", strerror(errno));
        return(NULL);
    }

    bam_a->bc_dat = kh_init(khab);
    if (bam_a->bc_dat == NULL){
        err_msg(-1, 0, "bam_atac_init: %s:", strerror(errno));
        return(NULL);
    }
    return(bam_a);
}

void bam_atac_free_dups(bam_atac_t *bam_a){
    if (bam_a == NULL) return;

    khint_t kbc;
    for (kbc = kh_begin(bam_a->bc_dat); kbc != kh_end(bam_a->bc_dat); ++kbc){
        if (!kh_exist(bam_a->bc_dat, kbc)) continue;
        bc_atac_t *bca = kh_val(bam_a->bc_dat, kbc);
        bc_atac_dstry_dups(bca);
    }
}

void bam_atac_free_pairs(bam_atac_t *bam_a){
    if (bam_a == NULL) return;

    khint_t kbc;
    for (kbc = kh_begin(bam_a->bc_dat); kbc != kh_end(bam_a->bc_dat); ++kbc){
        if (!kh_exist(bam_a->bc_dat, kbc)) continue;
        bc_atac_t *bca = kh_val(bam_a->bc_dat, kbc);
        bc_atac_dstry_pairs(bca);
    }
}

void bam_atac_dstry(bam_atac_t *bam_a){
    if (bam_a == NULL) return;

    khint_t kbc;
    for (kbc = kh_begin(bam_a->bc_dat); kbc != kh_end(bam_a->bc_dat); ++kbc){
        if (!kh_exist(bam_a->bc_dat, kbc)) continue;
        bc_atac_t *bca = kh_val(bam_a->bc_dat, kbc);
        char *bc = kh_key(bam_a->bc_dat, kbc);
        free(bc);
        bc_atac_dstry(bca);
    }
    kh_destroy(khab, bam_a->bc_dat);
    free(bam_a);
}

int bam_atac_add_read(bam_atac_t *bam_a, const char *bc, const atac_read1_t *ar, 
        qshort qname){
    if (bam_a == NULL || ar == NULL)
        return err_msg(-1, 0, "bam_atac_add_read: arguments cannot be NULL");

    int ret;
    khint_t kbc;
    kbc = kh_get(khab, bam_a->bc_dat, (char *)bc);
    // if barcode isn't present, add bc_atac
    if (kbc == kh_end(bam_a->bc_dat)){
        char *bc_cpy = strdup(bc);
        kbc = kh_put(khab, bam_a->bc_dat, bc_cpy, &ret);
        if (ret < 0){
            return err_msg(-1, 0, "bam_atac_add_reads: failed to add read to bcs");
        }
        kh_val(bam_a->bc_dat, kbc) = bc_atac_init();
    }
    bc_atac_t *bca = kh_val(bam_a->bc_dat, kbc);
    ret = bc_atac_add_read(bca, ar, qname);
    if (ret < 0)
        return err_msg(-1, 0, "bam_atac_add_reads: failed to add read to bam");
    return(ret);
}

int bam_atac_form_dups(bam_atac_t *bam_a){
    if (!bam_a)
        return err_msg(-1, 0, "bam_atac_form_dups: arguments must not be NULL");

    int ret;
    khint_t kbc;
    for (kbc = kh_begin(bam_a->bc_dat); kbc != kh_end(bam_a->bc_dat); ++kbc){
        if (!kh_exist(bam_a->bc_dat, kbc)) continue;
        bc_atac_t *bca = kh_val(bam_a->bc_dat, kbc);
        if (bca == NULL)
            return err_msg(-1, 0, "bam_atac_form_dups: barcode is NULL");
        ret = bc_atac_form_dups(bca);
        if (ret < 0) return(-1);
    }
    return(0);
}

int bam_atac_dedup(bam_atac_t *bam_a){
    if (!bam_a)
        return err_msg(-1, 0, "bam_atac_dedup: argument bam_a is NULL");

    int ret;
    khint_t kbc;
    for (kbc = kh_begin(bam_a->bc_dat); kbc != kh_end(bam_a->bc_dat); ++kbc){
        if (!kh_exist(bam_a->bc_dat, kbc)) continue;
        bc_atac_t *bca = kh_val(bam_a->bc_dat, kbc);
        if (bca == NULL)
            return err_msg(-1, 0, "bam_atac_dedup: barcode is NULL");

        ret = bc_atac_dedup(bca);
        if (ret < 0) return(-1);
    }
    return(0);
}

int bam_atac_var_call(bam_atac_t *bam_a, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual){
    if (!bam_a)
        return err_msg(-1, 0, "bam_atac_var_call: argument bam_a is NULL");

    int ret, n_add = 0;
    khint_t kbc;
    for (kbc = kh_begin(bam_a->bc_dat); kbc != kh_end(bam_a->bc_dat); ++kbc){
        if (!kh_exist(bam_a->bc_dat, kbc)) continue;
        bc_atac_t *bca = kh_val(bam_a->bc_dat, kbc);
        if (bca == NULL)
            return err_msg(-1, 0, "bam_atac_var_call: barcode is NULL");

        ret = bc_atac_var_call(bca, gv, cmap, min_qual);
        if (ret < 0) return(-1);
        n_add += ret;
    }
    return(n_add);
}

int bam_atac_peak_call(bam_atac_t *bam_a, iregs_t *pks, contig_map *cmap){
    if (!bam_a)
        return err_msg(-1, 0, "bam_atac_peak_call: argument bam_a is NULL");

    int ret, n_add = 0;
    khint_t kbc;
    for (kbc = kh_begin(bam_a->bc_dat); kbc != kh_end(bam_a->bc_dat); ++kbc){
        if (!kh_exist(bam_a->bc_dat, kbc)) continue;
        bc_atac_t *bca = kh_val(bam_a->bc_dat, kbc);
        if (bca == NULL)
            return err_msg(-1, 0, "bam_atac_var_call: barcode is NULL");

        ret = bc_atac_peak_call(bca, pks, cmap);
        if (ret < 0) return(-1);
        n_add += ret;
    }
    return(n_add);
}

/*******************************************************************************
 *
 *
 *
 *
 *
 *
 *
 *
 *
 *
 ******************************************************************************/





/*******************************************************************************
 * miscellaneous
 ******************************************************************************/

void count_frags(bam_atac_t *bam_a, int *n_frag, int *n_reads, int *n_bc){
    *n_frag = 0;
    *n_reads = 0;
    *n_bc = 0;
    khint_t k;
    bc_atac_t *b;

    for (k = kh_begin(bam_a->bc_dat); k != kh_end(bam_a->bc_dat); ++k){
        if (!kh_exist(bam_a->bc_dat, k)) continue;
        (*n_bc)++;

        b = kh_val(bam_a->bc_dat, k);

        *n_frag += kh_size(b->frags);

        if (b->pairs != NULL)
            *n_reads += kh_size(b->pairs);
    }
    printf("n_frag = %i; n_reads = %i; n_bc = %i\n", *n_frag, *n_reads, *n_bc);
}

uint64_t get_n_buckets(bam_atac_t *bam_a){
    if (bam_a == NULL || bam_a->bc_dat == NULL)
        return(0);

    uint64_t nb_frag = 0, nb_dup = 0, nb_pair = 0;
    uint64_t nb_bc = kh_n_buckets(bam_a->bc_dat);
    khint_t kb;
    for (kb = kh_begin(bam_a->bc_dat); kb != kh_end(bam_a->bc_dat); ++kb){
        if (!kh_exist(bam_a->bc_dat, kb)) continue;
        bc_atac_t *bc_atac = kh_val(bam_a->bc_dat, kb);
        if (bc_atac == NULL) continue;
        if (bc_atac->frags) nb_frag += (uint64_t)kh_n_buckets(bc_atac->frags);
        if (bc_atac->dups) nb_dup += (uint64_t)kh_n_buckets(bc_atac->dups);
        if (bc_atac->pairs) nb_pair += (uint64_t)kh_n_buckets(bc_atac->pairs);
    }
    printf("nb_frag = %" PRIu64 "; nb_dup = %" PRIu64 "; "
           "nb_pair = %" PRIu64 "; nb_barcodes = %" PRIu64 "\n", 
            nb_frag, nb_dup, nb_pair, nb_bc);
    return(nb_bc);
}

void print_vac_bam(bam_atac_t *b, GenomeVar *gv){
    if (b == NULL){
        err_msg(-1, 0, "print_vac_bam: arguments must not be NULL");
        return;
    }

    fprintf(stdout, "printing variant allele calls for atac fragments\n");

    khint_t k;
    for (k = kh_begin(b->bc_dat); k != kh_end(b->bc_dat); ++k){
        if (!kh_exist(b->bc_dat, k)) continue;

        // get frags from bam_atac
        char *barcode = kh_key(b->bc_dat, k);
        bc_atac_t *v = kh_val(b->bc_dat, k);
        if (v == NULL) continue;

        size_t n_f = kh_size(v->frags);
        // loop through frags
        khint_t kf;
        for (kf = kh_begin(v->frags); kf != kh_end(v->frags); ++kf){
            if (!kh_exist(v->frags, kf)) continue;
            g_reg_pair regp = kh_key(v->frags, kf);

            atac_frag_t *f_i = kh_val(v->frags, kf);
            if (f_i == NULL) continue;
            fprintf(stdout, "%s: (%zu frags) ", barcode, n_f);
            fprintf(stdout, "r1: %i:%"PRIhts_pos"-%"PRIhts_pos"", 
                    regp.r1.rid, regp.r1.start, regp.r1.end);
            fprintf(stdout, " r2: %i:%"PRIhts_pos"-%"PRIhts_pos"\n", 
                    regp.r2.rid, regp.r2.start, regp.r2.end);

            // print bases
            seq_base_t *sb = f_i->bases.head;
            fprintf(stdout, "bases:\n");
            while (sb != NULL){
                fprintf(stdout, "\t%i", sb->pos.rid);
                fprintf(stdout, "\t%"PRIhts_pos"", sb->pos.pos);
                fprintf(stdout, "\t%u", sb->base);
                fprintf(stdout, "\t%u\n", sb->qual);
                sb = sb->next;
            }

            fprintf(stdout, "variant alleles:\n");
            // print variant calls
            vac_t *va = f_i->vacs.head;
            size_t va_n = f_i->vacs.n;
            fprintf(stdout, "\t%zu variants\n", va_n);
            while (va != NULL){
                Var * var = gv_vari(gv, va->vix); 
                if (var == NULL)
                    fprintf(stdout, "\t%i is invalid\n", va->vix);
                bcf1_t *rec = var->b;
                fprintf(stdout, "\t%s", rec->d.id);
                fprintf(stdout, "\t%"PRIhts_pos"", rec->pos);
                fprintf(stdout, "\t%u\n", va->allele);
                va = va->next;
            }
        }
    }
}

void print_frag_dup(bam_atac_t *b){
    uint32_t counts[100] = {0}; // 1-5, 6-10, 11-15, ...
    int i, m = 1;
    khint_t k;
    for (k = kh_begin(b->bc_dat); k != kh_end(b->bc_dat); ++k){
        if (!kh_exist(b->bc_dat, k)) continue;

        // get frags from bam_atac
        bc_atac_t *v = kh_val(b->bc_dat, k);
        if (v == NULL) continue;

        // loop through frags
        khint_t kf;
        for (kf = kh_begin(v->frags); kf != kh_end(v->frags); ++kf){
            if (!kh_exist(v->frags, kf)) continue;

            atac_frag_t *f_i = kh_val(v->frags, kf);
            if (f_i == NULL) continue;

            for (i = 1; i <= 99; ++i){
                if (f_i->s <= (m * i)){
                    break;
                }
                ++counts[i-1];
            }
        }
    }

    for (i = 0; i < 100; ++i){
        fprintf(stdout, "%i-%i: %" PRIu32 "\n", i*m, (i+1)*m, counts[i]);
    }
    
}

void print_dups_n(bam_atac_t *b){
    // get number of duplicates
    uint32_t counts[100] = {0};
    int i;
    khint_t k;
    for (k = kh_begin(b->bc_dat); k != kh_end(b->bc_dat); ++k){
        if (!kh_exist(b->bc_dat, k)) continue;

        // get dups from bam_atac
        bc_atac_t *v = kh_val(b->bc_dat, k);
        if (v == NULL) continue;

        // loop through dups
        khint_t kd;
        for (kd = kh_begin(v->dups); kd != kh_end(v->dups); ++kd){
            if (!kh_exist(v->dups, kd)) continue;

            atac_dups_t *d_i = kh_val(v->dups, kd);
            if (d_i == NULL) continue;

            size_t s = 0;
            if (d_i->size > 100)
                s = 99;
            else
                s = d_i->size - 1;

            ++counts[s];
        }
    }

    for (i = 0; i < 100; ++i){
        fprintf(stdout, "%i: %" PRIu32 "\n", i+1, counts[i]);
    }
    
}

