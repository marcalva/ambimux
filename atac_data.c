
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "atac_data.h"
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
    seq_base_l_init(&ar->bl);
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

    seq_base_l_free(&ar->bl);
}

void atac_read_dstry(atac_read1_t *ar){
    if (ar == NULL) return;
    atac_read_free(ar);
    free(ar);
}

atac_read1_t *atac_read1_dup(const atac_read1_t *r, int *ret){
    *ret = 0;
    if (r == NULL) return(NULL);

    atac_read1_t *cpy = atac_read_init();
    if (cpy == NULL){
        *ret = -1;
        return(NULL);
    }

    cpy->loc = r->loc;

    if ( seq_base_l_cpy(&cpy->bl, &r->bl) < 0 ){
        *ret = -1;
        free(cpy);
        return(NULL);
    }

    return(cpy);
}

int atac_read1_cpy(atac_read1_t *dest, const atac_read1_t *src){
    if (src == NULL)
        return err_msg(-1, 0, "atac_read1_cpy: src is null");
    if (dest == NULL)
        return err_msg(-1, 0, "atac_read1_cpy: dest is null");

    dest->loc = src->loc;
    
    if ( seq_base_l_cpy(&dest->bl, &src->bl) < 0 )
        return(-1);

    return(0);
}

int atac_read_equal(atac_read1_t r1, atac_read1_t r2){

    /* compare region */
    int rcmp = regioncmp(r1.loc, r2.loc);
    if (rcmp != 0)
        return(0);

    /* compare bases, ignore qual */
    int bcmp = seq_base_l_equal(r1.bl, r2.bl, 0);
    if (bcmp != 1)
        return(0);
    return(1);
}

int atac_read1_add_base(atac_read1_t *r, seq_base_t base){
    if (r == NULL)
        return err_msg(-1, 0, "atac_read1_add_base: arguments are NULL");

    if (seq_base_l_insert(&r->bl, base, 1, 1) < 0) return(-1);

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

atac_rd_pair_t *atac_rd_pair_dup(const atac_rd_pair_t *rp, int *ret){
    *ret = 0;
    if (rp == NULL) return(NULL);

    atac_rd_pair_t *cp = (atac_rd_pair_t *)calloc(1, sizeof(atac_rd_pair_t));
    if (cp == NULL){
        *ret = err_msg(-1, 0, "atac_rd_pair_dup: %s", strerror(errno));
        return(NULL);
    }

    atac_read1_t *rtmp = NULL;

    int cret = 0;
    rtmp = atac_read1_dup(&rp->r1, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    cp->r1 = *rtmp;
    free(rtmp);

    rtmp = atac_read1_dup(&rp->r2, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    cp->r2 = *rtmp;
    free(rtmp);

    cp->s = rp->s;

    return(cp);
}

int atac_rd_pair_cpy(atac_rd_pair_t *dest, const atac_rd_pair_t *src){
    if (src == NULL)
        return err_msg(-1, 0, "atac_rd_pair_cpy: dest is null");
    if (dest == NULL)
        return err_msg(-1, 0, "atac_rd_pair_cpy: dest is null");

    if (atac_read1_cpy(&dest->r1, &src->r1) < 0)
        return(-1);
    if (atac_read1_cpy(&dest->r2, &src->r2) < 0) // TODO: free memory from r1 if error
        return(-1);

    dest->s = src->s;

    return(0);
}

int atac_rd_pair_add_read(atac_rd_pair_t *rp, const atac_read1_t *ar){
    if (rp == NULL || ar == NULL)
        return err_msg(-1, 0, "atac_rd_pair_add_read: arguments are null");

    int cret = 0;
    atac_read1_t *rtmp = NULL;
    rtmp = atac_read1_dup(ar, &cret);
    if (cret < 0) return(-1);

    if (rp->s == 0){
        rp->r1 = *rtmp;
        free(rtmp);
    } else if (rp->s == 1){
        rp->r2 = *rtmp;
        free(rtmp);

        /* switch order if r2 comes before r1 */
        int rcmp = regioncmp(rp->r2.loc, rp->r1.loc);
        if (rcmp < 0){
            atac_read1_t tmp = rp->r2;
            rp->r2 = rp->r1;
            rp->r1 = tmp;
        }
    } else {
        free(rtmp);
        return err_msg(-1, 0, "atac_rd_pair_add_read: trying to add read to full pair "
                "(query name was found a third time)");
    }
    ++rp->s;
    return(0);
}

int atac_rd_pair_equal(const atac_rd_pair_t *rp1, const atac_rd_pair_t *rp2){
    if (rp1 == NULL || rp2 == NULL)
        return err_msg(-1, 0, "atac_rd_pair_equal: one of read pairs is null");

    if (rp1->s != rp2->s)
        return(0);

    int rcmp1 = atac_read_equal(rp1->r1, rp2->r1); // compare read 1
    if (rcmp1 != 1) return(0);

    int rcmp2 = atac_read_equal(rp1->r2, rp2->r2); // compare read 2
    if (rcmp2 != 1) return(0);

    return(1);
}

int atac_rd_pair_match_qual(atac_rd_pair_t *rp, const atac_rd_pair_t *cmp){
    if (rp == NULL || cmp == NULL)
        return err_msg(-1, 0, "atac_rd_pair_match_qual: arguments are null");
    if (rp->s != cmp->s)
        return err_msg(-1, 0, "atac_rd_pair_match_qual: "
                "number of reads don't match");
    if (rp->s == 0)
        return(0);

    if (seq_base_l_match_qual(&rp->r1.bl, &cmp->r1.bl) < 0)
        return(-1);
    if (rp->s == 2 && seq_base_l_match_qual(&rp->r2.bl, &cmp->r2.bl) < 0)
        return(-1);

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
    size_t i;
    for (i = 0; i < d->size; ++i)
        atac_rd_pair_free(&d->dups[i].rd); // free underlying reads
    free(d->dups);
    d->dups = NULL;
    d->size = 0;
    d->m = 0;
}

void atac_dups_dstry(atac_dups_t *d){
    if (d == NULL) return;
    atac_dups_free(d);
    free(d);
}

int atac_dups_add_pair(atac_dups_t *d, const atac_rd_pair_t *rp){
    if (d == NULL || rp == NULL)
        return err_msg(-1, 0, "atac_dups_add_pair: arguments must not be NULL");

    size_t i;
    int tcmp = 0;
    for (i = 0; i < d->size; ++i){
        tcmp = atac_rd_pair_equal(&d->dups[i].rd, rp);
        if (tcmp < 0) return(-1);
        if (tcmp == 1){
            if (atac_rd_pair_match_qual(&d->dups[i].rd, rp) < 0) return(-1);
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
    if (atac_rd_pair_cpy(&d->dups[d->size].rd, rp) < 0)
        return(-1);

    d->dups[d->size].n_rd = 1;
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
    f->dups = atac_dups_init();
    if (f->dups == NULL) return(NULL);
    seq_base_l_init(&f->bl);
    seq_vac_l_init(&f->vl);
    mv_init(&f->pks);
    f->_dedup = 0;
    f->s = 0;
    return(f);
}

void atac_frag_dstry(atac_frag_t *f){
    if (f == NULL) return;

    atac_dups_dstry(f->dups);
    seq_base_l_free(&f->bl);
    seq_vac_l_free(&f->vl);
    mv_free(&f->pks);
    free(f);
}

int atac_frag_dedup(atac_frag_t *frag){
    if (frag == NULL)
        return err_msg(-1, 0, "atac_frag_dedup: frag is NULL");

    atac_dups_t *dups = frag->dups;
    assert(dups != NULL);

    // check for bugs
    if (dups->size == 0)
        return err_msg(-1, 0, "atac_frag_dedup: no duplicates present in frag");

    int ix_best = 0;
    size_t max_c = 0; // store read count
    size_t rp_w_max = 0; // number of read pairs with max_c read counts

    size_t i;
    for (i = 0; i < dups->size; ++i){
        assert(dups->dups[i].n_rd > 0);

        if (dups->dups[i].n_rd == max_c){
            ++rp_w_max;
        } else if (dups->dups[i].n_rd > max_c){
            max_c = dups->dups[i].n_rd;
            ix_best = i;
            rp_w_max = 1;
        }
    }

    assert(max_c > 0 && rp_w_max > 0);

    // skip if ambiguous
    if (rp_w_max > 1){
        atac_dups_dstry(dups);
        frag->dups = NULL;
        return(0);
    }

    atac_rd_pair_t rpb = dups->dups[ix_best].rd;

    ml_node_t(seq_base_l) *bn;

    assert(rpb.s == 2);

    /* add bases from read 1 */
    for (bn = ml_begin(&rpb.r1.bl); bn; bn = ml_node_next(bn)){
        seq_base_t base = ml_node_val(bn);
        if (seq_base_l_insert(&frag->bl, base, 1, 1) < 0) return(-1);
    }

    /* add bases from read 2 */
    for (bn = ml_begin(&rpb.r2.bl); bn; bn = ml_node_next(bn)){
        seq_base_t base = ml_node_val(bn);
        if (seq_base_l_insert(&frag->bl, base, 1, 1) < 0) return(-1);
    }

    /* add number of supporting reads */
    frag->s = dups->dups[ix_best].n_rd;

    // free dups
    atac_dups_dstry(dups);
    frag->dups = NULL;

    frag->_dedup = 1;

    return 0;
}

int atac_frag_var_call(atac_frag_t *f, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (f == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "atac_frag_var_call: arguments must not be NULL");
    }

    int ret = seq_vac_l_call_var(&f->bl, &f->vl, gv, cmap, min_qual);
    // TODO: remove this
    if (ret < 0){
        ml_node_t(seq_base_l) *n = ml_begin(&f->bl);
        while (n != NULL){
            seq_base_t b = ml_node_val(n);
            fprint_g_pos(stderr, b.pos);
            fprintf(stderr, "\n");
            n = ml_node_next(n);
        }
    }

    // free the bases since we don't need them anymore
    seq_base_l_free(&f->bl);

    return(ret);
}

int atac_frag_peak_call(atac_frag_t *f, g_reg_pair reg, iregs_t *pks, str_map *cmap){
    if (f == NULL)
        return err_msg(-1, 0, "atac_frag_peak_call: f is NULL");
    if (pks == NULL)
        return err_msg(-1, 0, "atac_frag_peak_call: pks is NULL");

    if (reg.r1.rid != reg.r2.rid)
        return err_msg(-1, 0, "atac_frag_peak_call: chromosomes in given region pair don't match");
    int rid = (int)reg.r1.rid;
    const char *chr = str_map_str(cmap, rid);

    int32_t beg = reg.r1.start;
    int32_t end = reg.r2.end - 1;
    if (end < beg)
        return err_msg(-1, 0, "atac_frag_peak_call: end %" PRIi32 " < beg %" PRIi32, end, beg);

    mv_t(int_vec) overlaps;
    mv_init(&overlaps);
    int ret = iregs_overlap(pks, chr, beg, end, &overlaps);
    if (ret < 0) return(-1);

    // check if peak overlaps by at least one base
    // if so, add to f->pks
    size_t i;
    for (i = 0; i < mv_size(&overlaps); ++i){
        int ix = mv_i(&overlaps, i);
        if (ix >= pks->n)
            return err_msg(-1, 0, "atac_frag_peak_call: %zu index >= num of peaks %zu", ix, pks->n);
        g_region pk_reg = pks->reg[ix];

        int64_t ovrlp;

        // check cut site 1 (+4 based on 10X)
        ovrlp = bp_overlap(beg+4, beg+5, '.', pk_reg.start, pk_reg.end, '.');
        if (ovrlp < 0) return(-1);
        if (ovrlp > 0 && mv_push(int_vec, &f->pks, ix) < 0)
            return(-1);

        // check cut site 2 (-5 based on 10X)
        ovrlp = bp_overlap(end-5, end-4, '.', pk_reg.start, pk_reg.end, '.');
        if (ovrlp < 0) return(-1);
        if (ovrlp > 0 && mv_push(int_vec, &f->pks, ix) < 0)
            return(-1);
    }

    mv_free(&overlaps);

    return(mv_size(&f->pks));
}

/*******************************************************************************
 * frags hash table
 ******************************************************************************/

int khaf_add_dup(khash_t(khaf) *frags, atac_rd_pair_t *rp, int skip_chim){
    if (frags == NULL || rp == NULL)
        return err_msg(-1, 0, "khaf_add_dup: 'frags' or 'rp' is null");

    // skip if chimeric
    if (skip_chim && rp->r1.loc.rid != rp->r2.loc.rid)
        return(0);

    int kret;
    khint_t kf;
    g_reg_pair reg_pair = get_reg_pair(rp->r1.loc, rp->r2.loc);
    kf = kh_get(khaf, frags, reg_pair);
    if (kf == kh_end(frags)){
        // add region key
        kf = kh_put(khaf, frags, reg_pair, &kret);// TODO: high heap alloc
        if (kret < 0){
            return err_msg(-1 , 0, "khaf_add_dup: failed to add "
                    "duplicate reg key to hash table");
        }
        // TODO: high heap alloc
        if ( (kh_val(frags, kf) = atac_frag_init()) == NULL) return(-1);
    }
    // add read pair to frag
    atac_frag_t *frag = kh_val(frags, kf);
    if (atac_dups_add_pair(frag->dups, rp) < 0)
        return(-1);

    return(0);
}

