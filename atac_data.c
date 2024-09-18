
#include <stdlib.h>
#include <string.h>
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
    atac_read_set0(ar);
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
        err_msg(-1, 0, "atac_read1_dup: %s", strerror(errno));
        return(NULL);
    }

    cpy->loc = r->loc;

    if ( seq_base_l_cpy(&cpy->bl, &r->bl) < 0 ){
        *ret = -1;
        err_msg(-1, 0, "atac_read1_dup: %s", strerror(errno));
        atac_read_dstry(cpy);
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
    
    if ( seq_base_l_cpy(&dest->bl, &src->bl) < 0 ) {
        atac_read_free(dest);
        return err_msg(-1, 0, "atac_read1_cpy: %s", strerror(errno));
    }

    return(0);
}

int atac_read1_cmp(atac_read1_t r1, atac_read1_t r2){

    /* compare region */
    int rcmp = regioncmp(r1.loc, r2.loc);
    if (rcmp != 0)
        return rcmp;

    /* compare bases, ignore qual */
    int bcmp = seq_base_l_cmp(r1.bl, r2.bl, 0);
    if (bcmp != 0)
        return bcmp;

    return 0;
}

int atac_read1_add_base(atac_read1_t *r, seq_base_t base){
    if (r == NULL)
        return err_msg(-1, 0, "atac_read1_add_base: arguments are NULL");

    if (seq_base_l_insert(&r->bl, base, 1, 1) < 0)
        return err_msg(-1, 0, "atac_read1_add_base: failed to add base to read");

    return(0);
}

/*******************************************************************************
 * atac_rd_pair_t
 ******************************************************************************/

void atac_rd_pair_set0(atac_rd_pair_t *rp){
    if (rp == NULL) return;
    rp->qname = 0;
    atac_read_set0(&rp->r1);
    atac_read_set0(&rp->r2);
    rp->s = 0;
    rp->n_reads = 0;
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

    rp->qname = 0;
    atac_read_free(&rp->r1);
    atac_read_set0(&rp->r1);
    atac_read_free(&rp->r2);
    atac_read_set0(&rp->r2);
    rp->s = 0;
    rp->n_reads = 0;
}

void atac_rd_pair_dstry(atac_rd_pair_t *rp){
    if (rp == NULL) return;

    atac_rd_pair_free(rp);
    free(rp);
}

int atac_rd_pair_is_chimeric(const atac_rd_pair_t *rp) {
    if (rp == NULL)
        return err_msg(-1, 0, "atac_rd_pair_is_chimeric: read pair is null");

    if (rp->s > 2)
        return err_msg(-1, 0, "atac_rd_pair_is_chimeric: read pair has %u reads", rp->s);

    if (rp->s < 2)
        return(1);

    if (rp->r1.loc.rid != rp->r2.loc.rid)
        return(1);

    return(0);
}

atac_rd_pair_t *atac_rd_pair_dup(const atac_rd_pair_t *rp, int *ret){
    *ret = 0;
    if (rp == NULL) {
        *ret = err_msg(-1, 0, "atac_rd_pair_dup: read pair is null");
        return(NULL);
    }

    atac_rd_pair_t *cp = (atac_rd_pair_t *)calloc(1, sizeof(atac_rd_pair_t));
    if (cp == NULL){
        *ret = err_msg(-1, 0, "atac_rd_pair_dup: %s", strerror(errno));
        return(NULL);
    }

    cp->qname = rp->qname;

    atac_read1_t *rtmp = NULL;

    int cret = 0;
    rtmp = atac_read1_dup(&rp->r1, &cret);
    if (cret < 0){
        atac_rd_pair_dstry(cp);
        *ret = err_msg(01, 0, "atac_rd_pair_dup: failed to duplicate read 1");
        return(NULL);
    }
    if (rtmp != NULL) {
        cp->r1 = *rtmp;
        free(rtmp);
    }

    rtmp = atac_read1_dup(&rp->r2, &cret);
    if (cret < 0){
        atac_rd_pair_dstry(cp);
        *ret = err_msg(01, 0, "atac_rd_pair_dup: failed to duplicate read 2");
        return(NULL);
    }
    if (rtmp != NULL) {
        cp->r2 = *rtmp;
        free(rtmp);
    }

    cp->s = rp->s;
    cp->n_reads = rp->n_reads;

    return(cp);
}

int atac_rd_pair_cpy(atac_rd_pair_t *dest, const atac_rd_pair_t *src){
    if (src == NULL)
        return err_msg(-1, 0, "atac_rd_pair_cpy: src is null");
    if (dest == NULL)
        return err_msg(-1, 0, "atac_rd_pair_cpy: dest is null");

    dest->qname = src->qname;

    if (atac_read1_cpy(&dest->r1, &src->r1) < 0) {
        atac_rd_pair_free(dest);
        return(-1);
    }
    if (atac_read1_cpy(&dest->r2, &src->r2) < 0) {
        atac_rd_pair_free(dest);
        return(-1);
    }

    dest->s = src->s;
    dest->n_reads = src->n_reads;

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
        /* switch order if r2 comes before r1 */
        int rcmp = regioncmp(rtmp->loc, rp->r1.loc);
        if (rcmp < 0){
            rp->r2 = rp->r1;
            rp->r1 = *rtmp;
        } else {
            rp->r2 = *rtmp;
        }
        free(rtmp);
    } else {
        atac_read_dstry(rtmp);
        return err_msg(0, 1, "atac_rd_pair_add_read: trying to add read to full pair "
                "(query name was found a third time)");
    }
    ++rp->s;
    return(0);
}

int atac_rd_pair_cmp(const atac_rd_pair_t *rp1, const atac_rd_pair_t *rp2) {
    if (rp1 == NULL || rp2 == NULL)
        return err_msg(-1, 0, "atac_rd_pair_cmp: one of read pairs is null");

    if (rp1->s < rp2->s)
        return -1;
    else if (rp1->s > rp2->s)
        return 1;

    int cmp1 = atac_read1_cmp(rp1->r1, rp2->r1);
    if (cmp1 != 0)
        return cmp1;

    int cmp2 = atac_read1_cmp(rp1->r2, rp2->r2);
    if (cmp2 != 0)
        return cmp2;

    return 0;
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

int atac_dups_init(atac_dups_t *d){
    if (d == NULL)
        return err_msg(-1, 0, "atac_dups_init: argument is null");

    init_reg_pair(&d->reg);

    mv_init(&d->dups);

    return 0;
}

atac_dups_t *atac_dups_alloc() {
    atac_dups_t *d = calloc(1, sizeof(atac_dups_t));
    if (d == NULL) {
        err_msg(-1, 0, "atac_dups_alloc: %s", strerror(errno));
        return NULL;
    }
    atac_dups_init(d);
    return d;
}

void atac_dups_free(atac_dups_t *d){
    if (d == NULL) return;

    init_reg_pair(&d->reg);

    size_t i;
    for (i = 0; i < mv_size(&d->dups); ++i){
        atac_rd_pair_t *rp = &mv_i(&d->dups, i);
        atac_rd_pair_free(rp);
    }
    mv_free(&d->dups);
}

void atac_dups_dstry(atac_dups_t *d){
    if (d == NULL) return;
    atac_dups_free(d);
    free(d);
}

int atac_dups_add_pair(atac_dups_t *d, const atac_rd_pair_t *rp){
    if (d == NULL || rp == NULL)
        return err_msg(-1, 0, "atac_dups_add_pair: arguments must not be NULL");

    // copy read pair
    int ret;
    atac_rd_pair_t *tmp = atac_rd_pair_dup(rp, &ret);
    if (ret < 0 || tmp == NULL)
        return err_msg(-1, 0, "atac_dups_add_pair: failed to duplicate read pair");

    size_t i;
    for (i = 0; i < mv_size(&d->dups); ++i){
        atac_rd_pair_t *rp_vec = &mv_i(&d->dups, i);
        // if read pair is found, set existing pair's base qualities to max
        //  of the two.
        if (atac_rd_pair_cmp(rp_vec, tmp) == 0){
            if (atac_rd_pair_match_qual(rp_vec, tmp) < 0){
                atac_rd_pair_dstry(tmp);
                return err_msg(-1, 0, "atac_dups_add_pair: failed to match quality of read pair");
            }
            rp_vec->n_reads += 1;
            atac_rd_pair_dstry(tmp);
            tmp = NULL;
            return(0);
        }
    }
    // if not found, add to vector
    if (mv_push(mv_arp, &d->dups, *tmp) < 0){
        atac_rd_pair_dstry(tmp);
        return err_msg(-1, 0, "atac_dups_add_pair: failed to add read pair to vector");
    }
    free(tmp);

    return 0;
}

/*******************************************************************************
 * atac_frag_t
 ******************************************************************************/

int atac_frag_init(atac_frag_t *f){
    if (f == NULL)
        return err_msg(-1, 0, "atac_frag_init: argument is null");
    init_reg_pair(&f->reg);
    seq_base_l_init(&f->bl);
    seq_vac_l_init(&f->vl);
    mv_init(&f->pks);
    f->s = 0;
    return 1;
}

atac_frag_t *atac_frag_alloc() {
    atac_frag_t *f = (atac_frag_t *)calloc(1, sizeof(atac_frag_t));
    if (f == NULL){
        err_msg(-1, 0, "atac_frag_alloc: %s", strerror(errno));
        return(NULL);
    }
    if (atac_frag_init(f) < 0) {
        atac_frag_dstry(f);
        return NULL;
    }
    return f;
}

void atac_frag_free(atac_frag_t *f) {
    if (f == NULL) return;

    init_reg_pair(&f->reg);
    seq_base_l_free(&f->bl);
    seq_vac_l_free(&f->vl);
    mv_free(&f->pks);
}

void atac_frag_dstry(atac_frag_t *f){
    if (f == NULL) return;
    atac_frag_free(f);
    free(f);
}

atac_frag_t *atac_dups_dedup(atac_dups_t *dups, int *ret) {
    *ret = 0;
    if (dups == NULL) {
        *ret = err_msg(-1, 0, "atac_dups_dedup: argument is null");
        return NULL;
    }

    if (mv_size(&dups->dups) == 0) {
        *ret = err_msg(-1, 0, "atac_dups_dedup: no read pairs in 'dups'");
        return NULL;
    }

    uint32_t max_c = 0; // store read count
    atac_rd_pair_t *rp_best = NULL;
    uint32_t rp_w_max = 0; // number of read pairs with max_c read counts

    size_t i;
    for (i = 0; i < mv_size(&dups->dups); ++i){
        atac_rd_pair_t *rp = &mv_i(&dups->dups, i);
        uint32_t rds_per_dup = rp->n_reads;
        if (rds_per_dup == 0) {
            *ret = err_msg(-1, 0, "atac_dups_dedup: read pair has no reads");
            return NULL;
        }
        if (rds_per_dup == max_c){
            ++rp_w_max;
        } else if (rds_per_dup > max_c){
            max_c = rds_per_dup;
            rp_best = rp;
            rp_w_max = 1;
        }
    }

    if (rp_w_max == 0) {
        *ret = err_msg(-1, 0, "atac_dups_dedup: no read pairs with reads");
        return NULL;
    }
    if (rp_best == NULL) {
        *ret = err_msg(-1, 0, "atac_dups_dedup: failed to find best read pair");
        return NULL;
    }

    // skip if ambiguous
    if (rp_w_max > 1)
        return NULL;

    // check that read pair is complete
    if (rp_best->s != 2) {
        *ret = err_msg(-1, 0, "atac_dups_dedup: read pair has %u reads", rp_best->s);
        return NULL;
    }

    atac_frag_t *frag = atac_frag_alloc();
    if (frag == NULL) {
        *ret = err_msg(-1, 0, "atac_dups_dedup: failed to allocate fragment");
        return NULL;
    }

    // add base calls from both reads of pair to fragment
    ml_node_t(seq_base_l) *bn;

    /* add bases from read 1 */
    for (bn = ml_begin(&rp_best->r1.bl); bn; bn = ml_node_next(bn)){
        seq_base_t base = ml_node_val(bn);
        if (seq_base_l_insert(&frag->bl, base, 1, 1) < 0) {
            atac_frag_dstry(frag);
            *ret = err_msg(-1, 0, "atac_dups_dedup: failed to insert base");
            return NULL;
        }
    }

    /* add bases from read 2 */
    for (bn = ml_begin(&rp_best->r2.bl); bn; bn = ml_node_next(bn)){
        seq_base_t base = ml_node_val(bn);
        if (seq_base_l_insert(&frag->bl, base, 1, 1) < 0) {
            atac_frag_dstry(frag);
            *ret = err_msg(-1, 0, "atac_dups_dedup: failed to insert base");
            return NULL;
        }
    }

    /* add number of supporting reads */
    frag->s = rp_w_max;
    frag->reg = dups->reg;

    return frag;
}

int atac_frag_var_call(atac_frag_t *f, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (f == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "atac_frag_var_call: arguments must not be NULL");
    }

    int ret = seq_vac_l_call_var(&f->bl, &f->vl, gv, cmap, min_qual);

    // free the bases since we don't need them anymore
    seq_base_l_free(&f->bl);

    return(ret);
}

int atac_frag_peak_call(atac_frag_t *f, g_reg_pair reg, iregs_t *pks, str_map *cmap){
    if (f == NULL)
        return err_msg(-1, 0, "atac_frag_peak_call: f is NULL");
    if (pks == NULL)
        return err_msg(-1, 0, "atac_frag_peak_call: pks is NULL");

    // check if chimeric
    if (reg.r1.rid != reg.r2.rid)
        return err_msg(-1, 0, "atac_frag_peak_call: chromosomes in given region pair don't match");

    int rid = (int)reg.r1.rid;
    const char *chr = str_map_str(cmap, rid);

    // if chr not present, no peaks to overlap
    if (chr == NULL)
        return 0;

    int32_t beg = reg.r1.start;
    int32_t end = reg.r2.end - 1;
    if (end < beg)
        return err_msg(-1, 0, "atac_frag_peak_call: end %" PRIi32 " < beg %" PRIi32, end, beg);

    mv_t(int_vec) overlaps;
    mv_init(&overlaps);
    int ret = iregs_overlap(pks, chr, beg, end, &overlaps);
    if (ret < 0) {
        mv_free(&overlaps);
        return err_msg(-1, 0, "atac_frag_peak_call: failed to overlap peaks");
    }

    // check if the individual cut sites overlap peaks
    // A fragment fully contained within a peak will contain two counts/instances
    //   for the peak, one for each cut site.
    size_t i;
    for (i = 0; i < mv_size(&overlaps); ++i){
        int ix = mv_i(&overlaps, i);
        if (ix >= pks->n || ix < 0) {
            mv_free(&overlaps);
            return err_msg(-1, 0, "atac_frag_peak_call: %zu index >= num of peaks %zu", ix, pks->n);
        }
        g_region pk_reg = pks->reg[ix];

        int64_t ovrlp;

        // check cut site 1 (+4 based on 10X)
        ovrlp = bp_overlap(beg+4, beg+5, '.', pk_reg.start, pk_reg.end, '.');
        if (ovrlp < 0) {
            mv_free(&overlaps);
            return err_msg(-1, 0, "atac_frag_peak_call: failed to overlap cut site 1");
        }
        if (ovrlp > 0 && mv_push(int_vec, &f->pks, ix) < 0) {
            mv_free(&overlaps);
            return err_msg(-1, 0, "atac_frag_peak_call: failed to push peak index");
        }

        // check cut site 2 (-5 based on 10X)
        ovrlp = bp_overlap(end-5, end-4, '.', pk_reg.start, pk_reg.end, '.');
        if (ovrlp < 0) {
            mv_free(&overlaps);
            return err_msg(-1, 0, "atac_frag_peak_call: failed to overlap cut site 2");
        }
        if (ovrlp > 0 && mv_push(int_vec, &f->pks, ix) < 0) {
            mv_free(&overlaps);
            return err_msg(-1, 0, "atac_frag_peak_call: failed to push peak index");
        }
    }

    mv_free(&overlaps);

    return(mv_size(&f->pks));
}

