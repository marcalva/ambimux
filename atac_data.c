
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
    f->dups = atac_dups_init();
    if (f->dups == NULL) return(NULL);
    init_seq_blist(&f->bases);
    init_vacs(&f->vacs);
    f->pks.ix = NULL;
    f->pks.n = 0;
    f->pks.m = 0;
    f->_dedup = 0;
    f->s = 0;
    return(f);
}

void atac_frag_dstry(atac_frag_t *f){
    if (f == NULL) return;

    atac_dups_dstry(f->dups);
    free_seq_blist(&f->bases);
    free_vacs(&f->vacs);
    free(f->pks.ix);
    free(f);
}

int atac_frag_dedup(atac_frag_t *frag){
    if (frag == NULL)
        return err_msg(-1, 0, "atac_frag_dedup: frag is NULL");

    atac_dups_t *dups = frag->dups;

    // check for bugs
    if (dups->size == 0)
        return err_msg(-1, 0, "atac_frag_dedup: no duplicates present in frag");

    int ix_best;
    size_t max_c = 0; // store read count
    size_t rp_w_max = 0; // number of read pairs with max_c read counts

    int i;
    for (i = 0; i < dups->size; ++i){
        // check for bugs
        if (dups->dups[i].n_rd == 0)
            return err_msg(-1, 0, "atac_frag_dedup: a dup read pair has 0 "
                    "supporting reads, there is a bug");

        if (dups->dups[i].n_rd == max_c){
            ++rp_w_max;
        } else if (dups->dups[i].n_rd > max_c){
            max_c = dups->dups[i].n_rd;
            ix_best = i;
            rp_w_max = 1;
        }
    }

    if (max_c == 0 || rp_w_max == 0)
        return err_msg(-1, 0, "atac_frag_dedup: no best duplicate found "
                ", there is a bug");

    // skip if ambiguous
    if (rp_w_max > 1)
        return(0);

    atac_rd_pair_t rpb = dups->dups[ix_best].rd;

    int bret;
    seq_base_t *base;

    // check for bugs
    if (rpb.s < 2)
        return err_msg(-1, 0, "atac_frag_dedup: dup read pair has less than 2 reads"
                ", there is a bug");

    /* add bases from read 1 */
    base = rpb.r1.bases.head;
    while (base){
        bret = blist_add_base(&frag->bases, base, 1, 1);
        if (bret < 0) return(-1);
        base = base->next;
    }

    /* add bases from read 2 */
    base = rpb.r2.bases.head;
    while (base){
        bret = blist_add_base(&frag->bases, base, 1, 1);
        if (bret < 0) return(-1);
        base = base->next;
    }

    /* add number of supporting reads */
    frag->s = dups->dups[ix_best].n_rd;

    // free dups
    atac_dups_dstry(dups);
    frag->dups = NULL;

    frag->_dedup = 1;

    return 0;
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

int atac_frag_var_call(atac_frag_t *f, g_var_t *gv, contig_map *cmap, 
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
        return err_msg(-1, 0, "atac_frag_peak_call: chromosomes in given region pair don't match");
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
        // add allocated frag // TODO: high heap alloc
        if ( (kh_val(frags, kf) = atac_frag_init()) == NULL) return(-1);
    }
    // add read pair to frag
    atac_frag_t *frag = kh_val(frags, kf);
    if (atac_dups_add_pair(frag->dups, rp) < 0)
        return(-1);

    return(0);
}

/*******************************************************************************
 * miscellaneous
 ******************************************************************************/

/*
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

void print_vac_bam(bam_atac_t *b, g_var_t *gv){
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
                var_t * var = gv_vari(gv, va->vix); 
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

*/

