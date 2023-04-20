
#include <stdlib.h>
#include "str_util.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "htslib/vcf.h"
#include "region.h"

/*******************************************************************************
 * g_region
 ******************************************************************************/

void init_g_region(g_region *reg){
    reg->rid = -1;
    reg->start = -1;
    reg->end = -1;
    reg->strand = '.';
}

void set_region(g_region *reg, int32_t rid, int32_t start, int32_t end, 
        char strand){
    if (reg == NULL) return;
    reg->rid = rid;
    reg->start = start;
    reg->end = end;
    reg->strand = strand;
}

int regioncmp(g_region r1, g_region r2){
    if (r1.rid != r2.rid){
        return(r1.rid - r2.rid);
    }
    if (r1.start != r2.start){
        return(r1.start - r2.start);
    }
    if (r1.end != r2.end){
        return(r1.end - r2.end);
    }
    if (r1.strand != r2.strand){
        return(r1.strand - r2.strand);
    }
    return(0);
}

void print_g_region(FILE *f, g_region g){
    fprintf(f, "%i:%" PRIi32 "-%" PRIi32 " (%c)\n", g.rid, g.start, g.end, g.strand);
}

/*******************************************************************************
 * g_pos
 ******************************************************************************/

void init_g_pos(g_pos *p){
    if (p == NULL) return;
    p->rid = -1;
    p->pos = -1;
    p->strand = '.';
}

int poscmp(g_pos p1, g_pos p2){
    if (p1.rid != p2.rid){
        return(p1.rid - p2.rid);
    }
    if (p1.pos != p2.pos){
        return(p1.pos - p2.pos);
    }
    if (p1.strand != p2.strand){
        return(p1.strand - p2.strand);
    }
    return(0);
}

char *str_g_pos(g_pos p){
    char *s = aloc_sprintf("(Chr ID):%" PRId32 " (Pos):%" PRIhts_pos " (Strand):%c", 
            p.rid, p.pos, p.strand);
    return(s);
}

void fprint_g_pos(FILE *f, g_pos p){
    char *s = str_g_pos(p);
    fprintf(f, "%s", s);
    free(s);
}

/*******************************************************************************
 * g_reg_pair
 ******************************************************************************/

void init_reg_pair(g_reg_pair *rp){
    init_g_region(&rp->r1);
    init_g_region(&rp->r2);
}

g_reg_pair get_reg_pair(g_region r1, g_region r2){
    g_reg_pair g;
    init_reg_pair(&g);
    g.r1 = r1;
    g.r2 = r2;
    return(g);
}

khint_t kh_reg_pair_hash(g_reg_pair p){
    khint_t h = 0;
    h = (h << 4) + (khint_t)(p.r1.rid);
    h = (h << 4) + (khint_t)(p.r2.rid);

    h = (h << 4) + (khint_t)(p.r1.start);
    h = (h << 4) + (khint_t)(p.r1.end);

    h = (h << 4) + (khint_t)(p.r2.start);
    h = (h << 4) + (khint_t)(p.r2.end);

    h = (h << 1) + (khint_t)(p.r1.strand);
    h = (h << 1) + (khint_t)(p.r2.strand);

    return(h);
}

int reg_pair_cmp(g_reg_pair p1, g_reg_pair p2) {
    int r1 = regioncmp(p1.r1, p2.r1);
    if (r1 != 0)
        return r1;

    int r2 = regioncmp(p1.r2, p2.r2);
    if (r2 != 0)
        return r2;

    return 0;
}

int kh_reg_pair_equal(g_reg_pair p1, g_reg_pair p2){

    int r1 = regioncmp(p1.r1, p2.r1);
    if (r1 != 0) return(0);

    int r2 = regioncmp(p1.r2, p2.r2);
    if (r2 != 0) return(0);

    return(1);
}

/*******************************************************************************
 * indexed region
 ******************************************************************************/

int iregn_add_ix(iregn_t *ireg, int ix){
    if (ireg == NULL) return(0);
    while (ireg->n >= ireg->m){
        ireg->m = ireg->m == 0 ? 1 : ireg->m + 1;
        ireg->ix = realloc(ireg->ix, ireg->m * sizeof(int));
        if (ireg->ix == NULL)
            return err_msg(-1, 0, "iregn_add_ix: %s", strerror(errno));
    }
    ireg->ix[ireg->n++] = ix;
    return(0);
}

iregs_t *iregs_init(){
    iregs_t *iregs = malloc(sizeof(iregs_t));
    if (!iregs){
        err_msg(0, 0, "iregs_init: %s", strerror(errno));
        return NULL;
    }
    iregs->idx = NULL;
    iregs->itr = NULL;
    iregs->reg = NULL;
    iregs->n = 0;
    iregs->m = 0;
    iregs->hash = kh_init(kh_reg);
    iregs->chr_map = init_str_map();
    if (iregs->chr_map == NULL || iregs->hash == NULL){
        err_msg(0, 0, "iregs_init: %s", strerror(errno));
        return NULL;
    }
    return(iregs);
}

void iregs_dstry(iregs_t *iregs){
    if (iregs == NULL) return;
    regidx_destroy(iregs->idx);
    if (iregs->itr) regitr_destroy(iregs->itr);
    iregs->itr = NULL;
    free(iregs->reg);
    kh_destroy(kh_reg, iregs->hash);
    destroy_str_map(iregs->chr_map);
    free(iregs);
}

int iregs_add2reghash(iregs_t *iregs, const char *chr, int32_t beg, int32_t end, char strand){
    if (iregs == NULL) return(0);
    if (iregs->reg == NULL){
        iregs->n = 0;
        iregs->m = 2;
        iregs->reg = calloc(iregs->m, sizeof(g_region));
        if (iregs->reg == NULL)
            return err_msg(-1, 0, "iregs_add2reghash:  %s", strerror(errno));
    }
    while (iregs->n >= iregs->m){
        iregs->m *= 2;
        iregs->reg = realloc(iregs->reg, iregs->m * sizeof(g_region));
        if (iregs->reg == NULL)
            return err_msg(-1, 0, "iregs_add2reghash:  %s", strerror(errno));
    }

    if (iregs->chr_map == NULL && (iregs->chr_map = init_str_map()) == NULL)
        return err_msg(-1, 0, "iregs_add2reghash: %s", strerror(errno));
    int found;
    int cix = add2str_map(iregs->chr_map, chr, &found);
    if (cix < 0) return(-1);

    // add to reg
    iregs->reg[iregs->n].rid = (int32_t)cix;
    iregs->reg[iregs->n].start = beg;
    iregs->reg[iregs->n].end = end;
    iregs->reg[iregs->n].strand = strand;
    
    // add to hash
    int ret;
    khint_t k = kh_put(kh_reg, iregs->hash, iregs->reg[iregs->n], &ret);
    if (ret < 0)
        return err_msg(-1, 0, "iregs_add2reghash: failed to add to iregs hash");
    kh_val(iregs->hash, k) = iregs->n;
    ++iregs->n;
    return(0);
}

int iregs_add_bed(iregs_t *iregs, const char *fn){
    if (iregs == NULL) return(0);
    iregs->idx = regidx_init(fn, regidx_parse_bed, NULL, 0, NULL);
    if (iregs->idx == NULL)
        return err_msg(-1, 0, "iregs_add_bed: failed to add bed file %s", fn);
    return(0);
}

// TODO: add strand capacity
int iregs_parse_bed(iregs_t *iregs){
    if (iregs == NULL || iregs->idx == NULL) return(0);
    regitr_t *itr = regitr_init(iregs->idx);
    if (itr == NULL)
        return err_msg(-1, 0, "iregs_parse_bed: failed to init regitr");
    while (regitr_loop(itr))
        if (iregs_add2reghash(iregs, itr->seq, (int32_t)itr->beg, 
                    (int32_t)itr->end, '.') < 0) return(-1);
    regitr_destroy(itr);
    iregs->itr = NULL;
    return(0);
}

// TODO: add strand capacity
int iregs_overlap(iregs_t *iregs, const char *chr, int32_t beg, int32_t end, 
        mv_t(int_vec) *overlaps){
    if (iregs == NULL) return(0);

    if (chr == NULL || overlaps == NULL)
        return err_msg(-1, 0, "iregs_overlap: argument is null");

    regitr_t *itr;
    if ( (itr = regitr_init(iregs->idx)) == NULL )
        return err_msg(-1, 0, "iregs_overlap: failed to initialize regitr");
    if (regidx_overlap(iregs->idx, chr, (hts_pos_t)beg, (hts_pos_t)end, itr) < 0)
        return err_msg(-1, 0, "iregs_overlap: failed to index overlap "
                "%s:%"PRIi32"-%"PRIi32, chr, beg+1, end+1);

    while ( regitr_overlap(itr) ){
        g_region tmp;
        tmp.rid = str_map_ix(iregs->chr_map, itr->seq);
        if (tmp.rid < 0)
            return err_msg(-1, 0, "iregs_overlap: could not find chromosome %s in iregs, "
                    "improper initialization", chr);
        tmp.start = (int32_t)itr->beg;
        tmp.end = (int32_t)itr->end;
        tmp.strand = '.';
        // get the index of the region in iregs->reg to store
        khint_t k_ix = kh_get(kh_reg, iregs->hash, tmp);
        if (k_ix == kh_end(iregs->hash))
            return err_msg(-1, 0, "iregs_overlap: could not find region in hash table, "
                    "there may be a bug");
        int ix       = kh_val(iregs->hash, k_ix);
        if (ix >= iregs->n)
            return err_msg(-1, 0, "iregs_overlap: ix %i > n %i", ix, iregs->n);
        if (mv_push(int_vec, overlaps, ix) < 0) return(-1);
    }
    regitr_destroy(itr);
    iregs->itr = NULL;
    return((int)mv_size(overlaps));
}

// TODO: write payload.
int iregs_write(iregs_t *iregs, BGZF *fp){
    if (iregs == NULL || fp == NULL) return(0);

    size_t strp_size = 20;
    char *strp = malloc(strp_size * sizeof(char));
    int i, len;
    for (i = 0; i < iregs->n; ++i){
        g_region reg = iregs->reg[i];
        const char *chr = str_map_str(iregs->chr_map, (int)(reg.rid));

        len = bgzf_write(fp, chr, strlen(chr));
        len = bgzf_write(fp, "\t", 1);

        len = htspos2strp(reg.start, &strp, &strp_size);
        if (len < 0) return(-1);
        len = bgzf_write(fp, strp, len);
        len = bgzf_write(fp, "\t", 1);

        len = htspos2strp(reg.end + 1, &strp, &strp_size);
        if (len < 0) return(-1);
        len = bgzf_write(fp, strp, len);
        len = bgzf_write(fp, "\n", 1);

        if (len < 0)
            return err_msg(-1, 0, "iregs_write: failed to write to file");
    }

    free(strp);

    return(0);
}

