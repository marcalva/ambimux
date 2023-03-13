
#include <stdlib.h>
#include "str_util.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "counts.h"
#include "variants.h"
#include "region.h"
#include "gtf_anno.h"
#include "g_list.h"

/*******************************************************************************
 * seq_base_t
 ******************************************************************************/

void seq_base_init(seq_base_t *sbase){
    init_g_pos(&sbase->pos);
    sbase->base = 15;
    sbase->qual = 0;
}

seq_base_t *seq_base_alloc(){
    seq_base_t *b = (seq_base_t *)calloc(1, sizeof(seq_base_t));
    if (b == NULL){
        err_msg(-1, 0, "seq_base_alloc: %s", strerror(errno));
        return(NULL);
    }

    seq_base_init(b);
    return(b);
}

void seq_base_dstry(seq_base_t *b){
    if (b == NULL) return;
    free(b);
}

seq_base_t *seq_base_dup(const seq_base_t *b, int *ret){
    *ret = 0;
    if (b == NULL) return(NULL);

    seq_base_t *c = (seq_base_t *)calloc(1, sizeof(seq_base_t));
    if (c == NULL){
        *ret = err_msg(-1, 0, "seq_base_dup: %s", strerror(errno));
        return(NULL);
    }

    c->pos = b->pos;
    c->base = b->base;
    c->qual = b->qual;

    return(c);
}

/*******************************************************************************
 * base list
 ******************************************************************************/

int seq_base_l_cmp(ml_t(seq_base_l) bl1, ml_t(seq_base_l) bl2, int cmp_qual){
    int sc = ml_size(&bl1) - ml_size(&bl2);
    if (sc != 0) return(sc);

    ml_node_t(seq_base_l) *n1 = ml_begin(&bl1);
    ml_node_t(seq_base_l) *n2 = ml_begin(&bl2);
    while (n1 != NULL && n2 != NULL){
        seq_base_t b1 = ml_node_val(n1);
        seq_base_t b2 = ml_node_val(n2);
        int bc = seq_base_cmp(b1, b2, cmp_qual);
        if (bc != 0) return(bc);
        n1 = ml_node_next(n1);
        n2 = ml_node_next(n2);
    }
    return(0);
}

int seq_base_l_equal(ml_t(seq_base_l) bl1, ml_t(seq_base_l) bl2, int cmp_qual){
    if (ml_size(&bl1) != ml_size(&bl2))
        return(0);

    ml_node_t(seq_base_l) *n1 = ml_begin(&bl1);
    ml_node_t(seq_base_l) *n2 = ml_begin(&bl2);
    while (n1 != NULL && n2 != NULL){
        seq_base_t b1 = ml_node_val(n1);
        seq_base_t b2 = ml_node_val(n2);
        if ( seq_base_equal(b1, b2, cmp_qual) != 1 )
            return(0);
        n1 = ml_node_next(n1);
        n2 = ml_node_next(n2);
    }
    return(1);
}

int seq_base_l_match_qual(ml_t(seq_base_l) *bl, const ml_t(seq_base_l) *cmp){
    if (bl == NULL || cmp == NULL)
        return err_msg(-1, 0, "seq_base_l_match_qual: argument is null");

    if (ml_size(bl) != ml_size(cmp))
        return err_msg(-1, 0, "seq_base_l_match_qual: "
                "number of bases don't match (%zu != %zu)", 
                ml_size(bl), ml_size(cmp));

    ml_node_t(seq_base_l) *n1 = ml_begin(bl), *n2 = ml_begin(cmp);
    while (n2 != NULL){
        if (seq_base_cmp(ml_node_val(n1), ml_node_val(n2), 0) != 0)
            return err_msg(-1, 0, "seq_base_l_match_qual: bases don't match");
        if ( ml_node_val(n2).qual > ml_node_val(n1).qual )
            ml_node_val(n1).qual = ml_node_val(n2).qual;

        n1 = ml_node_next(n1);
        n2 = ml_node_next(n2);
    }

    return(0);
}

/*******************************************************************************
 * seq_vac_t
 ******************************************************************************/

void seq_vac_init(seq_vac_t *v){
    v->vix = -1;
    v->allele = 0xf; // missing
    v->qual = 0;
}

seq_vac_t *seq_vac_alloc(){
    seq_vac_t *v = (seq_vac_t *)calloc(1, sizeof(seq_vac_t));
    if (v == NULL){
        err_msg(-1, 0, "seq_vac_alloc: %s", strerror(errno));
        return NULL;
    }
    seq_vac_init(v);
    return(v);
}

void seq_vac_dstry(seq_vac_t *v){
    if (v == NULL) return;
    free(v);
}

int seq_base_call_var(seq_base_t b, ml_t(seq_vac_l) *vl, g_var_t *gv, 
        str_map *cmap, uint8_t min_qual){
    if (vl == NULL || gv == NULL || cmap == NULL)
        return err_msg(-1, 0, "seq_base_call_var: argument is null");

    /* get chromosome and position */
    int32_t b_rid = b.pos.rid;
    const char *b_ref = str_map_str(cmap, b_rid);
    if (b_ref == NULL)
        return err_msg(-1, 0, "seq_base_call_var: cannot find chromosome ID %i", b_rid);
    int32_t b_beg = b.pos.pos, b_end = b_beg + 1;

    /* base call */
    uint8_t b_base = b.base;
    uint8_t b_qual = b.qual;

    // skip bases observed as N
    if (b_base == 0xf)
        return(0);

    if (b_qual < min_qual)
        return(0);

    /* Get overlapping variants */
    int n_added = 0;
    ml_t(vcfr_list) vars;
    ml_init(vcfr_list, &vars);
    int n_v = g_var_get_region_vars(gv, b_ref, b_beg, b_end, &vars);
    if (n_v < 0) return(-1);
    // debugging to check for mutliple snps at same site
    // if (1){
    if (n_v > 1 || ml_size(&vars) > 1){
        err_msg(0, 1, "seq_base_call_var: %i variants found at position %s:%"PRIi32"\n", 
                n_v, b_ref, b_beg);
    }
    // }

    /* loop through overlapping variants (should just be one) */
    ml_node_t(vcfr_list) *vn;
    for (vn = ml_begin(&vars); vn; vn = ml_node_next(vn)){
        var_t var = ml_node_val(vn);
        bcf1_t *rec = var.b;
        // skip if position doesn't overlap
        if (rec->pos != b_beg)
            continue;

        // skip indels
        if (!bcf_is_snp(rec))
            continue;

        // throw error if more alleles than present than can be stored.
        // last allele is reserved for missing/other
        if (rec->n_allele > (MAX_ALLELE-1))
            return err_msg(-1, 0, "seq_base_call_var: too many alleles: "
                    "%i > %i", rec->n_allele, MAX_ALLELE-1);

        // set up vac
        seq_vac_t v_i;
        v_i.vix = var.vix;
        assert(v_i.vix >= 0 && v_i.vix < (int32_t)mv_size(&gv->vix2var));
        v_i.allele = 0xf; // missing value
        v_i.qual = b_qual;

        char base = seq_nt16_str[b_base];
        int a_i;
        for (a_i = 0; a_i < rec->n_allele; ++a_i){
            if (strlen(rec->d.allele[a_i]) > 1){
                err_msg(0, 0, "seq_base_call_var: indel present %i allele is %s at variant %s\n", 
                        a_i, rec->d.allele[a_i], rec->d.id);
                break;
            }
            if (base == rec->d.allele[a_i][0]){
                v_i.allele = a_i;
                break;
            }
        }

        // add variant
        // TODO: change back to 0,0 
        if (seq_vac_l_insert(vl, v_i, 1, 1) < 0)
            return err_msg(-1, 0, "seq_base_call_var: error adding variant call, %i vars", n_v);
        ++n_added;
    }
    // free allocated var_t objects.
    ml_free(vcfr_list, &vars);

    /* debug */
    /*
    if (1){
        printf("var list:\n");
        ml_node_t(seq_vac_l) *nv;
        for (nv = ml_begin(vl); nv; nv = ml_node_next(nv)){
            seq_vac_t vi = ml_node_val(nv);
            printf("\tvix=%i,allele=%u,qual=%u\n", vi.vix, vi.allele, vi.qual);
        }
    }
    */
    return(n_added);
}

int seq_vac_l_call_var(ml_t(seq_base_l) *bl, ml_t(seq_vac_l) *vl, g_var_t *gv, 
        str_map *cmap, uint8_t min_qual){
    if (bl == NULL || vl == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "seq_vac_l_call_var: arguments cannot be NULL");
    }

    int n_added = 0;
    ml_node_t(seq_base_l) *vn;
    for (vn = ml_begin(bl); vn; vn = ml_node_next(vn)){
        seq_base_t b = ml_node_val(vn);
        int nadd = seq_base_call_var(b, vl, gv, cmap, min_qual);
        if (nadd < 0) return(-1);
        else n_added += nadd;
    }

    return n_added;
}

/*******************************************************************************
 * seq_gene_t
 ******************************************************************************/

void seq_gene_init(seq_gene_t *g){
    // if (g == NULL) return;
    g->gene_id = -1;
    g->splice = SNA; // NA/missing value
}

seq_gene_t *seq_gene_alloc(){
   seq_gene_t *g = malloc(sizeof(seq_gene_t));
   if (g == NULL){
       err_msg(-1, 0, "seq_gene_alloc: %s", strerror(errno));
       return(NULL);
   }
   seq_gene_init(g);
   return(g);
}

void seq_gene_dstry(seq_gene_t *g){
    if (g == NULL) return;
    free(g);
}

int seq_gene_cmp(seq_gene_t g1, seq_gene_t g2){
    int ic = g1.gene_id - g2.gene_id;
    if (ic != 0) return(ic);

    int sc = g1.splice - g2.splice;
    if (sc != 0) return(sc);

    return(0);
}

seq_gene_t *seq_gene_dup(const seq_gene_t *src, int *ret){
    *ret = 0;
    if (src == NULL)
        return(NULL);

    seq_gene_t *cpy = malloc(sizeof(seq_gene_t));
    if (cpy == NULL){
        *ret = err_msg(-1, 0, "seq_gene_dup: %s", strerror(errno));
        return(NULL);
    }

    *cpy = *src;

    return(cpy);
}

int seq_gene_l_cmp(ml_t(seq_gene_l) gl1, ml_t(seq_gene_l) gl2){
    int sc = ml_size(&gl1) - ml_size(&gl2);
    if (sc != 0) return(sc);

    ml_node_t(seq_gene_l) *n1 = ml_begin(&gl1);
    ml_node_t(seq_gene_l) *n2 = ml_begin(&gl1);
    while (n1 != NULL && n2 != NULL){
        seq_gene_t g1 = ml_node_val(n1);
        seq_gene_t g2 = ml_node_val(n2);
        int gc = seq_gene_cmp(g1, g2);
        if (gc != 0) return(gc);
        n1 = ml_node_next(n1);
        n2 = ml_node_next(n2);
    }

    return(0);
}

int seq_gene_l_equal(ml_t(seq_gene_l) gl1, ml_t(seq_gene_l) gl2){
    if (ml_size(&gl1) != ml_size(&gl2))
        return(0);

    ml_node_t(seq_gene_l) *n1 = ml_begin(&gl1);
    ml_node_t(seq_gene_l) *n2 = ml_begin(&gl1);
    while (n1 != NULL && n2 != NULL){
        seq_gene_t g1 = ml_node_val(n1);
        seq_gene_t g2 = ml_node_val(n2);
        if (g1.gene_id != g2.gene_id)
            return(0);
        if (g1.splice != g2.splice)
            return(0);
        n1 = ml_node_next(n1);
        n2 = ml_node_next(n2);
    }
    return(1);
}

