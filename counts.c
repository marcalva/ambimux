
#include <stdlib.h>
#include "str_util.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "counts.h"
#include "variants.h"
#include "region.h"
#include "gtf_anno.h"

/*******************************************************************************
 * seq_base_t
 ******************************************************************************/

void init_seq_base(seq_base_t *sbase){
    init_g_pos(&sbase->pos);
    sbase->base = 15;
    sbase->qual = 0;
    sbase->next = NULL;
}

seq_base_t *alloc_seq_base(){
    seq_base_t *b = (seq_base_t *)calloc(1, sizeof(seq_base_t));
    if (b == NULL){
        err_msg(-1, 0, "alloc_seq_base: %s", strerror(errno));
        return(NULL);
    }

    init_seq_base(b);
    return(b);
}

seq_base_t *dstry_seq_base(seq_base_t *b){
    seq_base_t *n = b->next;
    free(b);
    return(n);
}

seq_base_t *copy_seq_base(const seq_base_t *b, int *ret){
    *ret = 0;
    if (b == NULL) return(NULL);

    seq_base_t *c = (seq_base_t *)calloc(1, sizeof(seq_base_t));
    if (c == NULL){
        *ret = err_msg(-1, 0, "copy_seq_base: %s", strerror(errno));
        return(NULL);
    }

    *c = *b;
    return(c);
}

/*******************************************************************************
 * seq_blist_t
 ******************************************************************************/

void init_seq_blist(seq_blist_t *s){
    s->head = NULL;
    s->n = 0;
}

seq_blist_t *alloc_seq_blist(){
    seq_blist_t *s = (seq_blist_t *)calloc(1, sizeof(seq_blist_t));
    if (s == NULL){
        err_msg(-1, 0, "alloc_seq_blist: %s", strerror(errno));
        return(NULL);
    }

    init_seq_blist(s);
    return(s);
}

void free_seq_blist(seq_blist_t *s){
    if (s == NULL) return;

    seq_base_t *b = s->head;
    while (b != NULL)
        b = dstry_seq_base(b);
    s->head = NULL;
    s->n = 0;
}

int seq_blist_equal(seq_blist_t l1, seq_blist_t l2, int qual_cmp){
    if (l1.n != l2.n) return(0);

    seq_base_t *b1 = l1.head;
    seq_base_t *b2 = l2.head;
    while (b1 != NULL && b2 != NULL){
        if (b1->base != b2->base)
            return(0);
        if (qual_cmp && (b1->qual != b2->qual))
            return(0);
        b1 = b1->next;
        b2 = b2->next;
    }
    return(1);
}

int blist_add_base(seq_blist_t *s, const seq_base_t *b, int dup_ok, int skip_dup){
    if (s == NULL || b == NULL)
        return err_msg(-1, 0, "blist_add_base: arguments cannot be NULL");

    /* sorted insert based on genomic position */
    int cret;
    if (s->head == NULL || poscmp(b->pos, s->head->pos) < 0){
        seq_base_t *c = copy_seq_base(b, &cret);
        if (cret < 0)
            return(-1);
        c->next = s->head;
        s->head = c;
        ++s->n;
        return(2);
    }

    seq_base_t *p = s->head;
    /* b->pos always >= p->pos */
    int pcmp = poscmp(b->pos, p->pos);
    while (p->next != NULL && (pcmp = poscmp(b->pos, p->next->pos)) > 0)
        p = p->next;

    if (pcmp == 0){
        if (dup_ok == 0)
            return err_msg(-1, 0, "blist_add_base: trying to add base at "
                    "duplicate position");

        if (skip_dup == 1)
            return(0);
    }

    seq_base_t *c = copy_seq_base(b, &cret);
    if (cret < 0)
        return(-1);
    c->next = p->next;
    p->next = c;
    ++s->n;

    if (pcmp == 0) return(1);
    else return(2);
}

seq_blist_t *copy_seq_blist(const seq_blist_t *s, int *ret){
    *ret = 0;

    if (s == NULL)
        return(NULL);

    seq_blist_t *c = alloc_seq_blist();
    if (c == NULL){
        *ret = -1;
        return(NULL);
    }

    c->n = s->n;

    if (s->head == NULL) return(c);

    int cret;
    c->head = copy_seq_base(s->head, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    seq_base_t *c_p = c->head, *b_p = s->head->next;
    while (b_p != NULL){
        c_p->next = copy_seq_base(b_p, &cret);
        if (cret < 0){
            *ret = -1;
            return(NULL);
        }
        b_p = b_p->next;
        c_p = c_p->next;
    }

    return(c);
}

int seq_blist_match_qual(seq_blist_t *bl, const seq_blist_t *cmp){
    if (bl == NULL || cmp == NULL)
        return err_msg(-1, 0, "seq_blist_match_qual: arguments are NULL");

    if (bl->n != cmp->n)
        return err_msg(-1, 0, "seq_blist_match_qual: "
                "number of bases don't match");

    seq_base_t *b1, *b2;
    b1 = bl->head;
    b2 = cmp->head;
    while (b1 != NULL){
        if (b2->base != b1->base)
            return err_msg(-1, 0, "seq_blist_match_qual: "
                    "bases don't match between reads");
        if (b2->qual > b1->qual)
            b1->qual = b2->qual;
        b1 = b1->next;
        b2 = b2->next;
    }

    return(0);
}

/*******************************************************************************
 * vac_t
 ******************************************************************************/

void init_vac(vac_t *v){
    v->vix = -1;
    v->allele = 0xf; // missing
    v->qual = 0;
    v->next = NULL;
}

vac_t *alloc_vac(){
    vac_t *v = (vac_t *)calloc(1, sizeof(vac_t));
    if (v == NULL){
        err_msg(-1, 0, "alloc_vac: %s", strerror(errno));
        return NULL;
    }
    init_vac(v);
    return(v);
}

vac_t *dstry_vac(vac_t *v){
    vac_t *n = v->next;
    if (v) free(v); // The bcf record is freed elsewhere
    return(n);
}

void init_vacs(vacs_t *vacs){
    if (vacs == NULL) return;
    vacs->head = NULL;
    vacs->n = 0;
}

void free_vacs(vacs_t *vacs){
    if (vacs == NULL) return;
    vac_t *p = vacs->head;
    while (p != NULL){
        vac_t *n = p->next;
        free(p);
        p = n;
    }
    vacs->head = NULL;
    vacs->n = 0;
}

int vacs_add(vacs_t *vacs, vac_t *v){
    if (vacs == NULL || v == NULL){
        return err_msg(-1, 0, "vacs_add: arguments cannot be NULL");
    }

    if (vacs->head == NULL){
        vacs->head = v;
        vacs->n = 1;
    } else {
        vac_t *p = vacs->head;
        while (p->next != NULL)
            p = p->next;
        p->next = v;
        v->next = NULL;
        vacs->n++;
    }
    return(0);
}

int seq_base_call_var(seq_base_t *b, vacs_t *vacs, GenomeVar *gv, 
        contig_map *cmap, uint8_t min_qual){ 
    if (b == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "seq_base_call_var: arguments are NULL");
    }

    /* get chromosome and position */
    int32_t b_rid = b->pos.rid;
    const char *b_ref = cm_ix_to_chr(cmap, b_rid);
    if (b_ref == NULL)
        return err_msg(-1, 0, "seq_base_call_var: cannot find chromosome ID %i", b_rid);
    hts_pos_t b_beg = b->pos.pos, b_end = b_beg + 1;

    /* base call */
    uint8_t b_base = b->base;
    uint8_t b_qual = b->qual;

    // skip bases observed as N
    if (b_base == 15)
        return(0);

    if (b_qual < min_qual)
        return(0);

    /* Get overlapping variants */
    Var *v_a = NULL, *v_hd = NULL;
    int n_added = 0;
    int n_v = region_vars(gv, b_ref, b_beg, b_end, &v_hd);
    if (n_v < 0) return(-1);
    // debugging
    // if (1){
    if (n_v > 1){
        err_msg(0, 1, "%i variants found at position %s:%"PRIhts_pos"\n", 
                n_v, b_ref, b_beg);
    }
    // }

    /* loop through overlapping variants */
    for (v_a = v_hd; v_a != NULL; v_a = v_a->next){
        bcf1_t *rec = v_a->b;
        // skip if position doesn't overlap
        if (rec->pos != b_beg)
            continue;

        // skip indels
        int issnp = bcf_is_snp(rec);
        if (!issnp)
            continue;

        // throw error if more alleles than present than can be stored.
        // last allele is reserved for missing/other
        if (rec->n_allele > (MAX_ALLELE-1))
            return err_msg(-1, 0, "seq_base_call_var: variant has %u alleles, "
                    "more than %i that is allowed", rec->n_allele, MAX_ALLELE-1);

        // set up vac
        vac_t *vac_next = alloc_vac();
        if (vac_next == NULL) return(-1);

        vac_next->vix = v_a->vix;
        vac_next->allele = 0xf; // missing value
        vac_next->qual = b_qual;

        char base = seq_nt16_str[b_base];
        int a_i;
        for (a_i = 0; a_i < rec->n_allele; ++a_i){
            if (strlen(rec->d.allele[a_i]) > 1){
                err_msg(0, 0, "indel present %i allele is %s at variant %s\n", 
                        a_i, rec->d.allele[a_i], rec->d.id);
                break;
            }
            if (base == rec->d.allele[a_i][0]){
                vac_next->allele = a_i;
                break;
            }
        }

        // add variant
        if (vacs_add(vacs, vac_next) < 0) return(-1);
        ++n_added;
    }
    // free allocated Var objects.
    v_a = v_hd;
    while (v_a != NULL){
        Var *v_tmp = v_a->next;
        free(v_a);
        v_a = v_tmp;
    }
    return(n_added);
}

int seq_blist_call_var(seq_blist_t *s, vacs_t *vacs, GenomeVar *gv, 
        contig_map *cmap, uint8_t min_qual){
    if (s == NULL || vacs == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "pbcs_to_vacs: arguments cannot be NULL");
    }

    if (s->n == 0 || s->head == NULL)
        return(0);

    int n_added = 0;
    seq_base_t *b = s->head;
    for (; b != NULL; b = b->next){
        int nadd = seq_base_call_var(b, vacs, gv, cmap, min_qual);
        if (nadd < 0) return(-1);
        else n_added += nadd;
    }

    return n_added;
}

/*******************************************************************************
 * seq_gene_t
 ******************************************************************************/

void seq_gene_init(seq_gene_t *g){
    if (g == NULL) return;
    g->gene_id = -1;
    g->splice = SNA; // NA/missing value
    g->next = NULL;
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

seq_gene_t *seq_gene_dstry(seq_gene_t *g){
    if (g == NULL) return(NULL);

    seq_gene_t *n = g->next;
    free(g);
    return(n);
}

seq_gene_t *seq_gene_cpy(const seq_gene_t *src, int *ret){
    *ret = 0;
    if (src == NULL)
        return(NULL);

    seq_gene_t *cpy = malloc(sizeof(seq_gene_t));
    if (cpy == NULL){
        *ret = err_msg(-1, 0, "seq_gene_cpy: %s", strerror(errno));
        return(NULL);
    }

    *cpy = *src;
    cpy->next = NULL;

    return(cpy);
}

/*******************************************************************************
 * seq_glist_t
 ******************************************************************************/

void seq_glist_init(seq_glist_t *gl){
    if (gl == NULL) return;

    gl->head = NULL;
    gl->n = 0;
}

seq_glist_t *seq_glist_alloc(){
    seq_glist_t *gl = (seq_glist_t *)calloc(1, sizeof(seq_glist_t));
    if (gl == NULL){
        err_msg(-1, 0, "seq_glist_alloc: %s", strerror(errno));
        return(NULL);
    }
    seq_glist_init(gl);
    return(gl);
}

void seq_glist_free(seq_glist_t *gl){
    if (gl == NULL) return;

    seq_gene_t *g = gl->head;
    while (g != NULL)
        g = seq_gene_dstry(g);

    gl->head = NULL;
    gl->n = 0;
}

seq_glist_t *seq_glist_cpy(const seq_glist_t *src, int *ret){
    *ret = 0;
    if (src == NULL)
        return(NULL);

    seq_glist_t *cpy = seq_glist_alloc();
    if (cpy == NULL){
        *ret = -1;
        return(NULL);
    }

    cpy->n = src->n;
    if (cpy->n == 0)
        return(cpy);

    int cret = 0;
    cpy->head = seq_gene_cpy(src->head, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }

    seq_gene_t *p_src = src->head->next, *p_dst = cpy->head;
    while (p_src != NULL){
        p_dst->next = seq_gene_cpy(p_src, &cret);
        if (cret < 0){
            *ret = -1;
            return(NULL);
        }
        p_src = p_src->next;
        p_dst = p_dst->next;
    }

    return(cpy);
}

int seq_glist_add_gene(seq_glist_t *gl, const seq_gene_t *gene){
    if (gl == NULL || gene == NULL)
        return err_msg(-1, 0, "seq_gene_add_gene: arguments are NULL");

    // copy seq_gene_t
    int cret = 0;
    seq_gene_t *cpy = seq_gene_cpy(gene, &cret);
    if (cret < 0)
        return(-1);

    // add seq_gene
    if (gl->head == NULL || cpy->gene_id < gl->head->gene_id){
        cpy->next = gl->head;
        gl->head = cpy;
    } else {
        seq_gene_t *p = gl->head;
        while (p->next != NULL && cpy->gene_id >= p->next->gene_id)
            p = p->next;

        cpy->next = p->next;
        p->next = cpy;
    }
    ++gl->n;

    return(0);
}

