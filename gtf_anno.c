
#include "gtf_anno.h"
#include "bins.h"
#include "overlap.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>

void exon_init(exon_t *e){
    if (e == NULL) return;
    e->beg = -1;
    e->end = -1;
}

int isoform_init(isoform_t *iso){
    iso->id = NULL;
    iso->beg = 0;
    iso->end = 0;
    ml_init(exon_list, &iso->exons);
    ml_init(exon_list, &iso->utrs);
    return 0;
}

void isoform_free(isoform_t *iso){
    if (iso == NULL) return;
    free(iso->id);
    ml_free(exon_list, &iso->exons);
    ml_free(exon_list, &iso->utrs);
}

void gene_init(gene_t *gene){
    if (gene == NULL) return;

    gene->id = NULL;
    gene->name = NULL;
    gene->type = NULL;
    gene->beg = 0;
    gene->end = 0;
    gene->strand = 0;
    gene->chrm = -1;
    gene->bin = -1;
    gene->bt_isoforms = kb_init(kb_iso, KB_DEFAULT_SIZE);
    gene->isoforms_n = 0;
}

void gene_free(gene_t *gene){
    if (gene == NULL) return;

    free(gene->id);
    gene->id = NULL;
    free(gene->name);
    gene->name = NULL;
    free(gene->type);
    gene->type = NULL;

    kbtree_t(kb_iso) *bt = gene->bt_isoforms;
    kbitr_t itr;
    kb_itr_first(kb_iso, bt, &itr);
    for (; kb_itr_valid(&itr); kb_itr_next(kb_iso, bt, &itr)){
        isoform_t *iso = &kb_itr_key(isoform_t, &itr);
        isoform_free(iso);
    }
    kb_destroy(kb_iso, gene->bt_isoforms);
    gene->bt_isoforms = NULL;
}

void chr_genes_init(chr_genes_t *cg){
    if (cg == NULL) return;

    int i;
    for (i = 0; i < MAX_BIN; i++){
        ml_init(gl, &cg->bins[i]);
    }
}

chr_genes_t *chr_genes_alloc(){
    chr_genes_t *cg = (chr_genes_t*)calloc(1, sizeof(chr_genes_t));

    if (cg == NULL){
        err_msg(-1, 0, "init_chr_genes: %s", strerror(errno));
        return NULL;
    }

    int i;
    for (i = 0; i < MAX_BIN; i++) ml_init(gl, &cg->bins[i]);

    return cg;
}

void chr_genes_free(chr_genes_t *cg){
    if (cg == NULL) return;

    int i;
    for (i = 0; i < MAX_BIN; i++){
        ml_node_t(gl) *gn;
        for (gn = ml_begin(&cg->bins[i]); gn; gn = ml_node_next(gn)){
            gene_t *gene = ml_node_val(gn);
            gene_free(gene);
            free(gene);
        }
        ml_free(gl, &cg->bins[i]);
    }

}

int gene_anno_init(gene_anno_t *gene_anno){
    if (gene_anno == NULL) return(0);

    mv_init(&gene_anno->chrms);

    gene_anno->chrm_ix = init_str_map();
    if (gene_anno->chrm_ix == NULL)
        return err_msg(-1, 0, "gene_anno_init: %s", strerror(errno));

    gene_anno->gene_ix = init_str_map();
    if (gene_anno->gene_ix == NULL)
        return err_msg(-1, 0, "gene_anno_init: %s", strerror(errno));

    mv_init(&gene_anno->gix2gene);
    
    return(0);
}

gene_anno_t *gene_anno_alloc(){
    gene_anno_t *gene_anno = malloc(sizeof(gene_anno_t));
    if (gene_anno == NULL){
        err_msg(-1, 0, "gene_anno_alloc: %s", strerror(errno));
        return NULL;
    }
    if (gene_anno_init(gene_anno) < 0) return(NULL);
    return gene_anno;
}

void gene_anno_free(gene_anno_t *gene_anno){
    if (gene_anno == NULL) return;

    mv_free(&gene_anno->gix2gene);

    size_t i;
    for (i = 0; i < mv_size(&gene_anno->chrms); i++){
        chr_genes_t *c = mv_i(&gene_anno->chrms, i);
        chr_genes_free(c);
        free(c);
    }
    mv_free(&gene_anno->chrms);

    destroy_str_map(gene_anno->chrm_ix);
    gene_anno->chrm_ix = NULL;
    destroy_str_map(gene_anno->gene_ix);
    gene_anno->gene_ix = NULL;

}

void gene_anno_dstry(gene_anno_t *gene_anno){
    if (gene_anno == NULL) return;

    gene_anno_free(gene_anno);
    free(gene_anno);
}


int add_chrm(gene_anno_t *a, char *c){
    if (a == NULL || c == NULL)
        return err_msg(-1, 0, "add_chrm: argument is null");

    // add chromosome ID
    int chr_ix, found = 0;
    if ( (chr_ix = add2str_map(a->chrm_ix, (const char*)c, &found)) < 0 ) return -1;
    if (found == 0){
        chr_genes_t *p = chr_genes_alloc();
        if (p == NULL)
            return(-1);

        if (mv_push(cv, &a->chrms, p) < 0)
            return err_msg(-1, 0, "add_chrm: failed to add chromosome to index");
    }
    return chr_ix;
}

gene_t *gene_anno_get_gene(gene_anno_t *gene_anno, char *gene_id){
    if (gene_anno == NULL || gene_id == NULL){
        err_msg(-1, 0, "gene_anno_get_gene: argument is null");
        return(NULL);
    }

    int gix = str_map_ix(gene_anno->gene_ix, gene_id);
    if (gix < 0){
        err_msg(-1, 0, "gene_anno_get_gene: gene %s not found in index", 
                gene_id);
        return(NULL);
    }
    assert((size_t)gix < mv_size(&gene_anno->gix2gene));

    gene_t *p = mv_i(&gene_anno->gix2gene, gix);

    return(p);
}

int isoform_collapse_utrs(isoform_t *iso){
    if (iso == NULL)
        return err_msg(-1, 0, "isoform_collapse_utrs: argument is null");

    // Add UTRs as exons
    ml_node_t(exon_list) *n_utrs, *n_exons, *exon_next;
    for (n_utrs = ml_begin(&iso->utrs); n_utrs; n_utrs = ml_node_next(n_utrs)){
        exon_t utr = ml_node_val(n_utrs);
        int ex_hit = 0;
        for (n_exons = ml_begin(&iso->exons); n_exons; n_exons = ml_node_next(n_exons)){
            exon_t exon = ml_node_val(n_exons);
            if (utr.beg < exon.beg && utr.end >= exon.beg){
                ml_node_val(n_exons).beg = utr.beg;
                ex_hit = 1;
            }
            if (utr.beg <= exon.end && utr.end > exon.end){
                ml_node_val(n_exons).end = utr.end;
                ex_hit = 1;
            }
            // exons are sorted by pos, so break if exons passed UTR
            if (utr.end < exon.beg)
                break;
        }
        if (ex_hit == 0){
            if (ml_insert(exon_list, &iso->exons, utr, 0, 0) < 0)
                return err_msg(-1, 0, "isoform_collapse_utrs: failed to add UTR to exon list");
        }
    }

    //
    n_exons = ml_begin(&iso->exons);
    exon_next = ml_node_next(n_exons);
    while ((exon_next = ml_node_next(n_exons)) != NULL){
        exon_t exon1 = ml_node_val(n_exons);
        exon_t exon2 = ml_node_val(exon_next);
        if (exon1.end >= exon2.beg){
            int end = exon2.end > exon1.end ? exon2.end : exon1.end;
            ml_node_val(n_exons).end = end;
            ml_node_next(n_exons) = ml_node_next(exon_next);
            free(exon_next);
        } else {
            n_exons = ml_node_next(n_exons);
        }
    }
    return(0);
}

int gene_anno_collapse_utrs(gene_anno_t *a){
    if (a == NULL)
        return err_msg(-1, 0, "gene_anno_collapse_utrs: argument is null");

    int nc = a->chrm_ix->n;
    int i, j;
    for (i = 0; i < nc; i++){
        chr_genes_t *cg = mv_i(&a->chrms, i);
        for (j = 0; j < MAX_BIN; j++){
            ml_node_t(gl) *gn;
            for (gn = ml_begin(&cg->bins[j]); gn; gn = ml_node_next(gn)){
                gene_t *g = ml_node_val(gn);
                kbtree_t(kb_iso) *bt = g->bt_isoforms;
                kbitr_t itr;
                kb_itr_first(kb_iso, bt, &itr);
                for (; kb_itr_valid(&itr); kb_itr_next(kb_iso, bt, &itr)){
                    isoform_t *iso = &kb_itr_key(isoform_t, &itr);
                    if (isoform_collapse_utrs(iso) < 0)
                        return(-1);
                }
            }
        }
    }
    return 0;
}

/****************************
 * GTF processing
 ****************************/

gtf1_t *gtf1_alloc(){
    gtf1_t *g = (gtf1_t *)calloc(1, sizeof(gtf1_t));

    if (g == NULL){
        err_msg(-1, 0, "init_gtf1_t: %s", strerror(errno));
        return NULL;
    }

    ks_initialize( &(g->chrname) );
    ks_initialize( &(g->source) );
    ks_initialize( &(g->feature) );
    g->start = -1;
    g->end = -1;
    g->score = -1;
    g->strand = '.';
    g->frame = -1;
    ks_initialize( &(g->attribute) );
    g->attr_tag = NULL;
    g->attr_val = NULL;
    g->n_attr = 0;

    return g;
}

void gtf1_free(gtf1_t *g){
    if (g == NULL) return;

    ks_free(&(g->chrname));
    ks_free(&(g->source));
    ks_free(&(g->feature));

    /* free attributes data */
    ks_free(&(g->attribute));
    kstr_node *n;
    n = g->attr_tag;
    while (n) n = destroy_kstr_node(n);
    n = g->attr_val;
    while (n) n = destroy_kstr_node(n);
    /* */

    g->start = -1;
    g->end = -1;
    g->score = -1;
    g->strand = '.';
    g->frame = -1;
    g->attr_tag = NULL;
    g->attr_val = NULL;
    g->n_attr = 0;
}

void gtf1_dstry(gtf1_t *g){
    if (g == NULL) return;
    gtf1_free(g);
    free(g);
}

// Add gene from a GTF line
int gtf_gene2anno(gene_anno_t *a, gtf1_t *gl){
    if (a == NULL || gl == NULL)
        return err_msg(-1, 0, "gtf_gene2anno: argument is null");

    int ci = add_chrm(a, gl->chrname.s);
    if (ci < 0) return -1;

    char *gene_id = get_attr_val(gl, GENE_ID);
    if (gene_id == NULL)
        return err_msg(-1, 0, "gtf_gene2anno: gtf line must have %s attribute", GENE_ID);

    char *gene_name = get_attr_val(gl, GENE_NAME);
    if (gene_name == NULL)
        return err_msg(-1, 0, "gtf_gene2anno: gtf line must have %s attribute", GENE_NAME);

    char *gene_type = get_attr_val(gl, GENE_TYPE);
    if (gene_type == NULL)
        return err_msg(-1, 0, "gtf_gene2anno: gtf line must have %s attribute", GENE_TYPE);

    // sanity check
    if (gl->end <= gl->start){
        return err_msg(-1, 0, "gtf_gene2anno: end <= beg pos for gene "
                "%s", gene_name);
    }

    int gix, found = 0;
    if ( (gix = add2str_map(a->gene_ix, (const char*)gene_id, &found)) < 0 ) return -1;
    if (found)
        return err_msg(-1, 0, "gtf_gene2anno: gene %s found twice", gene_id);

    gene_t *gene = calloc(1, sizeof(gene_t));
    if (gene == NULL)
        return err_msg(-1, 0, "gtf_gene2anno: %s", strerror(errno));
    gene_init(gene);

    gene->id = strdup(gene_id);
    if (gene->id == NULL){
        gene_free(gene);
        return err_msg(-1, 0, "gtf_gene2anno: %s", strerror(errno));
    }

    gene->name = strdup(gene_name);
    if (gene->name == NULL){
        gene_free(gene);
        return err_msg(-1, 0, "gtf_gene2anno: %s", strerror(errno));
    }

    gene->type = strdup(gene_type);
    if (gene->type == NULL){
        gene_free(gene);
        return err_msg(-1, 0, "gtf_gene2anno: %s", strerror(errno));
    }

    gene->beg = gl->start;
    gene->end = gl->end;
    gene->strand = gl->strand;
    gene->chrm = ci;
    gene->bin = reg2bin(gene->beg, gene->end);

    // add to chr bins
    chr_genes_t *cg = mv_i(&a->chrms, ci);
    if (ml_insert(gl, &cg->bins[gene->bin], gene, 0, 0) < 0)
        return err_msg(-1, 0, "gtf_gene2anno: failed to add gene to bin list");

    // add gene pointer to gix2gene
    if (mv_push(gv, &a->gix2gene, gene) < 0)
        return err_msg(-1, 0, "gtf_gene2anno: failed to add gene to index");

    assert(mv_size(&a->gix2gene) == (size_t)a->gene_ix->n);

    return 0;
}

// add isoform from a GTF line
int gtf_iso2anno(gene_anno_t *a, gtf1_t *gl){
    int ci = add_chrm(a, gl->chrname.s);
    if (ci < 0) return -1;

    // get GENE ID
    char *gene_id = get_attr_val(gl, GENE_ID);
    if (!gene_id)
        return err_msg(-1, 0, "gtf_iso2anno: gtf line must have %s attribute", GENE_ID);

    // get TX ID
    char *tx_id = get_attr_val(gl, TX_ID);
    if (!tx_id)
        return err_msg(-1, 0, "gtf_iso2anno: gtf line for gene %s must have %s attribute", gene_id, TX_ID);
    
    // sanity check
    if (gl->end <= gl->start){
        return err_msg(-1, 0, "gtf_iso2anno: end <= beg pos for isoform "
                "%s %s", gene_id, tx_id);
    }

    isoform_t iso;
    if (isoform_init(&iso) < 0)
        return -1;
    
    iso.id = strdup(tx_id);
    iso.beg = gl->start;
    iso.end = gl->end;

    // get gene object
    gene_t *ag = gene_anno_get_gene(a, gene_id);
    if (ag == NULL){
        err_msg(-1, 0, "gtf_iso2anno: gene %s not found", gene_id);
        return -1;
    }

    // add isoform to btree
    isoform_t *p;
    p = kb_getp(kb_iso, ag->bt_isoforms, &iso);
    if (p == NULL)
        kb_putp(kb_iso, ag->bt_isoforms, &iso);
    else
        return err_msg(-1, 0, "gtf_iso2anno: trying to add duplicate isoform %s", tx_id);

    ++ag->isoforms_n;

    return 0;
}

// add exon from a GTF line
int gtf_exon2anno(gene_anno_t *a, gtf1_t *gl, int is_utr){
    int ci = add_chrm(a, gl->chrname.s);
    if (ci < 0) return -1;

    exon_t ex;
    exon_init(&ex);
    
    // get gene tx exon IDs
    char *gene_id = get_attr_val(gl, GENE_ID);
    if (gene_id == NULL){
        err_msg(-1, 0, "gtf_exon2anno: gtf line must have %s attribute", GENE_ID);
        return -1;
    }
    char *tx_id = get_attr_val(gl, TX_ID);
    if (tx_id == NULL){
        err_msg(-1, 0, "gtf_exon2anno: gtf line must have %s attribute", TX_ID);
        return -1;
    }
    char *exon_id = get_attr_val(gl, EXON_ID);
    if (exon_id == NULL){
        err_msg(-1, 0, "gtf_exon2anno: gtf line must have %s attribute", EXON_ID);
        return -1;
    }
    
    // sanity check
    if (gl->end <= gl->start){
        return err_msg(-1, 0, "gtf_exon2anno: end <= beg pos for exon "
                "%s %s %s", gene_id, tx_id, exon_id);
    }

    // fprintf(stderr, "start = %i\n", gl->start);
    ex.beg = gl->start;
    ex.end = gl->end;
    // fprintf(stderr, "ex start = %i\n", ex.beg);

    // get gene object
    gene_t *ag = gene_anno_get_gene(a, gene_id);
    if (ag == NULL){
        err_msg(-1, 0, "gtf_iso2anno: gene %s not found", gene_id);
        return -1;
    }

    // get isoform object
    isoform_t *iso, d_iso;
    d_iso.id = tx_id;
    iso = kb_getp(kb_iso, ag->bt_isoforms, &d_iso);
    if (iso == NULL)
        return err_msg(-1, 0, "gtf_exon2anno: could not find isoform %s to "
                "add exon %s to. GTF file is likely malformed", tx_id, exon_id);

    // make sure it's sorted
    if (is_utr == 0){
        if (ml_insert(exon_list, &iso->exons, ex, 0, 0) < 0)
            return err_msg(-1, 0, "gtf_exon2anno: could not add exon %s to isoform list", 
                    exon_id);
    } else {
        if (ml_insert(exon_list, &iso->utrs, ex, 0, 0) < 0)
            return err_msg(-1, 0, "gtf_exon2anno: could not add UTR %s to isoform list", 
                    exon_id);
    }

    // implement check to see if exons overlap

    return 0;
}

int gtf1_to_anno(gene_anno_t *a, gtf1_t *gl){

    int ret;
    char *type = gl->feature.s;

    if (strcmp(type, GENE) == 0)
        ret = gtf_gene2anno(a, gl);
    else if (strcmp(type, TX) == 0)
        ret = gtf_iso2anno(a, gl);
    else if (strcmp(type, EXON) == 0)
        ret = gtf_exon2anno(a, gl, 0);
    else if (strcmp(type, UTR) == 0)
        ret = gtf_exon2anno(a, gl, 1);
    else
        return 0;

    return ret;
}

int parse_gtf1(kstring_t *line, gtf1_t *g){
    char delims[] = "\t";
    char *token = NULL;
    char *rest = NULL;
    char *endptr = NULL;

    int save_errno = errno;
    errno = 0;

    if ( (token = strtok_r(line->s, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf1: GTF line is malformed\n%s", line);
    if (kputs(token, &(g->chrname)) == EOF) return -1;
    
    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf1: GTF line is malformed\n%s", line);
    if (kputs(token, &(g->source)) == EOF) return -1;

    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf1: GTF line is malformed\n%s", line);
    if (kputs(token, &(g->feature)) == EOF) return -1;
    
    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf1: GTF line is malformed\n%s", line);
    g->start = (int) strtol(token, &endptr, 10);
    if (g->start == 0 && errno > 0){
        return err_msg(-1, 0, "parse_gtf1: could not convert %s to int: %s", 
                token, strerror(errno));
    }
    g->start -= 1; // convert start to 0-based coordinate
    
    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf1: GTF line is malformed\n%s", line);
    g->end = strtol(token, &endptr, 10);
    if (g->end == 0 && errno > 0){
        return err_msg(-1, 0, "parse_gtf1: could not convert %s to int: %s", 
                token, strerror(errno));
    }

    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf1: GTF line is malformed\n%s", line);
    if (token[0] == '.')
        g->score = -1;
    else{
        g->score = strtol(token, &endptr, 10);
        if (g->score == 0 && errno > 0)
            return err_msg(-1, 0, "parse_gtf1: could not convert %s to int: %s", 
                    token, strerror(errno));
    }

    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf1: GTF line is malformed\n%s", line);
    g->strand = token[0];
    // can be '+' '-' '.':irrelevant '?':relevant but unknown

    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf1: GTF line is malformed\n%s", line);
    if (token[0] == '.')
        g->frame = -1;
    else{
        g->frame = strtol(token, &endptr, 10);
        if (g->frame == 0 && errno > 0)
            return err_msg(-1, 0, "parse_gtf1: could not convert %s to int: %s",
                    token, strerror(errno));
    }

    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf1: GTF line is malformed\n%s", line);;
    if (kputs(token, &(g->attribute)) == EOF) return -1;

    errno = save_errno;

    return 0;
}

int parse_gtf_attr(gtf1_t *g){
    kstring_t *s = &(g->attribute);
    char delim[4] = " ;\"";

    char *token;
    char *rest;
    token = strtok_r(s->s, delim, &rest);
    while (token){
        kstr_node *tag_node = init_kstr_node();
        kstr_node *val_node = init_kstr_node();
        if (tag_node == NULL || val_node == NULL)
            return err_msg(-1, 0, "parse_gtf_attr: %s");

        kputs(token, &(tag_node->str));
        token = strtok_r(NULL, delim, &rest);
        kputs(token, &(val_node->str));
        token = strtok_r(NULL, delim, &rest);

        if (g->n_attr == 0){
            g->attr_tag = tag_node;
            g->attr_val = val_node;
        }
        else {
            kstr_node *tag_next = g->attr_tag;
            while (tag_next->next) tag_next = tag_next->next;
            kstr_node *val_next = g->attr_val;
            while (val_next->next) val_next = val_next->next;
            tag_next->next = tag_node;
            val_next->next = val_node;
        }
        g->n_attr++;
    }
    return 0;
}

// returns char* pointer to value of first occurence of key.
char *get_attr_val(gtf1_t *g, char *key){
    kstr_node *tag_next = g->attr_tag;
    kstr_node *val_next = g->attr_val;
    int i;
    for (i = 0; i < g->n_attr; i++){
        if (strcmp(key, (tag_next->str).s) == 0)
            return (val_next->str).s;
        tag_next = tag_next->next;
        val_next = val_next->next;
    }
    return NULL;
}   

int has_key_val(gtf1_t *g, char *key, char *val){
    kstr_node *tag_next = g->attr_tag;
    kstr_node *val_next = g->attr_val;
    int i;
    for (i = 0; i < g->n_attr; i++){
        if ((strcmp(key, (tag_next->str).s) == 0) && 
                (strcmp(val, (val_next->str).s) == 0))
            return 1;
        tag_next = tag_next->next;
        val_next = val_next->next;
    }
    return 0;
}

gene_anno_t *read_from_gtf(char *file, int basic){

    gene_anno_t *a = gene_anno_alloc();

    if (a == NULL) return NULL;

    kstring_t line = KS_INITIALIZE;
    int len = 0;

    const char mode_read[3] = "r1\0";
    BGZF *fp = (BGZF *)bgzf_open(file, mode_read);
    if (fp == 0){
        err_msg(-1, 0, "read_from_gtf: could not open GTF file %s", file);
        gene_anno_free(a);
        free(a);
        bgzf_close(fp);
        return NULL;
    }

    // Read out header
    char pr;
    pr = bgzf_peek(fp);
    while (pr == '#'){
        len = bgzf_getline(fp, '\n', &line);
        pr = bgzf_peek(fp);
    }

    int nlines = 0;
    gtf1_t *gtf1 = gtf1_alloc();
    while ( (len = bgzf_getline(fp, '\n', &line)) > 0 ){
        nlines++;

        if (parse_gtf1(&line, gtf1) < 0){
            err_msg(-1, 0, "read_from_gtf: failed to parse GTF line\n%s", line.s);
            gtf1_dstry(gtf1);
            gene_anno_free(a); free(a);
            bgzf_close(fp);
            return NULL;
        }
        if (parse_gtf_attr(gtf1) < 0){
            err_msg(-1, 0, "read_from_gtf: failed to parse GTF line\n%s", line.s);
            gtf1_dstry(gtf1);
            gene_anno_free(a); free(a);
            bgzf_close(fp);
            return NULL;
        }

        // filter basic tag from tx or exon gtf line
        int line_is_basic = has_key_val(gtf1, "tag", "basic");
        int basic_flt = (!line_is_basic) && basic;
        if ((strcmp(gtf1->feature.s, GENE) != 0) && // not gene
            basic_flt){                           // filter for basic tx
            gtf1_free(gtf1);
            continue;
        }

        // add gtf line
        if (gtf1_to_anno(a, gtf1) < 0){
            err_msg(-1, 0, "read_from_gtf: failed to add GTF line\n%s", line.s);
            gtf1_dstry(gtf1);
            gene_anno_free(a);
            free(a);
            bgzf_close(fp);
            return NULL;
        }

        gtf1_free(gtf1);

        // printf("Feature: name = %s; id = %s; beg = %i; end = %i; strand = %c; chrm = %i; bin = %i\n", f.name.s, f.id.s, f.beg, f.end, f.strand, f.chrm, f.bin);
    }

    ks_free(&line);
    gtf1_dstry(gtf1);
    bgzf_close(fp);

    if (gene_anno_collapse_utrs(a) < 0){
        gene_anno_free(a);
        free(a);
        return(NULL);
    }

    return a;
}

int write_gene_data(BGZF *fp, gene_anno_t *anno, str_map *gene_ix){
    if (fp == NULL || anno == NULL || gene_ix == NULL)
        return err_msg(-1, 0, "write_gene_data: arguments are NULL");

    int il; // intstrp string length
    size_t intstrp_m = 1; // intstrp allocated size
    char *intstrp = malloc(sizeof(char) * intstrp_m);
    int ret, k;
    for (k = 0; k < gene_ix->n; ++k){
        char *gene_key = str_map_str(gene_ix, k);
        if (gene_key == NULL){
            fprintf(stderr, "gene index %i not found\n", k);
            continue;
        }

        gene_t *gene_obj = gene_anno_get_gene(anno, gene_key);

        char *chrm = str_map_str(anno->chrm_ix, gene_obj->chrm);
        if (chrm == NULL){
            fprintf(stderr, "chromosome index %i not found\n", gene_obj->chrm);
            continue;
        }
        ret = bgzf_write(fp, chrm, strlen(chrm));

        if ((il = int2strp(gene_obj->beg+1, &intstrp, &intstrp_m)) < 0) return -1;
        ret = bgzf_write(fp, "\t", 1);
        ret = bgzf_write(fp, intstrp, il);

        if ((il = int2strp(gene_obj->end, &intstrp, &intstrp_m)) < 0) return -1;
        ret = bgzf_write(fp, "\t", 1);
        ret = bgzf_write(fp, intstrp, il);

        ret = bgzf_write(fp, "\t", 1);
        ret = bgzf_write(fp, &gene_obj->strand, 1);

        ret = bgzf_write(fp, "\t", 1);
        ret = bgzf_write(fp, gene_obj->type, strlen(gene_obj->type));

        ret = bgzf_write(fp, "\t", 1);
        ret = bgzf_write(fp, gene_obj->name, strlen(gene_obj->name));

        ret = bgzf_write(fp, "\t", 1);
        ret = bgzf_write(fp, gene_key, strlen(gene_key));

        ret = bgzf_write(fp, "\n", 1);

        if (ret < 0)
            return err_msg(-1, 0, "write_gene_data: could not write to file");
    }
    free(intstrp);

    return(0);
}

/****************************
 * Region overlap
 *****************************/

int feats_from_region_p(const gene_anno_t *a, const char* ref, 
        int32_t beg, int32_t end, uint8_t stranded, char strand, 
        ml_t(gl) *genes, double p){
    if (a == NULL || ref == NULL || genes == NULL)
        return err_msg(-1, 0, "feats_from_region_p: argument is null");

    // overlap chromosomes
    int tid = str_map_ix(a->chrm_ix, (char *)ref);
    if (tid < 0) 
        return 0;

    // region length
    double reg_len = (double)end - (double)beg;
    if (reg_len < 0) 
        return err_msg(-1, 0, "feats_from_region_p: end (%li) < beg (%li)", (int32_t)end, (int32_t)beg);

    int ngenes = 0;

    // get bins for query region
    uint16_t list[MAX_BIN];
    int n_bin = reg2bins((int)beg, (int)end, list);

    ml_t(gl) *chr_bins = mv_i(&a->chrms, tid)->bins;

    int i;
    for (i = 0; i < n_bin; ++i){ // for each bin
        uint16_t bin = list[i];
        ml_node_t(gl) *gn;
        for (gn = ml_begin(&chr_bins[bin]); gn; gn = ml_node_next(gn)){
            gene_t *g = ml_node_val(gn);
            char reg_strand = '.';
            if (stranded) reg_strand = strand;

            int ovrlp = bp_overlap(beg, end, reg_strand, g->beg, g->end, g->strand);

            if (ovrlp < 0)
                return err_msg(-1, 0, "feats_from_region_p: failed to get bp_overlap");

            // fraction overlap of region
            double frac_ovrlp = (double)ovrlp / reg_len;
            if (frac_ovrlp < p) continue;

            if (g->id == NULL) // if gene has no ID
                return err_msg(-1, 0, "feats_from_region_p: overlapping genes must have IDs. Beg=%i End=%i", 
                        g->beg, g->end);

            if (ml_insert(gl, genes, g, 0, 0) < 0)
                return err_msg(-1, 0, "feats_from_region_p: could not add gene to list");

            ++ngenes;
        }
    }
    return ngenes;
}


/****************************
 * summary functions
 ****************************/

int n_feat(gene_anno_t *a, int *n_gene, int *n_iso, int *n_exon){
    *n_gene = 0;
    *n_iso = 0;
    *n_exon = 0;
    int nc = a->chrm_ix->n;
    int i, j;
    for (i = 0; i < nc; i++){
        chr_genes_t *cg = mv_i(&a->chrms, i);
        for (j = 0; j < MAX_BIN; j++){
            ml_node_t(gl) *gn;
            for (gn = ml_begin(&cg->bins[j]); gn; gn = ml_node_next(gn)){
                gene_t *g = ml_node_val(gn);
                *n_gene = *n_gene + 1;
                fprintf(stdout, "gene %s has %i isoforms\n", g->id, g->isoforms_n);
                kbtree_t(kb_iso) *bt = g->bt_isoforms;
                kbitr_t itr;
                kb_itr_first(kb_iso, bt, &itr);
                for (; kb_itr_valid(&itr); kb_itr_next(kb_iso, bt, &itr)){
                    isoform_t *iso = &kb_itr_key(isoform_t, &itr);
                    fprintf(stdout, "ID %s beg %i end %i n_exon %zu\n", iso->id, 
                            iso->beg, iso->end, ml_size(&iso->exons));
                    *n_iso += 1;
                    *n_exon += + ml_size(&iso->exons);
                }
            }
        }
    }
    return nc;
}

int test_anno(char *file){
    fprintf(stdout, "Reading from GTF\n");
    gene_anno_t *a = read_from_gtf(file, 1);
    fprintf(stdout, "finished\n");
    int n_gene = 0;
    int n_chr = 0;
    int n_iso = 0;
    int n_exon = 0;
    n_chr = n_feat(a, &n_gene, &n_iso, &n_exon);
    fprintf(stdout, "Number of chromosomes: %i\n", n_chr);
    fprintf(stdout, "Number of genes: %i\n", n_gene);
    fprintf(stdout, "Number of isoforms: %i\n", n_iso);
    fprintf(stdout, "Number of exons: %i\n", n_exon);
    gene_anno_free(a);
    free(a);
    return 0;
}

