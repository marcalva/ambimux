
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

Exon *init_exon(){
    Exon *e;
    
    if ((e = malloc(sizeof(Exon))) == NULL){
        err_msg(-1, 0, "init_exon: %s", strerror(errno));
    } else {
        e->beg = -1;
        e->end = -1;
        e->next = NULL;
    }

    return e;
}

Exon *destroy_exon(Exon *e){
    if (e == NULL) return NULL;
    Exon *n = e->next;
    free(e);
    return n;
}

Isoform *init_isoform(){
    Isoform *i;

    if ((i = malloc(sizeof(Isoform))) == NULL){
        err_msg(-1, 0, "init_isoform: %s", strerror(errno));
    } else {
        i->beg = 0;
        i->end = 0;
        i->exons = init_exon(); // dummy head node
        i->exons_n = 0;
    }

    return i;
}

void destroy_isoform(Isoform *iso){
    if (iso == NULL) return;
    Exon *e = iso->exons;
    while (e) e = destroy_exon(e);
    free(iso);
}

Gene *init_gene(){
    Gene *g;

    if ((g = malloc(sizeof(Gene))) == NULL){
        err_msg(-1, 0, "init_gene: %s", strerror(errno));
    } else {
        g->id = NULL;
        g->name = NULL;
        g->type = NULL;
        g->beg = 0;
        g->end = 0;
        g->strand = 0;
        g->chrm = -1;
        g->bin = -1;
        g->isoforms = kh_init(iso);
        g->isoforms_n = 0;
        g->next = NULL;
    }
    
    return g;
}

Gene *destroy_gene(Gene *g){
    khint_t k;
    for (k = kh_begin(g->isoforms); k != kh_end(g->isoforms); ++k){
        if (!kh_exist(g->isoforms, k)) continue;
        Isoform *iso = kh_val(g->isoforms, k);
        char *key = kh_key(g->isoforms, k);
        destroy_isoform(iso);
        free(key);
    }
    kh_destroy(iso, g->isoforms);
    free(g->name);
    free(g->type);
    Gene *n = g->next;
    free(g);
    return n;
}

Annotation *init_anno(){
    int init_n = 1<<8;
    Annotation *a;

    if ((a = malloc(sizeof(Annotation))) == NULL){
        err_msg(-1, 0, "init_anno: %s", strerror(errno));
        return NULL;
    } else {
        if ((a->chrms = malloc(init_n * sizeof(Chrm*))) == NULL){
            err_msg(-1, 0, "init_anno: %s", strerror(errno));
        }
        a->chrm_ix = init_str_map();
        a->gene_ix = init_str_map();
        a->gene_chrm = kh_init(str_int);
        a->gene_bin = kh_init(str_int);
        a->chrms_m = init_n;
    }

    return a;
}

Chrm *init_Chrm(){
    Chrm *chrm = (Chrm*)calloc(1, sizeof(Chrm));

    if (chrm == NULL){
        err_msg(-1, 0, "init_Chrm: %s", strerror(errno));
        return NULL;
    }

    int i;
    for (i = 0; i < MAX_BIN; i++){
        chrm->bins[i] = NULL;
        chrm->genes_n[i] = 0;
    }
    return chrm;
}

void destroy_anno(Annotation *a){

    if (a == NULL) return;

    int i;
    for (i = 0; i < a->chrm_ix->n; i++){
        Chrm *c = a->chrms[i];
        int j;
        // clear bins in chromosome
        for (j = 0; j < MAX_BIN; j++){
            Gene *g = c->bins[j];
            while (g) g = destroy_gene(g);
        }
        free(c);
    }
    free(a->chrms);

    destroy_str_map(a->chrm_ix);
    destroy_str_map(a->gene_ix);
    kh_destroy(str_int, a->gene_chrm);
    kh_destroy(str_int, a->gene_bin);

    free(a);
}

int add_chrm(Annotation *a, char *c){
    // add chromosome ID
    int chr_ix, found = 0;
    if ( (chr_ix = add2str_map(a->chrm_ix, (const char*)c, &found)) < 0 ) return -1;
    if (found == 0){
        while (a->chrm_ix->n >= a->chrms_m){
            a->chrms_m = (a->chrms_m)<<1;
            a->chrms = realloc(a->chrms, (a->chrms_m)*sizeof(Chrm*));

            if (a->chrms == NULL)
                return err_msg(-1, 0, "add_chrm: %s", strerror(errno));

        }
        a->chrms[chr_ix] = init_Chrm();
        if (a->chrms[chr_ix] == NULL) return -1;
    }
    return chr_ix;
}

Gene *gene_from_name_chrm(Annotation *a, int chrm_ix, char *gene_id){
    khint_t k;
    k = kh_get(str_int, a->gene_bin, gene_id);
    if (k == kh_end(a->gene_bin)){
        err_msg(-1, 0, "gene_from_name_chrm: gene %s not found in gene_bin", gene_id);
        return NULL;
    }
    int bin = kh_val(a->gene_bin, k);

    Gene *ag = a->chrms[chrm_ix]->bins[bin];
    for (; ag != NULL; ag = ag->next){
        if (strcmp(ag->id, gene_id) == 0) break;
    }
    return ag;
}

Gene *gene_from_name(Annotation *a, char *gene_id){
    khint_t k;
    k = kh_get(str_int, a->gene_bin, gene_id);
    if (k == kh_end(a->gene_bin)){
        err_msg(-1, 0, "gene_from_name: gene %s not found in gene_bin", gene_id);
        return NULL;
    }
    int bin = kh_val(a->gene_bin, k);

    k = kh_get(str_int, a->gene_chrm, gene_id);
    if (k == kh_end(a->gene_chrm)){
        err_msg(-1, 0, "gene_from_name: gene %s not found in gene_chrm", gene_id);
        return NULL;
    }
    int chrm_ix = kh_val(a->gene_chrm, k);

    Gene *ag = a->chrms[chrm_ix]->bins[bin];
    for (; ag != NULL; ag = ag->next){
        if (strcmp(ag->id, gene_id) == 0) break;
    }
    return ag;
}

/****************************
 * GTF processing
 ****************************/

gtf_line *init_gtf_line(){
    gtf_line *g = (gtf_line *)calloc(1, sizeof(gtf_line));

    if (g == NULL){
        err_msg(-1, 0, "init_gtf_line: %s", strerror(errno));
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

void clear_gtf_line(gtf_line *g){
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

void destroy_gtf_line(gtf_line *g){
    if (g == NULL) return;
    clear_gtf_line(g);
    free(g);
}

// Add gene from a GTF line
int gtf_gene2anno(Annotation *a, gtf_line *gl){
    int ci = add_chrm(a, gl->chrname.s);
    if (ci < 0) return -1;

    Gene *g = init_gene();
    if (g == NULL) return -1;
    
    char *gene_id = get_attr_val(gl, GENE_ID);
    if (gene_id == NULL){
        err_msg(-1, 0, "gtf_gene2anno: gtf line must have %s attribute", GENE_ID);
        return -1;
    }
    char *gene_name = get_attr_val(gl, GENE_NAME);
    char *gene_type = get_attr_val(gl, GENE_TYPE);
    if (gene_name == NULL || gene_type == NULL){
        return err_msg(-1, 0, "gtf_gene2anno: gtf line must have %s and %s "
                "attributes\n", GENE_NAME, GENE_TYPE);
    }

    int gix, found = 0;
    if ( (gix = add2str_map(a->gene_ix, (const char*)gene_id, &found)) < 0 ) return -1;

    // NOTE string object of g->id in g object is same as in gene_ix str_map. Free once
    g->id = str_map_str(a->gene_ix, gix);

    g->name = strdup(gene_name);
    if (g->name == NULL)
        return err_msg(-1, 0, "gtf_gene2anno: %s", strerror(errno));

    g->type = strdup(gene_type);
    if (g->type == NULL)
        return err_msg(-1, 0, "gtf_gene2anno: %s", strerror(errno));

    g->beg = gl->start;
    g->end = gl->end;
    g->strand = gl->strand;
    g->chrm = ci;
    g->bin = reg2bin(g->beg, g->end);

    // add to bins
    if (a->chrms[ci]->bins[g->bin]){
        Gene *ag = a->chrms[ci]->bins[g->bin];
        while (ag->next) ag = ag->next;
        ag->next = g;
    } else {
        a->chrms[ci]->bins[g->bin] = g;
    }
    (a->chrms[ci]->genes_n[g->bin])++;

    // add gene bin
    khint_t k;
    int ret;

    // printf("putting in gene_chrm\n");
    k = kh_put(str_int, a->gene_chrm, g->id, &ret);
    if (ret < 0)
        return err_msg(-1, 0, "gtf_gene2anno: failed to add %s to gene_chrm", gene_id);
    kh_val(a->gene_chrm, k) = ci;
    // printf("placed in gene_chrm\n");

    // add gene bin
    k = kh_put(str_int, a->gene_bin, g->id, &ret);
    if (ret < 0)
        return err_msg(-1, 0, "gtf_gene2anno: failed to add %s to gene_bin", gene_id);
    kh_val(a->gene_bin, k) = g->bin;

    return 0;
}

// add isoform from a GTF line
int gtf_iso2anno(Annotation *a, gtf_line *gl){
    int ci = add_chrm(a, gl->chrname.s);
    if (ci < 0) return -1;

    Isoform *iso = init_isoform();
    if (iso == NULL) return -1;
    
    // get GENE ID
    char *gene_id = get_attr_val(gl, GENE_ID);
    if (!gene_id){
        err_msg(-1, 0, "gtf_iso2anno: gtf line must have %s attribute", GENE_ID);
        return -1;
    }

    // get TX ID
    char *tx_id = get_attr_val(gl, TX_ID);
    if (!tx_id){
        err_msg(-1, 0, "gtf_iso2anno: gtf line for gene %s must have %s attribute", gene_id, TX_ID);
        return -1;
    }
    
    iso->beg = gl->start;
    iso->end = gl->end;

    // get gene object
    Gene *ag = gene_from_name_chrm(a, ci, gene_id);
    if (ag == NULL){
        err_msg(-1, 0, "gtf_iso2anno: gene %s not found", gene_id);
        return -1;
    }

    char *tx_id_cpy = calloc(strlen(tx_id)+1, sizeof(char));
    strcpy(tx_id_cpy, tx_id);
    
    int ret;
    khint_t k;

    // add key to hash
    k = kh_put(iso, ag->isoforms, tx_id_cpy, &ret);
    if (ret < 0){
        err_msg(-1, 0, "gtf_iso2anno: could not add %s to isoform hash table", tx_id);
        return -1;
    }

    // add isoform val to hash
    kh_val(ag->isoforms, k) = iso;
    ag->isoforms_n += 1;

    return 0;
}

// add exon from a GTF line
int gtf_exon2anno(Annotation *a, gtf_line *gl){
    int ci = add_chrm(a, gl->chrname.s);
    if (ci < 0) return -1;

    Exon *e = init_exon();
    if (e == NULL) return -1;
    
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
    
    e->beg = gl->start;
    e->end = gl->end;

    khint_t k;

    // get gene object
    Gene *ag = gene_from_name_chrm(a, ci, gene_id);
    if (ag == NULL){
        err_msg(-1, 0, "gtf_exon2anno: gene %s not found in bin index", gene_id);
        return -1;
    }

    // get isoform object
    k = kh_get(iso, ag->isoforms, tx_id);
    Isoform *iso = kh_val(ag->isoforms, k);

    // if no isoform
    if (k == kh_end(ag->isoforms) || !iso){
        err_msg(-1, 0, "gtf_exon2anno: gtf file is malformed. "
                "Found exon feature %s before transcript feature %s", exon_id, tx_id);
        return -1;
    }

    // sorted insert
    Exon *ge = iso->exons;
    while (ge->next != NULL && e->beg > ge->next->beg) ge = ge->next;

    if (ge->beg == e->beg){ // if exon found already
        err_msg(-1, 0, "gtf_exon2anno: duplicate exon positions found for %s %s %s", 
                gene_id, tx_id, exon_id);
        return -1;
    }

    e->next = ge->next;
    ge->next = e;
    iso->exons_n++;

    return 0;
}

int gtf_line2anno(Annotation *a, gtf_line *gl){

    int ret;
    char *type = gl->feature.s;

    if (strcmp(type, GENE) == 0)
        ret = gtf_gene2anno(a, gl);
    else if (strcmp(type, TX) == 0)
        ret = gtf_iso2anno(a, gl);
    else if (strcmp(type, EXON) == 0)
        ret = gtf_exon2anno(a, gl);
    else
        return 0;

    return ret;
}

int parse_gtf_line(kstring_t *line, gtf_line *g){
    char delims[] = "\t";
    char *token = NULL;
    char *rest = NULL;
    char *endptr = NULL;

    int save_errno = errno;
    errno = 0;

    if ( (token = strtok_r(line->s, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf_line: GTF line is malformed\n%s", line);
    if (kputs(token, &(g->chrname)) == EOF) return -1;
    
    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf_line: GTF line is malformed\n%s", line);
    if (kputs(token, &(g->source)) == EOF) return -1;

    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf_line: GTF line is malformed\n%s", line);
    if (kputs(token, &(g->feature)) == EOF) return -1;
    
    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf_line: GTF line is malformed\n%s", line);
    g->start = (int) strtol(token, &endptr, 10);
    if (g->start == 0 && errno > 0){
        return err_msg(-1, 0, "parse_gtf_line: could not convert %s to int: %s", 
                token, strerror(errno));
    }
    g->start -= 1; // convert start to 0-based coordinate
    
    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf_line: GTF line is malformed\n%s", line);
    g->end = strtol(token, &endptr, 10);
    if (g->end == 0 && errno > 0)
        return err_msg(-1, 0, "parse_gtf_line: could not convert %s to int: %s", 
                token, strerror(errno));
    
    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf_line: GTF line is malformed\n%s", line);
    if (token[0] == '.')
        g->score = -1;
    else{
        g->score = strtol(token, &endptr, 10);
        if (g->score == 0 && errno > 0)
            return err_msg(-1, 0, "parse_gtf_line: could not convert %s to int: %s", 
                    token, strerror(errno));
    }

    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf_line: GTF line is malformed\n%s", line);
    g->strand = token[0];
    // can be '+' '-' '.':irrelevant '?':relevant but unknown

    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf_line: GTF line is malformed\n%s", line);
    if (token[0] == '.')
        g->frame = -1;
    else{
        g->frame = strtol(token, &endptr, 10);
        if (g->frame == 0 && errno > 0)
            return err_msg(-1, 0, "parse_gtf_line: could not convert %s to int: %s",
                    token, strerror(errno));
    }

    if ( (token = strtok_r(NULL, delims, &rest)) == NULL)
        return err_msg(-1, 0, "parse_gtf_line: GTF line is malformed\n%s", line);;
    if (kputs(token, &(g->attribute)) == EOF) return -1;

    errno = save_errno;

    return 0;
}

int parse_gtf_attr(gtf_line *g){
    kstring_t *s = &(g->attribute);
    char delim[4] = " ;\"";

    char *token;
    char *rest;
    token = strtok_r(s->s, delim, &rest);
    while (token){
        kstr_node *tag_node = init_kstr_node();
        kstr_node *val_node = init_kstr_node();
        if (tag_node == NULL || val_node == NULL) return -1;

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
char *get_attr_val(gtf_line *g, char *key){
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

int has_key_val(gtf_line *g, char *key, char *val){
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

Annotation *read_from_gtf(char *file, int basic){

    Annotation *a = init_anno();

    if (a == NULL) return NULL;

    int ret;

    kstring_t line = KS_INITIALIZE;
    int len = 0;

    const char mode_read[3] = "r1\0";
    BGZF *fp = (BGZF *)bgzf_open(file, mode_read);
    if (fp == 0){
        err_msg(-1, 0, "read_from_gtf: could not open GTF file %s", file);
        destroy_anno(a);
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
    gtf_line *gl = init_gtf_line();
    while ( (len = bgzf_getline(fp, '\n', &line)) > 0 ){
        nlines++;

        if ( (ret = parse_gtf_line(&line, gl)) < 0){
            err_msg(-1, 0, "read_from_gtf: failed to parse GTF line\n%s", line.s);
            destroy_gtf_line(gl);
            destroy_anno(a);
            bgzf_close(fp);
            return NULL;
        }
        parse_gtf_attr(gl);

        // filter basic tag from tx or exon gtf line
        int line_is_basic = has_key_val(gl, "tag", "basic");
        if ((strcmp(gl->feature.s, GENE) != 0) && // not gene
            (basic)                            && // filter for basic tx
            (!line_is_basic)){ // line is not basic
            clear_gtf_line(gl);
            continue;
        }

        // add gtf line
        if ( (ret = gtf_line2anno(a, gl)) == -1 ){
            err_msg(-1, 0, "read_from_gtf: failed to add GTF line\n%s", line.s);
            destroy_gtf_line(gl);
            destroy_anno(a);
            bgzf_close(fp);
            return NULL;
        }

        clear_gtf_line(gl);

        // printf("Feature: name = %s; id = %s; beg = %i; end = %i; strand = %c; chrm = %i; bin = %i\n", f.name.s, f.id.s, f.beg, f.end, f.strand, f.chrm, f.bin);
    } // end read GTF

    ks_free(&line);
    destroy_gtf_line(gl);
    bgzf_close(fp);

    return a;
}

int write_gene_data(BGZF *fp, Annotation *anno, str_map *gene_ix){
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

        Gene *gene_obj = gene_from_name(anno, gene_key);

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


int feats_from_region_p(const Annotation *a, const char* ref, 
        int64_t beg, int64_t end, uint8_t stranded, char strand, 
        Gene ***genes, int *genes_len, double p){

    // overlap chromosomes
    int tid = str_map_ix(a->chrm_ix, (char *)ref);
    if (tid < 0) 
        return 0;

    // region length
    double reg_len = (double)end - (double)beg;
    if (reg_len < 0) 
        return err_msg(-1, 0, "feats_from_region_p: end (%lli) < beg (%lli)", (int64_t)end, (int64_t)beg);

    int ngenes = 0;

    // reallocate genes array
    if (*genes_len <= 0){
        *genes_len = 1;
        *genes = (Gene **)realloc(*genes, *genes_len * sizeof(Gene *));
        if (*genes == NULL)
            return err_msg(-1, 0, "feats_from_region_p: %s", strerror(errno));
    }

    // get bins for region
    uint16_t list[MAX_BIN];
    int n_bin = reg2bins((int)beg, (int)end, list);

    int i;
    for (i = 0; i < n_bin; ++i){ // for each bin
        uint16_t bin = list[i];
        Gene *g = a->chrms[tid]->bins[bin]; 
        for (; g; g = g->next){ // for each gene in bin
            char reg_strand = '.';
            if (stranded) reg_strand = strand;

            int ovrlp = bp_overlap(beg, end, reg_strand, g->beg, g->end, g->strand);

            if (ovrlp < 0) {
                err_msg(-1, 0, "feats_from_region_p: failed to get bp_overlap");
                return -1;
            }

            // fraction overlap of region
            double frac_ovrlp = (double)ovrlp / reg_len;
            if (frac_ovrlp < p) continue;

            // reallocate genes array
            while (ngenes >= *genes_len){
                *genes_len = (*genes_len) * 2;
                *genes = (Gene **)realloc(*genes, (*genes_len) * sizeof(Gene *));
                if (*genes == NULL)
                    return err_msg(-1, 0, "feats_from_region_p: %s", strerror(errno));
            }

            if (g->id == NULL){ // if gene has no ID
                err_msg(-1, 0, "feats_from_region_p: overlapping genes must have IDs. Beg=%i End=%i", 
                        g->beg, g->end);
                return -1;
            }

            // add to @p genes
            (*genes)[ngenes++] = g;
        }
    }
    return ngenes;
}


/****************************
 * summary functions
 ****************************/

int n_feat(Annotation *a, int *n_gene, int *n_iso, int *n_exon){
    *n_gene = 0;
    *n_iso = 0;
    *n_exon = 0;
    int nc = a->chrm_ix->n;
    int i, j;
    for (i = 0; i < nc; i++){
        Chrm *chrm = a->chrms[i];
        for (j = 0; j < MAX_BIN; j++){
            Gene *g = chrm->bins[j];
            while (g){
                *n_gene = *n_gene + 1;
                khint_t k_iso;
                fprintf(stdout, "gene %s has %i isoforms\n", g->id, g->isoforms_n);
                for (k_iso = kh_begin(g->isoforms); k_iso != kh_end(g->isoforms); ++k_iso){
                    if (!kh_exist(g->isoforms, k_iso)) continue;
                    Isoform *iso = kh_val(g->isoforms, k_iso);
                    char *key = kh_key(g->isoforms, k_iso);
                    fprintf(stdout, "ID %s beg %i end %i n_exon %i\n", key, iso->beg, iso->end, iso->exons_n);
                    *n_iso = *n_iso + 1;
                    *n_exon = *n_exon + iso->exons_n;
                }
                g = g->next;
            }
        }
    }
    return nc;
}

int test_anno(char *file){
    fprintf(stdout, "Reading from GTF\n");
    Annotation *a = read_from_gtf(file, 1);
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
    destroy_anno(a);
    return 0;
}


