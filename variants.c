
#include "variants.h"
#include "bins.h"
#include "gtf_anno.h"
#include "overlap.h"
#include "str_util.h"
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include <math.h>
#include "g_list.h"

/************************************************************************
 * VCF file
 ***********************************************************************/

int load_vcf(const char *vcffn, const char *region, int region_set, 
        bcf_srs_t **sr, bcf_hdr_t **vcf_hdr){
    *sr = bcf_sr_init();
    bcf_sr_set_opt(*sr, BCF_SR_REQUIRE_IDX);
    if ( region_set ){
        int tret = bcf_sr_set_regions(*sr, region, 0);
        if ( tret == -1 )
            return err_msg(-1, 0, "load_vcf: could not set region %s in VCF file", region);
    }
    if ( bcf_sr_add_reader(*sr, vcffn) == 0 )
        return err_msg(-1, 0, "load_vcf: could not read VCF file %s", vcffn);
    *vcf_hdr = (*sr)->readers[0].header;

    return(0);
}

int sub_vcf_samples(bcf_hdr_t **vcf_hdr, const char *samplefn){
    if (samplefn == NULL) return(0);

    int rets = bcf_hdr_set_samples(*vcf_hdr, (const char*)samplefn, 1);
    if (rets < 0)
        return err_msg(-1, 0, "sub_vcf_samples: failed to subset samples in VCF");
    if (rets > 0)
        return err_msg(-1, 0, "sub_vcf_samples: failed to subset samples in VCF: sample %i not present", rets);
    return(0);
}

int bcf_hdr_chr_ix(const bcf_hdr_t *hdr, str_map *cm){
    if (hdr == NULL || cm == NULL)
        return err_msg(-1, 0, "bcf_hdr_to_sm: hdr or sm is NULL");

    int ret, i, n_chr = hdr->n[BCF_DT_CTG], found = 0;
    for (i = 0; i < n_chr; ++i){
        const char *hdr_chr = bcf_hdr_id2name(hdr, i);

        // check that chr ix matches between hdr and i
        int tmp_id = bcf_hdr_name2id(hdr, hdr_chr);
        if (i != tmp_id)
            return err_msg(-1, 0, "bcf_hdr_to_sm: index %i does not match "
                    "vcf header rid %i", i, tmp_id);

        ret = add2str_map(cm, hdr_chr, &found);
        if (ret < 0) return(-1);
        if (found)
            return err_msg(-1, 0, "bcf_hdr_to_sm: chromosome %s found more "
                    "than once", hdr_chr);
    }
    return(0);
}

/************************************************************************
 * 
 ***********************************************************************/

int get_var_len(bcf1_t *b){
    int len = 1;
    int i;
    for (i = 0; i < b->n_allele; i++){
        int ilen = strlen(b->d.allele[i]);
        if (ilen > len) 
            len = ilen;
    }
    return len;
}

int get_var_bin(bcf1_t *b){
    int beg = b->pos;
    int len = get_var_len(b);
    int end = beg + len;
    int bin = reg2bin(beg, end);
    return bin;
}

char *var_id(const bcf_hdr_t *h, bcf1_t *b, char delim){

    // chromosome
    const char *chrm = bcf_seqname(h, b);
    size_t nchrm = strlen(chrm);

    // pos
    size_t npos = 0;
    char *pos = int2str(b->pos + 1, &npos);
    if (pos == NULL) return NULL;

    // ID
    size_t nid = strlen(b->d.id);

    // alleles
    size_t allele_len = 0;
    int *allele_lens = (int *)malloc(sizeof(int) * b->n_allele);
    int i;
    for (i = 0; i < b->n_allele; ++i){
        allele_lens[i] = strlen(b->d.allele[i]);
        if (i > 1) allele_lens[i] += 1;
        allele_len += allele_lens[i];
    }

    // total
    size_t ntotal = nchrm + npos + nid + allele_len + 5;
    char *out = (char *)calloc(ntotal, sizeof(char));
    if (out == NULL){
        err_msg(-1, 0, "var_id: %s", strerror(errno));
        return(NULL);
    }
    char *tmp = out;

    memcpy(tmp, chrm, nchrm * sizeof(char));
    tmp = tmp + nchrm;
    *(tmp++) = delim;

    memcpy(tmp, pos, npos * sizeof(char));
    tmp = tmp + npos;
    *(tmp++) = delim;

    memcpy(tmp, b->d.id, nid * sizeof(char));
    tmp = tmp + nid;
    *(tmp++) = delim;

    memcpy(tmp, b->d.allele[0], allele_lens[0] * sizeof(char));
    tmp = tmp + allele_lens[0];
    *(tmp++) = delim;

    for (i = 1; i < b->n_allele; ++i){
        if (i > 1) *(tmp++) = ';';
        memcpy(tmp, b->d.allele[i], allele_lens[i] * sizeof(char));
        tmp = tmp + allele_lens[i];
    }

    free(pos);
    free(allele_lens);

    return out;
}

uint8_t base_ref_alt(bcf1_t *b, char base){
    if ( base == b->d.allele[0][0] ) return ((uint8_t)REF);
    else if ( base == b->d.allele[1][0] ) return ((uint8_t)ALT);
    else return((uint8_t)OTHER);
}

/************************************************************************
 * var_t
 ***********************************************************************/

var_t *var_alloc(){
    var_t *v = (var_t*)calloc(1, sizeof(var_t));
    if (v == NULL){
        err_msg(-1, 0, "var_alloc: %s", strerror(errno));
        return NULL;
    }
    v->vix = -1;
    v->next = NULL;
    return v;
}

g_var_t *init_genomevar(){
    int init_n = 1<<8;
    g_var_t *gv = (g_var_t*)calloc(1, sizeof(g_var_t));
    
    if (gv == NULL){
        err_msg(-1, 0, "init_genomevar: %s", strerror(errno));
        return NULL;
    }

    mv_init(&gv->chr_bins);

    gv->chrm_ix = init_str_map();
    gv->chrms = (chr_var_t**)calloc(init_n, sizeof(chr_var_t*));

    if (gv->chrms == NULL){
        err_msg(-1, 0, "init_genomevar: %s", strerror(errno));
        return NULL;
    }
    gv->chrms_m = init_n;
    
    mv_init(&gv->vix2var);

    gv->n_v = 0;
    gv->n_e = 0;
    gv->n_a = 1<<8;
    gv->ix2var = (var_t **)calloc(gv->n_a, sizeof(var_t *));
    if (gv->ix2var == NULL){
        err_msg(-1, 0, "init_genomevar: %s", strerror(errno));
        return NULL;
    }

    return gv;
}

chr_bins_t *chr_bins_alloc(){
    chr_bins_t *bins = calloc(1, sizeof(chr_bins_t));
    if (bins == NULL){
        err_msg(-1, 0, "chr_bins_init: %s", strerror(errno));
        return NULL;
    }

    int i;
    for (i = 0; i < MAX_BIN; ++i){
        if (vcfr_list_init(&bins->bins[i]) < 0){
            err_msg(-1, 0, "chr_bins_init: failed to initialize list");
            return(NULL);
        }
    }
    return(bins);
}

void chr_bins_free(chr_bins_t *bins){
    if (bins == NULL) return;

    int i;

    // free the bcf records, then the list
    for (i = 0; i < MAX_BIN; ++i){
        ml_node_t(vcfr_list) *rn;
        for (rn = ml_begin(&bins->bins[i]); rn; rn = ml_node_next(rn)){
            var_t var = ml_node_val(rn);
            if (var.b) bcf_destroy(var.b);
        }
        vcfr_list_free(&bins->bins[i]);
    }
}

chr_var_t *init_chr_var(){
    chr_var_t *chrm = (chr_var_t*)calloc(1, sizeof(chr_var_t));

    if (chrm == NULL){
        err_msg(-1, 0, "init_chr_var: %s", strerror(errno));
        return NULL;
    }

    int i;
    for (i = 0; i < MAX_BIN; i++){
        chrm->bins[i] = NULL;
        chrm->vars_n[i] = 0;
    }
    return chrm;
}

int gv_add_hdr(g_var_t *gv, bcf_hdr_t *hdr){
    if (gv == NULL || hdr == NULL)
        return err_msg(-1, 0, "gv_add_hdr: arguments are NULL");

    gv->vcf_hdr = hdr;
    return(0);
}

int add_var(g_var_t *gv, bcf1_t *b, const bcf_hdr_t *hdr){
    var_t *tv = var_alloc();
    if (tv == NULL) return -1;
    tv->b = bcf_dup(b);
    tv->vix = gv->n_e;
    bcf_unpack(tv->b, BCF_UN_FLT); // unpack up to and including info field

    int bin = get_var_bin(b);
    int found;
    int rid = tv->b->rid;
    const char *seqname = bcf_hdr_id2name(hdr, rid);

    // add chromosome ID
    int chr_ix, bs = (int)mv_size(&gv->chr_bins);
    if ( (chr_ix = add2str_map(gv->chrm_ix, (const char*)seqname, &found)) < 0 ) return -1;

    if (chr_ix >= bs || mv_i(&gv->chr_bins, chr_ix) == NULL){
        chr_bins_t *p = chr_bins_alloc();
        if (p == NULL)
            return(-1);

        int ret = mv_insert(cbin_vec, &gv->chr_bins, p, chr_ix);
        if (ret < 0)
            return err_msg(-1, 0, "add_var: failed to add chr bin");
    }

    printf("size=%zu chr_ix=%i\n", mv_size(&gv->chr_bins), chr_ix);

    if (found == 0){
        while (gv->chrm_ix->n >= gv->chrms_m){
            gv->chrms_m = (gv->chrms_m)<<1;
            gv->chrms = realloc(gv->chrms, (gv->chrms_m)*sizeof(chr_var_t*));
            if (gv->chrms == NULL)
                return err_msg(-1, 0, "add_var: %s", strerror(errno));
        }
        gv->chrms[chr_ix] = init_chr_var();
        if (gv->chrms[chr_ix] == NULL) return -1;
    }

    // add to bins
    ml_vcfr_list_t *ll = &mv_i(&gv->chr_bins, chr_ix)->bins[bin];
    if ( vcfr_list_insert(ll, *tv, 0, 0) < 0 )
        return err_msg(-1, 0, "add_var: failed to insert to list");

    if (gv->chrms[chr_ix]->bins[bin]){
        var_t *gvv = gv->chrms[chr_ix]->bins[bin];
        while (gvv->next){
            gvv = gvv->next;
        }
        gvv->next = tv;
    }
    else {
        gv->chrms[chr_ix]->bins[bin] = tv;
    }
    (gv->chrms[chr_ix]->vars_n[bin])++;

    // add to ix2var
    if (mv_push(vcfr_vec, &gv->vix2var, *tv) < 0)
        return err_msg(-1, 0, "add_var: failed to add var to vix2var");

    while (gv->n_e >= gv->n_a){
        gv->n_a <<= 1;
        gv->ix2var = realloc(gv->ix2var, (gv->n_a) * sizeof(var_t));
        if (gv->ix2var == NULL)
            return err_msg(-1, 0, "add_var: %s", strerror(errno));
    }
    gv->ix2var[gv->n_e] = tv;
    ++gv->n_e;
    ++gv->n_v;

    return 0;
}

g_var_t *vcf2gv(bcf_srs_t *sr, bcf_hdr_t *vcf_hdr, int max_miss, double maf_cut){
    g_var_t *gv = init_genomevar();
    if (gv == NULL) return NULL;

    gv->vcf_hdr = vcf_hdr;

    size_t i;

    // add chromosomes from header
    if (bcf_hdr_chr_ix(vcf_hdr, gv->chrm_ix) < 0) return(NULL);
    size_t n_chr = (size_t)gv->chrm_ix->n;

    // create chromosome bin vector
    if (mv_resize(cbin_vec, &gv->chr_bins, n_chr) < 0) return(NULL);

    size_t cix;
    for (cix = 0; cix < n_chr; ++cix){
        chr_bins_t *p = chr_bins_alloc();
        if (p == NULL)
            return(NULL);

        int ret = mv_insert(cbin_vec, &gv->chr_bins, p, cix);
        if (ret < 0){
            err_msg(-1, 0, "vcf2gv: failed to add chr bin");
            return(NULL);
        }
    }

    // initialize chromosomes read from header
    if (n_chr > 0)
        gv->chrms = (chr_var_t **)realloc(gv->chrms, n_chr * sizeof(chr_var_t *));
    for (i = 0; i < n_chr; ++i){
        gv->chrms[i] = init_chr_var();
        if (gv->chrms[i] == NULL){
            err_msg(-1, 0, "vcf2gv: failed to init chr_var_t");
            return(NULL);
        }
    }

    while ( bcf_sr_next_line(sr) ){
        bcf1_t *vcf_r = bcf_sr_get_line(sr, 0);
        if (vcf_r == NULL) fprintf(stderr, "warning: a VCF record is NULL\n");

        if (is_biallelic_snp(vcf_r) < 0) continue;

        int na = 0, nmiss=0, ntotal=1;
        if (bcf1_t_miss_maf(vcf_hdr, vcf_r, &nmiss, &na, &ntotal) < 0)
            return(NULL);

        if (max_miss >= 0 && nmiss > max_miss) continue;
        
        double maf = ((double)na) / ((double)ntotal);
        if (maf_cut >= 0 && (maf <= maf_cut || maf >= (1 - maf_cut))) continue;

        if ( add_var(gv, vcf_r, vcf_hdr) < 0 ) return NULL;
    }

    return gv;
}

void destroy_gv(g_var_t *gv){

    if (gv == NULL) return;

    // chrms
    int i, j;
    for (i = 0; i < gv->chrm_ix->n; i++){
        for (j = 0; j < MAX_BIN; j++){
            var_t *v = gv->chrms[i]->bins[j];
            while (v){
                bcf1_t *b = v->b;
                if (b) bcf_destroy(b);
                var_t *vn = v->next;
                free(v);
                v = vn;
            }
        }
        free(gv->chrms[i]);
    }
    free(gv->chrms);
    destroy_str_map(gv->chrm_ix);

    mv_free(&gv->chr_bins);
    mv_free(&gv->vix2var);

    free(gv->ix2var);

    free(gv);
}

void free_bcf1_t(g_var_t *gv){
    int i, j;
    for (i = 0; i < gv->chrm_ix->n; i++){
        for (j = 0; j < MAX_BIN; j++){
            var_t *v = gv->chrms[i]->bins[j];
            while (v){
                bcf1_t *b = v->b;
                if (b) bcf_destroy(b);
                v->b = NULL;
                v = v->next;
            }
        }
    }
}

int is_biallelic_snp(bcf1_t *b){
    if (b->n_allele != 2 || 
        bcf_get_variant_type(b, 1) != VCF_SNP) 
        return -1;
    return 0;
}

int bcf1_t_miss_maf(bcf_hdr_t *vcf_hdr, bcf1_t *b, int *nmiss, int *na, int *ntotal){
    bcf_unpack(b, BCF_UN_FMT);

    // check if fmt tag is present
    bcf_fmt_t *fmt = bcf_get_fmt(vcf_hdr, b, "GT");
    if (!fmt)
        return err_msg(-1, 0, "bcf1_t_miss_maf: format not present for %s", b->d.id);

    int n_samples = bcf_hdr_nsamples(vcf_hdr);
    int s;
    void *gt_arr = NULL;
    int ndst = 0, ngt;
    ngt = bcf_get_format_values(vcf_hdr, b, "GT", &gt_arr, &ndst, BCF_HT_INT); // GT must be given as BCF_HT_INT

    // check if tag is present again
    if (ngt < 0){
        free(gt_arr);
        return err_msg(-1, 0, "bcf1_t_miss_maf: format not present for %s", b->d.id);
    }
    int n_val = ngt / n_samples;

    if (n_val != fmt->n){
        fprintf(stderr, "bcf1_t_miss_maf: nval doesn't match: %i vs %i\n", n_val, fmt->n);
    }

    // check missing
    void *p = gt_arr;
    int n_miss = 0, n_a = 0, n_total = 0;
    for (s = 0; s < n_samples; s++){ // each sample
        int t;
        for (t = 0; t < n_val; t++){ // each allele
            int32_t *a = (int32_t *)p;
            a = a + (s * n_val) + t;

            // if end, sample has smaller ploidy, break
            if ( *a == bcf_int32_vector_end ) break;

            if ( bcf_gt_is_missing(*a) ){
                n_miss++;
            } else {
                if (bcf_gt_allele(*a) > 0) n_a++;
            }
            n_total++;
        }
    }
    *nmiss = n_miss;
    *na = n_a;
    *ntotal = n_total;
    free(gt_arr);

    return 0;
}

int is_gt_valid(bcf_hdr_t *vcf_hdr, bcf1_t *b){
    bcf_unpack(b, BCF_UN_FMT);

    // variant must be bi-allelic
    if (is_biallelic_snp(b) < 0) return -2;
    if (b->n_allele != 2) return 0;

    // check if fmt tag is present
    bcf_fmt_t *fmt = bcf_get_fmt(vcf_hdr, b, "GT");
    if (!fmt) return -1;

    int n_miss = 0, n_a = 0, n_total = 0;
    if (bcf1_t_miss_maf(vcf_hdr, b, &n_miss, &n_a, &n_total) < 0){
        return -3;
    }

    // if (n_miss == n_samples) return -1;

    return 0;
}
        
int is_gp_valid(bcf_hdr_t *vcf_hdr, bcf1_t *b){
    bcf_unpack(b, BCF_UN_FMT);

    // variant must be bi-allelic
    if (is_biallelic_snp(b) < 0) return -2;
    if (b->n_allele != 2) return 0;

    // check if fmt tag is present
    bcf_fmt_t *fmt = bcf_get_fmt(vcf_hdr, b, "GP");
    if (!fmt) return -1;

    return 0;
}
        
float *bcf1_ap_gt(bcf_hdr_t *vcf_hdr, bcf1_t *b, int extra){

    bcf_unpack(b, BCF_UN_FMT);

    int n_samples = bcf_hdr_nsamples(vcf_hdr);
    if (n_samples == 0) return NULL;

    int s;
    void *gt_arr = NULL;
    int ndst = 0, ngt;
    ngt = bcf_get_format_values(vcf_hdr, b, "GT", &gt_arr, &ndst, BCF_HT_INT); // GT must be given as BCF_HT_INT
    if (ngt < 0){
        err_msg(-1, 0, "bcf1_ap_gt: failed to get genotypes");
        free(gt_arr);
        return NULL;
    }
    int n_val = ngt / n_samples;

    // allocate dose array
    float *dose = malloc((n_samples + extra) * sizeof(float));
    if (dose == NULL){
        err_msg(-1, 0, "bcf1_ap_gt: %s", strerror(errno));
        return NULL;
    }

    // fill dose array
    void *p = gt_arr;
    int n_miss = 0;
    for (s = 0; s < n_samples; s++){ // each sample
        float gt = 0;
        int t, total = 0;
        for (t = 0; t < n_val; t++){ // each allele
            int32_t *a = (int32_t *)p;
            a = a + (s * n_val) + t;

            // if end, sample has smaller ploidy, break
            if ( *a == bcf_int32_vector_end ) break;

            // if any missing
            if ( bcf_gt_is_missing(*a) ) {
                n_miss++;
                gt = -1.0;
                total = 1;
                break;
            }

            gt += (float)bcf_gt_allele(*a);
            total++;
        }
        
        if (total == 0) dose[s] = -1;
        else dose[s] = gt / total;
    }
    free(gt_arr);

    // if all missing, return NULL
    // if (n_miss == n_samples) return NULL;

    return dose;
}

float *bcf1_ap_gp(bcf_hdr_t *vcf_hdr, bcf1_t *b, int extra){

    bcf_unpack(b, BCF_UN_FMT);

    int n_samples = bcf_hdr_nsamples(vcf_hdr);
    if (n_samples == 0) return NULL;

    int s;
    void *gt_arr = NULL;
    int ndst = 0, ngt;
    ngt = bcf_get_format_values(vcf_hdr, b, "GP", &gt_arr, &ndst, BCF_HT_REAL);
    if (ngt < 0){
        err_msg(-1, 0, "bcf1_dose_gp: failed to get genotypes %s", b->d.id);
        free(gt_arr);
        return NULL;
    }
    int n_val = ngt / n_samples;

    bcf_fmt_t *fmt = bcf_get_fmt(vcf_hdr, b, "GP");
    if (fmt->n != n_val)
        fprintf(stderr, "fmt->n %i doesn't match n_val %i for SNP %s\n", fmt->n, n_val, b->d.id);

    // allocate dose array
    float *dose = malloc((n_samples + extra) * sizeof(float));
    if (dose == NULL){
        err_msg(-1, 0, "bcf1_dose_gp: %s", strerror(errno));
        free(gt_arr);
        return NULL;
    }

    // fill dose array
    void *p = gt_arr;
    int n_miss = 0;
    for (s = 0; s < n_samples; s++){ // each sample
        float gt = 0;
        int t, total = 0;
        for (t = 0; t < n_val; t++){ // each allele
            float *a = (float *)p;
            a = a + (s * n_val) + t;

            // if end, sample has smaller ploidy, break
             if ( (uint32_t)*a == bcf_float_vector_end ) break;

            // if any are missing, set to -1
            if ( (uint32_t)*a == bcf_float_missing ){
                fprintf(stdout, "missing\n");
                n_miss++;
                gt = -1.0;
                total = 1;
                break;
            }

            gt += t * *a;
            total++;
        }

        switch ( total ) {
            case 0: 
                dose[s] = -1;
                break;
            case 1: // when prob. is given directly
                dose[s] = gt;
                break;
            case 2:
                dose[s] = -1;
                break;
            case 3:
                dose[s] = gt / 2;
                break;
            default:
                err_msg(-1, 0, "bcf1_dose_gp: total number of genotypes for %s is inconsistent: %i", 
                        b->d.id, total);
                return NULL;
        }
    }
    free(gt_arr);

    // if all missing, return NULL
    if (n_miss == n_samples) return NULL;

    return dose;
}

float **ap_array_gt(g_var_t *gv, bcf_hdr_t *vcf_hdr, int32_t *ids, int ni, char *field){
    int j;

    // check input
    if ( strcmp(field, "GP") != 0 && strcmp(field, "GT") != 0){
        err_msg(-1, 0, "ap_array_gt: %s: genotype field must be one of GT or GP\n", field);
        return NULL;
    }

    // check if field is present in header
    int tag_id = bcf_hdr_id2int(vcf_hdr, BCF_DT_ID, field);

    if (tag_id < 0){
        err_msg(-1, 0, "ap_array_gt: %s is not present in the VCF header\n", field);
        return NULL;
    }

    // allocate ap array
    float **ap = malloc(ni * sizeof(float *));
    if (ap == NULL){
        err_msg(-1, 0, "ap_array_gt: %s", strerror(errno));
        return NULL;
    }
    int i;
    for (i = 0; i < ni; ++i) ap[i] = NULL;

    // loop over each variant
    for (i = 0; i < ni; ++i){
        int32_t tix = ids[i];
        var_t *var = gv_vari(gv, tix);
        // If var is NULL, let ap[i] = NULL
        // TODO: might need a better way to handle missing variants
        if (var == NULL){
            continue;
            // err_msg(-1, 0, "ap_array_gt: index %i is NULL\n", tix);
            // return(NULL);
        }
        char *vid = var_id(vcf_hdr, var->b, '\t');

        // check if genotype fields are valid
        int is_valid;
        if (strcmp(field, "GP") == 0) is_valid = is_gp_valid(vcf_hdr, var->b);
        else is_valid = is_gt_valid(vcf_hdr, var->b);
        if (is_valid <= -2)
            err_msg(-1, 0, "ap_array_gt: %s is not bi-allelic\n", vid);
        if (is_valid == -1)
            err_msg(-1, 0, "ap_array_gt: %s does not have format data for %s\n", vid, field);
        if (is_valid < 0) goto cleanup;

        // get alt allele prob.
        if ( strcmp(field, "GP") == 0 )
            ap[i] = bcf1_ap_gp(vcf_hdr, var->b, 0);
        else
            ap[i] = bcf1_ap_gt(vcf_hdr, var->b, 0);

        // if error getting genotypes
        if ( ap[i] == NULL ){
            err_msg(-1, 0, "ap_array_gt: error getting genotypes from variant %s", vid);
            goto cleanup;
        }
        free(vid);
    }
    return ap;
cleanup:
    for (j = 0; j < i; ++j) free(ap[j]);
    free(ap);
    return NULL;
}

int region_vars(g_var_t *gv, const char* ref, int32_t beg, 
        int32_t end, var_t **vars){
    if (gv == NULL)
        return err_msg(-2, 0, "region_vars: gv must not be NULL");

    int tid = str_map_ix(gv->chrm_ix, (char *)ref);
    if (tid < 0)
        return err_msg(-1, 1, "region_vars: Chromosome %s not found in VCF", ref);

    double reg_len = (double)end - (double)beg;
    if (reg_len < 0)
        return err_msg(-2, 0, "region_vars: end (%i) < beg (%i)", beg, end);

    int n_vars = 0;

    var_t *pv = NULL, *begin = NULL;
    uint16_t list[MAX_BIN];
    int n_bin = reg2bins((int)beg, (int)end, list);
    int i;
    for (i = 0; i < n_bin; ++i){
        uint16_t bin = list[i];
        ml_vcfr_list_t *ll = &mv_i(&gv->chr_bins, tid)->bins[bin];
        ml_node_t(vcfr_list) *vn;
        for (vn = ml_begin(ll); vn; vn = ml_node_next(vn)){
            var_t var = ml_node_val(vn);
            bcf1_t *b = var.b;
            int ovrlp = bp_overlap((int)beg, (int)end, '.', b->pos, b->pos + b->rlen, '.');
            if (ovrlp < 0)
                return err_msg(-2, 0, "region_vars: failed to get bp_overlap");
            if (ovrlp == 0) continue;

            // Add var_t to list
            ++n_vars;
            var_t *cpy = (var_t *)calloc(1, sizeof(var_t));
            if (cpy == NULL)
                return err_msg(-2, 0, "region_vars: %s", strerror(errno));

            *cpy = var;
            cpy->next = NULL;
            if ( pv == NULL ){
                pv = cpy;
                begin = pv;
            }
            else {
                pv->next = cpy;
                pv = pv->next;
            }
        }
    }

    // If var was NULL
    if (*vars == NULL) *vars = begin;
    else{
        pv = *vars;
        while (pv->next != NULL) pv = pv->next;
        pv->next = begin;
    }

    return(n_vars);
}

int n_snp(g_var_t *gv, int *n_snp){
    *n_snp = 0;
    int nc = gv->chrm_ix->n;
    int i, j;
    for (i = 0; i < nc; i++){
        for (j = 0; j < MAX_BIN; j++){
            *n_snp = *n_snp + gv->chrms[i]->vars_n[j];
        }
    }
    return nc;
}

var_t *gv_vari(g_var_t *gv, int32_t ix){
    if (gv == NULL){
        err_msg(-1, 0, "gv_vari: gv is NULL");
        return(NULL);
    }
    if (ix < 0 || ix >= gv->n_e){
        err_msg(-1, 0, "gv_vari: ix=%i is invalid, only [%i,%i) are valid", 
                ix, 0, gv->n_e);
        return(NULL);
    }
    var_t *var = gv->ix2var[ix];
    return(var);
}

