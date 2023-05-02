
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
 * bcf functions
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

int get_bcf1_bin(bcf1_t *b){
    int beg = b->pos;
    int len = get_var_len(b);
    int end = beg + len;
    int bin = reg2bin(beg, end);
    return bin;
}

int is_biallelic_snp(bcf1_t *b){
    if (b->n_allele != 2 || 
        bcf_get_variant_type(b, 1) != VCF_SNP) 
        return -1;
    return 0;
}

int bcf1_t_miss_maf(bcf_hdr_t *vcf_hdr, bcf1_t *b, int *nmiss, int *n_allele, 
        int *ntotal){
    bcf_unpack(b, BCF_UN_FMT);

    // only bi-allelic variants are supported
    if (b->n_allele > 2)
        return err_msg(-1, 0, "bcf1_t_miss_maf: only biallelic variants, "
                "check %s", b->d.id);

    // check if fmt tag "GT" is present before access
    bcf_fmt_t *fmt = bcf_get_fmt(vcf_hdr, b, "GT");
    if (!fmt)
        return err_msg(-1, 0, "bcf1_t_miss_maf: format not present for %s", b->d.id);

    void *gt_arr = NULL;
    int ndst = 0, ngt, n_samples = bcf_hdr_nsamples(vcf_hdr);
    ngt = bcf_get_format_values(vcf_hdr, b, "GT", &gt_arr, &ndst, BCF_HT_INT); // GT must be given as BCF_HT_INT

    // check for error, if GT was not present
    if (ngt < 0){
        free(gt_arr);
        return err_msg(-1, 0, "bcf1_t_miss_maf: format not present for %s", b->d.id);
    }

    int n_val = ngt / n_samples;
    assert(n_val == fmt->n);

    // check missing
    int32_t *p = (int32_t *)gt_arr;
    int s, t, n_miss = 0, n_a = 0, n_total = 0;
    for (s = 0; s < n_samples; ++s){ // each sample
        int32_t *a = p + s*n_val;
        for (t = 0; t < n_val; ++t){ // each allele
            // if end, sample has smaller ploidy, continue
            if ( a[t] == bcf_int32_vector_end ) continue;

            int32_t is_miss = bcf_gt_is_missing(a[t]);
            int32_t sgt = bcf_gt_allele(a[t]);
            if (sgt >= 2)
                return err_msg(-1, 0, "bcf1_t_miss_maf: there are multiple "
                        "alleles for variant '%s' and sample '%s'. Only bi-allelic variants "
                        "are supported", b->d.id, vcf_hdr->samples[s]);
            // if this fails, need new way to check for bi-allelic
            assert(sgt < 2);

            if (is_miss)
                ++n_miss;
            else if (sgt > 0)
                ++n_a;
            
            ++n_total;
        }
    }
    *nmiss = n_miss;
    *n_allele = n_a;
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
        
/************************************************************************
 * var_t functions
 ***********************************************************************/

var_t *var_alloc(){
    var_t *v = (var_t*)calloc(1, sizeof(var_t));
    if (v == NULL){
        err_msg(-1, 0, "var_alloc: %s", strerror(errno));
        return NULL;
    }
    v->vix = -1;
    return v;
}

void var_free(var_t *var){
    if (var == NULL) return;
    if (var->b) bcf_destroy(var->b);
}

/************************************************************************
 * chr_bins_t functions
 ***********************************************************************/

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
    for (i = 0; i < MAX_BIN; ++i)
        vcfr_list_free(&bins->bins[i]);
}

/************************************************************************
 * g_var_t functions
 ***********************************************************************/

g_var_t *g_var_alloc(){
    g_var_t *gv = (g_var_t*)calloc(1, sizeof(g_var_t));
    
    if (gv == NULL){
        err_msg(-1, 0, "g_var_alloc: %s", strerror(errno));
        return NULL;
    }

    gv->chrm_ix = init_str_map();
    if (gv->chrm_ix == NULL)
        return(NULL);

    mv_init(&gv->chr_bins);

    mv_init(&gv->vix2var);

    return gv;
}

int g_var_add_chr(g_var_t *gv, const char *chr){
    if (gv == NULL || chr == NULL)
        return err_msg(-1, 0, "g_var_add_chr: argument is null");

    // add to the index
    int found = 0, chr_ix, bs = (int)mv_size(&gv->chr_bins);
    chr_ix = add2str_map(gv->chrm_ix, chr, &found);
    if (chr_ix < 0)
        return(-1);

    // create chr_bins if not present
    if (chr_ix >= bs){
        chr_bins_t *p = chr_bins_alloc();
        if (p == NULL)
            return(-1);

        int ret = mv_insert(cbin_vec, &gv->chr_bins, p, chr_ix);
        if (ret < 0)
            return err_msg(-1, 0, "g_var_add_chr: failed to add chr bin");
    }

    return(chr_ix);
}

int g_var_add_var(g_var_t *gv, bcf1_t *b, const bcf_hdr_t *hdr){
    // create var_t object from bcf1_t
    var_t tv;
    tv.b = bcf_dup(b);
    if (tv.b == NULL)
        return err_msg(-1, 0, "g_var_add_var: %s", strerror(errno));
    tv.vix = mv_size(&gv->vix2var);
    bcf_unpack(tv.b, BCF_UN_FLT); // unpack up to and including info field

    int bin = get_bcf1_bin(b);
    int rid = tv.b->rid;
    const char *seqname = bcf_hdr_id2name(hdr, rid);

    // add chromosome ID
    int chr_ix;
    if ( (chr_ix = g_var_add_chr(gv, seqname)) < 0)
        return(-1);

    // add to chr_bins
    // insert with 0, 1 allows and adds duplicate variants.
    // throw a warning if duplicates are found
    ml_vcfr_list_t *ll = &mv_i(&gv->chr_bins, chr_ix)->bins[bin];
    int dup;
    if ( (dup = vcfr_list_insert(ll, tv, 0, 1)) < 0 )
        return err_msg(-1, 0, "g_var_add_var: failed to insert to list");
    if (dup == 1)
        err_msg(0, 1, "variant %s is in a duplicate position. "
                "multiallelic variants should be removed.", b->d.id);

    // add to vix2var
    if (mv_push(vcfr_vec, &gv->vix2var, tv) < 0)
        return err_msg(-1, 0, "g_var_add_var: failed to add var to vix2var");

    return 0;
}

g_var_t *g_var_read_vcf(bcf_srs_t *sr, bcf_hdr_t *vcf_hdr, int max_miss, double maf_cut){
    g_var_t *gv = g_var_alloc();
    if (gv == NULL) return NULL;

    gv->vcf_hdr = vcf_hdr;

    // add chromosomes from header to str map
    if (bcf_hdr_chr_ix(vcf_hdr, gv->chrm_ix) < 0) return(NULL);
    assert(gv->chrm_ix->n >= 0);
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
            err_msg(-1, 0, "g_var_read_vcf: failed to add chr bin");
            return(NULL);
        }
    }

    while ( bcf_sr_next_line(sr) ){
        bcf1_t *vcf_r = bcf_sr_get_line(sr, 0);
        if (vcf_r == NULL) fprintf(stderr, "warning: a VCF record is NULL\n");

        // skip multiallelic snps and indels
        if (is_biallelic_snp(vcf_r) < 0) continue;

        int na = 0, nmiss=0, ntotal=1;
        if (bcf1_t_miss_maf(vcf_hdr, vcf_r, &nmiss, &na, &ntotal) < 0)
            return(NULL);

        if (max_miss >= 0 && nmiss > max_miss) continue;
        
        double maf = ((double)na) / ((double)ntotal);
        if (maf_cut >= 0.0 && (maf <= maf_cut || maf >= (1.0 - maf_cut))) continue;

        if ( g_var_add_var(gv, vcf_r, vcf_hdr) < 0 ) return NULL;
    }

    return gv;
}

void g_var_free_bcfr(g_var_t *gv){
    if (gv == NULL) return;
    size_t vix, nvar = mv_size(&gv->vix2var);
    for (vix = 0; vix < nvar; ++vix){
        var_t var = mv_i(&gv->vix2var, vix);
        if (var.b) bcf_destroy(var.b);
    }
}

void g_var_free(g_var_t *gv){
    if (gv == NULL) return;

    destroy_str_map(gv->chrm_ix);

    g_var_free_bcfr(gv);

    size_t cix, nchr = mv_size(&gv->chr_bins);
    for (cix = 0; cix < nchr; ++cix){
        chr_bins_t *p = mv_i(&gv->chr_bins, cix);
        chr_bins_free(p);
        free(p);
    }
    mv_free(&gv->chr_bins);

    mv_free(&gv->vix2var);
}

void g_var_dstry(g_var_t *gv){
    if (gv == NULL) return;
    g_var_free(gv);
    free(gv);
}

float *bcf1_ap_gt(bcf_hdr_t *vcf_hdr, bcf1_t *b, int extra){

    bcf_unpack(b, BCF_UN_FMT);

    if (b->n_allele > 2){
        err_msg(-1, 0, "bcf1_ap_gt: only biallelic variants supported, %s", b->d.id);
        return NULL;
    }

    int n_samples = bcf_hdr_nsamples(vcf_hdr);
    if (n_samples == 0) return NULL;

    // call htslib function 
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
    int32_t *p = (int32_t *)gt_arr;
    int n_miss = 0;
    for (s = 0; s < n_samples; s++){ // each sample
        float gt = 0;
        int t, total = 0;
        int32_t *a = p + s*n_val;
        for (t = 0; t < n_val; t++){ // each allele

            // if end, sample has smaller ploidy, so continue
            if ( a[t] == bcf_int32_vector_end ) continue;

            // if any missing
            if ( bcf_gt_is_missing(a[t]) ) {
                n_miss++;
                gt = -1.0;
                total = 1;
                break;
            }

            gt += (float)bcf_gt_allele(a[t]);
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

    if (b->n_allele > 2){
        err_msg(-1, 0, "bcf1_ap_gp: only biallelic variants supported, %s", b->d.id);
        return NULL;
    }

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

int g_var_get_region_vars(g_var_t *gv, const char* ref, int32_t beg, 
        int32_t end, ml_t(vcfr_list) *vars){
    if (gv == NULL || ref == NULL || vars == NULL)
        return err_msg(-2, 0, "g_var_get_region_vars: "
                "argument is null");

    int tid = str_map_ix(gv->chrm_ix, (char *)ref);
    if (tid < 0)
        return err_msg(-1, 1, "region_vars: Chromosome %s not found", ref);

    if (end <= beg)
    if (beg < 0 || end < 0 || end <= beg)
        return err_msg(-2, 0, "region_vars: end (%i) < beg (%i)", beg, end);

    // get the bins that the region overlaps
    // region must be [beg, end) 0-based.
    uint16_t list[MAX_BIN];
    int n_bin = reg2bins((int)beg, (int)end, list);

    int n_vars = 0;

    int bin_ix;
    for (bin_ix = 0; bin_ix < n_bin; ++bin_ix){
        uint16_t bin = list[bin_ix];
        ml_vcfr_list_t *ll = &mv_i(&gv->chr_bins, tid)->bins[bin];

        ml_node_t(vcfr_list) *vn;
        for (vn = ml_begin(ll); vn; vn = ml_node_next(vn)){
            var_t var = ml_node_val(vn);
            bcf1_t *b = var.b;
            if (b == NULL){
                return err_msg(-2, 0, "g_var_get_region_vars: "
                        "bcf record in variant is null");
            }

            // variants are sorted, so skip if past query region.
            if (b->pos > end)
                break;

            // see if variant overlaps region
            int ovrlp = bp_overlap(beg, end, '.', b->pos, b->pos + b->rlen, '.');
            if (ovrlp < 0){
                return err_msg(-2, 0, "g_var_get_region_vars: "
                        "failed to get bp_overlap");
            }
            if (ovrlp == 0) continue;

            // TODO how to handle multiallelic variants
            // Add var_t to list
            int added;
            if ( (added = ml_insert(vcfr_list, vars, var, 0, 1)) < 0){
                return err_msg(-2, 0, "g_var_get_region_vars: "
                        "failed to add variant to return list");
            }
            ++n_vars;
        }
    }

    return n_vars;
}

var_t *gv_vari(g_var_t *gv, int32_t ix){
    if (gv == NULL){
        err_msg(-1, 0, "gv_vari: gv is NULL");
        return(NULL);
    }
    size_t n_var = mv_size(&gv->vix2var);
    if (ix < 0){
        err_msg(-1, 0, "gv_vari: ix=%i must be within [0,%zu)", ix, n_var);
        return(NULL);
    }
    size_t s_ix = (size_t)ix;
    if (s_ix >= n_var){
        err_msg(-1, 0, "gv_vari: ix=%i must be within [0,%zu)", s_ix, n_var);
        return(NULL);
    }
    var_t *var = &mv_i(&gv->vix2var, s_ix);
    return(var);
}

