
#include <stdlib.h>
#include <string.h>
#include "clopts.h"

cl_opts *init_cl_opts(){
    cl_opts *opts = (cl_opts *)calloc(1, sizeof(cl_opts));
    if (opts == NULL){
        err_msg(-1, 0, "init_cl_opts: %s", strerror(errno));
        return(NULL);
    }

    // file names
    opts->rna_bam_fn = NULL;
    opts->atac_bam_fn = NULL;
    opts->vcf_fn = NULL;
    opts->gtf_fn = NULL;
    opts->exclude_fn = NULL;
    opts->peaks_fn = NULL;
    opts->out_fn = strdup("ambimux");
    opts->sample_fn = NULL;

    opts->wl_bc_fn = NULL;

    // parameters
    opts->out_min = 100;
    opts->min_phred = 30;
    opts->max_nh = 0;
    opts->k = 1;

    strcpy(opts->rna_bc_tag, "CB");
    strcpy(opts->atac_bc_tag, "CB");
    strcpy(opts->rna_umi_tag, "UB");
    strcpy(opts->rna_nh_tag, "NH");
    strcpy(opts->atac_nh_tag, "NH");

    opts->rna_mapq = 30;
    opts->atac_mapq = 30;
    opts->tx_basic = 0;
    opts->counts_only = 0;
    opts->no_counts_o = 1;
    // opts->alpha_vars = 1;

    opts->region = strdup(".");
    opts->region_set = 0;

    opts->eps = 1e-6;
    opts->max_iter = 100;
    // opts->alpha_max = 1;

    opts->mdl_intra_reads = 1;
    opts->mdl_inter_reads = 0;

    opts->alpha_prior_w = 1.0;

    opts->threads = 1;

    opts->verbose = 0;

    return(opts);
}

void dstry_cl_opts(cl_opts *opts){
    free(opts->rna_bam_fn);
    free(opts->atac_bam_fn);
    free(opts->vcf_fn);
    free(opts->gtf_fn);
    free(opts->exclude_fn);
    free(opts->peaks_fn);
    free(opts->out_fn);
    free(opts->wl_bc_fn);
    free(opts->sample_fn);
    free(opts->region);

    free(opts);
}

/*******************************************************************************
 * structure to hold initialized objects for main functions
 ******************************************************************************/

obj_pars *init_obj_pars(){
    obj_pars *p = (obj_pars *)calloc(1, sizeof(obj_pars));
    if (p == NULL){
        err_msg(-1, 0, "init_obj_pars: %s", strerror(errno));
        return(NULL);
    }

    p->rna_bam = NULL;
    p->rna_bam_hdr = NULL;
    p->rna_bam_idx = NULL;

    p->atac_bam = NULL;
    p->atac_bam_hdr = NULL;
    p->atac_bam_idx = NULL;

    p->sr = NULL;
    p->vcf_hdr = NULL;
    p->gv = NULL;
    p->cmap = NULL;

    p->pks = NULL;

    p->exclude = NULL;

    p->anno = NULL;

    p->wl_bcs = NULL;

    p->out_min = 100;
    p->max_nh = 0;
    p->k = 1;

    strcpy(p->rna_bc_tag, "CB");
    strcpy(p->atac_bc_tag, "CB");
    strcpy(p->rna_umi_tag, "UB");
    strcpy(p->rna_nh_tag, "NH");
    strcpy(p->atac_nh_tag, "NH");

    p->region = strdup(".");
    p->region_set = 0;

    p->min_phred = 30;
    p->rna_mapq = 30;
    p->atac_mapq = 30;
    // p->alpha_vars = 1;

    p->eps = 1e-6;
    p->max_iter = 100;

    p->mdl_intra_reads = 1;
    p->mdl_inter_reads = 0;

    p->alpha_prior_w = 1.0;

    p->threads = 1;

    p->out_fn = NULL;

    p->verbose = 0;

    p->samples = NULL;

    return(p);
}

void dstry_obj_pars(obj_pars *objs){
    if (objs == NULL) return;

    sam_hdr_destroy(objs->rna_bam_hdr);
    if (objs->rna_bam) sam_close(objs->rna_bam);
    if (objs->rna_bam_idx) hts_idx_destroy(objs->rna_bam_idx);

    sam_hdr_destroy(objs->atac_bam_hdr);
    if (objs->atac_bam) sam_close(objs->atac_bam);
    if (objs->atac_bam_idx) hts_idx_destroy(objs->atac_bam_idx);

    if (objs->sr) bcf_sr_destroy(objs->sr);
    // if (objs->vcf_hdr) bcf_hdr_destroy(objs->vcf_hdr);

    if (objs->anno) gene_anno_dstry(objs->anno);

    if (objs->pks) iregs_dstry(objs->pks);

    if (objs->exclude) iregs_dstry(objs->exclude);

    if (objs->wl_bcs) destroy_str_map(objs->wl_bcs);

    if (objs->region) free(objs->region);

    g_var_dstry(objs->gv);
    
    destroy_str_map(objs->cmap);

    if (objs->out_fn) free(objs->out_fn);

    if (objs->samples) destroy_str_map(objs->samples);

    free(objs);
}

int load_rna_bam(cl_opts *opts, obj_pars *objs){
    if (opts == NULL || objs == NULL)
        return err_msg(-1, 0, "load_rna_bam: 'opts' or 'objs' is null");
    if (opts->rna_bam_fn == NULL) return(0);

    int ret = load_bam(opts->rna_bam_fn, &objs->rna_bam, &objs->rna_bam_hdr, &objs->rna_bam_idx);
    if (ret < 0) return(-1);
    else return(1);
}

int load_atac_bam(cl_opts *opts, obj_pars *objs){
    if (opts == NULL || objs == NULL)
        return err_msg(-1, 0, "load_atac_bam: 'opts' or 'objs' is null");
    if (opts->atac_bam_fn == NULL) return(0);

    int ret = load_bam(opts->atac_bam_fn, &objs->atac_bam, &objs->atac_bam_hdr, &objs->atac_bam_idx);
    if (ret < 0) return(-1);
    else return(1);
}

int load_vars(cl_opts *opts, obj_pars *objs){
    if (opts == NULL || objs == NULL)
        return err_msg(-1, 0, "load_vars: 'opts' or 'objs' is null");
    if (opts->vcf_fn == NULL) return(0);

    int ret;

    if (opts->verbose)
        log_msg("reading VCF file");

    /* Load in VCF file */
    ret = load_vcf(opts->vcf_fn, opts->region, opts->region_set, &objs->sr, &objs->vcf_hdr);
    if (ret < 0){
        return err_msg(-1, 0, "load_vars: failed to load vcf file %s", opts->vcf_fn);
    }

    /* subset to samples */
    ret = sub_vcf_samples(&objs->vcf_hdr, opts->sample_fn);
    if (ret < 0){
        return err_msg(-1, 0, "load_vars: failed to subset samples in vcf file %s "
                "from %s", opts->vcf_fn, opts->sample_fn);
    }

    objs->gv = g_var_read_vcf(objs->sr, objs->vcf_hdr, 0, 0);
    if (objs->gv == NULL) return(-1);

    // chr map
    objs->cmap = init_str_map();
    if (objs->cmap == NULL) return(-1);

    ret = bcf_hdr_chr_ix(objs->vcf_hdr, objs->cmap);
    if (ret < 0) return(-1);

    /*
    int cmi = 0;
    for (cmi = 0; cmi < cmap->chr_map->n; ++cmi){
        fprintf(stdout, "%i=%s\n", cmi, str_map_str(cmap->chr_map, cmi));
    }
    */

    return(1);
}

int load_gtf(cl_opts *opts, obj_pars *objs){
    if (opts == NULL || objs == NULL)
        return err_msg(-1, 0, "load_gtf: 'opts' or 'objs' is null");

    if (opts->gtf_fn == NULL)
        return(0);

    if (opts->verbose)
        log_msg("reading GTF file");

    objs->anno = read_from_gtf(opts->gtf_fn, opts->tx_basic);
    if (objs->anno == NULL)
        return err_msg(-1, 0, "load_gtf: failed to read GTF file");

    return(0);
}

int load_samples(cl_opts *opts, obj_pars *objs){
    if (opts == NULL || objs == NULL)
        return err_msg(-1, 0, "load_samples: 'opts' or 'objs' is null");

    if (opts->sample_fn == NULL)
        return(0);

    objs->samples = read_str_map(opts->sample_fn);
    if (objs->samples == NULL)
        return err_msg(-1, 0, "load_samples: failed to read samples "
                "from %s", opts->sample_fn);

    return(0);
}

int load_bcs(cl_opts *opts, obj_pars *objs){
    if (opts == NULL || objs == NULL)
        return err_msg(-1, 0, "load_bcs: 'opts' or 'objs' is null");

    if (opts->verbose)
        log_msg("reading barcode whitelist");

    if (opts->wl_bc_fn != NULL){
        objs->wl_bcs = read_str_map(opts->wl_bc_fn);
        if (objs->wl_bcs == NULL)
            return err_msg(-1, 0, "load_bcs: failed to read wl_bcs "
                    "from %s", opts->wl_bc_fn);
    }

    return(0);
}

int load_peaks(cl_opts *opts, obj_pars *objs){
    if (opts == NULL || objs == NULL)
        return err_msg(-1, 0, "load_peaks: 'opts' or 'objs' is null");

    if (opts->peaks_fn == NULL) return(0);

    if (opts->verbose) log_msg("reading peaks");

    if ((objs->pks = iregs_init()) == NULL) return(-1);
    if (iregs_add_bed(objs->pks, opts->peaks_fn) < 0) return(-1);
    if (iregs_parse_bed(objs->pks) < 0) return(-1);

    return(0);
}

int load_exclude(cl_opts *opts, obj_pars *objs){
    if (opts == NULL || objs == NULL)
        return err_msg(-1, 0, "load_exclude: 'opts' or 'objs' is null");

    if (opts->exclude_fn == NULL) return(0);

    if (opts->verbose) log_msg("reading exclusion regions");

    if ((objs->exclude = iregs_init()) == NULL) return(-1);
    if (iregs_add_bed(objs->exclude, opts->exclude_fn) < 0) return(-1);
    if (iregs_parse_bed(objs->exclude) < 0) return(-1);

    return(0);
}

int copy_options(cl_opts *opts, obj_pars *objs){
    if (opts == NULL || objs == NULL)
        return err_msg(-1, 0, "copy_options: 'opts' or 'objs' is null");
    objs->out_min = opts->out_min;
    objs->max_nh = opts->max_nh;
    strcpy(objs->rna_bc_tag, opts->rna_bc_tag);
    strcpy(objs->atac_bc_tag, opts->atac_bc_tag);
    strcpy(objs->rna_umi_tag, opts->rna_umi_tag);
    strcpy(objs->rna_nh_tag, opts->rna_nh_tag);
    free(objs->region);
    if (opts->region) objs->region = strdup(opts->region);
    objs->region_set = opts->region_set;
    objs->min_phred = opts->min_phred;
    objs->rna_mapq = opts->rna_mapq;
    objs->atac_mapq = opts->atac_mapq;
    // objs->alpha_vars = opts->alpha_vars;
    objs->eps = opts->eps;
    objs->k = opts->k;
    objs->max_iter = opts->max_iter;
    objs->alpha_prior_w = opts->alpha_prior_w;
    objs->threads = opts->threads;
    objs->mdl_intra_reads = opts->mdl_intra_reads;
    objs->mdl_inter_reads = opts->mdl_inter_reads;
    if (opts->out_fn) objs->out_fn = strdup(opts->out_fn);
    objs->verbose = opts->verbose;

    return(0);
}
