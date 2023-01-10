
#ifndef CLOPTS_H
#define CLOPTS_H

#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "sam_read.h"
#include "sam_read.h"
#include "variants.h"
#include "region.h"
#include "gtf_anno.h"

/*******************************************************************************
 * structure to hold options passed
 ******************************************************************************/

typedef struct {
    // RNA BAM
    char *rna_bam_fn;

    // ATAC BAM
    char *atac_bam_fn;

    // VCF
    char *vcf_fn;

    // GTF
    char *gtf_fn;

    // peaks
    char *peaks_fn;

    // exclusion regions
    char *exclude_fn;

    // out
    char *out_fn;

    // samples
    char *sample_fn;

    // filtered barcodes, if any
    uint32_t flt_n_bcs;
    char *flt_bc_fn; // flt barcoes
    char *wl_bc_fn; // whitelist of barcodes

    // options
    int out_min;
    int min_phred;
    int max_nh;
    char rna_bc_tag[3];
    char atac_bc_tag[3];
    char rna_umi_tag[3];
    char rna_nh_tag[3];
    char atac_nh_tag[3];
    int rna_mapq;
    int atac_mapq;
    int tx_basic;
    int counts_only;
    char *region;
    int region_set;

    int verbose;
} cl_opts;

/* Initialize and destroy cl_opts */
cl_opts *init_cl_opts();
void dstry_cl_opts(cl_opts *opts);

/*******************************************************************************
 * structure to hold initialized objects for main functions
 ******************************************************************************/

typedef struct {
    // RNA BAM
    samFile *rna_bam;
    sam_hdr_t *rna_bam_hdr;
    hts_idx_t *rna_bam_idx;

    // ATAC BAM
    samFile *atac_bam;
    sam_hdr_t *atac_bam_hdr;
    hts_idx_t *atac_bam_idx;

    // VCF
    bcf_srs_t *sr;
    bcf_hdr_t *vcf_hdr;
    g_var_t *gv;
    contig_map *cmap;

    // GTF
    gene_anno_t *anno;

    // ATAC peaks
    iregs_t *pks;

    // Blacklist regions
    iregs_t *exclude;

    // Barcodes
    uint32_t flt_n_bcs;
    str_map *flt_bcs;
    str_map *wl_bcs;

    // options
    int out_min;
    int max_nh;
    char rna_bc_tag[3];
    char atac_bc_tag[3];
    char rna_umi_tag[3];
    char rna_nh_tag[3];
    char atac_nh_tag[3];
    char *region;
    int region_set;
    int min_phred;
    int rna_mapq;
    int atac_mapq;

    char *out_fn;

    int verbose;

    str_map *samples;

} obj_pars;

/* Set default obj_pars
 *
 * Allocate and obj_pars object and initialize to default values.
 *
 * @return Pointer to obj_pars, or NULL on error.
 */
obj_pars *init_obj_pars();
void dstry_obj_pars(obj_pars *objs);

/* Load RNA BAM file
 *
 * If the file name is not given, don't do anything and return 0.
 * If there was an error loading the BAM file, return -1.
 * If loaded successfully, return 1.
 *
 * The bam, bam_hdr, and bam_idx fields are updated in @p op
 *
 * @return -1 on error, 0 if nothing happened, 1 if loaded successfully.
 */
int load_rna_bam(cl_opts *opts, obj_pars *objs);

/* Load ATAC BAM file
 *
 * If the file name is not given, don't do anything and return 0.
 * If there was an error loading the BAM file, return -1.
 * If loaded successfully, return 1.
 *
 * The bam, bam_hdr, and bam_idx fields are updated in @p op
 *
 * @return -1 on error, 0 if nothing happened, 1 if loaded successfully.
 */
int load_atac_bam(cl_opts *opts, obj_pars *objs);

/* Load the VCF file into g_var_t object.
 *
 * Loads the vcf file, subsets samples, creates a g_var_t object, and 
 * creates a cmap object.
 *
 * @return -1 on error, 0 if nothing happened, 1 if loaded successfully.
 */
int load_vars(cl_opts *opts, obj_pars *objs);

/* Load samples from file.
 *
 */
int load_samples(cl_opts *opts, obj_pars *objs);

/* Load barcodes
 *
 */
int load_bcs(cl_opts *opts, obj_pars *objs);

/* Load the GTF file into the gene_anno_t object.
 *
 * Calls read_from_gtf
 *
 * @return 0 on success, -1 on error.
 */
int load_gtf(cl_opts *opts, obj_pars *objs);

/* Load BED file into peaks.
 *
 * If the bed file is not provided (peaks_fn is null), do nothing and 
 * return 0.
 *
 * @return 0 on success, -1 on error.
 */
int load_peaks(cl_opts *opts, obj_pars *objs);

/* Load BED file of regions to exclude
 *
 * If the bed file is not provided (peaks_fn is null), do nothing and 
 * return 0.
 *
 * @return 0 on success, -1 on error.
 */
int load_exclude(cl_opts *opts, obj_pars *objs);

/* copy options from opts to obj */
int copy_options(cl_opts *opts, obj_pars *objs);

#endif // CLOPTS_H
