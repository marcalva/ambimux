
/*
#include <string.h>
*/
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <getopt.h>
#include <errno.h>
#include "mod.h"
#include "bam_proc.h"

#define HOPT 1000

#define ER(ret) \
    if ((ret) < 0){ \
        ret = EXIT_FAILURE; \
        goto cleanup; \
    } \

char am_version_str[] = "v0.4.0";

static void usage(FILE *fp, int exit_status){
    fprintf(fp, 
            "\n"
            "ambimux %s: single-cell demultiplexing\n"
            "\n"
            "Options:\n"
            "\n"
            "Input options\n"
            "\n"
            "  -a, --atac-bam      Indexed ATAC BAM file.\n"
            "  -r, --rna-bam       Indexed RNA BAM file.\n"
            "  -v, --vcf           Indexed VCF file.\n" 
            "  -g, --gtf           GTF file.\n"
            "  -p, --peaks         BED file containing peaks.\n"
            "  -e, --exclude       BED file containing ATAC regions to exclude.\n"
            "  -s, --samples       File listing sample IDs in VCF file to de-multiplex.\n"
            "\n"
            "Output options\n"
            "\n"
            "  -o, --out           Output file prefix. File names are appended with '.' delimiter [ambimux].\n"
            "  -C, --counts-only   Flag argument to produce counts only and not fit a demultiplexing model.\n"
            "  -n, --no-counts-o   Do not output variant allele, gene, and peak counts. by default, ambimux\n"
            "                      output UMI/fragment counts.\n"
            "\n"
            "Alignment options\n"
            "\n"
            "  -w, --bc-wl         Optional file containing a whitelist of barcodes. One barcode per line.\n"
            "  -u, --rna-umi-tag   RNA BAM tag for UMI [UB].\n"
            "  -b, --rna-bc-tag    RNA BAM tag for cell barcode [CB].\n"
            "  -c, --atac-bc-tag   ATAC BAM tag for cell barcode [CB].\n"
            "  -H, --nh-tag [type] BAM tag for the number of alignments of a read, where type is one \n"
            "                      of 'rna' or 'atac'  [NH].\n"
            "\n"
            "Mapping thresholds\n"
            "\n"
            "  -m, --max-nh        Only process reads with a maximum of this many alignments. \n"
            "                      Set to 0 to ignore (default) [0].\n"
            "  -P, --min-phredq    Minimum base phred quality score in read to consider for variant calling [30].\n"
            "  -z, --rna-mapq      Minimum MAPQ (mapping quality) of an RNA alignment [30].\n"
            "  -Z, --atac-mapq     Minimum MAPQ (mapping quality) of an ATAC alignment [30].\n"
            "  -R, --region        Region (hts format), for example 'chr21,chr21:10-,chr21-10-20'.\n"
            "\n"
            "Model options\n"
            "\n"
            "  -x, --out-min       Calculate the likelihood and demultiplex barcodes that have at least this many \n"
            "                      RNA UMIs or ATAC fragments. If there both the UMIs and fragments are below this, skip [100].\n"
            "  -h, --eps           Convergence threshold, where the percent change in parameters\n"
            "                      must be less than this value [1e-6].\n"
            "  -j, --amb-eps       Convergence threshold for droplet ambient contamination parameter alpha,\n"
            "                      the percent change in value must be less than this argument [1e-6].\n"
            "  -q, --max-iter      Maximum number of iterations to perform for EM [100].\n"
            "  -d, --amb-prior-w   Weight of the ambient fraction prior. Effectively gives the number\n"
            "                      of reads added to the likelihood as a prior. Must be > 0 [1].\n"
            "  -i, --mdl-reads     The type of reads to use for demultiplexing and ambient estimation.\n"
            "                      One of 'intra', 'inter', or 'all'. 'intra' (default) specifies the model\n"
            "                      should use only reads contained in peaks or genes.\n"
            "                      'inter' specifies to use only inter-peak or inter-gene reads that\n"
            "                      fall outside peaks or genes. 'all' specifies to use all reads including\n"
            "                      inter- and intra- gene and peak reads [intra].\n"
            "  -T, --threads       Optional number of threads to use [1].\n"
            "\n"
            "GTF options\n"
            "\n"
            "  -t, --tx-basic      Read only transcripts tagged as 'basic' in the GTF file.\n"
            "\n"
            "  -V, --verbose       Write status on output.\n"
            "      --help          Print this help screen.\n"
            "\n",
            am_version_str);
    exit(exit_status);
}

int main(int argc, char *argv[]){

    if (argc == 1) usage(stderr, EXIT_FAILURE);

    static const struct option loptions[] =
    {
        {"atac-bam", required_argument, NULL, 'a'}, 
        {"rna-bam", required_argument, NULL, 'r'}, 
        {"vcf", required_argument, NULL, 'v'}, 
        {"gtf", required_argument, NULL, 'g'},
        {"peaks", required_argument, NULL, 'p'},
        {"exclude", required_argument, NULL, 'e'},
        {"samples", required_argument, NULL, 's'}, 
        {"out", required_argument, NULL, 'o'},
        {"counts-only", no_argument, NULL, 'C'}, 
        {"no-counts-o", no_argument, NULL, 'n'}, 
        {"out-min", required_argument, NULL, 'x'}, 
        {"bc-wl", required_argument, NULL, 'w'}, 
        {"rna-umi-tag", required_argument, NULL, 'u'}, 
        {"rna-bc-tag", required_argument, NULL, 'b'}, 
        {"atac-bc-tag", required_argument, NULL, 'c'}, 
        {"nh-tag", required_argument, NULL, 'H'}, 
        {"max-nh", required_argument, NULL, 'm'}, 
        {"min-phredq", required_argument, NULL, 'P'},
        {"rna-mapq", required_argument, NULL, 'z'},
        {"atac-mapq", required_argument, NULL, 'Z'},
        {"region", required_argument, NULL, 'R'},
        {"eps", required_argument, NULL, 'h'},
        {"amb-eps", required_argument, NULL, 'j'},
        {"max-iter", required_argument, NULL, 'q'},
        {"amb-prior-w", required_argument, NULL, 'd'},
        {"mdl-reads", required_argument, NULL, 'i'},
        {"threads", required_argument, NULL, 'T'},
        {"tx-basic", no_argument, NULL, 't'}, 
        {"verbose", no_argument, NULL, 'V'},
        {"help", no_argument, NULL, HOPT},
        {NULL, 0, NULL, 0}
    };

    cl_opts *opts = init_cl_opts();
    obj_pars *objs = init_obj_pars();

    // parameters
    int ret = EXIT_SUCCESS;

    // local variables
    bam_data_t *bam_dat = NULL;
    mdl_t *mdl = NULL;

    char *p_end = NULL;
    int option_index = 0;
    int cm, eno;
    while ((cm = getopt_long_only(argc, argv, "a:r:v:g:p:e:s:o:Cnx:w:u:b:c:H:m:P:z:Z:R:h:j:q:d:i:T:tV", loptions, &option_index)) != -1){
        switch(cm){
            case 'a':
                opts->atac_bam_fn = strdup(optarg);
                if (opts->atac_bam_fn == NULL){
                    ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                    goto cleanup;
                }
                break; 
            case 'r':
                opts->rna_bam_fn = strdup(optarg);
                if (opts->rna_bam_fn == NULL){
                    ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                    goto cleanup;
                }
                break; 
            case 'v':
                opts->vcf_fn = strdup(optarg); 
                if (opts->vcf_fn == NULL){
                    ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                    goto cleanup;
                }
                break; 
            case 'g':
                opts->gtf_fn = strdup(optarg); 
                if (opts->vcf_fn == NULL){
                    ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                    goto cleanup;
                }
                break; 
            case 'p':
                opts->peaks_fn = strdup(optarg); 
                if (opts->peaks_fn == NULL){
                    ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                    goto cleanup;
                }
                break; 
            case 'e':
                opts->exclude_fn = strdup(optarg); 
                if (opts->exclude_fn == NULL){
                    ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                    goto cleanup;
                }
                break; 
            case 's':
                opts->sample_fn = strdup(optarg);
                if (opts->sample_fn == NULL){
                    ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                    goto cleanup;
                }
                break; 
            case 'o':
                free(opts->out_fn);
                opts->out_fn = strdup(optarg); 
                if (opts->out_fn == NULL){
                    ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                    goto cleanup;
                }
                break; 
            case 'C':
                opts->counts_only = 1;
                break;
            case 'n':
                opts->no_counts_o = 1;
                break;
            case 'x':
                errno = 0;
                opts->out_min = (int)strtol(optarg, &p_end, 10);
                if (opts->out_min == 0 && errno > 0){
                    ret =  err_msg(EXIT_FAILURE, 0, 
                                   "could not convert --out-min %s to int: %s", 
                                   optarg, strerror(errno));
                    goto cleanup;
                }
                if (opts->out_min <= 0){
                    ret = err_msg(EXIT_FAILURE, 0, "--out-min must be greater than 0"); 
                    goto cleanup;
                }
                break;
            case 'w':
                opts->wl_bc_fn = strdup(optarg);
                if (opts->wl_bc_fn == NULL){
                    ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                    goto cleanup;
                }
                break; 
            case 'u':
                if (strlen(optarg) != 2){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "rna-umi-tag must be a 2 character string: %s given", optarg);
                    goto cleanup;
                }
                strcpy(opts->rna_umi_tag, optarg);
                break;
            case 'b':
                if (strlen(optarg) != 2){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "rna-bc-tag must be a 2 character string: %s given", optarg);
                    goto cleanup;
                }
                strcpy(opts->rna_bc_tag, optarg);
                break;
            case 'c':
                if (strlen(optarg) != 2){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "atac-bc-tag must be a 2 character string: %s given", optarg);
                    goto cleanup;
                }
                strcpy(opts->atac_bc_tag, optarg);
                break;
            case 'H': 
                if (strlen(optarg) != 2){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "nh-tag must be a 2 character string: %s given", optarg);
                    goto cleanup;
                }
                strcpy(opts->rna_nh_tag, optarg);
                break;
            case 'm':
                errno = 0;
                opts->max_nh = (int)strtol(optarg, &p_end, 10);
                if (opts->max_nh == 0 && errno > 0){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "could not convert --max-nh %s to int: %s", 
                                  optarg, strerror(errno));
                    goto cleanup;
                }
                break;
            case 'P':
                errno = 0;
                opts->min_phred = (int)strtol(optarg, &p_end, 10);
                if (opts->min_phred == 0 && errno > 0){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "could not convert --min-phred %s to int: %s", 
                                  optarg, strerror(errno));
                    goto cleanup;
                }
                break;
            case 'z':
                errno = 0;
                opts->rna_mapq = (int)strtol(optarg, &p_end, 10);
                if (opts->rna_mapq == 0 && errno > 0){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "could not convert --rna-mapq %s to int: %s", 
                                  optarg, strerror(errno));
                    goto cleanup;
                }
                if (opts->rna_mapq < 0 || opts->rna_mapq > 255){
                    ret = err_msg(EXIT_FAILURE, 0, "--rna-mapq must be between 0 and 255"); 
                    goto cleanup;
                }
                break;
            case 'Z':
                errno = 0;
                opts->atac_mapq = (int)strtol(optarg, &p_end, 10);
                if (opts->atac_mapq == 0 && errno > 0){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "could not convert --atac-mapq %s to int: %s", 
                                  optarg, strerror(errno));
                    goto cleanup;
                }
                if (opts->atac_mapq < 0 || opts->atac_mapq > 255){
                    ret = err_msg(EXIT_FAILURE, 0, "--atac-mapq must be between 0 and 255"); 
                    goto cleanup;
                }
                break;
            case 'R': free(opts->region);
                opts->region = strdup(optarg);
                opts->region_set = 1;
                break;
            case 'h':
                errno = 0;
                float eps = strtof(optarg, &p_end);
                if (errno == ERANGE){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "could not convert --eps %s to int: %s", 
                                  optarg, strerror(errno));
                    goto cleanup;
                }
                if (eps <= 0 || eps >= 1){
                    ret = err_msg(EXIT_FAILURE, 0, "--eps must be between 0 and 1"); 
                    goto cleanup;
                }
                opts->eps = eps;
                break;
            case 'j':
                errno = 0;
                float alpha_eps = strtof(optarg, &p_end);
                if (errno == ERANGE){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "could not convert --amb-eps %s to int: %s", 
                                  optarg, strerror(errno));
                    goto cleanup;
                }
                if (alpha_eps <= 0 || alpha_eps >= 1){
                    ret = err_msg(EXIT_FAILURE, 0, "--amb-eps must be between 0 and 1"); 
                    goto cleanup;
                }
                opts->alpha_eps = alpha_eps;
                break;
            case 'q':
                errno = 0;
                opts->max_iter = (uint16_t)strtoul(optarg, &p_end, 10);
                eno = errno;
                if (eno == EINVAL || eno == ERANGE){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "could not convert --max-iter %s to int: %s", 
                                  optarg, strerror(errno));
                    goto cleanup;
                }
                break;
            case 'd':
                errno = 0;
                opts->alpha_prior_w = strtof(optarg, &p_end);
                if (errno == ERANGE){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "could not convert --amb-prior-w '%s' to float: %s", 
                                  optarg, strerror(errno));
                    goto cleanup;
                }
                if (opts->alpha_prior_w <= 0) {
                    ret = err_msg(EXIT_FAILURE, 0, "--amb-prior-w '%s' must > 0", optarg);
                    goto cleanup;
                }
                break;
            case 'i':
                errno = 0;
                if ( strcmp(optarg, "intra") == 0 ) {
                    opts->mdl_intra_reads = 1;
                    opts->mdl_inter_reads = 0;
                } else if ( strcmp(optarg, "inter") == 0 ) {
                    opts->mdl_intra_reads = 0;
                    opts->mdl_inter_reads = 1;
                } else if ( strcmp(optarg, "all") == 0 ) {
                    opts->mdl_intra_reads = 1;
                    opts->mdl_inter_reads = 1;
                } else {
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "invalid argument for --mdl-reads: %s", 
                                  optarg);
                    goto cleanup;
                }
                break;
            case 'T':
                errno = 0;
                opts->threads = (uint16_t)strtoul(optarg, &p_end, 10);
                eno = errno;
                if (eno == EINVAL || eno == ERANGE){
                    ret = err_msg(EXIT_FAILURE, 0, 
                                  "could not convert --threads %s to int: %s", 
                                  optarg, strerror(errno));
                    goto cleanup;
                }
                break;
            case 't':
                opts->tx_basic = 1;
                break;
            case 'V':
                opts->verbose = 1;
                break;
            case HOPT:
                usage(stdout, EXIT_SUCCESS);
                break;
            default: 
                // err_msg(EXIT_FAILURE, 0, "unrecognized option: %s", loptions[option_index].name);
                usage(stdout, EXIT_FAILURE);
                break;
        }
    }

    copy_options(opts, objs);

    if (objs->verbose) log_msg("running ambimux %s", am_version_str);

    /* Load in BAM files */
    ret = load_atac_bam(opts, objs);
    ER(ret);

    ret = load_rna_bam(opts, objs);
    ER(ret);

    /* Load samples */
    ret = load_samples(opts, objs);
    ER(ret);

    /* Load in VCF file */
    ret = load_vars(opts, objs);
    ER(ret);

    /* Load GTF file */
    ret = load_gtf(opts, objs);
    ER(ret);

    /* Load peaks */
    ret = load_peaks(opts, objs);
    ER(ret);

    /* load exclusion */
    ret = load_exclude(opts, objs);
    ER(ret);
    
    /* Load/set barcodes */
    ret = load_bcs(opts, objs);
    ER(ret);

    // init bam data object.
    if ( (bam_dat = bam_data_init()) == NULL ){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    // Pileup for RNA and ATAC:
    // Run ATAC
    if (objs->atac_bam != NULL){
        ret = run_atac(objs, bam_dat);
        ER(ret);
    }

    // Run RNA
    if (objs->rna_bam != NULL){
        ret = run_rna(objs, bam_dat);
        ER(ret);
    }

    if (objs->verbose) log_msg("getting barcode data from BAM records ");
    ret = bam_data_fill_bcs(bam_dat, objs->wl_bcs);
    ER(ret);

    ret = bam_data_fill_stats(bam_dat, objs->rna_bam_hdr, objs->atac_bam_hdr);
    ER(ret);

    // generate gene counts
    if (!opts->no_counts_o) {
        ret = bam_count(bam_dat, objs, objs->out_fn);
        ER(ret);
    }

    if (opts->counts_only)
        goto cleanup;

    // fit model
    ret = mdl_fit(bam_dat, objs);
    ER(ret);

    if ( ret == EXIT_SUCCESS && opts->verbose) log_msg("done");

cleanup:

    bam_data_dstry(bam_dat);
    dstry_obj_pars(objs);
    dstry_cl_opts(opts);
    mdl_dstry(mdl);

    return(ret);
}

