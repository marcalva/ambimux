
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <getopt.h>
#include <errno.h>
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/klist.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "htslib/kstring.h"
#include "htslib/khash_str2int.h"
#include "sam_read.h"
#include "str_util.h"
#include "variants.h"
#include "gtf_anno.h"
#include "overlap.h"
#include "atac_data.h"
#include "rna_data.h"
#include "counts.h"
#include "region.h"
#include "clopts.h"
#include "bc_stats.h"
#include "mod.h"
#include "bam_dat.h"
#include "r_count.h"
#include "bam_proc.h"

#define Er(x, r) \
    if ((x) < 0){ \
        r = EXIT_FAILURE; \
        goto cleanup; \
    }

static void usage(FILE *fp, int exit_status){
    fprintf(fp, 
            "\n"
            "ambimux v0.1.0: single-cell demultiplexing\n"
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
            "  -e, --exclude       BED file containing regions to exclude.\n"
            "  -s, --samples       File listing sample IDs in VCF file to de-multiplex.\n"
            "\n"
            "Output options\n"
            "\n"
            "  -o, --out           Output file prefix [ambimux]. File names are appended with '.' delimiter\n"
            "  -C, --counts-only   Flag argument to produce counts only and not fit a demultiplexing model.\n"
            "  -x, --out-min       Calculate the likelihood and demultiplex barcodes that have at least this many \n"
            "                      RNA UMIs or ATAC fragments. If there both the UMIs and fragments are below this, skip. [100]\n"
            "\n"
            "Alignment options\n"
            "\n"
            "  -f, --flt-bcs       File containing list of filtered barcodes for alpha estimation\n"
            "                      These correspond to nonempty droplets. If set, --flt-n is ignored.\n"
            "                      Likelihood and summary results will be output only for these filtered barcodes,\n"
            "                      whereas counts will be generated for all barcodes.\n"
            "  -n, --flt-n         Alternative to flt-bcs, estimate alpha by assuming the top n barcodes\n"
            "                      ranked by number or RNA+ATAC reads are non-empty [2000].\n"
            "  -w, --bc-wl         Optional file containing a whitelist of barcodes.\n"
            "  -u, --rna-umi-tag   RNA BAM tag for UMI [UB].\n"
            "  -b, --rna-bc-tag    RNA BAM tag for cell barcode [CB].\n"
            "  -c, --atac-bc-tag   ATAC BAM tag for cell barcode [CB].\n"
            "  -H, --nh-tag [type] BAM tag for the number of alignments of a read, where type is one \n"
            "                      of 'rna' or 'atac'  [NH].\n"
            "\n"
            "Mapping thresholds\n"
            "\n"
            "  -m, --max-nh        Only process reads with a maximum of this many alignments [0]. \n"
            "                      Set to 0 to ignore (default).\n"
            "  -P, --min-phredq    Minimum base phred quality score in read to consider for variant calling [30].\n"
            "  -z, --rna-mapq      Minimum MAPQ (mapping quality) of an RNA alignment [30].\n"
            "  -Z, --atac-mapq     Minimum MAPQ (mapping quality) of an ATAC alignment [30].\n"
            "  -R, --region        Region (hts format), for example 'chr21,chr21:10-,chr21-10-20'.\n"
            "\n"
            "GTF options\n"
            "\n"
            "  -t, --tx-basic      Read only transcripts tagged as 'basic' in the GTF file.\n"
            "\n"
            "  -V, --verbose       Write status on output.\n"
            "\n");
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
        {"out-min", required_argument, NULL, 'x'}, 
        {"flt-bcs", required_argument, NULL, 'f'}, 
        {"flt-n-bcs", required_argument, NULL, 'n'}, 
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
        {"tx-basic", no_argument, NULL, 't'}, 
        {"verbose", no_argument, NULL, 'V'},
        {NULL, 0, NULL, 0}
    };

    cl_opts *opts = init_cl_opts();
    obj_pars *objs = init_obj_pars();

    /* parameters */
    int ret = EXIT_SUCCESS;
    /* */

    /* local variables */
    bam_data_t *bam_dat = NULL;
    mdl_t *mdl = NULL;
    /* */

    char *p_end = NULL;
    int option_index = 0;
    int cm;
    while ((cm = getopt_long_only(argc, argv, "a:r:v:g:p:e:s:oC:x:f:n:w:u:b:c:H:m:P:z:Z:tR:V", loptions, &option_index)) != -1){
        switch(cm){
            case 'a': opts->atac_bam_fn = strdup(optarg);
                      if (opts->atac_bam_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'r': opts->rna_bam_fn = strdup(optarg);
                      if (opts->rna_bam_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'v': opts->vcf_fn = strdup(optarg); 
                      if (opts->vcf_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'g': opts->gtf_fn = strdup(optarg); 
                      if (opts->vcf_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'p': opts->peaks_fn = strdup(optarg); 
                      if (opts->peaks_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'e': opts->exclude_fn = strdup(optarg); 
                      if (opts->exclude_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 's': opts->sample_fn = strdup(optarg);
                      if (opts->sample_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'o': free(opts->out_fn);
                      opts->out_fn = strdup(optarg); 
                      if (opts->out_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'C': opts->counts_only = 1;
                      break;
            case 'x': errno = 0;
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
            case 'f': opts->flt_bc_fn = strdup(optarg);
                      if (opts->flt_bc_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'n': errno = 0;
                      opts->flt_n_bcs = (uint32_t)strtoul(optarg, &p_end, 10);
                      int eno = errno;
                      if (eno == EINVAL || eno == ERANGE){
                          ret = err_msg(EXIT_FAILURE, 0, 
                                  "could not convert --flt_n_bcs %s to int: %s", 
                                  optarg, strerror(errno));
                          goto cleanup;
                      }
                      break;
            case 'w': opts->wl_bc_fn = strdup(optarg);
                      if (opts->wl_bc_fn == NULL){
                          ret = err_msg(EXIT_FAILURE, 0, "%s", strerror(errno));
                          goto cleanup;
                      }
                      break; 
            case 'u': if (strlen(optarg) != 2){
                          ret = err_msg(EXIT_FAILURE, 0, 
                                  "rna-umi-tag must be a 2 character string: %s given", optarg);
                          goto cleanup;
                      }
                      strcpy(opts->rna_umi_tag, optarg);
                      break;
            case 'b': if (strlen(optarg) != 2){
                          ret = err_msg(EXIT_FAILURE, 0, 
                                  "rna-bc-tag must be a 2 character string: %s given", optarg);
                          goto cleanup;
                      }
                      strcpy(opts->rna_bc_tag, optarg);
                      break;
            case 'c': if (strlen(optarg) != 2){
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
            case 'm': errno = 0;
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
            case 't':
                      opts->tx_basic = 1;
                      break;
            case 'R': free(opts->region);
                      opts->region = strdup(optarg);
                      opts->region_set = 1;
                      break;
            case 'V': opts->verbose = 1;
                      break;
            default: 
                      err_msg(EXIT_FAILURE, 0, 
                              "unrecognized option: %s", loptions[option_index].name);
                      usage(stdout, EXIT_FAILURE);
        }
    }

    copy_options(opts, objs);

    /* Load in BAM files */
    ret = load_atac_bam(opts, objs);
    if (ret < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    ret = load_rna_bam(opts, objs);
    if (ret < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* Load samples */
    ret = load_samples(opts, objs);
    if (ret < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* Load in VCF file */
    ret = load_vars(opts, objs);
    if (ret < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* Load GTF file */
    ret = load_gtf(opts, objs);
    if (ret < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* Load peaks */
    if (load_peaks(opts, objs) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    /* load exclusion */
    if (load_exclude(opts, objs) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }
    
    // Load/set barcodes
    ret = load_bcs(opts, objs);
    if (ret < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    // init bam data object.
    if ( (bam_dat = bam_data_init()) == NULL ){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    // Pileup for RNA and ATAC:
    // Run ATAC
    if (objs->atac_bam != NULL){
        ret = run_atac(objs, bam_dat);
        if (ret < 0){
            ret = EXIT_FAILURE;
            goto cleanup;
        }
    }

    // Run RNA
    if (objs->rna_bam != NULL){
        ret = run_rna(objs, bam_dat);
        if (ret < 0){
            ret = EXIT_FAILURE;
            goto cleanup;
        }
    }

    if (bam_data_fill_bcs(bam_dat, objs->wl_bcs) < 0)
        goto cleanup;
    if (bam_data_fill_stats(bam_dat) < 0)
        goto cleanup;

    // generate gene counts
    if (bam_count(bam_dat, objs, objs->out_fn) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

    if (opts->counts_only)
        goto cleanup;

    // fit model
    if (fit_mdl_pars(bam_dat, objs) < 0){
        ret = EXIT_FAILURE;
        goto cleanup;
    }

cleanup:

    if ( ret == EXIT_SUCCESS && opts->verbose) log_msg("done");

    bam_data_dstry(bam_dat);
    dstry_obj_pars(objs);
    dstry_cl_opts(opts);
    mdl_dstry(mdl);

    return(ret);
}

