
#include <stdio.h>
#include <pthread.h>
#include "htslib/sam.h"
#include "clopts.h"
#include "atac_data.h"
#include "htslib/hts.h"
#include "rna_data.h"
#include "bam_dat.h"
#include "str_util.h"
#include "math_util.h"
#include "overlap.h"
#include "r_count.h"
#include "bam_proc.h"

int proc_atac1(bam1_t *bam_r, obj_pars *objs, bam_data_t *bam_data){
    if (bam_r == NULL)
        return err_msg(-1, 0, "proc_atac1: 'bam_r' is null");
    if (objs == NULL)
        return err_msg(-1, 0, "proc_atac1: 'objs' is null");
    if (bam_data == NULL)
        return err_msg(-1, 0, "proc_atac1: 'bam' is null");

    if (bam1_unmapped(bam_r))
        return 0; // unmapped
    if ( ((bam_r->core.flag)&(BAM_FSECONDARY)) != 0 )
        return 0; // secondary alignment
    if ( ((bam_r->core.flag)&(BAM_FSUPPLEMENTARY)) != 0 )
        return 0; // supplemetary alignment
    if ( (bam_r->core.qual) < (uint8_t)(objs->atac_mapq))
        return 0; // low qual

    // coordinates are begin and end of alignment, excluding clipped bases.
    int32_t bam_rid;
    int32_t bam_beg, bam_end;
    if (get_rcoord_bam(bam_r, &bam_rid, &bam_beg, &bam_end, 0) < 0)
        return err_msg(-1, 0, "proc_atac1: failed get coordinates");

    const char *chr = sam_hdr_tid2name(objs->atac_bam_hdr, bam_r->core.tid);
    if (chr == NULL)
        return err_msg(-1, 0, "proc_atac1: could not get chromosome name %d",
                       bam_r->core.tid);

    // skip read if in exclusion regions
    if (objs->exclude && 
            iregs_has_overlap(objs->exclude, chr, bam_beg, bam_end-1) > 0)
        return 0;

    const char *qname = bam_get_qname(bam_r);

    // error if no chromosome
    if (bam_rid == -1)
        return err_msg(-1, 0, "proc_atac1: bam rid is -1 for '%s'", qname);

    uint8_t *tag_ptr;

    // skip if no barcode
    const char *b_cb = get_tag(bam_r, objs->atac_bc_tag);
    if (b_cb == NULL || strlen(b_cb) < 4)
        return 0; // missing can be represented as dash

    // skip if number of hits of query is above threshold
    if (objs->max_nh > 0){
        tag_ptr = bam_aux_get(bam_r, objs->atac_nh_tag);
        if (tag_ptr == NULL){
            return err_msg(-1, 0, "proc_atac1: NH set but could not find tag '%s' in "
                    "read '%s'\n", objs->atac_nh_tag, qname);
        }
        int nh = (int)bam_aux2i(tag_ptr);
        if (nh > objs->max_nh)
            return 0;
    }

    // Set short query name of read
    qshort qs = qname2qshort(qname);

    // init atac read
    atac_read1_t *atac_read = atac_read_init();
    if (atac_read == NULL)
        return err_msg(-1, 0, "proc_atac1: failed to init atac read");

    atac_read->loc.rid = bam_rid;
    atac_read->loc.start = bam_beg;
    atac_read->loc.end = bam_end;
    atac_read->loc.strand = '.'; // unstranded data

    // Get overlapping variants, base calls, and base qualities.
    if ( bam1_seq_base(objs->atac_bam_hdr, bam_r, objs->gv, &atac_read->bl) < 0) {
        atac_read_dstry(atac_read);
        return err_msg(-1, 0, "proc_atac1: could not overlap variants and bases in bam record");
    }

    // add read to bam_data
    if (bam_data_atac_add_read(bam_data, b_cb, atac_read, qs) < 0) {
        atac_read_dstry(atac_read);
        return err_msg(-1, 0, "proc_atac1: failed add read");
    }

    atac_read_dstry(atac_read); // destroy since atac read is copied into bam_data

    return(0);
}

int run_atac(obj_pars *objs, bam_data_t *bam_data){
    if (objs == NULL)
        return err_msg(-1, 0, "run_atac: 'objs' is null");
    if (bam_data == NULL)
        return err_msg(-1, 0, "run_atac: 'bam_data' is null");
    if (objs->atac_bam == NULL)
        return err_msg(-1, 0, "run_atac: ATAC BAM file not given");
    if (objs->gv == NULL)
        return err_msg(-1, 0, "run_atac: VCF file not given");

    hts_itr_t *bam_itr = NULL;
    bam1_t *bam_r = NULL;

    uint64_t n_reads = 0;
    
    /* Set BAM iterator to region */
    bam_itr = sam_itr_querys(objs->atac_bam_idx, objs->atac_bam_hdr, objs->region);
    if (bam_itr == NULL){
        return err_msg(-1, 0, "run_atac: could not set region '%s' in BAM file", objs->region);
    }
    /* */

    /* Loop over BAM records */
    if (objs->verbose) log_msg("processing ATAC alignments");
    int iter_ret;
    bam_r = bam_init1();
    if (bam_r == NULL) {
        hts_itr_destroy(bam_itr);
        return err_msg(-1, 0, "run_atac: failed to init bam record");
    }
    while ( (iter_ret = sam_itr_next(objs->atac_bam, bam_itr, bam_r)) >= -1){
        if ( (objs->verbose) && (n_reads > 0) && (n_reads % (uint64_t)10e6 == 0)){
            if (bam1_unmapped(bam_r)){
                log_msg("processed %" PRIu64 " alignments (unmapped)", n_reads/(uint64_t)1e6);
            } else {
                const char *chr = sam_hdr_tid2name(objs->atac_bam_hdr, bam_r->core.tid);
                print_status_plp("million alignments", n_reads/(uint64_t)1e6, chr, 
                        (int)(bam_r->core.pos+1));
            }
        }
        n_reads++;

        if (proc_atac1(bam_r, objs, bam_data) < 0) {
            hts_itr_destroy(bam_itr);
            bam_destroy1(bam_r);
            return err_msg(-1, 0, "run_atac: failed proc_atac1");
        }
        // no more data (should this be before proc_atac1?)
        if (iter_ret == -1)
            break;
    }

    hts_itr_destroy(bam_itr);
    bam_destroy1(bam_r);

    // free dangling pairs
    bam_data_atac_free_pairs(bam_data);

    if (iter_ret < -1)
        return err_msg(-1, 0, "run_atac: failed sam_iter_next, check bam file and index");

    if (objs->verbose)
        log_msg("processed %" PRIu64 " alignments", n_reads);

    if (objs->verbose)
        log_msg("deduplicating ATAC reads");
    if (bam_data_atac_dedup(bam_data) < 0)
        return err_msg(-1, 0, "run_atac: failed to deduplicate ATAC reads");

    // go from base pair to ref/alt/other allele at variants
    if (objs->verbose)
        log_msg("running ATAC pileup at SNP sites");
    if (bam_data_atac_var_call(bam_data, objs->gv, objs->cmap, objs->min_phred) < 0)
        return err_msg(-1, 0, "run_atac: failed to call variants");

    // call peaks at fragments
    if (objs->pks){
        if (objs->verbose)
            log_msg("overlapping ATAC fragment peaks");
        if (bam_data_atac_peak_call(bam_data, objs->pks, objs->cmap) < 0)
            return err_msg(-1, 0, "run_atac: failed to call peaks");
    }

    return 0;
}

int proc_rna1(bam1_t *bam_r, obj_pars *objs, bam_data_t *bam_data){
    if (objs == NULL)
        return err_msg(-1, 0, "proc_rna1: 'objs' is null");
    if (bam_r == NULL)
        return err_msg(-1, 0, "proc_rna1: 'bam_r' is null");
    if (bam_data == NULL)
        return err_msg(-1, 0, "proc_rna1: 'bam_data' is null");

    if (bam1_unmapped(bam_r))
        return 0; // unmapped
    if ( ((bam_r->core.flag)&(BAM_FSECONDARY)) != 0 )
        return 0; // secondary alignment
    if ( ((bam_r->core.flag)&(BAM_FSUPPLEMENTARY)) != 0 )
        return 0; // supplemetary alignment
    if ( (bam_r->core.qual) < (uint8_t)(objs->rna_mapq))
        return 0; // low qual

    // coordinates are begin and end of alignment, excluding clipped bases.
    int32_t bam_rid;
    int32_t bam_beg, bam_end;
    if (get_rcoord_bam(bam_r, &bam_rid, &bam_beg, &bam_end, 0) < 0)
        return err_msg(-1, 0, "proc_rna1: failed get coordinates");

    // print_bam1_t(bam_r);
    const char *qname = bam_get_qname(bam_r);

    // error if no chromosome
    if (bam_rid == -1)
        return err_msg(-1, 0, "proc_rna1: bam rid is -1 for '%s'", qname);

    // get umi, skip if none
    const char *b_umi = get_tag(bam_r, objs->rna_umi_tag);
    if (b_umi == NULL || strlen(b_umi) < 4)
        return 0;

    // get barcode, skip if none
    const char *b_cb = get_tag(bam_r, objs->rna_bc_tag);
    if (b_cb == NULL || strlen(b_cb) < 4)
        return 0; // missing can be represented as dash

    const char *chr = sam_hdr_tid2name(objs->rna_bam_hdr, bam_r->core.tid);
    if (chr == NULL)
        return err_msg(-1, 0, "proc_rna1: could not get chromosome name %d",
                       bam_r->core.tid);

    // skip read if in exclusion regions
    /*
    const char *chr = sam_hdr_tid2name(objs->rna_bam_hdr, bam_r->core.tid);
    if (objs->exclude && 
            iregs_has_overlap(objs->exclude, chr, bam_beg, bam_end-1) > 0)
    */

    uint8_t *tag_ptr;

    // check number of hits of query, filter if above threshold
    if (objs->max_nh > 0){
        tag_ptr = bam_aux_get(bam_r, objs->rna_nh_tag);
        if (tag_ptr == NULL){
            return err_msg(-1, 0, "proc_rna1: NH set but could not find "
                    "tag '%s' in read '%s'\n", objs->rna_nh_tag, qname);
        }
        int nh = (int)bam_aux2i(tag_ptr);
        if (nh > objs->max_nh)
            return 0;
    }

    int fret, vret;

    umishort umih = umi2umishort(b_umi);

    // create rna_read_t object.
    rna_read1_t *rna_read = rna_read1_alloc();
    if (rna_read == NULL)
        return err_msg(-1, 0, "proc_rna1: failed to allocate rna_read");

    rna_read->loc.rid = bam_rid;
    rna_read->loc.start = bam_beg;
    rna_read->loc.end = bam_end;
    rna_read->loc.strand = bam_is_rev(bam_r) ? '-' : '+'; // rna is stranded
    rna_read->umi = umih;
    rna_read->n_reads = 1;

    // Get overlapping features
    if ( (fret = bam1_feat_overlap(objs->rna_bam_hdr, bam_r, objs->anno, &rna_read->gl)) < 0) {
        rna_read1_dstry(rna_read);
        return err_msg(-1, 0, "proc_rna1: could not overlap features in bam record");
    }

    // Get overlapping variants, base calls, and base qualities.
    if ( (vret = bam1_seq_base(objs->rna_bam_hdr, bam_r, objs->gv, &rna_read->bl)) < 0) {
        rna_read1_dstry(rna_read);
        return err_msg(-1, 0, "proc_rna1: could not overlap variants and bases in bam record");
    }

    // skip if no overlapping features or variants
    if (vret == 0 && fret == 0){
        rna_read1_dstry(rna_read);
        return 0;
    }

    // add read to rna
    if (bam_data_rna_add_read(bam_data, b_cb, rna_read, umih) < 0) {
        rna_read1_dstry(rna_read);
        return err_msg(-1, 0, "proc_rna1: failed add read");
    }

    rna_read1_dstry(rna_read);

    return 0;
}

int run_rna(obj_pars *objs, bam_data_t *bam_data){
    if (objs == NULL || bam_data == NULL)
        return err_msg(-1, 0, "run_rna:'objs' or 'bam_data' is null");
    hts_itr_t *bam_itr = NULL;
    bam1_t *bam_r = NULL;

    if (objs->anno == NULL)
        return err_msg(-1, 0, "run_rna: GTF file not given");
    if (objs->rna_bam == NULL)
        return err_msg(-1, 0, "run_rna: RNA BAM file not given");
    if (objs->gv == NULL)
        return err_msg(-1, 0, "run_rna: VCF file not given");

    uint64_t n_reads = 0;
    
    /* Set BAM iterator to region */
    bam_itr = sam_itr_querys(objs->rna_bam_idx, objs->rna_bam_hdr, objs->region);
    if (bam_itr == NULL)
        return err_msg(-1, 0, "run_rna: could not set region '%s' in BAM file", objs->region);

    /* Loop over BAM records */
    if (objs->verbose)
        log_msg("processing RNA alignments");
    int iter_ret;
    bam_r = bam_init1();
    if (bam_r == NULL) {
        hts_itr_destroy(bam_itr);
        return err_msg(-1, 0, "run_rna: failed to init bam record");
    }
    while ( (iter_ret = sam_itr_next(objs->rna_bam, bam_itr, bam_r)) >= -1){
        if ( (objs->verbose) && (n_reads > 0) && (n_reads % (uint64_t)10e6 == 0)){
            if (bam1_unmapped(bam_r)){
                log_msg("processed %" PRIu64 " alignments (unmapped)", n_reads/(uint64_t)1e6);
            } else {
                const char *chr = sam_hdr_tid2name(objs->rna_bam_hdr, bam_r->core.tid);
                print_status_plp("million alignments", n_reads/(uint64_t)1e6, chr, 
                        (int)(bam_r->core.pos+1));
            }
        }
        n_reads++;

        if (proc_rna1(bam_r, objs, bam_data) < 0){
            hts_itr_destroy(bam_itr);
            bam_destroy1(bam_r);
            return err_msg(-1, 0, "run_rna: failed proc_rna1");
        }
        if (iter_ret == -1)
            break;
    }

    hts_itr_destroy(bam_itr);
    bam_destroy1(bam_r);

    if (iter_ret < -1)
        return err_msg(-1, 0, "run_rna: failed sam_iter_next");

    if (objs->verbose)
        log_msg("processed %" PRIu64 " alignments", n_reads);

    if (objs->verbose)
        log_msg("deduplicating RNA reads");
    if (bam_data_rna_dedup(bam_data) < 0)
        return err_msg(-1, 0, "run_rna: failed to deduplicate RNA reads");

    if (objs->verbose)
        log_msg("running RNA pileup at SNP sites");
    if (bam_data_rna_var_call(bam_data, objs->gv, objs->cmap, objs->min_phred) < 0)
        return err_msg(-1, 0, "run_rna: failed to call variants");

    return 0;
}

int bam_count(bam_data_t *bam_dat, obj_pars *objs, char *filename){
    if (bam_dat == NULL || objs == NULL || filename == NULL)
        return err_msg(-1, 0, "bam_count: arguments are NULL");

    bam_counts_t *agc = bam_counts_init();
    if (objs->gv != NULL && (bam_counts_add_gv(agc, objs->gv) < 0)) {
        bam_counts_dstry(agc);
        return err_msg(-1, 0, "bam_count: failed to add GV");
    }

    if (objs->anno != NULL && (bam_counts_add_gene_map(agc, objs->anno->gene_ix) < 0)) {
        bam_counts_dstry(agc);
        return err_msg(-1, 0, "bam_count: failed to add gene map");
    }

    if (objs->pks != NULL && (bam_counts_add_peaks(agc, objs->pks) < 0)) {
        bam_counts_dstry(agc);
        return err_msg(-1, 0, "bam_count: failed to add peaks");
    }

    if (bam_counts_add_bc_map(agc, bam_dat->bcs) < 0) {
        bam_counts_dstry(agc);
        return err_msg(-1, 0, "bam_count: failed to add barcode map");
    }

    if (objs->verbose)
        log_msg("generating feature counts");
    if (bam_counts_count(agc, bam_dat) < 0) {
        bam_counts_dstry(agc);
        return err_msg(-1, 0, "bam_count: failed to count features");
    }

    if (objs->verbose)
        log_msg("writing feature counts");
    if (bam_counts_write(agc, objs->anno, objs->gv, filename, 0) < 0) {
        bam_counts_dstry(agc);
        return err_msg(-1, 0, "bam_count: failed to write counts");
    }

    bam_counts_dstry(agc);
    return 0;
}

void *atac_thrd_fx(void *arg){
    bam_thrd_args *t = (bam_thrd_args *)arg;

    uint32_t i = 0;
    for (i = 0; i < t->n; ++i){
        if (proc_atac1(t->a[i], t->objs, t->bam_data) < 0){
            t->ret = -1;
            return(NULL);
        }
        if (t->a[i]) bam_destroy1(t->a[i]);
        t->a[i] = NULL;
    }
    return(NULL);
}

void *rna_thrd_fx(void *arg){
    bam_thrd_args *t = (bam_thrd_args *)arg;

    uint32_t i = 0;
    for (i = 0; i < t->n; ++i){
        if (proc_rna1(t->a[i], t->objs, t->bam_data) < 0){
            t->ret = -1;
            return(NULL);
        }
        if (t->a[i]) bam_destroy1(t->a[i]);
        t->a[i] = NULL;
    }
    return(NULL);
}

int rna_atac_thrd_call(bam1_t **a, uint32_t n, bam_data_t *bam_data, obj_pars *objs, 
        int rna_atac){
    if (objs == NULL)
        return err_msg(-1, 0, "rna_atac_thrd_call: 'objs' is null");
    uint32_t i;
    if (objs->threads < 1)
        return err_msg(-1, 0, "rna_atac_thrd_call: threads=%u is less than 1", objs->threads);

    if (rna_atac != 0 && rna_atac != 1)
        return err_msg(-1, 0, "rna_atac_thrd_call: rna_atac must be 0 or 1");

    uint32_t mt1 = objs->threads - 1;
    double step = (double)n / (double)objs->threads;
    step = ceil(step);
    uint32_t stepi = (uint32_t)ROUND_2_INT(step);
    uint32_t step_last = n - ((mt1) * stepi);

    bam_thrd_args *targs = malloc(objs->threads * sizeof(bam_thrd_args));
    if (targs == NULL)
        return err_msg(-1, 0, "rna_atac_thrd_call: could not malloc thread arguments");
    for (i = 0; i < objs->threads; ++i){
        targs[i].a = a + (stepi * i);
        targs[i].n = i == (mt1) ? step_last : stepi;
        targs[i].bam_data = bam_data;
        targs[i].objs = objs;
        targs[i].ret = 0;
        if (stepi * i + targs[i].n > n) {
            free(targs);
            return err_msg(-1, 0, "rna_atac_thrd_call: improper step calculated");
        }
    }

    pthread_t *ids = malloc(objs->threads * sizeof(pthread_t));
    if (ids == NULL) {
        free(targs);
        return err_msg(-1, 0, "rna_atac_thrd_call: could not malloc thread ids");
    }

    int err;
    for (i = 0; i < objs->threads; ++i){
        if (rna_atac == 0){
            err = pthread_create(ids + i, NULL, rna_thrd_fx, &targs[i]);
        } else if (rna_atac == 1){
            err = pthread_create(ids + i, NULL, atac_thrd_fx, &targs[i]);
        } else {
            free(targs);
            free(ids);
            return err_msg(-1, 0, "rna_atac_thrd_call: rna_atac must be 0 or 1");
        }

        if (err != 0) {
            uint32_t j;
            for (j = 0; j < i; ++j)
                pthread_join(ids[j], NULL);

            free(targs);
            free(ids);
            return err_msg(-1, 0, "rna_atac_thrd_call: could not create thread '%u': '%d'", i, err);
        }
    }

    for (i = 0; i < objs->threads; ++i){
        err = pthread_join(ids[i], NULL);
        if (err != 0) {
            free(targs);
            free(ids);
            return err_msg(-1, 0, "rna_atac_thrd_call: could not join thread '%u': '%d'", i, err);
        }
        if (targs[i].ret < 0){
            err_msg(-1, 0, "rna_atac_thrd_call: thread '%u' returned error '%d'", i, targs[i].ret);
            free(targs);
            free(ids);
            return -1;
        }
    }

    free(targs);
    free(ids);

    return 0;
}

