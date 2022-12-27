
#include <stdio.h>
#include "clopts.h"
#include "atac_data.h"
#include "rna_data.h"
#include "bam_dat.h"
#include "str_util.h"
#include "overlap.h"
#include "r_count.h"


int run_atac(obj_pars *objs, bam_atac_t **bam_atac){
    bam_atac_t *bama = NULL;
    hts_itr_t *bam_itr = NULL;
    bam1_t *bam_r = NULL;

    if (objs->atac_bam == NULL)
        return err_msg(-1, 0, "run_atac: ATAC BAM file not given");
    if (objs->gv == NULL)
        return err_msg(-1, 0, "run_atac: VCF file not given");
    bama = bam_atac_init();
    uint64_t n_reads = 0;
    
    /* Set BAM iterator to region */
    bam_itr = sam_itr_querys(objs->atac_bam_idx, objs->atac_bam_hdr, objs->region);
    if (bam_itr == NULL){
        return err_msg(-1, 0, "run_atac: could not set region %s in BAM file", objs->region);
    }
    /* */

    /* Loop over BAM records */
    if (objs->verbose) log_msg("processing ATAC alignments");
    int iter_ret;
    bam_r = bam_init1();
    while ( (iter_ret = sam_itr_next(objs->atac_bam, bam_itr, bam_r)) >= 0){
        const char *chr = sam_hdr_tid2name(objs->atac_bam_hdr, bam_r->core.tid);
        if ( (objs->verbose) && (n_reads > 0) && (n_reads % (uint64_t)10e6 == 0)){
            if (bam1_unmapped(bam_r)){
                log_msg("processed %" PRIu64 " alignments (unmapped)", n_reads/(uint64_t)1e6);
            } else {
                print_status_plp("million alignments", n_reads/(uint64_t)1e6, chr, 
                        (int)(bam_r->core.pos+1));
            }
        }
        n_reads++;

        if (bam1_unmapped(bam_r)) continue; // unmapped
        if ( ((bam_r->core.flag)&(BAM_FSECONDARY)) != 0 ) continue; // secondary alignment
        if ( ((bam_r->core.flag)&(BAM_FSUPPLEMENTARY)) != 0 ) continue; // supplemetary alignment
        if ( (bam_r->core.qual) < (uint8_t)(objs->atac_mapq)) continue; // low qual
        
        // coordinates are begin and end of alignment, excluding clipped bases.
        int32_t bam_rid;
        hts_pos_t bam_beg, bam_end;
        if (get_rcoord_bam(bam_r, &bam_rid, &bam_beg, &bam_end, 0) < 0)
            return err_msg(-1, 0, "run_atac: failed get coordinates");

        // skip read if in exclusion regions
        if (objs->exclude && 
            iregs_has_overlap(objs->exclude, chr, bam_beg, bam_end-1) > 0)
            continue;

        // print_bam1_t(bam_r);
        const char *qname = bam_get_qname(bam_r);

        uint8_t *tag_ptr;

        const char *b_cb = get_tag(bam_r, objs->atac_bc_tag);
        if (b_cb == NULL || strlen(b_cb) < 4) continue; // missing can be represented as dash

        if (objs->max_nh > 0){
            tag_ptr = bam_aux_get(bam_r, objs->rna_nh_tag);
            if (tag_ptr == NULL){
                return err_msg(-1, 0, "NH set but could not find tag %s in "
                        "read %s\n", objs->rna_nh_tag, qname);
            }
            int nh = (int)bam_aux2i(tag_ptr);
            if (nh > objs->max_nh) continue;
        }

        // Set region and short query name of read
        qshort qs = qname2qshort(qname);

        // init atac read
        atac_read1_t *ar = atac_read_init();
        if (ar == NULL) return(-1);

        // coordinates are begin and end of alignment, excluding clipped bases.
        ar->loc.rid = bam_rid;
        ar->loc.start = bam_beg;
        ar->loc.end = bam_end;
        ar->loc.strand = '.'; // unstranded data

        // Get overlapping variants, base calls, and base qualities.
        int vret;
        seq_blist_t *bl = NULL;
        if ( (vret = bam1_seq_base(objs->atac_bam_hdr, bam_r, objs->gv, &bl)) < 0)
            return(-1);

        // create seq_base_t
        seq_base_t *base = bl->head;
        while (base){
            if (atac_read1_add_base(ar, base) < 0) return(-1);
            base = base->next;
        }
        free_seq_blist(bl);
        free(bl);

        // add read to bama
        if (bam_atac_add_read(bama, b_cb, ar, qs) < 0)
            return err_msg(-1, 0, "run_atac: failed add read");

        atac_read_dstry(ar); // destroy since atac read is copied into bama
    }
    if (objs->verbose) log_msg("processed %" PRIu64 " alignments", n_reads);

    hts_itr_destroy(bam_itr);
    bam_destroy1(bam_r);

    bam_atac_free_pairs(bama); // free the rest of the reads that weren't paired to dups.

    if (objs->verbose) log_msg("deduplicating reads");
    if (bam_atac_dedup(bama) < 0) return(-1);

    bam_atac_free_dups(bama);

    // go from base pair to ref/alt/other allele at variants
    if (objs->verbose) log_msg("calling variants");
    if (bam_atac_var_call(bama, objs->gv, objs->cmap, objs->min_phred) < 0) return(-1);

    // call peaks at fragments
    if (objs->pks){
        if (objs->verbose) log_msg("calling peaks");
        if (bam_atac_peak_call(bama, objs->pks, objs->cmap) < 0) return(-1);
    }

    *bam_atac = bama;

    return(0);
}

int run_rna(obj_pars *objs, bam_rna_t **bam_rna){
    bam_rna_t *br = NULL;
    hts_itr_t *bam_itr = NULL;
    bam1_t *bam_r = NULL;

    if (objs->anno == NULL)
        return err_msg(-1, 0, "run_rna: GTF file not given");
    if (objs->rna_bam == NULL)
        return err_msg(-1, 0, "run_rna: RNA BAM file not given");
    if (objs->gv == NULL)
        return err_msg(-1, 0, "run_rna: VCF file not given");

    br = bam_rna_alloc();
    uint64_t n_reads = 0;
    
    /*
    uint8_t n_feat = 0, m_feat = 4;
    char **feat = (char **)calloc(m_feat, sizeof(char *));
    uint8_t *splice = (uint8_t *)calloc(m_feat, sizeof(uint8_t));;
    if (feat == NULL || splice == NULL)
        return err_msg(-1, 0, "run_rna: %s", strerror(errno));
        */

    /* Set BAM iterator to region */
    bam_itr = sam_itr_querys(objs->rna_bam_idx, objs->rna_bam_hdr, objs->region);
    if (bam_itr == NULL){
        return err_msg(-1, 0, "run_rna: could not set region %s in BAM file", objs->region);
    }
    /* */

    /* Loop over BAM records */
    if (objs->verbose) log_msg("processing RNA alignments");
    int iter_ret;
    bam_r = bam_init1();
    while ( (iter_ret = sam_itr_next(objs->rna_bam, bam_itr, bam_r)) >= 0){
        const char *chr = sam_hdr_tid2name(objs->rna_bam_hdr, bam_r->core.tid);
        if ( (objs->verbose) && (n_reads > 0) && (n_reads % (uint64_t)10e6 == 0)){
            if (bam1_unmapped(bam_r)){
                log_msg("processed %" PRIu64 " alignments (unmapped)", n_reads/(uint64_t)1e6);
            } else {
                print_status_plp("million alignments", n_reads/(uint64_t)1e6, chr, 
                        (int)(bam_r->core.pos+1));
            }
        }
        n_reads++;

        if (bam1_unmapped(bam_r)) continue; // unmapped
        if ( ((bam_r->core.flag)&(BAM_FSECONDARY)) != 0 ) continue; // secondary alignment
        if ( ((bam_r->core.flag)&(BAM_FSUPPLEMENTARY)) != 0 ) continue; // supplemetary alignment
        if ( (bam_r->core.qual) < (uint8_t)(objs->atac_mapq)) continue; // low qual

        // coordinates are begin and end of alignment, excluding clipped bases.
        int32_t bam_rid;
        hts_pos_t bam_beg, bam_end;
        if (get_rcoord_bam(bam_r, &bam_rid, &bam_beg, &bam_end, 0) < 0)
            return err_msg(-1, 0, "run_atac: failed get coordinates");

        // skip read if in exclusion regions
        if (objs->exclude && 
            iregs_has_overlap(objs->exclude, chr, bam_beg, bam_end-1) > 0)
            continue;

        // print_bam1_t(bam_r);
        const char *qname = bam_get_qname(bam_r);

        uint8_t *tag_ptr;

        const char *b_umi = get_tag(bam_r, objs->rna_umi_tag);
        if (b_umi == NULL) continue;
        if (strlen(b_umi) < 4) continue; // missing can be represented as dash

        const char *b_cb = get_tag(bam_r, objs->rna_bc_tag);
        if (b_cb == NULL || strlen(b_cb) < 4) continue; // missing can be represented as dash

        if (objs->max_nh > 0){
            tag_ptr = bam_aux_get(bam_r, objs->rna_nh_tag);
            if (tag_ptr == NULL){
                return err_msg(-1, 0, "run_rna: NH set but could not find "
                        "tag %s in read %s\n", objs->rna_nh_tag, qname);
            }
            int nh = (int)bam_aux2i(tag_ptr);
            if (nh > objs->max_nh) continue;
        }

        int fret, vret;

        // Get overlapping features
        seq_glist_t *gl = NULL;
        if ( (fret = bam1_feat_overlap(objs->rna_bam_hdr, bam_r, objs->anno, &gl)) < 0)
            return(-1);
        if (fret < 0)
            return err_msg(-1, 0, "run_rna: could not overlap features in bam record");
        
        // Get overlapping variants, base calls, and base qualities.
        seq_blist_t *bl = NULL;
        if ( (vret = bam1_seq_base(objs->rna_bam_hdr, bam_r, objs->gv, &bl)) < 0)
            return(-1);

        // skip if no overlapping features or variants
        if (vret == 0 && fret == 0){
            seq_glist_free(gl);
            free(gl);
            free_seq_blist(bl);
            free(bl);
            continue;
        }

        // create rna_read_t object.
        rna_read1_t *rna_read = rna_read1_alloc();
        if (rna_read == NULL) return(-1);

        rna_read->loc.rid = bam_rid;
        rna_read->loc.start = bam_beg;
        rna_read->loc.end = bam_end;
        rna_read->loc.strand = bam_is_rev(bam_r) ? '-' : '+'; // rna is stranded


        // create seq_gene_t
        seq_gene_t *g = gl->head;
        while (g){
            if (rna_read1_add_gene(rna_read, g) < 0) return(-1);
            g = g->next;
        }
        seq_glist_free(gl);
        free(gl);
        
        // create seq_base_t
        seq_base_t *base = bl->head;
        while (base){
            if (rna_read1_add_base(rna_read, base) < 0) return(-1);
            base = base->next;
        }
        free_seq_blist(bl);
        free(bl);

        // add read to rna
        if (bam_rna_add_read(br, b_cb, rna_read, b_umi) < 0)
            return(-1);

        rna_read1_dstry(rna_read);
    }
    if (objs->verbose) log_msg("processed %" PRIu64 " alignments", n_reads);

    hts_itr_destroy(bam_itr);
    bam_destroy1(bam_r);
    /*
    free(feat);
    free(splice);
    */

    *bam_rna = br;

    if (objs->verbose) log_msg("deduplicating reads");
    if (bam_rna_dedup(br) < 0)
        return(-1);

    bam_rna_free_dups(br);

    if (objs->verbose) log_msg("calling variants");
    if (bam_rna_var_call(br, objs->gv, objs->cmap, objs->min_phred) < 0)
        return(-1);

    return(0);
}

int bam_count(bam_data_t *bam_dat, obj_pars *objs, char *filename){
    if (bam_dat == NULL || objs == NULL)
        return err_msg(-1, 0, "bam_count: arguments are NULL");

    bam_ag_t *agc = bam_ag_init();
    if (objs->gv != NULL && (bam_ag_add_gv(agc, objs->gv) < 0))
        return -1;

    if (objs->anno != NULL && (bam_ag_add_gene_map(agc, objs->anno->gene_ix) < 0))
        return -1;

    if (objs->pks != NULL && (bam_ag_add_peaks(agc, objs->pks) < 0))
        return -1;

    if (bam_ag_add_bc_map(agc, bam_dat->bcs) < 0)
        return -1;

    if (objs->verbose) log_msg("generating gene counts from RNA");
    if (bam_rna_gc_count(agc, bam_dat) < 0)
        return -1;

    if (objs->pks){
        if (objs->verbose) log_msg("generating peak counts from ATAC");
        if (bam_atac_pc_count(agc, bam_dat) < 0)
            return -1;
    }

    if (objs->verbose) log_msg("generating allele counts from ATAC");
    if (bam_atac_ac_count(agc, bam_dat) < 0)
        return -1;

    if (objs->verbose) log_msg("generating allele counts from RNA");
    if (bam_rna_ac_count(agc, bam_dat) < 0)
        return -1;

    if (objs->verbose) log_msg("writing counts");
    if (bam_ag_write(agc, objs->anno, objs->gv, filename) < 0)
        return -1;

    bam_ag_dstry(agc);
    return 0;
}

