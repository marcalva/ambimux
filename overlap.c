
#include "overlap.h"
#include "sam_read.h"
#include "str_util.h"
#include "bins.h"
#include "htslib/khash.h"
#include "kbtree.h"
#include <string.h>

int bam1_feat_overlap(const sam_hdr_t *h, bam1_t *b, const gene_anno_t *a, 
        seq_glist_t **gl){
    if (h == NULL || b == NULL || a == NULL)
        return err_msg(-1, 0, "bam1_feat_overlap: argument is null");

    int32_t tid = b->core.tid;
    const char *ref = sam_hdr_tid2name(h, (int)tid);
    
    hts_pos_t b_beg = b->core.pos;
    hts_pos_t b_end = bam_endpos(b);
    if (b_end == b_beg + 1) return(-1); // unmapped;

    char strand = '+';
    if (bam_is_rev(b)) strand = '-';

    // get features that contain bam1_t read b.
    int genes_len = 1;
    gene_t **genes = malloc(genes_len * sizeof(gene_t *)); // must be freed.

    if (genes == NULL)
        return err_msg(-1, 0, "bam1_feat_overlap: %s", strerror(errno));

    int ngenes = feats_from_region(a, ref, b_beg, b_end, 1, strand, &genes, &genes_len);
    if (ngenes < 0) return -1;

    if (*gl == NULL){
        *gl = seq_glist_alloc();
        if (*gl == NULL) return(-1);
    }

    // add gene ID and splice status
    seq_gene_t *sg = seq_gene_alloc();
    if (sg == NULL) return(-1);
    int i;
    for (i = 0; i < ngenes; ++i){
        // check for bugs
        if (genes[i] == NULL)
            return err_msg(-1, 0, "bam1_feat_overlap: gene is null");

        char *id = genes[i]->id;
        if (id == NULL)
            return err_msg(-1, 0, "bam1_feat_overlap: gene name not found");
        int32_t gid = str_map_ix(a->gene_ix, id);
        if (gid < 0) return(-1);
        sg->gene_id = gid;

        int sret = 0;
        sg->splice = bam1_spliced(b, genes[i], &sret);
        if (sret < 0)
            return err_msg(-1, 0, "bam1_feat_overlap: splice failed");
        
        if (seq_glist_add_gene(*gl, sg) < 0) return(-1);
    }
    seq_gene_dstry(sg); // sg is copied into gl
    free(genes);

    int ret = (int)((*gl)->n);
    return ret;
}

/* 
 * Detect whether a bam record is spliced. 
 * Assumes the bam alignment and gene lie in the same chromosome.
 */
uint8_t bam1_spliced(bam1_t *b, gene_t *g, int *ret){
    *ret = 0;

    if (g == NULL){
        *ret = -1;
        return 0;
    }

    int n_iso = g->isoforms_n;

    if (n_iso <= 0 || g->bt_isoforms == NULL){
        err_msg(-1, 0, "bam1_spliced: %s has no isoforms", g->id);
        *ret = -1;
        return 0;
    }

    // overlap of read segment that consumes query and reference (CQR segment)
    double *iro = calloc(n_iso, sizeof(double)); // length of CQR segment that overlaps intron
    double *ero = calloc(n_iso, sizeof(double)); // length of CQR segment that overlaps exon
    double *rl = calloc(n_iso, sizeof(double)); // length of CQR segment that overlaps the isoform
    double *pi = calloc(n_iso, sizeof(double)); // percent of CQR segment that overlaps intron
    double *pe = calloc(n_iso, sizeof(double)); // percent of CQR segment that overlaps exon
    double *pt = calloc(n_iso, sizeof(double)); // percent of CQR segment that overlaps the isoform
    uint8_t *sj = calloc(n_iso, sizeof(uint8_t)); // splice junctions. Set to 1 if alignment overlaps splice juction. 0 otherwise

    if (iro == NULL || ero == NULL || rl == NULL || pe == NULL || sj == NULL){
        err_msg(-1, 0, "bam1_spliced: %s", strerror(errno));
        *ret = -1;
        return 0;
    }

    int i;
    for (i = 0; i < n_iso; ++i){
        iro[i] = 0;
        ero[i] = 0;
        rl[i] = 0;
        pe[i] = 0;
        pi[i] = 0;
        pt[i] = 0;
        sj[i] = 0;
    }

    char strand = '+';
    if (bam_is_rev(b)) strand = '-';

    uint32_t *cigar_raw = bam_get_cigar(b);
    uint32_t n_cigar = b->core.n_cigar;

    // pos is 0-based leftmost base of first CIGAR op that consumes reference.
    hts_pos_t pos = b->core.pos;

    // get length of cigar length that consumes reference and query
    uint32_t qr_len = 0;
    for (i = 0; i < n_cigar; ++i){
        uint32_t c_op = bam_cigar_op(cigar_raw[i]); // lower 4 bits is cigar op
        uint32_t c_len = bam_cigar_oplen(cigar_raw[i]); // higher 24 bits is length

        int cr = bam_cigar_type(c_op)&2; // consumes reference
        int cq = bam_cigar_type(c_op)&1; // consumes query

        if (cr && cq) qr_len += c_len;
    }

    int iso_i = 0; // isoform index
    kbtree_t(kb_iso) *bt = g->bt_isoforms;
    kbitr_t itr;
    kb_itr_first(kb_iso, bt, &itr);
    for (; kb_itr_valid(&itr); kb_itr_next(kb_iso, bt, &itr)){
        isoform_t *iso = &kb_itr_key(isoform_t, &itr);

        if (iso == NULL){
            err_msg(-1, 0, "bam1_spliced: isoform not stored properly");
            *ret = -1;
            return 0;
        }

        exon_t *exon = iso->exons;
        if (exon == NULL){
            iso_i++;
            continue;
        }

        // position of CIGAR segment in reference sequence
        hts_pos_t r_beg = pos, r_end = pos;

        for (i = 0; i < n_cigar; ++i){ // for each CIGAR op
            uint32_t c_op = bam_cigar_op(cigar_raw[i]); // lower 4 bits is cigar op
            uint32_t c_len = bam_cigar_oplen(cigar_raw[i]); // higher 24 bits is length

            int cr = bam_cigar_type(c_op)&2; // consumes reference
            int cq = bam_cigar_type(c_op)&1; // consumes query

            if (cr == 0 && cq == 0)
                continue;

            // set begin ref position of CIGAR seg to previous seg ref end position
            r_beg = r_end;

            // if consumes reference, add cigar length to get end position 
            // of CIGAR seg in reference seq.
            // if doesn't, e.g. insertion, r_end is still equal to r_beg
            if (cr)
                r_end = r_beg + c_len;

            /*********************************************************
             * check splice junction N
             *********************************************************/

            if (c_op == BAM_CREF_SKIP){
                exon = iso->exons;
                while (exon != NULL && exon->next != NULL && exon->end < r_beg){
                    exon = exon->next;
                }
                if (exon != NULL && exon->next != NULL && 
                    exon->end == r_beg && exon->next->beg == r_end){
                    sj[iso_i] = 1;
                }
            }

            /*********************************************************
             * CQR length
             * overlap with isoform [beg,end)
             *********************************************************/

            if (cr && cq){
                rl[iso_i] += (double)bp_overlap(r_beg, r_end, strand, 
                        iso->beg, iso->end, g->strand);
            }

            /*********************************************************
             * check intron overlap
             *********************************************************/

            if (cr && cq){
                int64_t ovrlp = 0;
                exon = iso->exons->next; // first exon is dummy node
                while (exon != NULL && exon->next != NULL){
                    int64_t tmp_ovrlp = bp_overlap(r_beg, r_end, strand, 
                            exon->end, exon->next->beg, g->strand);
                    if (tmp_ovrlp < 0){
                        *ret = -1;
                        return 0;
                    }
                    ovrlp += tmp_ovrlp;
                    exon = exon->next;
                }
                iro[iso_i] += (double)ovrlp; // add overlap of CIGAR segment
            }

            /*********************************************************
             * check exon overlap
             *********************************************************/

            if (cr && cq){
                int64_t ovrlp = 0;
                exon = iso->exons->next;
                while (exon != NULL && exon->beg < r_end){
                    int64_t tmp_ovrlp = bp_overlap(r_beg, r_end, strand, 
                            exon->beg, exon->end, g->strand);
                    if (tmp_ovrlp < 0){
                        *ret = -1;
                        return 0;
                    }
                    ovrlp += tmp_ovrlp;
                    exon = exon->next;
                }

                ero[iso_i] += (double)ovrlp; // add overlap of CIGAR segment
            }
        }

        pt[iso_i] = (double)rl[iso_i] / (double)qr_len;
        if (rl[iso_i] > 0){
            pi[iso_i] = iro[iso_i] / rl[iso_i];
            pe[iso_i] = ero[iso_i] / rl[iso_i];
        }
        iso_i++;
    }

    /*********************************************************
     * get spliced or unspliced status of record for each isoform
     *********************************************************/

    int n_o_iso = 0; // number isoforms with full overlap with UMI
    int is_spl = 0, is_unspl = 0;
    for (i = 0; i < n_iso; ++i){
        if (pt[i] < 1) continue;
        n_o_iso++;
        if (sj[i] == 1) // if overlaps splice junction
            is_spl++;
        if (pi[i] > 0) // if at least one base pair overlaps intron
            is_unspl++;
    }

    /*********************************************************
     * if read overlaps a splice junction from an isoform, then it is 
     *  considered spliced
     * if read does not overlap any splice junction and at least one 
     *  base pair overlaps an intron for all isoforms, it is unspliced
     * otherwise it is ambiguous.
     *********************************************************/

    uint8_t spl_stat;
    if (is_spl > 0) spl_stat = SPLICE;
    else if ((is_unspl > 0) && (is_unspl == n_o_iso)) spl_stat = UNSPLICE;
    else spl_stat = AMBIG;

    free(pt);
    free(pi);
    free(pe);
    free(iro);
    free(ero);
    free(rl);
    free(sj);

    return spl_stat;
}

/* new bam1 vars overlap. */
int bam1_vars_overlap(const sam_hdr_t *h, bam1_t *b, g_var_t *gv, 
        var_t **vars){
    if (b == NULL || h == NULL || gv == NULL)
        return err_msg(-1, 0, "bam1_vars_overlap: argument is null");

    uint32_t *cigar_raw;
    cigar_raw = bam_get_cigar(b);
    uint32_t n_cigar = b->core.n_cigar;

    int32_t tid = b->core.tid;
    const char *ref = sam_hdr_tid2name(h, (int)tid);
    hts_pos_t left_pos = b->core.pos; // position of first base that consumes the reference.

    int n_vars = 0;
    /* r_pos stores position in ref. sequence 0-based */
    int64_t r_pos_beg = left_pos;
    int64_t r_pos_end = left_pos;
    int i;
    for (i = 0; i < n_cigar; i++){
        
        uint32_t cigar_op = bam_cigar_op(cigar_raw[i]);
        uint32_t cigar_oplen = bam_cigar_oplen(cigar_raw[i]);

        int cr = bam_cigar_type(cigar_op)&2; // consumes ref seq.
        int cq = bam_cigar_type(cigar_op)&1; // consumes query seq.
        if (cr){
            r_pos_end += cigar_oplen;
        }
        // if cigar doesn't consume query and reference, we can't match the query to the 
        // reference and there is no overlapping base.

        // get overlapping variants from r_pos region.
        if ( cr && cq ){
            int ret = region_vars(gv, ref, r_pos_beg, r_pos_end, vars);
            if (ret < 0)
                return err_msg(-1, 0, "bam1_vars_overlap: failed region overlap at "
                        "%s:%"PRIhts_pos"-%"PRIhts_pos"\n", ref, r_pos_beg, r_pos_end);
            else n_vars += ret;
        }

        r_pos_beg = r_pos_end;
    }
    return(n_vars);
}

int bam1_seq_base(const sam_hdr_t *h, bam1_t *b, g_var_t *gv, seq_blist_t **bl){
    if (h == NULL || b == NULL || gv == NULL)
        return err_msg(-1, 0, "bam1_seq_base: argument is null");

    int vret;
    var_t *ovars = NULL;
    if ( (vret = bam1_vars_overlap(h, b, gv, &ovars)) < 0 )
        return(-1);

    if (*bl == NULL){
        *bl = alloc_seq_blist();
        if (*bl == NULL) return(-1);
    }

    seq_base_t *sb = alloc_seq_base();
    if (sb == NULL) return(-1);
    var_t *tvar;
    for (tvar = ovars; tvar != NULL; tvar = tvar->next){
        // Get position of the variant
        int32_t v_rid = tvar->b->rid;
        hts_pos_t v_pos = tvar->b->pos;
        // const char *v_rstr = bcf_hdr_id2name(objs->vcf_hdr, v_rid);

        // Get base pair at overlapping site
        uint8_t base, qual;
        if (bam1_site_base(b, b->core.tid, v_pos, &base, &qual) < 0)
            return(-1);

        // create seq_base_t object
        sb->base = base;
        sb->qual = qual;
        sb->pos.rid = v_rid;
        sb->pos.pos = v_pos;
        sb->pos.strand = '.';
        sb->next = NULL;

        // TODO: better handle multi-allelic SNPs
        if (blist_add_base(*bl, sb, 1, 1) < 0)
            return(-1);

    }
    dstry_seq_base(sb); // sb is copied to bl

    while (ovars != NULL){
        tvar = ovars->next;
        free(ovars);
        ovars = tvar;
    }

    int ret = (int)((*bl)->n);
    return ret;
}

int64_t bp_overlap(int64_t a1, int64_t a2, char a_strand, int64_t b1, int64_t b2, char b_strand){
    if (a2 < a1 || b2 < b1){
        return err_msg(-1, 0, "bp_overlap: incorrect parameters A=[%i, %i) B=[%i, %i)", 
                a1, a2, b1, b2);
    }

    if (a2 <= b1 || b2 <= a1)
        return 0;

    if ( (a_strand == '+') && (b_strand == '-') )
        return 0;
    if ( (a_strand == '-') && (b_strand == '+') )
        return 0;

    // c1 is greatest pos of a1 and b1
    // c2 is least pos
    int64_t c1 = a1 > b1 ? a1 : b1; // max of a1, b1
    int64_t c2 = a2 < b2 ? a2 : b2; // min of a2, b2
    
    int64_t o = c2 - c1;
    if (o < 0){
        return err_msg(-1, 0, "bp_overlap: overlap is negative");
    }

    return o;
}

