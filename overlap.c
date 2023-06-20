
#include "overlap.h"
#include "sam_read.h"
#include "str_util.h"
#include "bins.h"
#include "htslib/khash.h"
#include "kbtree.h"
#include <string.h>

int bam1_feat_overlap(const sam_hdr_t *h, bam1_t *b, const gene_anno_t *a, 
        ml_t(seq_gene_l) *gl){
    if (h == NULL || b == NULL || a == NULL)
        return err_msg(-1, 0, "bam1_feat_overlap: argument is null");

    if (gl == NULL)
        return err_msg(-1, 0, "bam1_feat_overlap: gl is null");
    if (ml_size(gl) != 0)
        return err_msg(-1, 0, "bam1_feat_overlap: gl is not empty");

    int32_t tid = b->core.tid;
    const char *ref = sam_hdr_tid2name(h, (int)tid);
    
    int32_t b_beg = (int32_t)b->core.pos;
    int32_t b_end = (int32_t)bam_endpos(b);
    if (b_end == b_beg + 1) return(-1); // unmapped;

    char strand = '+';
    if (bam_is_rev(b)) strand = '-';

    // get features that contain bam1_t read b.
    ml_t(gl) genes;
    ml_init(gl, &genes);

    int ngenes = feats_from_region(a, ref, b_beg, b_end, 1, strand, &genes);
    if (ngenes < 0) return -1;

    // add gene ID and splice status
    seq_gene_t sg;
    seq_gene_init(&sg);
    ml_node_t(gl) *gn;
    for (gn = ml_begin(&genes); gn; gn = ml_node_next(gn)){
        gene_t *gene = ml_node_val(gn);

        char *id = gene->id;
        assert(id != NULL);
        int32_t gid = str_map_ix(a->gene_ix, id);
        if (gid < 0)
            return err_msg(-1, 0, "bam1_feat_overlap: failed to get gene index from gtf anno");
        sg.gene_id = gid;

        int sret = 0;
        sg.splice = bam1_spliced(b, gene, &sret);
        if (sret < 0)
            return err_msg(-1, 0, "bam1_feat_overlap: splice failed");

        if (ml_insert(seq_gene_l, gl, sg, 0, 0) < 0)
            return(-1);
        
    }

    ml_free(gl, &genes);

    int ret = (int)(ml_size(gl));
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

    int ix;
    for (ix = 0; ix < n_iso; ++ix){
        iro[ix] = 0.0;
        ero[ix] = 0.0;
        rl[ix] = 0.0;
        pe[ix] = 0.0;
        pi[ix] = 0.0;
        pt[ix] = 0.0;
        sj[ix] = 0.0;
    }

    char strand = '+';
    if (bam_is_rev(b)) strand = '-';

    uint32_t *cigar_raw = bam_get_cigar(b);
    uint32_t n_cigar = b->core.n_cigar;

    // pos is 0-based leftmost base of first CIGAR op that consumes reference.
    int32_t pos = (int32_t)b->core.pos;

    // get length of cigar ops that consumes reference and query
    uint32_t qr_len = 0;
    uint32_t ci;
    for (ci = 0; ci < n_cigar; ++ci){
        uint32_t c_op = (uint32_t)bam_cigar_op(cigar_raw[ci]); // lower 4 bits is cigar op
        uint32_t c_len = (uint32_t)bam_cigar_oplen(cigar_raw[ci]); // higher 24 bits is length

        int crq = bam_cigar_type(c_op)&3; // consumes both

        if (crq) qr_len += c_len;
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

        ml_t(exon_list) *exons = &iso->exons;
        size_t n_exons = ml_size(exons);
        if (n_exons == 0){
            iso_i++;
            continue;
        }

        // position of CIGAR segment in reference sequence
        int32_t r_beg = pos, r_end = pos;

        for (ci = 0; ci < n_cigar; ++ci){ // for each CIGAR op
            uint32_t c_op = (uint32_t)bam_cigar_op(cigar_raw[ci]); // lower 4 bits is cigar op
            uint32_t c_len = (uint32_t)bam_cigar_oplen(cigar_raw[ci]); // higher 24 bits is length

            int cr = bam_cigar_type(c_op)&2; // consumes reference
            int cq = bam_cigar_type(c_op)&1; // consumes query
            int crq = bam_cigar_type(c_op)&3; // consumes both

            if (cr == 0 && cq == 0)
                continue;

            // set begin ref position of CIGAR seg to previous seg ref end position
            r_beg = r_end;

            // if consumes reference, add cigar length to get end position 
            // of CIGAR seg in reference seq.
            // if doesn't, e.g. insertion, r_end is still equal to r_beg
            if (cr)
                r_end = r_beg + c_len;

            int64_t e_ovrlp = 0, i_ovrlp = 0;
            ml_node_t(exon_list) *en;
            for (en = ml_begin(exons); en; en = ml_node_next(en)){
                exon_t ex1 = ml_node_val(en);

                /*********************************************************
                 * check exon overlap
                 *********************************************************/
                if (crq){
                    int64_t e_tmp_ovrlp = bp_overlap(r_beg, r_end, strand, 
                            ex1.beg, ex1.end, g->strand);
                    if (e_tmp_ovrlp < 0){
                        *ret = -1;
                        return 0;
                    }
                    e_ovrlp += e_tmp_ovrlp;
                }

                ml_node_t(exon_list) *en_next = ml_node_next(en);
                if (en_next){
                    exon_t ex2 = ml_node_val(en_next);
                    assert(ex1.end <= ex2.beg);

                    /*********************************************************
                     * check splice junction N
                     *********************************************************/
                    if (c_op == BAM_CREF_SKIP && ex1.end == r_beg && ex2.beg == r_end){
                        sj[iso_i] = 1;
                        break;
                    }

                    /*********************************************************
                     * check intron overlap
                     *********************************************************/
                    if (crq){
                        int64_t i_tmp_ovrlp = bp_overlap(r_beg, r_end, strand, 
                                ex1.end, ex2.beg, g->strand);
                        if (i_tmp_ovrlp < 0){
                            *ret = -1;
                            return 0;
                        }
                        i_ovrlp += i_tmp_ovrlp;
                    }
                }

                if (ex1.beg >= r_end)
                    break;
            }
            ero[iso_i] += (double)e_ovrlp; // add overlap of CIGAR segment
            iro[iso_i] += (double)i_ovrlp; // add overlap of CIGAR segment

            // isoform overlap
            if (crq)
                rl[iso_i] += (double)bp_overlap(r_beg, r_end, strand, 
                        iso->beg, iso->end, g->strand);
        }

        pt[iso_i] = qr_len ? (double)rl[iso_i] / (double)qr_len : 0;
        assert(iro[iso_i] <= rl[iso_i]);
        assert(ero[iso_i] <= rl[iso_i]);
        if (rl[iso_i] > 0){
            pi[iso_i] = rl[iso_i] >= 1 ? iro[iso_i] / rl[iso_i] : 0;
            pe[iso_i] = rl[iso_i] >= 1 ? ero[iso_i] / rl[iso_i] : 0;
        }
        iso_i++;
    }

    /*********************************************************
     * get spliced or unspliced status of record for each isoform
     *********************************************************/

    // TODO: check how to handle reads that partially overal an isoform.
    int n_o_iso = 0; // number isoforms with full overlap with UMI
    int is_spl = 0, is_unspl = 0;
    for (ix = 0; ix < n_iso; ++ix){
        if (pt[ix] < 1) continue;
        ++n_o_iso;
        if (sj[ix] == 1) // if overlaps splice junction
            ++is_spl;
        if (pi[ix] > 0) // if at least one base pair overlaps intron
            ++is_unspl;
    }

    /*********************************************************
     * if read overlaps a splice junction from any isoform, then it is 
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
        ml_t(vcfr_list) *vars){
    if (b == NULL || h == NULL || gv == NULL || vars == NULL)
        return err_msg(-1, 0, "bam1_vars_overlap: argument is null");

    uint32_t *cigar_raw = (uint32_t *)bam_get_cigar(b);
    uint32_t n_cigar = (uint32_t)b->core.n_cigar;

    int32_t tid = b->core.tid;
    const char *ref = sam_hdr_tid2name(h, (int)tid);
    int32_t left_pos = (int32_t)b->core.pos; // position of first base that consumes the reference.

    int n_vars = 0;

    // if chromosome not present in g_var_t, no variants to overlap
    if (str_map_ix(gv->chrm_ix, (char *)ref) < 0)
        return 0;

    /* r_pos stores position in ref. sequence 0-based */
    int32_t r_pos_beg = left_pos;
    int32_t r_pos_end = left_pos;
    uint32_t ci;
    for (ci = 0; ci < n_cigar; ci++){
        
        uint32_t cigar_op = (uint32_t)bam_cigar_op(cigar_raw[ci]);
        uint32_t cigar_oplen = (uint32_t)bam_cigar_oplen(cigar_raw[ci]);

        int cr = bam_cigar_type(cigar_op)&2; // consumes ref seq.
        int cq = bam_cigar_type(cigar_op)&1; // consumes query seq.
        if (cr){
            r_pos_end += cigar_oplen;
        }
        // if cigar doesn't consume query and reference, we can't match the query to the 
        // reference and there is no overlapping base.

        // get overlapping variants from r_pos region.
        if ( cr && cq ){
            int ret = g_var_get_region_vars(gv, ref, r_pos_beg, r_pos_end, vars);
            if (ret < -1)
                return err_msg(-1, 0, "bam1_vars_overlap: failed region overlap at "
                        "%s:%"PRIi32"-%"PRIi32"\n", ref, r_pos_beg, r_pos_end);
            else n_vars += ret;
        }

        r_pos_beg = r_pos_end;
    }
    return(n_vars);
}

int bam1_seq_base(const sam_hdr_t *h, bam1_t *b, g_var_t *gv, ml_t(seq_base_l) *bl){
    if (h == NULL || b == NULL || gv == NULL || bl == NULL)
        return err_msg(-1, 0, "bam1_seq_base: argument is null");

    int vret;
    ml_t(vcfr_list) ovars;
    ml_init(vcfr_list, &ovars);
    if ( (vret = bam1_vars_overlap(h, b, gv, &ovars)) < 0 )
        return(-1);

    int n_base = 0;
    ml_node_t(vcfr_list) *vn;
    for (vn = ml_begin(&ovars); vn; vn = ml_node_next(vn)){
        var_t tvar = ml_node_val(vn);
        // Get position of the variant
        int32_t v_rid = tvar.b->rid;
        int32_t v_pos = tvar.b->pos;

        // Get base pair at overlapping site
        uint8_t base, qual;
        if (bam1_site_base(b, b->core.tid, v_pos, &base, &qual) < 0)
            return(-1);

        // create seq_base_t object
        seq_base_t sb;
        seq_base_init(&sb);
        sb.base = base;
        sb.qual = qual;
        sb.pos.rid = v_rid;
        sb.pos.pos = v_pos;
        sb.pos.strand = '.';

        // TODO: better handle multi-allelic SNPs. Currently just skips adding 
        // duplicate bases if there are two or more variants at the same position.
        if ( seq_base_l_insert(bl, sb, 1, 1) < 0 )
            return(-1);

        ++n_base;
    }

    ml_free(vcfr_list, &ovars);

    return n_base;
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
    assert(o >= 0);

    return o;
}

