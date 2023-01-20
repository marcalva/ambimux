
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "htslib/sam.h"
#include "htslib/vcf.h"
#include "htslib/vcfutils.h"
#include "htslib/synced_bcf_reader.h"
#include "htslib/hts.h"
#include "sam_read.h"
#include "str_util.h"

int get_rcoord_bam(const bam1_t *b, int32_t *tid, int32_t *start, int32_t *end, 
        int adj_soft_clip){
    if (b == NULL || start == NULL || end == NULL){
        return err_msg(-1, 0, "get_rcoord_bam: one of arguments are NULL pointers");
    }
    *tid = b->core.tid;
    *start = (int32_t)b->core.pos;
    *end = (int32_t)bam_endpos(b);
    // if seq is soft clipped at beg or end, subtract or add the length to start or end
    if (adj_soft_clip){
        uint32_t *cigar_raw;
        cigar_raw = bam_get_cigar(b);
        uint32_t n_cigar = b->core.n_cigar;
        if (n_cigar < 1)
            return(0);
        uint32_t op1, op2, op1len, op2len;
        op1 = bam_cigar_op(cigar_raw[0]);
        op1len = bam_cigar_oplen(cigar_raw[0]);
        if (op1 == BAM_CSOFT_CLIP)
            start -= (int32_t)op1len;
        if (n_cigar > 1){
            op2 = bam_cigar_op(cigar_raw[n_cigar-1]);
            op2len = bam_cigar_oplen(cigar_raw[n_cigar-1]);
            if (op2 == BAM_CSOFT_CLIP)
                end += op2len;
        }
    }
    return(0);
}

int load_bam(const char *bamfn, samFile **bam, sam_hdr_t **bam_hdr, hts_idx_t **bam_idx){

    *bam = sam_open(bamfn, "r");
    if (*bam == NULL)
        return err_msg(-1, 0, "load_bam: could not open BAM file %s", bamfn);

    *bam_hdr = sam_hdr_read(*bam);
    if (*bam_hdr == NULL)
        return err_msg(-1, 0, "load_bam: could not read header for BAM file %s", bamfn);

    *bam_idx = sam_index_load(*bam, bamfn);
    if (*bam_idx == NULL)
        return err_msg(-1, 0, "load_bam: could not load index for BAM file %s", bamfn);

    return(0);
}


// Overlap tid
// t1 and t2 must be freed.
int ovrlp_tid(sam_hdr_t *sam_hdr, bcf_hdr_t *bcf_hdr, int **t1, int **t2){
    int n_sam = sam_hdr->n_targets;
    char **sam_names = sam_hdr->target_name;
    int n_bcf;
    const char **bcf_names = bcf_hdr_seqnames(bcf_hdr, &n_bcf);

    *t1 = (int*)malloc(n_bcf * n_sam * sizeof(int));
    if (*t1 == NULL) return err_msg(-1, 0, "ovrlp_tid: %s", strerror(errno));
    *t2 = (int*)malloc(n_bcf * n_sam * sizeof(int));
    if (*t2 == NULL) return err_msg(-1, 0, "ovrlp_tid: %s", strerror(errno));

    int i,j;
    int n = 0;
    for (i = 0; i < n_sam; i++){
        for (j = 0; j < n_bcf; j++){
            if (strcmp(sam_names[i], bcf_names[j]) != 0)
                continue;
            int ix1 = sam_hdr_name2tid(sam_hdr, sam_names[i]);
            int ix2 = bcf_hdr_name2id(bcf_hdr, bcf_names[j]);
            (*t1)[n] = ix1;
            (*t2)[n] = ix2;
            n++;
            break;
        }
    }
    free(bcf_names);

    (*t1) = realloc((*t1), n * sizeof(int));
    (*t2) = realloc((*t2), n * sizeof(int));
    return n;
}

// Return 0 on no overlapping chromosome
// Return 1 on success
// overlapping base stored in base. Stores 'N' if no overlap.
// overlapping qual stored in qual. Stores 0 if no overlap
int get_ovrlp_base(const bam1_t *b, int32_t base_ref, int32_t base_pos, 
        char *base, uint8_t *qual){
    *base = 'N';
    *qual = 0;
    uint32_t *cigar_raw;
    cigar_raw = bam_get_cigar(b);
    uint32_t n_cigar = b->core.n_cigar;

    int32_t tid = b->core.tid;
    if (tid != base_ref)
        return(0);
    int32_t left_pos = (int32_t)b->core.pos;
    uint32_t i;
    uint32_t prev_q_index = 0;
    uint32_t next_q_index = 0;
    int64_t qlen = bam_cigar2qlen(n_cigar, cigar_raw);
    int32_t prev_r_index = left_pos;
    int32_t next_r_index = left_pos;
    for (i = 0; i < n_cigar; i++){
        if (base_pos < prev_r_index)
            break;
        uint32_t op = bam_cigar_op(cigar_raw[i]);
        uint32_t cigar_oplen = bam_cigar_oplen(cigar_raw[i]);

        int cr = bam_cigar_type(op)&2;
        int cq = bam_cigar_type(op)&1;
        // if consumes reference
        if (cr){
            next_r_index += cigar_oplen;
        }
        // if consumes query
        if (cq){
            next_q_index += cigar_oplen;
        }
        if ( !cq ){
            prev_q_index = next_q_index;
            prev_r_index = next_r_index;
            continue;
        }

        if ( (base_pos >= prev_r_index) && (base_pos < next_r_index)){
            int32_t add_to = base_pos - prev_r_index;
            uint32_t q_ix = prev_q_index + add_to;
            if (q_ix >= qlen)
                return err_msg(-1, 0, "get_ovrlp_base: q index (%llu) is greater than query length (%i), "
                        "probable bug in program", q_ix, qlen);
            uint8_t *p = bam_get_seq(b);
            uint8_t *q = bam_get_qual(b);
            *base = seq_nt16_str[bam_seqi(p, q_ix)];
            *qual = q[q_ix];
        }
        else {
            prev_q_index = next_q_index;
            prev_r_index = next_r_index;
        }
    }
    return(1);
}

int bam1_site_base(const bam1_t *b, int32_t ref, int32_t pos, 
        uint8_t *base, uint8_t *qual){

    uint32_t *cigar_raw = bam_get_cigar(b);
    uint32_t n_cigar = b->core.n_cigar;

    int32_t tid = b->core.tid;
    int32_t left_pos = (int32_t)b->core.pos; // position of first base that consumes the reference.
    int32_t right_pos = (int32_t)bam_endpos(b); // position after last base that consumes reference.
    
    if (tid != ref)
        return(0);
    if ((pos < left_pos) || (pos >= right_pos))
        return(0);
    uint32_t i;

    /* q_index stores position in query sequence
     * r_index stores position in ref. sequence
     * both are 0-based, [beg, end) */
    uint64_t q_pos_beg = 0;
    uint64_t q_pos_end = 0;
    int qlen = bam_cigar2qlen(n_cigar, cigar_raw);
    int32_t r_pos_beg = left_pos;
    int32_t r_pos_end = left_pos;
    for (i = 0; i < n_cigar; i++){
        // If we passed the seeked position
        if (pos < r_pos_beg)
            break;
        
        uint32_t cigar_op = bam_cigar_op(cigar_raw[i]);
        uint32_t cigar_oplen = bam_cigar_oplen(cigar_raw[i]);

        int cr = bam_cigar_type(cigar_op)&2; // consumes ref seq.
        int cq = bam_cigar_type(cigar_op)&1; // consumes query seq.
        if (cr){
            r_pos_end += cigar_oplen;
        }
        if (cq){
            q_pos_end += cigar_oplen;
        }
        // if cigar doesn't consume query and reference, we can't match the query to the 
        // reference and there is no overlapping base.
        if ( !(cr && cq) ){
            q_pos_beg = q_pos_end;
            r_pos_beg = r_pos_end;
            continue;
        }

        if ( (pos >= r_pos_beg) && (pos < r_pos_end)){
            int32_t add_to = pos - r_pos_beg;
            uint64_t q_ix = q_pos_beg + add_to;
            if (q_ix >= qlen)
                return err_msg(-1, 0, "bam1_site_base: q index (%llu) is greater than query length (%i), "
                        "probable bug in program", q_ix, qlen);
            uint8_t *p = bam_get_seq(b);
            uint8_t *q = bam_get_qual(b);
            *base = bam_seqi(p, q_ix); // returns 4-bit integer encoding of base
            *qual = q[q_ix];
            return(1);
        }
        else {
            q_pos_beg = q_pos_end;
            r_pos_beg = r_pos_end;
        }
    }
    return(0);
}

// Returns string from tag, NULL if not found.
char *get_tag(const bam1_t *b, char tag[2]){
    uint8_t *ptr = bam_aux_get(b, tag);
    if (ptr == NULL)
        return NULL;
    char *bc = bam_aux2Z(ptr);
    return bc;
}

// Returns integer from tag, NULL if not found.
int64_t get_itag(const bam1_t *b, char tag[2]){
    uint8_t *ptr = bam_aux_get(b, tag);
    int64_t ret = bam_aux2i(ptr);
    return ret;
}

void print_bam1_t(const bam1_t *b){
    char *qname = bam_get_qname(b);
    int32_t pos = (int32_t)b->core.pos;

    fprintf(stdout, "qname=%s; pos=%i; seq=", qname, (int)pos);
    int i;
    for (i = 0; i < b->core.l_qseq; ++i){
        uint8_t *p = bam_get_seq(b);
        char base = seq_nt16_str[bam_seqi(p, i)];
        fprintf(stdout, "%c", base);
    }
    fprintf(stdout, "\n");
}

