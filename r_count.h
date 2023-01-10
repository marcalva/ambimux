
#ifndef R_COUNT_H
#define R_COUNT_H

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <inttypes.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/bgzf.h"
#include "str_util.h"
#include "variants.h"
#include "bam_dat.h"
#include "counts.h"
#include "kbtree.h"

/*******************************************************************************
 * count node
 ******************************************************************************/

typedef struct {
    int ix;
    uint32_t counts[MAX_ALLELE];
} cnt_node_t;

#define cnt_node_cmp(p, q) ((p).ix - (q).ix)
KBTREE_INIT(kh_cnode, cnt_node_t, cnt_node_cmp);

typedef struct {
    kbtree_t(kh_cnode) *b;
} bc_btree_t;

KHASH_INIT(kh_bt, char*, bc_btree_t, 1, kh_str_hash_func, kh_str_hash_equal);

/*******************************************************************************
 * Barcode counts
 ******************************************************************************/

typedef struct  bc_counts_t {
    kbtree_t(kh_cnode) *atac_ac;
    kbtree_t(kh_cnode) *atac_pc;
    kbtree_t(kh_cnode) *rna_ac;
    kbtree_t(kh_cnode) *rna_gc;
} bc_counts_t;

KHASH_INIT(kh_bc_cnode, char *, bc_counts_t *, 1, kh_str_hash_func, kh_str_hash_equal);

/*******************************************************************************
 * BAM counts
 ******************************************************************************/

typedef struct {
    g_var_t *gv; // variants
    str_map *gene_ix; // genes
    str_map *bc_ix; // barcodes
    iregs_t *pks; // ATAC peaks. Don't free

    // key-val for hash tables:
    // key is barcode ID string. Key is same memory as bc_ix. 
    // val is bc_counts_t
    khash_t(kh_bc_cnode) *bc_counts;

    // flags
    uint8_t has_atac_ac;
    uint8_t has_atac_pc;
    uint8_t has_rna_ac;
    uint8_t has_rna_gc;

    int rna_gcs_nz; // number of non-zero elements in rna_gcs
    int rna_acs_nz; // number of non-zero elements in rna_acs
    int atac_acs_nz; // number of non-zero elements in atac_acs
    int atac_pcs_nz; // number of non-zero elements in atac_pcs
} bam_counts_t;

bam_counts_t *bam_counts_init();
void bam_counts_dstry(bam_counts_t *agc);

int bam_counts_add_gv(bam_counts_t *agc, g_var_t *gv);
int bam_counts_add_gene_map(bam_counts_t *agc, str_map *sm);
int bam_counts_add_bc_map(bam_counts_t *agc, str_map *sm);
int bam_counts_add_peaks(bam_counts_t *agc, iregs_t *pks);

int bam_counts_count(bam_counts_t *agc, bam_data_t *bam_data);

int write_num(BGZF *fp, int n, char *c, char **intstrp, size_t *intstrp_m);

int bam_counts_write_atac_ac(bam_counts_t *agc, char *fn);
int bam_counts_write_rna_gc(bam_counts_t *agc, gene_anno_t *anno, char *fn);
int bam_counts_write_rna_ac(bam_counts_t *agc, gene_anno_t *anno, char *fn);
int bam_counts_write_atac_pc(bam_counts_t *agc, char *fn);
int bam_counts_write(bam_counts_t *agc, gene_anno_t *anno, g_var_t *gv, char *fn);

int mtx_write_hdr(BGZF *fp, char *len1, char *len2, char *len3, char *delim);

#endif // R_COUNT_H
