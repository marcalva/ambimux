
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
KBTREE_INIT(str, cnt_node_t, cnt_node_cmp);

typedef struct {
    kbtree_t(str) *b;
} bc_btree_t;

KHASH_INIT(kh_bt, char*, bc_btree_t, 1, kh_str_hash_func, kh_str_hash_equal);

/*******************************************************************************
 * allele counts
 ******************************************************************************/

/* hold allele counts for a variant */
typedef struct ac_node {
    int ix; // variant index from gv
    uint32_t counts[MAX_ALLELE]; // ref, alt1, alt2, ..., N
    struct ac_node *next; // pointer to next ac_node, NULL if none
} ac_node;

// value in hash is a dummy head node with index=-1
KHASH_INIT(ac, char*, ac_node, 1, kh_str_hash_func, kh_str_hash_equal);

ac_node *init_ac_node();
ac_node *destroy_ac_node(ac_node *n);
void ac_node_set_zero(ac_node *n);

/*******************************************************************************
 * gene counts
 ******************************************************************************/

/* hold counts for a gene */
typedef struct gc_node {
    int ix; // gene index from gene_ix
    uint32_t counts[N_SPL];
    struct gc_node *next; // pointer to next gc_node, NULL if none
} gc_node;

// value in hash is a dummy head node with index=-1
KHASH_INIT(gc, char*, gc_node, 1, kh_str_hash_func, kh_str_hash_equal);

gc_node *init_gc_node();
gc_node *destroy_gc_node(gc_node *n);
void gc_node_set_zero(gc_node *n);

/*******************************************************************************
 * barcode counts
 ******************************************************************************/

typedef struct {
    GenomeVar *gv; // variants
    str_map *gene_ix; // genes
    str_map *bc_ix; // barcodes
    iregs_t *pks; // ATAC peaks. Don't free

    // key-val for hash tables:
    // key is barcode ID string. Key is same memory as bc_ix. 
    // val is gc_node
    // b-tree
    khash_t(kh_bt) *bt_atac_ac; // allele counts per SNP (ATAC)
    khash_t(kh_bt) *bt_atac_pc; // frag counts per peak (ATAC)
    khash_t(kh_bt) *bt_rna_ac; // allele counts per SNP (RNA)
    khash_t(kh_bt) *bt_rna_gc; // allele counts per SNP (ATAC)

    int rna_gcs_nz; // number of non-zero elements in rna_gcs
    int rna_acs_nz; // number of non-zero elements in rna_acs
    int atac_acs_nz; // number of non-zero elements in atac_acs
    int atac_pcs_nz; // number of non-zero elements in atac_pcs
} bam_ag_t;

bam_ag_t *bam_ag_init();
void bam_ag_dstry(bam_ag_t *agc);

int bam_ag_add_gv(bam_ag_t *agc, GenomeVar *gv);
int bam_ag_add_gene_map(bam_ag_t *agc, str_map *sm);
int bam_ag_add_bc_map(bam_ag_t *agc, str_map *sm);
int bam_ag_add_peaks(bam_ag_t *agc, iregs_t *pks);

int write_num(BGZF *fp, int n, char *c, char **intstrp, size_t *intstrp_m);
int bam_rna_gc_count(bam_ag_t *agc, bam_data_t *bam_dat);
int bam_atac_ac_count(bam_ag_t *agc, bam_data_t *bam_dat);
int bam_rna_ac_count(bam_ag_t *agc, bam_data_t *bam_dat);
int bam_atac_pc_count(bam_ag_t *agc, bam_data_t *bam_dat);

int bam_ag_write_atac_ac(bam_ag_t *agc, char *fn);
int bam_ag_write_rna_gc(bam_ag_t *agc, Annotation *anno, char *fn);
int bam_ag_write_rna_ac(bam_ag_t *agc, Annotation *anno, char *fn);
int bam_ag_write_atac_pc(bam_ag_t *agc, char *fn);
int bam_ag_write(bam_ag_t *agc, Annotation *anno, GenomeVar *gv, char *fn);

int mtx_write_hdr(BGZF *fp, char *len1, char *len2, char *len3, char *delim);

#endif // R_COUNT_H
