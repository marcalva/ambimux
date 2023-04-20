
#include "rna_data.h"
#include <stdlib.h>
#include <string.h>
#include "str_util.h"
#include "variants.h"
#include "region.h"
#include "counts.h"


/*******************************************************************************
 * rna_read1_t
 ******************************************************************************/

void rna_read1_init(rna_read1_t *r){
    if (r == NULL) return;
    init_g_region(&r->loc);
    seq_base_l_init(&r->bl);
    seq_gene_l_init(&r->gl);
}

rna_read1_t *rna_read1_alloc(){
   rna_read1_t *r = malloc(sizeof(rna_read1_t));
   if (r == NULL){
       err_msg(-1, 0, "rna_read1_alloc: %s", strerror(errno));
       return(NULL);
   }
   rna_read1_init(r);
   return(r);
}

void rna_read1_free(rna_read1_t *r){
    if (r == NULL) return;
    seq_base_l_free(&r->bl);
    seq_gene_l_free(&r->gl);
}

void rna_read1_dstry(rna_read1_t *r){
    if (r == NULL) return;
    rna_read1_free(r);
    free(r);
}

rna_read1_t *rna_read1_cpy(const rna_read1_t *r, int *ret){
    *ret = 0;
    if (r == NULL)
        return(NULL);

    rna_read1_t *cpy = rna_read1_alloc();
    if (cpy == NULL){
        *ret = -1;
        return(NULL);
    }

    cpy->loc = r->loc;

    if (seq_base_l_cpy(&cpy->bl, &r->bl) < 0){
        *ret = -1;
        free(cpy);
        return(NULL);
    }

    if (seq_gene_l_cpy(&cpy->gl, &r->gl) < 0){
        *ret = -1;
        free(cpy);
        return(NULL);
    }

    return(cpy);
}

int rna_read1_cmp(const rna_read1_t *r1, const rna_read1_t *r2, int cmp_qual){
    int rc = regioncmp(r1->loc, r2->loc);
    if (rc != 0) return(rc);

    int bsc = ml_size(&r1->bl) - ml_size(&r2->bl);
    if (bsc != 0) return(bsc);

    int gsc = ml_size(&r1->gl) - ml_size(&r2->gl);
    if (gsc != 0) return(gsc);

    int bc = seq_base_l_cmp(r1->bl, r2->bl, cmp_qual);
    if (bc != 0) return(bc);

    int gc = seq_gene_l_cmp(r1->gl, r2->gl);
    if (gc != 0) return(gc);

    return(0);
}

int rna_read1_equal(const rna_read1_t *r1, const rna_read1_t *r2, int cmp_qual){
    if (r1 == NULL || r2 == NULL)
        return err_msg(-1, 0, "rna_read1_equal: r1 or r2 are null");

    if (regioncmp(r1->loc, r2->loc) != 0)
        return(0);

    if (ml_size(&r1->bl) != ml_size(&r2->bl))
        return(0);

    if (ml_size(&r1->gl) != ml_size(&r2->gl))
        return(0);

    if ( seq_base_l_equal(r1->bl, r2->bl, cmp_qual) != 1 )
        return(0);

    if ( seq_gene_l_equal(r1->gl, r2->gl) != 1 )
        return(0);

    return(1);
}

int rna_read1_add_gene(rna_read1_t *r, const seq_gene_t *gene){
    if (r == NULL || gene == NULL)
        return err_msg(-1, 0, "rna_read1_add_gene: argument is null");

    if ( seq_gene_l_insert(&r->gl, *gene, 0, 0) < 0 ) return(-1);

    return(0);
}

int rna_read1_add_base(rna_read1_t *r, seq_base_t base){
    if (r == NULL)
        return err_msg(-1, 0, "rna_read1_add_base: argument is null");

    if (seq_base_l_insert(&r->bl, base, 1, 1) < 0) return(-1);

    return(0);
}

int rna_read1_match_qual(rna_read1_t *r, const rna_read1_t *cmp){
    if (r == NULL || cmp == NULL)
        return err_msg(-1, 0, "rna_read1_match_qual: argument is null");

    if (seq_base_l_match_qual(&r->bl, &cmp->bl) < 0)
        return(-1);

    return(0);
}

/*******************************************************************************
 * rna_dups_t
 ******************************************************************************/

int rna_dups_init(rna_dups_t *rd){
    if (rd == NULL) return(0);

    rd->size = 0;
    rd->m = 1;
    rd->dups = (rna_read1_t *)calloc(rd->m, sizeof(rna_read1_t));
    rd->rds_per_dup = (uint32_t *)calloc(rd->m, sizeof(uint32_t));
    if (rd->dups == NULL || rd->rds_per_dup == NULL)
        return err_msg(-1 , 0, "rna_dups_init: %s", strerror(errno));

    return(0);
}

rna_dups_t *rna_dups_alloc(){
    rna_dups_t *rd = (rna_dups_t *)calloc(1, sizeof(rna_dups_t));
    if (rd == NULL){
        err_msg(-1, 0, "rna_dups_alloc: %s", strerror(errno));
        return(NULL);
    }

    if (rna_dups_init(rd) < 0)
        return(NULL);
    
    return(rd);
}

void rna_dups_free(rna_dups_t *rd){
    if (rd == NULL) return;
    uint32_t i;
    for (i = 0; i < rd->size; ++i)
        rna_read1_free(&rd->dups[i]);

    free(rd->dups);
    rd->dups = NULL;
    free(rd->rds_per_dup);
    rd->rds_per_dup = NULL;
    rd->size = 0;
    rd->m = 0;
}

void rna_dups_dstry(rna_dups_t *rd){
    if (rd == NULL) return;

    rna_dups_free(rd);

    free(rd);
}

int rna_dups_add_read(rna_dups_t *rd, const rna_read1_t *r){
    if (rd == NULL || r == NULL)
        return err_msg(-1, 0, "rna_dups_add_read: arguments are null");

    if (rd->dups == NULL || rd->rds_per_dup == NULL)
        return err_msg(-1, 0, "rna_dups_add_read: rd->dups or rds_per_dup is null");

    uint32_t i;
    int ret;
    for (i = 0; i < rd->size; ++i){
        int rcmp = rna_read1_cmp(&rd->dups[i], r, 0);
        if (rcmp == 0){
            if (rna_read1_match_qual(&rd->dups[i], r) < 0)
                return(-1);
            ++rd->rds_per_dup[i];
            return(0);
        }
    }
    while (rd->size >= rd->m){
        rd->m = rd->size + 1;
        rd->dups = realloc(rd->dups, rd->m * sizeof(rna_read1_t));
        rd->rds_per_dup = realloc(rd->rds_per_dup, rd->m * sizeof(uint32_t));
        if (rd->dups == NULL || rd->rds_per_dup == NULL)
            return err_msg(-1, 0, "rna_dups_add_read: %s", strerror(errno));
    }

    // copy read pair
    rna_read1_t *tmp = rna_read1_cpy(r, &ret);
    if (ret < 0)
        return(-1);

    assert(tmp->loc.rid >= 0);

    rd->dups[rd->size] = *tmp;
    rd->rds_per_dup[rd->size] = 1;
    free(tmp);

    ++rd->size;
    return(0);
}

rna_mol_t *rna_dups_dedup(rna_dups_t *dups, int *ret){
    *ret = 0;
    if (dups == NULL) {
        *ret = err_msg(-1, 0, "rna_dups_dedup: mol is null");
        return NULL;
    }

    // if no dups, return 0 successfully
    if (dups->size == 0)
        return NULL;

    int ix_best = 0;
    uint32_t max_c = 0; // store read count
    uint32_t rd_w_max = 0; // number of read pairs with max_c read counts

    uint32_t i;
    for (i = 0; i < dups->size; ++i){
        // check for bugs
        assert(dups->rds_per_dup[i] > 0);

        if (dups->rds_per_dup[i] == max_c){
            ++rd_w_max;
        } else if (dups->rds_per_dup[i] > max_c){
            max_c = dups->rds_per_dup[i];
            ix_best = i;
            rd_w_max = 1;
        }
    }

    // check for bugs
    assert(max_c > 0 && rd_w_max > 0);

    // initialize molecule
    rna_read1_t rb = dups->dups[ix_best];

    rna_mol_t *mol = calloc(1, sizeof(rna_mol_t));
    // rna_mol_t *mol = rna_mol_alloc();
    if (mol == NULL) {
        *ret = -1;
        return NULL;
    }

    // rna_dups_free(mol->dups);
    mol->loc = rb.loc;

    seq_base_l_init(&mol->bl);
    if ( seq_base_l_cpy(&mol->bl, &rb.bl) < 0 ) {
        *ret = err_msg(-1, 0, "rna_dups_dedup: failed to copy base list");
        free(mol);
        return NULL;
    }

    seq_gene_l_init(&mol->gl);
    if ( seq_gene_l_cpy(&mol->gl, &rb.gl) < 0 ) {
        *ret = err_msg(-1, 0, "rna_dups_dedup: failed to copy gene list");
        free(mol);
        return NULL;
    }

    seq_vac_l_init(&mol->vl);
    mol->n_reads = dups->rds_per_dup[ix_best];

    return mol;
}

/*******************************************************************************
 * rna_dups_bag_t
 ******************************************************************************/

int rna_dup_node_init(rna_dup_node *node) {
    if (node == NULL)
        return err_msg(-1, 0, "rna_dup_node_init: argument is null");

    node->key = 0;
    if ( rna_dups_init(&node->rna_dup) < 0 )
        return -1;

     return 0;
}

rna_dup_node *rna_dup_node_alloc() {
    rna_dup_node *node = calloc(1, sizeof(rna_dup_node));
    if (node == NULL) {
        err_msg(-1, 0, "rna_dup_node_alloc: %s", strerror(errno));
        return NULL;
    }
    if (rna_dup_node_init(node) < 0)
        return NULL;
    return node;
}

void rna_dup_node_free(rna_dup_node *node) {
    if (node == NULL) return;
    rna_dups_free(&node->rna_dup);
}

int rna_dups_bag_init(rna_dups_bag_t *bag) {
    if (bag == NULL) return -1;
    bag->rna_dups = NULL;
    return 0;
}

void rna_dups_bag_free(rna_dups_bag_t *bag) {
    if (bag == NULL) return;
    if (bag->rna_dups == NULL) return;

    kavl_itr_t(bt_umi_dup) itr;
    kavl_itr_first(bt_umi_dup, bag->rna_dups, &itr);  // place at first
    do {                             // traverse
        rna_dup_node *p = (rna_dup_node *)kavl_at(&itr);
        rna_dup_node_free(p);
        free((void*)p);                // free node
    } while (kavl_itr_next(bt_umi_dup, &itr));
    bag->rna_dups = NULL;
}

int rna_dups_bag_add_read(rna_dups_bag_t *bag, const rna_read1_t *r,
        umishort umih) {
    if (bag == NULL || r == NULL)
        return err_msg(-1, 0, "rna_dups_bag_add_read: argument is null");

    rna_dup_node *node_i, *node_j = rna_dup_node_alloc();
    if (node_j == NULL)
        return -1;
    node_j->key = umih;
    node_i = kavl_insert(bt_umi_dup, &bag->rna_dups, node_j, NULL);
    if (node_i != node_j) {
        rna_dup_node_free(node_j);
        free(node_j);
    }

    if (rna_dups_add_read(&node_i->rna_dup, r) < 0)
        return -1;

    return 0;
}

void rna_dups_bag_itr_first(rna_dups_bag_itr *itr, rna_dups_bag_t *bag) {
    assert(itr != NULL || bag != NULL);
    if (bag->rna_dups == NULL) {
        itr->next = 0;
        return;
    }
    itr->next = 1;
    kavl_itr_first(bt_umi_dup, bag->rna_dups, &itr->itr);
}

int rna_dups_bag_itr_next(rna_dups_bag_itr *itr) {
    assert(itr != NULL);
    if (itr->next == 0) return 0;
    itr->next = kavl_itr_next(bt_umi_dup, &itr->itr);
    return(itr->next);
}

int rna_dups_bag_itr_alive(rna_dups_bag_itr *itr) {
    assert(itr != NULL);
    return(itr->next);
}

umishort *rna_dups_bag_itr_key(rna_dups_bag_itr *itr) {
    assert(itr != NULL);
    rna_dup_node *p = (rna_dup_node *)kavl_at(&itr->itr);
    if (p == NULL)
        return NULL;
    return(&p->key);
}

rna_dups_t *rna_dups_bag_itr_val(rna_dups_bag_itr *itr) {
    assert(itr != NULL);
    rna_dup_node *p = (rna_dup_node *)kavl_at(&itr->itr);
    if (p == NULL)
        return NULL;
    return(&p->rna_dup);
}

/*******************************************************************************
 * rna_mlcl_bag_t
 ******************************************************************************/

int rna_mlc_node_init(rna_mlc_node *node) {
    if (node == NULL)
        return err_msg(-1, 0, "rna_mlc_node_init: argument is null");

    node->key = 0;
    if ( rna_mol_init(&node->rna_mol) < 0 )
        return -1;

     return 0;
}

rna_mlc_node *rna_mlc_node_alloc() {
    rna_mlc_node *node = calloc(1, sizeof(rna_mlc_node));
    if (node == NULL) {
        err_msg(-1, 0, "rna_mlc_node_alloc: %s", strerror(errno));
        return NULL;
    }
    if (rna_mlc_node_init(node) < 0)
        return NULL;
    return node;
}

void rna_mlc_node_free(rna_mlc_node *node) {
    if (node == NULL) return;
    rna_mol_free(&node->rna_mol);
}

int rna_mlc_bag_init(rna_mlc_bag_t *bag) {
    if (bag == NULL) return -1;
    bag->rna_mlcs = NULL;
    return 0;
}

void rna_mlc_bag_free(rna_mlc_bag_t *bag) {
    if (bag == NULL) return;
    if (bag->rna_mlcs == NULL) return;

    kavl_itr_t(bt_umi_mlc) itr;
    kavl_itr_first(bt_umi_mlc, bag->rna_mlcs, &itr);  // place at first
    do {                             // traverse
        rna_mlc_node *p = (rna_mlc_node *)kavl_at(&itr);
        rna_mlc_node_free(p);
        free((void*)p);                // free node
    } while (kavl_itr_next(bt_umi_mlc, &itr));
    bag->rna_mlcs = NULL;
}

int rna_mlc_bag_add(rna_mlc_bag_t *bag, rna_mol_t *rna_mol, 
        umishort umih) {
    if (bag == NULL || rna_mol == NULL)
        return err_msg(-2, 0, "rna_mlc_bag_add: argument is null");

    rna_mlc_node *node_i, *node_j = calloc(1, sizeof(rna_mlc_node));
    if (node_j == NULL)
        return err_msg(-2, 0, "rna_mlc_bag_add: %s", strerror(errno));
    node_j->key = umih;
    node_j->rna_mol = *rna_mol;
    node_i = kavl_insert(bt_umi_mlc, &bag->rna_mlcs, node_j, NULL);
    if (node_i != node_j) {
        rna_mlc_node_free(node_j);
        free(node_j);
        return -1;
    }

    return 0;
}

void rna_mlc_bag_itr_first(rna_mlc_bag_itr *itr, rna_mlc_bag_t *bag) {
    assert(itr != NULL || bag != NULL);
    if (bag->rna_mlcs == NULL) {
        itr->next = 0;
        return;
    }
    itr->next = 1;
    kavl_itr_first(bt_umi_mlc, bag->rna_mlcs, &itr->itr);
}

int rna_mlc_bag_itr_next(rna_mlc_bag_itr *itr) {
    assert(itr != NULL);
    if (itr->next == 0) return 0;
    itr->next = kavl_itr_next(bt_umi_mlc, &itr->itr);
    return(itr->next);
}

int rna_mlc_bag_itr_alive(rna_mlc_bag_itr *itr) {
    assert(itr != NULL);
    return(itr->next);
}

umishort *rna_mlc_bag_itr_key(rna_mlc_bag_itr *itr) {
    assert(itr != NULL);
    rna_mlc_node *p = (rna_mlc_node *)kavl_at(&itr->itr);
    if (p == NULL)
        return NULL;
    return(&p->key);
}

rna_mol_t *rna_mlc_bag_itr_val(rna_mlc_bag_itr *itr) {
    assert(itr != NULL);
    rna_mlc_node *p = (rna_mlc_node *)kavl_at(&itr->itr);
    if (p == NULL)
        return NULL;
    return(&p->rna_mol);
}

/*******************************************************************************
 * rna_mol_t
 ******************************************************************************/

int rna_mol_init(rna_mol_t *rmol){
    if (rmol == NULL) return(-1);
    seq_base_l_init(&rmol->bl);
    seq_vac_l_init(&rmol->vl);
    seq_gene_l_init(&rmol->gl);
    rmol->n_reads = 0;
    return(0);
}

rna_mol_t *rna_mol_alloc(){
    rna_mol_t *rmol = (rna_mol_t *)calloc(1, sizeof(rna_mol_t));
    if (rmol == NULL){
        err_msg(-1, 0, "rna_mol_alloc: %s", strerror(errno));
        return(NULL);
    }
    if (rna_mol_init(rmol) < 0)
        return(NULL);
    return(rmol);
}

void rna_mol_free(rna_mol_t *rmol){
    if (rmol == NULL) return;

    seq_base_l_free(&rmol->bl);
    seq_vac_l_free(&rmol->vl);
    seq_gene_l_free(&rmol->gl);
    rmol->n_reads = 0;
}

void rna_mol_dstry(rna_mol_t *rmol){
    if (rmol == NULL) return;
    rna_mol_free(rmol);
    free(rmol);
}

rna_mol_t *rna_mol_cpy(const rna_mol_t *m, int *ret){
    *ret = 0;
    if (m == NULL)
        return(NULL);

    rna_mol_t *cpy = rna_mol_alloc();
    if (cpy == NULL){
        *ret = -1;
        return(NULL);
    }

    cpy->loc = m->loc;

    if (seq_base_l_cpy(&cpy->bl, &m->bl) < 0){
        *ret = -1;
        rna_mol_dstry(cpy);
        return(NULL);
    }

    if (seq_vac_l_cpy(&cpy->vl, &m->vl) < 0){
        *ret = -1;
        rna_mol_dstry(cpy);
        return(NULL);
    }

    if (seq_gene_l_cpy(&cpy->gl, &m->gl) < 0){
        *ret = -1;
        rna_mol_dstry(cpy);
        return(NULL);
    }

    cpy->n_reads = m->n_reads;

    return(cpy);
}

int rna_mol_var_call(rna_mol_t *m, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (m == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "rna_mol_var_call: argument is null");
    }

    int ret = seq_vac_l_call_var(&m->bl, &m->vl, gv, cmap, min_qual);

    // free the bases
    seq_base_l_free(&m->bl);

    return(ret);
}

