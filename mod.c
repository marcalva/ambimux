
#include "mod.h"
#include "str_util.h"
#include "rna_data.h"
#include "clopts.h"
#include "math_util.h"
#include "array_util.h"
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include <pthread.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include "bits.h"
#include "g_list.h"
#include "rna_data.h"
#include "atac_data.h"

#define ATAC_IX 0
#define RNA_IX 1

/*******************************************************************************
 * num check
 ******************************************************************************/

static int num_invalid(f_t x) {
    if (isnan(x))
        return 1;
    if (isinf(x))
        return 1;
    return 0;
}

static int prob_invalid(f_t x) {
    if (num_invalid(x))
        return 1;
    if (x < 0.0 && x != -0.0 || x > 1.0)
        return 1;
    return 0;
}

static int psum_invalid(f_t x, f_t lb, f_t ub) {
    if (num_invalid(x))
        return 1;
    if (x < lb || x > ub)
        return 1;
    return 0;
}

static int iszero(f_t x) {
    if (x <= 0.0 && x >= 0.0)
        return 1;
    return 0;
}

/*******************************************************************************
 * math
 ******************************************************************************/

int f_t_lt_cmp(const void *a, const void *b) {
    f_t arg1 = *(const f_t *)a;
    f_t arg2 = *(const f_t *)b;
    return (arg1 > arg2) - (arg1 < arg2);
}

int proj_splx(f_t *x, f_t *xp, size_t n) {
    if (x == NULL || xp == NULL)
        return err_msg(-1, 0, "proj_splx: argument is null");

    if (n <= 1)
        return err_msg(-1, 0, "proj_splx: n<=1 makes no sense");

    size_t i;
    f_t *y = (f_t *)malloc(n * sizeof(f_t));
    if (y == NULL) {
        return err_msg(-1, 0, "proj_splx: %s", strerror(errno));
    }
    for (i = 0; i < n; ++i) {
        y[i] = x[i];
    }

    qsort(y, n, sizeof(f_t), f_t_lt_cmp);

    f_t sum = 0.0;
    f_t t_i, t_hat;
    for (i = n; i > 0; --i) {
        sum += y[i - 1];
        t_i = (sum - 1) / (n - (i - 1));
        if (i == 1 || t_i >= y[i - 2]) {
            t_hat = t_i;
            break;
        }
    }

    for (i = 0; i < n; ++i) {
        xp[i] = (x[i] - t_hat) > 0 ? (x[i] - t_hat) : 0;
    }

    free(y);

    return 0;
}

/*******************************************************************************
 * seq bases
 ******************************************************************************/

f_t phred_to_perr(uint8_t phred) {
    return pow(10.0, -(f_t)phred / 10.0);
}

/*******************************************************************************
 * mdl_mlcl_t
 ******************************************************************************/

int mdl_mlcl_cmp(mdl_mlcl_t m1, mdl_mlcl_t m2){
    size_t i;
    if (mv_size(&m1.feat_ixs) != mv_size(&m2.feat_ixs))
        return( mv_size(&m1.feat_ixs) - mv_size(&m2.feat_ixs));

    for (i = 0; i < mv_size(&m1.feat_ixs); ++i){
        if (mv_i(&m1.feat_ixs, i) != mv_i(&m2.feat_ixs, i))
            return( mv_i(&m1.feat_ixs, i) - mv_i(&m2.feat_ixs, i) );
    }

    if (mv_size(&m1.var_ixs) != mv_size(&m2.var_ixs))
        return( mv_size(&m1.var_ixs) - mv_size(&m2.var_ixs));

    for (i = 0; i < mv_size(&m1.var_ixs); ++i){
        if (mv_i(&m1.var_ixs, i) != mv_i(&m2.var_ixs, i))
            return( mv_i(&m1.var_ixs, i) - mv_i(&m2.var_ixs, i) );
    }

    if (mv_size(&m1.bquals) != mv_size(&m2.bquals))
        return( mv_size(&m1.bquals) - mv_size(&m2.bquals));

    for (i = 0; i < mv_size(&m1.bquals); ++i){
        if (mv_i(&m1.bquals, i) != mv_i(&m2.bquals, i))
            return( mv_i(&m1.bquals, i) - mv_i(&m2.bquals, i) );
    }

    return(0);
}

void mdl_mlcl_init(mdl_mlcl_t *mlcl){
    mv_init(&mlcl->feat_ixs);
    mv_init(&mlcl->var_ixs);
    mv_init(&mlcl->bquals);
    mlcl->counts = 0;
}
void mdl_mlcl_free(mdl_mlcl_t *mlcl){
    if (mlcl == NULL) return;
    mv_free(&mlcl->feat_ixs);
    mv_free(&mlcl->var_ixs);
    mv_free(&mlcl->bquals);
    mlcl->counts = 0;
}

void mdl_mlcl_print(FILE *f, mdl_mlcl_t *mlcl){
    fprintf(f, "\tmolecule: ");
    size_t i;
    fprintf(f, "features:");
    for (i = 0; i < mv_size(&mlcl->feat_ixs); ++i){
        fprintf(f, ", %i", mv_i(&mlcl->feat_ixs, i));
    }
    fprintf(f, " variants:");
    for (i = 0; i < mv_size(&mlcl->var_ixs); ++i){
        fprintf(f, ", %i (%u)", mv_i(&mlcl->var_ixs, i), mv_i(&mlcl->bquals, i));
    }
    fprintf(f, ", count=%u", mlcl->counts);
    fprintf(f, "\n");
}

uint32_t mdl_mlcl_tot_count(kbtree_t(kb_mdl_mlcl) *bt){
    kbitr_t itr;
    uint32_t counts = 0;
    kb_itr_first(kb_mdl_mlcl, bt, &itr); 
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, bt, &itr)){
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        counts += mlcl->counts;
    }
    return(counts);
}

int mdl_mlcl_info_count(kbtree_t(kb_mdl_mlcl) *bt, khash_t(iset) *var_ixs, uint32_t *counts) {
    if (bt == NULL || var_ixs == NULL || counts == NULL)
        return err_msg(-1, 0, "mdl_mlcl_var_flg: argument is null");

    *counts = 0;
    int khret;
    kbitr_t itr;
    kb_itr_first(kb_mdl_mlcl, bt, &itr); 
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, bt, &itr)){
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        size_t i, n_var = mv_size(&mlcl->var_ixs);
        for (i = 0; i < n_var; ++i) {
            uint32_t pack = mv_i(&mlcl->var_ixs, i);
            uint32_t v_ix;
            uint8_t allele;
            mdl_mlcl_unpack_var(pack, &v_ix, &allele);
            kh_put(iset, var_ixs, v_ix, &khret);
            if (khret < 0)
                return err_msg(-1, 0, "mdl_mlcl_var_flg: failed to push '%i' to set", v_ix);
        }
        if (n_var > 0)
            *counts += mlcl->counts;
    }
    return 0;
}

int mdl_mlcl_pack_var(uint32_t var_ix, uint8_t allele, uint32_t *pack) {
    uint32_t max_var = 0x1<<28;
    if (var_ix >= max_var)
        return err_msg(-1, 0, "mdl_mlcl_pack_var: cannot pack var_ix=%u", var_ix);

    *pack = (var_ix << 4) | allele;

    return 0;
}

void mdl_mlcl_unpack_var(uint32_t pack, uint32_t *var_ix, uint8_t *allele) {
    *var_ix = pack >> 4;
    *allele = pack & 0xf;
}

int mdl_mlcl_bc_info_count(mdl_mlcl_bc_t *mlcl_bc, size_t n_var,
                           uint32_t *rna_count, uint32_t *atac_count,
                           uint32_t *rna_var_count, uint32_t *atac_var_count) {
    if (mlcl_bc == NULL)
        return err_msg(-1, 0, "mdl_mlcl_bc_info: argument is null");

    if (rna_count == NULL || atac_count == NULL || 
        rna_var_count == NULL || atac_var_count == NULL)
        return err_msg(-1, 0, "mdl_mlcl_bc_info: argument is null");

    khash_t(iset) *rna_vars = kh_init(iset);
    if (rna_vars == NULL)
        return err_msg(-1, 0, "mdl_mlcl_bc_info: failed to initialize rna_vars");
    khash_t(iset) *atac_vars = kh_init(iset);
    if (atac_vars == NULL)
        return err_msg(-1, 0, "mdl_mlcl_bc_info: failed to initialize atac_vars");

    if (kh_resize(iset, rna_vars, n_var) < 0)
        return err_msg(-1, 0, "mdl_mlcl_bc_info: failed resize");
    if (kh_resize(iset, atac_vars, n_var) < 0)
        return err_msg(-1, 0, "mdl_mlcl_bc_info: failed resize");
    
    if (mdl_mlcl_info_count(mlcl_bc->rna, rna_vars, rna_count) < 0)
        return err_msg(-1, 0, "mdl_mlcl_info_count: failed to get rna counts");
    if (mdl_mlcl_info_count(mlcl_bc->atac, atac_vars, atac_count) < 0)
        return err_msg(-1, 0, "mdl_mlcl_info_count: failed to get atac counts");

    *rna_var_count = kh_size(rna_vars);
    *atac_var_count = kh_size(atac_vars);

    kh_destroy(iset, rna_vars);
    rna_vars = NULL;
    kh_destroy(iset, atac_vars);
    atac_vars = NULL;

    return 0;
}

int mdl_mlcl_add_rna(mdl_mlcl_t *mlcl, rna_mol_t *mol,
        int n_genes) {
    if (mlcl == NULL || mol == NULL)
        return err_msg(-1, 0, "mdl_mlcl_add_rna: argument is null");

    size_t mol_n_genes = ml_size(&mol->gl);
    // add gene if only 1 is present, skip multi-gene RNAs
    if (mol_n_genes == 1){
        ml_node_t(seq_gene_l) *g_n;
        for (g_n = ml_begin(&mol->gl); g_n; g_n = ml_node_next(g_n)){
            seq_gene_t gene = ml_node_val(g_n);
            if (gene.gene_id < 0)
                return err_msg(-1, 0, "mdl_mlcl_add_rna: invalid gene id: %d", gene.gene_id);
            if (gene.splice >= 3)
                return err_msg(-1, 0, "mdl_mlcl_add_rna: invalid splice: %u", gene.splice);
            int32_t f_ix = gene.gene_id + (n_genes * gene.splice);
            if (mv_push(i32, &mlcl->feat_ixs, f_ix) < 0)
                return err_msg(-1, 0, "mdl_mlcl_add_rna: failed to push feature index");
        }
    }

    ml_node_t(seq_vac_l) *v_n;
    for (v_n = ml_begin(&mol->vl); v_n; v_n = ml_node_next(v_n)){
        seq_vac_t vac = ml_node_val(v_n);
        if (vac.vix < 0)
            return err_msg(-1, 0, "mdl_mlcl_add_rna: invalid variant index: %d", vac.vix);
        // set allele to 0:ref, 1:alt, 2:any other
        uint8_t allele = vac.allele < 2 ? vac.allele : 2;
        uint32_t var_ix = vac.vix;
        uint32_t v_ix = -1;
        if (mdl_mlcl_pack_var(var_ix, allele, &v_ix) < 0)
            return err_msg(-1, 0, "mdl_mlcl_add_rna: failed to pack variant");
        if (mv_push(u32, &mlcl->var_ixs, v_ix) < 0)
            return err_msg(-1, 0, "mdl_mlcl_add_rna: failed to push variant index");
        uint8_t bqual = vac.qual;
        // cap phred at 40
        if (bqual > 93 && bqual < 0xff)
            bqual = 93;
        if (mv_push(u8, &mlcl->bquals, bqual) < 0)
            return err_msg(-1, 0, "mdl_mlcl_add_rna: failed to push base quality");
    }
    return 0;
}

int mdl_mlcl_add_atac(mdl_mlcl_t *mlcl, atac_frag_t *frag) {
    if (mlcl == NULL || frag == NULL)
        return err_msg(-1, 0, "mdl_mlcl_add_atac: argument is null");

    // add peak
    // index 0 is outside peak
    // index 1+ is peak.
    int32_t in_pk = mv_size(&frag->pks) > 0;
    int32_t pk_ix = in_pk ? mv_i(&frag->pks, 0) + 1 : 0;
    if (mv_push(i32, &mlcl->feat_ixs, pk_ix) < 0)
        return err_msg(-1, 0, "mdl_mlcl_add_atac: failed to push peak index");

    // variant(s)
    ml_node_t(seq_vac_l) *v_n;
    for (v_n = ml_begin(&frag->vl); v_n; v_n = ml_node_next(v_n)){
        seq_vac_t vac = ml_node_val(v_n);
        if (vac.vix < 0)
            return err_msg(-1, 0, "mdl_mlcl_add_atac: invalid variant index: %d", vac.vix);

        // set allele to 0:ref, 1:alt, 2:any other
        uint8_t allele = vac.allele < 2 ? vac.allele : 2;
        uint32_t var_ix = vac.vix;
        uint32_t v_ix = -1;
        if (mdl_mlcl_pack_var(var_ix, allele, &v_ix) < 0)
            return err_msg(-1, 0, "mdl_mlcl_add_atac: failed to pack variant");
        if (mv_push(u32, &mlcl->var_ixs, v_ix) < 0)
            return err_msg(-1, 0, "mdl_mlcl_add_atac: failed to push variant index");
        uint8_t bqual = vac.qual;
        // cap phred at 40
        if (bqual > 93 && bqual < 0xff)
            bqual = 93;
        if (mv_push(u8, &mlcl->bquals, bqual) < 0)
            return err_msg(-1, 0, "mdl_mlcl_add_atac: failed to push base quality");
    }
    return 0;
}

f_t mdl_bc_frip(mdl_mlcl_bc_t *mdl_bc){
    uint32_t counts = 0, pks = 0;
    kbtree_t(kb_mdl_mlcl) *bt = mdl_bc->atac;
    kbitr_t itr;
    kb_itr_first(kb_mdl_mlcl, bt, &itr); 
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, bt, &itr)){
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        if ( mv_i(&mlcl->feat_ixs, 0) > 0 ) pks += mlcl->counts;
        counts += mlcl->counts;
    }
    if (counts == 0)
        return 0.0;
    f_t ret = (f_t)pks / (f_t)counts;
    return(ret);
}

void mdl_bc_counts(mdl_mlcl_bc_t *mdl_bc, uint32_t *rna, uint32_t *atac){
    if (mdl_bc == NULL || rna == NULL || atac == NULL)
        return;
    *rna = mdl_mlcl_tot_count(mdl_bc->rna);
    *atac = mdl_mlcl_tot_count(mdl_bc->atac);
}

void mdl_bc_print(FILE *f, mdl_mlcl_bc_t *mdl_bc){

    kbtree_t(kb_mdl_mlcl) *bt;
    kbitr_t itr;

    size_t counts;

    counts = 0;
    // fprintf(f, "RNA: \n");
    bt = mdl_bc->rna;
    kb_itr_first(kb_mdl_mlcl, bt, &itr); 
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, bt, &itr)){
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        counts += mlcl->counts;
        // fprintf(f, "counts=%u\n", mlcl->counts);
        // mdl_mlcl_print(f, mlcl);
    }
    fprintf(f, "\t%zu", counts);

    counts = 0;
    // fprintf(f, "ATAC: \n");
    bt = mdl_bc->atac;
    kb_itr_first(kb_mdl_mlcl, bt, &itr); 
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, bt, &itr)){
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        counts += mlcl->counts;
        // fprintf(f, "counts=%u\n", mlcl->counts);
        // mdl_mlcl_print(f, mlcl);
    }
    fprintf(f, "\t%zu", counts);

    f_t frip = mdl_bc_frip(mdl_bc);
    fprintf(f, "\t%f\n", frip);
}

void mdl_print_count(FILE *f, mdl_t *mdl){

    size_t i;
    for (i = 0; i < mv_size(&mdl->mdl_bc_dat->bc_mlcl); ++i){
        char *bc_key = str_map_str(mdl->mdl_bc_dat->all_bcs, i);
        fprintf(f, "%s", bc_key);
        mdl_bc_print(f, &mv_i(&mdl->mdl_bc_dat->bc_mlcl, i));
    }
}

int mdl_mlcl_bc_init(mdl_mlcl_bc_t *mdl_bc){
    if (mdl_bc == NULL)
        return err_msg(-1, 0, "mdl_mlcl_bc_init: argument is null");

    mdl_bc->rna = kb_init(kb_mdl_mlcl, KB_DEFAULT_SIZE);
    if (mdl_bc->rna == NULL){
        err_msg(-1, 0, "mdl_mlcl_bc_init: failed to initialize kbtree");
        return err_msg(-1, 0, "mdl_mlcl_bc_init: %s", strerror(errno));
    }

    mdl_bc->atac = kb_init(kb_mdl_mlcl, KB_DEFAULT_SIZE);
    if (mdl_bc->atac == NULL){
        err_msg(-1, 0, "mdl_mlcl_bc_init: failed to initialize kbtree");
        return err_msg(-1, 0, "mdl_mlcl_bc_init: %s", strerror(errno));
    }

    mdl_bc->n_bc = 0;

    return 0;
}

void mdl_mlcl_bc_free(mdl_mlcl_bc_t *mdl_bc){
    if (mdl_bc == NULL) return;
    kbitr_t itr;
    kbtree_t(kb_mdl_mlcl) *bt;

    bt = mdl_bc->rna;
    if (bt != NULL) {
        kb_itr_first(kb_mdl_mlcl, bt, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, bt, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            mdl_mlcl_free(mlcl);
        }
        kb_destroy(kb_mdl_mlcl, mdl_bc->rna);
        mdl_bc->rna = NULL;
    }

    bt = mdl_bc->atac;
    if (bt != NULL) {
        kb_itr_first(kb_mdl_mlcl, bt, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, bt, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            mdl_mlcl_free(mlcl);
        }
        kb_destroy(kb_mdl_mlcl, mdl_bc->atac);
        mdl_bc->atac = NULL;
    }

    mdl_bc->n_bc = 0;
}

/*******************************************************************************
 * index structs
 ******************************************************************************/

void par_ix_init(par_ix_t *par_ix) {
    if (par_ix == NULL)
        return;
    par_ix->hs_ix = -1;
    par_ix->hd = -1;
    par_ix->s1 = -1;
    par_ix->s2 = -1;
    int i;
    for (i = 0; i < 3; ++i)
        par_ix->t_ix[i] = -1;
    par_ix->t_n = -1;
}

void hs_ix_init(hs_ix_t *hs_ix) {
    if (hs_ix == NULL)
        return;
    hs_ix->n_sam = 0;
    hs_ix->n_hs = 0;
    hs_ix->ix2pars = NULL;
}

hs_ix_t *hs_ix_alloc() {
    hs_ix_t *hs_ix = calloc(1, sizeof(hs_ix_t));
    if (hs_ix == NULL) {
        err_msg(-1, 0, "hs_ix_alloc: %s", strerror(errno));
        return NULL;
    }
    hs_ix_init(hs_ix);
    return hs_ix;
}

void hs_ix_free(hs_ix_t *hs_ix) {
    if (hs_ix == NULL)
        return;
    hs_ix->n_sam = 0;
    hs_ix->n_hs = 0;
    free(hs_ix->ix2pars);
    hs_ix->ix2pars = NULL;
}

void hs_ix_dstry(hs_ix_t *hs_ix) {
    if (hs_ix == NULL)
        return;
    hs_ix_free(hs_ix);
    free(hs_ix);
}

void hs_ix_set(hs_ix_t *hs_ix, uint16_t n_samples) {
    if (hs_ix == NULL)
        return;

    // indices
    int n_par = 3;
    hs_ix->n_sam = n_samples;
    hs_ix->n_hs = 1 + n_samples + (n_samples * (n_samples - 1) / 2);
    hs_ix->ix2pars = calloc(n_par * hs_ix->n_hs, sizeof(int));

    int s1 = 0, s2 = 1;
    uint32_t ixi = 0;
    hs_ix->ix2pars[CMI(0, ixi, n_par)] = 0;
    hs_ix->ix2pars[CMI(1, ixi, n_par)] = -1;
    hs_ix->ix2pars[CMI(2, ixi, n_par)] = -1;
    ++ixi;

    for (s1 = 0; s1 < n_samples; ++s1){
        hs_ix->ix2pars[CMI(0, ixi, n_par)] = 1;
        hs_ix->ix2pars[CMI(1, ixi, n_par)] = s1;
        hs_ix->ix2pars[CMI(2, ixi, n_par)] = -1;
        ++ixi;
    }
    for (s1 = 0; s1 < n_samples - 1; ++s1){
        for (s2 = s1 + 1; s2 < n_samples; ++s2){
            hs_ix->ix2pars[CMI(0, ixi, n_par)] = 2;
            hs_ix->ix2pars[CMI(1, ixi, n_par)] = s1;
            hs_ix->ix2pars[CMI(2, ixi, n_par)] = s2;
            ++ixi;
        }
    }
    if (ixi != hs_ix->n_hs)
        err_msg(0, 1, "hs_ix_set: likely bug: ixi=%u != n_hs=%u", ixi, hs_ix->n_hs);
}

int hs_ix_get_pars(hs_ix_t *hs_ix, unsigned ix,
        par_ix_t *par_ix) {

    int n_par = 3;
    if (ix >= hs_ix->n_hs)
        return err_msg(-1, 0, "hs_ix_get_pars: ix=%i > n_ixs=%u", ix, hs_ix->n_hs);

    // -1 for if NA/invalid
    par_ix->hs_ix = ix;
    uint32_t divby = hs_ix->n_hs - 1;
    if (ix == 0)
        par_ix->hs_ix = ix;
    else
        par_ix->hs_ix = ((ix - 1) % divby) + 1;

    if ((uint32_t)par_ix->hs_ix >= hs_ix->n_hs)
        return err_msg(-1, 0, "hs_ix_get_pars: par_ix->hs_ix=%i >= n_ixs=%u", par_ix->hs_ix, hs_ix->n_hs);

    par_ix->hd = hs_ix->ix2pars[CMI(0, ix, n_par)];
    par_ix->s1 = hs_ix->ix2pars[CMI(1, ix, n_par)];
    par_ix->s2 = hs_ix->ix2pars[CMI(2, ix, n_par)];

    switch (par_ix->hd) {
        case 0:
            par_ix->t_ix[0] = hs_ix->n_sam;
            par_ix->t_ix[1] = -1;
            par_ix->t_ix[2] = -1;
            par_ix->t_n = 1;
            break;
        case 1:
            par_ix->t_ix[0] = hs_ix->n_sam;
            par_ix->t_ix[1] = par_ix->s1;
            par_ix->t_ix[2] = -1;
            par_ix->t_n = 2;
            break;
        case 2:
            par_ix->t_ix[0] = hs_ix->n_sam;
            par_ix->t_ix[1] = par_ix->s1;
            par_ix->t_ix[2] = par_ix->s2;
            par_ix->t_n = 3;
            break;
        default:
            return err_msg(-1, 0, "hs_ix_get_pars: "
                    "hd=%i is invalid, there is a bug", par_ix->hd);
    }
    return(0);
}

/*******************************************************************************
 * mdl_pars_t
 ******************************************************************************/

int mdl_pars_init(mdl_pars_t *mp){
    if (mp == NULL)
        return err_msg(-1, 0, "mdl_pars_init: argument is null");

    mp->D = 0;
    mp->T = 0;
    mp->V = 0;
    mp->M = 0;
    mp->n_hs = 0;
    
    mp->pi = NULL;
    mp->pi_amb = NULL;
    mp->alpha_rna1 = NULL;
    mp->alpha_atac1 = NULL;
    mp->alpha_rna_se = NULL;
    mp->alpha_atac_se = NULL;
    mp->alpha_rna_info = NULL;
    mp->alpha_atac_info = NULL;
    mp->gamma = NULL;
    mp->_pi_sum = NULL;
    mp->max_par_delta = 0;

    int i;
    for (i = 0; i < 3; ++i) {
        mp->lambda[i] = 0;
        mp->_lambda_sum[i] = 0;
    }

    int err;
    if ((err = pthread_mutex_init(&mp->sum_lock, NULL)) != 0)
        return err_msg(-1, 0, "mdl_pars_init: failed to initialize mutex: %i", err);
    return 0;
}

mdl_pars_t *mdl_pars_alloc() {
    mdl_pars_t *mp = calloc(1, sizeof(mdl_pars_t));
    if (mp == NULL) {
        err_msg(-1, 0, "mdl_pars_alloc: %s", strerror(errno));
        return NULL;
    }
    if (mdl_pars_init(mp) < 0) {
        free(mp);
        return NULL;
    }

    return(mp);
}

void mdl_pars_free(mdl_pars_t *mp) {
    if (mp == NULL) return;


    free(mp->pi);
    free(mp->pi_amb);
    free(mp->alpha_rna1);
    free(mp->alpha_atac1);
    free(mp->alpha_rna_se);
    free(mp->alpha_atac_se);
    free(mp->alpha_rna_info);
    free(mp->alpha_atac_info);
    free(mp->gamma);
    free(mp->_pi_sum);
    
    mp->D = 0;
    mp->T = 0;
    mp->V = 0;
    mp->M = 0;
    mp->n_hs = 0;
    
    mp->pi = NULL;
    mp->pi_amb = NULL;
    mp->alpha_rna1 = NULL;
    mp->alpha_rna1 = NULL;
    mp->alpha_atac1 = NULL;
    mp->alpha_rna_se = NULL;
    mp->alpha_atac_se = NULL;
    mp->alpha_rna_info = NULL;
    mp->alpha_atac_info = NULL;
    mp->gamma = NULL;
    mp->_pi_sum = NULL;
    mp->max_par_delta = 0;

    int i;
    for (i = 0; i < 3; ++i) {
        mp->lambda[i] = 0;
        mp->_lambda_sum[i] = 0;
    }

    pthread_mutex_destroy(&mp->sum_lock);
}

void mdl_pars_dstry(mdl_pars_t *mp){
    if (mp == NULL) return;
    mdl_pars_free(mp);
    free(mp);
}

int mld_pars_set_num_alloc(mdl_pars_t *mp, uint32_t D, uint32_t T,
        uint32_t V, uint16_t M) {
    if (mp == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: argument is null");

    mp->D = D;
    mp->T = T;
    mp->V = V;
    mp->M = M;

    mp->n_hs = 1 + (mp->M) + (mp->M * (mp->M - 1) / 2);

    // allocate parameter fields
    mp->pi = calloc(mp->M, sizeof(f_t));
    if (mp->pi == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->pi_amb = calloc(mp->M, sizeof(f_t));
    if (mp->pi_amb == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->alpha_rna1 = calloc(D * mp->n_hs, sizeof(f_t));
    if (mp->alpha_rna1 == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->alpha_atac1 = calloc(D * mp->n_hs, sizeof(f_t));
    if (mp->alpha_atac1 == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->alpha_rna_se = calloc(D * mp->n_hs, sizeof(f_t));
    if (mp->alpha_rna_se == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));
    mp->alpha_atac_se = calloc(D * mp->n_hs, sizeof(f_t));
    if (mp->alpha_atac_se == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->alpha_rna_info = calloc(D * mp->n_hs, sizeof(f_t));
    if (mp->alpha_rna_info == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));
    mp->alpha_atac_info = calloc(D * mp->n_hs, sizeof(f_t));
    if (mp->alpha_atac_info == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->gamma = calloc((mp->M + 1) * mp->V, sizeof(f_t));
    if (mp->gamma == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    // allocate sum fields
    mp->_pi_sum = calloc(mp->M, sizeof(f_t));
    if (mp->_pi_sum == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->max_par_delta = 0;

    return 0;
}

int mdl_pars_reset_sums(mdl_pars_t *mp, f_t psc) {
    if (mp == NULL)
        return err_msg(-1, 0, "mdl_pars_reset_sums: argument is null");

    if (mp->D == 0)
        return 0;

    uint32_t i;
    for (i = 0; i < 3; ++i)
        mp->_lambda_sum[i] = psc;

    for (i = 0; i < mp->M; ++i)
        mp->_pi_sum[i] = psc;

    mp->max_par_delta = 0;

    return 0;
}

int mdl_pars_pi_fix(mdl_pars_t *mp){
    if (mp == NULL)
        return err_msg(-1, 0, "mdl_pars_pi_fix: argument is null");
    if (mp->pi == NULL)
        return err_msg(-1, 0, "mdl_pars_pi_fix: pi is null");

    uint32_t i, j, M = mp->M;

    // set pi uniform
    f_t pi_upd = 1.0 / M;
    for (i = 0; i < M; ++i)
        mp->pi[i] = pi_upd;

    f_t pi_sum = 0;
    for (i = 0; i < M; ++i){
        for (j = i+1; j < M; ++j){
            pi_sum += (mp->pi[i] * mp->pi[j]);
        }
    }
    mp->pi_d_sum = pi_sum;

    return(0);
}

int mdl_pars_pi_amb_fix(mdl_pars_t *mp) {
    if (mp == NULL)
        return err_msg(-1, 0, "mdl_pars_pi_amb_fix: argument is null");
    if (mp->pi_amb == NULL)
        return err_msg(-1, 0, "mdl_pars_pi_amb_fix: pi_amb is null");

    uint32_t i, M = mp->M;

    // set pi_amb uniform
    f_t pi_upd = 1.0 / M;
    for (i = 0; i < M; ++i)
        mp->pi_amb[i] = pi_upd;

    return(0);
}

int mdl_pars_add_gamma(mdl_pars_t *mp, float **a, int nv, int ns){
    if (mp == NULL || a == NULL)
        return err_msg(-1, 0, "mdl_pars_add_gamma: a is NULL");

    if (mp->M != (uint16_t)ns)
        return err_msg(-1, 0, "mdl_pars_add_gamma: "
                "ns and number of samples in a don't match");

    if (nv < 1 || ns < 1)
        return err_msg(-1, 0, "mdl_pars_add_gamma: nv=%i and ns=%i "
                "must be > 0", nv, ns);

    uint16_t ns1 = (uint16_t)ns + 1;

    free(mp->gamma);
    mp->gamma = calloc(nv * ns1, sizeof(f_t));
    if (mp->gamma == NULL)
        return err_msg(-1, 0, "mdl_pars_add_gamma: %s", strerror(errno));
    
    int v;
    for (v = 0; v < nv; ++v){
        int n_miss = 0;
        f_t d_sum = 0.0;
        int s;
        for (s = 0; s < ns; ++s){
            // if variant is missing, set to -1
            f_t av = a[v] == NULL ? -1 : a[v][s];
            mp->gamma[CMI(s, v, ns1)] = av;
            if (av >= 0)
                d_sum += av;
            else
                n_miss++;
        }
        f_t ambp;
        if (ns == n_miss){
            ambp = -1;
        } else {
            ambp = d_sum / (f_t)(ns - n_miss);
        }
        mp->gamma[CMI(s, v, ns1)] = ambp;
    }

    return 0;
}

int mdl_pars_set_gamma_amb(mdl_pars_t *mp){
    if (mp == NULL)
        return err_msg(-1, 0, "mdl_pars_set_gamma_amb: argument is null");
    if (mp->pi_amb == NULL)
        return err_msg(-1, 0, "mdl_pars_set_gamma_amb: pi_amb is null");
    if (mp->gamma == NULL)
        return err_msg(-1, 0, "mdl_pars_set_gamma_amb: gamma is null");

    uint32_t i, M = mp->M, V = mp->V, M1 = M + 1;

    // ensure that pi_amb is in simplex
    if (mdl_check_pi_amb(mp) < 0)
        return err_msg(-1, 0, "mdl_pars_set_gamma_amb: pi_amb is invalid");

    // update gamma
    for (i = 0; i < V; ++i){
        uint16_t st, n_miss = 0;
        f_t ap = 0;
        f_t pi_tot = 0;
        for (st = 0; st < M; ++st){
            f_t gprob = mp->gamma[CMI(st, i, M1)];
            if (gprob < 0) {
                ++n_miss; // if missing, skip
            } else {
                ap += mp->pi_amb[st] * gprob;
                pi_tot += mp->pi_amb[st];
            }
        }
        if (n_miss == M || pi_tot <= 0) {
            mp->gamma[CMI(st, i, M1)] = -1; // if all missing, set amb to missing
        } else {
            mp->gamma[CMI(st, i, M1)] = ap / pi_tot;
            // ensure valid probabilities
            if (mp->gamma[CMI(st, i, M1)] > 1)
                mp->gamma[CMI(st, i, M1)] = 1;
        }
    }
    return(0);
}

int mdl_pars_set_dat(mdl_pars_t *mp, mdl_bc_dat_t *bd, obj_pars *objs,
        uint16_t n_samples){
    if (mp == NULL || bd == NULL || objs == NULL)
        return err_msg(-1, 0, "mdl_pars_set_dat: argument is null");
    if (objs->gv == NULL)
        return err_msg(-1, 0, "mdl_pars_set_dat: no genotype data present");

    if (bd->all_bcs->n < 1)
        return err_msg(-1, 0, "mdl_pars_set_dat: no barcodes in bc_dat");
    if (bd->test_bcs->n < 1)
        return err_msg(-1, 0, "mdl_pars_set_dat: no test barcodes in bc_dat");

    // rerror if no variants are present
    if (bd->V < 1)
        return err_msg(-1, 0, "mdl_pars_set_dat: no variants found");

    // set numbers and allocate arrays
    if (mld_pars_set_num_alloc(mp, bd->all_bcs->n, bd->test_bcs->n,
                bd->V, n_samples) < 0)
        return -1;

    // set pi and pi_amb uniform
    if (mdl_pars_pi_fix(mp) < 0)
        return(-1);
    if (mdl_pars_pi_amb_fix(mp) < 0)
        return(-1);

    // prior of 1 for lambda
    f_t n_amb_bcs = 1;
    f_t n_nuc_bcs = 1;

    uint32_t i, j;
    for (i = 0; i < mp->D; ++i){
        int fl;

        // if no reads
        fl = bflg_get(&bd->absent_bc, i);
        if (fl == 1) continue;

        // if ambient
        fl = bflg_get(&bd->amb_flag, i);

        // bam data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, i);

        if (fl)
            n_amb_bcs += bc_dat.n_bc;
        else
            n_nuc_bcs += bc_dat.n_bc;

        // set alpha
        for (j = 0; j < mp->n_hs; ++j) {
            uint32_t aix = CMI(i, j, mp->D);
            if (fl) {
                mp->alpha_rna1[aix] = 1.0;
                mp->alpha_atac1[aix] = 1.0;
            } else {
                mp->alpha_rna1[aix] = NAN;
                mp->alpha_atac1[aix] = NAN;
            }
            mp->alpha_rna_se[aix] = NAN;
            mp->alpha_atac_se[aix] = NAN;

            mp->alpha_rna_info[aix] = NAN;
            mp->alpha_atac_info[aix] = NAN;
        }
    }

    // set lambda
    f_t n_tot_bcs = n_amb_bcs + n_nuc_bcs;
    f_t n_sng = n_nuc_bcs * 0.9;
    f_t n_dbl = n_nuc_bcs * 0.1;
    mp->lambda[0] = n_amb_bcs / n_tot_bcs;
    mp->lambda[1] = n_sng / n_tot_bcs;
    mp->lambda[2] = n_dbl / n_tot_bcs;

    // create array of allele probabilities from VCF
    int32_t *vixs = calloc(bd->V, sizeof(int32_t));
    if (vixs == NULL)
        return err_msg(-1, 0, "mdl_pars_set_dat: %s", strerror(errno));
    for (i = 0; i < bd->V; ++i) vixs[i] = i;
    float **gm = ap_array_gt(objs->gv, objs->vcf_hdr, vixs, bd->V, "GT");
    free(vixs);
    if (gm == NULL)
        return(-1);

    // set gamma from allele probabilities array
    if (mdl_pars_add_gamma(mp, gm, bd->V, n_samples) < 0) {
        for (i = 0; i < bd->V; ++i){
            free(gm[i]);
        }
        free(gm);
        return -1;
    }

    // set tau
    mp->tau = TAU;
    if (prob_invalid(mp->tau))
        return err_msg(-1, 0, "mdl_pars_set_dat: invalid prob tau = %e", mp->tau);

    // alpha prior
    mp->alpha_prior_w = objs->alpha_prior_w;
    mp->alpha_prior_rna = 0.1;
    mp->alpha_prior_atac = 0.1;

    // free
    for (i = 0; i < bd->V; ++i){
        free(gm[i]);
    }
    free(gm);
    
    return(0);
}

int is_in_simplex(f_t *x, size_t len) {
    if (x == NULL)
        return err_msg(-1, 0, "is_in_simplex: argument is null");

    f_t lt1 = 1 - 1e-2, ut1 = 1 + 1e-2;
    f_t tot = 0.0;
    unsigned int i;
    for (i = 0; i < len; ++i) {
        if (prob_invalid(x[i]))
            return err_msg(-1, 0, "is_in_simplex: x[%i] = %e", i, x[i]);
        tot += x[i];
    }
    if (psum_invalid(tot, lt1, ut1))
        return err_msg(-1, 0, "is_in_simplex: x sum = %e", tot);

    return 0;
}

int mdl_check_lambda(mdl_pars_t *mp) {
    if (mp == NULL)
        return err_msg(-1, 0, "mdl_check_lambda: argument is null");

    if (is_in_simplex(mp->lambda, 3) < 0)
        return err_msg(-1, 0, "mdl_check_lambda: lambda is invalid");

    return 0;
}

int mdl_check_pi(mdl_pars_t *mp) {
    if (mp == NULL)
        return err_msg(-1, 0, "mdl_check_pi: argument is null");

    if (mp->pi == NULL)
        return err_msg(-1, 0, "mdl_check_pi: pi is null");

    if (is_in_simplex(mp->pi, mp->M) < 0)
        return err_msg(-1, 0, "mdl_check_pi: pi is invalid");

    return 0;
}

int mdl_check_pi_amb(mdl_pars_t *mp) {
    if (mp == NULL)
        return err_msg(-1, 0, "mdl_check_pi_amb: argument is null");

    if (mp->pi_amb == NULL)
        return err_msg(-1, 0, "mdl_check_pi_amb: pi_amb is null");

    // check pi amb
    if (is_in_simplex(mp->pi_amb, mp->M) < 0)
        return err_msg(-1, 0, "mdl_check_pi_amb: pi_amb is invalid");

    return 0;
}

int mdl_pars_check(mdl_pars_t *mp){

    uint32_t i, j, M1 = mp->M + 1;

    // check lambda
    if (mdl_check_lambda(mp) < 0)
        return -1;

    // check pi
    if (mdl_check_pi(mp) < 0)
        return -1;

    // check pi_amb
    if (mdl_check_pi_amb(mp) < 0)
        return -1;

    // alpha
    for (i = 0; i < mp->D; ++i){
        for (j = 0; j < mp->n_hs; ++j) {
            f_t a1;
            a1 = mp->alpha_rna1[CMI(i, j, mp->D)];
            if (!isnan(a1) && (a1 < 0 || a1 > 1))
                return err_msg(-1, 0, "mdl_pars_check: alpha_rna1[%i,%i] = %e", i, j, a1);
            a1 = mp->alpha_atac1[CMI(i, j, mp->D)];
            if (!isnan(a1) && (a1 < 0 || a1 > 1))
                return err_msg(-1, 0, "mdl_pars_check: alpha_atac1[%i,%i] = %e", i, j, a1);
        }
    }

    /*
    for (i = 0; i < mp->D; ++i){
        for (j = 0; j < mp->n_hs; ++j) {
            f_t a1, a2;

            a1 = mp->alpha_rna1[CMI(i, j, mp->D)];
            if (num_invalid(a1) || a1 < 0)
                return err_msg(-1, 0, "mdl_pars_check: alpha_rna1[%i,%i] = %f", i, j, a1);
            if (num_invalid(a1 + a2) || a1 + a2 < 0)
                return err_msg(-1, 0, "mdl_pars_check: alpha_rna12[%i,%i] = %f", i, j, a1 + a2);

            a1 = mp->alpha_atac1[CMI(i, j, mp->D)];
            if (num_invalid(a1) || a1 < 0)
                return err_msg(-1, 0, "mdl_pars_check: alpha_atac1[%i,%i] = %f", i, j, a1);
            if (num_invalid(a1 + a2) || a1 + a2 < 0)
                return err_msg(-1, 0, "mdl_pars_check: alpha_atac12[%i,%i] = %f", i, j, a1 + a2);
            a1 = mp->alpha_rna_se[CMI(i, j, mp->D)];
            if (a1 < 0)
                return err_msg(-1, 0, "mdl_pars_check: alpha_rna_se[%i,%i] = %f", i, j, a1);
            a1 = mp->alpha_atac_se[CMI(i, j, mp->D)];
            if (a1 < 0)
                return err_msg(-1, 0, "mdl_pars_check: alpha_atac_se[%i,%i] = %f", i, j, a1);
        }
    }
    */

    // gamma
    for (i = 0; i < M1; ++i) {
        for (j = 0; j < mp->V; ++j) {
            f_t p = mp->gamma[CMI(i, j, M1)];
            // if missing, don't check
            if (p > (-1 - 1e-2) && p < (-1 + 1e-2))
                continue;
            if (prob_invalid(p))
                return err_msg(-1, 0, "mdl_pars_check: gamma[%i, %i] = %e",
                        i, j, p);
        }
    }

    return(0);
}

/*******************************************************************************
 * mdl_bc_dat_t
 ******************************************************************************/

int mdl_bc_dat_init(mdl_bc_dat_t *mdl_bc_dat) {
    if (mdl_bc_dat == NULL)
        return err_msg(-1, 0, "mdl_bc_dat_init: argument is null");
    mdl_bc_dat->all_bcs = init_str_map();
    mdl_bc_dat->test_bcs = init_str_map();

    bflg_init_empty(&mdl_bc_dat->amb_flag);

    bflg_init_empty(&mdl_bc_dat->absent_bc);

    mdl_bc_dat->V = 0;
    mv_init(&mdl_bc_dat->bc_mlcl);
    return 0;
}

mdl_bc_dat_t *mdl_bc_dat_alloc() {
    mdl_bc_dat_t *mdl_bc_dat = calloc(1, sizeof(mdl_bc_dat_t));
    if (mdl_bc_dat == NULL) {
        err_msg(-1, 0, "mdl_bc_dat_alloc: %s", strerror(errno));
        return NULL;
    }

    if (mdl_bc_dat_init(mdl_bc_dat) < 0) {
        free(mdl_bc_dat);
        return NULL;
    }
    return mdl_bc_dat;
}

void mdl_bc_dat_free(mdl_bc_dat_t *mdl_bc_dat) {
    if (mdl_bc_dat == NULL)
        return;

    uint32_t i;
    for (i = 0; i < mv_size(&mdl_bc_dat->bc_mlcl); ++i)
        mdl_mlcl_bc_free(&mv_i(&mdl_bc_dat->bc_mlcl, i));

    destroy_str_map(mdl_bc_dat->all_bcs);
    destroy_str_map(mdl_bc_dat->test_bcs);

    bflg_free(&mdl_bc_dat->amb_flag);
    bflg_free(&mdl_bc_dat->absent_bc);

    mdl_bc_dat->V = 0;

    mv_free(&mdl_bc_dat->bc_mlcl);
}

void mdl_bc_dat_dstry(mdl_bc_dat_t *mdl_bc_dat) {
    if (mdl_bc_dat == NULL)
        return;

    mdl_bc_dat_free(mdl_bc_dat);

    free(mdl_bc_dat);
}

int mdl_bc_dat_add_bc(mdl_bc_dat_t *mdl_bc_dat, char *bc_key, int absent, int ambient,
        int dup_ok, int *found){
    if (mdl_bc_dat == NULL || bc_key == NULL)
        return err_msg(-1, 0, "mdl_bc_dat_add_bc: argument is null");
    int bc_ix;
    if ( (bc_ix = add2str_map(mdl_bc_dat->all_bcs, bc_key, found)) < 0)
        return(-1);
    if (bc_ix < 0)
        return(-1);

    if (!dup_ok && *found)
        return err_msg(-1, 0, "mdl_bc_dat_add_bc: barcode '%s' already added", bc_key);

    // add bc_mlcl element
    if (*found) {
        mv_i(&mdl_bc_dat->bc_mlcl, bc_ix).n_bc += 1;
        return 0;
    } else {
        mdl_mlcl_bc_t mdl_mlcl;
        mdl_mlcl_bc_init(&mdl_mlcl);
        mdl_mlcl.n_bc = 1;
        if (mv_insert(mdl_bc_v, &mdl_bc_dat->bc_mlcl, mdl_mlcl, bc_ix) < 0)
            return(-1);
    }

    if (bflg_resize(&mdl_bc_dat->absent_bc, mdl_bc_dat->all_bcs->n) < 0)
        return(-1);
    if (bflg_resize(&mdl_bc_dat->amb_flag, mdl_bc_dat->all_bcs->n) < 0)
        return(-1);

    // set to absent if no reads/empty
    if (absent) {
        bflg_set(&mdl_bc_dat->absent_bc, bc_ix);
    } else {
        bflg_unset(&mdl_bc_dat->absent_bc, bc_ix);
    }

    // set ambient flag and add to test_bcs
    if (ambient) {
        bflg_set(&mdl_bc_dat->amb_flag, bc_ix);
    } else {
        bflg_unset(&mdl_bc_dat->amb_flag, bc_ix);
    }
    
    // add to test bcs if not ambient and not empty
    if (!ambient && !absent) {
        if (add2str_map(mdl_bc_dat->test_bcs, bc_key, found) < 0)
            return -1;
    }
    
    return(bc_ix);
}

int mdl_bc_dat_bam_data(mdl_bc_dat_t *mdl_bc_dat, bam_data_t *bam_data, obj_pars *objs){
    if (mdl_bc_dat == NULL || bam_data == NULL || objs == NULL)
        return err_msg(-1, 0, "mdl_bc_dat_bam_data: argument is null");

    khint_t k_bc;

    // not modeling genes or peaks currently
    uint32_t n_genes = 0;
    uint32_t n_vars = 0;
    if (objs->gv)
        n_vars = mv_size(&objs->gv->vix2var);
    mdl_bc_dat->V = n_vars;
    if (mdl_bc_dat->V < 1)
        return err_msg(-1, 0, "mdl_bc_dat_bam_data: no variants found");

    uint32_t c_thresh = objs->out_min;
    int bc_ix, found = 0;

    // add the empty droplets barcode
    char bc_empty[] = "Empty";
    bc_ix = mdl_bc_dat_add_bc(mdl_bc_dat, bc_empty, 0, 1, 1, &found);
    if (bc_ix < 0)
        return(-1);

    khash_t(kh_bc_dat) *bam_bcs = bam_data->bc_data;

    for (k_bc = kh_begin(bam_bcs); k_bc != kh_end(bam_bcs); ++k_bc){
        if (!kh_exist(bam_bcs, k_bc)) continue;

        char *bc_key = kh_key(bam_bcs, k_bc);

        bc_data_t *bam_bc = kh_val(bam_bcs, k_bc);
        if (bam_bc == NULL)
            return err_msg(-1, 0, "mdl_bc_dat_bam_data: barcode data is NULL");

        bc_stats_t *bc_stat = bam_bc->bc_stats;
        if (bc_stat == NULL)
            return err_msg(-1, 0, "mdl_bc_dat_bam_data: barcode stats is NULL");

        // barcode can be both absent and empty
        // low_count=1 when num. atac and num. rna are both below threshold
        int low_count = (bc_stat->atac_counts < c_thresh) && (bc_stat->rna_counts < c_thresh);
        int dup_ok = low_count; // only low count empty barcode can be duplicate

        // if droplet is low count, set barcode to empty (also sets empty)
        if (low_count)
            bc_key = bc_empty;

        found = 0;
        bc_ix = mdl_bc_dat_add_bc(mdl_bc_dat, bc_key, 0, low_count, dup_ok, &found);
        if (bc_ix < 0)
            return(-1);

        // add reads to mdl_mlcl
        mdl_mlcl_bc_t *mdl_bc = &mv_i(&mdl_bc_dat->bc_mlcl, bc_ix);

        // loop through RNA
        ml_node_t(ml_rm) *rna_node = ml_begin(bam_bc->rna_mols);
        for (; rna_node != NULL; rna_node = ml_node_next(rna_node)) {
            rna_mol_t *mol = &ml_node_val(rna_node);

            // skip RNA molecules without variants
            if (ml_size(&mol->vl) == 0)
                continue;

            // filter intra/inter gene reads
            uint32_t n_feats = ml_size(&mol->gl);
            if (n_feats == 0 && objs->mdl_inter_reads == 0) {
                continue;
            } else if (n_feats > 0 && objs->mdl_intra_reads == 0) {
                continue;
            }

            // initialize mlcl to 0
            mdl_mlcl_t mlcl;
            mdl_mlcl_init(&mlcl);
            mlcl.counts = 1;

            if (mdl_mlcl_add_rna(&mlcl, mol, n_genes) < 0) {
                mdl_mlcl_free(&mlcl);
                return -1;
            }

            mdl_mlcl_t *p; // pointer for btree add
            p = kb_getp(kb_mdl_mlcl, mdl_bc->rna, &mlcl);
            if (!p) {
                kb_putp(kb_mdl_mlcl, mdl_bc->rna, &mlcl);
            } else {
                p->counts += 1;
                mdl_mlcl_free(&mlcl);
            }
        }

        // loop through atac
        ml_node_t(ml_af) *atac_node = ml_begin(bam_bc->atac_frgs);
        for (; atac_node != NULL; atac_node = ml_node_next(atac_node)) {
            atac_frag_t *frag = &ml_node_val(atac_node);

            // skip ATAC read if no variants
            if (ml_size(&frag->vl) == 0)
                continue;

            // filter intra/inter peak reads
            uint32_t n_feats = mv_size(&frag->pks);
            if (n_feats == 0 && objs->mdl_inter_reads == 0) {
                continue;
            } else if (n_feats > 0 && objs->mdl_intra_reads == 0) {
                continue;
            }

            // initialize mlcl to all 0
            mdl_mlcl_t mlcl;
            mdl_mlcl_init(&mlcl);
            mlcl.counts = 1;

            if (mdl_mlcl_add_atac(&mlcl, frag) < 0) {
                mdl_mlcl_free(&mlcl);
                return -1;
            }

            mdl_mlcl_t *p;
            p = kb_getp(kb_mdl_mlcl, mdl_bc->atac, &mlcl);
            if (!p) {
                kb_putp(kb_mdl_mlcl, mdl_bc->atac, &mlcl);
            } else {
                p->counts += 1;
                mdl_mlcl_free(&mlcl);
            }
        }

        // note the read data in `bam_data' is freed from the barcode here
        bc_data_free_reads(bam_bc);
    }

    if ((uint32_t)mdl_bc_dat->all_bcs->n != mv_size(&mdl_bc_dat->bc_mlcl))
        return err_msg(-1, 0, "mdl_bc_dat_bam_data: number of barcodes in bc_dat "
                       "and bc_mlcl do not match");

    return(0);
}

/*******************************************************************************
 * mdl_t
 ******************************************************************************/

mdl_t *mdl_alloc(){
    mdl_t *mdl = (mdl_t *)calloc(1, sizeof(mdl_t));
    if (mdl == NULL){
        err_msg(-1, 0, "mdl_alloc: %s", strerror(errno));
        return(NULL);
    }

    mdl->samples = init_str_map();
    if (mdl->samples == NULL) {
        free(mdl);
        return NULL;
    }

    hs_ix_init(mdl->hs_ix);
    mdl->hs_ix = hs_ix_alloc();
    if (mdl->hs_ix == NULL) {
        free(mdl);
        return NULL;
    }

    mdl->n_hs = -1;

    mdl->sub_lp_x = NULL;
    mdl->sub_lp_hs = NULL;
    mdl->sub_lp_h = NULL;
    mdl->sub_best_hs = NULL;
    mdl->sub_best_h1 = NULL;
    mdl->sub_best_h2 = NULL;
    mdl->sub_pp_h = NULL;

    mdl->alpha_eps = 1e-6;
    mdl->eps = 1e-6;
    mdl->max_iter = 100;

    mdl->has_rna = 0;
    mdl->has_atac = 0;

    mdl->threads = 1;

    mdl->mdl_bc_dat = mdl_bc_dat_alloc();
    if (mdl->mdl_bc_dat == NULL) {
        hs_ix_dstry(mdl->hs_ix);
        free(mdl);
        return NULL;
    }

    mdl->mp = mdl_pars_alloc();
    if (mdl->mp == NULL) {
        mdl_bc_dat_dstry(mdl->mdl_bc_dat);
        hs_ix_dstry(mdl->hs_ix);
        free(mdl);
        return(NULL);
    }

    return(mdl);
}

void mdl_dstry(mdl_t *m){
    if (m == NULL) return;

    destroy_str_map(m->samples);
    hs_ix_dstry(m->hs_ix);

    free(m->sub_lp_hs);
    free(m->sub_lp_h);
    free(m->sub_lp_x);
    free(m->sub_best_hs);
    free(m->sub_best_h1);
    free(m->sub_best_h2);
    free(m->sub_pp_h);

    mdl_bc_dat_dstry(m->mdl_bc_dat);
    mdl_pars_dstry(m->mp);
    free(m);
}

int mdl_set_samples(mdl_t *mdl, bcf_hdr_t *vcf_hdr){
    if (mdl == NULL || vcf_hdr == NULL)
        return err_msg(-1, 0, "mdl_set_samples: arguments are NULL");

    int n_samples = bcf_hdr_nsamples(vcf_hdr);
    destroy_str_map(mdl->samples);
    mdl->samples = init_str_map_array(vcf_hdr->samples, n_samples);
    if (mdl->samples == NULL)
        return err_msg(-1, 0, "mdl_set_samples: failed to copy samples from vcf header");


    return(0);
}

int mdl_set_hs_ix(mdl_t *mdl) {
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_set_hs_ix: argument is null");

    if (mdl->samples == NULL)
        return err_msg(-1, 0, "mdl_set_hs_ix: run mdl_set_samples before");
    if (mdl->samples->n < 1)
        return err_msg(-1, 0, "mdl_set_hs_ix: run mdl_set_samples before");

    uint16_t n_sam = mdl->samples->n;
    hs_ix_set(mdl->hs_ix, n_sam);
    return 0;
}

int mdl_set_rna_atac(mdl_t *mdl, uint8_t has_rna, uint8_t has_atac) {
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_set_rna_atac: argument is null");
    mdl->has_rna = has_rna;
    mdl->has_atac = has_atac;
    return 0;
}

int mdl_alloc_probs(mdl_t *mdl) {
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_alloc_probs: argument is null");

    if (mdl->mdl_bc_dat == NULL || mdl->hs_ix == NULL)
        return err_msg(-1, 0, "mdl_alloc_probs: 'mdl_bc_dat' and 'hs_ix' "
                "must be initialized");
    if (mdl->mdl_bc_dat->all_bcs->n < 1)
        return err_msg(-1, 0, "mdl_alloc_probs: no barcodes, cannot alloc probs");

    if (mdl->hs_ix->n_hs < 1)
        return err_msg(-1, 0, "mdl_alloc_probs: indices must be set with "
                "'mdl_set_hs_ix'");

    uint32_t D = mdl->mdl_bc_dat->all_bcs->n;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    mdl->n_hs = n_hs;

    mdl->sub_lp_hs = calloc(D * n_hs, sizeof(f_t));
    mdl->sub_lp_h = calloc(D * 3, sizeof(f_t));
    mdl->sub_lp_x = calloc(D, sizeof(f_t));
    mdl->sub_best_hs = calloc(D, sizeof(int32_t));
    mdl->sub_best_h1 = calloc(D, sizeof(int32_t));
    mdl->sub_best_h2 = calloc(D, sizeof(int32_t));
    mdl->sub_pp_h = calloc(D * 3, sizeof(f_t));

    if (mdl->sub_pp_h == NULL)
        return err_msg(-1, 0, "mdl_alloc_probs: %s", strerror(errno));

    return 0;
}

/*******************************************************************************
 * Probability functions
 ******************************************************************************/
// all functions output p(x), not log p(x)

void pr_hd(mdl_pars_t *mp, par_ix_t *par_ix, f_t *prob){
    *prob = mp->lambda[par_ix->hd];
}

void pr_sd(mdl_pars_t *mp, par_ix_t *par_ix, f_t *prob) {
    f_t p_sd;
    f_t *pi = mp->pi;
    int s1 = par_ix->s1, s2 = par_ix->s2;
    f_t off_diag = mp->pi_d_sum;
    switch (par_ix->hd) {
        case 0:
            p_sd = 1;
            break;
        case 1:
            p_sd = pi[s1];
            break;
        case 2:
            p_sd = (pi[s1] * pi[s2]) / off_diag;
            break;
        default:
            p_sd = -1;
    }
    *prob = p_sd;
}

int pr_tdm(mdl_pars_t *mp, int mol_type, int bc_ix, par_ix_t *par_ix,
        int t_ix, f_t *prob) {

    int s_ix = par_ix->t_ix[t_ix];
    uint16_t M = mp->M;
    uint8_t c_ix = s_ix == M ? 0 : 1; // 0 is ambient, 1 is cell
    f_t al;
    if (mol_type == RNA_IX) {
        al = mp->alpha_rna1[CMI(bc_ix, par_ix->hs_ix, mp->D)];
    } else if (mol_type == ATAC_IX) {
        al = mp->alpha_atac1[CMI(bc_ix, par_ix->hs_ix, mp->D)];
    } else {
        return err_msg(-1, 0, "pr_tdm: invalid `mol_type` index '%i'", mol_type);
    }
    if (num_invalid(al))
        return err_msg(-1.0, 0, "pr_tdm: alpha=%e is invalid", al);
    f_t pr = 0;
    switch (par_ix->hd) {
        case 0:
            if (c_ix == 0) pr = 1.0;
            else
                return err_msg(-1, 0, "pr_tdm: hd = 0 but s_ix != M", par_ix->hd);
            break;
        case 1:
            if (c_ix == 0) pr = al;
            else pr = 1.0 - al;
            break;
        case 2:
            if (c_ix == 0) pr = al;
            else pr = (1.0-al)/2.0;
            break;
        default:
            return err_msg(-1, 0, "pr_tdm: hd=%i is invalid, there is a bug", par_ix->hd);
    }
    *prob = pr;
    return(0);
}

f_t pr_gamma_var(f_t *gamma, uint32_t v_ix, uint8_t allele, 
        uint32_t s_ix, uint32_t gamma_nrow, f_t tau){
    if (allele > 1) return(1.0); // return 1 if allele is missing
    uint32_t eix = CMI(s_ix, v_ix, gamma_nrow);
    f_t ap = gamma[eix];
    if (ap < 0) return 1.0; // if allele is missing
    if (allele == 0) ap = 1.0 - ap;
    f_t p_be0 = (1.0 - tau) * ap;
    f_t p_be1 = tau * 0.5;
    f_t p_bdm = p_be0 + p_be1;
    return(p_bdm);
}

f_t mol_info(mdl_pars_t *mp, mdl_mlcl_t *mlcl, par_ix_t *par_ix) {
    // if empty, nothing to check
    if (par_ix->hd == 0)
        return 0;

    uint16_t M = mp->M;
    uint32_t gamma_nrow = M + 1;

    f_t info = 0;

    size_t i, n_vars = mv_size(&mlcl->var_ixs);

    for (i = 0; i < n_vars; ++i){
        // get base call for variant allele
        uint32_t pack = mv_i(&mlcl->var_ixs, i);
        uint32_t v_ix;
        uint8_t allele;
        mdl_mlcl_unpack_var(pack, &v_ix, &allele);
        if (allele > 1) continue; // if missing

        int t_im;
        f_t t_info = 0;
        for (t_im = 1; t_im < par_ix->t_n; ++t_im){
            int s_ix = par_ix->t_ix[t_im];
            uint32_t g1ix = CMI(s_ix, v_ix, gamma_nrow);
            uint32_t g2ix = CMI(M, v_ix, gamma_nrow);
            f_t ap1 = mp->gamma[g1ix], ap2 = mp->gamma[g2ix];
            if (ap1 < 0 || ap2 < 0)
                continue;
            t_info += fabs(ap1 - ap2) / (par_ix->t_n - 1);
        }
        // t_info = t_info * t_info;
        info += t_info;
    }
    if (n_vars)
        info /= (f_t)n_vars;

    return info;
}

int bc_info(mdl_pars_t *mp, kbtree_t(kb_mdl_mlcl) *mols, par_ix_t *par_ix,
            f_t *info) {
    if (par_ix->hd == 0)
        return 0;

    *info = 0;

    kbitr_t itr;
    kb_itr_first(kb_mdl_mlcl, mols, &itr); 
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols, &itr)){
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        uint32_t n_vars = mv_size(&mlcl->var_ixs);
        if (n_vars < 1)
            continue;

        *info += mol_info(mp, mlcl, par_ix);
    }

    return 0;
}

int p_var(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int s_ix, f_t *prob){
    uint16_t M = mp->M;
    uint32_t gamma_nrow = M + 1;

    // probability terms
    f_t pbm = 1;

    // Pr(B_dm | T_dm, \rho)
    size_t i, n_vars = mv_size(&mlcl->var_ixs);
    for (i = 0; i < n_vars; ++i){
        // get base call for variant allele
        uint32_t pack = mv_i(&mlcl->var_ixs, i);
        uint32_t v_ix;
        uint8_t allele;
        mdl_mlcl_unpack_var(pack, &v_ix, &allele);
        if (allele > 1) continue; // if missing

        // get base quality score
        // set to default tau if missing
        uint8_t bqual = mv_i(&mlcl->bquals, i);
        f_t perr = bqual != 0xff ? phred_to_perr(bqual) : mp->tau;

        f_t pp = pr_gamma_var(mp->gamma, v_ix, allele, s_ix, gamma_nrow, perr);
        if (prob_invalid(pp))
            return err_msg(-1, 0, "p_var: invalid base probability value '%e'", pp);
        pbm *= pp;
    }
    if (prob_invalid(pbm))
        return err_msg(-1, 0, "p_var: invalid variant probability value '%e'", pbm);

    *prob = pbm;

    return 0;
}

int p_bd(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int mol_type, int bc_ix, par_ix_t *par_ix,
        f_t *p_b, f_t *psum) {
    *psum = 0;

    int t_im, ret;
    for (t_im = 0; t_im < par_ix->t_n; ++t_im){
        p_b[t_im] = 1;

        f_t p_t = -1;
        ret = pr_tdm(mp, mol_type, bc_ix, par_ix, t_im, &p_t);
        if (ret < 0)
            return -1;

        f_t p_v = -1;
        ret = p_var(mp, mlcl, par_ix->t_ix[t_im], &p_v);
        if (ret < 0)
            return -1;

        p_b[t_im] = p_t * p_v;
        //if (prob_invalid(p_b[t_im]))
        //    return err_msg(-1, 0, "p_bd: prob=%f is invalid", p_b[t_im]);
        *psum += p_b[t_im];
    }

    return 0;
}

int d_p_bd_d_alpha(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int mol_type, int bc_ix,
                   par_ix_t *par_ix, f_t *d1, f_t *d2, f_t *p_ba, f_t *p_bs) {

    f_t x = 0.0, psum = 0.0;
    *p_ba = 0.0;
    *p_bs = 0.0;

    int t_im, ret;
    for (t_im = 0; t_im < par_ix->t_n; ++t_im){
        int s_ix = par_ix->t_ix[t_im];
        if (s_ix < 0)
            return err_msg(-1, 0, "d_p_bd_d_alpha: invalid s_ix=%i", s_ix);

        f_t p_t = -1;
        ret = pr_tdm(mp, mol_type, bc_ix, par_ix, t_im, &p_t);
        if (ret < 0)
            return -1;

        f_t p_v = -1;
        ret = p_var(mp, mlcl, s_ix, &p_v);
        if (ret < 0)
            return -1;

        psum += p_t * p_v;
        if (t_im == 0) {
            x += p_v;
            *p_ba += p_v;
       } else if (t_im > 0 && par_ix->hd == 1) {
            x -= p_v;
            *p_bs += p_v;
        } else if (t_im > 0 && par_ix->hd == 2) {
            x -= 0.5 * p_v;
            *p_bs += 0.5 * p_v;
        } else {
            return err_msg(-1, 0, "d_p_bd_d_alpha: invalid internal parameters");
        }
    }
    if (iszero(x) && iszero(psum)) {
        *d1 = 0;
        *d2 = 0;
        return 0;
    }
    *d1 = x / psum;
    *d2 = ((-1.0) * (x * x)) / (psum * psum);

    if (*d2 > 0)
        return err_msg(-1, 0, "d_p_bd_d_alpha: result d2/da2=%e is invalid", *d2);

    return 0;
}

int d_lpd_dlogita(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int mol_type, int bc_ix,
                  par_ix_t *par_ix, f_t *d1, f_t *d2) {

    f_t pb = 0.0, p1sum = 0.0, p2sum = 0.0;

    // get alpha parameter
    f_t al;
    if (mol_type == RNA_IX) {
        al = mp->alpha_rna1[CMI(bc_ix, par_ix->hs_ix, mp->D)];
    } else if (mol_type == ATAC_IX) {
        al = mp->alpha_atac1[CMI(bc_ix, par_ix->hs_ix, mp->D)];
    } else {
        return err_msg(-1, 0, "d_lpd_dlogita: invalid `mol_type` index '%i'", mol_type);
    }
    f_t pbs_norm = par_ix->hd == 2 ? 0.5 : 1.0;

    // derivatives of logistic function w.r.t. alpha
    f_t fd1 = al * (1 - al);
    f_t fd2 = al * (1 - al) * (1 - (2*al));

    int t_im, ret;
    for (t_im = 0; t_im < par_ix->t_n; ++t_im){
        int s_ix = par_ix->t_ix[t_im];
        if (s_ix < 0)
            return err_msg(-1, 0, "d_lpd_dlogita: invalid s_ix=%i", s_ix);

        f_t p_t = -1;
        ret = pr_tdm(mp, mol_type, bc_ix, par_ix, t_im, &p_t);
        if (ret < 0)
            return -1;

        f_t p_v = -1;
        ret = p_var(mp, mlcl, s_ix, &p_v);
        if (ret < 0)
            return -1;

        pb += p_t * p_v;

        if (t_im == 0) {
            p1sum += fd1 * p_v;
            p2sum += fd2 * p_v;
        } else {
            p1sum -= pbs_norm * fd1 * p_v;
            p2sum -= pbs_norm * fd2 * p_v;
        }
    }
    if (iszero(p1sum) && iszero(pb)) {
        err_msg(0, 1, "d_lpd_dlogita: p1sum and pb are 0");
        *d1 = 0;
        return 0;
    }
    if (iszero(p2sum) && iszero(pb)) {
        err_msg(0, 1, "d_lpd_dlogita: p2sum and pb are 0");
        *d2 = 0;
        return 0;
    }
    *d1 = p1sum / pb;
    *d2 = (p2sum / pb) - (p1sum * p1sum) / (pb * pb);

    //if (*d2 > 0)
    //    return err_msg(-1, 0, "d_lpd_dlogita: result d2/da2=%f is invalid", *d2);

    return 0;
}

/*******************************************************************************
 * Expectation
 ******************************************************************************/

int mdl_init_alpha(mdl_t *mdl, f_t val) {
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_init_alpha: argument is null");

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;

    uint8_t n_mod = 2;
    f_t *alpha1[2] = {NULL, NULL};
    alpha1[ATAC_IX] = mdl->mp->alpha_atac1;
    alpha1[RNA_IX] = mdl->mp->alpha_rna1;
    int has_mod[2] = {0, 0};
    has_mod[ATAC_IX] = mdl->has_atac;
    has_mod[RNA_IX] = mdl->has_rna;

    if (val <= 0 || val >= 1)
        return err_msg(-1, 0, "mdl_init_alpha: val='%e' must be within (0, 1)",
                       val);

    uint8_t mod_ix;
    for (mod_ix = 0; mod_ix < n_mod; ++mod_ix) {
        if (!has_mod[mod_ix])
            continue;

        uint32_t bc_ix;
        for (bc_ix = 0; bc_ix < mdl->mp->D; ++bc_ix){
            // if no barcode
            int efl = bflg_get(&bd->absent_bc, bc_ix);
            if (efl == 1)
                continue;

            int hs_ix;
            for (hs_ix = 0; hs_ix < mdl->n_hs; ++hs_ix) {
                uint32_t aix = CMI(bc_ix, hs_ix, mdl->mp->D);
                // if fixed ambient, skip
                int afl = bflg_get(&bd->amb_flag, bc_ix);
                if (hs_ix == 0) {
                    alpha1[mod_ix][aix] = 1;
                } else if (afl == 1) {
                    alpha1[mod_ix][aix] = 1;
                } else {
                    alpha1[mod_ix][aix] = val;
                }
            }
        }
    }
    return 0;
}

int mdl_get_mmle_alpha(mdl_t *mdl, int *ixs, uint32_t ix_len) {
    if (mdl == NULL || ixs == NULL)
        return err_msg(-1, 0, "mdl_get_mmle_alpha: arguments are null");

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "mdl_get_mmle_alpha: model parameters are missing");
    if (mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "mdl_get_mmle_alpha: barcode count data is missing");
    if (mdl->mdl_bc_dat->all_bcs->n < 1)
        return err_msg(-1, 0, "mdl_get_mmle_alpha: no barcodes found");

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    uint32_t hs_ix;
    par_ix_t par_ix;
    par_ix_init(&par_ix);


    // store data
    uint8_t n_mod = 2;
    f_t *alpha1[2] = {NULL, NULL};
    alpha1[ATAC_IX] = mdl->mp->alpha_atac1;
    alpha1[RNA_IX] = mdl->mp->alpha_rna1;
    f_t *alpha_se[2] = {NULL, NULL};
    alpha_se[ATAC_IX] = mdl->mp->alpha_atac_se;
    alpha_se[RNA_IX] = mdl->mp->alpha_rna_se;
    f_t *alpha_info[2] = {NULL, NULL};
    alpha_info[ATAC_IX] = mdl->mp->alpha_atac_info;
    alpha_info[RNA_IX] = mdl->mp->alpha_rna_info;
    int has_mod[2] = {0, 0};
    has_mod[ATAC_IX] = mdl->has_atac;
    has_mod[RNA_IX] = mdl->has_rna;
    char mod_str[2][20];
    strcpy(mod_str[ATAC_IX], "ATAC");
    strcpy(mod_str[RNA_IX], "RNA");
    f_t alpha_prior[2];
    alpha_prior[ATAC_IX] = mdl->mp->alpha_prior_atac;
    alpha_prior[RNA_IX] = mdl->mp->alpha_prior_rna;

    // prior for alpha
    f_t al_prior_w = mdl->mp->alpha_prior_w;

    // priors needed
    if (alpha_prior[0] <= 0 || alpha_prior[0] >= 1)
        return err_msg(-1, 0, "mdl_get_mmle_alpha: prior='%e' needs to be within (0,1)",
                       alpha_prior[0]);
    if (alpha_prior[1] <= 0 || alpha_prior[1] >= 1)
        return err_msg(-1, 0, "mdl_get_mmle_alpha: prior='%e' needs to be within (0,1)",
                       alpha_prior[1]);
    if (al_prior_w <= 0)
        return err_msg(-1, 0, "mwl_get_mmle_alpha: prior weight='%e' needs to be greater than 0",
                       al_prior_w);

    // upper and lower bounds for initial alpha
    f_t lbnd_eps = 1e-12, ubnd_eps = 1.0 - 1e-10;

    uint32_t i;
    for (i = 0; i < ix_len; ++i){
        int bc_ix = ixs[i];

        // if no barcode
        int efl = bflg_get(&bd->absent_bc, bc_ix);
        if (efl == 1)
            continue;

        // if fixed ambient, skip
        int afl = bflg_get(&bd->amb_flag, bc_ix);
        if (afl == 1)
            continue;

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);
        kbtree_t(kb_mdl_mlcl) *mols[2] = {NULL, NULL};
        mols[RNA_IX] = bc_dat.rna;
        mols[ATAC_IX] = bc_dat.atac;

        // loop over RNA and ATAC
        uint8_t mod_ix;
        for (mod_ix = 0; mod_ix < n_mod; ++mod_ix) {
            if (!has_mod[mod_ix])
                continue;

            // loop over singlet possibilities
            for (hs_ix = 0; hs_ix < n_hs; ++hs_ix) {
                // get (h,s,t) from index, store in par_ix
                if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                    return -1;

                uint32_t aix = CMI(bc_ix, hs_ix, mdl->mp->D);

                // skip if empty
                if (par_ix.hd == 0) {
                    alpha1[mod_ix][aix] = 1.0;
                    alpha_se[mod_ix][aix] = NAN;
                    alpha_info[mod_ix][aix] = NAN;
                    continue;
                }

                f_t delta = 0.0, delta_thresh = mdl->eps;
                int alpha_iter, iter_max = 500;

                //if (kb_size(mols_rna) > 100)
                //    printf("BC hs=%i\n", hs_ix);

                f_t best_alpha = alpha1[mod_ix][aix];
                for (alpha_iter = 0; alpha_iter < iter_max; ++alpha_iter) {
                    // initialize with smoothing
                    f_t d1_sum = 0, d2_sum = 0;
                    d1_sum = 0.0;
                    d2_sum = 0.0;
                    d1_sum += (alpha_prior[mod_ix] * al_prior_w) * (1 / best_alpha);
                    d1_sum += ((1 - alpha_prior[mod_ix]) * al_prior_w) * (-1 / (1 - best_alpha));
                    d2_sum -= (alpha_prior[mod_ix] * al_prior_w) * (1 / (best_alpha*best_alpha));
                    d2_sum -= ((1 - alpha_prior[mod_ix]) * al_prior_w) * (1 / ((1 - best_alpha)*(1 - best_alpha)));

                    if (kb_size(mols[mod_ix]) == 0) {
                        alpha1[mod_ix][aix] = alpha_prior[mod_ix];
                        alpha_se[mod_ix][aix] = NAN;
                        alpha_info[mod_ix][aix] = NAN;
                        break;
                    }

                    kbitr_t itr;
                    kb_itr_first(kb_mdl_mlcl, mols[mod_ix], &itr); 
                    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols[mod_ix], &itr)){
                        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
                        uint32_t n_vars = mv_size(&mlcl->var_ixs);
                        if (n_vars < 1)
                            continue;

                        // calculate derivative at current alpha
                        f_t d1, d2, p_ba, p_bs;
                        int pret = d_p_bd_d_alpha(mdl->mp, mlcl, mod_ix, bc_ix,
                                                  &par_ix, &d1, &d2, &p_ba, &p_bs);
                        if (pret < 0)
                            return -1;
                        if (isinf(d2))
                            return err_msg(-1, 0, "mdl_get_mmle_alpha: d2/da2='%e', likelihood collapse. "
                                           "There may be a bug.", d2);
                        if (isnan(d2))
                            return err_msg(-1, 0, "mdl_get_mmle_alpha: d2/da2='nan'. There may be a bug.");

                        d1_sum += mlcl->counts * d1;
                        d2_sum += mlcl->counts * d2;
                    }

                    // sanity check
                    if (isnan(d1_sum))
                        return err_msg(-1, 0, "mdl_get_mmle_alpha: barcode 1st derivative is nan (a=%e)", best_alpha);
                    if (isnan(d2_sum))
                        return err_msg(-1, 0, "mdl_get_mmle_alpha: barcode 2nd derivative is nan (a=%e)", best_alpha);
                    if (isinf(d1_sum))
                        return err_msg(-1, 0, "mdl_get_mmle_alpha: barcode 1st derivative is inf (a=%e)", best_alpha);
                    if (isinf(d2_sum))
                        return err_msg(-1, 0, "mdl_get_mmle_alpha: barcode 2nd derivative is inf (a=%e)", best_alpha);

                    // update alpha, ensure update is valid
                    f_t alpha_t1;
                    if (iszero(d2_sum) && iszero(d1_sum)) {
                        alpha_t1 = best_alpha;
                    } else {
                        alpha_t1 = best_alpha - (d1_sum / d2_sum);
                    }

                    // bound alpha if update is outside domain
                    if (alpha_t1 < lbnd_eps)
                        alpha_t1 = lbnd_eps;
                    if (alpha_t1 > ubnd_eps)
                        alpha_t1 = ubnd_eps;

                    // get delta and store new alpha
                    delta = fabs(alpha_t1 - best_alpha) / fabs(best_alpha);
                    best_alpha = alpha_t1;
                    alpha1[mod_ix][aix] = best_alpha;

                    // update std. error estimate
                    if (iszero(d2_sum))
                        alpha_se[mod_ix][aix] = NAN;
                    else
                        alpha_se[mod_ix][aix] = 1.0 / sqrt(-d2_sum);
                    if (delta < delta_thresh)
                        break;
                    //if (isnan(delta) || isinf(delta) || fabs(delta) < delta_thresh)
                    //    break;
                }
                //if (par_ix.hd == 1 && kb_size(mols[mod_ix]) > 100) {
                //    printf("%s ML [iters=%i] [BC %i] [HS %i] alpha=%e, delta=%e\n",
                //           mod_str[mod_ix], alpha_iter+1, bc_ix, hs_ix, best_alpha, delta);
                //}
                // if (alpha_iter == iter_max)
                //     err_msg(0, 1,
                //             "mdl_get_mmle_alpha: failed to converge for RNA alpha of BC %i "
                //             "[%i reads] [alpha=%e] (delta=%e) alpha_range=[%e,%e]",
                //             bc_ix, kb_size(mols[mod_ix]), best_alpha, delta, alpha_min, alpha_max);

                // effective number of reads
                if (bc_info(mdl->mp, mols[mod_ix], &par_ix, alpha_info[mod_ix] + aix) < 0)
                    return -1;
            }
        }
    }

    return 0;
}

int mdl_em_alpha(mdl_t *mdl, int *ixs, uint32_t ix_len) {
    if (mdl == NULL || ixs == NULL)
        return err_msg(-1, 0, "mdl_get_mmle_alpha: arguments are null");

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "mdl_get_mmle_alpha: model parameters are missing");
    if (mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "mdl_get_mmle_alpha: barcode count data is missing");
    if (mdl->mdl_bc_dat->all_bcs->n < 1)
        return err_msg(-1, 0, "mdl_get_mmle_alpha: no barcodes found");

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    uint32_t hs_ix;

    // store data
    uint8_t n_mod = 2;
    f_t *alpha1[2] = {NULL, NULL};
    alpha1[ATAC_IX] = mdl->mp->alpha_atac1;
    alpha1[RNA_IX] = mdl->mp->alpha_rna1;
    f_t *alpha_info[2] = {NULL, NULL};
    alpha_info[ATAC_IX] = mdl->mp->alpha_atac_info;
    alpha_info[RNA_IX] = mdl->mp->alpha_rna_info;
    char mod_str[2][20];
    strcpy(mod_str[ATAC_IX], "ATAC");
    strcpy(mod_str[RNA_IX], "RNA");

    // prior for alpha
    f_t al_prior_w = mdl->mp->alpha_prior_w;
    f_t alpha_prior[2];
    alpha_prior[ATAC_IX] = mdl->mp->alpha_prior_atac;
    alpha_prior[RNA_IX] = mdl->mp->alpha_prior_rna;

    // priors needed
    if (alpha_prior[0] <= 0 || alpha_prior[0] >= 1)
        return err_msg(-1, 0, "mdl_em_alpha: prior=%e needs to be within (0,1)",
                       alpha_prior[0]);
    if (alpha_prior[1] <= 0 || alpha_prior[1] >= 1)
        return err_msg(-1, 0, "mdl_em_alpha: prior=%e needs to be within (0,1)",
                       alpha_prior[1]);
    if (al_prior_w <= 0)
        return err_msg(-1, 0, "mwl_em_alpha: prior weight=%e needs to be greater than 0",
                       al_prior_w);

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    uint32_t i;
    for (i = 0; i < ix_len; ++i){
        int bc_ix = ixs[i];

        // if no barcode
        int efl = bflg_get(&bd->absent_bc, bc_ix);
        if (efl == 1)
            continue;

        // if fixed ambient, skip
        int afl = bflg_get(&bd->amb_flag, bc_ix);
        if (afl == 1)
            continue;

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);
        kbtree_t(kb_mdl_mlcl) *mols[2] = {NULL, NULL};
        mols[RNA_IX] = bc_dat.rna;
        mols[ATAC_IX] = bc_dat.atac;

        // loop over singlet possibilities
        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix) {
            // get (h,s,t) from index, store in par_ix
            if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                return -1;

            // skip if empty
            if (par_ix.hd == 0)
                continue;

            uint32_t aix = CMI(bc_ix, hs_ix, mdl->mp->D);

            // loop over RNA and ATAC
            uint8_t mod_ix;
            for (mod_ix = 0; mod_ix < n_mod; ++mod_ix) {

                f_t delta = 0.0, delta_thresh = 1e-15;
                int alpha_iter, iter_max = 500;

                //if (kb_size(mols_rna) > 100)
                //    printf("BC hs=%i\n", hs_ix);
                f_t lp_prev = -INFINITY;
                f_t best_alpha = alpha1[mod_ix][aix], init_alpha = best_alpha;
                // alpha1[mod_ix][aix] = best_alpha;
                for (alpha_iter = 0; alpha_iter < iter_max; ++alpha_iter) {
                    f_t asum[2] = {alpha_prior[mod_ix] * al_prior_w, (1 - alpha_prior[mod_ix]) * al_prior_w}; // 1st index is ambient
                    f_t lp = 0.0;
                    lp += (alpha_prior[mod_ix] * al_prior_w) * log(best_alpha);
                    lp += ((1 - alpha_prior[mod_ix]) * al_prior_w) * log(1 - best_alpha);

                    kbitr_t itr;
                    kb_itr_first(kb_mdl_mlcl, mols[mod_ix], &itr); 
                    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols[mod_ix], &itr)){
                        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
                        uint32_t n_vars = mv_size(&mlcl->var_ixs);
                        if (n_vars < 1)
                            continue;

                        f_t psum = 0, p_tgpb[3] = {1,1,1};
                        int pret = p_bd(mdl->mp, mlcl, mod_ix, bc_ix, &par_ix, p_tgpb, &psum);
                        if (pret < 0)
                            return -1;
                        if (psum < 0 || num_invalid(psum))
                            return err_msg(-1, 1, "mdl_em_alpha: invalid read likelihood '%e'", psum);
                        lp += mlcl->counts * log(psum);

                        int t_im;
                        for (t_im = 0; t_im < par_ix.t_n; ++t_im){
                            int s_ix = par_ix.t_ix[t_im];
                            uint8_t c_ix = s_ix == mdl->mp->M ? 0 : 1; // 0 is ambient, 1 is cell
                            f_t pt = p_tgpb[t_im] / psum;

                            // alpha
                            asum[c_ix] += mlcl->counts * pt;
                        }

                    }
                    f_t a_tot = asum[0] + asum[1];
                    f_t alpha_t1 = asum[0] / a_tot;
                    delta = alpha_t1 - alpha1[mod_ix][aix];
                    alpha1[mod_ix][aix] = alpha_t1;
                    best_alpha = alpha_t1;
                    // check for llk decrease
                    if (lp < lp_prev)
                        err_msg(-1, 1, "mdl_em_alpha: [iter %i] llk='%e' decreased [%s] [alpha=%e] [delta=%e]",
                                alpha_iter, lp, mod_str[mod_ix], alpha_t1, delta);
                    lp_prev = lp;
                    // if (par_ix.hd == 1 && kb_size(mols[mod_ix]) > 100) {
                    //     printf("%s EM [iter %i] [BC %i] [HS %i] alpha=%e, 1-alpha=%e, delta=%e\n",
                    //            mod_str[mod_ix], alpha_iter+1, bc_ix, hs_ix, best_alpha, 1-best_alpha, delta);
                    // }
                    // f_t f_i = (ar[0] / (best_rna_alpha * best_rna_alpha)) + (ar[1] / ((1 - best_rna_alpha)*(1 - best_rna_alpha)));
                    // mdl->mp->alpha_rna_se[aix] = 1.0 / sqrt(f_i);
                    if (fabs(delta) < delta_thresh)
                        break;
                }
                if (alpha_iter > 0) {
                    printf("%s EM [iters=%i] [BC %i] [HS %i] alpha=%e, init=%e, delta=%e\n",
                           mod_str[mod_ix], alpha_iter+1, bc_ix, hs_ix, best_alpha, init_alpha, delta);
                }
                // if (par_ix.hd == 1 && kb_size(mols[mod_ix]) > 100) {
                //     printf("%s EM [iters=%i] [BC %i] [HS %i] alpha=%e, init=%e, delta=%e\n",
                //            mod_str[mod_ix], alpha_iter+1, bc_ix, hs_ix, best_alpha, init_alpha, delta);
                // }
                // effective number of reads
                if (bc_info(mdl->mp, mols[mod_ix], &par_ix, alpha_info[mod_ix] + aix) < 0)
                    return -1;
            }
        }
    }

    return 0;
}

int mdl_set_lp(mdl_t *mdl, f_t *lp_htd, uint32_t bc_ix) {
    if (mdl == NULL || lp_htd == NULL)
        return err_msg(-1, 0, "mdl_set_lp: argument is null");

    uint32_t n_hs = mdl->n_hs;
    par_ix_t par_ix;
    par_ix_init(&par_ix);

    // set log Pr(H_d, S_d, X_d) to mdl
    memcpy(&mdl->sub_lp_hs[CMI(0, bc_ix, n_hs)], lp_htd, n_hs * sizeof(f_t));

    // set log Pr(H_d, X_d) to mdl
    mdl->sub_lp_h[CMI(0, bc_ix, 3)] = -INFINITY;
    mdl->sub_lp_h[CMI(1, bc_ix, 3)] = -INFINITY;
    mdl->sub_lp_h[CMI(2, bc_ix, 3)] = -INFINITY;
    
    uint32_t hs_ix;
    for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
        if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
            return -1;
        f_t prev_lp_h = mdl->sub_lp_h[CMI(par_ix.hd, bc_ix, 3)];
        mdl->sub_lp_h[CMI(par_ix.hd, bc_ix, 3)] = logsum2expd(prev_lp_h, lp_htd[hs_ix]);
    }

    int lsret = 0;
    // Pr(X_d | \Theta)
    mdl->sub_lp_x[bc_ix] = logsumexpd(mdl->sub_lp_h + (3 * bc_ix), 3, &lsret);
    if (lsret < 0)
        return err_msg(-1, 0, "mdl_set_lp: could not logsumexp");
    if (num_invalid(mdl->sub_lp_x[bc_ix]))
        return err_msg(-1, 0, "mdl_set_lp: sub_lp_x[%i]=%.6e is invalid\n",
                       bc_ix, mdl->sub_lp_x[bc_ix]);

    return 0;
}

// sub model
int mdl_sub_e(mdl_t *mdl, int *ixs, uint32_t ix_len) {
    if (mdl == NULL || ixs == NULL)
        return err_msg(-1, 0, "mdl_sub_e: arguments are null");

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "mdl_sub_e: model parameters are missing");
    if (mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "mdl_sub_e: barcode count data is missing");
    if (mdl->mdl_bc_dat->all_bcs->n < 1)
        return err_msg(-1, 0, "mdl_sub_e: no barcodes found");

    uint32_t M = mdl->mp->M;
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    uint32_t hs_ix;

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    // modality and alpha data
    uint8_t n_mod = 2;
    int has_mod[2] = {0, 0};
    has_mod[ATAC_IX] = mdl->has_atac;
    has_mod[RNA_IX] = mdl->has_rna;
    f_t *alpha1[2] = {NULL, NULL};
    alpha1[ATAC_IX] = mdl->mp->alpha_atac1;
    alpha1[RNA_IX] = mdl->mp->alpha_rna1;
    f_t alpha_prior[2];
    alpha_prior[ATAC_IX] = mdl->mp->alpha_prior_atac;
    alpha_prior[RNA_IX] = mdl->mp->alpha_prior_rna;

    // prior for alpha
    f_t al_prior_w = mdl->mp->alpha_prior_w;


    // pre-calculate P(H_d, S_d). Same for all droplets
    f_t *lp_hs_v = calloc(n_hs, sizeof(f_t));
    for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
        if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
            return -1;
        f_t p_hd = -1, p_sd = -1;
        pr_hd(mdl->mp, &par_ix, &p_hd);
        pr_sd(mdl->mp, &par_ix, &p_sd);
        if (prob_invalid(p_hd) || prob_invalid(p_sd))
            return err_msg(-1, 0, "mdl_sub_e: failed to get pr_hd or pr_sd");
        lp_hs_v[hs_ix] = log(p_hd) + log(p_sd);
    }

    // store P(H_d, S_d, X_d) during loop
    f_t *lp_htd = calloc(n_hs, sizeof(f_t));
    if (lp_htd == NULL)
        return err_msg(-1, 0, "mdl_sub_e: %s", strerror(errno));

    // lambda counter
    f_t lambda_sums[3] = {0, 0, 0};
    f_t *pi_sums = calloc(M, sizeof(f_t));

    uint32_t i;
    for (i = 0; i < ix_len; ++i){
        // get barcode and bam data
        int bc_ix = ixs[i];

        // if no barcode
        int efl = bflg_get(&bd->absent_bc, bc_ix);
        if (efl == 1)
            continue;

        // if fixed, set Pr(empty) = 1
        int afl = bflg_get(&bd->amb_flag, bc_ix);
        uint32_t bc_n_ix = afl ? 1 : n_hs;

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);

        kbtree_t(kb_mdl_mlcl) *mols[2] = {NULL, NULL};
        mols[RNA_IX] = bc_dat.rna;
        mols[ATAC_IX] = bc_dat.atac;

        // initialize log Pr(H_d, S_d, X_d)
        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix) {
            if (hs_ix < bc_n_ix)
                lp_htd[hs_ix] = lp_hs_v[hs_ix];
            else
                lp_htd[hs_ix] = -INFINITY;
        }

        // Get P(H_d, S_d, X_d) by marginalizing T_dm
        for (hs_ix = 0; hs_ix < bc_n_ix; ++hs_ix){
            // get (h,s,t) from index, store in par_ix
            if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                return -1;

            uint32_t aix = CMI(bc_ix, hs_ix, mdl->mp->D);

            // loop over RNA and ATAC
            uint8_t mod_ix;
            for (mod_ix = 0; mod_ix < n_mod; ++mod_ix) {
                if (!has_mod[mod_ix])
                    continue;

                f_t best_alpha = alpha1[mod_ix][aix];
                if (par_ix.hd > 0) {
                    lp_htd[hs_ix] += (alpha_prior[mod_ix] * al_prior_w) * log(best_alpha);
                    lp_htd[hs_ix] += ((1 - alpha_prior[mod_ix]) * al_prior_w) * log(1 - best_alpha);
                }

                kbitr_t itr;
                kb_itr_first(kb_mdl_mlcl, mols[mod_ix], &itr); 
                for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols[mod_ix], &itr)){
                    mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
                    uint32_t n_vars = mv_size(&mlcl->var_ixs);
                    if (n_vars < 1)
                        continue;

                    f_t psum = 0, p_tgpb[3] = {1,1,1};
                    int pret = p_bd(mdl->mp, mlcl, mod_ix, bc_ix, &par_ix, p_tgpb, &psum);
                    if (pret < 0)
                        return -1;
                    lp_htd[hs_ix] += mlcl->counts * log(psum);
                }
            }
        }

        if (mdl_set_lp(mdl, lp_htd, bc_ix) < 0)
            return -1;

        // add lambda and pi parameters
        for (hs_ix = 0; hs_ix < bc_n_ix; ++hs_ix){
            // get (h,s,t) from index
            if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                return(-1);

            // get P(H,S|X)
            uint32_t mix = CMI(hs_ix, bc_ix, n_hs);
            f_t cp_hsx = exp(mdl->sub_lp_hs[mix] - mdl->sub_lp_x[bc_ix]);
            if (prob_invalid(cp_hsx))
                return err_msg(-1, 0, "mdl_sub_e: P(Z|X)=%.6e is an "
                               "invalid probability", cp_hsx);

            // lambda
            lambda_sums[par_ix.hd] += bc_dat.n_bc * cp_hsx;

            // pi
            if (par_ix.hd == 1)
                pi_sums[par_ix.s1] += bc_dat.n_bc * cp_hsx;
            if (par_ix.hd == 2) {
                pi_sums[par_ix.s1] += bc_dat.n_bc * cp_hsx;
                pi_sums[par_ix.s2] += bc_dat.n_bc * cp_hsx;
            }

        }
    }
    free(lp_htd);
    free(lp_hs_v);

    // add to tmp sum variables
    pthread_mutex_lock(&mdl->mp->sum_lock);

    // add lambda
    for (i = 0; i < 3; ++i){
        mdl->mp->_lambda_sum[i] += lambda_sums[i];
    }

    // add pi
    for (i = 0; i < M; ++i) {
        mdl->mp->_pi_sum[i] += pi_sums[i];
    }

    pthread_mutex_unlock(&mdl->mp->sum_lock);

    free(pi_sums);

    return 0;
}

/*******************************************************************************
 * Maximization
 ******************************************************************************/

int mdl_m_lambda(mdl_t *mdl) {
    f_t new_par = 0;
    uint32_t i;

    f_t lambda0_norm = 0, lambdad_norm = 0;
    // maximize lambda
    f_t lambda_tot = 0;
    for (i = 0; i < 3; ++i)
        lambda_tot += mdl->mp->_lambda_sum[i];
    for (i = 0; i < 3; ++i){
        lambda0_norm += mdl->mp->lambda[i] * mdl->mp->lambda[i];
        new_par = mdl->mp->_lambda_sum[i] / lambda_tot;
        f_t ldiff = mdl->mp->lambda[i] - new_par;
        lambdad_norm += ldiff * ldiff;
        // f_t delta = fabs(mdl->mp->lambda[i] - new_par) / fabs(mdl->mp->lambda[i]);
        // if (delta > mdl->mp->max_par_delta)
        //     mdl->mp->max_par_delta = delta;
        mdl->mp->lambda[i] = new_par;
    }
    lambda0_norm = sqrt(lambda0_norm);
    lambdad_norm = sqrt(lambdad_norm);
    f_t delta = lambdad_norm / lambda0_norm;
    if (delta > mdl->mp->max_par_delta)
        mdl->mp->max_par_delta = delta;

    return 0;
}

int mdl_m_pi(mdl_t *mdl) {
    f_t new_par = 0;
    uint16_t M = mdl->mp->M;
    f_t pi0_norm = 0, pid_norm = 0;
    f_t pi_tot = 0;
    unsigned int i, j;
    for (i = 0; i < M; ++i) {
        pi_tot += mdl->mp->_pi_sum[i];
    }
    for (i = 0; i < M; ++i) {
        pi0_norm = mdl->mp->pi[i] * mdl->mp->pi[i];
        new_par = mdl->mp->_pi_sum[i] / pi_tot;
        f_t ldiff = mdl->mp->pi[i] - new_par;
        pid_norm += ldiff * ldiff;
        // f_t delta = fabs(mdl->mp->pi[i] - new_par) / fabs(mdl->mp->pi[i]);
        // if (delta > mdl->mp->max_par_delta)
        //     mdl->mp->max_par_delta = delta;
        mdl->mp->pi[i] = new_par;
    }
    pi0_norm = sqrt(pi0_norm);
    pid_norm = sqrt(pid_norm);
    f_t delta = pid_norm / pi0_norm;
    if (delta > mdl->mp->max_par_delta)
        mdl->mp->max_par_delta = delta;

    f_t pi_sum = 0;
    for (i = 0; i < M; ++i){
        for (j = i+1; j < M; ++j){
            pi_sum += (mdl->mp->pi[i] * mdl->mp->pi[j]);
        }
    }
    mdl->mp->pi_d_sum = pi_sum;
    return 0;
}

/**
 * @brief Process molecules to update ambient sample fraction estimates.
 *
 * This function processes individual molecules (RNA or ATAC) to update the
 * estimates of ambient sample fractions. It handles missing values and ensures
 * that only valid probability values are used in calculations.
 *
 * @param mdl Pointer to the model structure.
 * @param mols Pointer to the kbtree containing molecule data.
 * @param M Number of samples.
 * @param gamma_nrow Number of rows in the gamma matrix.
 * @param d_pi_amb Array to store the updates for pi_amb.
 * @param n_tot Pointer to the total count of processed molecules.
 * @return 0 if successful, or a negative value if an error occurred.
 */
static int pi_amb_process_molecules(mdl_t *mdl,
                                    kbtree_t(kb_mdl_mlcl) *mols,
                                    uint16_t M,
                                    uint32_t gamma_nrow,
                                    f_t *d_pi_amb,
                                    f_t *n_tot) {
    kbitr_t itr;
    kb_itr_first(kb_mdl_mlcl, mols, &itr); 
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols, &itr)){
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        uint32_t n_vars = mv_size(&mlcl->var_ixs);

        size_t j;
        for (j = 0; j < n_vars; ++j){
            // get base call for variant allele
            uint32_t pack = mv_i(&mlcl->var_ixs, j);
            uint32_t v_ix;
            uint8_t allele;
            mdl_mlcl_unpack_var(pack, &v_ix, &allele);
            if (allele > 1) continue; // if missing

            // get base quality score
            // set to default tau if missing
            uint8_t bqual = mv_i(&mlcl->bquals, j);
            f_t perr = bqual != 0xff ? phred_to_perr(bqual) : mdl->mp->tau;

            f_t amb_alt_freq = mdl->mp->gamma[CMI(M, v_ix, gamma_nrow)];
            if (amb_alt_freq < 0)
                continue;

            if (prob_invalid(amb_alt_freq)) {
                char estr[] = "pi_amb_process_molecules: invalid ambient gamma value '%e'";
                return err_msg(-1, 0, estr, amb_alt_freq);
            }

            uint16_t s_ix;
            for (s_ix = 0; s_ix < M; ++s_ix) {
                f_t s_alt_freq = mdl->mp->gamma[CMI(s_ix, v_ix, gamma_nrow)];
                if (s_alt_freq < 0)
                    continue;

                if (prob_invalid(s_alt_freq)) {
                    char estr[] = "mdl_m_pi_amb: invalid sample gamma value '%e'";
                    return err_msg(-1, 0, estr, s_alt_freq);
                }

                f_t dp1 = allele * s_alt_freq / (amb_alt_freq);
                f_t dp2 = (1.0 - allele) * s_alt_freq / (1.0 - amb_alt_freq);
                f_t dp = dp1 - dp2;

                if (num_invalid(dp))
                    continue;

                d_pi_amb[s_ix] += (1 - perr) * mlcl->counts * dp;
            }
            *n_tot += (1 - perr) * mlcl->counts;
        }
    }
    return 0;
}

/**
 * @brief Process molecules using the likelihood-gradient update for pi_amb.
 *
 * Accumulate d / d pi_as log p(B | T=a) using model
 * A(x) = x^B (1-x)^(1-B)
 * p(B | x, \tau) = (1 − \tau) A(x) + \tau/2,
 * x = \gamma_av = \sum_s \pi_as \gamma_sv
 * dx/d\pi_as = \gamma_sv
 *
 * For each observed variant, adds
 * ((1 − \tau) A(x) / p) * ((B − x) / (x (1 − x))) * \gamma_sv
 * to d_pi_amb[s], weighted by mlcl->counts, skipping missing allele/genotype.
 *
 * @param mdl Pointer to the model structure.
 * @param mols Pointer to the kbtree containing molecule data.
 * @param M Number of samples.
 * @param gamma_nrow Number of rows in the gamma matrix.
 * @param d_pi_amb Array to store the updates for pi_amb.
 * @param n_tot Pointer to the total count of processed molecules.
 * @return 0 if successful, or a negative value if an error occurred.
 */
static int pi_amb_process_molecules_p(mdl_t *mdl,
                                      kbtree_t(kb_mdl_mlcl) *mols,
                                      uint16_t M,
                                      uint32_t gamma_nrow,
                                      f_t *d_pi_amb,
                                      f_t *n_tot) {
    const f_t eps = 1e-8;

    kbitr_t itr;
    kb_itr_first(kb_mdl_mlcl, mols, &itr);
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols, &itr)) {
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        uint32_t n_vars = mv_size(&mlcl->var_ixs);

        size_t j;
        for (j = 0; j < n_vars; ++j) {
            uint32_t pack = mv_i(&mlcl->var_ixs, j);
            uint32_t v_ix; uint8_t allele;
            mdl_mlcl_unpack_var(pack, &v_ix, &allele);
            if (allele > 1) continue; // skip missing

            // per-base error
            uint8_t bqual = mv_i(&mlcl->bquals, j);
            f_t tau = (bqual != 0xff) ? phred_to_perr(bqual) : mdl->mp->tau;

            // ambient alt allele frequency x
            f_t x = mdl->mp->gamma[CMI(M, v_ix, gamma_nrow)];
            if (x < 0) continue; // ambient missing for this variant
            if (prob_invalid(x)) {
                return err_msg(-1, 0, "pi_amb_process_molecules_p: invalid ambient gamma value '%e'", x);
            }
            // clip x
            if (x <= eps) x = eps;
            if (x >= 1.0 - eps) x = 1.0 - eps;

            // B in [0,1];  may replace w/ imputed
            f_t B = (f_t)allele;

            // A(x) and denominator p(x)
            f_t A = pow(x, B) * pow(1.0 - x, 1.0 - B);
            f_t denom = (1.0 - tau) * A + 0.5 * tau;
            if (!(denom > 0) || num_invalid(A) || num_invalid(denom))
                continue;

            // d log p / dx
            f_t dlogp_dx = ((1.0 - tau) * A / denom) * ((B - x) / (x * (1.0 - x)));
            if (num_invalid(dlogp_dx)) continue;

            // accumulate gradient over non-missing sample genotypes
            uint16_t s_ix;
            for (s_ix = 0; s_ix < M; ++s_ix) {
                f_t gsv = mdl->mp->gamma[CMI(s_ix, v_ix, gamma_nrow)];
                if (gsv < 0) continue;
                f_t contrib = mlcl->counts * gsv * dlogp_dx;
                if (num_invalid(contrib)) continue;
                d_pi_amb[s_ix] += contrib;
            }

            *n_tot += (1.0 - tau) * mlcl->counts;
        }
    }
    return 0;

}

/**
 * @brief Calculate the delta between old and new pi_amb values and update pi_amb.
 *
 * This function calculates the difference (delta) between the old and new pi_amb
 * values, updates the pi_amb values in the model, and returns the calculated delta.
 *
 * @param mdl Pointer to the model structure.
 * @param M Number of samples.
 * @param new_pi_ambp Array containing the new pi_amb values.
 * @return The calculated delta value, or a negative value if an error occurred.
 */
static f_t calculate_delta_and_update_pi_amb(mdl_t *mdl, uint16_t M, const f_t *new_pi_ambp) {
    f_t pi0_norm = 0.0, pid_norm = 0.0;
    uint16_t s_ix;
    for (s_ix = 0; s_ix < M; ++s_ix) {
        pi0_norm += mdl->mp->pi_amb[s_ix] * mdl->mp->pi_amb[s_ix];
        f_t ldiff = new_pi_ambp[s_ix] - mdl->mp->pi_amb[s_ix];
        pid_norm += ldiff * ldiff;
        if (prob_invalid(new_pi_ambp[s_ix])) {
            char estr[] = "calculate_delta_and_update_pi_amb: Invalid pi_amb value %e";
            return err_msg(-1, 0, estr, new_pi_ambp[s_ix]);
        }
        mdl->mp->pi_amb[s_ix] = new_pi_ambp[s_ix];
    }
    if (pi0_norm <= 1e-15)
        return err_msg(-1, 0, "calculate_delta_and_update_pi_amb: Invalid pi0_norm");
    
    return sqrt(pid_norm) / sqrt(pi0_norm);
}

int mdl_m_pi_amb(mdl_t *mdl) {
    if (mdl == NULL || mdl->mp == NULL || mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "mdl_m_pi_amb: invalid input: parameters are null");

    uint32_t i;
    uint16_t M = mdl->mp->M;
    uint32_t D = mdl->mp->D;
    uint32_t gamma_nrow = M + 1;
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;

    int has_mod[2] = {0, 0};
    has_mod[ATAC_IX] = mdl->has_atac;
    has_mod[RNA_IX] = mdl->has_rna;

    f_t *d_pi_amb = calloc(M, sizeof(f_t));
    f_t *new_pi_amb = calloc(M, sizeof(f_t));
    f_t *new_pi_ambp = calloc(M, sizeof(f_t));
    if (d_pi_amb == NULL || new_pi_amb == NULL || new_pi_ambp == NULL)
        return err_msg(-1, 0, "mdl_m_pi_amb: %s", strerror(errno));

    // get ambient barcodes
    mv_t(i32) amb_bcs;
    mv_init(&amb_bcs);
    for (i = 0; i < D; ++i) {
        int afl = bflg_get(&bd->amb_flag, i);
        if (afl)
            mv_push(i32, &amb_bcs, i);
    }
    uint32_t E = mv_size(&amb_bcs);
    if (E == 0) {
        mv_free(&amb_bcs);
        free(d_pi_amb);
        free(new_pi_amb);
        free(new_pi_ambp);
        return err_msg(0, 1, "mdl_m_pi_amb: no ambient barcodes found");
    }


    int max_iter = 10000, iter = 0;
    f_t delta_thresh = mdl->eps, delta = delta_thresh + 1;
    if (delta_thresh < 1e-15)
        delta_thresh = 1e-15;
    while (delta >= delta_thresh && iter < max_iter) {
        memset(d_pi_amb, 0, M * sizeof(f_t));
        memset(new_pi_amb, 0, M * sizeof(f_t));
        memset(new_pi_ambp, 0, M * sizeof(f_t));

        f_t n_tot = 0.0;
        for (i = 0; i < E; ++i){
            // get barcode and bam data
            uint32_t bc_ix = mv_i(&amb_bcs, i);

            // barcode count data
            mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);
            kbtree_t(kb_mdl_mlcl) *mols[2] = {NULL, NULL};
            mols[RNA_IX] = bc_dat.rna;
            mols[ATAC_IX] = bc_dat.atac;

            int ret_rna = pi_amb_process_molecules_p(
                mdl, bc_dat.rna, M, gamma_nrow, d_pi_amb, &n_tot
            );
            int ret_atac = pi_amb_process_molecules_p(
                mdl, bc_dat.atac, M, gamma_nrow, d_pi_amb, &n_tot
            );
            if (ret_rna < 0 || ret_atac < 0) {
                mv_free(&amb_bcs);
                free(d_pi_amb);
                free(new_pi_amb);
                free(new_pi_ambp);
                char estr[] = "mdl_m_pi_amb: failed to process molecules";
                return err_msg(-1, 0, estr);
            }
        }

        if (n_tot <= 1e-15) {
            err_msg(0, 1, "mdl_m_pi_amb: no reads found for ambient barcodes");
            break;
        }

        f_t vareps = 1.0 / n_tot;

        // Update pi_amb
        uint16_t s_ix;
        for (s_ix = 0; s_ix < M; ++s_ix) {
            new_pi_amb[s_ix] = mdl->mp->pi_amb[s_ix] + (vareps * d_pi_amb[s_ix]);
        }

        if (proj_splx(new_pi_amb, new_pi_ambp, M) < 0) {
            mv_free(&amb_bcs);
            free(d_pi_amb);
            free(new_pi_amb);
            free(new_pi_ambp);
            return err_msg(0, 1, "mdl_m_pi_amb: failed simplex projection");
        }

        // Calculate delta and update pi_amb
        delta = calculate_delta_and_update_pi_amb(mdl, M, new_pi_ambp);

        if (mdl_pars_set_gamma_amb(mdl->mp) < 0) {
            mv_free(&amb_bcs);
            free(d_pi_amb);
            free(new_pi_amb);
            free(new_pi_ambp);
            return -1;
        }

        ++iter;
    }

    if (iter == max_iter)
        err_msg(0, 1, "failed to converge calculating ambient sample fractions");

    mv_free(&amb_bcs);
    free(d_pi_amb);
    free(new_pi_amb);
    free(new_pi_ambp);

    return iter;
}

int mdl_sub_m(mdl_t *mdl) {
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_sub_m: mdl is null");

    if (mdl_m_lambda(mdl) < 0)
        return -1;

    if (mdl_m_pi(mdl) < 0)
        return -1;

    // update gamma
    if (mdl_pars_set_gamma_amb(mdl->mp) < 0) return(-1);

    return 0;
}

f_t mdl_get_avg_sng_alpha(mdl_t *mdl) {
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_get_avg_sng_alpha: argument is null");

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    uint8_t n_mod = 2;
    f_t *alpha1[2] = {NULL, NULL};
    alpha1[ATAC_IX] = mdl->mp->alpha_atac1;
    alpha1[RNA_IX] = mdl->mp->alpha_rna1;
    f_t *alpha_prior[2] = {NULL, NULL};
    alpha_prior[ATAC_IX] = &mdl->mp->alpha_prior_atac;
    alpha_prior[RNA_IX] = &mdl->mp->alpha_prior_rna;
    int has_mod[2] = {0, 0};
    has_mod[ATAC_IX] = mdl->has_atac;
    has_mod[RNA_IX] = mdl->has_rna;

    f_t max_delta = 0;

    int mod_ix;
    for (mod_ix = 0; mod_ix < n_mod; ++mod_ix) {
        if (!has_mod[mod_ix])
            continue;
        f_t alpha_avg = 0, cp_tot = 0;
        uint32_t bc_ix;
        for (bc_ix = 0; bc_ix < mdl->mp->D; ++bc_ix) {
            // if no barcode
            int efl = bflg_get(&mdl->mdl_bc_dat->absent_bc, bc_ix);
            if (efl == 1)
                continue;

            // if fixed
            int afl = bflg_get(&mdl->mdl_bc_dat->amb_flag, bc_ix);
            if (afl == 1)
                continue;

            int hs_ix;
            for (hs_ix = 0; hs_ix < mdl->n_hs; ++hs_ix) {
                if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                    return(-1);
                if (par_ix.hd != 1)
                    continue;

                uint32_t aix = CMI(hs_ix, bc_ix, mdl->n_hs);
                f_t cp_hs = exp(mdl->sub_lp_hs[aix] - mdl->sub_lp_x[bc_ix]);
                f_t alpha = alpha1[mod_ix][CMI(bc_ix, hs_ix, mdl->mp->D)];
                alpha_avg += cp_hs * alpha;
                cp_tot += cp_hs;
            }
        }
        if (iszero(cp_tot)) {
            err_msg(0, 1, "mdl_get_avg_sng_alpha: no singlets found, can't set alpha prior");
            break;
        }
        alpha_avg /= cp_tot;
        // bound alpha_avg to avoid collapsing likelihood to 0
        if (alpha_avg < 1e-4)
            alpha_avg = 1e-4;
        if (alpha_avg > 1 - 1e-4)
            alpha_avg = 1 - 1e-4;
        f_t t_delta = fabs(*alpha_prior[mod_ix] - alpha_avg) / *alpha_prior[mod_ix];
        if (t_delta > max_delta)
            max_delta = t_delta;
        *alpha_prior[mod_ix] = alpha_avg;
    }
    return max_delta;
}

int mdl_delta_q(f_t q1, f_t q2, f_t *q_delta) {
    *q_delta = (q2 - q1) / fabs(q1);
    return 0;
}

/*******************************************************************************
 * EM functions
 ******************************************************************************/

typedef struct thrd_args {
    mdl_t *mdl;
    int *ixs;
    uint32_t ix_len;
} thrd_args;

void *mdl_sub_thrd_fx(void *arg){
    thrd_args *t = (thrd_args *)arg;
    if (mdl_sub_e(t->mdl, t->ixs, t->ix_len) < 0) return((void *)-1);
    return((void *)0);
}

void *mdl_est_alpha_thrd_fx(void *arg) {
    thrd_args *t = (thrd_args *)arg;
    if (mdl_get_mmle_alpha(t->mdl, t->ixs, t->ix_len) < 0) return((void *)-1);
    // if (mdl_em_alpha(t->mdl, t->ixs, t->ix_len) < 0) return((void *)-1);
    return((void *)0);
}

int mdl_thrd_call(mdl_t *mdl, int *ixs, uint32_t ix_len){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_thrd_call: arguments are NULL");

    uint32_t i;
    if (mdl->threads < 1)
        return err_msg(-1, 0, "mdl_thrd_call: threads=%u is less than 1", mdl->threads);

    if (ix_len == 0)
        return err_msg(-1, 0, "mdl_thrd_call: invalid num. barcodes %u", ix_len);

    int step = (int)ix_len / (int)(mdl->threads);
    int rem = (int)ix_len % (int)(mdl->threads);
    int stepl = step + rem;
    /*
       div_t di = div((int)ix_len, (int)(mdl->threads));
       int step = di.quot;
       int stepl = di.rem == 0 ? step : di.rem;
       */

    uint32_t mt1 = mdl->threads - 1;
    thrd_args *targs = malloc(mdl->threads * sizeof(thrd_args));
    for (i = 0; i < mdl->threads; ++i){
        targs[i].mdl = mdl;
        targs[i].ixs = ixs + (step * i);
        targs[i].ix_len = i == (mt1) ? stepl : step;
    }

    pthread_t *ids = malloc(mdl->threads * sizeof(pthread_t));

    int err;
    for (i = 0; i < mdl->threads; ++i){
        err = pthread_create(ids + i, NULL, mdl_sub_thrd_fx, &targs[i]);
        if (err != 0)
            return err_msg(-1, 0, "mdl_thrd_call: could not create thread");
    }

    for (i = 0; i < mdl->threads; ++i){
        void *ret;
        err = pthread_join(ids[i], &ret);
        if (err != 0)
            return err_msg(-1, 0, "mdl_thrd_call: could not join thread");
        if ((long)ret < 0)
            return(-1);
    }

    free(targs);
    free(ids);

    return 0;
}

int mdl_thrd_est_alpha(mdl_t *mdl, int *ixs, uint32_t ix_len) {
    int step = (int)ix_len / (int)(mdl->threads);
    int rem = (int)ix_len % (int)(mdl->threads);
    int stepl = step + rem;

    uint32_t mt1 = mdl->threads - 1;
    thrd_args *targs = malloc(mdl->threads * sizeof(thrd_args));
    uint32_t i;
    for (i = 0; i < mdl->threads; ++i){
        targs[i].mdl = mdl;
        targs[i].ixs = ixs + (step * i);
        targs[i].ix_len = i == (mt1) ? stepl : step;
    }

    pthread_t *ids = malloc(mdl->threads * sizeof(pthread_t));

    int err;
    for (i = 0; i < mdl->threads; ++i){
        err = pthread_create(ids + i, NULL, mdl_est_alpha_thrd_fx, &targs[i]);
        if (err != 0)
            return err_msg(-1, 0, "mdl_thrd_est_alpha: could not create thread");
    }

    for (i = 0; i < mdl->threads; ++i){
        void *ret;
        err = pthread_join(ids[i], &ret);
        if (err != 0)
            return err_msg(-1, 0, "mdl_thrd_est_alpha: could not join thread");
        if ((long)ret < 0)
            return(-1);
    }

    free(targs);
    free(ids);

    return 0;
}

int mdl_em(mdl_t *mdl, obj_pars *objs){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_em_step: arguments are NULL");

    //f_t llk_s1, llk_s2;
    //f_t llk_delta;
    f_t prior = 1e-8;
    uint32_t i, j;
    uint32_t T = mdl->mp->T;
    uint32_t D = mdl->mp->D;
    int *ixs = malloc(D * sizeof(int));
    for (i = 0; i < D; ++i)
        ixs[i] = i;
    permute(ixs, D);

    int fret; // store function return

    fret = mdl_init_alpha(mdl, 0.5);
    if (fret < 0)
        return -1;

    if (objs->verbose)
        log_msg("fitting %u test droplets", T);

    for (i = 0; i < mdl->max_iter; ++i) {

        //if (i == 0) {
        //    llk_s1 = -INFINITY;
        //} else {
        //    fret = mdl_sub_llk(mdl, &llk_s1);
        //    if (fret < 0)
        //        return -1;
        //}

        fret = mdl_thrd_est_alpha(mdl, ixs, D);
        if (fret < 0)
            return -1;

        //fret = mdl_sub_llk(mdl, &llk_s2);
        //if (fret < 0)
        //    return -1;
        //if (mdl_delta_q(llk_s1, llk_s2, &llk_delta) < 0)
        //    return -1;
        //if (llk_delta < 0)
        //    err_msg(0, 1, "[alpha] llk decreased: [iter %u] llk1=%.6e llk2=%.6e", i, llk_s1, llk_s2);

        for (j = 0; j < mdl->max_iter; ++j){
            if (mdl_pars_check(mdl->mp) < 0)
                return -1;

            //if (i == 0 || j == 0) {
            //    llk_s1 = -INFINITY;
            //} else {
            //    fret = mdl_sub_llk(mdl, &llk_s1);
            //    if (fret < 0)
            //        return -1;
            //}

            // initialize sums to prior pseudocount
            fret = mdl_pars_reset_sums(mdl->mp, prior);
            if (fret < 0)
                return -1;

            fret = mdl_thrd_call(mdl, ixs, D);
            if (fret < 0)
                return -1;

            fret = mdl_sub_m(mdl);
            if (fret < 0)
                return -1;

            //fret = mdl_sub_llk(mdl, &llk_s2);
            //if (fret < 0)
            //    return -1;
            //if (mdl_delta_q(llk_s1, llk_s2, &llk_delta) < 0)
            //    return -1;
            //if (llk_delta < 0)
            //    err_msg(0, 1, "[H, S] llk decreased: [iter %u, %u] llk1=%.6e llk2=%.6e", i, j, llk_s1, llk_s2);

            if (mdl->mp->max_par_delta < mdl->eps)
                break;
        }

        f_t max_alpha_avg_delta = mdl_get_avg_sng_alpha(mdl);
        if (max_alpha_avg_delta < 0)
            return -1;

        if (mdl_pars_check(mdl->mp) < 0)
            return -1;

        if (objs->verbose)
            log_msg("iteration %u: delta=%.6e", i+1, max_alpha_avg_delta);

        if (max_alpha_avg_delta < mdl->eps)
            break;
    }

    if (objs->verbose) log_msg("finished fitting parameters");

    free(ixs);

    return(0);
}

int mdl_sub_llk(mdl_t *mdl, f_t *llk) {
    if (mdl == NULL || llk == NULL)
        return err_msg(-1, 0, "mdl_sub_llk: argument is null");

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    if (bd == NULL)
        return err_msg(-1, 0, "mdl_sub_llk: barcode data is not initialized");

    uint32_t i, n_bcs = (uint32_t)bd->all_bcs->n;
    f_t llkt = 0;
    for (i = 0; i < n_bcs; ++i){
        if (bflg_get(&bd->absent_bc, i))
            continue;
        llkt += mdl->sub_lp_x[i];
        if (num_invalid(llkt))
            return err_msg(-1, 0, "mdl_sub_llk: invalid log likelihood %e", llkt);
    }
    *llk = llkt;
    return(0);
}

int mdl_sub_best_hs(mdl_t *mdl) {
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_sub_best_hs: argument is null");

    if (mdl->sub_lp_hs == NULL || mdl->mp == NULL || mdl->mp->D == 0)
        return err_msg(-1, 0, "mdl_sub_best_hs: model must be initialized");

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    if (bd == NULL)
        return err_msg(-1, 0, "mdl_sub_best_hs: barcode data is not initialized");

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    uint32_t n_hs = mdl->n_hs;
    uint32_t D = mdl->mp->D;
    uint32_t M = mdl->mp->M;
    if (mdl->sub_best_hs == NULL)
        return err_msg(-1, 0, "mdl_sub_best_hs: model not initialized");

    uint32_t bc_i;
    for (bc_i = 0; bc_i < D; ++bc_i){
        if (bflg_get(&bd->absent_bc, bc_i))
            continue;
        // get best (H,S) index
        uint32_t b_hs_ix = 0;
        uint32_t t_hs_ix = 0;
        uint32_t b_h1_ix = 1;
        uint32_t b_h2_ix = M + 1;
        for (t_hs_ix = 1; t_hs_ix < n_hs; ++t_hs_ix) {
            if (hs_ix_get_pars(mdl->hs_ix, t_hs_ix, &par_ix) < 0)
                return(-1);

            uint32_t b_mix = CMI(b_hs_ix, bc_i, n_hs);
            uint32_t t_mix = CMI(t_hs_ix, bc_i, n_hs);
            uint32_t b1_mix = CMI(b_h1_ix, bc_i, n_hs);
            uint32_t b2_mix = CMI(b_h2_ix, bc_i, n_hs);
            if (mdl->sub_lp_hs[t_mix] > mdl->sub_lp_hs[b_mix])
                b_hs_ix = t_hs_ix;
            if (par_ix.hd == 1 &&
                mdl->sub_lp_hs[t_mix] > mdl->sub_lp_hs[b1_mix])
                b_h1_ix = t_hs_ix;
            if (par_ix.hd == 2 &&
                mdl->sub_lp_hs[t_mix] > mdl->sub_lp_hs[b2_mix])
                b_h2_ix = t_hs_ix;
        }
        mdl->sub_best_hs[bc_i] = b_hs_ix;
        mdl->sub_best_h1[bc_i] = b_h1_ix;
        mdl->sub_best_h2[bc_i] = b_h2_ix;
    }

    return 0;
}

int mdl_sub_pp_h(mdl_t *mdl) {
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_sub_pp_h: argument is null");

    if (mdl->sub_lp_h == NULL)
        return err_msg(-1, 0, "mdl_sub_pp_h: model must be initialized");

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    if (bd == NULL)
        return err_msg(-1, 0, "mdl_sub_pp_h: barcode data is not initialized");

    uint32_t D = mdl->mp->D;
    free(mdl->sub_pp_h);
    mdl->sub_pp_h = malloc(3 * D * sizeof(f_t));
    if (mdl->sub_pp_h == NULL)
        return err_msg(-1, 0, "mdl_sub_pp_h: %s", strerror(errno));

    uint32_t bc_i;
    for (bc_i = 0; bc_i < D; ++bc_i){
        if (bflg_get(&bd->absent_bc, bc_i))
            continue;
        // calculate posterior prob
        unsigned h;
        for (h = 0; h < 3; ++h) {
            f_t pp_h = exp(mdl->sub_lp_h[CMI(h, bc_i, 3)] - mdl->sub_lp_x[bc_i]);
            mdl->sub_pp_h[CMI(h, bc_i, 3)] = pp_h;
        }
    }

    return 0;
}

int mdl_fit(bam_data_t *bam_dat, obj_pars *objs){
    if (bam_dat == NULL || objs == NULL)
        return err_msg(-1, 0, "mdl_fit: arguments are NULL");

    if (bam_dat->bc_data == NULL)
        return err_msg(-1, 0, "mdl_fit: rna and atac are NULL");

    if (objs->gv == NULL || objs->vcf_hdr == NULL)
        return err_msg(-1, 0, "mdl_fit: variant data  is NULL");

    mdl_t *mdl = mdl_alloc();
    if (mdl == NULL)
        return -1;

    mdl->eps = objs->eps;
    mdl->alpha_eps = objs->alpha_eps;
    mdl->max_iter = objs->max_iter;
    mdl->threads = objs->threads;

    if (mdl_set_rna_atac(mdl, bam_dat->has_rna, bam_dat->has_atac) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    if (objs->verbose)
        log_msg("initializing data");
    if (mdl_bc_dat_bam_data(mdl->mdl_bc_dat, bam_dat, objs) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    if (mdl_set_samples(mdl, objs->vcf_hdr) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    if (mdl_set_hs_ix(mdl) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    if (mdl_alloc_probs(mdl) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    // initialize parameters
    if (objs->verbose)
        log_msg("initializing parameters");
    if (mdl_pars_set_dat(mdl->mp, mdl->mdl_bc_dat, objs,
                mdl->samples->n) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    int pi_amb_iter;
    if (objs->verbose)
        log_msg("calculating ambient sample fractions");
    if ((pi_amb_iter = mdl_m_pi_amb(mdl)) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    // run sub model
    if (mdl_em(mdl, objs) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    // get best (H,S) index
    if (mdl_sub_best_hs(mdl) < 0) {
        mdl_dstry(mdl);
        return -1;
    }
    // get posterior prob for (H)
    if (mdl_sub_pp_h(mdl) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    if (objs->verbose)
        log_msg("getting barcode likelihoods");

    // save output
    if (objs->verbose)
        log_msg("writing out likelihoods");
    if (write_sub_llk(mdl, objs->out_fn) < 0) {
        mdl_dstry(mdl);
        return -1;
    }
    if (write_samples(mdl, objs->out_fn) < 0) {
        mdl_dstry(mdl);
        return -1;
    }
    if (objs->verbose)
        log_msg("writing out parameters");
    if (write_lambda(mdl, objs->out_fn) < 0) {
        mdl_dstry(mdl);
        return -1;
    } 
    if (write_pi(mdl, objs->out_fn) < 0) {
        mdl_dstry(mdl);
        return -1;
    }
    if (bam_dat->has_rna && write_alpha_rna(mdl, objs->out_fn) < 0) {
        mdl_dstry(mdl);
        return -1;
    }
    if (bam_dat->has_atac && write_alpha_atac(mdl, objs->out_fn) < 0) {
        mdl_dstry(mdl);
        return -1;
    }
    if (write_alpha_prior(mdl, objs->out_fn) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    if (objs->verbose)
        log_msg("writing out summary results");
    if (write_res(mdl, bam_dat, objs->out_fn) < 0) {
        mdl_dstry(mdl);
        return -1;
    }

    mdl_dstry(mdl);

    return(0);
}

int write_lambda(mdl_t *mdl, char *fn){
    if (mdl == NULL || fn == NULL)
        return err_msg(-1, 0, "write_lambda: arguments are NULL");

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "write_lambda: model hasn't been initialized");

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    char *lambda_fn = ".lambda.txt.gz";
    char *out_lambda_fn = strcat2((const char*)fn, (const char*)lambda_fn);

    // row names
    char **row_names = malloc(3 * sizeof(char *));
    row_names[0] = strdup("0");
    row_names[1] = strdup("1");
    row_names[2] = strdup("2");
    if (row_names[0] == NULL || row_names[1] == NULL || row_names[2] == NULL)
        return err_msg(-1, 0, "write_lambda: %s", strerror(errno));

    // col names
    char **col_names = malloc(sizeof(char *));
    col_names[0] = strdup("Lambda");
    if (col_names[0] == NULL)
        return err_msg(-1, 0, "write_lambda: %s", strerror(errno));
    
    // write matrix
    ret = write_matrix_double(out_lambda_fn, mdl->mp->lambda, NULL, NULL, NULL, 
            row_names, 3, col_names, 1, 
            delim, nl, decp);
    if (ret < 0)
        return err_msg(-1, 0, "write_lambda: failed to write matrix to file");

    free(col_names[0]);
    free(col_names);
    free(row_names[0]);
    free(row_names[1]);
    free(row_names[2]);
    free(row_names);
    free(out_lambda_fn);

    return(0);
}

int write_pi(mdl_t *mdl, char *fn) {
    if (mdl == NULL || fn == NULL)
        return err_msg(-1, 0, "write_pi: arguments are NULL");

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "write_pi: model hasn't been initialized");

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;
    uint16_t M = mdl->mp->M;

    char *pi_fn = ".pi.txt.gz";
    char *out_pi_fn = strcat2((const char*)fn, (const char*)pi_fn);

    // row names
    char **row_names = malloc(M * sizeof(char *));
    unsigned int i;
    for (i = 0; i < M; ++i) {
        row_names[i] = str_map_str(mdl->samples, i);
    }

    // col names
    char **col_names = malloc(2 * sizeof(char *));
    col_names[0] = strdup("Pi");
    col_names[1] = strdup("Pi_amb");
    if (col_names[0] == NULL || col_names[1] == NULL)
        return err_msg(-1, 0, "write_pi: %s", strerror(errno));

    // pi array
    f_t *pi_all = malloc(sizeof(f_t) * M * 2);
    memcpy(pi_all, mdl->mp->pi, sizeof(f_t) * M);
    memcpy(pi_all + M, mdl->mp->pi_amb, sizeof(f_t) * M);
    
    // write matrix
    ret = write_matrix_double(out_pi_fn, pi_all, NULL, NULL, NULL, 
            row_names, M, col_names, 2, 
            delim, nl, decp);
    if (ret < 0)
        return err_msg(-1, 0, "write_pi: failed to write matrix to file");

    free(col_names[0]);
    free(col_names[1]);
    free(col_names);
    free(row_names);
    free(out_pi_fn);
    free(pi_all);

    return(0);
}

int write_alpha_rna(mdl_t *mdl, char *fn){
    if (mdl == NULL)
        return err_msg(-1, 0, "write_alpha_rna: arguments are NULL");

    if (mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "write_alpha_rna: no barcodes found");
    if (mdl->mp == NULL || mdl->mp->alpha_rna1 == NULL)
        return err_msg(-1, 0, "write_alpha_rna: model hasn't been initialized");
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_ixs = mdl->hs_ix->n_hs;
    int n_test_bc = mdl->mp->T;
    uint32_t n_all_bc = mdl->mp->D;

    if (n_test_bc != bd->test_bcs->n)
        return err_msg(-1, 0, "write_alpha_rna: number of test barcodes doesn't match");

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    char *alpha_fn = ".alpha_rna.txt.gz";
    char *out_alpha_fn = strcat2((const char*)fn, (const char*)alpha_fn);

    // row names
    char **bc_row_names = str_map_ca(bd->test_bcs);

    // alpha array for test barcodes
    f_t *al = malloc(n_test_bc * n_ixs * sizeof(f_t));
    int i;
    uint32_t j;
    for (i = 0; i < n_test_bc; ++i){
        char *bc = str_map_str(bd->test_bcs, i);
        int bc_ix = str_map_ix(bd->all_bcs, bc);
        for (j = 0; j < n_ixs; ++j) {
            f_t a = mdl->mp->alpha_rna1[CMI(bc_ix, j, n_all_bc)];
            al[CMI(i, j, n_test_bc)] = a;
        }
    }

    // write matrix
    ret = write_matrix_double(out_alpha_fn, al, NULL, NULL, NULL, 
            bc_row_names, n_test_bc, NULL, n_ixs, 
            delim, nl, decp);
    free(al);
    if (ret < 0)
        return err_msg(-1, 0, "write_alpha_rna: failed to write matrix to file");

    for (i = 0; i < n_test_bc; ++i)
        free(bc_row_names[i]);
    free(bc_row_names);
    free(out_alpha_fn);

    return 0;
}

int write_alpha_atac(mdl_t *mdl, char *fn){
    if (mdl == NULL)
        return err_msg(-1, 0, "write_alpha_atac: arguments are NULL");

    if (mdl->mp == NULL || mdl->mp->alpha_atac1 == NULL)
        return err_msg(-1, 0, "write_alpha_atac: model hasn't been initialized");
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_ixs = mdl->hs_ix->n_hs;
    int n_test_bc = mdl->mp->T;
    uint32_t n_all_bc = mdl->mp->D;

    if (n_test_bc != bd->test_bcs->n)
        return err_msg(-1, 0, "write_alpha_atac: number of test barcodes doesn't match");

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    char *alpha_fn = ".alpha_atac.txt.gz";
    char *out_alpha_fn = strcat2((const char*)fn, (const char*)alpha_fn);

    // row names
    char **bc_row_names = str_map_ca(bd->test_bcs);

    // alpha array for test barcodes
    f_t *al = malloc(n_test_bc * n_ixs * sizeof(f_t));
    int i;
    uint32_t j;
    for (i = 0; i < n_test_bc; ++i){
        char *bc = str_map_str(bd->test_bcs, i);
        int bc_ix = str_map_ix(bd->all_bcs, bc);
        for (j = 0; j < n_ixs; ++j) {
            f_t a = mdl->mp->alpha_atac1[CMI(bc_ix, j, n_all_bc)];
            al[CMI(i, j, n_test_bc)] = a;
        }
    }

    // write matrix
    ret = write_matrix_double(out_alpha_fn, al, NULL, NULL, NULL, 
            bc_row_names, n_test_bc, NULL, n_ixs, 
            delim, nl, decp);
    free(al);
    if (ret < 0)
        return err_msg(-1, 0, "write_alpha_atac: failed to write matrix to file");

    for (i = 0; i < n_test_bc; ++i)
        free(bc_row_names[i]);
    free(bc_row_names);
    free(out_alpha_fn);

    return 0;
}

int write_alpha_prior(mdl_t *mdl, char *fn) {
    if (mdl == NULL)
        return err_msg(-1, 0, "write_alpha_prior: arguments are NULL");

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    char *alpha_fn = ".alpha_prior.txt.gz";
    char *out_alpha_fn = strcat2((const char*)fn, (const char*)alpha_fn);
    if (out_alpha_fn == NULL)
        return -1;

    char *row_names[2] = {"prior_value", "prior_weight"};
    char *col_names[2] = {"RNA", "ATAC"};

    // alpha priors array
    f_t al[4] = {
        mdl->mp->alpha_prior_rna,
        mdl->mp->alpha_prior_w,
        mdl->mp->alpha_prior_atac,
        mdl->mp->alpha_prior_w
    };

    // write matrix
    ret = write_matrix_double(out_alpha_fn, al, NULL, NULL, NULL, 
                              row_names, 2, col_names, 2, 
                              delim, nl, decp);
    if (ret < 0) {
        free(out_alpha_fn);
        return err_msg(-1, 0, "write_alpha_prior: failed to write matrix to file");
    }

    free(out_alpha_fn);

    return 0;
}

int write_sub_llk(mdl_t *mdl, char *fn) {
    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    uint32_t i, j;
    int ret = 0;
    char *llk_fn = ".sample_llk.txt.gz";

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    uint32_t T = mdl->mp->T;

    char *out_llk_fn = strcat2((const char*)fn, (const char*)llk_fn);
    if (out_llk_fn == NULL)
        return err_msg(-1, 0, "write_sub_llk: %s", strerror(errno));

    // row names
    char **bc_row_names = str_map_ca(bd->test_bcs);
    if (bc_row_names == NULL)
        return err_msg(-1, 0, "write_sub_llk: failed to write to file");

    // get llk matrix
    f_t *llks = malloc(T * n_hs * sizeof(f_t));
    for (i = 0; i < T; ++i) {
        char *bc_name = str_map_str(bd->test_bcs, i);
        int all_bc_ix = str_map_ix(bd->all_bcs, bc_name);
        if (all_bc_ix < 0)
            return err_msg(-1, 0, "write_sub_llk: invalid bc index %d", all_bc_ix);

        for (j = 0; j < n_hs; ++j) {
            llks[CMI(i, j, T)] = mdl->sub_lp_hs[CMI(j, all_bc_ix, n_hs)];
        }
    }

    // write matrix
    ret = write_matrix_double(out_llk_fn, llks, NULL, NULL, NULL, 
            bc_row_names, T, NULL, n_hs, 
            delim, nl, decp);
    if (ret < 0)
        return err_msg(-1, 0, "write_sub_llk: failed to write matrix to file");

    for (i = 0; i < T; ++i)
        free(bc_row_names[i]);
    free(bc_row_names);
    free(llks);
    free(out_llk_fn);

    return 0;
}

int write_samples(mdl_t *mdl, char *fn){
    char *s_fn = ".samples.txt";
    int delim = '\t';
    int nl = '\n';

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "write_samples: failed to create output directory for %s", fn);

    char *out_s_fn = strcat2((const char*)fn, (const char*)s_fn);
    if (out_s_fn == NULL)
        return -1;
    FILE *fp = fopen(out_s_fn, "w");
    if (fp == NULL){
        err_msg(-1, 0, "write_samples: Could not open file %s", out_s_fn);
        free(out_s_fn);
        return -1;
    }

    const char *na_str = "NA";
    uint32_t n_hs = mdl->hs_ix->n_hs;
    par_ix_t par_ix;
    par_ix_init(&par_ix);
    int fret;
    uint32_t hs_ix;
    for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
        if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
            return -1;

        // if singlet
        char *h_str = aloc_sprintf("%i", par_ix.hd);
        if (h_str == NULL)
            return -1;
        const char *s1_str = NULL, *s2_str = NULL;
        switch (par_ix.hd) {
            case 0:
                s1_str = na_str;
                s2_str = na_str;
                break;
            case 1:
                s1_str = mdl->samples->strs[par_ix.s1];
                s2_str = na_str;
                break;
            case 2:
                s1_str = mdl->samples->strs[par_ix.s1];
                s2_str = mdl->samples->strs[par_ix.s2];
                break;
            default:
                return err_msg(-1, 0, "write_samples: invalid par H=%i", par_ix.hd);
        }

        fret = fputs(h_str, fp);
        fret = fputc(delim, fp);
        fret = fputs(s1_str, fp);
        fret = fputc(delim, fp);
        fret = fputs(s2_str, fp);
        fret = fputc(nl, fp);
        free(h_str);
    }

    // close file
    if ( (fret = fclose(fp)) != 0 )
        return err_msg(-1, 0, "write_samples: could not close file %s: %s", 
                out_s_fn, strerror(errno));

    free(out_s_fn);

    return 0;
}

int write_res(mdl_t *mdl, bam_data_t *bam_dat, char *fn){
    char *r_fn = ".summary.txt";
    unsigned int decp = 8;
    int fret;
    int delim = '\t';
    int nl = '\n';
    const char *nastr = "NA";

    if (mdl == NULL || bam_dat == NULL)
        return err_msg(-1, 0, "write_res: argument is null");

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_test_bc = mdl->mp->T;
    uint32_t n_all_bc = mdl->mp->D;
    uint32_t n_vars = bd->V;
    str_map *samples = mdl->samples;
    str_map *test_bcs = bd->test_bcs;
    par_ix_t par_ix;
    par_ix_init(&par_ix);

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "write_res: failed to create output directory for %s", fn);

    char *out_r_fn = strcat2((const char*)fn, (const char*)r_fn);
    if (out_r_fn == NULL)
        return -1;
    FILE *fp = fopen(out_r_fn, "w");
    if (fp == NULL){
        err_msg(-1, 0, "write_res: Could not open file %s", out_r_fn);
        free(out_r_fn);
        return -1;
    }

    char hdr[] = "Barcode\tn_rna_molecules\tn_atac_molecules\tn_features\tn_peaks"
        "\tn_rna_info\tn_atac_info"
        "\tn_rna_variants\tn_atac_variants\trna_pct_mt\tatac_pct_mt\tFRIG\tFRIP"
        "\tbest_type\tbest_sample"
        "\tbest_rna_ambient\trna_ambient_info"
        "\tbest_atac_ambient\tatac_ambient_info"
        "\tbest_singlet\tbest_doublet"
        "\tPP0\tPP1\tPP2"
        "\tLLK0\tLLK1\tLLK2\n";

    fputs(hdr, fp);

    size_t buf_size = decp + 1000;
    char *pstr = (char *)malloc(buf_size * sizeof(char));
    if (pstr == NULL)
        return err_msg(-1, 0, "write_res: %s", strerror(errno));

    int pstr_len;
    uint32_t bc_i, i;
    for (bc_i = 0; bc_i < n_test_bc; ++bc_i){
        char *bc_name = str_map_str(test_bcs, bc_i);
        fret = fputs(bc_name, fp);
        int all_bc_ix = str_map_ix(bd->all_bcs, bc_name);
        if (all_bc_ix < 0)
            return err_msg(-1, 0, "write_res: invalid bc index %d", all_bc_ix);

        double par;

        khint_t k_bc = kh_get(kh_bc_dat, bam_dat->bc_data, bc_name);
        if (k_bc == kh_end(bam_dat->bc_data))
            return err_msg(-1, 0, "write_res: barcode %s not found", bc_name);

        bc_data_t *bc_data = kh_val(bam_dat->bc_data, k_bc);
        bc_stats_t *bcc = bc_data->bc_stats;
        mdl_mlcl_bc_t *mlcl_bc = &mv_i(&bd->bc_mlcl, all_bc_ix);

        // get number of informative reads and num. variants
        uint32_t n_rna_info = 0, n_atac_info = 0, n_rna_var = 0, n_atac_var = 0;
        if (mdl_mlcl_bc_info_count(mlcl_bc, n_vars, &n_rna_info, &n_atac_info,
                                   &n_rna_var, &n_atac_var) < 0)
            return -1;

        // write bc stats
        fputc(delim, fp);
        int2strp(bcc->rna_counts, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        int2strp(bcc->atac_counts, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        int2strp(bcc->n_gene, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        int2strp(bcc->n_peak, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        int2strp(n_rna_info, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        int2strp(n_atac_info, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        int2strp(n_rna_var, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        int2strp(n_atac_var, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        double2str_in((double)bcc->rna_mt, &pstr, &buf_size, 4);
        fputs(pstr, fp);

        fputc(delim, fp);
        double2str_in((double)bcc->atac_mt, &pstr, &buf_size, 4);
        fputs(pstr, fp);

        fputc(delim, fp);
        double2str_in((double)bcc->frig, &pstr, &buf_size, 4);
        fputs(pstr, fp);

        fputc(delim, fp);
        double2str_in((double)bcc->frip, &pstr, &buf_size, 4);
        fputs(pstr, fp);

        // get best sample
        f_t *pp_h = mdl->sub_pp_h + (3 * all_bc_ix);
        uint8_t best_h = 0;
        if (pp_h[1] > pp_h[2] && pp_h[1] > pp_h[0])
            best_h = 1;
        if (pp_h[2] > pp_h[1] && pp_h[2] > pp_h[0])
            best_h = 2;
        int32_t best_h1s = mdl->sub_best_h1[all_bc_ix];
        int32_t best_h2s = mdl->sub_best_h2[all_bc_ix];

        // write best H_d string
        fret = fputc(delim, fp);
        if (best_h == 0) {
            fret = fputs("Empty", fp);
        } else if (best_h == 1) {
            fret = fputs("Singlet", fp);
        } else {
            fret = fputs("Doublet", fp);
        }

        // get best S_d strings
        fret = fputc(delim, fp);
        if (best_h == 0) {
            fputs("Empty", fp);
        } else if (best_h == 1) {
            if (hs_ix_get_pars(mdl->hs_ix, best_h1s, &par_ix) < 0)
                return -1;
            fputs(str_map_str(samples, par_ix.s1), fp);
        } else {
            if (hs_ix_get_pars(mdl->hs_ix, best_h2s, &par_ix) < 0)
                return -1;
            fputs(str_map_str(samples, par_ix.s1), fp);
            fputs(":", fp);
            fputs(str_map_str(samples, par_ix.s2), fp);
        }

        // write alphas
        int mod_ix;
        f_t *mod_alpha[2] = {NULL, NULL};
        mod_alpha[0] = mdl->mp->alpha_rna1;
        mod_alpha[1] = mdl->mp->alpha_atac1;
        f_t *mod_alpha_info[2] = {NULL, NULL};
        mod_alpha_info[0] = mdl->mp->alpha_rna_info;
        mod_alpha_info[1] = mdl->mp->alpha_atac_info;
        for (mod_ix = 0; mod_ix < 2; ++mod_ix) {
            f_t alpha_v, alpha_info;
            if (best_h == 0) {
                alpha_v = NAN;
                alpha_info = NAN;
            } else if (best_h == 1) {
                alpha_v = mod_alpha[mod_ix][CMI(all_bc_ix, best_h1s, n_all_bc)];
                alpha_info = mod_alpha_info[mod_ix][CMI(all_bc_ix, best_h1s, n_all_bc)];
            } else {
                alpha_v = mod_alpha[mod_ix][CMI(all_bc_ix, best_h2s, n_all_bc)];
                alpha_info = mod_alpha_info[mod_ix][CMI(all_bc_ix, best_h2s, n_all_bc)];
            }

            // alpha estimate
            if (isnan(alpha_v)) {
                fret = fputc(delim, fp);
                fret = fputs(nastr, fp);
            } else {
                par = alpha_v;
                if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
                    return err_msg(-1, 0, "write_res: failed to convert %e to string", par);
                fret = fputc(delim, fp);
                fret = fputs(pstr, fp);
            }
            // alpha info
            if (isnan(alpha_info)) {
                fret = fputc(delim, fp);
                fret = fputs(nastr, fp);
            } else {
                par = alpha_info;
                if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
                    return err_msg(-1, 0, "write_res: failed to convert %e to string", par);
                fret = fputc(delim, fp);
                fret = fputs(pstr, fp);
            }
        }

        // write best singlet
        if (hs_ix_get_pars(mdl->hs_ix, best_h1s, &par_ix) < 0)
            return -1;
        fret = fputc(delim, fp);
        fret = fputs(str_map_str(samples, par_ix.s1), fp);

        // write best doublet
        if (hs_ix_get_pars(mdl->hs_ix, best_h2s, &par_ix) < 0)
            return -1;
        fret = fputc(delim, fp);
        fret = fputs(str_map_str(samples, par_ix.s1), fp);
        fret = fputs(":", fp);
        fret = fputs(str_map_str(samples, par_ix.s2), fp);

        // write posterior probs of H_d
        for (i = 0; i < 3; ++i){
            par = mdl->sub_pp_h[CMI(i, all_bc_ix, 3)];
            if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
                return err_msg(-1, 0, "write_res: failed to convert %e to string", par);
            fret = fputc(delim, fp);
            fret = fputs(pstr, fp);
        }

        // write llks of H_d
        for (i = 0; i < 3; ++i){
            par = mdl->sub_lp_h[CMI(i, all_bc_ix, 3)];
            if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
                return err_msg(-1, 0, "write_res: failed to convert %e to string", par);
            fret = fputc(delim, fp);
            fret = fputs(pstr, fp);
        }

        fret = fputc(nl, fp);
    }

    free(pstr);

    // close file
    if ( (fret = fclose(fp)) != 0 )
        return err_msg(-1, 0, "write_res: could not close file %s: %s", 
                out_r_fn, strerror(errno));
    free(out_r_fn);

    return 0;
}

