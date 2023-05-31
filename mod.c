
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
#include <assert.h>
#include "bits.h"

#define feps (FLT_EPSILON * 100e3)
#define ATAC_IX 0
#define RNA_IX 1

/*******************************************************************************
 * num check
 ******************************************************************************/

int num_invalid(f_t x) {
    if (isnan(x))
        return 1;
    if (isinf(x))
        return 1;
    return 0;
}

int prob_invalid(f_t x) {
    if (num_invalid(x))
        return 1;
    if (x < 0 || x > 1)
        return 1;
    return 0;
}

int psum_invalid(f_t x, f_t lb, f_t ub) {
    if (num_invalid(x))
        return 1;
    if (x < lb || x > ub)
        return 1;
    return 0;
}

/*******************************************************************************
 * mdl_mlcl_t
 ******************************************************************************/

void mdl_mlcl_init(mdl_mlcl_t *mlcl){
    mv_init(&mlcl->feat_ixs);
    mv_init(&mlcl->var_ixs);
    mlcl->counts = 0;
}
void mdl_mlcl_free(mdl_mlcl_t *mlcl){
    if (mlcl == NULL) return;
    mv_free(&mlcl->feat_ixs);
    mv_free(&mlcl->var_ixs);
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
        fprintf(f, ", %i", mv_i(&mlcl->var_ixs, i));
    }
    fprintf(f, ", count=%u", mlcl->counts);
    fprintf(f, "\n");
}

uint32_t mdl_mlcl_feat_count(kbtree_t(kb_mdl_mlcl) *bt){
    kbitr_t itr;
    uint32_t counts = 0;
    kb_itr_first(kb_mdl_mlcl, bt, &itr); 
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, bt, &itr)){
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        counts += mlcl->counts;
    }
    return(counts);
}

int mdl_mlcl_add_rna(mdl_mlcl_t *mlcl, rna_mol_t *mol,
        int n_genes, int n_vars) {
    size_t mol_n_genes = ml_size(&mol->gl);
    // add gene if only 1 is present, skip multi-gene RNAs
    if (mol_n_genes == 1){
        ml_node_t(seq_gene_l) *g_n;
        for (g_n = ml_begin(&mol->gl); g_n; g_n = ml_node_next(g_n)){
            seq_gene_t gene = ml_node_val(g_n);
            assert(gene.gene_id >= 0);
            assert(gene.splice < 3);
            int32_t f_ix = gene.gene_id + (n_genes * gene.splice);
            if (mv_push(i32, &mlcl->feat_ixs, f_ix) < 0)
                return(-1);
        }
    }

    ml_node_t(seq_vac_l) *v_n;
    for (v_n = ml_begin(&mol->vl); v_n; v_n = ml_node_next(v_n)){
        seq_vac_t vac = ml_node_val(v_n);
        assert(vac.vix >= 0);
        // set allele to 0:ref, 1:alt, 2:any other
        uint8_t allele = vac.allele < 2 ? vac.allele : 2;
        int32_t v_ix = vac.vix + (n_vars * allele);
        if (mv_push(i32, &mlcl->var_ixs, v_ix) < 0)
            return(-1);
    }
    return 0;
}

int mdl_mlcl_add_atac(mdl_mlcl_t *mlcl, atac_frag_t *frag, int n_vars) {
    // add peak
    int32_t pk_ix = mv_size(&frag->pks) > 0 ? 1 : 0;
    if (mv_push(i32, &mlcl->feat_ixs, pk_ix) < 0)
        return(-1);

    // variant(s)
    ml_node_t(seq_vac_l) *v_n;
    for (v_n = ml_begin(&frag->vl); v_n; v_n = ml_node_next(v_n)){
        seq_vac_t vac = ml_node_val(v_n);
        assert(vac.vix >= 0);
        // set allele to 0:ref, 1:alt, 2:any other
        uint8_t allele = vac.allele < 2 ? vac.allele : 2;
        int32_t v_ix = vac.vix + (n_vars * allele);
        if (mv_push(i32, &mlcl->var_ixs, v_ix) < 0)
            return(-1);
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
    f_t ret = (f_t)pks / (f_t)counts;
    return(ret);
}

void mdl_bc_counts(mdl_mlcl_bc_t *mdl_bc, uint32_t *rna, uint32_t *atac){
    *rna = mdl_mlcl_feat_count(mdl_bc->rna);
    *atac = mdl_mlcl_feat_count(mdl_bc->atac);
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
    kbtree_t(kb_mdl_mlcl) *p;
    p = kb_init(kb_mdl_mlcl, KB_DEFAULT_SIZE);
    if (p == NULL){
        return err_msg(-1, 0, "mdl_mlcl_bc_init: %s", strerror(errno));
    }
    mdl_bc->rna = p;

    p = kb_init(kb_mdl_mlcl, KB_DEFAULT_SIZE);
    if (p == NULL){
        return err_msg(-1, 0, "mdl_mlcl_bc_init: %s", strerror(errno));
    }
    mdl_bc->atac = p;

    mdl_bc->n_bc = 0;

    return 0;
}

void mdl_mlcl_bc_free(mdl_mlcl_bc_t *mdl_bc){
    if (mdl_bc == NULL) return;
    kbitr_t itr;
    kbtree_t(kb_mdl_mlcl) *bt;

    bt = mdl_bc->rna;
    kb_itr_first(kb_mdl_mlcl, bt, &itr); 
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, bt, &itr)){
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        mdl_mlcl_free(mlcl);
    }
    kb_destroy(kb_mdl_mlcl, mdl_bc->rna);
    mdl_bc->rna = NULL;

    bt = mdl_bc->atac;
    kb_itr_first(kb_mdl_mlcl, bt, &itr); 
    for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, bt, &itr)){
        mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
        mdl_mlcl_free(mlcl);
    }
    kb_destroy(kb_mdl_mlcl, mdl_bc->atac);
    mdl_bc->atac = NULL;

    mdl_bc->n_bc = 0;
}

/*******************************************************************************
 * index structs
 ******************************************************************************/

void par_ix_init(par_ix_t *par_ix) {
    if (par_ix == NULL)
        return;
    par_ix->ix = -1;
    par_ix->hs_ix = -1;
    par_ix->hd = -1;
    par_ix->s1 = -1;
    par_ix->s2 = -1;
    par_ix->k = -1;
    int i;
    for (i = 0; i < 3; ++i)
        par_ix->t_ix[i] = -1;
    par_ix->t_n = -1;
}

void hs_ix_init(hs_ix_t *hs_ix) {
    if (hs_ix == NULL)
        return;
    hs_ix->n_sam = 0;
    hs_ix->n_ixs = 0;
    hs_ix->n_hs = 0;
    hs_ix->k = 0;
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
    hs_ix->n_ixs = 0;
    hs_ix->n_hs = 0;
    hs_ix->k = 0;
    free(hs_ix->ix2pars);
    hs_ix->ix2pars = NULL;
}

void hs_ix_dstry(hs_ix_t *hs_ix) {
    if (hs_ix == NULL)
        return;
    hs_ix_free(hs_ix);
    free(hs_ix);
}

void hs_ix_set(hs_ix_t *hs_ix, uint16_t n_samples, uint16_t n_clust) {
    if (hs_ix == NULL)
        return;

    // indices
    hs_ix->n_sam = n_samples;
    hs_ix->k = n_clust;
    hs_ix->n_hs = 1 + n_samples + (n_samples * (n_samples - 1) / 2);
    hs_ix->n_ixs = 1 + (n_clust * n_samples) +
        ( n_clust * n_samples * (n_samples - 1) / 2);
    hs_ix->ix2pars = calloc(4 * hs_ix->n_ixs, sizeof(int));

    int s1 = 0, s2 = 1;
    uint32_t ixi = 0;
    hs_ix->ix2pars[CMI(ixi, 0, hs_ix->n_ixs)] = 0;
    hs_ix->ix2pars[CMI(ixi, 1, hs_ix->n_ixs)] = -1;
    hs_ix->ix2pars[CMI(ixi, 2, hs_ix->n_ixs)] = -1;
    hs_ix->ix2pars[CMI(ixi, 3, hs_ix->n_ixs)] = -1;
    ++ixi;

    uint16_t ki;
    for (ki = 0; ki < n_clust; ++ki) {
        for (s1 = 0; s1 < n_samples; ++s1){
            hs_ix->ix2pars[CMI(ixi, 0, hs_ix->n_ixs)] = 1;
            hs_ix->ix2pars[CMI(ixi, 1, hs_ix->n_ixs)] = s1;
            hs_ix->ix2pars[CMI(ixi, 2, hs_ix->n_ixs)] = -1;
            hs_ix->ix2pars[CMI(ixi, 3, hs_ix->n_ixs)] = ki;
            ++ixi;
        }
        for (s1 = 0; s1 < n_samples - 1; ++s1){
            for (s2 = s1 + 1; s2 < n_samples; ++s2){
                hs_ix->ix2pars[CMI(ixi, 0, hs_ix->n_ixs)] = 2;
                hs_ix->ix2pars[CMI(ixi, 1, hs_ix->n_ixs)] = s1;
                hs_ix->ix2pars[CMI(ixi, 2, hs_ix->n_ixs)] = s2;
                hs_ix->ix2pars[CMI(ixi, 3, hs_ix->n_ixs)] = ki;
                ++ixi;
            }
        }
    }
    assert(ixi == hs_ix->n_ixs);
}

int hs_ix_get_pars(hs_ix_t *hs_ix, int ix,
        par_ix_t *par_ix) {

    if (ix < 0)
        return err_msg(-1, 0, "ix=%i < 0", ix);
    if ((uint16_t)ix >= hs_ix->n_ixs)
        return err_msg(-1, 0, "ix=%i > n_ixs=%u", ix, hs_ix->n_ixs);

    // -1 for if NA/invalid
    par_ix->ix = ix;
    uint32_t divby = hs_ix->n_hs - 1;
    if (ix == 0)
        par_ix->hs_ix = ix;
    else
        par_ix->hs_ix = ((ix - 1) % divby) + 1;
    assert((uint32_t)par_ix->hs_ix < hs_ix->n_hs);
    par_ix->hd = hs_ix->ix2pars[CMI(ix, 0, hs_ix->n_ixs)];
    par_ix->s1 = hs_ix->ix2pars[CMI(ix, 1, hs_ix->n_ixs)];
    par_ix->s2 = hs_ix->ix2pars[CMI(ix, 2, hs_ix->n_ixs)];
    par_ix->k = hs_ix->ix2pars[CMI(ix, 3, hs_ix->n_ixs)];

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
    mp->G = 0;
    mp->V = 0;
    mp->M = 0;
    mp->M = 0;
    mp->n_hs = 0;
    mp->n_ix = 0;

    int err;
    if ((err = pthread_mutex_init(&mp->sum_lock, NULL)) != 0)
        return err_msg(-1, 0, "mdl_pars_init: failed to initialize mutex %i", err);
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
    free(mp->kappa);
    free(mp->alpha_rna);
    free(mp->alpha_atac);
    free(mp->rho);
    free(mp->gamma);
    
    pthread_mutex_destroy(&mp->sum_lock);
    free(mp->_pi_sum);
    free(mp->_kappa_sum);
    free(mp->_rho_sum);
}

void mdl_pars_dstry(mdl_pars_t *mp){
    if (mp == NULL) return;
    mdl_pars_free(mp);
    free(mp);
}

int mld_pars_set_num_alloc(mdl_pars_t *mp, uint32_t D, uint32_t T,
        uint32_t G, uint32_t V, uint16_t M, uint16_t K) {
    if (mp == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: argument is null");

    mp->D = D;
    mp->T = T;
    mp->G = G;
    mp->V = V;
    mp->M = M;
    mp->K = K;
    mp->n_ix = 1 + (mp->K * mp->M) +
        ( mp->K * mp->M * (mp->M - 1) / 2);

    mp->n_hs = 1 + (mp->M) + (mp->M * (mp->M - 1) / 2);

    // allocate parameter fields
    mp->pi = calloc(mp->M, sizeof(f_t));
    if (mp->pi == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->kappa = calloc(mp->K, sizeof(f_t));
    if (mp->kappa == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->alpha_rna = calloc(D * mp->n_hs, sizeof(f_t));
    if (mp->alpha_rna == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->alpha_atac = calloc(D * mp->n_hs, sizeof(f_t));
    if (mp->alpha_atac == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->rho = calloc((mp->G * 3) * (mp->K + 1), sizeof(f_t));
    if (mp->rho == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->gamma = calloc((mp->M + 1) * mp->V, sizeof(f_t));
    if (mp->gamma == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    // allocate sum fields
    mp->_pi_sum = calloc(mp->M, sizeof(f_t));
    if (mp->_pi_sum == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->_kappa_sum = calloc(mp->K, sizeof(f_t));
    if (mp->_kappa_sum == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->_rho_sum = calloc((mp->G * 3) * (mp->K + 1), sizeof(f_t));
    if (mp->_rho_sum == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->_par_diff = 0;

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

    for (i = 0; i < mp->K; ++i)
        mp->_kappa_sum[i] = psc;

    uint32_t rho_ne = (mp->G * 3) * (mp->K + 1);
    for (i = 0; i < rho_ne; ++i)
        mp->_rho_sum[i] = psc;

    for (i = 0; i < 4; ++i)
        mp->_sigma_sum[i] = psc;

    return 0;
}

int mdl_pars_pi_fix(mdl_pars_t *mp){

    uint32_t i, j;
    uint16_t M = mp->M;

    // maximize pi (keep fixed uniform for now)
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

int mdl_pars_kappa_uni(mdl_pars_t *mp) {
    if (mp == NULL)
        return err_msg(-1, 0, "mdl_pars_kappa_uni: argument is null");

    uint16_t i;
    f_t ku = 1.0 / (f_t)(mp->K);
    for (i = 0; i < mp->K; ++i)
        mp->kappa[i] = ku;
    return 0;
}

int mdl_pars_add_gamma(mdl_pars_t *mp, float **a, int nv, int ns){
    if (mp == NULL || a == NULL)
        return err_msg(-1, 0, "mdl_pars_add_gamma: a is NULL");

    if (mp->M != (uint16_t)ns)
        return err_msg(-1, 0, "mdl_pars_add_gamma: "
                "ns and number of samples in a don't match");

    if (nv < 1 || ns < 1)
        return err_msg(-1, 0, "mdl_pars_add_gamma: nv=%i and ns=%i "
                "must be non-negative", nv, ns);

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
            mp->gamma[CMI(s, v, ns1)] = (f_t)av;
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
        mp->gamma[CMI(s, v, ns1)] = (f_t)ambp;
    }

    return 0;
}

int mdl_pars_set_gamma_amb(mdl_pars_t *mp){

    uint32_t i, M = mp->M, V = mp->V, M1 = M + 1;

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
                ap += mp->pi[st] * mp->gamma[CMI(st, i, M1)];
                pi_tot += mp->pi[st];
            }
        }
        if (n_miss == M)
            mp->gamma[CMI(st, i, M1)] = -1; // if all missing, set missing
        else
            mp->gamma[CMI(st, i, M1)] = ap / pi_tot;
    }
    return(0);
}

int mdl_pars_set_dat(mdl_pars_t *mp, mdl_bc_dat_t *bd, obj_pars *objs,
        uint16_t n_samples, uint16_t K){
    if (mp == NULL || bd == NULL)
        return err_msg(-1, 0, "mdl_pars_set_dat: argument is null");

    if (bd->all_bcs->n < 1)
        return err_msg(-1, 0, "mdl_pars_set_dat: no barcodes in bc_dat");
    if (bd->test_bcs->n < 1)
        return err_msg(-1, 0, "mdl_pars_set_dat: no test barcodes in bc_dat");

    // warn if no variants are present
    if (bd->V < 1)
        err_msg(0, 1, "mdl_pars_set_dat: no variants present (%u) "
                "initializing parameters", bd->V);

    // set numbers and allocate arrays
    if (mld_pars_set_num_alloc(mp, bd->all_bcs->n, bd->test_bcs->n,
                bd->G, bd->V, n_samples, K) < 0)
        return -1;

    uint32_t rho_nrow = 3 * mp->G;
    uint32_t rho_ncol = 1 + mp->K;
    uint32_t sig_nrow = 2;

    // keep pi uniform
    if (mdl_pars_pi_fix(mp) < 0)
        return(-1);

    // set kappa uniform
    if (mdl_pars_kappa_uni(mp) < 0)
        return -1;

    // arrays to sum reads
    uint32_t rho_ne = rho_nrow * rho_ncol;
    f_t *rho_sum = calloc(rho_ne, sizeof(f_t));
    uint32_t sig_ne = sig_nrow * 2;
    f_t *sig_sum = calloc(sig_ne, sizeof(f_t)); // peak by droptype matrix
    if (rho_sum == NULL || sig_sum == NULL)
        return err_msg(-1, 0, "mdl_pars_set_dat: %s", strerror(errno));

    mp->tau = TAU;
    uint32_t i, j;

    f_t psc = 1;
    for (i = 0; i < rho_ne; ++i)
        rho_sum[i] = psc;
    for (i = 0; i < sig_ne; ++i)
        sig_sum[i] = psc;

    f_t *rho_tot = calloc(rho_ncol, sizeof(f_t));
    for (i = 0; i < rho_ncol; ++i)
        rho_tot[i] = psc * rho_nrow;
    f_t sig_tot[2] = {psc * sig_nrow, psc * sig_nrow};

    f_t n_amb_bcs = 1; // prior of 1 for lambda
    f_t n_nuc_bcs = 1;

    for (i = 0; i < mp->D; ++i){
        int fl;

        // if no reads
        fl = bflg_get(bd->absent_bc, i);
        if (fl == 1) continue;

        // if ambient
        fl = bflg_get(bd->amb_flag, i);
        uint32_t c_ix = fl ? 0 : 1;

        // sample cell type
        int rand_k = rand() % mp->K;
        uint32_t rho_col = fl ? 0 : rand_k + 1;

        // bam data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, i);
        kbtree_t(kb_mdl_mlcl) *mols = bc_dat.rna;
        kbtree_t(kb_mdl_mlcl) *frags = bc_dat.atac;

        if (fl)
            n_amb_bcs += bc_dat.n_bc;
        else
            n_nuc_bcs += bc_dat.n_bc;

        // loop over RNA
        kbitr_t itr;
        kb_itr_first(kb_mdl_mlcl, mols, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            assert(mlcl);
            size_t f_ix;
            for (f_ix = 0; f_ix < mv_size(&mlcl->feat_ixs); ++f_ix){
                uint32_t rho_row = mv_i(&mlcl->feat_ixs, f_ix);
                uint32_t rho_ix = CMI(rho_row, rho_col, rho_nrow);
                rho_sum[rho_ix] += (f_t)mlcl->counts;
                rho_tot[rho_col] += (f_t)mlcl->counts;
            }
        }

        // loop over ATAC
        kb_itr_first(kb_mdl_mlcl, frags, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, frags, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            assert(mlcl);
            size_t f_ix;
            assert(mv_size(&mlcl->feat_ixs) == 1);
            for (f_ix = 0; f_ix < mv_size(&mlcl->feat_ixs); ++f_ix){
                uint32_t pk = (uint32_t)(mv_i(&mlcl->feat_ixs, f_ix));
                uint32_t six = CMI(pk, c_ix, sig_nrow);
                sig_sum[six] += (f_t)mlcl->counts;
                sig_tot[c_ix] += (f_t)mlcl->counts;
            }
        }

        // set alpha
        for (j = 0; j < mp->n_hs; ++j) {
            uint32_t aix = CMI(i, j, mp->D);
            if (fl) {
                mp->alpha_rna[aix] = 1;
                mp->alpha_atac[aix] = 1;
            } else {
                mp->alpha_rna[aix] = 0.5;
                mp->alpha_atac[aix] = 0.5;
            }
        }
    }

    // set lambda
    f_t n_tot_bcs = n_amb_bcs + n_nuc_bcs;
    f_t n_sng = n_nuc_bcs * 0.9;
    f_t n_dbl = n_nuc_bcs * 0.1;
    mp->lambda[0] = n_amb_bcs / n_tot_bcs;
    mp->lambda[1] = n_sng / n_tot_bcs;
    mp->lambda[2] = n_dbl / n_tot_bcs;

    // set rho
    for (i = 0; i < rho_ncol; ++i) {
        for (j = 0; j < rho_nrow; ++j) {
            uint32_t rho_ix = CMI(j, i, rho_nrow);
            mp->rho[rho_ix] = rho_sum[rho_ix] / rho_tot[i];
        }
    }

    // set sigma
    uint32_t six0 = CMI(0, 0, sig_nrow);
    uint32_t six1 = CMI(0, 1, sig_nrow);
    mp->sigma[0] = sig_sum[six0] / sig_tot[0];
    mp->sigma[1] = sig_sum[six1] / sig_tot[1];

    // set gamma from VCF
    // first index: variant, second index: sample
    int32_t *vixs = calloc(bd->V, sizeof(int32_t));
    for (i = 0; i < bd->V; ++i) vixs[i] = i;
    float **gm = ap_array_gt(objs->gv, objs->vcf_hdr, vixs, bd->V, "GT");
    free(vixs);
    if (gm == NULL)
        return(-1);

    if (mdl_pars_add_gamma(mp, gm, bd->V, n_samples) < 0)
        return -1;

    // set tau
    mp->tau = TAU;

    // free
    for (i = 0; i < bd->V; ++i){
        free(gm[i]);
    }
    free(gm);
    free(rho_sum);
    free(rho_tot);
    free(sig_sum);
    
    return(0);
}

// TODO: check pi_sum
int mdl_pars_check(mdl_pars_t *mp){

    f_t lt1 = 1 - 1e-8;
    f_t ut1 = 1 + 1e-8;

    uint32_t i, j, M1 = mp->M + 1;

    // check lambda
    f_t lambda_tot = 0;
    for (i = 0; i < 3; ++i) {
        if (prob_invalid(mp->lambda[i]))
            return err_msg(-1, 0, "mdl_pars_check: lambda[%i] = %f",
                    i, mp->lambda[i]);
        lambda_tot += mp->lambda[i];
    }
    if (psum_invalid(lambda_tot, lt1, ut1))
        return err_msg(-1, 0, "mdl_pars_check: lambda sum = %f", lambda_tot);

    // check pi
    f_t pi_sum = 0;
    for (i = 0; i < mp->M; ++i) {
        if (prob_invalid(mp->pi[i]))
            return err_msg(-1, 0, "mdl_pars_check: pi[%i] = %f",
                    i, mp->pi[i]);
        pi_sum += mp->pi[i];
    }
    if (psum_invalid(pi_sum, lt1, ut1))
        return err_msg(-1, 0, "mdl_pars_check: pi sum = %f", pi_sum);

    // kappa
    f_t kappa_sum = 0;
    for (i = 0; i < mp->K; ++i) {
        if (prob_invalid(mp->kappa[i]))
            return err_msg(-1, 0, "mdl_pars_check: kappa[%i] = %f",
                    i, mp->kappa[i]);
        kappa_sum += mp->kappa[i];
    }
    if (psum_invalid(kappa_sum, lt1, ut1))
        return err_msg(-1, 0, "mdl_pars_check: kappa sum = %f", kappa_sum);

    // alpha
    for (i = 0; i < mp->D; ++i){
        for (j = 0; j < mp->n_hs; ++j) {
            f_t a;
            a = mp->alpha_rna[CMI(i, j, mp->D)];
            if (prob_invalid(a))
                return err_msg(-1, 0, "mdl_pars_check: alpha_rna[%i,%i] = %f",
                        i, j, a);
            a = mp->alpha_atac[CMI(i, j, mp->D)];
            if (prob_invalid(a))
                return err_msg(-1, 0, "mdl_pars_check: alpha_atac[%i,%i] = %f",
                        i, j, a);
        }
    }

    // rho
    uint32_t rho_nrow = mp->G * 3;
    uint32_t rho_ncol = 1 + mp->K;
    f_t *rho_tot = calloc(rho_ncol, sizeof(f_t));
    for (i = 0; i < rho_ncol; ++i)
        rho_tot[i] = 0;
    for (i = 0; i < rho_ncol; ++i){
        for (j = 0; j < rho_nrow; ++j){
            rho_tot[i] += mp->rho[CMI(j, i, rho_nrow)];
            f_t p = mp->rho[CMI(j, i, rho_nrow)];
            if (prob_invalid(p))
                return err_msg(-1, 0, "mdl_pars_check: rho[%i,%i] = %f",
                        i, j, p);
        }
    }
    for (i = 0; i < 2; ++i){
        if (psum_invalid(rho_tot[i], lt1, ut1))
            return err_msg(-1, 0, "mdl_pars_check: rho sum %i = %f",
                    i, rho_tot[i]);
    }

    // sigma
    for (i = 0; i < 2; ++i){
        if (prob_invalid(mp->sigma[i]))
            return err_msg(-1, 0, "mdl_pars_check: sigma[%i] = %f",
                    i, mp->sigma[i]);
    }

    // gamma
    for (i = 0; i < M1; ++i) {
        for (j = 0; j < mp->V; ++j) {
            f_t p = mp->gamma[CMI(i, j, M1)];
            if (p > (-1 - 1e-8) && p < (1 + 1e-8))
                continue;
            if (prob_invalid(p))
                return err_msg(-1, 0, "mdl_pars_check: gamma[%i, %i] = %f",
                        i, j, p);
        }
    }

    // tau
    if (prob_invalid(mp->tau))
        return err_msg(-1, 0, "mdl_pars_check: tau = %f",
                mp->tau);

    free(rho_tot);

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

    mdl_bc_dat->amb_flag = calloc(1, sizeof(bflg_t));
    mdl_bc_dat->absent_bc = calloc(1, sizeof(bflg_t));
    if (mdl_bc_dat->absent_bc == NULL)
        return err_msg(-1, 0, "mdl_bc_dat_init: %s", strerror(errno));
    bflg_init_empty(mdl_bc_dat->amb_flag);
    bflg_init_empty(mdl_bc_dat->absent_bc);
    mdl_bc_dat->G = 0;
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

    int i;
    for (i = 0; i < mdl_bc_dat->all_bcs->n; ++i)
        mdl_mlcl_bc_free(&mv_i(&mdl_bc_dat->bc_mlcl, i));

    destroy_str_map(mdl_bc_dat->all_bcs);
    destroy_str_map(mdl_bc_dat->test_bcs);

    bflg_free(mdl_bc_dat->amb_flag);
    bflg_free(mdl_bc_dat->absent_bc);
    free(mdl_bc_dat->amb_flag);
    free(mdl_bc_dat->absent_bc);
    mdl_bc_dat->G = 0;
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

    if ((!dup_ok) && *found)
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

    if (bflg_resize(mdl_bc_dat->absent_bc, mdl_bc_dat->all_bcs->n) < 0)
        return(-1);
    if (bflg_resize(mdl_bc_dat->amb_flag, mdl_bc_dat->all_bcs->n) < 0)
        return(-1);

    // set to absent if no reads/empty
    if (absent) {
        bflg_set(mdl_bc_dat->absent_bc, bc_ix);
    } else {
        bflg_unset(mdl_bc_dat->absent_bc, bc_ix);
    }

    // set ambient flag and add to test_bcs
    if (ambient) {
        bflg_set(mdl_bc_dat->amb_flag, bc_ix);
    } else {
        bflg_unset(mdl_bc_dat->amb_flag, bc_ix);
    }
    
    // add to test bcs if not ambient and not empty
    if ((!ambient) && (!absent)) {
        if (add2str_map(mdl_bc_dat->test_bcs, bc_key, found) < 0)
            return -1;
    }
    
    return(bc_ix);
}

int mdl_bc_dat_bam_data(mdl_bc_dat_t *mdl_bc_dat, bam_data_t *bam_data, obj_pars *objs){
    if (mdl_bc_dat == NULL || bam_data == NULL || objs == NULL)
        return err_msg(-1, 0, "mdl_bc_dat_bam_data: argument is null");

    khint_t k_bc;

    uint32_t n_genes = 0;
    if (objs->anno)
        n_genes = objs->anno->gene_ix->n;
    mdl_bc_dat->G = n_genes;
    uint32_t n_vars = 0;
    if (objs->gv)
        n_vars = mv_size(&objs->gv->vix2var);
    mdl_bc_dat->V = n_vars;
    if (mdl_bc_dat->V < 1)
        return err_msg(-1, 0, "mdl_bc_dat_bam_data: no variants found");

    // flag to include only reads with >= variant bases overlapping
    uint8_t flg_rd_w_vars = 0;

    uint32_t c_thresh = (objs->out_min);
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
        assert(bam_bc != NULL);

        bc_stats_t *bc_stat = bam_bc->bc_stats;
        assert(bc_stat != NULL);

        int absent = (bc_stat->atac_counts == 0) && (bc_stat->rna_counts == 0);
        int low_count = (bc_stat->atac_counts < c_thresh) && (bc_stat->rna_counts < c_thresh);
        int dup_ok = low_count; // only low count empty barcode can be duplicate

        // if barcode is low count, set to empty
        if (low_count)
            bc_key = bc_empty;

        found = 0;
        bc_ix = mdl_bc_dat_add_bc(mdl_bc_dat, bc_key, absent, low_count, dup_ok, &found);
        if (bc_ix < 0)
            return(-1);

        // add reads to mdl_mlcl
        mdl_mlcl_bc_t *mdl_bc = &mv_i(&mdl_bc_dat->bc_mlcl, bc_ix);

        // loop through RNA
        rna_mlc_bag_itr itr, *ritrp = &itr;
        rna_mlc_bag_itr_first(ritrp, &bam_bc->rna_mlcs);
        for (; rna_mlc_bag_itr_alive(ritrp); rna_mlc_bag_itr_next(ritrp)) {
            rna_mol_t *mol = rna_mlc_bag_itr_val(ritrp);

            // skip RNA molecules with no data, and only consider unique feature
            size_t mol_n_genes = ml_size(&mol->gl), mol_n_vars = ml_size(&mol->vl);
            int skip_gene = !( mol_n_genes == 1 || mol_n_vars > 0 );
            if (skip_gene)
                continue;

            // skip RNA read if no variants and flag set
            if (flg_rd_w_vars > 0 && mol_n_vars == 0)
                continue;

            // initialize mlcl to 0
            mdl_mlcl_t mlcl;
            mdl_mlcl_init(&mlcl);
            mlcl.counts = 1;

            if (mdl_mlcl_add_rna(&mlcl, mol, n_genes, n_vars) < 0)
                return -1;

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
        atac_frag_bag_itr aitr, *aitrp = &aitr;
        atac_frag_bag_itr_first(aitrp, &bam_bc->atac_frgs);
        for (; atac_frag_bag_itr_alive(aitrp); atac_frag_bag_itr_next(aitrp)) {
            atac_frag_t *frag = atac_frag_bag_itr_val(aitrp);
            size_t frag_n_vars = ml_size(&frag->vl);

            // skip ATAC read if no variants and flag set
            if (flg_rd_w_vars > 0 && frag_n_vars == 0)
                continue;

            // initialize mlcl to all 0
            mdl_mlcl_t mlcl;
            mdl_mlcl_init(&mlcl);
            mlcl.counts = 1;

            if (mdl_mlcl_add_atac(&mlcl, frag, n_vars) < 0)
                return -1;

            mdl_mlcl_t *p;
            p = kb_getp(kb_mdl_mlcl, mdl_bc->atac, &mlcl);
            if (!p) {
                kb_putp(kb_mdl_mlcl, mdl_bc->atac, &mlcl);
            } else {
                p->counts += 1;
                mdl_mlcl_free(&mlcl);
            }
        }
        // note the read data is freed from the barcode here
        // TODO: implement an argument flag to control this
        bc_data_free_reads(bam_bc);
    }

    assert((uint32_t)mdl_bc_dat->all_bcs->n == mv_size(&mdl_bc_dat->bc_mlcl));
    log_msg("initialized %zu barcodes\n", mv_size(&mdl_bc_dat->bc_mlcl));

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

    mdl->n_ix = -1;
    mdl->n_hs = -1;
    mdl->k = -1;
    mdl->lp_hsk = NULL;
    mdl->sub_lp_hs = NULL;
    mdl->p_x = NULL;

    mdl->eps = 1e-5;
    mdl->max_iter = 100;
    mdl->alpha_max = 1;

    mdl->has_rna = 0;
    mdl->has_atac = 0;

    mdl->threads = 1;
    mdl->alpha_vars = 1;

    mdl->mdl_bc_dat = mdl_bc_dat_alloc();
    if (mdl->mdl_bc_dat == NULL) {
        hs_ix_dstry(mdl->hs_ix);
        free(mdl);
        return NULL;
    }

    mdl->sub_best_hs = NULL;
    mdl->best_hsk_ix = NULL;
    mdl->best_sng_ix = NULL;
    mdl->sec_sng_ix = NULL;
    mdl->best_dbl_ix = NULL;
    mdl->pp = NULL;
    mdl->lp_h = NULL;

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

    free(m->lp_hsk);
    free(m->sub_lp_hs);
    free(m->p_x);
    free(m->sub_best_hs);
    free(m->best_hsk_ix);
    free(m->best_sng_ix);
    free(m->sec_sng_ix);
    free(m->best_dbl_ix);
    free(m->pp);
    free(m->lp_h);

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
    if (mdl->k < 0)
        return err_msg(-1, 0, "mdl_set_hs_ix: run mdl_set_k before");

    uint16_t n_sam = mdl->samples->n;
    uint16_t n_clust = mdl->k;
    hs_ix_set(mdl->hs_ix, n_sam, n_clust);
    return 0;
}

int mdl_set_k(mdl_t *mdl, int k) {
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_set_k: argument is null");
    if (k < 1)
        return err_msg(-1, 0, "mdl_set_k: k=%i must be greater than 0", k);
    mdl->k = k;
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

    if (mdl->hs_ix->n_ixs < 1)
        return err_msg(-1, 0, "mdl_alloc_probs: indices must be set with "
                "'mdl_set_hs_ix'");
    if (mdl->k < 1)
        return err_msg(-1, 0, "mdl_alloc_probs: cell type k=%i must be >= 1", mdl->k);

    uint32_t D = mdl->mdl_bc_dat->all_bcs->n;
    uint32_t n_ix = mdl->hs_ix->n_ixs;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    mdl->n_ix = n_ix;
    mdl->n_hs = n_hs;
    mdl->lp_hsk = calloc(n_ix * D, sizeof(f_t));
    mdl->sub_lp_hs = calloc(D * n_hs, sizeof(f_t));
    mdl->p_x = calloc(D , sizeof(f_t));
    if (mdl->lp_hsk == NULL || mdl->p_x == NULL)
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

int pr_kd(mdl_pars_t *mp, par_ix_t *par_ix, f_t *prob) {
    int k = par_ix->k;
    if (k >= mp->K)
        return err_msg(-1, 0, "pr_kd:k=%u must be < K=%u", k, mp->K);
    if (k < 0) {
        *prob = 1;
        return 0;
    }
    if (prob_invalid(mp->kappa[k]) < 0)
        return -1;
    *prob = mp->kappa[k];
    return 0;
}

int pr_tdm(mdl_pars_t *mp, int rna, int bc_ix, par_ix_t *par_ix,
        int t_ix, f_t *prob) {

    int hs_ix = par_ix->hs_ix;
    int hsk_ix = par_ix->ix;
    int hd = par_ix->hd;
    int s_ix = par_ix->t_ix[t_ix];
    uint16_t M = mp->M;
    uint8_t c_ix = s_ix == M ? 0 : 1; // 0 is ambient, 1 is cell
    // TODO: how many droplets are in alpha?
    f_t al;
    if (rna) {
        al = mp->alpha_rna[CMI(bc_ix, hs_ix, mp->D)];
    } else {
        al = mp->alpha_atac[CMI(bc_ix, hs_ix, mp->D)];
    }
    f_t pr = 0;
    switch (hd) {
        case 0:
            if (c_ix == 0) pr = 1.0;
            else
                return err_msg(-1, 0, "pr_tdm: hd = 0 but s_ix != M", hd);
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
            return err_msg(-1, 0, "pr_tdm: hd=%i is invalid, there is a bug", hd);
    }
    *prob = pr;
    return(0);
}

f_t pr_rho_gene(f_t *rho, seq_gene_t seq_gene, uint32_t col, uint32_t G, uint32_t rho_nrow){
    uint32_t g_ix = (uint32_t)seq_gene.gene_id;
    uint32_t spl = (uint32_t)seq_gene.splice;
    uint32_t f_ix = g_ix + (spl * G);
    f_t p_gdm = rho[CMI(f_ix, col, rho_nrow)];
    return(p_gdm);
}

f_t pr_gamma_var(f_t *gamma, uint32_t v_ix, uint8_t allele, 
        uint32_t s_ix, uint32_t gamma_nrow, f_t tau){
    if (allele > 1) return(1.0); // return 1 if allele is missing
    uint32_t eix = CMI(s_ix, v_ix, gamma_nrow);
    f_t ap = gamma[eix];
    if (ap < 0) return 1.0; // if allele is missing
    if (allele == 0) ap = 1.0 - ap;
    f_t p_be0 = (1.0 - tau) * ap;
    f_t p_be1 = tau * 0.25;
    f_t p_bdm = p_be0 + p_be1;
    return(p_bdm);
}

int p_var(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int s_ix, f_t *prob){
    uint32_t V = mp->V;
    uint16_t M = mp->M;
    uint32_t gamma_nrow = M + 1;
    assert(s_ix < gamma_nrow);

    // probability terms
    f_t pbm = 1;

    // Pr(B_dm | T_dm, \rho)
    size_t i, n_vars = mv_size(&mlcl->var_ixs);
    for (i = 0; i < n_vars; ++i){
        int32_t vi = mv_i(&mlcl->var_ixs, i);
        if (vi < 0)
            return err_msg(-1, 0, "p_var: variant index=%i is negative", vi);
        uint32_t v = (uint32_t)vi;
        div_t di = div((int)v, (int)V);
        assert(di.quot < 3); // TODO: remove after testing
        uint8_t allele = di.quot;
        uint32_t v_ix = di.rem;
        if (allele > 1) continue; // if missing

        f_t pp = pr_gamma_var(mp->gamma, v_ix, allele, s_ix,
                gamma_nrow, mp->tau);
        if (prob_invalid(pp))
            return err_msg(-1, 0, "p_var: invalid base probability value '%e'", pp);
        pbm *= pp;
    }

    *prob = pbm;

    return 0;
}

int p_rna(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int s_ix, uint16_t k, f_t *prob){
    uint16_t M = mp->M;
    uint32_t G = mp->G;
    uint32_t V = mp->V;
    uint32_t gamma_nrow = M + 1;
    uint32_t rho_nrow = 3 * G;
    uint32_t rho_ncol = 1 + mp->K;
    uint8_t c_ix = s_ix == M ? 0 : k + 1; // 0 is ambient, 1 is cell

    // probability terms
    f_t pgbm = 1;

    // probability of variant
    if (p_var(mp, mlcl, s_ix, &pgbm) < 0)
        return -1;

    // Pr(G_dm | T_dm, \rho)
    size_t i;

    assert(pgbm > 0 && pgbm <= 1);

    // gene probability
    size_t n_feat = mv_size(&mlcl->feat_ixs);
    if (n_feat > 1)
        n_feat = 0;
    for (i = 0; i < n_feat; ++i){
        uint32_t row_ix = mv_i(&mlcl->feat_ixs, i);
        assert(row_ix < (rho_nrow)); // TODO: remove after testing
        uint32_t rho_ix = CMI(row_ix, c_ix, rho_nrow);
        assert(rho_ix < (rho_nrow * rho_ncol)); // TODO: remove after testing
        f_t rp = mp->rho[rho_ix];
        assert(rp > 0 && rp <= 1);
        pgbm *= mp->rho[rho_ix];
    }

    if (prob_invalid(pgbm))
        return err_msg(-1, 0, "p_rna: prob=%f is invalid", pgbm);

    *prob = pgbm;
    return(0);
}

int p_atac(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int s_ix, f_t *prob){
    uint16_t M = mp->M;
    uint32_t V = mp->V;
    uint32_t gamma_nrow = M + 1;
    uint8_t c_ix = s_ix == M ? 0 : 1; // 0 is ambient, 1 is cell

    // prob. value
    f_t ppbm = 1;

    // probability of variant
    if (p_var(mp, mlcl, s_ix, &ppbm) < 0)
        return -1;

    // peak probability
    int32_t pk = mv_i(&mlcl->feat_ixs, 0); // 0 for no peak, 1 for in peak

    f_t p_pdm;
    if (pk == 0)
        p_pdm = 1.0 - mp->sigma[c_ix];
    else
        p_pdm = mp->sigma[c_ix];

    ppbm *= p_pdm;

    if (prob_invalid(ppbm))
        return err_msg(-1, 0, "p_atac: prob=%f is invalid", ppbm);

    *prob = ppbm;
    return 0;
}

int p_bd(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int rna, int bc_ix, par_ix_t *par_ix,
        f_t *p_b, f_t *psum) {
    *psum = 0;

    int t_im, ret;
    for (t_im = 0; t_im < par_ix->t_n; ++t_im){
        p_b[t_im] = 1;
        int s_ix = par_ix->t_ix[t_im];
        assert(s_ix >= 0);

        f_t p_t = -1;
        ret = pr_tdm(mp, rna, bc_ix, par_ix, t_im, &p_t);
        if (ret < 0)
            return -1;

        f_t p_v = -1;
        ret = p_var(mp, mlcl, s_ix, &p_v);
        if (ret < 0)
            return -1;

        p_b[t_im] = p_t * p_v;
        if (prob_invalid(p_b[t_im]))
            return err_msg(-1, 0, "p_bd: prob=%f is invalid", p_b[t_im]);
        *psum += p_b[t_im];
    }

    return 0;
}

int p_f_v(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int rna, int bc_ix, par_ix_t *par_ix,
        f_t *p_tgb, f_t *psum) {
    *psum = 0;

    int t_im, ret;
    for (t_im = 0; t_im < par_ix->t_n; ++t_im){
        p_tgb[t_im] = 1;
        int s_ix = par_ix->t_ix[t_im];
        int k_ix = par_ix->k;
        assert(s_ix >= 0);
        f_t p_t;
        ret = pr_tdm(mp, rna, bc_ix, par_ix, t_im, &p_t);
        if (ret < 0)
            return -1;
        f_t p_fv = -1;
        if (rna != 0) {
            ret = p_rna(mp, mlcl, s_ix, k_ix, &p_fv);
        } else {
            ret = p_atac(mp, mlcl, s_ix, &p_fv);
        }
        if (ret < 0)
            return -1;
        p_tgb[t_im] = p_t * p_fv;
        if (prob_invalid(p_tgb[t_im]))
            return err_msg(-1, 0, "p_f_v: Pr(t, f, v)=%f is invalid", p_tgb[t_im]);
        *psum += p_tgb[t_im];
    }
    if (prob_invalid(*psum))
        return err_msg(-1, 0, "p_f_v: Pr(f, v)=%f is invalid", psum);
    return 0;
}

/*******************************************************************************
 * Expectation
 ******************************************************************************/

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

    int lsret = 0;

    uint32_t D = mdl->mp->D;
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    uint32_t hs_ix;

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    // pre-calculate P(H_d, S_d). Same for all droplets
    f_t *lp_hs_v = calloc(n_hs, sizeof(f_t));
    for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
        if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
            return -1;
        f_t p_hd = -1, p_sd = -1;
        pr_hd(mdl->mp, &par_ix, &p_hd);
        pr_sd(mdl->mp, &par_ix, &p_sd);
        if (p_hd < 0 || p_sd < 0)
            return err_msg(-1, 0, "mdl_sub_e: failed to get pr_hd or pr_sd");
        lp_hs_v[hs_ix] = log(p_hd) + log(p_sd);
    }

    // store probs
    f_t *lp_htd = calloc(n_hs, sizeof(f_t));
    if (lp_htd == NULL)
        return err_msg(-1, 0, "mdl_sub_e: %s", strerror(errno));

    uint32_t fixed = 0, unfixed = 0;
    uint32_t i;
    for (i = 0; i < ix_len; ++i){
        // get barcode and bam data
        int bc_ix = ixs[i];

        // if no barcode
        int efl = bflg_get(bd->absent_bc, bc_ix);
        if (efl == 1) continue;

        // if fixed, set Pr(empty) = 1
        int afl = bflg_get(bd->amb_flag, bc_ix);
        uint32_t bc_n_ix = afl == 1 ? 1 : n_hs;

        if (afl) ++fixed;
        else ++unfixed;

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);
        kbtree_t(kb_mdl_mlcl) *mols = bc_dat.rna;
        kbtree_t(kb_mdl_mlcl) *frags = bc_dat.atac;

        kbitr_t itr;

        // initialize log Pr(H_d, S_d, K_d)
        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix) {
            if (hs_ix < bc_n_ix)
                lp_htd[hs_ix] = bc_dat.n_bc * lp_hs_v[hs_ix];
            else
                lp_htd[hs_ix] = -INFINITY;
        }

        // loop over RNA
        kb_itr_first(kb_mdl_mlcl, mols, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            uint32_t n_vars = mv_size(&mlcl->var_ixs);
            if (n_vars < 1)
                continue;
            for (hs_ix = 0; hs_ix < bc_n_ix; ++hs_ix){
                // get (h,s,t) from index, store in par_ix
                if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                    return -1;

                f_t psum = 0, p_tgpb[3] = {1,1,1};
                int pret = p_bd(mdl->mp, mlcl, RNA_IX, bc_ix, &par_ix, p_tgpb, &psum);
                if (pret < 0)
                    return -1;
                assert(psum > 0);
                lp_htd[hs_ix] += mlcl->counts * log(psum);
            }
        }
        // loop over ATAC
        kb_itr_first(kb_mdl_mlcl, frags, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, frags, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            uint32_t n_vars = mv_size(&mlcl->var_ixs);
            if (n_vars < 1)
                continue;
            for (hs_ix = 0; hs_ix < bc_n_ix; ++hs_ix){
                // get (h,s,t) from index, store in par_ix
                if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                    return -1;

                f_t psum = 0, p_tgpb[3] = {1,1,1};
                int pret = p_bd(mdl->mp, mlcl, ATAC_IX, bc_ix, &par_ix, p_tgpb, &psum);
                if (pret < 0)
                    return -1;
                lp_htd[hs_ix] += mlcl->counts * log(psum);
            }
        }
        // add log Pr(H_d, S_d, K_d) to mdl
        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix) {
            uint32_t mix = CMI(bc_ix, hs_ix, D);
            mdl->sub_lp_hs[mix] = lp_htd[hs_ix];
        }

        // Pr(X_d | \Theta)
        mdl->p_x[bc_ix] = logsumexpd(lp_htd, bc_n_ix, &lsret);
        if (lsret < 0)
            return err_msg(-1, 0, "mdl_sub_e: could not logsumexp");
        if (mdl->p_x[bc_ix] > 1)
            fprintf(stdout, "mdl_sub_e: p_x[%i]=%.6e\n", bc_ix, mdl->p_x[bc_ix]);
        if (num_invalid(mdl->p_x[bc_ix]))
            return err_msg(-1, 0, "mdl_sub_e: p_x[%i]=%.6e is invalid\n",
                    bc_ix, mdl->p_x[bc_ix]);
    }
    free(lp_htd);
    free(lp_hs_v);

    return 0;
}

int mdl_full_e(mdl_t *mdl, int *ixs, uint32_t ix_len){
    if (mdl == NULL || ixs == NULL)
        return err_msg(-1, 0, "mdl_full_e: arguments are null");

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "mdl_full_e: model parameters are missing");
    if (mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "mdl_full_e: barcode count data is missing");
    if (mdl->mdl_bc_dat->all_bcs->n < 1)
        return err_msg(-1, 0, "mdl_full_e: no barcodes found");

    int lsret = 0;

    uint32_t D = mdl->mp->D;
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_hsk = mdl->hs_ix->n_ixs;
    uint32_t hs_ix;

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    // pre-calculate P(H_d, S_d, K_d). Same for all droplets
    f_t *lp_hsk_v = calloc(n_hsk, sizeof(f_t));
    for (hs_ix = 0; hs_ix < n_hsk; ++hs_ix){
        if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
            return -1;
        f_t p_hd = -1, p_sd = -1;
        pr_hd(mdl->mp, &par_ix, &p_hd);
        pr_sd(mdl->mp, &par_ix, &p_sd);
        if (p_hd < 0 || p_sd < 0)
            return err_msg(-1, 0, "mdl_full_e: failed to get pr_hd or pr_sd");
        lp_hsk_v[hs_ix] = log(p_hd) + log(p_sd);
        if (hs_ix > 0) {
            f_t p_kd;
            if (pr_kd(mdl->mp, &par_ix, &p_kd) < 0)
                return -1;
            assert(p_kd > 0);
            lp_hsk_v[hs_ix] += log(p_kd);
        }
    }

    // store probs
    f_t *lp_htd = calloc(n_hsk, sizeof(f_t));
    if (lp_htd == NULL)
        return err_msg(-1, 0, "mdl_full_e: %s", strerror(errno));

    f_t *cp_hskd = calloc(n_hsk, sizeof(f_t));
    if (cp_hskd == NULL)
        return err_msg(-1, 0, "mdl_full_e: %s", strerror(errno));

    uint32_t i, n_hs;
    for (i = 0; i < ix_len; ++i){
        // get barcode and bam data
        int bc_ix = ixs[i];

        // if no barcode
        int efl = bflg_get(bd->absent_bc, bc_ix);
        if (efl == 1) continue;

        // if fixed, set Pr(empty) = 1
        int afl = bflg_get(bd->amb_flag, bc_ix);
        uint32_t bc_n_ix = afl == 1 ? 1 : n_hsk;

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);
        kbtree_t(kb_mdl_mlcl) *mols = bc_dat.rna;
        kbtree_t(kb_mdl_mlcl) *frags = bc_dat.atac;

        // best (H,S) from sub model
        int32_t sub_best_hs = mdl->sub_best_hs[bc_ix];

        kbitr_t itr;

        // initialize log Pr(H_d, S_d, K_d)
        for (hs_ix = 0; hs_ix < n_hsk; ++hs_ix) {
            if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                return -1;
            if (par_ix.hs_ix == sub_best_hs)
                lp_htd[hs_ix] = bc_dat.n_bc * lp_hsk_v[hs_ix];
            else
                lp_htd[hs_ix] = -INFINITY;
        }

        // loop over RNA
        kb_itr_first(kb_mdl_mlcl, mols, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            for (hs_ix = 0; hs_ix < bc_n_ix; ++hs_ix){
                // get (h,s,t) from index, store in par_ix
                if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                    return -1;
                if (par_ix.hs_ix != sub_best_hs)
                    continue;

                f_t psum = 0, p_tgpb[3] = {1,1,1};
                int pret = p_f_v(mdl->mp, mlcl, RNA_IX, bc_ix, &par_ix, p_tgpb, &psum);
                if (pret < 0)
                    return -1;
                if (psum <= 0)
                    fprintf(stdout, "fl=%i (h=%i,s1=%i,s2=%i,k=%i) psum=%f\n",
                            afl, par_ix.hd, par_ix.s1, par_ix.s2, par_ix.k, psum);
                assert(psum > 0 && psum <= 1);
                lp_htd[hs_ix] += mlcl->counts * log(psum);
            }
        }
        // loop over ATAC
        kb_itr_first(kb_mdl_mlcl, frags, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, frags, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            for (hs_ix = 0; hs_ix < bc_n_ix; ++hs_ix){
                // get (h,s,t) from index, store in par_ix
                if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                    return -1;
                if (par_ix.hs_ix != sub_best_hs)
                    continue;

                f_t psum = 0, p_tgpb[3] = {1,1,1};
                int pret = p_f_v(mdl->mp, mlcl, ATAC_IX, bc_ix, &par_ix, p_tgpb, &psum);
                if (pret < 0)
                    return -1;
                if (psum <= 0)
                    fprintf(stdout, "fl=%i (h=%i,s1=%i,s2=%i,k=%i) psum=%f\n",
                            afl, par_ix.hd, par_ix.s1, par_ix.s2, par_ix.k, psum);
                lp_htd[hs_ix] += mlcl->counts * log(psum);
            }
        }
        // add log Pr(H_d, S_d, K_d) to mdl
        for (hs_ix = 0; hs_ix < n_hsk; ++hs_ix) {
            uint32_t mix = CMI(bc_ix, hs_ix, D);
            mdl->lp_hsk[mix] = lp_htd[hs_ix];
        }

        // Pr(X_d | \Theta)
        mdl->p_x[bc_ix] = logsumexpd(lp_htd, bc_n_ix, &lsret);
        if (lsret < 0)
            return err_msg(-1, 0, "mdl_full_e: could not logsumexp");
        if (mdl->p_x[bc_ix] > 1)
            fprintf(stdout, "mdl_full_e: p_x[%i]=%.6e\n", bc_ix, mdl->p_x[bc_ix]);
        if (num_invalid(mdl->p_x[bc_ix]))
            return err_msg(-1, 0, "mdl_full_e: p_x[%i]=%.6e is invalid\n",
                    bc_ix, mdl->p_x[bc_ix]);
    }
    free(lp_htd);
    free(lp_hsk_v);
    return 0;
}

/*******************************************************************************
 * Maximization
 ******************************************************************************/

int mdl_sub_m(mdl_t *mdl, int *ixs, uint32_t ix_len) {
    if (mdl == NULL || ixs == NULL)
        return err_msg(-1, 0, "mdl_sub_m: arguments are null");

    uint32_t i;
    uint16_t M = mdl->mp->M;
    uint32_t G = mdl->mp->G;
    uint32_t D = mdl->mp->D;
    uint16_t K = mdl->mp->K;
    uint32_t rho_ncol = mdl->mp->K + 1;
    uint32_t rho_nrow = 3 * G;
    uint32_t rho_ne = rho_ncol * rho_nrow;
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_hs = mdl->hs_ix->n_hs;

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    // lambda counter
    f_t lambda_sums[3] = {0,0,0};

    uint32_t hs_ix;
    for (i = 0; i < ix_len; ++i){
        // get barcode and bam data
        uint32_t bc_ix = ixs[i];

        // if no barcode
        int efl = bflg_get(bd->absent_bc, bc_ix);
        if (efl == 1) continue;

        // if fixed, set Pr(empty) = 1
        int afl = bflg_get(bd->amb_flag, bc_ix);

        uint32_t n_ix = n_hs;

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);

        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
            // get (h,s,t) from index
            if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                return(-1);

            if (afl && hs_ix > 0)
                continue;

            // for alpha
            f_t a_tot, new_par, pdiff;
            f_t psc = 1e-8;
            f_t ar[2] = {psc, psc}; // 1st index is ambient
            f_t aa[2] = {psc, psc};

            // get P(H,S|X)
            uint32_t mix = CMI(bc_ix, hs_ix, D);
            f_t cp_hsx = exp(mdl->sub_lp_hs[mix] - mdl->p_x[bc_ix]);
            if (prob_invalid(cp_hsx))
                return err_msg(-1, 0, "mdl_sub_m: "
                        "P(Z|X)=%.6e is an invalid probability",
                        cp_hsx);

            // skip over very small prob values
            if (cp_hsx < 1e-100)
                continue;

            // lambda
            lambda_sums[par_ix.hd] += bc_dat.n_bc * cp_hsx;

            // Get P(H,S,X)
            // Divide by \sum P(T,G,P,B,X)
            // Multiple by P(T,G,P,B,X)

            // bam data
            kbtree_t(kb_mdl_mlcl) *mols = bc_dat.rna;
            kbtree_t(kb_mdl_mlcl) *frags = bc_dat.atac;

            kbitr_t itr_rna, itr_atac;

            // loop over RNA
            kb_itr_first(kb_mdl_mlcl, mols, &itr_rna); 
            for (; kb_itr_valid(&itr_rna); kb_itr_next(kb_mdl_mlcl, mols, &itr_rna)){
                mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr_rna);
                uint32_t n_vars = mv_size(&mlcl->var_ixs);
                if (n_vars < 1)
                    continue;

                // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
                f_t psum = 0, p_tfb[3] = {1,1,1};
                int pret = p_bd(mdl->mp, mlcl, RNA_IX, bc_ix, &par_ix, p_tfb, &psum);
                if (pret < 0)
                    return -1;
                assert(psum > 0);

                // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
                int t_im;
                for (t_im = 0; t_im < par_ix.t_n; ++t_im){
                    int s_ix = par_ix.t_ix[t_im];
                    uint8_t c_ix = s_ix == M ? 0 : 1; // 0 is ambient, 1 is cell
                    f_t pt = p_tfb[t_im] / psum;
                    if (prob_invalid(pt)) {
                        printf("BC RNA with fix=%i has cp=%f psum=%f, "
                                "pt=%f, hs=%i\n", afl, cp_hsx, psum, pt, hs_ix);
                        return err_msg(-1, 0, "mdl_sub_m: pt=%f is invalid prob", pt);
                    }

                    // alpha
                    uint32_t coln = c_ix == 0 ? 0 : 1;
                    ar[coln] += mlcl->counts * pt;
                    // mix = CMI(hs_ix, coln, n_hs);
                    // alpha_sums[mix] += mlcl->counts * pt;
                }
            }

            // RNA alpha
            a_tot = ar[0] + ar[1];
            new_par = ar[0] / a_tot;
            pdiff = new_par - mdl->mp->alpha_rna[CMI(bc_ix, hs_ix, D)];
            mdl->mp->_par_diff += fabs(pdiff);
            mdl->mp->alpha_rna[CMI(bc_ix, hs_ix, D)] = new_par;

            // loop over ATAC
            kb_itr_first(kb_mdl_mlcl, frags, &itr_atac); 
            for (; kb_itr_valid(&itr_atac); kb_itr_next(kb_mdl_mlcl, frags, &itr_atac)){
                mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr_atac);
                uint32_t n_vars = mv_size(&mlcl->var_ixs);
                if (n_vars < 1)
                    continue;

                // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
                f_t psum = 0, p_tfb[3] = {1,1,1};
                int pret = p_bd(mdl->mp, mlcl, ATAC_IX, bc_ix, &par_ix, p_tfb, &psum);
                if (pret < 0)
                    return -1;
                assert(psum > 0);

                // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
                int t_im;
                for (t_im = 0; t_im < par_ix.t_n; ++t_im){
                    int s_ix = par_ix.t_ix[t_im];
                    uint32_t c_ix = s_ix == M ? 0 : 1; // 0 for ambient, 1 for cell
                    f_t pt = p_tfb[t_im] / psum; // removes cp_hsx constant
                    if (prob_invalid(pt)) {
                        printf("BC ATAC has cp=%f psum=%f, "
                                "pt=%f, hs=%i\n", cp_hsx, psum, pt, hs_ix);
                        return err_msg(-1, 0, "mdl_sub_m: pt=%f is invalid prob", pt);
                    }

                    // alpha
                    uint32_t coln = c_ix == 0 ? 0 : 1;
                    aa[coln] += mlcl->counts * pt;
                    // uint32_t coln = c_ix == 0 ? 2 : 3;
                    // mix = CMI(hs_ix, coln, n_hs);
                    // alpha_sums[mix] += mlcl->counts * pt;
                }
            }

            // ATAC alpha
            a_tot = aa[0] + aa[1];
            new_par = aa[0] / a_tot;
            pdiff = new_par - mdl->mp->alpha_atac[CMI(bc_ix, hs_ix, D)];
            mdl->mp->_par_diff += fabs(pdiff);
            mdl->mp->alpha_atac[CMI(bc_ix, hs_ix, D)] = new_par;
        }
    }

    // add to tmp sum variables
    pthread_mutex_lock(&mdl->mp->sum_lock);
    // add lambda
    for (i = 0; i < 3; ++i){
        mdl->mp->_lambda_sum[i] += lambda_sums[i];
    }

    pthread_mutex_unlock(&mdl->mp->sum_lock);

    return(0);
}

int mdl_full_m(mdl_t *mdl, int *ixs, uint32_t ix_len) {
    if (mdl == NULL || ixs == NULL)
        return err_msg(-1, 0, "mdl_full_m: arguments are null");


    uint32_t i;
    uint16_t M = mdl->mp->M;
    uint32_t G = mdl->mp->G;
    uint32_t D = mdl->mp->D;
    uint16_t K = mdl->mp->K;
    uint32_t rho_ncol = mdl->mp->K + 1;
    uint32_t rho_nrow = 3 * G;
    uint32_t rho_ne = rho_ncol * rho_nrow;
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    uint32_t n_hsk = mdl->hs_ix->n_ixs;

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    // rho counter
    f_t *rho_sum = malloc(rho_ne * sizeof(f_t));
    for (i = 0; i < rho_ne; ++i) rho_sum[i] = 0;

    // sigma counter
    uint32_t sig_ncol = 2;
    uint32_t sig_nrow = 2;
    uint32_t sig_ne = sig_ncol * sig_nrow;
    // sig_sum col 0: ambient, col 1: cell, row 0: outside peak, row 1: inside peak
    f_t *sig_sum = malloc(sig_ne * sizeof(f_t)); // outside/inside peak by ambient/nuclear
    for (i = 0; i < sig_ne; ++i)
        sig_sum[i] = 0;

    // kappa counter
    f_t *kappa_sum = malloc(K * sizeof(f_t));
    for (i = 0; i < K; ++i)
        kappa_sum[i] = 0;

    uint32_t hs_ix, hsk_ix;
    for (i = 0; i < ix_len; ++i){
        // get barcode and bam data
        uint32_t bc_ix = ixs[i];

        int fl;

        // if no barcode
        fl = bflg_get(bd->absent_bc, bc_ix);
        if (fl == 1) continue;

        // if fixed, set Pr(empty) = 1
        fl = bflg_get(bd->amb_flag, bc_ix);
        // n_hsk = fl == 1 ? 1 : n_ixs;
        uint32_t n_ix = n_hsk;

        // best (H,S) from sub model
        int32_t sub_best_hs = mdl->sub_best_hs[bc_ix];

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);

        for (hsk_ix = 0; hsk_ix < n_hsk; ++hsk_ix){
            // get (h,s,t) from index
            if (hs_ix_get_pars(mdl->hs_ix, hsk_ix, &par_ix) < 0)
                return(-1);
            hs_ix = par_ix.hs_ix;
            if (par_ix.hs_ix != sub_best_hs)
                continue;

            if (fl && hsk_ix > 0)
                continue;

            int n_k = par_ix.k;

            // get P(H,S|X)
            uint32_t mix = CMI(bc_ix, hsk_ix, D);
            f_t cp_hsx = exp(mdl->lp_hsk[mix] - mdl->p_x[bc_ix]);
            if (prob_invalid(cp_hsx))
                return err_msg(-1, 0, "mdl_full_m: "
                        "P(Z|X)=%.6e is an invalid probability",
                        cp_hsx);
            /*
               if (bc_ix % 1000 == 0)
               fprintf(stdout, "ix=%u, cp_hsx[%u,%u]=%f\n", bc_ix, hsk_ix, i_k, cp_hsx);
               */

            // skip over very small prob values
            if (cp_hsx < 1e-100)
                continue;

            // kappa
            if (n_k >= 0)
                kappa_sum[n_k] += bc_dat.n_bc * cp_hsx;

            // Get P(H,S,X)
            // Divide by \sum P(T,G,P,B,X)
            // Multiple by P(T,G,P,B,X)

            // bam data
            kbtree_t(kb_mdl_mlcl) *mols = bc_dat.rna;
            kbtree_t(kb_mdl_mlcl) *frags = bc_dat.atac;

            kbitr_t itr_rna, itr_atac;

            // loop over RNA
            kb_itr_first(kb_mdl_mlcl, mols, &itr_rna); 
            for (; kb_itr_valid(&itr_rna); kb_itr_next(kb_mdl_mlcl, mols, &itr_rna)){
                mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr_rna);

                // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
                f_t psum = 0, p_tfb[3] = {1,1,1};
                int pret = p_f_v(mdl->mp, mlcl, RNA_IX, bc_ix, &par_ix, p_tfb, &psum);
                if (pret < 0)
                    return -1;
                assert(psum > 0);

                // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
                int t_im;
                for (t_im = 0; t_im < par_ix.t_n; ++t_im){
                    int s_ix = par_ix.t_ix[t_im];
                    uint8_t c_ix = s_ix == M ? 0 : n_k + 1; // 0 is ambient, otherwise cell
                    f_t pt = p_tfb[t_im] / psum;
                    if (prob_invalid(pt)) {
                        printf("BC RNA with fix=%i has cp=%f psum=%f, "
                                "pt=%f, hs=%i\n", fl, cp_hsx, psum, pt, hsk_ix);
                        return err_msg(-1, 0, "mdl_full_m: pt=%f is invalid prob", pt);
                    }

                    // rho
                    assert(mv_size(&mlcl->feat_ixs) < 2);
                    uint32_t l;
                    for (l = 0; l < mv_size(&mlcl->feat_ixs); ++l){
                        uint32_t row_ix = (uint32_t)(mv_i(&mlcl->feat_ixs, l));
                        assert(row_ix < rho_nrow);
                        uint32_t rho_ix = CMI(row_ix, c_ix, rho_nrow);
                        rho_sum[rho_ix] += mlcl->counts * cp_hsx * pt;
                    }
                }
            }

            // loop over ATAC
            kb_itr_first(kb_mdl_mlcl, frags, &itr_atac); 
            for (; kb_itr_valid(&itr_atac); kb_itr_next(kb_mdl_mlcl, frags, &itr_atac)){
                mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr_atac);

                // assert(mlcl->counts == 1);

                // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
                f_t psum = 0, p_tfb[3] = {1,1,1};
                int pret = p_f_v(mdl->mp, mlcl, ATAC_IX, bc_ix, &par_ix, p_tfb, &psum);
                if (pret < 0)
                    return err_msg(-1, 0, "failed for bc=%u hs=%u k=%u", bc_ix,
                            hsk_ix, n_k);
                assert(psum > 0 && psum <= 1);

                assert(mv_size(&mlcl->feat_ixs) == 1);
                int32_t pk32 = mv_i(&mlcl->feat_ixs, 0);
                if (pk32 != 0 && pk32 != 1)
                    fprintf(stderr, "pk32=%i\n", pk32);
                assert(pk32 == 0 || pk32 == 1);
                uint32_t pk = pk32 == 0 ? 0 : 1; // 0: outside peak, 1: inside peak

                // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
                int t_im;
                for (t_im = 0; t_im < par_ix.t_n; ++t_im){
                    int s_ix = par_ix.t_ix[t_im];
                    uint32_t a_type = s_ix == M ? 0 : 1; // 0 for ambient, 1 for cell
                                                         // f_t pt = phs1 * p_tfb[t_im];
                    f_t pt = p_tfb[t_im] / psum; // removes cp_hsx constant
                    if (prob_invalid(pt)) {
                        printf("BC ATAC has cp=%f psum=%f, "
                                "pt=%f, hs=%i\n", cp_hsx, psum, pt, hsk_ix);
                        return err_msg(-1, 0, "mdl_full_m: pt=%f is invalid prob", pt);
                    }

                    // sigma
                    uint32_t sig_ix = CMI(pk, a_type, sig_nrow);
                    sig_sum[sig_ix] += mlcl->counts * cp_hsx * pt;
                }
            }

        }
    }

    // add to tmp sum variables
    pthread_mutex_lock(&mdl->mp->sum_lock);
    // add kappa
    for (i = 0; i < K; ++i){
        mdl->mp->_kappa_sum[i] += kappa_sum[i];
    }

    // add rho
    if (mdl->has_rna){
        for (i = 0; i < rho_ne; ++i){
            mdl->mp->_rho_sum[i] += rho_sum[i];
        }
    }

    // add sigma
    if (mdl->has_atac){
        for (i = 0; i < sig_ne; ++i){
            mdl->mp->_sigma_sum[i] += sig_sum[i];
        }
    }
    pthread_mutex_unlock(&mdl->mp->sum_lock);

    free(kappa_sum);
    free(rho_sum);
    free(sig_sum);

    return(0);
}

int mdl_m_lambda(mdl_t *mdl) {
    f_t new_par = 0;
    // mdl->mp->_par_diff = 0;
    uint32_t i;

    // maximize lambda
    f_t lambda_tot = 0;
    for (i = 0; i < 3; ++i)
        lambda_tot += mdl->mp->_lambda_sum[i];
    for (i = 0; i < 3; ++i){
        new_par = mdl->mp->_lambda_sum[i] / lambda_tot;
        mdl->mp->_par_diff += fabs(new_par - mdl->mp->lambda[i]);
        mdl->mp->lambda[i] = new_par;
    }

    return 0;
}

int mdl_m_kappa(mdl_t *mdl) {
    // maximize kappa
    f_t kappa_tot = 0;
    uint32_t i;
    for (i = 0; i < mdl->mp->K; ++i)
        kappa_tot += mdl->mp->_kappa_sum[i];
    for (i = 0; i < mdl->mp->K; ++i)
        mdl->mp->kappa[i] = mdl->mp->_kappa_sum[i] / kappa_tot;

    return 0;
}

int mdl_m_rho(mdl_t *mdl) {
    f_t new_par = 0;

    uint32_t G = mdl->mp->G;
    uint32_t rho_nrow = 3 * G;
    uint32_t rho_ncol = 1 + mdl->mp->K;

    uint32_t i, j;
    if (mdl->has_rna){
        f_t *rho_tot = calloc(rho_ncol, sizeof(f_t));
        for (i = 0; i < rho_ncol; ++i)
            rho_tot[i] = 0;
        for (i = 0; i < rho_nrow; ++i){
            for (j = 0; j < rho_ncol; ++j){
                uint32_t rix = CMI(i, j, rho_nrow);
                rho_tot[j] += mdl->mp->_rho_sum[rix];
            }
        }

        // check for errors
        for (j = 0; j < rho_ncol; ++j){
            if (num_invalid(rho_tot[j]))
                return err_msg(-1, 0, "mdl_m_rho: rho prob. sum is invalid: [%f]", 
                        rho_tot[j]);
        }

        for (i = 0; i < rho_nrow; ++i){
            for (j = 0; j < rho_ncol; ++j){
                uint32_t rix = CMI(i, j, rho_nrow);
                new_par = mdl->mp->_rho_sum[rix] / rho_tot[j];
                mdl->mp->_par_diff += fabs(new_par - mdl->mp->rho[rix]);
                mdl->mp->rho[rix] = new_par;
            }
        }
        free(rho_tot);
    }

    return 0;
}

int mdl_m_sigma(mdl_t *mdl) {
    f_t new_par = 0;
    uint32_t i;
    if (mdl->has_atac){
        uint32_t sig_nrow = 2;
        f_t sig_tot[2];
        for (i = 0; i < 2; ++i){
            sig_tot[i] = mdl->mp->_sigma_sum[CMI(0,i,sig_nrow)] + 
                mdl->mp->_sigma_sum[CMI(1,i,sig_nrow)];
        }
        if (num_invalid(sig_tot[0]) || num_invalid(sig_tot[1]))
            return err_msg(-1, 0, "mdl_m_rho: sigma prob. sum is invalid: [%f,%f]", 
                    sig_tot[0], sig_tot[1]);

        new_par = mdl->mp->_sigma_sum[ CMI(1, 0, sig_nrow) ] / sig_tot[0];
        mdl->mp->_par_diff += fabs(new_par - mdl->mp->sigma[0]);
        mdl->mp->sigma[0] = new_par;

        new_par = mdl->mp->_sigma_sum[ CMI(1, 1, sig_nrow) ] / sig_tot[1];
        mdl->mp->_par_diff += fabs(new_par - mdl->mp->sigma[1]);
        mdl->mp->sigma[1] = new_par;
    }
    return 0;
}

int mdl_sub_est(mdl_t *mdl) {
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_sub_est: mdl is null");

    if (mdl_m_lambda(mdl) < 0)
        return -1;

    // maximize pi (keep fixed uniform for now)
    if (mdl_pars_pi_fix(mdl->mp) < 0)
        return -1;

    // update gamma
    if (mdl_pars_set_gamma_amb(mdl->mp) < 0) return(-1);

    // alpha estimated in m_sum

    return 0;
}

int mdl_full_est(mdl_t *mdl){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_full_est: mdl is null");

    // maximize kappa
    if (mdl_m_kappa(mdl) < 0)
        return -1;

    // rho
    if (mdl_m_rho(mdl) < 0)
        return -1;

    // sigma
    if (mdl_m_sigma(mdl) < 0)
        return -1;

    return(0);
}

int mdl_delta_pars(mdl_t *mdl, f_t *delta){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_delta_pars: arguments are NULL");

    f_t n_pars = 0.0;
    n_pars += 3.0; // lambda
    n_pars += mdl->mp->M; // pi
    n_pars += mdl->mp->K; // kappa
    if (mdl->has_rna) {
        n_pars += (f_t)mdl->mp->D * (f_t)mdl->mp->n_ix; // alpha (for unfixed droplets)
        n_pars += 3.0 * ((f_t)mdl->mp->K + 1) * (f_t)mdl->mp->G; // rho
    }
    if (mdl->has_atac) {
        n_pars += (f_t)mdl->mp->D * (f_t)mdl->mp->n_ix; // alpha (for unfixed droplets)
        n_pars += 2.0; // sigma
    }

    *delta = mdl->mp->_par_diff / n_pars;
    mdl->mp->_par_diff = 0;

    return(0);
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
    if (mdl_sub_m(t->mdl, t->ixs, t->ix_len) < 0) return((void *)-1);
    return((void *)0);
}

void *mdl_full_thrd_fx(void *arg){
    thrd_args *t = (thrd_args *)arg;
    if (mdl_full_e(t->mdl, t->ixs, t->ix_len) < 0) return((void *)-1);
    if (mdl_full_m(t->mdl, t->ixs, t->ix_len) < 0) return((void *)-1);
    return((void *)0);
}

int mdl_thrd_call(mdl_t *mdl, int *ixs, uint32_t ix_len, int sub){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_thrd_call: arguments are NULL");

    uint32_t i;
    if (mdl->threads < 1)
        return err_msg(-1, 0, "mdl_thrd_call: threads=%u is less than 1", mdl->threads);

    if (ix_len == 0)
        return err_msg(-1, 0, "mdl_thrd_call: no barcodes present");

    assert(ix_len > 0);
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
        if (sub)
            err = pthread_create(ids + i, NULL, mdl_sub_thrd_fx, &targs[i]);
        else
            err = pthread_create(ids + i, NULL, mdl_full_thrd_fx, &targs[i]);
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

int mdl_em(mdl_t *mdl, obj_pars *objs, int sub){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_em_step: arguments are NULL");

    f_t q = 0, q_prev = -INFINITY;
    f_t par_delta, q_delta;
    f_t prior = 1e-8;
    uint32_t i;
    uint32_t D = mdl->mp->D;
    int *ixs = malloc(D * sizeof(int));
    for (i = 0; i < D; ++i)
        ixs[i] = i;
    permute(ixs, D);

    int fret; // store function return
    if (objs->verbose) {
        log_msg("%u test droplets", mdl->mp->T);
        log_msg("running EM with %u droplets", D);
    }
    for (i = 0; i < mdl->max_iter; ++i){
        if (mdl_pars_check(mdl->mp) < 0)
            return -1;

        // initialize sums to prior pseudocount
        fret = mdl_pars_reset_sums(mdl->mp, prior);
        if (fret < 0)
            return -1;

        log_msg("E step");
        if (sub)
            fret = mdl_thrd_call(mdl, ixs, D, 1);
        else
            fret = mdl_thrd_call(mdl, ixs, D, 0);
        if (fret < 0)
            return -1;

        log_msg("M step");
        if (sub)
            fret = mdl_sub_est(mdl);
        else
            fret = mdl_full_est(mdl);
        if (fret < 0)
            return -1;

        log_msg("llk");
        fret = mdl_data_llk(mdl, &q);
        if (fret < 0)
            return -1;

        if (mdl_delta_q(q_prev, q, &q_delta) < 0)
            return -1;
        if (objs->verbose)
            log_msg("Q=%.6e; delta=%.6e", q, q_delta);
        if (q_delta < 0)
            err_msg(0, 1, "llk decreased: q1=%.6e q2=%.6e", q_prev, q);
        q_prev = q;

        if (mdl_delta_pars(mdl, &par_delta) < 0)
            return(-1);
        log_msg("iteration %u: delta=%f", i+1, par_delta);
        if (q_delta < mdl->eps) break;
    }
    if (objs->verbose) log_msg("finished EM");
    if (mdl_pars_check(mdl->mp) < 0)
        return -1;

    free(ixs);

    return(0);
}

int mdl_data_llk(mdl_t *mdl, f_t *llk){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_llk: mdl is null");
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    if (bd == NULL)
        return err_msg(-1, 0, "mdl_llk: barcode data is not initialized");
    if (bd->all_bcs == NULL)
        return err_msg(-1, 0, "mdl_llk: no barcodes added");
    if (bd->all_bcs->n < 1)
        return err_msg(-1, 0, "mdl_llk: no barcodes present");

    uint32_t i, n_bcs = (uint32_t)bd->all_bcs->n;
    f_t llkt = 0;
    int fl1;
    for (i = 0; i < n_bcs; ++i){
        fl1 = bflg_get(bd->absent_bc, i);
        if (fl1)
            continue;
        if (mdl->p_x[i] > 1)
            fprintf(stdout, "p_x[%i]=%.6e\n", i, mdl->p_x[i]);
        llkt += mdl->p_x[i];
        assert(!num_invalid(llkt));
    }
    *llk = llkt;
    return(0);
}

f_t *mdl_get_lp_hs(mdl_t *mdl) {
    if (mdl == NULL)
        return NULL;

    uint32_t i, j;
    uint32_t D = mdl->mp->D;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    uint32_t n_hsk = mdl->mp->n_ix;
    uint32_t n_e = n_hs * D;
    f_t *lp_hs = malloc(n_e * sizeof(f_t));
    if (lp_hs == NULL)
        return NULL;
    for (i = 0; i < n_e; ++i)
        lp_hs[i] = -INFINITY;

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    uint32_t bc_i;
    for (bc_i = 0; bc_i < D; ++bc_i){
        // calculate P(H_d, S_d | X_d)
        for (j = 0; j < n_hsk; ++j) {
            if (hs_ix_get_pars(mdl->hs_ix, j, &par_ix) < 0)
                return NULL;
            uint32_t hs_ix = par_ix.hs_ix;
            assert(hs_ix <= n_hs);
            uint32_t mix = CMI(bc_i, j, D);
            lp_hs[CMI(bc_i, hs_ix, D)] =
                logsum2expd(lp_hs[CMI(bc_i, hs_ix, D)], mdl->lp_hsk[mix]);
        }
    }
    return lp_hs;
}

f_t *mdl_get_cp_hsk(mdl_t *mdl) {
    if (mdl == NULL)
        return NULL;

    if (mdl->mp == NULL) {
        err_msg(-1, 0, "mdl_get_cp_hsk: no parameters initialized");
        return NULL;
    }
    if (mdl->lp_hsk == NULL || mdl->p_x == NULL) {
        err_msg(-1, 0, "mdl_get_cp_hsk: no probabilities found");
        return NULL;
    }

    uint32_t D = mdl->mp->D;
    uint32_t T = mdl->mp->T;
    uint32_t n_hsk = mdl->mp->n_ix;

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    f_t *cp_hs = malloc(T * n_hsk * sizeof(f_t));
    if (cp_hs == NULL)
        return NULL;
    uint32_t bc_i, hsk_i;
    for (bc_i = 0; bc_i < T; ++bc_i){
        // get barcode indices
        char *bc = str_map_str(bd->test_bcs, bc_i);
        int all_bc_ix = str_map_ix(bd->all_bcs, bc);

        for (hsk_i = 0; hsk_i < n_hsk; ++hsk_i) {
            cp_hs[CMI(bc_i, hsk_i, T)] = exp(mdl->lp_hsk[CMI(all_bc_ix, hsk_i, D)] - mdl->p_x[all_bc_ix]);
        }
    }

    return cp_hs;
}

int mdl_get_prob_h(mdl_t *mdl) {
    if (mdl == NULL)
        return 0;

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "mdl_get_prob_h: no parameters initialized");
    if (mdl->lp_hsk == NULL || mdl->p_x == NULL)
        return err_msg(-1, 0, "mdl_get_prob_h: no probabilities found");
    if (mdl->samples == NULL)
        return err_msg(-1, 0, "mdl_get_prob_h: no samples found");
    if (mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "mdl_get_prob_h: no barcodes initialized");
    if (mdl->mdl_bc_dat->test_bcs->n < 1)
        return err_msg(-1, 0, "mdl_get_prob_h: no test barcodes present");

    int lsret = 0;
    uint32_t i;
    uint32_t D = mdl->mp->D;
    uint32_t T = mdl->mp->T;
    uint32_t n_hsk = mdl->mp->n_ix;
    mdl->pp = malloc(T * 3 * sizeof(f_t));
    mdl->lp_h = malloc(T * 3 * sizeof(f_t));
    if (mdl->pp == NULL || mdl->lp_h == NULL)
        return err_msg(-1, 0, "mdl_get_prob_h: %s", strerror(errno));

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    f_t lsum[3];

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t bc_i;
    for (bc_i = 0; bc_i < T; ++bc_i){
        // get barcode indices
        char *bc = str_map_str(bd->test_bcs, bc_i);
        int all_bc_ix = str_map_ix(bd->all_bcs, bc);

        for (i = 0; i < 3; ++i)
            lsum[i] = -INFINITY;

        // get posterior probabilities of H_d
        for (i = 0; i < n_hsk; ++i) {
            if (hs_ix_get_pars(mdl->hs_ix, i, &par_ix) < 0)
                return -1;
            uint32_t mix = CMI(all_bc_ix, i, D);
            lsum[par_ix.hd] =
                logsum2expd(lsum[par_ix.hd], mdl->lp_hsk[mix]);
        }
        f_t lsumall = logsumexpd(lsum, 3, &lsret);
        if (lsret < 0)
            return -1;
        for (i = 0; i < 3; ++i) {
            mdl->lp_h[CMI(bc_i, i, T)] = lsum[i];
            mdl->pp[CMI(bc_i, i, T)] = exp(lsum[i] - lsumall);
        }
    }

    return 0;
}

int mdl_sub_best_hs(mdl_t *mdl) {
    if (mdl == NULL)
        return 0;

    if (mdl->sub_lp_hs == NULL || mdl->mp == NULL || mdl->mp->D == 0)
        return err_msg(-1, 0, "mdl_sub_best_hs: model must be initialized");

    uint32_t n_hs = mdl->n_hs;
    uint32_t D = mdl->mp->D;
    mdl->sub_best_hs = malloc(D * sizeof(uint32_t));
    if (mdl->sub_best_hs == NULL)
        return err_msg(-1, 0, "mdl_sub_best_hs: %s", strerror(errno));

    uint32_t hs_ix = 0;
    uint32_t bc_i;
    for (bc_i = 0; bc_i < D; ++bc_i){
        // get best (H,S) index
        uint32_t b_hs_ix = 0;
        uint32_t t_hs_ix = 0;
        for (t_hs_ix = 1; t_hs_ix < n_hs; ++t_hs_ix) {
            uint32_t b_mix = CMI(bc_i, b_hs_ix, D);
            uint32_t t_mix = CMI(bc_i, t_hs_ix, D);
            if (mdl->sub_lp_hs[t_mix] > mdl->sub_lp_hs[b_mix])
                b_hs_ix = t_hs_ix;
        }
        mdl->sub_best_hs[bc_i] = b_hs_ix;
    }

    return 0;
}

int mdl_get_best(mdl_t *mdl){
    if (mdl == NULL)
        return 0;

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "mdl_get_best: no parameters initialized");
    if (mdl->lp_hsk == NULL || mdl->p_x == NULL)
        return err_msg(-1, 0, "mdl_get_best: no probabilities found");
    if (mdl->pp == NULL)
        return err_msg(-1, 0, "mdl_get_best: mdl_get_prob_h must be run");
    if (mdl->samples == NULL)
        return err_msg(-1, 0, "mdl_get_best: no samples found");
    if (mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "mdl_get_best: no barcodes initialized");
    if (mdl->mdl_bc_dat->test_bcs->n < 1)
        return err_msg(-1, 0, "mdl_get_best: no test barcodes present");

    uint32_t M = mdl->samples->n;
    uint32_t D = mdl->mp->D;
    uint32_t T = mdl->mp->T;
    uint32_t n_hsk = mdl->hs_ix->n_ixs;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    mdl->best_hsk_ix = malloc(T * sizeof(uint32_t));
    mdl->best_sng_ix = malloc(T * sizeof(uint32_t));
    mdl->sec_sng_ix = malloc(T * sizeof(uint32_t));
    mdl->best_dbl_ix = malloc(T * sizeof(uint32_t));

    f_t *lp_hs = mdl_get_lp_hs(mdl);
    if (lp_hs == NULL)
        return err_msg(-1, 0, "mdl_get_best: failed to get lp_hs");

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t hs_ix = 0;
    uint32_t bc_i;
    for (bc_i = 0; bc_i < T; ++bc_i){
        // get barcode indices
        char *bc = str_map_str(bd->test_bcs, bc_i);
        if (bc == NULL)
            return err_msg(-1, 0, "mdl_get_best: test bc %u not found", bc_i);
        int all_bc_ix = str_map_ix(bd->all_bcs, bc);
        if (all_bc_ix < 0)
            return err_msg(-1, 0, "mdl_get_best: test bc %s not found in all bcs", bc);

        // get best (H,S,K) index
        uint32_t best_hsk = 0;
        for (hs_ix = 0; hs_ix < n_hsk; ++hs_ix) {
            uint32_t mix1 = CMI(all_bc_ix, hs_ix, D);
            uint32_t mix2 = CMI(all_bc_ix, best_hsk, D);
            if (mdl->lp_hsk[mix1] >= mdl->lp_hsk[mix2])
                best_hsk = hs_ix;
        }

        // get singlet index
        uint32_t si[2] = {1, 2};
        f_t sllk[2] = {-INFINITY, -INFINITY};
        for (hs_ix = 1; hs_ix < M + 1; ++hs_ix){
            uint32_t mix = CMI(all_bc_ix, hs_ix, D);
            if (lp_hs[mix] >= sllk[0]) {
                sllk[0] = lp_hs[mix];
                si[0] = hs_ix;
                continue;
            }
            if (lp_hs[mix] > sllk[1]) {
                sllk[1] = lp_hs[mix];
                si[1] = hs_ix;
            }
        }

        // get doublet index
        uint32_t di[2] = {hs_ix, hs_ix + 1};
        f_t dllk[2] = {-INFINITY, -INFINITY};
        for (; hs_ix < n_hs; ++hs_ix){
            uint32_t mix = CMI(all_bc_ix, hs_ix, D);
            if (lp_hs[mix] >= dllk[0]) {
                dllk[0] = lp_hs[mix];
                di[0] = hs_ix;
                continue;
            }
            if (lp_hs[mix] > dllk[1]) {
                dllk[1] = lp_hs[mix];
                di[1] = hs_ix;
            }
        }
        mdl->best_sng_ix[bc_i] = si[0];
        mdl->sec_sng_ix[bc_i] = si[1];
        mdl->best_dbl_ix[bc_i] = di[0];
    }

    free(lp_hs);
    return 0;
}

int mdl_fit(bam_data_t *bam_dat, obj_pars *objs){
    if (bam_dat == NULL || objs == NULL)
        return err_msg(-1, 0, "mdl_fit: arguments are NULL");

    if (bam_dat->bc_data == NULL)
        return err_msg(-1, 0, "mdl_fit: rna and atac are NULL");

    if (objs->gv == NULL || objs->vcf_hdr == NULL)
        return err_msg(-1, 0, "mdl_fit: variant data  is NULL");

    if (objs->verbose) log_msg("initializing model");

    mdl_t *mdl = mdl_alloc();
    if (mdl == NULL)
        return -1;

    mdl->eps = objs->eps;
    mdl->max_iter = objs->max_iter;
    mdl->threads = objs->threads;
    mdl->alpha_vars = objs->alpha_vars > 0;
    mdl->alpha_max = objs->alpha_max;

    if (mdl_set_rna_atac(mdl, bam_dat->has_rna, bam_dat->has_atac) < 0)
        return -1;

    if (mdl_bc_dat_bam_data(mdl->mdl_bc_dat, bam_dat, objs) < 0)
        return -1;

    if (mdl_set_samples(mdl, objs->vcf_hdr) < 0)
        return -1;

    if (mdl_set_k(mdl, objs->k) < 0)
        return -1;

    if (mdl_set_hs_ix(mdl) < 0)
        return -1;

    if (mdl_alloc_probs(mdl) < 0)
        return -1;

    // initialize parameters
    if (mdl_pars_set_dat(mdl->mp, mdl->mdl_bc_dat, objs,
                mdl->samples->n, mdl->k) < 0)
        return -1;

    // run sub model
    if (mdl_em(mdl, objs, 1) < 0)
        return -1;

    // get best (H,S) index
    if (mdl_sub_best_hs(mdl) < 0)
        return -1;

    // run full model
    if (mdl_em(mdl, objs, 0) < 0)
        return -1;

    if (mdl_get_prob_h(mdl) < 0)
        return -1;

    if (objs->verbose)
        log_msg("getting barcode likelihoods");
    if (mdl_get_best(mdl) < 0)
        return -1;

    // save output
    if (objs->verbose)
        log_msg("writing out likelihoods");
    if (write_llk(mdl, objs->out_fn) < 0)
        return -1;
    if (write_samples(mdl, objs->out_fn) < 0)
        return -1;
    if (objs->verbose)
        log_msg("writing out parameters");
    if (write_lambda(mdl, objs->out_fn) < 0)
        return -1;
    if (bam_dat->has_rna && write_alpha_rna(mdl, objs->out_fn) < 0)
        return -1;
    if (bam_dat->has_atac && write_alpha_atac(mdl, objs->out_fn) < 0)
        return -1;
    if (bam_dat->has_rna && write_rho(mdl, objs, objs->out_fn) < 0)
        return -1;
    if (bam_dat->has_atac && write_sigma(mdl, objs, objs->out_fn) < 0)
        return -1;

    if (objs->verbose)
        log_msg("writing out summary results");
    if (write_res(mdl, bam_dat, objs->out_fn) < 0)
        return -1;

    // correct counts
    if (objs->verbose)
        log_msg("correcting counts");
    iregs_t *pks = NULL;
    if (objs->pks)
        pks = objs->pks;
    char *corr_fn = ".corrected";
    char *corr_pre = strcat2((const char*)(objs->out_fn), (const char*)corr_fn);
    if (mod_correct_counts(mdl, objs->gv, objs->anno, pks, corr_pre) < 0)
        return -1;
    free(corr_pre);
    if (objs->verbose)
        log_msg("model fit finished");

    mdl_dstry(mdl);

    return(0);
}

char **mdl_s_names(mdl_t *mdl) {
    if (mdl == NULL) {
        err_msg(-1, 0, "mdl_s_names: argument is null");
        return NULL;
    }

    uint32_t n_ixs = mdl->hs_ix->n_ixs;

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    char **col_names = malloc(n_ixs * sizeof(char *));
    if (col_names == NULL) {
        err_msg(-1, 0, "mdl_s_names: %s", strerror(errno));
        return NULL;
    }
    char ssep[] = "-";
    uint32_t hs_ix;
    for (hs_ix = 0; hs_ix < n_ixs; ++hs_ix){
        if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
            return NULL;
        char *s1, *s2;
        
        if (par_ix.hd == 0){
            if ( (col_names[hs_ix] = strdup("Empty")) == NULL ) {
                err_msg(-1, 0, "mdl_s_names: %s", strerror(errno));
                return NULL;
            }
        } else if (par_ix.hd == 1) {
            s1 = str_map_str(mdl->samples, par_ix.s1);
            if (s1 == NULL) {
                err_msg(-1, 0, "mdl_s_names: failed to get sample for index %i", par_ix.s1);
                return NULL;
            }
            if ( (col_names[hs_ix] = strdup(s1)) == NULL ) {
                err_msg(-1, 0, "mdl_s_names: %s", strerror(errno));
                return NULL;
            }
        } else {
            s1 = str_map_str(mdl->samples, par_ix.s1);
            if (s1 == NULL) {
                err_msg(-1, 0, "mdl_s_names: failed to get sample for index %i", par_ix.s1);
                return NULL;
            }
            s2 = str_map_str(mdl->samples, par_ix.s2);
            if (s2 == NULL) {
                err_msg(-1, 0, "mdl_s_names: failed to get sample for index %i", par_ix.s2);
                return NULL;
            }
            const char *sa[2] = {s1, s2};
            if ( (col_names[hs_ix] = cat_strs(sa, 2, ssep)) == NULL ) {
                err_msg(-1, 0, "mdl_s_names: failed to cat samples '%s' and '%s'", s1, s2);
                return NULL;
            }
        }
    }

    return col_names;
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

int write_alpha_rna(mdl_t *mdl, char *fn){
    if (mdl == NULL)
        return err_msg(-1, 0, "write_alpha_rna: arguments are NULL");

    if (mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "write_alpha_rna: no barcodes found");
    if (mdl->mp == NULL || mdl->mp->alpha_rna == NULL)
        return err_msg(-1, 0, "write_alpha_rna: model hasn't been initialized");
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_ixs = mdl->hs_ix->n_hs;
    int n_test_bc = mdl->mp->T;
    uint32_t n_all_bc = mdl->mp->D;

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    char *alpha_fn = ".alpha_rna.txt.gz";
    char *out_alpha_fn = strcat2((const char*)fn, (const char*)alpha_fn);

    // row names
    char **bc_row_names = str_map_ca(bd->test_bcs);
    assert(n_test_bc == bd->test_bcs->n);

    // alpha array for test barcodes
    f_t *al = malloc(n_test_bc * n_ixs * sizeof(f_t));
    int i;
    uint32_t j;
    for (i = 0; i < n_test_bc; ++i){
        char *bc = str_map_str(bd->test_bcs, i);
        int bc_ix = str_map_ix(bd->all_bcs, bc);
        for (j = 0; j < n_ixs; ++j)
            al[CMI(i, j, n_test_bc)] = mdl->mp->alpha_rna[CMI(bc_ix, j, n_all_bc)];
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

    if (mdl->mp == NULL || mdl->mp->alpha_atac == NULL)
        return err_msg(-1, 0, "write_alpha_atac: model hasn't been initialized");
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_ixs = mdl->hs_ix->n_hs;
    int n_test_bc = mdl->mp->T;
    uint32_t n_all_bc = mdl->mp->D;

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    char *alpha_fn = ".alpha_atac.txt.gz";
    char *out_alpha_fn = strcat2((const char*)fn, (const char*)alpha_fn);

    // row names
    char **bc_row_names = str_map_ca(bd->test_bcs);
    assert(n_test_bc == bd->test_bcs->n);

    // alpha array for test barcodes
    f_t *al = malloc(n_test_bc * n_ixs * sizeof(f_t));
    int i;
    uint32_t j;
    for (i = 0; i < n_test_bc; ++i){
        char *bc = str_map_str(bd->test_bcs, i);
        int bc_ix = str_map_ix(bd->all_bcs, bc);
        for (j = 0; j < n_ixs; ++j)
            al[CMI(i, j, n_test_bc)] = mdl->mp->alpha_atac[CMI(bc_ix, j, n_all_bc)];
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

int write_rho(mdl_t *mdl, obj_pars *objs, char *fn){
    if (mdl == NULL || objs == NULL || fn == NULL)
        return err_msg(-1, 0, "write_rho: arguments are NULL");

    if (mdl->mp == NULL || mdl->mp->rho == NULL)
        return err_msg(-1, 0, "write_rho: model hasn't been initialized");

    if (objs->anno == NULL)
        return err_msg(-1, 0, "write_rho: annotation not provided");

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    char *rho_fn = ".rho.txt.gz";
    char *out_rho_fn = strcat2((const char*)fn, (const char*)rho_fn);

    // create row names for genes
    uint32_t i, gs, g, G = mdl->mp->G, GA = 3 * G;
    uint32_t rho_ncol = 1 + mdl->mp->K;
    char r_apnd[3][20] = {"_SPLICE", "_UNSPLICE", "_AMBIG"};
    char **gene_row_names = malloc(GA * sizeof(char *));
    if (gene_row_names == NULL)
        return err_msg(-1, 0, "write_rho: %s", strerror(errno));
    for (gs = 0; gs < 3; ++gs){
        for (g = 0; g < G; ++g){
            const char *g_str = str_map_str(objs->anno->gene_ix, g);
            char *gspl_str = strcat2(g_str, (const char*)r_apnd[gs]);
            if (gspl_str == NULL)
                return err_msg(-1, 0, "write_rho: %s", strerror(errno));
            gene_row_names[g + (gs * G)] = gspl_str;
        }
    }

    // column names
    char **col_names = malloc(rho_ncol * sizeof(char*));
    col_names[0] = strdup("Ambient");
    if (col_names[0] == NULL)
        return err_msg(-1, 0, "write_rho: %s", strerror(errno));
    size_t buf_size = 100;
    for (i = 1; i < rho_ncol; ++i) {
        col_names[i] = malloc(buf_size * sizeof(char));
        int2strp(i, col_names + i, &buf_size);
    }

    // write matrix
    ret = write_matrix_double(out_rho_fn, mdl->mp->rho, NULL, NULL, NULL, 
            gene_row_names, GA, col_names, rho_ncol, 
            delim, nl, decp);
    if (ret < 0)
        return err_msg(-1, 0, "write_rho: failed to write matrix to file");

    for (i = 0; i < rho_ncol; ++i)
        free(col_names[i]);
    free(col_names);
    for (g = 0; g < GA; ++g)
        free(gene_row_names[g]);
    free(gene_row_names);
    free(out_rho_fn);

    return(0);
}

int write_sigma(mdl_t *mdl, obj_pars *objs, char *fn){
    if (mdl == NULL || objs == NULL || fn == NULL)
        return err_msg(-1, 0, "write_sigma: arguments are NULL");

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "write_sigma: model hasn't been initialized");

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    char *sigma_fn = ".sigma.txt.gz";
    char *out_sigma_fn = strcat2((const char*)fn, (const char*)sigma_fn);

    // row names
    char **row_names = malloc(2 * sizeof(char *));
    row_names[0] = strdup("Ambient");
    if (row_names[0] == NULL)
        return err_msg(-1, 0, "write_sigma: %s", strerror(errno));
    row_names[1] = strdup("Cell");
    if (row_names[1] == NULL)
        return err_msg(-1, 0, "write_sigma: %s", strerror(errno));

    // col names
    char **col_names = malloc(sizeof(char *));
    col_names[0] = strdup("Sigma");
    if (col_names[0] == NULL)
        return err_msg(-1, 0, "write_sigma: %s", strerror(errno));
    
    // write matrix
    ret = write_matrix_double(out_sigma_fn, mdl->mp->sigma, NULL, NULL, NULL, 
            row_names, 2, col_names, 1, 
            delim, nl, decp);
    if (ret < 0)
        return err_msg(-1, 0, "write_sigma: failed to write matrix to file");

    free(col_names[0]);
    free(col_names);
    free(row_names[0]);
    free(row_names[1]);
    free(row_names);
    free(out_sigma_fn);

    return(0);
}

int write_llk(mdl_t *mdl, char *fn){
    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    uint32_t i, j;
    int ret = 0;
    char *llk_fn = ".llk.txt.gz";

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_hsk = mdl->hs_ix->n_ixs;
    uint32_t D = mdl->mp->D;
    uint32_t T = mdl->mp->T;

    char *out_llk_fn = strcat2((const char*)fn, (const char*)llk_fn);
    if (out_llk_fn == NULL)
        return err_msg(-1, 0, "write_llk: %s", strerror(errno));

    // row names
    char **bc_row_names = str_map_ca(bd->test_bcs);
    if (bc_row_names == NULL)
        return err_msg(-1, 0, "write_llk: failed to write to file");

    // get llk matrix
    f_t *llks = malloc(T * n_hsk * sizeof(f_t));
    for (i = 0; i < T; ++i) {
        char *bc_name = str_map_str(bd->test_bcs, i);
        int all_bc_ix = str_map_ix(bd->all_bcs, bc_name);
        assert(all_bc_ix >= 0);
        for (j = 0; j < n_hsk; ++j) {
            llks[CMI(i, j, T)] = mdl->lp_hsk[CMI(all_bc_ix, j, D)];
        }
    }

    // write matrix
    ret = write_matrix_double(out_llk_fn, llks, NULL, NULL, NULL, 
            bc_row_names, T, NULL, n_hsk, 
            delim, nl, decp);
    if (ret < 0)
        return err_msg(-1, 0, "write_llk: failed to write matrix to file");

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

    uint32_t n_hsk = mdl->hs_ix->n_ixs;
    par_ix_t par_ix;
    par_ix_init(&par_ix);
    int fret;
    uint32_t hs_ix;
    for (hs_ix = 1; hs_ix < n_hsk; ++hs_ix){
        if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
            return -1;

        // if singlet
        char *h_str = aloc_sprintf("%i", par_ix.hd);
        if (h_str == NULL)
            return -1;
        char *k_str = aloc_sprintf("%i", par_ix.k);
        if (k_str == NULL)
            return -1;
        if (par_ix.hd == 0) {
            fret = fputs(h_str, fp);
            fret = fputc(delim, fp);
            fret = fputs("NA", fp);
            fret = fputc(delim, fp);
            fret = fputs("NA", fp);
            fret = fputc(delim, fp);
            fret = fputs("Empty", fp);
        }
        if (par_ix.hd == 1){ // if singlet
            fret = fputs(h_str, fp);
            fret = fputc(delim, fp);
            fret = fputs(mdl->samples->strs[par_ix.s1], fp);
            fret = fputc(delim, fp);
            fret = fputs("NA", fp);
            fret = fputc(delim, fp);
            fret = fputs(k_str, fp);
        } else { // if doublet
            fret = fputs(h_str, fp);
            fret = fputc(delim, fp);
            fret = fputs(mdl->samples->strs[par_ix.s1], fp);
            fret = fputc(delim, fp);
            fret = fputs(mdl->samples->strs[par_ix.s2], fp);
            fret = fputc(delim, fp);
            fret = fputs(k_str, fp);
        }
        fret = fputc(nl, fp);
        free(h_str);
        free(k_str);
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

    if (mdl == NULL || bam_dat == NULL)
        return err_msg(-1, 0, "write_res: argument is null");
    if (mdl->best_sng_ix == NULL)
        return err_msg(-1, 0, "write_res: barcodes not classified yet, "
                "'mdl_get_best' must be run");


    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_hsk = mdl->hs_ix->n_ixs;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    uint32_t n_test_bc = mdl->mp->T;
    uint32_t D = mdl->mp->D;
    uint32_t n_all_bc = mdl->mp->D;
    char **sam_ids = mdl_s_names(mdl);
    str_map *samples = mdl->samples;
    str_map *test_bcs = bd->test_bcs;
    par_ix_t par_ix;
    par_ix_init(&par_ix);

    f_t *lp_hs = mdl_get_lp_hs(mdl);
    if (lp_hs == NULL)
        return err_msg(-1, 0, "write_res: failed to get lp_hs");

    if (sam_ids == NULL)
        return err_msg(-1, 0, "write_res: failed to get sample IDs");

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
        "\tn_rna_variants\tn_atac_variants\tatac_pct_mt\trna_pct_mt\tFRIP"
        "\tbest_type\tbest_sample\tbest_cluster"
        "\tbest_rna_alpha\tbest_atac_alpha"
        "\tPP0\tPP1\tPP2"
        "\tLLK0\tLLK1\tLLK2\n";

    fputs(hdr, fp);

    size_t buf_size = decp + 1000;
    char *pstr = (char *)malloc(buf_size * sizeof(char));
    if (pstr == NULL)
        return err_msg(-1, 0, "write_res: %s", strerror(errno));

    int pstr_len;
    uint32_t bc_i, i, j;
    for (bc_i = 0; bc_i < n_test_bc; ++bc_i){
        char *bc_name = str_map_str(test_bcs, bc_i);
        fret = fputs(bc_name, fp);
        int all_bc_ix = str_map_ix(bd->all_bcs, bc_name);
        assert(all_bc_ix >= 0);

        double par;

        khint_t k_bc = kh_get(kh_bc_dat, bam_dat->bc_data, bc_name);
        if (k_bc == kh_end(bam_dat->bc_data))
            return err_msg(-1, 0, "write_res: barcode %s not found", bc_name);

        bc_data_t *bc_data = kh_val(bam_dat->bc_data, k_bc);
        bc_stats_t *bcc = bc_data->bc_stats;

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
        int2strp(bcc->n_rna_vars, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        int2strp(bcc->n_atac_vars, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        double2str_in((double)bcc->atac_mt, &pstr, &buf_size, 4);
        fputs(pstr, fp);

        fputc(delim, fp);
        double2str_in((double)bcc->rna_mt, &pstr, &buf_size, 4);
        fputs(pstr, fp);

        fputc(delim, fp);
        double2str_in((double)bcc->frip, &pstr, &buf_size, 4);
        fputs(pstr, fp);

        // get best assignment
        uint32_t best_hsk = 0;
        for (j = 1; j < n_hsk; ++j) {
            if (mdl->lp_hsk[CMI(all_bc_ix, j, D)] > mdl->lp_hsk[CMI(all_bc_ix, best_hsk, D)])
                best_hsk = j;
        }
        uint32_t best_hs = 0;
        for (j = 1; j < n_hs; ++j) {
            if (lp_hs[CMI(all_bc_ix, j, D)] > lp_hs[CMI(all_bc_ix, best_hs, D)])
                best_hs = j;
        }
        if (hs_ix_get_pars(mdl->hs_ix, best_hsk, &par_ix) < 0)
            return -1;

        // write best H_d
        fret = fputc(delim, fp);
        if (par_ix.hd == 0) {
            fret = fputs("Empty", fp);
        } else if (par_ix.hd == 1) {
            fret = fputs("Singlet", fp);
        } else {
            fret = fputs("Doublet", fp);
        }

        // get best S_d
        fret = fputc(delim, fp);
        if (best_hsk == 0) {
            fputs("Empty", fp);
        } else {
            fputs(str_map_str(samples, par_ix.s1), fp);
            if (par_ix.s2 >= 0) {
                fputs(":", fp);
                fputs(str_map_str(samples, par_ix.s2), fp);
            }
        }

        // get best K_d
        fret = fputc(delim, fp);
        if (best_hsk == 0) {
            fputs("Empty", fp);
        } else {
            if (int2strp(par_ix.k, &pstr, &buf_size) < 0)
                return -1;
            fputs(pstr, fp);
        }

        // write best singlet alpha RNA
        f_t alpha_rna = mdl->mp->alpha_rna[CMI(all_bc_ix, best_hs, n_all_bc)];
        par = alpha_rna;
        if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
            return err_msg(-1, 0, "write_res: failed to convert %f to string", par);
        fret = fputc(delim, fp);
        fret = fputs(pstr, fp);

        // write best singlet alpha ATAC
        f_t alpha_atac = mdl->mp->alpha_atac[CMI(all_bc_ix, best_hs, n_all_bc)];
        par = alpha_atac;
        if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
            return err_msg(-1, 0, "write_res: failed to convert %f to string", par);
        fret = fputc(delim, fp);
        fret = fputs(pstr, fp);

        // write posterior probs of H_d
        for (i = 0; i < 3; ++i){
            par = mdl->pp[CMI(bc_i, i, n_test_bc)];
            if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
                return err_msg(-1, 0, "write_res: failed to convert %f to string", par);
            fret = fputc(delim, fp);
            fret = fputs(pstr, fp);
        }

        // write llks of H_d
        for (i = 0; i < 3; ++i){
            par = mdl->lp_h[CMI(bc_i, i, n_test_bc)];
            if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
                return err_msg(-1, 0, "write_res: failed to convert %f to string", par);
            fret = fputc(delim, fp);
            fret = fputs(pstr, fp);
        }

        // write log likelihood ratio of sng to doublet
        fret = fputc(nl, fp);
    }

    free(pstr);
    uint32_t hs_ix;
    for (hs_ix = 0; hs_ix < n_hsk; ++hs_ix)
        free(sam_ids[hs_ix]);
    free(sam_ids);
    free(lp_hs);

    // close file
    if ( (fret = fclose(fp)) != 0 )
        return err_msg(-1, 0, "write_res: could not close file %s: %s", 
                out_r_fn, strerror(errno));
    free(out_r_fn);

    return 0;
}

/*******************************************************************************
 * count correction
 ******************************************************************************/

int mod_correct_counts(mdl_t *mdl, g_var_t *gv, gene_anno_t *gene_anno,
        iregs_t *pks, char *out_fn) {

    if (mdl == NULL)
        return err_msg(-1, 0, "mod_correct_counts: argument is null");

    if (mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "mod_correct_counts: barcode data not present");
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;

    if (gv == NULL)
        return err_msg(-1, 0, "mod_correct_counts: 'gv' variants are required");

    bam_counts_t *bam_counts = bam_counts_init();
    if (bam_counts == NULL)
        return -1;

    if (bam_counts_add_gv(bam_counts, gv) < 0)
        return -1;


    int32_t n_genes = 0;
    if (gene_anno) {
        if (bam_counts_add_gene_map(bam_counts, gene_anno->gene_ix) < 0)
            return -1;
        n_genes = gene_anno->gene_ix->n;
    }

    int32_t n_vars = mv_size(&gv->vix2var);

    str_map *all_bcs = bd->all_bcs;
    str_map *test_bcs = bd->test_bcs;
    if (all_bcs->n < 1)
        return err_msg(-1, 0, "mod_correct_counts: no barcodes present");
    if (test_bcs->n < 1)
        return err_msg(-1, 0, "mod_correct_counts: no test barcodes present");
    if (bam_counts_add_bc_map(bam_counts, test_bcs) < 0)
        return -1;

    if (pks) {
        if (bam_counts_add_peaks(bam_counts, pks) < 0)
            return -1;
    }

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    uint32_t D = all_bcs->n;
    uint32_t n_hsk = mdl->hs_ix->n_ixs;

    // loop through test bcs
    int bc_i, n_test_bc = test_bcs->n;
    for (bc_i = 0; bc_i < n_test_bc; ++bc_i) {
        char *bc_key = str_map_str(test_bcs, bc_i);
        assert(bc_key);
        int d_bc_i = str_map_ix(all_bcs, bc_key);
        assert(d_bc_i >= 0);

        // add barcode to kbtree
        int ret;
        khint_t k_bt = kh_put(kh_bc_cnode, bam_counts->bc_counts, bc_key, &ret);
        if (ret == -1){
            return err_msg(-1, 0, "mod_correct_counts: failed to add barcode %s to "
                    "hash table", bc_key);
        } else if (ret == 0){
            return err_msg(-1, 0, "mod_correct_counts: barcode %s found twice\n", bc_key);
        }
        bc_counts_t *bcc = bc_counts_init();
        if (bcc == NULL) return(-1);

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, d_bc_i);

        // bam data
        kbtree_t(kb_mdl_mlcl) *mols = bc_dat.rna;
        kbitr_t itr;

        // loop over RNA
        kb_itr_first(kb_mdl_mlcl, mols, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);

            // skip UMIs with no gene or variant
            if (mv_size(&mlcl->feat_ixs) == 0 && mv_size(&mlcl->var_ixs) == 0)
                continue;

            // skip multi-gene UMIs
            if (mv_size(&mlcl->feat_ixs) > 1)
                continue;

            // hold the amb posterior prob.
            f_t amb_pp = 0.0;

            uint32_t hsk_ix;
            for (hsk_ix = 0; hsk_ix < n_hsk; ++hsk_ix) {
                // Pr(H_d, S_d, K_d | X_d)
                f_t cp_hsk_i = exp(mdl->lp_hsk[CMI(d_bc_i, hsk_ix, D)] - mdl->p_x[d_bc_i]);

                // skip over low prob assignments
                if (cp_hsk_i < 0.01)
                    continue;

                if (hs_ix_get_pars(mdl->hs_ix, hsk_ix, &par_ix) < 0)
                    return -1;
                f_t psum = 0, p_tfb[3] = {1,1,1};
                int pret = p_f_v(mdl->mp, mlcl, RNA_IX, d_bc_i, &par_ix, p_tfb, &psum);
                if (pret < 0)
                    return -1;
                if (prob_invalid(psum))
                    return err_msg(-1, 0, "mdl_correct_counts: p_f_v returned sum %f", psum);
                assert(psum > 0 && psum <= 1);

                // remove ambient portion of best assignment
                int t_im = 0;
                f_t pt = p_tfb[t_im] / psum; // Pr(T_dm = a | X_d)
                f_t pp = cp_hsk_i * pt; // Pr(H_d, S_d, K_d, T_dm = a | X_d)
                assert(!prob_invalid(pp));
                amb_pp += pp;
            }
            assert(!prob_invalid(amb_pp));
            f_t all_reads = (f_t)mlcl->counts;
            f_t nuc_reads = all_reads * (1.0 - amb_pp);

            uint32_t mv_ix;
            // add feature read count
            for (mv_ix = 0; mv_ix < mv_size(&mlcl->feat_ixs); ++mv_ix) {
                int32_t fs_ix, s_ix, f_ix;
                fs_ix = mv_i(&mlcl->feat_ixs, 0);
                s_ix = fs_ix / n_genes;
                f_ix = fs_ix % n_genes;

                cnt_node_t *p, t;
                memset(&t, 0, sizeof(cnt_node_t));
                t.ix = (int)f_ix;
                t.counts[s_ix] = nuc_reads;
                p = kb_getp(kh_cnode, bcc->rna_gc, &t);
                if (!p){
                    kb_putp(kh_cnode, bcc->rna_gc, &t);
                    ++bam_counts->rna_gcs_nz;
                } else {
                    p->counts[s_ix] += nuc_reads;
                }
                bam_counts->has_rna_gc = 1;
            }

            // add variant read count
            for (mv_ix = 0; mv_ix < mv_size(&mlcl->var_ixs); ++mv_ix) {
                int32_t va_ix, v_ix, a_ix;
                va_ix = mv_i(&mlcl->var_ixs, mv_ix);
                v_ix = va_ix % n_vars;
                a_ix = va_ix / n_vars;
                if (a_ix > MAX_ALLELE || a_ix < 0)
                    return err_msg(-1, 0, "mod_correct_counts: a_ix=%i", a_ix);

                cnt_node_t *p, t;
                memset(&t, 0, sizeof(cnt_node_t));
                t.ix = (int)v_ix;
                t.counts[a_ix] = nuc_reads;
                p = kb_getp(kh_cnode, bcc->rna_ac, &t);
                if (!p){
                    kb_putp(kh_cnode, bcc->rna_ac, &t);
                    ++bam_counts->rna_acs_nz;
                } else {
                    p->counts[a_ix] += nuc_reads;
                }
                bam_counts->has_rna_ac = 1;
            }
        }

        /*
        // loop over ATAC
        kbtree_t(kb_mdl_mlcl) *frags = bc_dat.atac;
        kb_itr_first(kb_mdl_mlcl, frags, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, frags, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            
            // skip fragments with no peak or variant
            if (mv_size(&mlcl->feat_ixs) == 0 && mv_size(&mlcl->var_ixs) == 0)
                continue;

            // skip multi-peak fragments
            if (mv_size(&mlcl->feat_ixs) > 1)
                continue;

            // hold the amb posterior prob.
            f_t amb_pp = 0.0;

            uint32_t hsk_ix;
            for (hsk_ix = 0; hsk_ix < n_hsk; ++hsk_ix) {
                // Pr(H_d, S_d, K_d | X_d)
                f_t cp_hsk_i = exp(mdl->lp_hsk[CMI(d_bc_i, hsk_ix, D)] - mdl->p_x[d_bc_i]);

                // skip over low prob assignments
                if (cp_hsk_i < 0.01)
                    continue;

                if (hs_ix_get_pars(mdl->hs_ix, hsk_ix, &par_ix) < 0)
                    return -1;
                f_t psum = 0, p_tfb[3] = {1,1,1};
                int pret = p_f_v(mdl->mp, mlcl, 0, d_bc_i, &par_ix, p_tfb, &psum);
                if (pret < 0)
                    return -1;
                if (prob_invalid(psum))
                    return err_msg(-1, 0, "mdl_full_e: p_f_v returned sum %f", psum);
                assert(psum > 0 && psum <= 1);

                // remove ambient portion of best assignment
                int t_im = 0;
                f_t pt = p_tfb[t_im] / psum; // Pr(T_dm = a | X_d)
                f_t pp = cp_hsk_i * pt; // Pr(H_d, S_d, K_d, T_dm = a | X_d)
                assert(!prob_invalid(pp));
                amb_pp += pp;
            }
            assert(!prob_invalid(amb_pp));
            f_t all_reads = (f_t)mlcl->counts;
            f_t nuc_reads = all_reads * (1.0 - amb_pp);

            uint32_t mv_ix;
            // add feature read count
            for (mv_ix = 0; mv_ix < mv_size(&mlcl->feat_ixs); ++mv_ix) {
                int32_t fs_ix, s_ix, f_ix;
                fs_ix = mv_i(&mlcl->feat_ixs, 0);
                s_ix = fs_ix / n_genes;
                f_ix = fs_ix % n_genes;

                cnt_node_t *p, t;
                memset(&t, 0, sizeof(cnt_node_t));
                t.ix = (int)f_ix;
                t.counts[s_ix] = nuc_reads;
                p = kb_getp(kh_cnode, bcc->rna_gc, &t);
                if (!p){
                    kb_putp(kh_cnode, bcc->rna_gc, &t);
                    ++bam_counts->rna_gcs_nz;
                } else {
                    p->counts[s_ix] += nuc_reads;
                }
                bam_counts->has_rna_gc = 1;
            }

            // add variant read count
            for (mv_ix = 0; mv_ix < mv_size(&mlcl->var_ixs); ++mv_ix) {
                int32_t va_ix, v_ix, a_ix;
                va_ix = mv_i(&mlcl->var_ixs, mv_ix);
                v_ix = va_ix / n_vars;
                a_ix = va_ix % n_vars;

                cnt_node_t *p, t;
                memset(&t, 0, sizeof(cnt_node_t));
                t.ix = (int)v_ix;
                t.counts[a_ix] = nuc_reads;
                p = kb_getp(kh_cnode, bcc->rna_ac, &t);
                if (!p){
                    kb_putp(kh_cnode, bcc->rna_ac, &t);
                    ++bam_counts->rna_acs_nz;
                } else {
                    p->counts[a_ix] += nuc_reads;
                }
                bam_counts->has_rna_ac = 1;
            }
        }
        */

        kh_val(bam_counts->bc_counts, k_bt) = bcc;
    }

    if (bam_counts_write(bam_counts, gene_anno, gv, out_fn, 1) < 0)
        return -1;

    bam_counts_dstry(bam_counts);
    return 0;
}
