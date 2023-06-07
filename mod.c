
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
    int32_t in_pk = mv_size(&frag->pks) > 0;
    int32_t pk_ix = in_pk ? mv_i(&frag->pks, 0) + 1 : 0;
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
    assert(ixi == hs_ix->n_hs);
}

int hs_ix_get_pars(hs_ix_t *hs_ix, unsigned ix,
        par_ix_t *par_ix) {

    int n_par = 3;
    if (ix >= hs_ix->n_hs)
        return err_msg(-1, 0, "ix=%i > n_ixs=%u", ix, hs_ix->n_hs);

    // -1 for if NA/invalid
    par_ix->hs_ix = ix;
    uint32_t divby = hs_ix->n_hs - 1;
    if (ix == 0)
        par_ix->hs_ix = ix;
    else
        par_ix->hs_ix = ((ix - 1) % divby) + 1;
    assert((uint32_t)par_ix->hs_ix < hs_ix->n_hs);
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
    mp->G = 0;
    mp->V = 0;
    mp->M = 0;
    mp->M = 0;
    mp->n_hs = 0;

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
    free(mp->sigma);
    free(mp->gamma);
    
    pthread_mutex_destroy(&mp->sum_lock);
    free(mp->_pi_sum);
    free(mp->_rho_sum);
    free(mp->_sigma_sum);
}

void mdl_pars_dstry(mdl_pars_t *mp){
    if (mp == NULL) return;
    mdl_pars_free(mp);
    free(mp);
}

int mld_pars_set_num_alloc(mdl_pars_t *mp, uint32_t D, uint32_t T,
        uint32_t G, uint32_t P, uint32_t V, uint16_t M, uint16_t K) {
    if (mp == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: argument is null");

    mp->D = D;
    mp->T = T;
    mp->G = G;
    mp->P = P;
    mp->V = V;
    mp->M = M;
    mp->K = K;

    log_msg("D=%u, T=%u, G=%u, P=%u, V=%u, M=%u, K=%u",
            D, T, G, P, V, M, K);

    mp->n_hs = 1 + (mp->M) + (mp->M * (mp->M - 1) / 2);

    // allocate parameter fields
    mp->pi = calloc(mp->M, sizeof(f_t));
    if (mp->pi == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->kappa = calloc(mp->D * mp->K, sizeof(f_t));
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

    mp->sigma = calloc(mp->P * (mp->K + 1), sizeof(f_t));
    if (mp->sigma == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->gamma = calloc((mp->M + 1) * mp->V, sizeof(f_t));
    if (mp->gamma == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    // allocate sum fields
    mp->_pi_sum = calloc(mp->M, sizeof(f_t));
    if (mp->_pi_sum == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->_rho_sum = calloc((mp->G * 3) * (mp->K + 1), sizeof(f_t));
    if (mp->_rho_sum == NULL)
        return err_msg(-1, 0, "mld_pars_set_num_alloc: %s", strerror(errno));

    mp->_sigma_sum = calloc(mp->P * (mp->K + 1), sizeof(f_t));
    if (mp->_sigma_sum == NULL)
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

    uint32_t rho_ne = (mp->G * 3) * (mp->K + 1);
    for (i = 0; i < rho_ne; ++i)
        mp->_rho_sum[i] = 0;

    uint32_t sigma_ne = mp->P * (mp->K + 1);
    for (i = 0; i < sigma_ne; ++i)
        mp->_sigma_sum[i] = 0;

    return 0;
}

int mdl_pars_pi_fix(mdl_pars_t *mp){

    uint32_t i, j;
    uint32_t M = mp->M;

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

    uint32_t i;
    f_t ku = 1.0 / (f_t)(mp->K);
    uint32_t kappa_ne = (mp->D * mp->K);
    for (i = 0; i < kappa_ne; ++i)
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
                bd->G, bd->P, bd->V, n_samples, K) < 0)
        return -1;

    uint32_t rho_nrow = 3 * mp->G;
    uint32_t rho_ncol = 1 + mp->K;
    uint32_t sig_nrow = mp->P;
    uint32_t sig_ncol = 1 + mp->K;

    // keep pi uniform
    if (mdl_pars_pi_fix(mp) < 0)
        return(-1);

    // set kappa uniform
    if (mdl_pars_kappa_uni(mp) < 0)
        return -1;

    // arrays to sum reads
    uint32_t rho_ne = rho_nrow * rho_ncol;
    f_t *rho_sum = calloc(rho_ne, sizeof(f_t));

    uint32_t sig_ne = sig_nrow * sig_ncol;
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
    f_t *sig_tot = calloc(sig_ncol, sizeof(f_t));
    for (i = 0; i < sig_ncol; ++i)
        sig_tot[i] = psc * sig_nrow;

    f_t n_amb_bcs = 1; // prior of 1 for lambda
    f_t n_nuc_bcs = 1;

    for (i = 0; i < mp->D; ++i){
        int fl;

        // if no reads
        fl = bflg_get(bd->absent_bc, i);
        if (fl == 1) continue;

        // if ambient
        fl = bflg_get(bd->amb_flag, i);

        // sample cell type
        int rand_k = rand() % mp->K;
        uint32_t ct_col = fl ? 0 : rand_k + 1;

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
                uint32_t rho_ix = CMI(rho_row, ct_col, rho_nrow);
                rho_sum[rho_ix] += (f_t)mlcl->counts;
                rho_tot[ct_col] += (f_t)mlcl->counts;
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
                uint32_t six = CMI(pk, ct_col, sig_nrow);
                sig_sum[six] += (f_t)mlcl->counts;
                sig_tot[ct_col] += (f_t)mlcl->counts;
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
    for (i = 0; i < sig_ncol; ++i) {
        for (j = 0; j < sig_nrow; ++j) {
            uint32_t sig_ix = CMI(j, i, sig_nrow);
            mp->sigma[sig_ix] = sig_sum[sig_ix] / sig_tot[i];
        }
    }

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
    free(sig_tot);
    
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
    for (i = 0; i < mp->D; ++i){
        f_t kappa_sum = 0;
        for (j = 0; j < mp->K; ++j) {
            uint32_t mix = CMI(j, i, mp->K);
            kappa_sum += mp->kappa[mix];
            if (prob_invalid(mp->kappa[mix]))
                return err_msg(-1, 0, "mdl_pars_check: kappa[%i] = %f",
                        i, mp->kappa[mix]);
        }
        if (psum_invalid(kappa_sum, lt1, ut1))
            return err_msg(-1, 0, "mdl_pars_check: kappa sum %i = %f",
                    i, kappa_sum);
    }

    // if (psum_invalid(kappa_sum, lt1, ut1))
    //     return err_msg(-1, 0, "mdl_pars_check: kappa sum = %f", kappa_sum);

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
    for (i = 0; i < rho_ncol; ++i){
        if (psum_invalid(rho_tot[i], lt1, ut1))
            return err_msg(-1, 0, "mdl_pars_check: rho sum %i = %f",
                    i, rho_tot[i]);
    }

    // sigma
    uint32_t sig_nrow = mp->P;
    uint32_t sig_ncol = 1 + mp->K;
    f_t *sig_tot = calloc(sig_ncol, sizeof(f_t));
    for (i = 0; i < sig_ncol; ++i)
        sig_tot[i] = 0;
    for (i = 0; i < sig_ncol; ++i){
        for (j = 0; j < sig_nrow; ++j){
            sig_tot[i] += mp->sigma[CMI(j, i, sig_nrow)];
            f_t p = mp->sigma[CMI(j, i, sig_nrow)];
            if (prob_invalid(p))
                return err_msg(-1, 0, "mdl_pars_check: sig[%i,%i] = %f",
                        i, j, p);
        }
    }
    for (i = 0; i < sig_ncol; ++i){
        if (psum_invalid(sig_tot[i], lt1, ut1))
            return err_msg(-1, 0, "mdl_pars_check: sig sum %i = %f",
                    i, sig_tot[i]);
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
    free(sig_tot);

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
    mdl_bc_dat->P = 0;
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
    mdl_bc_dat->P = 0;
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
    uint32_t n_peaks = 0;
    if (objs->pks)
        n_peaks = objs->pks->n + 1; // add outside peak
    mdl_bc_dat->P = n_peaks;
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
        rna_mlc_bag_itr ritr;
        rna_mlc_bag_itr_first(&ritr, &bam_bc->rna_mlcs);
        for (; rna_mlc_bag_itr_alive(&ritr); rna_mlc_bag_itr_next(&ritr)) {
            rna_mol_t *mol = rna_mlc_bag_itr_val(&ritr);

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
    mdl->k = -1;

    mdl->sub_lp_x = NULL;
    mdl->sub_lp_hs = NULL;
    mdl->sub_lp_h = NULL;
    mdl->sub_best_hs = NULL;
    mdl->sub_pp_h = NULL;

    mdl->full_lp_x = NULL;
    mdl->full_best_k = NULL;

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
    free(m->sub_pp_h);

    free(m->full_lp_x);
    free(m->full_best_k);

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

    if (mdl->hs_ix->n_hs < 1)
        return err_msg(-1, 0, "mdl_alloc_probs: indices must be set with "
                "'mdl_set_hs_ix'");
    if (mdl->k < 1)
        return err_msg(-1, 0, "mdl_alloc_probs: cell type k=%i must be >= 1", mdl->k);

    uint32_t D = mdl->mdl_bc_dat->all_bcs->n;
    uint32_t n_hs = mdl->hs_ix->n_hs;
    mdl->n_hs = n_hs;

    mdl->sub_lp_hs = calloc(D * n_hs, sizeof(f_t));
    mdl->sub_lp_h = calloc(D * 3, sizeof(f_t));
    mdl->sub_lp_x = calloc(D, sizeof(f_t));
    mdl->sub_best_hs = calloc(D, sizeof(int32_t));
    mdl->sub_pp_h = calloc(D * 3, sizeof(f_t));

    mdl->full_lp_x = calloc(D, sizeof(f_t));
    mdl->full_best_k = calloc(D, sizeof(int32_t));

    if (mdl->full_best_k == NULL)
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
        al = mp->alpha_rna[CMI(bc_ix, par_ix->hs_ix, mp->D)];
    } else if (mol_type == ATAC_IX) {
        al = mp->alpha_atac[CMI(bc_ix, par_ix->hs_ix, mp->D)];
    } else {
        return err_msg(-1, 0, "pr_tdm: invalid `mol_type` index '%i'", mol_type);
    }
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
    uint32_t G = mp->G;
    uint32_t rho_nrow = 3 * G;
    uint32_t rho_ncol = 1 + mp->K;
    if (k >= rho_ncol)
        return err_msg(-1, 0, "p_rna: invalid k=%u", k);

    // probability terms
    f_t pgbm = 1;

    // probability of variant
    if (p_var(mp, mlcl, s_ix, &pgbm) < 0)
        return -1;
    assert(pgbm > 0 && pgbm <= 1);

    size_t n_feat = mv_size(&mlcl->feat_ixs);
    if (n_feat > 1) // skip multi-gene UMIs
        n_feat = 0;

    // Pr(G_dm | T_dm, \rho)
    size_t i;
    for (i = 0; i < n_feat; ++i){
        uint32_t row_ix = mv_i(&mlcl->feat_ixs, i);
        if (row_ix >= rho_nrow)
            return err_msg(-1, 0, "p_rna: invalid gene ix=%u", row_ix);
        uint32_t rho_ix = CMI(row_ix, k, rho_nrow);
        f_t rp = mp->rho[rho_ix];
        if (prob_invalid(rp))
            return err_msg(-1, 0, "p_rna: invalid prob=%f", rp);
        pgbm *= rp;
    }

    if (prob_invalid(pgbm))
        return err_msg(-1, 0, "p_rna: prob=%f is invalid", pgbm);

    *prob = pgbm;
    return(0);
}

int p_atac(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int s_ix, uint16_t k, f_t *prob){
    // prob. value
    f_t ppbm = 1;

    // probability of variant
    if (p_var(mp, mlcl, s_ix, &ppbm) < 0)
        return -1;

    // peak probability
    uint32_t sig_nrow = mp->P;
    uint32_t row_ix = mv_i(&mlcl->feat_ixs, 0);
    if (row_ix >= sig_nrow)
        return err_msg(-1, 0, "p_atac: invalid peak ix=%u", row_ix);
    uint32_t sig_ix = CMI(row_ix, k, sig_nrow);
    f_t sp = mp->sigma[sig_ix];
    if (prob_invalid(sp))
        return err_msg(-1, 0, "p_atac: invalid prob=%f", sp);
    ppbm *= sp;

    if (prob_invalid(ppbm))
        return err_msg(-1, 0, "p_atac: prob=%f is invalid", ppbm);

    *prob = ppbm;
    return 0;
}

int p_bd(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int mol_type, int bc_ix, par_ix_t *par_ix,
        f_t *p_b, f_t *psum) {
    *psum = 0;

    int t_im, ret;
    for (t_im = 0; t_im < par_ix->t_n; ++t_im){
        p_b[t_im] = 1;
        int s_ix = par_ix->t_ix[t_im];
        assert(s_ix >= 0);

        f_t p_t = -1;
        ret = pr_tdm(mp, mol_type, bc_ix, par_ix, t_im, &p_t);
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

        // Pr(T_dm)
        f_t p_t;
        ret = pr_tdm(mp, rna, bc_ix, par_ix, t_im, &p_t);
        if (ret < 0)
            return -1;
        
        // Pr(G_dm, B_dm)
        f_t p_fv = -1;
        if (rna == RNA_IX) {
            ret = p_rna(mp, mlcl, s_ix, k_ix, &p_fv);
        } else if (rna == ATAC_IX) {
            ret = p_atac(mp, mlcl, s_ix, k_ix, &p_fv);
        } else {
            return err_msg(-1, 0, "p_fv: invalid index %i", rna);
        }
        if (ret < 0)
            return -1;

        // Pr(T_dm, G_dm, B_dm)
        p_tgb[t_im] = p_t * p_fv;
        if (prob_invalid(p_tgb[t_im]))
            return err_msg(-1, 0, "p_f_v: Pr(t, f, v)=%f is invalid", p_tgb[t_im]);

        *psum += p_tgb[t_im];
    }
    if (prob_invalid(*psum))
        return err_msg(-1, 0, "p_f_v: Pr(f, v)=%f is invalid", psum);
    return 0;
}

int p_kf_v(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int rna, int bc_ix, par_ix_t *par_ix,
        f_t **p_tgb, uint32_t p_len, f_t *psum) {
    *psum = 0.0;

    uint32_t K = mp->K;
    uint32_t i;
    for (i = 0; i < p_len; ++i)
        (*p_tgb)[i] = 0.0;

    int ret;
    uint32_t k_ix;
    int t_im;
    for (t_im = 0; t_im < par_ix->t_n; ++t_im){
        int s_ix = par_ix->t_ix[t_im];
        assert(s_ix >= 0);

        uint32_t tk_n = t_im == 0 ? 1 : K;
        for (k_ix = 0; k_ix < tk_n; ++k_ix) {
            uint32_t ct_col = t_im == 0 ? 0 : k_ix + 1;
            // Pr(T_dm)
            f_t p_t;
            ret = pr_tdm(mp, rna, bc_ix, par_ix, t_im, &p_t);
            if (ret < 0)
                return -1;

            if (t_im > 0)
                p_t *= mp->kappa[CMI(k_ix, bc_ix, K)];

            // Pr(G_dm, B_dm)
            f_t p_fv = -1;
            if (rna == RNA_IX) {
                ret = p_rna(mp, mlcl, s_ix, ct_col, &p_fv);
            } else if (rna == ATAC_IX) {
                ret = p_atac(mp, mlcl, s_ix, ct_col, &p_fv);
            } else {
                return err_msg(-1, 0, "p_fv: invalid index %i", rna);
            }
            if (ret < 0)
                return -1;

            // Pr(T_dm, G_dm, B_dm)
            if (ct_col >= p_len)
                return err_msg(-1, 0, "p_kf_v: p ix=%u > p_len=%u", ct_col, p_len);

            f_t ptfv = p_t * p_fv;
            if (prob_invalid(ptfv))
                return err_msg(-1, 0, "p_f_v: Pr(tk, f, v)=%f is invalid", ptfv);

            (*p_tgb)[ct_col] += ptfv;
            *psum += ptfv;
        }
    }
    if (num_invalid(*psum)) {
        fprintf(stderr, "ix=%i\n", par_ix->hs_ix);
        return err_msg(-1, 0, "sum p_f_v: Pr(f, v)=%e is invalid", psum);
    }
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
        if (prob_invalid(p_hd) || prob_invalid(p_sd))
            return err_msg(-1, 0, "mdl_sub_e: failed to get pr_hd or pr_sd");
        lp_hs_v[hs_ix] = log(p_hd) + log(p_sd);
    }

    // store P(H_d, S_d, X_d) during loop
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
        uint32_t bc_n_ix = afl ? 1 : n_hs;

        if (afl) ++fixed;
        else ++unfixed;

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);
        kbtree_t(kb_mdl_mlcl) *mols = bc_dat.rna;
        kbtree_t(kb_mdl_mlcl) *frags = bc_dat.atac;

        kbitr_t itr;

        // initialize log Pr(H_d, S_d, X_d)
        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix) {
            if (hs_ix < bc_n_ix)
                lp_htd[hs_ix] = bc_dat.n_bc * lp_hs_v[hs_ix];
            else
                lp_htd[hs_ix] = -INFINITY;
        }

        // only test molecules overlapping variants
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
        // add log Pr(H_d, S_d, X_d) to mdl
        memcpy(&mdl->sub_lp_hs[CMI(0, bc_ix, n_hs)], lp_htd, n_hs * sizeof(f_t));

        // set log Pr(H_d, X_d) to mdl
        mdl->sub_lp_h[CMI(0, bc_ix, 3)] = -INFINITY;
        mdl->sub_lp_h[CMI(1, bc_ix, 3)] = -INFINITY;
        mdl->sub_lp_h[CMI(2, bc_ix, 3)] = -INFINITY;
        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
            if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                return -1;
            mdl->sub_lp_h[CMI(par_ix.hd, bc_ix, 3)] = 
                logsum2expd(mdl->sub_lp_h[CMI(par_ix.hd, bc_ix, 3)],
                        lp_htd[hs_ix]);
        }

        // Pr(X_d | \Theta)
        mdl->sub_lp_x[bc_ix] = logsumexpd(lp_htd, bc_n_ix, &lsret);
        if (lsret < 0)
            return err_msg(-1, 0, "mdl_sub_e: could not logsumexp");
        if (num_invalid(mdl->sub_lp_x[bc_ix]))
            return err_msg(-1, 0, "mdl_sub_e: sub_lp_x[%i]=%.6e is invalid\n",
                    bc_ix, mdl->sub_lp_x[bc_ix]);
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

    uint32_t K = mdl->k;
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
        if (prob_invalid(p_hd) || prob_invalid(p_sd))
            return err_msg(-1, 0, "mdl_sub_e: failed to get pr_hd or pr_sd");
        lp_hs_v[hs_ix] = log(p_hd) + log(p_sd);
    }

    uint32_t p_len = K + 1;
    f_t *p_tgpb = calloc(p_len, sizeof(f_t));

    uint32_t i;
    for (i = 0; i < ix_len; ++i){
        // get barcode and bam data
        int bc_ix = ixs[i];

        // if no barcode
        int efl = bflg_get(bd->absent_bc, bc_ix);
        if (efl == 1) continue;

        // if fixed, set Pr(K_d = 0) = 1
        int afl = bflg_get(bd->amb_flag, bc_ix);

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);
        kbtree_t(kb_mdl_mlcl) *mols = bc_dat.rna;
        kbtree_t(kb_mdl_mlcl) *frags = bc_dat.atac;

        // best (H,S) from sub model
        int32_t sub_best_hs = mdl->sub_best_hs[bc_ix];
        hs_ix_get_pars(mdl->hs_ix, sub_best_hs, &par_ix);

        kbitr_t itr;

        // remove
        // initialize log Pr(H_d, S_d, K_d)
        f_t lp_kd;
        lp_kd = bc_dat.n_bc * lp_hs_v[sub_best_hs];

        // loop over RNA
        kb_itr_first(kb_mdl_mlcl, mols, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, mols, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            f_t psum = 0;
            // get probs P(T_dm, G_dm, B_dm, X_d)
            int pret = p_kf_v(mdl->mp, mlcl, RNA_IX, bc_ix, &par_ix, &p_tgpb, p_len, &psum);
            if (pret < 0)
                return -1;
            if (psum <= 0)
                fprintf(stdout, "fl=%i (h=%i,s1=%i,s2=%i,k=%i) psum=%f\n",
                        afl, par_ix.hd, par_ix.s1, par_ix.s2, par_ix.k, psum);
            assert(psum > 0 && psum <= 1);
            lp_kd += mlcl->counts * log(psum);
        }

        // loop over ATAC
        kb_itr_first(kb_mdl_mlcl, frags, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, frags, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);
            f_t psum = 0;
            // get probs P(T_dm, G_dm, B_dm, X_d)
            int pret = p_kf_v(mdl->mp, mlcl, ATAC_IX, bc_ix, &par_ix, &p_tgpb, p_len, &psum);
            if (pret < 0)
                return -1;
            if (psum <= 0)
                fprintf(stdout, "fl=%i (h=%i,s1=%i,s2=%i,k=%i) psum=%f\n",
                        afl, par_ix.hd, par_ix.s1, par_ix.s2, par_ix.k, psum);
            assert(psum > 0 && psum <= 1);
            lp_kd += mlcl->counts * log(psum);
        }
        // Pr(X_d | \Theta)
        mdl->full_lp_x[bc_ix] = lp_kd;
        if (lsret < 0)
            return err_msg(-1, 0, "mdl_full_e: could not logsumexp");
        if (num_invalid(mdl->full_lp_x[bc_ix]))
            return err_msg(-1, 0, "mdl_full_e: p_x[%i]=%.6e is invalid\n",
                    bc_ix, mdl->full_lp_x[bc_ix]);
    }
    free(p_tgpb);
    free(lp_hs_v);
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
    uint32_t D = mdl->mp->D;
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

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);

        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
            // get (h,s,t) from index
            if (hs_ix_get_pars(mdl->hs_ix, hs_ix, &par_ix) < 0)
                return(-1);

            if (afl && hs_ix > 0)
                break;

            // for alpha
            f_t a_tot, new_par, pdiff;
            f_t psc = 1e-8;
            f_t ar[2] = {psc, psc}; // 1st index is ambient
            f_t aa[2] = {psc, psc};

            // get P(H,S|X)
            uint32_t mix = CMI(hs_ix, bc_ix, n_hs);
            f_t cp_hsx = exp(mdl->sub_lp_hs[mix] - mdl->sub_lp_x[bc_ix]);
            if (prob_invalid(cp_hsx))
                return err_msg(-1, 0, "mdl_sub_m: "
                        "P(Z|X)=%.6e is an invalid probability", cp_hsx);

            // skip over very small prob values
            // if (cp_hsx < 1e-20)
            //     continue;

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
                    ar[c_ix] += mlcl->counts * pt;
                }
            }

            // Directly maximize RNA alpha
            a_tot = ar[0] + ar[1];
            new_par = a_tot <= 0.0 ? 0.5 : ar[0] / a_tot;
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
                    aa[c_ix] += mlcl->counts * pt;
                }
            }

            // ATAC alpha
            a_tot = aa[0] + aa[1];
            new_par = a_tot <= 0.0 ? 0.5 : aa[0] / a_tot;
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
    uint16_t K = mdl->mp->K;
    uint32_t G = mdl->mp->G;
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;

    par_ix_t par_ix;
    par_ix_init(&par_ix);

    // rho counter
    uint32_t rho_ncol = mdl->mp->K + 1;
    uint32_t rho_nrow = 3 * G;
    uint32_t rho_ne = rho_ncol * rho_nrow;
    f_t *rho_sum = malloc(rho_ne * sizeof(f_t));
    for (i = 0; i < rho_ne; ++i) rho_sum[i] = 0;

    // sigma counter
    uint32_t sig_ncol = mdl->mp->K + 1;
    uint32_t sig_nrow = mdl->mp->P;
    uint32_t sig_ne = sig_ncol * sig_nrow;
    // sig_sum col 0: ambient, col 1: cell, row 0: outside peak, row 1: inside peak
    f_t *sig_sum = malloc(sig_ne * sizeof(f_t)); // outside/inside peak by ambient/nuclear
    for (i = 0; i < sig_ne; ++i)
        sig_sum[i] = 0;

    // kappa counter
    uint32_t kappa_ne = K;
    f_t *kappa_sum = malloc(kappa_ne * sizeof(f_t));
    if (kappa_sum == NULL)
        return err_msg(-1, 0, "mdl_full_m: %s", strerror(errno));

    uint32_t p_len = K + 1;
    f_t *p_tgpb = calloc(p_len, sizeof(f_t));

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

        // best (H,S) from sub model
        int32_t sub_best_hs = mdl->sub_best_hs[bc_ix];
        hs_ix_get_pars(mdl->hs_ix, sub_best_hs, &par_ix);
        uint32_t tk_n = par_ix.hs_ix == 0 ? 1 : p_len;

        // kappa counter
        uint32_t k_ix;
        for (k_ix = 0; k_ix < K; ++k_ix)
            kappa_sum[k_ix] = 1e-8;
        f_t kappa_tot = K * 1e-8;

        // barcode count data
        mdl_mlcl_bc_t bc_dat = mv_i(&bd->bc_mlcl, bc_ix);
        kbtree_t(kb_mdl_mlcl) *mols = bc_dat.rna;
        kbtree_t(kb_mdl_mlcl) *frags = bc_dat.atac;
        kbitr_t itr_rna, itr_atac;

        // Get P(H,S,X)
        // Divide by \sum P(T,G,P,B,X)
        // Multiple by P(T,G,P,B,X)

        // loop over RNA
        kb_itr_first(kb_mdl_mlcl, mols, &itr_rna); 
        for (; kb_itr_valid(&itr_rna); kb_itr_next(kb_mdl_mlcl, mols, &itr_rna)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr_rna);

            // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
            f_t psum = 0;
            int pret = p_kf_v(mdl->mp, mlcl, RNA_IX, bc_ix, &par_ix, &p_tgpb, tk_n, &psum);
            if (pret < 0)
                return -1;
            assert(psum > 0);
            f_t kct_sum = psum - p_tgpb[0]; // sum of kappas for cell type
            if (kct_sum <= 0)
                kct_sum = 1e-8;

            // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
            uint32_t t_im;
            for (t_im = 0; t_im < tk_n; ++t_im) {
                k_ix = t_im - 1;
                uint8_t c_ix = t_im; // 0 is ambient, otherwise cell
                f_t pt = p_tgpb[t_im] / psum;
                if (prob_invalid(pt))
                    return err_msg(-1, 0, "mdl_full_m: invalid prob pt=%f", pt);

                // kappa
                if (t_im > 0) {
                    f_t kappa_add = mlcl->counts * p_tgpb[t_im] / kct_sum;
                    kappa_sum[k_ix] += kappa_add;
                    kappa_tot += kappa_add;
                }

                // rho
                uint32_t n_feats = mv_size(&mlcl->feat_ixs);
                if (n_feats > 1) n_feats = 0;
                uint32_t l;
                for (l = 0; l < n_feats; ++l){
                    uint32_t row_ix = (uint32_t)(mv_i(&mlcl->feat_ixs, l));
                    assert(row_ix < rho_nrow);
                    uint32_t rho_ix = CMI(row_ix, c_ix, rho_nrow);
                    rho_sum[rho_ix] += mlcl->counts * pt;
                }
            }
        }

        // loop over ATAC
        kb_itr_first(kb_mdl_mlcl, frags, &itr_atac); 
        for (; kb_itr_valid(&itr_atac); kb_itr_next(kb_mdl_mlcl, frags, &itr_atac)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr_atac);

            // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
            f_t psum = 0;
            int pret = p_kf_v(mdl->mp, mlcl, ATAC_IX, bc_ix, &par_ix, &p_tgpb, tk_n, &psum);
            if (pret < 0)
                return -1;
            assert(psum > 0);
            f_t kct_sum = psum - p_tgpb[0]; // sum of kappas for cell type
            if (kct_sum <= 0)
                kct_sum = 1e-8;

            // get peak ix
            assert(mv_size(&mlcl->feat_ixs) == 1);
            int32_t pk32 = mv_i(&mlcl->feat_ixs, 0);
            assert(pk32 >= 0);
            uint32_t pk = pk32; // 0: outside peak, 1: peak index

            // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
            uint32_t t_im;
            for (t_im = 0; t_im < tk_n; ++t_im) {
                k_ix = t_im - 1;
                uint8_t c_ix = t_im; // 0 is ambient, otherwise cell
                f_t pt = p_tgpb[t_im] / psum;
                if (prob_invalid(pt))
                    return err_msg(-1, 0, "mdl_full_m: invalid prob pt=%f", pt);
                
                // kappa
                if (t_im > 0) {
                    f_t kappa_add = mlcl->counts * p_tgpb[t_im] / kct_sum;
                    kappa_sum[k_ix] += kappa_add;
                    kappa_tot += kappa_add;
                }

                // sigma
                uint32_t sig_ix = CMI(pk, c_ix, sig_nrow);
                sig_sum[sig_ix] += mlcl->counts * pt;
            }
        }

        // maximize kappa directly
        for (k_ix = 0; k_ix < mdl->mp->K; ++k_ix) {
            uint32_t mix = CMI(k_ix, bc_ix, mdl->mp->K);
            mdl->mp->kappa[mix] = kappa_sum[k_ix] / kappa_tot;
        }
    }

    // add to tmp sum variables
    pthread_mutex_lock(&mdl->mp->sum_lock);

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

    free(p_tgpb);
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
                if (num_invalid(mdl->mp->_rho_sum[rix]))
                    return err_msg(-1, 0, "mdl_m_rho: invalid %f", 
                            mdl->mp->_rho_sum[rix]);
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

    uint32_t P = mdl->mp->P;
    uint32_t sig_nrow = P;
    uint32_t sig_ncol = 1 + mdl->mp->K;
    uint32_t i, j;
    if (mdl->has_atac){
        f_t *sig_tot = calloc(sig_ncol, sizeof(f_t));
        for (i = 0; i < sig_ncol; ++i)
            sig_tot[i] = 0;
        for (i = 0; i < sig_nrow; ++i) {
            for (j = 0; j < sig_ncol; ++j) {
                uint32_t six = CMI(i, j, sig_nrow);
                if (num_invalid(mdl->mp->_sigma_sum[six]))
                    return err_msg(-1, 0, "mdl_m_sigma: invalid %f", 
                            mdl->mp->_sigma_sum[six]);
                sig_tot[j] += mdl->mp->_sigma_sum[six];
            }
        }
        // check for errors
        for (j = 0; j < sig_ncol; ++j) {
            if (num_invalid(sig_tot[j]))
                return err_msg(-1, 0, "mdl_m_sig: sig prob. sum is invalid: [%f]", 
                        sig_tot[j]);
        }
        for (i = 0; i < sig_nrow; ++i){
            for (j = 0; j < sig_ncol; ++j){
                uint32_t six = CMI(i, j, sig_nrow);
                new_par = mdl->mp->_sigma_sum[six] / sig_tot[j];
                mdl->mp->_par_diff += fabs(new_par - mdl->mp->sigma[six]);
                mdl->mp->sigma[six] = new_par;
            }
        }
        free(sig_tot);
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

    // rho
    if (mdl_m_rho(mdl) < 0)
        return -1;

    // sigma
    if (mdl_m_sigma(mdl) < 0)
        return -1;

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
    f_t q_delta;
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

        if (sub)
            fret = mdl_thrd_call(mdl, ixs, D, 1);
        else
            fret = mdl_thrd_call(mdl, ixs, D, 0);
        if (fret < 0)
            return -1;

        if (sub)
            fret = mdl_sub_est(mdl);
        else
            fret = mdl_full_est(mdl);
        if (fret < 0)
            return -1;

        if (sub)
            fret = mdl_sub_llk(mdl, &q);
        else
            fret = mdl_full_llk(mdl, &q);
        if (fret < 0)
            return -1;

        if (mdl_delta_q(q_prev, q, &q_delta) < 0)
            return -1;
        if (objs->verbose)
            log_msg("iteration %u: Q=%.6e, delta=%.6e", i+1, q, q_delta);
        if (q_delta < 0)
            err_msg(0, 1, "llk decreased: q1=%.6e q2=%.6e", q_prev, q);
        q_prev = q;

        if (q_delta < mdl->eps) break;
    }
    if (objs->verbose) log_msg("finished EM");
    if (mdl_pars_check(mdl->mp) < 0)
        return -1;

    free(ixs);

    return(0);
}

int mdl_full_llk(mdl_t *mdl, f_t *llk) {
    if (mdl == NULL || llk == NULL)
        return err_msg(-1, 0, "mdl_full_llk: argument is null");

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    if (bd == NULL)
        return err_msg(-1, 0, "mdl_full_llk: barcode data is not initialized");

    uint32_t i, n_bcs = (uint32_t)bd->all_bcs->n;
    f_t llkt = 0;
    for (i = 0; i < n_bcs; ++i){
        if (bflg_get(bd->absent_bc, i))
            continue;
        llkt += mdl->full_lp_x[i];
        assert(!num_invalid(llkt));
    }
    *llk = llkt;
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
        if (bflg_get(bd->absent_bc, i))
            continue;
        llkt += mdl->sub_lp_x[i];
        assert(!num_invalid(llkt));
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

    uint32_t n_hs = mdl->n_hs;
    uint32_t D = mdl->mp->D;
    free(mdl->sub_best_hs);
    mdl->sub_best_hs = malloc(D * sizeof(uint32_t));
    if (mdl->sub_best_hs == NULL)
        return err_msg(-1, 0, "mdl_sub_best_hs: %s", strerror(errno));

    uint32_t bc_i;
    for (bc_i = 0; bc_i < D; ++bc_i){
        if (bflg_get(bd->absent_bc, bc_i))
            continue;
        // get best (H,S) index
        uint32_t b_hs_ix = 0;
        uint32_t t_hs_ix = 0;
        for (t_hs_ix = 1; t_hs_ix < n_hs; ++t_hs_ix) {
            uint32_t b_mix = CMI(b_hs_ix, bc_i, n_hs);
            uint32_t t_mix = CMI(t_hs_ix, bc_i, n_hs);
            if (mdl->sub_lp_hs[t_mix] > mdl->sub_lp_hs[b_mix])
                b_hs_ix = t_hs_ix;
        }
        mdl->sub_best_hs[bc_i] = b_hs_ix;
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
        if (bflg_get(bd->absent_bc, bc_i))
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

int mdl_full_best_k(mdl_t *mdl){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_full_best_k: argument is null");

    if (mdl->mp == NULL || mdl->mp->kappa == NULL || mdl->mp->D == 0)
        return err_msg(-1, 0, "mdl_full_best_k: model must be initialized");

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    if (bd == NULL)
        return err_msg(-1, 0, "mdl_full_best_k: barcode data is not initialized");

    int32_t K = mdl->k;
    uint32_t D = mdl->mp->D;
    free(mdl->full_best_k);
    mdl->full_best_k = malloc(D * sizeof(int32_t));
    if (mdl->full_best_k == NULL)
        return err_msg(-1, 0, "mdl_full_best_k: %s", strerror(errno));

    uint32_t bc_i;
    for (bc_i = 0; bc_i < D; ++bc_i){
        if (bflg_get(bd->absent_bc, bc_i))
            continue;
        if (mdl->sub_best_hs[bc_i] == 0) {
            mdl->full_best_k[bc_i] = 0;
            continue;
        }

        // get best (K) index
        int32_t b_k_ix = 0;
        int32_t t_k_ix = 0;
        for (t_k_ix = 1; t_k_ix < K; ++t_k_ix) {
            int32_t b_mix = CMI(b_k_ix, bc_i, K);
            int32_t t_mix = CMI(t_k_ix, bc_i, K);
            if (mdl->mp->kappa[t_mix] > mdl->mp->kappa[b_mix])
                b_k_ix = t_k_ix + 1;
        }
        mdl->full_best_k[bc_i] = b_k_ix;
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
    mdl->max_iter = objs->max_iter;
    mdl->threads = objs->threads;
    mdl->alpha_vars = objs->alpha_vars > 0;
    mdl->alpha_max = objs->alpha_max;

    if (mdl_set_rna_atac(mdl, bam_dat->has_rna, bam_dat->has_atac) < 0)
        return -1;

    if (objs->verbose)
        log_msg("initializing data", mv_size(&mdl->mdl_bc_dat->bc_mlcl));
    if (mdl_bc_dat_bam_data(mdl->mdl_bc_dat, bam_dat, objs) < 0)
        return -1;
    if (objs->verbose)
        log_msg("initialized %zu barcodes", mv_size(&mdl->mdl_bc_dat->bc_mlcl));

    if (mdl_set_samples(mdl, objs->vcf_hdr) < 0)
        return -1;

    if (mdl_set_k(mdl, objs->k) < 0)
        return -1;

    if (mdl_set_hs_ix(mdl) < 0)
        return -1;

    if (objs->verbose)
        log_msg("allocating arrays");
    if (mdl_alloc_probs(mdl) < 0)
        return -1;

    // initialize parameters
    if (objs->verbose)
        log_msg("initializing parameters");
    if (mdl_pars_set_dat(mdl->mp, mdl->mdl_bc_dat, objs,
                mdl->samples->n, mdl->k) < 0)
        return -1;
    if (objs->verbose)
        log_msg("running EM");

    // run sub model
    if (mdl_em(mdl, objs, 1) < 0)
        return -1;

    // get best (H,S) index
    if (mdl_sub_best_hs(mdl) < 0)
        return -1;
    // get posterior prob for (H)
    if (mdl_sub_pp_h(mdl) < 0)
        return -1;

    // run full model
    if (mdl_em(mdl, objs, 0) < 0)
        return -1;

    // get best (K) index
    if (mdl_full_best_k(mdl) < 0)
        return -1;

    if (objs->verbose)
        log_msg("getting barcode likelihoods");

    // save output
    if (objs->verbose)
        log_msg("writing out likelihoods");
    if (write_sub_llk(mdl, objs->out_fn) < 0)
        return -1;
    // if (write_full_llk(mdl, objs->out_fn) < 0)
    //     return -1;
    if (write_samples(mdl, objs->out_fn) < 0)
        return -1;
    if (objs->verbose)
        log_msg("writing out parameters");
    if (write_lambda(mdl, objs->out_fn) < 0)
        return -1;
    if (write_kappa(mdl, objs->out_fn) < 0)
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

int write_kappa(mdl_t *mdl, char *fn){
    if (mdl == NULL)
        return err_msg(-1, 0, "write_kappa: arguments are NULL");

    if (mdl->mdl_bc_dat == NULL)
        return err_msg(-1, 0, "write_kappa: no barcodes found");
    if (mdl->mp == NULL || mdl->mp->kappa == NULL)
        return err_msg(-1, 0, "write_kappa: model hasn't been initialized");
    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    int n_test_bc = mdl->mp->T;
    uint32_t K = mdl->mp->K;

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    char *kappa_fn = ".kappa.txt.gz";
    char *out_kappa_fn = strcat2((const char*)fn, (const char*)kappa_fn);

    // row names
    char **bc_row_names = str_map_ca(bd->test_bcs);
    assert(n_test_bc == bd->test_bcs->n);

    // kappa array for test barcodes
    f_t *ka = malloc(n_test_bc * K * sizeof(f_t));
    int i;
    uint32_t j;
    for (i = 0; i < n_test_bc; ++i){
        char *bc = str_map_str(bd->test_bcs, i);
        int bc_ix = str_map_ix(bd->all_bcs, bc);
        for (j = 0; j < K; ++j)
            ka[CMI(i, j, n_test_bc)] = mdl->mp->kappa[CMI(j, bc_ix, K)];
    }

    // write matrix
    ret = write_matrix_double(out_kappa_fn, ka, NULL, NULL, NULL, 
            bc_row_names, n_test_bc, NULL, K, 
            delim, nl, decp);
    free(ka);
    if (ret < 0)
        return err_msg(-1, 0, "write_kappa: failed to write matrix to file");

    for (i = 0; i < n_test_bc; ++i)
        free(bc_row_names[i]);
    free(bc_row_names);
    free(out_kappa_fn);

    return 0;
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
    /*
    char **row_names = malloc(2 * sizeof(char *));
    row_names[0] = strdup("Ambient");
    if (row_names[0] == NULL)
        return err_msg(-1, 0, "write_sigma: %s", strerror(errno));
    row_names[1] = strdup("Cell");
    if (row_names[1] == NULL)
        return err_msg(-1, 0, "write_sigma: %s", strerror(errno));
    */

    // col names
    uint32_t sig_nrow = mdl->mp->P;
    uint32_t sig_ncol = mdl->mp->K + 1;
    char **col_names = malloc(sig_ncol * sizeof(char*));
    col_names[0] = strdup("Ambient");
    if (col_names[0] == NULL)
        return err_msg(-1, 0, "write_sig: %s", strerror(errno));
    size_t buf_size = 100;
    uint32_t i;
    for (i = 1; i < sig_ncol; ++i) {
        col_names[i] = malloc(buf_size * sizeof(char));
        int2strp(i, col_names + i, &buf_size);
    }

    // write matrix
    ret = write_matrix_double(out_sigma_fn, mdl->mp->sigma, NULL, NULL, NULL, 
            NULL, sig_nrow, col_names, sig_ncol, 
            delim, nl, decp);
    if (ret < 0)
        return err_msg(-1, 0, "write_sigma: failed to write matrix to file");

    for (i = 0; i < sig_ncol; ++i)
        free(col_names[i]);
    free(col_names);
    free(out_sigma_fn);

    return(0);
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
        assert(all_bc_ix >= 0);
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

    if (mdl == NULL || bam_dat == NULL)
        return err_msg(-1, 0, "write_res: argument is null");

    mdl_bc_dat_t *bd = mdl->mdl_bc_dat;
    uint32_t n_test_bc = mdl->mp->T;
    uint32_t n_all_bc = mdl->mp->D;
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
    uint32_t bc_i, i;
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

        // get best (H,S) assignment
        uint32_t best_hs = mdl->sub_best_hs[all_bc_ix];
        if (hs_ix_get_pars(mdl->hs_ix, best_hs, &par_ix) < 0)
            return -1;

        // write best H_d string
        fret = fputc(delim, fp);
        if (par_ix.hd == 0) {
            fret = fputs("Empty", fp);
        } else if (par_ix.hd == 1) {
            fret = fputs("Singlet", fp);
        } else {
            fret = fputs("Doublet", fp);
        }

        // get best S_d strings
        fret = fputc(delim, fp);
        if (best_hs == 0) {
            fputs("Empty", fp);
        } else {
            fputs(str_map_str(samples, par_ix.s1), fp);
            if (par_ix.s2 >= 0) {
                fputs(":", fp);
                fputs(str_map_str(samples, par_ix.s2), fp);
            }
        }

        // get best K_d
        int32_t best_k = mdl->full_best_k[all_bc_ix];
        fret = fputc(delim, fp);
        if (int2strp(best_k, &pstr, &buf_size) < 0)
            return -1;
        fputs(pstr, fp);

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
            par = mdl->sub_pp_h[CMI(i, all_bc_ix, 3)];
            if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
                return err_msg(-1, 0, "write_res: failed to convert %f to string", par);
            fret = fputc(delim, fp);
            fret = fputs(pstr, fp);
        }

        // write llks of H_d
        for (i = 0; i < 3; ++i){
            par = mdl->sub_lp_h[CMI(i, all_bc_ix, 3)];
            if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
                return err_msg(-1, 0, "write_res: failed to convert %f to string", par);
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

    uint32_t p_len = mdl->mp->K + 1;
    f_t *p_tgpb = calloc(p_len, sizeof(f_t));

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

        // get best (H, S) and (K)
        uint32_t best_hs = mdl->sub_best_hs[d_bc_i];
        if (hs_ix_get_pars(mdl->hs_ix, best_hs, &par_ix) < 0)
            return -1;

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

            // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
            f_t psum = 0;
            int pret = p_kf_v(mdl->mp, mlcl, RNA_IX, d_bc_i, &par_ix, &p_tgpb, p_len, &psum);
            if (pret < 0)
                return -1;
            assert(psum > 0);

            // hold the amb posterior prob.
            f_t amb_pp = p_tgpb[0] / psum;

            assert(!prob_invalid(amb_pp));
            f_t all_reads = (f_t)mlcl->counts;
            f_t nuc_reads = all_reads * (1.0 - amb_pp);

            uint32_t mv_ix;
            // add feature read count
            for (mv_ix = 0; mv_ix < mv_size(&mlcl->feat_ixs); ++mv_ix) {
                int32_t fs_ix, s_ix, f_ix;
                fs_ix = mv_i(&mlcl->feat_ixs, mv_ix);
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

        // loop over ATAC
        kbtree_t(kb_mdl_mlcl) *frags = bc_dat.atac;
        kb_itr_first(kb_mdl_mlcl, frags, &itr); 
        for (; kb_itr_valid(&itr); kb_itr_next(kb_mdl_mlcl, frags, &itr)){
            mdl_mlcl_t *mlcl = &kb_itr_key(mdl_mlcl_t, &itr);

            // skip UMIs with no gene or variant
            if (mv_size(&mlcl->feat_ixs) == 0 && mv_size(&mlcl->var_ixs) == 0)
                continue;

            // skip multi-peak frags
            if (mv_size(&mlcl->feat_ixs) > 1)
                continue;

            // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
            f_t psum = 0;
            int pret = p_kf_v(mdl->mp, mlcl, ATAC_IX, d_bc_i, &par_ix, &p_tgpb, p_len, &psum);
            if (pret < 0)
                return -1;
            assert(psum > 0);

            // hold the amb posterior prob.
            f_t amb_pp = p_tgpb[0] / psum;

            assert(!prob_invalid(amb_pp));
            f_t all_reads = (f_t)mlcl->counts;
            f_t nuc_reads = all_reads * (1.0 - amb_pp);

            uint32_t mv_ix;
            // add feature read count
            for (mv_ix = 0; mv_ix < mv_size(&mlcl->feat_ixs); ++mv_ix) {
                int32_t p_ix;
                p_ix = mv_i(&mlcl->feat_ixs, 0);
                if (p_ix == 0)
                    continue;
                --p_ix;
                assert(p_ix >= 0);

                cnt_node_t *p, t;
                memset(&t, 0, sizeof(cnt_node_t));
                t.ix = (int)p_ix;
                t.counts[0] = nuc_reads;
                p = kb_getp(kh_cnode, bcc->atac_pc, &t);
                if (!p){
                    kb_putp(kh_cnode, bcc->atac_pc, &t);
                    ++bam_counts->atac_pcs_nz;
                } else {
                    p->counts[0] += nuc_reads;
                }
                bam_counts->has_atac_pc = 1;
            }

            // add variant read count
            for (mv_ix = 0; mv_ix < mv_size(&mlcl->var_ixs); ++mv_ix) {
                int32_t va_ix, v_ix, a_ix;
                va_ix = mv_i(&mlcl->var_ixs, mv_ix);
                v_ix = va_ix % n_vars;
                a_ix = va_ix / n_vars;
                assert(a_ix < MAX_ALLELE);

                cnt_node_t *p, t;
                memset(&t, 0, sizeof(cnt_node_t));
                t.ix = (int)v_ix;
                t.counts[a_ix] = nuc_reads;
                p = kb_getp(kh_cnode, bcc->atac_ac, &t);
                if (!p){
                    kb_putp(kh_cnode, bcc->atac_ac, &t);
                    ++bam_counts->atac_acs_nz;
                } else {
                    p->counts[a_ix] += nuc_reads;
                }
                bam_counts->has_atac_ac = 1;
            }
        }

        kh_val(bam_counts->bc_counts, k_bt) = bcc;
    }

    if (bam_counts_write(bam_counts, gene_anno, gv, out_fn, 1) < 0)
        return err_msg(-1, 0, "failed write corrected counts");

    bam_counts_dstry(bam_counts);
    free(p_tgpb);
    return 0;
}
