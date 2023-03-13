
#ifndef MOD_H
#define MOD_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include "htslib/khash.h"
#include "gtf_anno.h"
#include "array_util.h"
#include "str_util.h"
#include "variants.h"
#include "clopts.h"
#include "rna_data.h"
#include "atac_data.h"
#include "bam_dat.h"
#include "bc_stats.h"
#include "bits.h"

#define f_t double

// return a random number in [0,1]
#define frand() ((f_t) rand() / (RAND_MAX))
// return a random number in [0,1)
#define frandlt1() ((f_t) rand() / (RAND_MAX+1.0))
// return a random number in (0,1)
#define frandin() ((f_t) (rand()+1.0) / (RAND_MAX+2.0))
// return a random number in [0.1, 0.9]
#define frandin2() (0.1 + (0.8 * (f_t)rand() / (RAND_MAX)))

#define TAU 0.01

/*******************************************************************************
 * hold bam barcode counts
 ******************************************************************************/

typedef struct mdl_mlcl_t {
    mv_t(i32) feat_ixs; // vector of feature indices (length 1 for now)
                      // For RNA, it is the gene index.
                      // For ATAC, it is peak (1) or not (0)
    mv_t(i32) var_ixs; // vector of variant indices
    uint32_t counts; // number of molecules with this feature-variant comb.
} mdl_mlcl_t;

static inline int mdl_mlcl_cmp(mdl_mlcl_t m1, mdl_mlcl_t m2){
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

    return(0);
}

// btree to hold molecules per barcode
KBTREE_INIT(kb_mdl_mlcl, mdl_mlcl_t, mdl_mlcl_cmp);

typedef struct mdl_bc_t {
    kbtree_t(kb_mdl_mlcl) *rna;
    kbtree_t(kb_mdl_mlcl) *atac;
} mdl_bc_t;

/*********************
 * parameters
 *********************/

/*! @typedef mdl_pars_t
 * All matrices are column major indexed.
 */
typedef struct {
    // Parameters
    uint32_t D; // number of droplets
    uint32_t G; // number of genes
    uint32_t V; // number of variants
    uint16_t M; // number of samples

    f_t lambda[3]; // empty, singlet, doublet prop (3 x 1 array)
    f_t *pi; // sample prop (M x 1 array)
    f_t pi_d_sum;
    f_t *alpha; // droplet contamination prob. (D x 1 array)
    f_t *rho; // CM expression probs (3G x 2 array) col 0 is ambient col 1 is cell
    f_t sigma[2]; // prob. in peak (2 x 1 array), 0: ambient, 1: cell
    f_t *gamma; // CM fixed genotypes prob. {0,1} (M+1 x V array).
    f_t tau; // probability of a sequencing error (fixed at 0.01)

    // temporary counts for the maximization step (log values)
    pthread_mutex_t sum_lock; // lock for the sum variables
    f_t _lambda_sum[3];
    f_t *_pi_sum; // sample prop (M x 1 array)
    f_t _pi_d_sum;
    f_t *_alpha_sum; // CM droplet contamination prob. (D x 2 array) rows droplets, col 0: ambient col 1: cell
    f_t *_rho_sum; // CM expression probs (3G x 2 array) col 0 is ambient col 1 is cell
    f_t _sigma_sum[4]; // CM open chromatin peak (2 x 2 array), col: ambient,cell; row: outside,inside peak

    f_t _par_diff;
} mdl_pars_t;

/*! @typedef mdl_fit_t
 * The indexes are the hs_ix indices
 */
typedef struct {
    f_t *bc_llks; // BC x ix array. The indices follow hs_ix

    // All vectors are length BC
    f_t *best_sng_llk;
    f_t *sec_sng_llk;
    f_t *best_dbl_llk;
    uint32_t *best_sng_ix;
    uint32_t *sec_sng_ix;
    uint32_t *best_dbl_ix;

    f_t *pp; // BC x 3 array; posterior probs for empty, singlet, doublet
} mdl_fit_t;

/*! @typedef c_probs_t
 * These are log probabilities
 * p(x) is P(X | D) = \sum_Z P(X, Z | D).
 * cp_hs is length 1 + M + M(M-1)/2
 */
typedef struct {
    f_t p_x; // P(X_d | Theta) data probability
    f_t *lp_hs; // log P(H_d, S_d | Theta) length: _nrow_hs
    f_t *cp_hs; // P(H_d, S_d | X_d, Theta) length: _nrow_hs
} c_probs_t;

/*! @typedef mdl_t
 */
typedef struct {
    str_map *all_bcs; // whitelist of barcodes (may not be in same order as bam_data)
    str_map *test_bcs; // barcodes to calculate llk (those that are unfixed in fix_flag)
    str_map *samples;

    // order of barcodes in flags is the order of barcodes in all_bcs
    // all length of all_bcs
    bflg_t *fix_flag; // 0: unfixed, 1: fixed (ambient)
    bflg_t *absent_bc; // 0: unfixed, 1: fixed (absent)

    uint32_t c_probs_len;
    c_probs_t *c_probs; // length of all_bcs

    // H values are 0, 1, 2
    // S values are 2 values specifying sources in droplet (0, ..., M-1, M)
    // source values are 0, ..., M-1 for samples, and M for ambient. -1 for invalid
    uint32_t _nrow_hs, _nrow_hst;
    int *hs_ix; // (1 + M + M(M-1)/2) by 3 array, columns are H, S (sample 1), S (sample 2)

    f_t eps;
    uint16_t max_iter;

    uint8_t has_rna, has_atac;

    uint16_t threads;

    mdl_pars_t *mp;
    mdl_fit_t *mf;
} mdl_t;

// initialization functions return NULL on error
/* initialize/destroy mdl_pars struct */
mdl_pars_t *init_mdl_pars();
void destroy_mdl_pars(mdl_pars_t *gp);

/* initialize/destroy mdl_fit_t */
mdl_fit_t *init_mdl_fit();
void destroy_mdl_fit(mdl_fit_t *mf);

/* initialize/destroy mdl_t struct */
mdl_t *mdl_alloc();
void mdl_dstry(mdl_t *m);

/* Set filtered and all barcodes in the mdl struct
 *
 * @return 0 on success, -1 on error.
 */
int mdl_set_bcs(mdl_t *mdl, bam_data_t *bam_dat, obj_pars *objs);
int mdl_set_samples(mdl_t *mdl, obj_pars *objs);

/* Initialize data members and parameters */
int mdl_init_all(mdl_t *mdl, bam_data_t *bam_data, obj_pars *objs);

int mdl_get_hst(mdl_t *mdl, int hs_ix, int *hd, int *s1, int *s2, int t_ix[3], int *t_n);

/*******************************************************************************
 * Probability functions
 ******************************************************************************/

// get rho
static inline f_t pr_rho_gene(f_t *rho, seq_gene_t seq_gene, uint32_t col, uint32_t G, uint32_t rho_nrow){
    uint32_t g_ix = (uint32_t)seq_gene.gene_id;
    uint32_t spl = (uint32_t)seq_gene.splice;
    uint32_t f_ix = g_ix + (spl * G);
    f_t p_gdm = rho[CMI(f_ix, col, rho_nrow)];
    return(p_gdm);
}

// get gamma
static inline f_t pr_gamma_var(f_t *gamma, seq_vac_t vac, uint32_t s_ix, uint32_t nrow, f_t tau){
    if (vac.allele > 2) return 1.0; // other alleles are considered missing
    uint32_t vix = (uint32_t)vac.vix;
    f_t ap = gamma[ CMI(s_ix, vix, nrow) ];
    if (ap < 0) return 1.0; // if allele is missing
    if (vac.allele == 0) ap = 1.0 - ap;
    f_t p_be0 = (1.0 - tau) * ap;
    f_t p_be1 = tau * 0.25;
    f_t p_bdm = p_be0 + p_be1;
    return(p_bdm);
}

// Pr(H_d)
void pr_hd(mdl_t *mdl, int hd, f_t *prob);
// Pr(S_d)
void pr_sd(mdl_t *mdl, int hd, int s1, int s2, f_t *prob);
// Pr(T_d)
int pr_tdm(mdl_t *mdl, int bc_ix, int hd, int s_ix, f_t *prob);
// Pr(G_dm, B_dm)
int p_rna(mdl_t *mdl, rna_mol_t *mol, int s_ix, f_t *prob);
// Pr(P_dm, B_dm)
int p_atac(mdl_t *mdl, atac_frag_t *frag, int s_ix, f_t *prob);

/*******************************************************************************
 * Expectation
 ******************************************************************************/

int mdl_e_hs(mdl_t *mdl, bam_data_t *bam_data, int *ixs, uint32_t ix_len);

/*******************************************************************************
 * Maximization
 ******************************************************************************/

/* Maximize parameters for model.
 *
 * @param ixs Barcode indices to run maximization over. Currently just use all
 *  barcodes.
 * @param ix_len Length of ixs vector.
 * @param prior Prior for each element of the maximization.
 */
int mdl_m_init(mdl_t *mdl, f_t prior);
int mdl_m_sum(mdl_t *mdl, bam_data_t *bam_data, int *ixs, uint32_t ix_len);
int mdl_m_pi_fix(mdl_t *mdl);

/* Add gamma parameters
 *
 * The gamma parameters are the genotypes (probabilities of alternate allele).
 * Adds gamma (M+1) rows and V columns stored in an array in column major 
 * format.
 * 
 * @return 0 if successful, 1 on error.
 * 
 */
int mdl_pars_add_gamma(mdl_pars_t *gp, float **a, int nv, int ns);

int mdl_m_set_gamma_amb(mdl_t *mdl);

int mdl_m_est(mdl_t *mdl);
int mdl_delta(mdl_t *mdl, f_t *delta);

int mdl_init_par_dat(mdl_t *mdl, bam_data_t *bam_data);
int mdl_init_par_uni(mdl_t *mdl);
int mdl_init_par_rand(mdl_t *mdl, unsigned int seed);

void *mdl_thrd_fx(void *arg);
int mdl_thrd_call(mdl_t *mdl, bam_data_t *bam_data, int *ixs, uint32_t ix_len);
int mdl_em(mdl_t *mdl, bam_data_t *bam_data, obj_pars *objs);

int mdl_data_llk(mdl_t *mdl, f_t *llk);

int mdl_check_pars(mdl_t *mdl);

int mdl_get_llk(mdl_t *mdl);

/* main function to fit model.
 *
 * Call this function to set parameters and fit the model.
 *
 */
int mdl_fit(bam_data_t *bam_dat, obj_pars *objs);
int fit_mdl_pars(bam_data_t *bam_dat, obj_pars *objs);

// output functions
int write_lambda(mdl_t *mdl, char *fn);
int write_alpha(mdl_t *mdl, char *fn);
int write_rho(mdl_t *mdl, obj_pars *objs, char *fn);
int write_sigma(mdl_t *mdl, obj_pars *objs, char *fn);
int write_llk(mdl_t *mdl, char *fn);
int write_samples(mdl_t *mdl, char *fn);
int write_res(mdl_t *mdl, bam_data_t *bam_dat, char *fn);

#endif // MOD_H

