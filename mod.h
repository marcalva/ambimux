
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
#include "r_count.h"
#include "region.h"

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
 * mdl_mlcl_t
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

    // return(-1);
    return(0);
}

// btree to hold molecules per barcode
KBTREE_INIT(kb_mdl_mlcl, mdl_mlcl_t, mdl_mlcl_cmp);

typedef struct {
    kbtree_t(kb_mdl_mlcl) *rna;
    kbtree_t(kb_mdl_mlcl) *atac;
    uint32_t n_bc; // number of barcodes
} mdl_mlcl_bc_t;

mv_declare(mdl_bc_v, mdl_mlcl_bc_t);

void mdl_mlcl_init(mdl_mlcl_t *mlcl);
void mdl_mlcl_free(mdl_mlcl_t *mlcl);

/*******************************************************************************
 * index structs
 ******************************************************************************/

/*@typedef struct par_ix_t
 * Store the indexes for nuc. num, sample 1, sample 2, T_i, N_T
 */
typedef struct {
    int hs_ix; // (H,S) index
    int hd; // num. of nuclei
    int s1; // sample 1 index
    int s2; // sample 2 index
    int k; // cell type index
    int t_ix[3]; // droplet source index
    int t_n; // num. of droplet sources
} par_ix_t;

void par_ix_init(par_ix_t *par_ix);

/* hs_ix_t struct
 * Specify the possible indices for num of nuclei, sample 1, sample 2,
 * and droplet source.
 * H values are 0, 1, 2
 * S values are 2 values specifying sources in droplet (0, ..., M-1, M)
 * source values are 0, ..., M-1 for samples, and M for ambient. -1 for invalid
 *
 * @field hs_ix a (1 + M + M(M-1)/2) by 3 array. Columns correspond to
 *  H, S1 (sample 1), and S2 (sample2)
 */
typedef struct {
    uint16_t n_sam;
    uint32_t n_hs; // Total number of (H,S) indices.
    int *ix2pars; // 3 by (1 + M + M(M-1)/2) by 3 array, rows are H, S (sample 1), S (sample 2)
                  // Array is stored in column major vector
} hs_ix_t;

void hs_ix_init(hs_ix_t *hs_ix);
hs_ix_t *hs_ix_alloc();
void hs_ix_free(hs_ix_t *hs_ix);
void hs_ix_dstry(hs_ix_t *hs_ix);
void hs_ix_set(hs_ix_t *hs_ix, uint16_t n_samples);

/* hs index to parameters
 * Given an hs index, return num nuclei, sample 1, and sample 2
 * t_ix stores the possible sources in the droplet
 * t_n stores the number of possible sources in droplet
 *
 * @return 0 if successfull, -1 on error
 */
int hs_ix_get_pars(hs_ix_t *hs_ix, unsigned ix, par_ix_t *par_ix);

/*******************************************************************************
 * mdl_pars_t
 ******************************************************************************/

/*! @typedef mdl_pars_t
 * All matrices are column major indexed.
 */
typedef struct {
    // Parameters
    uint32_t D; // number of droplets
    uint32_t T; // number of test droplets
    uint32_t G; // number of genes
    uint32_t V; // number of variants
    uint16_t M; // number of samples
    uint16_t K; // number of cell types (not including ambient)
    uint32_t n_hs;

    f_t lambda[3]; // empty, singlet, doublet prop (3 x 1 array)
    f_t *pi; // sample prop (M x 1 array)
    f_t pi_d_sum;
    f_t *kappa; // cluster probability (K x D array) 
    f_t *alpha_rna; // droplet contamination prob. (D x _nrow_hs array)
    f_t *alpha_atac; // droplet contamination prob. (D x _nrow_hs array)
    f_t *rho; // CM expression probs (3G x (K+1) array) col 0 is ambient col, then cell types
    f_t sigma[2]; // prob. in peak (2 x 1 array), 0: ambient, 1: cell
    f_t *gamma; // CM fixed genotypes prob. {0,1} (M+1 x V array).
    f_t tau; // probability of a sequencing error (fixed at 0.01)

    // temporary counts for the maximization step (log values)
    pthread_mutex_t sum_lock; // lock for the sum variables
    f_t _lambda_sum[3];
    f_t *_pi_sum; // sample prop (M x 1 array)
    f_t *_rho_sum; // CM expression probs (3G x 2 array) col 0 is ambient col 1 is cell
    f_t _sigma_sum[4]; // CM open chromatin peak (2 x 2 array), col: ambient,cell; row: outside,inside peak

    f_t _par_diff;
} mdl_pars_t;

/*******************************************************************************
 * mdl_bc_dat_t
 ******************************************************************************/

typedef struct {
    // do not add directly to str_map, use provided functions mdl_bc_add and mdl_set_samples
    str_map *all_bcs; // barcode ID strings of all barcodes in the data. Don't add
    str_map *test_bcs; // barcodes to calculate llk (those that are unfixed in amb_flag)

    // order of barcodes in flags is the order of barcodes in all_bcs
    // all length of all_bcs
    bflg_t *amb_flag; // 0: unfixed, 1: fixed (ambient)
    bflg_t *absent_bc; // 0: unfixed, 1: fixed (absent)

    // store the barcode counts, order corresponds to all_bcs
    uint32_t G; // num of genes
    uint32_t P; // num of peaks
    uint32_t V; // num of variants
    mv_t(mdl_bc_v) bc_mlcl;

} mdl_bc_dat_t;

/*******************************************************************************
 * mdl_t
 ******************************************************************************/

/*! @typedef mdl_t
 */
typedef struct {
    // sample IDs
    str_map *samples;

    // H values are 0, 1, 2
    // S values are 2 values specifying sources in droplet (0, ..., M-1, M)
    // source values are 0, ..., M-1 for samples, and M for ambient. -1 for invalid
    hs_ix_t *hs_ix;

    int n_hs; // Number of (H,S) indices
    f_t *sub_lp_hs; // log Pr(H_d, S_d, X_d | \Theta) n_hs x D (column major)
    f_t *sub_lp_h; // log Pr(H_d, X_d | \Theta) 3 x D (column major)
    f_t *sub_lp_x; // log Pr(X_d | \Theta) D x 1 (column major)
    int32_t *sub_best_hs; // D x 1 array. Best (H,S) index for each droplet d under 
                           // the sub model.
    f_t *sub_pp_h; // 3 x BC array; posterior probabilities `H_d` for each droplet.

    int k; // number of cell types
    f_t *full_lp_x; // full model w cell types log Pr(X_d | \Theta) D x 1 (column major)
                    // calculated as logsumexp during E step
    int32_t *full_best_k; // D x 1 array. Best (K) index for each droplet d under full model

    f_t eps;
    uint16_t max_iter;
    f_t alpha_max;
    uint8_t has_rna, has_atac;
    uint16_t threads;
    uint8_t alpha_vars;

    mdl_bc_dat_t *mdl_bc_dat;
    mdl_pars_t *mp;
} mdl_t;

/*******************************************************************************
 * mdl_mlcl_bc_t
 ******************************************************************************/

int mdl_mlcl_bc_init(mdl_mlcl_bc_t *mdl_bc);
void mdl_mlcl_bc_free(mdl_mlcl_bc_t *mdl_bc);

/*******************************************************************************
 * mdl_pars_t
 ******************************************************************************/

int mdl_pars_init(mdl_pars_t *mp);
mdl_pars_t *mdl_pars_alloc();
void mdl_pars_free(mdl_pars_t *mp);
void mdl_pars_dstry(mdl_pars_t *mp);

/* Set the parameter number fields (D, T, G, V, M).
 * Allocate memory for the parameter fields in @p mp.
 * Allocate memory for the _sum fields in @p mp.
 */
int mld_pars_set_num_alloc(mdl_pars_t *mp, uint32_t D, uint32_t T,
        uint32_t G, uint32_t V, uint16_t M, uint16_t K);

/* Set the _sum variables in the mdl_pars_t object to psc value.
 */
int mdl_pars_reset_sums(mdl_pars_t *mp, f_t psc);

/* Set pi value to uniform across samples/indices
 */
int mdl_pars_pi_fix(mdl_pars_t *mp);

/* Add gamma parameters
 *
 * The gamma parameters are the genotypes (probabilities of alternate allele).
 * Adds gamma (M+1) rows and V columns stored in an array in column major 
 * format.
 * 
 * @return 0 if successful, -1 on error.
 */
int mdl_pars_add_gamma(mdl_pars_t *gp, float **a, int nv, int ns);

/* Set the ambient gamma parameter as a linear function of the 
 * sample proportions from parameter pi
 * @return 0 if successful, -1 on error.
 */
int mdl_pars_set_gamma_amb(mdl_pars_t *mp);

/* Initialize parameters from data
 *
 * Initializes the parameters in @mp from the data in @p bd.
 */
int mdl_pars_set_dat(mdl_pars_t *mp, mdl_bc_dat_t *bd, obj_pars *objs,
        uint16_t n_samples, uint16_t K);

/* Check validity of parameters
 * @return -1 if any invalid, 0 if OK
 */
int mdl_pars_check(mdl_pars_t *mp);

/*******************************************************************************
 * mdl_bc_dat
 ******************************************************************************/

int mdl_bc_dat_init(mdl_bc_dat_t *mdl_bc_dat);
mdl_bc_dat_t *mdl_bc_dat_alloc();
void mdl_bc_dat_free(mdl_bc_dat_t *mdl_bc_dat);
void mdl_bc_dat_dstry(mdl_bc_dat_t *mdl_bc_dat);

int mdl_bc_dat_add_bc(mdl_bc_dat_t *mdl_bc_dat, char *bc_key, int absent, int ambient,
        int dup_ok, int *found);

int mdl_bc_dat_bam_data(mdl_bc_dat_t *mdl_bc_dat, bam_data_t *bam_data, obj_pars *objs);

/*******************************************************************************
 * mdl_t
 ******************************************************************************/

/* initialize/destroy mdl_t struct */
mdl_t *mdl_alloc();
void mdl_dstry(mdl_t *m);

int mdl_set_samples(mdl_t *mdl, bcf_hdr_t *vcf_hdr);
int mdl_set_hs_ix(mdl_t *mdl);

/* Set the number of cell types, not including the ambient cluster. */
int mdl_set_k(mdl_t *mdl, int k);
int mdl_set_rna_atac(mdl_t *mdl, uint8_t has_rna, uint8_t has_atac);
int mdl_alloc_probs(mdl_t *mdl);

/*******************************************************************************
 * Probability functions
 ******************************************************************************/

// Pr(H_d)
void pr_hd(mdl_pars_t *mp, par_ix_t *par_ix, f_t *prob);
// Pr(S_d)
void pr_sd(mdl_pars_t *mp, par_ix_t *par_ix, f_t *prob);
// Pr(T_d)
int pr_tdm(mdl_pars_t *mp, int mol_type, int bc_ix, par_ix_t *par_ix,
        int t_ix, f_t *prob);
// get rho
f_t pr_rho_gene(f_t *rho, seq_gene_t seq_gene, uint32_t col,
        uint32_t G, uint32_t rho_nrow);
// get gamma
f_t pr_gamma_var(f_t *gamma, uint32_t v_ix, uint8_t allele, 
        uint32_t s_ix, uint32_t gamma_nrow, f_t tau);

/* Get probability of observed base calls in a molecule.
 * mp gives the model parameters and should have gamma and tau set.
 * mlcl gives the observed base calls at overlapping variants.
 * s_ix gives the sample index.
 * The resulting probability value is stored in the argument pointed to by
 * prob. If there are no overlapping variants and base calls, then 1 is
 * stored in prob.
 * If any pointer arguments are null, behaviour is undefined.
 * Returns 0 on success, -1 on error.
 */
int p_var(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int s_ix, f_t *prob);

/* Get Pr(G_dm, B_dm, X_d | \Theta) for an RNA molecule.
 * Calculate probability of observed bases and probability of observed gene.
 * `s_ix` gives the sample index, must be within range [0, M+1). Last index
 * is the ambient cluster.
 * `k` gives the cluster index, must be within range [0, K+1). The first index
 * is the ambient cluster, the rest are cell types.
 * Stores the probability in `prob`.
 * Returns 0 on success, -1 on error.
 */
int p_rna(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int s_ix, uint16_t k, f_t *prob);

/* Get Pr(P_dm, B_dm, X_d | \Theta) for an ATAC molecule.
 * Calculate probability of observed bases and probability of observed peak.
 * `s_ix` gives the sample index, must be within range [0, M+1). Last index
 * is the ambient cluster.
 * `k` gives the cluster index, must be within range [0, K+1). The first index
 * is the ambient cluster, the rest are cell types.
 * Stores the probability in `prob`.
 * Returns 0 on success, -1 on error.
 */
int p_atac(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int s_ix, f_t *prob);

/* Get Pr(B_dm) = \sum_t=1^tn Pr(T_dm = t, B_dm).
 * This calculates the probability of observing a base call B_dm = b after
 * marginalizing out T_dm.
 * mp gives the model parameters.
 * mlcl gives the observed bases.
 * rna specifies whether the molecule is RNA (=1) or ATAC (=0).
 * bc_ix gives the barcode index in D (the index across all barcodes).
 * par_ix gives the index parameters, including t_n, s_n.
 * The resulting probabilities for each T_i are stored in the array p_b, which
 * must have at least 3 elements allocated.
 * The sum of the probabilities Pr(B_dm) is stored in psum.
 * If any pointer arguments are null, behaviour is undefined.
 * Returns 0 on success, -1 on error.
 */
int p_bd(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int mol_type, int bc_ix, par_ix_t *par_ix,
        f_t *p_b, f_t *psum);

/* Calculate Pr(T_dm, F_dm, B_dm) = \sum_t=1^tn Pr(T_dm = t, F_dm, B_dm)
 * where F_dm is gene or peak probability.
 * Set rna != 0 if mlcl is RNA, and rna = 0 if mlcl is ATAC
 * @return 0 on success, -1 on error.
 */
int p_f_v(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int rna, int bc_ix, par_ix_t *par_ix,
        f_t *p_tgb, f_t *psum);

int p_kf_v(mdl_pars_t *mp, mdl_mlcl_t *mlcl, int rna, int bc_ix, par_ix_t *par_ix,
        f_t **p_tgb, uint32_t p_len, f_t *psum);

/*******************************************************************************
 * Expectation
 ******************************************************************************/

/* Calculate conditional probabilities for sub (H,S) model.
 */
int mdl_sub_e(mdl_t *mdl, int *ixs, uint32_t ix_len);

/* this calculates the conditional probabilities.
 */
int mdl_full_e(mdl_t *mdl, int *ixs, uint32_t ix_len);

/*******************************************************************************
 * Maximization
 ******************************************************************************/

/* Maximization step for sub model.
 * Maximizes parameters lambda, alpha, pi.
 */
int mdl_sub_m(mdl_t *mdl, int *ixs, uint32_t ix_len);

/* Maximization for full model.
 * Maximizes kappa, rho, and sigma.
 */
int mdl_full_m(mdl_t *mdl, int *ixs, uint32_t ix_len);

int mdl_sub_est(mdl_t *mdl);
int mdl_full_est(mdl_t *mdl);

int mdl_delta_q(f_t q1, f_t q2, f_t *q_delta);

/*******************************************************************************
 * EM functions
 ******************************************************************************/

void *mdl_sub_thrd_fx(void *arg);
int mdl_thrd_call(mdl_t *mdl, int *ixs, uint32_t ix_len, int sub);
int mdl_em(mdl_t *mdl, obj_pars *objs, int sub);

int mdl_full_llk(mdl_t *mdl, f_t *llk);
int mdl_sub_llk(mdl_t *mdl, f_t *llk);

/* Get the best (H,S) index for each droplet.
 * This function gets the index of the best `mdl.sub_lp_hs` and stores in
 * `mdl.sub_best_hs` for each of `D` droplets. 
 * If the model isn't initialized, returns -1 on error.
 * Returns 0 on success, -1 on error.
 */
int mdl_sub_best_hs(mdl_t *mdl);

/* Get the posterior probabilities Pr(H_d | X_d) for each droplet.
 * Takes the log probabilities from `mdl.sub_lp_h` and 
 * calculates each H_d = 0, 1, 2 for all droplets `D`.
 * Stores the results in `mdl.sub_pp_h`.
 * Returns 0 on success, -1 on error.
 */
int mdl_sub_pp_h(mdl_t *mdl);

/* get the best (K) index for each droplet `d` in `D`.
 * Gets the index of highest `mdl.mp.kappa` and stores in `mdl.full_best_k`.
 * Index 0 is empty, while 1, ..., K give the clusters.
 * If the model isn't initialized, returns -1 on error.
 * Returns 0 on success, -1 on error.
 */
int mdl_full_best_k(mdl_t *mdl);

/* main function to fit model.
 *
 * Call this function to set parameters and fit the model.
 *
 */
int mdl_fit(bam_data_t *bam_dat, obj_pars *objs);

// output functions
int write_lambda(mdl_t *mdl, char *fn);
int write_kappa(mdl_t *mdl, char *fn);
int write_alpha_rna(mdl_t *mdl, char *fn);
int write_alpha_atac(mdl_t *mdl, char *fn);
int write_rho(mdl_t *mdl, obj_pars *objs, char *fn);
int write_sigma(mdl_t *mdl, obj_pars *objs, char *fn);

/* Write the log likelihoods `log Pr(H_d, S_d, X_d | \Theta)` for test
 * droplets under the sub model. Write to file `fn.sample_llk.txt.gz`.
 * The test barcodes are in the rows, and columns correspond to the (H,S)
 * indices. The row names are written but the header is excluded.
 * Returns 0 on success, -1 on error.
 */
int write_sub_llk(mdl_t *mdl, char *fn);

/* Write the log likelihoods `log Pr(K_d, X_d | \Theta)` for test
 * droplets under the full model. Write to file `fn.cluster_llk.txt.gz`.
 * The test barcodes are in the rows, and columns correspond to the (K)
 * indices, ambient + cell types. The row names are written but the header
 * is excluded.
 * Returns 0 on success, -1 on error.
 */
int write_full_llk(mdl_t *mdl, char *fn);
int write_samples(mdl_t *mdl, char *fn);
int write_res(mdl_t *mdl, bam_data_t *bam_dat, char *fn);

/*******************************************************************************
 * count correction
 ******************************************************************************/

int mod_correct_counts(mdl_t *mdl, g_var_t *gv, gene_anno_t *gene_anno,
        iregs_t *pks, char *out_fn);

#endif // MOD_H

