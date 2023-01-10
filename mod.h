
#ifndef MOD_H
#define MOD_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdint.h>
#include "htslib/khash.h"
#include "gtf_anno.h"
#include "str_util.h"
#include "variants.h"
#include "clopts.h"
#include "rna_data.h"
#include "atac_data.h"
#include "bam_dat.h"
#include "bc_stats.h"

// return a random number in [0,1]
#define frand() ((double) rand() / (RAND_MAX))
// return a random number in [0,1)
#define frandlt1() ((double) rand() / (RAND_MAX+1.0))
// return a random number in (0,1)
#define frandin() ((double) (rand()+1) / (RAND_MAX+2.0))
// return a random number in [0.1, 0.9]
#define frandin2() (0.1 + (0.8 * (double)rand() / (RAND_MAX)))

#define TAU 0.01

#define GROUPN 2
enum GROUP{UNFX, AMBN};

#define JMAX 3
#define LMAX 3

#define C_AMBN 0
#define C_CELL 1

// [i,j] is [row,column]. I is number of rows, J is number of columns
// column-major indexing
#define CMI(i,j,I) ((I*j)+i) // column-major indexing
#define RMI(i,j,J) ((J*i)+j) // row-major indexing

/* j_ix indexing.
 * The single integer j_ix gives information on the k cell type and n sample.
 * ix takes integer values in 0, 1, ..., NK
 * ix corresponds to (n,k,a) values of [ (0,0,0), (0,1,0), ..., (N-1,K-2,0), (N-1,K-1,0), (N,K,1) ].
 * n and k consist of the two-dimensional indices. N and K give the max values 
 * for the indixes
 * a=1 indicates it is the NK+1 (ambient) index.
 */
void ix2indices(int ix, int M, int *s1, int *s2);

static inline void indices2ix(int *ix, int K, int k, int N, int n, int a){
    int NK = N*K;
    if (a == AMBN) *ix = NK;
    else *ix = (K * n) + k;
}

/*********************
 * parameters
 *********************/

/*
 *
 * An example workflow:
 *
 * mdl_pars_t *gp = init_mdl_pars();
 * set_mdl_pars(gp, n_bcs, n_genes, n_vars, n_samples);
 * mdl_pars_set_bcs(gp, bcs, 1);
 * mdl_pars_est_alpha(gp, bam_rna, flt_bcs, 1.0);
 * mdl_pars_set_samples(gp, samples, 1);
 * mdl_pars_add_gamma(gp, dose_array, n_vars, n_samples);
 * mdl_pars_set_beta(gp);
 * mdl_pars_set_tau(gp, 0.5);
 * mdl_pars_set_eps(gp, 0.01);
 * mdl_llk(bam_rna, bam_atac, gp, &bc_llks);
 *
 */

typedef struct {
    // Parameters
    uint32_t D; // number of droplets
    uint32_t G; // number of genes
    uint32_t V; // number of variants
    uint16_t M; // number of samples

    double *alpha; // expression prob. (G x 2 array).
    float *gamma; // fixed genotypes prob. {0,1} (M+1 x V array).
    double *beta; // prior single doublet sample prob. ( M + M*(M-1) / 2 vector).
    double tau; // prior UMI contamination prob.
    double eps; // prior for error rate.
} mdl_pars_t;

typedef struct {
    double *bc_llks; // BC x ix array

    double *best_sng_llk; // BC-length vector
    double *sec_sng_llk; // BC-length vector
    double *best_dbl_llk; // BC-length vector
    uint32_t *best_sng_ix; // BC-length vector
    uint32_t *sec_sng_ix; // BC-length vector
    uint32_t *best_dbl_ix; // BC-length vector

    double *llr; // BC-length vector
} mdl_fit_t;

typedef struct {
    str_map *all_bcs; // whitelist of barcodes
    str_map *flt_bcs; // barcodes to estimate alpha
    str_map *test_bcs; // barcodes to calculate llk
    str_map *samples;

    mdl_pars_t *mp;
    mdl_fit_t *mf;
} mdl_t;

/* Allocate and initialize and mdl_t struct
 *
 * @return Pointer to mdl_t object, or NULL on failure.
 */
mdl_t *mdl_alloc();
void mdl_dstry(mdl_t *m);

/* initialize mdl_pars struct */
mdl_pars_t *init_mdl_pars();

/* initialize mdl_pars struct
 * @param D Number of droplets
 * @param G Number of genes
 * @param V Number of variants
 * @param M Number of samples
 * @return
 * @return pointer to initialized mdl_pars struct, NULL on failure.
 */
int set_mdl_pars(mdl_pars_t *gp, uint32_t D, uint32_t G, uint32_t V, uint16_t M);

/* destroy mdl_pars struct
 * @param gp Pointer to mdl_pars struct.
 */
void destroy_mdl_pars(mdl_pars_t *gp);

/* initialize mdl_fit_t */
mdl_fit_t *init_mdl_fit();
int set_mdl_fit(mdl_fit_t *mf, str_map *bcs, str_map *samples);

/* destroy mdl_fit_t */
void destroy_mdl_fit(mdl_fit_t *mf);

/* Estimate alpha parameter.
 *
 * The alpha parameter gives the gene expression probabilities for the 
 * background and cell/nucleus groups.
 *
 * Assumes the barcodes in @p flt_bcs contain cells/nuclei, and are used 
 * to estimate alpha in the cell group. Barcodes not in flt_bcs are used 
 * to estimate alpha in the ambient group.
 *
 * @param gp Pointer to mdl_pars struct.
 * @param br Pointer to bam_rna_t object that contains aligned gene for each
 *  molecule in each barcode.
 * @param flt_bcs A str_map containing a filtered set of barcodes to use for 
 *  estimating alpha.
 * @param smooth add to estimate.
 * @return -1 on error, 0 if success.
 */
int mdl_pars_est_alpha(mdl_pars_t *gp, bam_data_t *bam_data, str_map *flt_bcs, 
        double smooth);

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

/* Set beta parameter
 *
 * Beta is the prior probability of singlet or doublet. This sets a 
 * uniform prior, so equal probability for any singlet or doublet configuration.
 *
 * @return 0 on success, -1 on error.
 */
int mdl_pars_set_beta(mdl_pars_t *gp);

/* Set tau parameter
 * 
 * Tau is prior probability of contamination.
 */
int mdl_pars_set_tau(mdl_pars_t *gp, double tau);

/* Set eps parameter
 *
 * Eps is probability of sequencing error.
 */
int mdl_pars_set_eps(mdl_pars_t *gp, double eps);

int p_brnli(int x, double p, double *ret);

/* calculate P(B | gamma, S, E, C)
 *
 * If genotype of the variant is missing, then a probability of 1 is 
 * returned. This applies to singlets, and when at least one of the 
 * samples of the doublet pair has a missing genotype, then a 
 * probability of 1 is also returned. This is equivalent to marginalizing 
 * out the allele variable.
 *
 * @param allele int 0 or 1 only.
 * @param v variant index.
 * @param gp pointer to mdl_pars struct.
 * @param s sample index.
 * @param e 0 or 1 error variable.
 * @param c 0 or 1 contamination variable.
 * @param prob pointer to double where prob is updated after call.
 * @return 0 on success, 1 on error.
 */
int p_b_gce(int allele, int v, mdl_pars_t *gp, int s, int e, int c, double *prob);

/* calculate P(E | epsilon)
 *
 * @param e error 0 or 1.
 * @param eps prior probability parameter for error (e = 1).
 * @param prob pointer to double where prob is updated after call.
 * @return 0 on success, 1 on error.
 */
int p_e(int e, double eps, double *prob);


int p_f_c(int f, int spl, mdl_pars_t *gp, int c, double *prob);


int p_c(int c, double tau, double *prob);

int mdl_llk(mdl_t *mdl, bam_data_t *bam_dat, int verbose);
int mdl_get_best_llk(mdl_t *mdl);

/* Set filtered and all barcodes in the mdl struct
 *
 * @return 0 on success, -1 on error.
 */
int mdl_set_bcs(mdl_t *mdl, bam_data_t *bam_dat, obj_pars *objs);

/* Set samples in mdl struct
 *
 * @return 0 on success, -1 on error
 */
int mdl_set_samples(mdl_t *mdl, obj_pars *objs);

/* Set parameters for model
 *
 * Set numbers, estimate alpha, set famma, set beta, set tau, and set eps.
 */
int mdl_est_pars(mdl_t *mdl, bam_data_t *bam_dat, obj_pars *objs);

/* main function to fit model.
 *
 * Call this function to set parameters and fit the model.
 *
 */
int fit_mdl_pars(bam_data_t *bam_dat, obj_pars *objs);

int write_alpha(mdl_t *mdl, obj_pars *objs, char *fn);

int write_llk(mdl_t *mdl, char *fn);

int write_samples(mdl_t *mdl, char *fn);

int write_res(mdl_t *mdl, bam_data_t *bam_dat, char *fn);

#endif // MOD_H

