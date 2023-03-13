
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

mdl_pars_t *mdl_pars_init(){
    mdl_pars_t *gp = calloc(1, sizeof(mdl_pars_t));
    int err;
    if ((err = pthread_mutex_init(&gp->sum_lock, NULL)) != 0){
        err_msg(-1, 0, "mdl_pars_init: failed to initialize mutex %i", err);
        return(NULL);
    }
    return gp;
}

void destroy_mdl_pars(mdl_pars_t *gp){
    if (gp == NULL) return;

    free(gp->pi);
    free(gp->alpha);
    free(gp->rho);
    free(gp->gamma);
    
    pthread_mutex_destroy(&gp->sum_lock);
    free(gp->_pi_sum);
    free(gp->_alpha_sum);
    free(gp->_rho_sum);
    free(gp);
}

mdl_fit_t *init_mdl_fit(){
    mdl_fit_t *mf = calloc(1, sizeof(mdl_fit_t));
    return mf;
}

void destroy_mdl_fit(mdl_fit_t *mf){
    if (mf == NULL) return;
    free(mf->bc_llks);
    free(mf->best_sng_llk);
    free(mf->sec_sng_llk);
    free(mf->best_dbl_llk);
    free(mf->best_sng_ix);
    free(mf->sec_sng_ix);
    free(mf->best_dbl_ix);
    free(mf->pp);
    free(mf);
}

mdl_t *mdl_alloc(){
    mdl_t *mdl = (mdl_t *)calloc(1, sizeof(mdl_t));
    if (mdl == NULL){
        err_msg(-1, 0, "mdl_alloc: %s", strerror(errno));
        return(NULL);
    }

    mdl->all_bcs = NULL;
    mdl->test_bcs = NULL;
    mdl->samples = NULL;

    mdl->fix_flag = NULL;

    mdl->c_probs_len = 0;
    mdl->c_probs = NULL;

    mdl->_nrow_hs = 0;
    mdl->_nrow_hst = 0;
    mdl->hs_ix = NULL;

    mdl->eps = 1e-5;
    mdl->max_iter = 100;
    mdl->has_rna = 0;
    mdl->has_atac = 0;
    mdl->threads = 1;

    mdl->mp = calloc(1, sizeof(mdl_pars_t));
    mdl->mf = calloc(1, sizeof(mdl_fit_t));
    if (mdl->mp == NULL || mdl->mf == NULL)
        return(NULL);

    return(mdl);
}

void mdl_dstry(mdl_t *m){
    if (m == NULL) return;

    if (m->all_bcs) destroy_str_map(m->all_bcs);
    if (m->test_bcs) destroy_str_map(m->test_bcs);
    if (m->samples) destroy_str_map(m->samples);

    bflg_free(m->fix_flag);
    bflg_free(m->absent_bc);
    free(m->fix_flag);
    free(m->absent_bc);

    uint32_t i;
    for (i = 0; i < m->c_probs_len; ++i){
        free(m->c_probs[i].cp_hs);
        free(m->c_probs[i].lp_hs);
    }
    free(m->c_probs);

    free(m->hs_ix);

    destroy_mdl_pars(m->mp);
    destroy_mdl_fit(m->mf);
    free(m);
}

int mdl_set_bcs(mdl_t *mdl, bam_data_t *bam_data, obj_pars *objs){
    if (mdl == NULL || bam_data == NULL || objs == NULL)
        return err_msg(-1, 0, "mdl_set_bcs: arguments are null");

    if (bam_data->has_stats == 0)
        return err_msg(-1, 0, "mdl_set_bcs: 'bam_data' stats are required");

    khint_t k;

    // set all barcodes for model
    int found;
    if (objs->wl_bcs != NULL){
        // set to wl_bc if given
        mdl->all_bcs = str_map_copy(objs->wl_bcs);
        if (mdl->all_bcs == NULL)
            return err_msg(-1, 0, "mdl_set_bcs: failed to copy all_bcs");
    } else {
        if ( (mdl->all_bcs = init_str_map()) == NULL )
            return(-1);

        // if wl_bc not given, set to all barcodes present in data
        for (k = kh_begin(bam_data->bc_data); k != kh_end(bam_data->bc_data); ++k){
            if (!kh_exist(bam_data->bc_data, k)) continue;

            char *bc_key = kh_key(bam_data->bc_data, k);

            if (add2str_map(mdl->all_bcs, bc_key, &found) < 0)
                return(-1);
            assert(found != 1);
        }
    }

    // set the output barcodes to calculate llk for
    int i, bc_n = mdl->all_bcs->n;

    mdl->fix_flag = calloc(1, sizeof(bflg_t));
    mdl->absent_bc = calloc(1, sizeof(bflg_t));
    if (mdl->fix_flag == NULL || mdl->absent_bc == NULL)
        return err_msg(-1, 0, "mdl_set_bcs: %s", strerror(errno));

    if (bflg_init(mdl->fix_flag, bc_n) < 0)
        return err_msg(-1, 0, "mdl_set_bcs: failed to init fix flag");
    if (bflg_init(mdl->absent_bc, bc_n) < 0)
        return err_msg(-1, 0, "mdl_set_bcs: failed to init absent flag");

    if (objs->out_min < 0)
        return err_msg(-1, 0, "mdl_set_bcs: out-min=%i is less than 0", objs->out_min);
    uint32_t c_thresh = (uint32_t)(objs->out_min);
    mdl->test_bcs = init_str_map();
    for (i = 0; i < bc_n; ++i){
        char *bc_key = str_map_str(mdl->all_bcs, i);
        assert(bc_key != NULL);

        khint_t k = kh_get(kh_bc_dat, bam_data->bc_data, bc_key);

        // if barcode is not present in data
        if (k == kh_end(bam_data->bc_data)){
            bflg_set(mdl->fix_flag, i);
            continue;
        }

        bc_data_t *bc_data = kh_val(bam_data->bc_data, k);
        assert(bc_data != NULL);
        bc_stats_t *bc_stat = bc_data->bc_stats;
        assert(bc_stat != NULL);
        if (bc_stat->atac_counts < c_thresh && bc_stat->rna_counts < c_thresh){
            bflg_set(mdl->fix_flag, i);
        } else {
            if (add2str_map(mdl->test_bcs, bc_key, &found) < 0) return(-1);
            assert(found != 1);
        }
    }

    return(0);
}

int mdl_set_samples(mdl_t *mdl, obj_pars *objs){
    if (mdl == NULL || objs == NULL)
        return err_msg(-1, 0, "mdl_set_samples: arguments are NULL");

    if (objs->vcf_hdr == NULL)
        return err_msg(-1, 0, "mdl_set_samples: vcf header is not present");

    int n_samples = bcf_hdr_nsamples(objs->vcf_hdr);
    mdl->samples = init_str_map_array(objs->vcf_hdr->samples, n_samples);
    if (mdl->samples == NULL)
        return err_msg(-1, 0, "mdl_set_samples: failed to copy samples from vcf header");

    return(0);
}

int mdl_init_all(mdl_t *mdl, bam_data_t *bam_data, obj_pars *objs){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_init_all: mdl is null");

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "mdl_init_all: mdl->mp is null");
    mdl_pars_t *pars = mdl->mp;

    if (objs == NULL)
        return err_msg(-1, 0, "mdl_init_all: objs is null");

    /* */
    if (mdl->all_bcs == NULL || mdl->test_bcs == NULL)
        return err_msg(-1, 0, "mdl_init_all: mdl->all_bcs or mdl->test_bcs is null, call mdl_set_bcs");

    uint32_t i, j;

    uint32_t n_bcs = (uint32_t)mdl->all_bcs->n;
    if (n_bcs == 0)
        return err_msg(-1, 0, "mdl_init_all: no barcodes");
    uint32_t n_genes = 0;
    if (objs->anno)
        n_genes = (uint32_t)objs->anno->gene_ix->n;
    uint16_t n_samples = (uint16_t)bcf_hdr_nsamples(objs->vcf_hdr);
    if (n_samples == 0)
        return err_msg(-1, 0, "mdl_init_all: no samples");

    uint32_t n_vars = (uint32_t)mv_size(&objs->gv->vix2var);

    pars->D = n_bcs;
    pars->G = n_genes;
    pars->M = n_samples;
    pars->V = n_vars;
    uint32_t gene_max = N_SPL * pars->G;

    /* */

    // indices
    mdl->_nrow_hs = 1 + n_samples + ( n_samples * (n_samples - 1) / 2);
    mdl->hs_ix = calloc(3 * mdl->_nrow_hs, sizeof(int));
    mdl->_nrow_hst = 1 + (2 * n_samples) + (3 * n_samples * (n_samples - 1) / 2);

    int s1 = 0, s2 = 1;
    uint32_t ixi = 0;
    mdl->hs_ix[CMI(ixi, 0, mdl->_nrow_hs)] = 0;
    mdl->hs_ix[CMI(ixi, 1, mdl->_nrow_hs)] = -1;
    mdl->hs_ix[CMI(ixi, 2, mdl->_nrow_hs)] = -1;
    ++ixi;

    for (s1 = 0; s1 < n_samples; ++s1){
        mdl->hs_ix[CMI(ixi, 0, mdl->_nrow_hs)] = 1;
        mdl->hs_ix[CMI(ixi, 1, mdl->_nrow_hs)] = s1;
        mdl->hs_ix[CMI(ixi, 2, mdl->_nrow_hs)] = -1;
        ++ixi;
    }
    for (s1 = 0; s1 < n_samples - 1; ++s1){
        for (s2 = s1 + 1; s2 < n_samples; ++s2){
            mdl->hs_ix[CMI(ixi, 0, mdl->_nrow_hs)] = 2;
            mdl->hs_ix[CMI(ixi, 1, mdl->_nrow_hs)] = s1;
            mdl->hs_ix[CMI(ixi, 2, mdl->_nrow_hs)] = s2;
            ++ixi;
        }
    }
    assert(ixi == mdl->_nrow_hs);

    pars->pi = calloc(pars->M, sizeof(f_t));
    pars->alpha = calloc(pars->D, sizeof(f_t));
    pars->rho = calloc(gene_max * 2, sizeof(f_t));

    pars->_pi_sum = calloc(pars->M, sizeof(f_t));
    pars->_alpha_sum = calloc(2 * pars->D, sizeof(f_t));
    pars->_rho_sum = calloc(gene_max * 2, sizeof(f_t));
    if (pars->_rho_sum == NULL)
        return err_msg(-1, 0, "mdl_pars_init_all: %s", strerror(errno));

    // conditional probabilities per barcode
    mdl->c_probs_len = n_bcs;
    mdl->c_probs = calloc(n_bcs, sizeof(c_probs_t));
    if (mdl->c_probs == NULL)
        return err_msg(-1, 0, "mdl_pars_init_all: %s", strerror(errno));

    for (i = 0; i < (uint32_t)mdl->all_bcs->n; ++i){
        mdl->c_probs[i].p_x = 0;
        int fl;
        fl = bflg_get(mdl->fix_flag, i);
        uint32_t n_hs = fl == 1 ? 1 : mdl->_nrow_hs;
        mdl->c_probs[i].lp_hs = calloc(n_hs, sizeof(f_t));
        mdl->c_probs[i].cp_hs = calloc(n_hs, sizeof(f_t));
        if (mdl->c_probs[i].cp_hs == NULL)
            return err_msg(-1, 0, "mdl_pars_init_all: %s", strerror(errno));
        mdl->c_probs[i].lp_hs[0] = 0;
        for (j = 1; j < n_hs; ++j) mdl->c_probs[i].lp_hs[j] = -INFINITY;
        mdl->c_probs[i].cp_hs[0] = 1;
        for (j = 1; j < n_hs; ++j) mdl->c_probs[i].cp_hs[j] = 0;
    }

    mdl->eps = objs->eps;
    mdl->max_iter = objs->max_iter;

    mdl->has_rna = bam_data->has_rna;
    mdl->has_atac = bam_data->has_atac;

    // Initialize/set gamma
    // first index: variant, second index: sample
    int32_t *vixs = calloc(n_vars, sizeof(int32_t));
    for (i = 0; i < n_vars; ++i) vixs[i] = i;
    float **gm = ap_array_gt(objs->gv, objs->vcf_hdr, vixs, n_vars, "GT");
    free(vixs);
    if (gm == NULL)
        return(-1);

    if (mdl_pars_add_gamma(pars, gm, n_vars, n_samples) < 0)
        return(-1);

    // free gm
    for (i = 0; i < (uint32_t)n_vars; ++i){
        free(gm[i]);
    }
    free(gm);
    
    if (mdl_init_par_dat(mdl, bam_data) < 0) return(-1);

    mdl->threads = objs->threads;

    return(0);
}

int mdl_get_hst(mdl_t *mdl, int hs_ix, int *hd, int *s1, int *s2, int t_ix[3], int *t_n){
    // get (h,s,t) from index

    uint16_t M = mdl->mp->M;

    // -1 for if NA/invalid
    *hd = mdl->hs_ix[CMI(hs_ix, 0, mdl->_nrow_hs)];
    *s1 = mdl->hs_ix[CMI(hs_ix, 1, mdl->_nrow_hs)];
    *s2 = mdl->hs_ix[CMI(hs_ix, 2, mdl->_nrow_hs)];

    switch (*hd) {
        case 0:
            t_ix[0] = M;
            t_ix[1] = -1;
            t_ix[2] = -1;
            *t_n = 1;
            break;
        case 1:
            t_ix[0] = M;
            t_ix[1] = *s1;
            t_ix[2] = -1;
            *t_n = 2;
            break;
        case 2:
            t_ix[0] = M;
            t_ix[1] = *s1;
            t_ix[2] = *s2;
            *t_n = 3;
            break;
        default:
            return err_msg(-1, 0, "mdl_get_hst: hd=%i is invalid, there is a bug", hd);
    }
    return(0);
}

/*******************************************************************************
 * Probability functions
 ******************************************************************************/

void pr_hd(mdl_t *mdl, int hd, f_t *prob){
    *prob = mdl->mp->lambda[hd];
}

// not in log
void pr_sd(mdl_t *mdl, int hd, int s1, int s2, f_t *prob){
    f_t p_sd, p_hd = mdl->mp->lambda[hd];
    switch (hd) {
        case 0:
            p_sd = 1;
            break;
        case 1:
            p_sd = mdl->mp->pi[s1];
            break;
        case 2:
            p_sd = (mdl->mp->pi[s1] * mdl->mp->pi[s2]) / mdl->mp->pi_d_sum;
            break;
        default:
            p_sd = -1;
    }
    *prob = p_sd;
}

// not in log
int pr_tdm(mdl_t *mdl, int bc_ix, int hd, int s_ix, f_t *prob){

    uint16_t M = mdl->mp->M;
    uint8_t c_ix = s_ix == M ? 0 : 1; // 0 is ambient, 1 is cell
    f_t al = mdl->mp->alpha[bc_ix];
    f_t pr = 0;
    switch (hd) {
        case 0:
            if (c_ix == 0) pr = 1.0;
            else
                return err_msg(-1, 0, "pr_tdm: hd = 0 but s_ix != M", hd);
            break;
        case 1:
            if (c_ix == 0) pr = al;
            else pr = (1.0 - al);
            break;
        case 2:
            if (c_ix == 0) pr = al;
            else pr = ((1.0-al)/2.0);
            break;
        default:
            return err_msg(-1, 0, "pr_tdm: hd=%i is invalid, there is a bug", hd);
    }
    *prob = pr;
    return(0);
}

// not in log
int p_rna(mdl_t *mdl, rna_mol_t *mol, int s_ix, f_t *prob){
    uint16_t M = mdl->mp->M, M1 = M + 1;
    uint32_t G = mdl->mp->G;
    uint32_t rho_nrow = 3 * G;
    uint8_t c_ix = s_ix == M ? 0 : 1; // 0 is ambient, 1 is cell

    f_t lp = 1; // prob

    // probability terms
    f_t p_gdm = 1, p_bdm = 1;

    // Pr(G_dm | T_dm, \rho)
    if (ml_size(&mol->gl) == 1){
        ml_node_t(seq_gene_l) *gn = ml_begin(&mol->gl);
        seq_gene_t seq_gene = ml_node_val(gn);
        p_gdm = pr_rho_gene(mdl->mp->rho, seq_gene, (uint32_t)c_ix, G, rho_nrow);
    }
    lp *= p_gdm;

    // Pr(B_dm | T_dm, \gamma, pi)
    ml_node_t(seq_vac_l) *vn;
    for (vn = ml_begin(&mol->vl); vn; vn = ml_node_next(vn)){
        seq_vac_t vac = ml_node_val(vn);
        p_bdm = pr_gamma_var(mdl->mp->gamma, vac, (uint32_t)s_ix, M1, mdl->mp->tau);
        lp *= p_bdm;
    }

    *prob = lp;
    return(0);
}

// not in log
int p_atac(mdl_t *mdl, atac_frag_t *frag, int s_ix, f_t *prob){

    size_t n_peaks = mv_size(&frag->pks);

    uint16_t M = mdl->mp->M, M1 = M + 1;
    uint8_t c_ix = s_ix == M ? 0 : 1; // 0 is ambient, 1 is cell

    f_t lp = 1; // prob

    // probability terms
    f_t p_pdm = 1, p_bdm = 1;

    // Pr(P_dm | T_dm, \sigma)
    p_pdm = mdl->mp->sigma[c_ix]; // frag in peak
    if (n_peaks == 0) p_pdm = 1 - p_pdm; // frag out peak
    lp *= p_pdm;

    // Pr(B_dm | T_dm, \gamma, pi)
    ml_node_t(seq_vac_l) *vn;
    for (vn = ml_begin(&frag->vl); vn; vn = ml_node_next(vn)){
        seq_vac_t vac = ml_node_val(vn);
        if (vac.allele > 2) continue; // other alleles are considered missing
        f_t ap = mdl->mp->gamma[ CMI(s_ix, vac.vix, M1) ];
        if (ap < 0) continue; // if allele is missing
        if (vac.allele == 0) ap = 1.0 - ap;
        f_t p_be0 = (1.0 - mdl->mp->tau) * ap;
        f_t p_be1 = mdl->mp->tau * 0.25;
        p_bdm = p_be0 + p_be1;
        lp *= p_bdm;
    }
    *prob = lp;
    return(0);
}

/*******************************************************************************
 * Expectation
 ******************************************************************************/

int mdl_e_hs(mdl_t *mdl, bam_data_t *bam_data, int *ixs, uint32_t ix_len){
    if (mdl == NULL || bam_data == NULL || ixs == NULL)
        return err_msg(-1, 0, "mdl_e_hs: arguments are null");

    assert(mdl->mp != NULL);
    assert(mdl->all_bcs != NULL);
    assert(mdl->c_probs != NULL);
    assert(bam_data->bc_data != NULL);

    uint16_t M = mdl->mp->M, M1 = M + 1;

    int lsret = 0;

    // pre-calculate P(H_d, S_d)
    f_t *lp_hs_v = calloc(mdl->_nrow_hs, sizeof(f_t));
    uint32_t hs_i = 0;
    for (hs_i = 0; hs_i < mdl->_nrow_hs; ++hs_i){
        int hd, s1, s2, t_ix[3], t_n;
        if (mdl_get_hst(mdl, hs_i, &hd, &s1, &s2, t_ix, &t_n) < 0)
            return(-1);
        f_t p_hd = -1, p_sd = -1;
        pr_hd(mdl, hd, &p_hd);
        pr_sd(mdl, hd, s1, s2, &p_sd);
        lp_hs_v[hs_i] = log(p_hd) + log(p_sd);
    }

    f_t *p_tm_v[3];
    uint32_t h_i = 0, t_i = 0;
    for (h_i = 0; h_i < 3; ++h_i) p_tm_v[h_i] = calloc(M1, sizeof(f_t));

    uint32_t i, n_hs;
    khint_t k;
    for (i = 0; i < ix_len; ++i){
        // get barcode and bam data
        int bc_ix = ixs[i];

        uint32_t hs_ix;
        int fl;

        // if no barcode
        fl = bflg_get(mdl->absent_bc, bc_ix);
        if (fl == 1) continue;

        // if fixed, set Pr(empty) = 1
        fl = bflg_get(mdl->fix_flag, bc_ix);
        n_hs = fl == 1 ? 1 : mdl->_nrow_hs;

        mdl->c_probs[bc_ix].cp_hs[0] = 1;
        for (hs_ix = 1; hs_ix < n_hs; ++hs_ix)
            mdl->c_probs[bc_ix].cp_hs[hs_ix] = 0;

        char *bc_key = str_map_str(mdl->all_bcs, bc_ix);
        k = kh_get(kh_bc_dat, bam_data->bc_data, bc_key);
        if (k == kh_end(bam_data->bc_data)) continue;

        // get barcode data
        bc_data_t *bc_dat = kh_val(bam_data->bc_data, k);
        assert(bc_dat != NULL);
        khash_t(khrmn) *mols = bc_dat->rna_mols;
        assert(mols != NULL);
        khash_t(khaf) *frags = bc_dat->atac_frags;
        assert(frags != NULL);

        // pre-calculate P(T_dm) for this barcode
        h_i = 0;
        if (pr_tdm(mdl, bc_ix, h_i, M, p_tm_v[h_i] + M) < 0)
            return(-1);
        for (h_i = 1; h_i < 3; ++h_i){
            for (t_i = 0; t_i < M1; ++t_i){
                if (pr_tdm(mdl, bc_ix, h_i, t_i, p_tm_v[h_i] + t_i) < 0)
                    return(-1);
            }
        }

        int moli = 0;
        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
            // get (h,s,t) from index
            // hd: num nuc. in droplet; s1: sample 1; s2: sample 2; t: mol sample
            int hd, s1, s2, t_ix[3], t_n;
            if (mdl_get_hst(mdl, hs_ix, &hd, &s1, &s2, t_ix, &t_n) < 0)
                return(-1);

            // hold Pr(H_d, S_d, X_d)
            f_t lp_htd = 0;

            // Pr(H_d, S_d | \lambda, \pi)
            lp_htd += lp_hs_v[hs_ix];
            moli=0;

            ml_node_t(mlar) *rn;
            for (rn = ml_begin(&bc_dat->mols_l); rn; rn = ml_node_next(rn)){
                rna_mol_t *mol = ml_node_val(rn);

                f_t psum = 0, p_tgpb[3] = {1,1,1};
                int t_i;
                for (t_i = 0; t_i < t_n; ++t_i){
                    int s_ix = t_ix[t_i];
                    f_t p_t = p_tm_v[hd][s_ix];
                    f_t p_gv = 0;
                    if (p_rna(mdl, mol, s_ix, &p_gv) < 0)
                        return(-1);
                    p_tgpb[t_i] = p_t * p_gv;
                    psum += p_tgpb[t_i];
                }
                lp_htd += log(psum);
                // ++moli;
            }
            ml_node_t(mlaf) *fn;
            for (fn = ml_begin(&bc_dat->frags_l); fn; fn = ml_node_next(fn)){
                atac_frag_t *frag = ml_node_val(fn);

                f_t psum = 0, p_tgpb[3] = {1,1,1};
                int t_i;
                for (t_i = 0; t_i < t_n; ++t_i){
                    int s_ix = t_ix[t_i];
                    f_t p_t = p_tm_v[hd][s_ix];
                    f_t p_pv = 0;
                    if (p_atac(mdl, frag, s_ix, &p_pv) < 0)
                        return(-1);
                    p_tgpb[t_i] = p_t * p_pv;
                    psum += p_tgpb[t_i];
                }
                // lp_htd += logsumexpd2(lp_tgpb, t_n);
                lp_htd += log(psum);
                ++moli;
            }
            mdl->c_probs[bc_ix].lp_hs[hs_ix] = lp_htd;
        }
        // Pr(X_d | \Theta)
        mdl->c_probs[bc_ix].p_x = logsumexpd2(mdl->c_probs[bc_ix].lp_hs, n_hs);
        assert(!isnan(mdl->c_probs[bc_ix].p_x));
        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
            f_t ldiv = mdl->c_probs[bc_ix].lp_hs[hs_ix] - 
                mdl->c_probs[bc_ix].p_x;
            mdl->c_probs[bc_ix].cp_hs[hs_ix] = exp(ldiv);
        }
    }
    free(lp_hs_v);
    for (t_i = 0; t_i < 3; ++t_i) free(p_tm_v[t_i]);
    return(0);
}

/*******************************************************************************
 * Maximization
 ******************************************************************************/

int mdl_m_init(mdl_t *mdl, f_t prior){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_m_init: mdl is null");

    uint32_t D = mdl->mp->D, G = mdl->mp->G;
    uint32_t rho_ne = 3 * G * 2;
    uint32_t alpha_ne = 2 * D;
    uint16_t M = mdl->mp->M;
    uint32_t i;
    /*
    for (i = 0; i < 3; ++i) mdl->mp->_lambda_sum[i] = log(prior / 3.0);
    for (i = 0; i < M; ++i) mdl->mp->_pi_sum[i] = log(prior / M);
    for (i = 0; i < alpha_ne; ++i) mdl->mp->_alpha_sum[i] = log(prior / 2.0);
    for (i = 0; i < rho_ne; ++i) mdl->mp->_rho_sum[i] = log(prior / (G * 3));
    for (i = 0; i < 4; ++i) mdl->mp->_sigma_sum[i] = log(prior / 2.0);
    */
    for (i = 0; i < 3; ++i) mdl->mp->_lambda_sum[i] = prior / 3.0;
    for (i = 0; i < M; ++i) mdl->mp->_pi_sum[i] = prior / M;
    for (i = 0; i < alpha_ne; ++i) mdl->mp->_alpha_sum[i] = prior / 2.0;
    for (i = 0; i < rho_ne; ++i) mdl->mp->_rho_sum[i] = prior / (G * 3);
    for (i = 0; i < 4; ++i) mdl->mp->_sigma_sum[i] = prior / 2.0;

    return(0);
}

int mdl_m_sum(mdl_t *mdl, bam_data_t *bam_data, int *ixs, uint32_t ix_len){
    if (mdl == NULL || ixs == NULL)
        return err_msg(-1, 0, "mdl_max_lambda: arguments are null");

    uint32_t i;
    uint16_t M = mdl->mp->M, M1 = M + 1;
    uint32_t G = mdl->mp->G;
    uint32_t D = mdl->mp->D;
    uint32_t rho_ncol = 2;
    uint32_t rho_nrow = 3 * G;
    uint32_t rho_ne = rho_ncol * rho_nrow;
    f_t *rho_sum = malloc(rho_ne * sizeof(f_t));
    // for (i = 0; i < rho_ne; ++i) rho_sum[i] = -INFINITY;
    for (i = 0; i < rho_ne; ++i) rho_sum[i] = 0;

    uint32_t sig_ncol = 2;
    uint32_t sig_nrow = 2;
    uint32_t sig_ne = sig_ncol * sig_nrow;
    // sig_sum col 0: ambient, col 1: cell, row 0: outside peak, row 1: inside peak
    f_t *sig_sum = malloc(sig_ne * sizeof(f_t)); // outside/inside peak by ambient/nuclear
    // for (i = 0; i < sig_ne; ++i) sig_sum[i] = -INFINITY;
    for (i = 0; i < sig_ne; ++i) sig_sum[i] = 0;

    f_t *p_tm_v[3];
    uint32_t h_i = 0, t_i = 0;
    for (h_i = 0; h_i < 3; ++h_i) p_tm_v[h_i] = calloc(M1, sizeof(f_t));

    uint32_t hs_ix;
    int lsret = 0;
    khint_t k;
    // f_t lambda_sums[3] = {-INFINITY, -INFINITY, -INFINITY};
    f_t lambda_sums[3] = {0,0,0};
    for (i = 0; i < ix_len; ++i){
        uint32_t bc_ix = ixs[i];

        int fl;

        // if no barcode
        fl = bflg_get(mdl->absent_bc, bc_ix);
        if (fl == 1) continue;

        c_probs_t c_p = mdl->c_probs[bc_ix];

        char *bc_key = str_map_str(mdl->all_bcs, bc_ix);
        k = kh_get(kh_bc_dat, bam_data->bc_data, bc_key);
        if (k == kh_end(bam_data->bc_data)) continue;

        fl = bflg_get(mdl->fix_flag, bc_ix);
        uint32_t n_hs = fl == 1 ? 1 : mdl->_nrow_hs;

        // if (!fl){
        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
            int h_d = mdl->hs_ix[CMI(hs_ix, 0, mdl->_nrow_hs)];
            lambda_sums[h_d] += c_p.cp_hs[hs_ix];
            // lambda_sums[h_d] = logsum2expd(lambda_sums[h_d], c_p.cp_hs[hs_ix]);
        }
        // }

        // f_t alpha_sums[2] = {-INFINITY,-INFINITY}; // index 0: ambient, index 1: cell
        f_t alpha_sums[2] = {0, 0}; // index 0: ambient, index 1: cell

        // bam data
        bc_data_t *bc_dat = kh_val(bam_data->bc_data, k);

        // pre-calculate P(T_dm) for this barcode
        h_i = 0;
        if (pr_tdm(mdl, bc_ix, h_i, M, p_tm_v[h_i] + M) < 0)
            return(-1);
        for (h_i = 1; h_i < 3; ++h_i){
            for (t_i = 0; t_i < M1; ++t_i){
                if (pr_tdm(mdl, bc_ix, h_i, t_i, p_tm_v[h_i] + t_i) < 0)
                    return(-1);
            }
        }

        for (hs_ix = 0; hs_ix < n_hs; ++hs_ix){
            // get (h,s,t) from index
            int hd, s1, s2, t_ix[3], t_n;
            if (mdl_get_hst(mdl, hs_ix, &hd, &s1, &s2, t_ix, &t_n) < 0)
                return(-1);

            // Get P(H,S,X)
            // Divide by \sum P(T,G,P,B,X)
            // Multiple by P(T,G,P,B,X)
            f_t lp_hsx = mdl->c_probs[bc_ix].cp_hs[hs_ix];

            ml_node_t(mlar) *rn;
            for (rn = ml_begin(&bc_dat->mols_l); rn; rn = ml_node_next(rn)){
                rna_mol_t *mol = ml_node_val(rn);

                size_t n_genes = ml_size(&mol->gl);
                if (n_genes > 1) n_genes = 0; // skip multigene UMIs

                // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
                f_t psum = 0, p_tgpb[3] = {0,0,0};
                int t_i;
                for (t_i = 0; t_i < t_n; ++t_i){
                    int s_ix = t_ix[t_i];
                    f_t p_t = p_tm_v[hd][s_ix];
                    f_t p_gv = 0;
                    if (p_rna(mdl, mol, s_ix, &p_gv) < 0)
                        return(-1);
                    p_tgpb[t_i] = p_t * p_gv;
                    psum += p_tgpb[t_i];
                }
                // f_t lp_tgbx = log(psum);
                // f_t lp1 = lp_hsx - lp_tgbx;
                f_t phs1 = lp_hsx / psum;
                for (t_i = 0; t_i < t_n; ++t_i){
                    int s_ix = t_ix[t_i];
                    int a_type = s_ix == M ? 0 : 1; // 0 for ambient, 1 for cell
                    // f_t lpt = lp1 + log(p_tgpb[t_i]);
                    f_t pt = phs1 * p_tgpb[t_i];
                    assert(!isnan(pt));

                    // alpha
                    if (hd >= 1){
                        // alpha_sums[a_type] = logsum2expd(alpha_sums[a_type], lpt);;
                        alpha_sums[a_type] += pt;
                    }

                    // rho
                    if (n_genes == 1){
                        ml_node_t(seq_gene_l) *gn = ml_begin(&mol->gl);
                        seq_gene_t seq_gene = ml_node_val(gn);
                        int32_t g_ix = seq_gene.gene_id;
                        uint8_t spl = seq_gene.splice;
                        uint32_t f_ix = g_ix + (spl * mdl->mp->G);
                        // rho_sum[ CMI(f_ix, a_type, rho_nrow) ] = 
                        //     logsum2expd(rho_sum[ CMI(f_ix, a_type, rho_nrow) ], lpt);
                        rho_sum[ CMI(f_ix, a_type, rho_nrow) ] += pt;
                    }
                }
            }
            ml_node_t(mlaf) *fn;
            for (fn = ml_begin(&bc_dat->frags_l); fn; fn = ml_node_next(fn)){
                atac_frag_t *frag = ml_node_val(fn);


                size_t n_peaks = mv_size(&frag->pks);
                uint8_t pk = n_peaks == 0 ? 0 : 1; // 0: outside peak, 1: inside peak

                // Recalculate Pr(T_dm, P_dm, B_dm | \Theta)
                f_t psum = 0, p_tgpb[3] = {0,0,0};
                int t_i;
                for (t_i = 0; t_i < t_n; ++t_i){
                    int s_ix = t_ix[t_i];
                    f_t p_t = p_tm_v[hd][s_ix];
                    f_t p_pv = 0;
                    if (p_atac(mdl, frag, s_ix, &p_pv) < 0)
                        return(-1);
                    p_tgpb[t_i] = p_t * p_pv;
                    psum += p_tgpb[t_i];
                }
                // f_t lp_tgbx = log(psum);
                // f_t lp1 = lp_hsx - lp_tgbx;
                f_t phs1 = lp_hsx / psum;
                for (t_i = 0; t_i < t_n; ++t_i){
                    int s_ix = t_ix[t_i];
                    int a_type = s_ix == M ? 0 : 1; // 0 for ambient, 1 for cell
                    // f_t lpt = lp1 + log(p_tgpb[t_i]);
                    f_t pt = phs1 * p_tgpb[t_i];

                    // alpha
                    if (hd >= 1){
                        // alpha_sums[a_type] = logsum2expd(alpha_sums[a_type], lpt);;
                        alpha_sums[a_type] += pt;
                    }

                    // sig_sum[ CMI(pk, a_type, sig_nrow) ] = 
                    //     logsum2expd(sig_sum[ CMI(pk, a_type, sig_nrow) ], lpt);
                    sig_sum[ CMI(pk, a_type, sig_nrow) ] += pt;
                }
            }
        }

        mdl->mp->_alpha_sum[CMI(bc_ix, 0, D)] += alpha_sums[0];
        mdl->mp->_alpha_sum[CMI(bc_ix, 1, D)] += alpha_sums[1];
        /*
        mdl->mp->_alpha_sum[CMI(bc_ix, 0, D)] = 
            logsum2expd(mdl->mp->_alpha_sum[CMI(bc_ix, 0, D)], alpha_sums[0]);
        mdl->mp->_alpha_sum[CMI(bc_ix, 1, D)] = 
            logsum2expd(mdl->mp->_alpha_sum[CMI(bc_ix, 1, D)], alpha_sums[1]);
        */
    }

    // add to tmp sum variables
    // add lambda
    pthread_mutex_lock(&mdl->mp->sum_lock);
    for (i = 0; i < 3; ++i){
        // mdl->mp->_lambda_sum[i] = logsum2expd(mdl->mp->_lambda_sum[i], lambda_sums[i]);
        mdl->mp->_lambda_sum[i] += lambda_sums[i];
    }

    // add rho
    if (bam_data->has_rna){
        for (i = 0; i < rho_ne; ++i){
            // mdl->mp->_rho_sum[i] = logsum2expd(mdl->mp->_rho_sum[i], rho_sum[i]);
            mdl->mp->_rho_sum[i] += rho_sum[i];
        }
    }

    // add sigma
    if (bam_data->has_atac){
        for (i = 0; i < sig_ne; ++i){
            // mdl->mp->_sigma_sum[i] = logsum2expd(mdl->mp->_sigma_sum[i], sig_sum[i]);
            mdl->mp->_sigma_sum[i] += sig_sum[i];
        }
    }
    pthread_mutex_unlock(&mdl->mp->sum_lock);

    free(rho_sum);
    free(sig_sum);

    for (t_i = 0; t_i < 3; ++t_i) free(p_tm_v[t_i]);

    return(0);
}

int mdl_m_pi_fix(mdl_t *mdl){

    uint32_t i, j, M = mdl->mp->M;

    // maximize pi (keep fixed uniform for now)
    f_t pi_upd = 1.0 / M;
    for (i = 0; i < M; ++i) mdl->mp->pi[i] = pi_upd;
    f_t pi_sum = 0;
    for (i = 0; i < M; ++i){
        for (j = i+1; j < M; ++j){
            pi_sum += (mdl->mp->pi[i] * mdl->mp->pi[j]);
        }
    }
    mdl->mp->pi_d_sum = pi_sum;

    return(0);
}

int mdl_pars_add_gamma(mdl_pars_t *gp, float **a, int nv, int ns){
    if (a == NULL)
        return err_msg(-1, 0, "mdl_pars_add_gamma: a is NULL");

    if (gp->M != (uint16_t)ns)
        return err_msg(-1, 0, "mdl_pars_add_gamma: ns and number of samples in a don't match");

    uint16_t ns1 = (uint16_t)ns + 1;

    gp->gamma = calloc(nv * ns1, sizeof(f_t));
    
    int v;
    for (v = 0; v < nv; ++v){
        int n_miss = 0;
        f_t d_sum = 0.0;
        int s;
        for (s = 0; s < ns; ++s){
            // if variant is missing, set to -1
            f_t av = a[v] == NULL ? -1 : a[v][s];
            gp->gamma[CMI(s, v, ns1)] = (f_t)av;
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
        gp->gamma[CMI(s, v, ns1)] = (f_t)ambp;
    }

    return 0;
}

int mdl_m_set_gamma_amb(mdl_t *mdl){

    uint32_t i, M = mdl->mp->M, V = mdl->mp->V, M1 = M + 1;

    // update gamma
    for (i = 0; i < V; ++i){
        uint16_t st, n_miss = 0;
        f_t ap = 0;
        for (st = 0; st < M; ++st){
            f_t gprob = mdl->mp->gamma[CMI(st, i, M1)];
            if (gprob < 0)
                ++n_miss; // if missing, skip
            else
                ap += mdl->mp->pi[st] * mdl->mp->gamma[CMI(st, i, M1)];
        }
        if (n_miss == M)
            mdl->mp->gamma[CMI(st, i, M1)] = -1; // if all missing, set missing
        else
            mdl->mp->gamma[CMI(st, i, M1)] = ap;
    }
    return(0);
}

int mdl_m_est(mdl_t *mdl){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_m_est: mdl is null");

    f_t new_par = 0;
    mdl->mp->_par_diff = 0;
    int lsret;
    uint32_t i, j;

    uint32_t D = mdl->mp->D, G = mdl->mp->G;
    uint32_t rho_nrow = 3 * G;

    // maximize lambda
    // f_t lambda_tot = logsumexpd2(mdl->mp->_lambda_sum, 3);
    f_t lambda_tot = 0;
    for (i = 0; i < 3; ++i) lambda_tot += mdl->mp->_lambda_sum[i];
    assert( (!isnan(lambda_tot)) && (!isinf(lambda_tot)) );
    for (i = 0; i < 3; ++i){
        // new_par = exp(mdl->mp->_lambda_sum[i] - lambda_tot);
        new_par = mdl->mp->_lambda_sum[i] / lambda_tot;
        mdl->mp->_par_diff += fabs(new_par - mdl->mp->lambda[i]);
        mdl->mp->lambda[i] = new_par;
    }

    // maximize pi (keep fixed uniform for now)
    if (mdl_m_pi_fix(mdl) < 0) return(-1);

    // update gamma
    if (mdl_m_set_gamma_amb(mdl) < 0) return(-1);

    // alpha
    f_t a_total = 0, a_sum[2] = {0,0};
    for (i = 0; i < D; ++i){
        int fl;
        fl = bflg_get(mdl->absent_bc, i);
        if (fl){
            mdl->mp->alpha[i] = 0.5;
            continue;
        }

        a_sum[0] = mdl->mp->_alpha_sum[CMI(i, 0, D)];
        a_sum[1] = mdl->mp->_alpha_sum[CMI(i, 1, D)];
        // a_total = logsumexpd2(a_sum, 2);
        a_total = a_sum[0] + a_sum[1];

        if (isfinite(a_total)){
            // new_par = exp(a_sum[0] - a_total);
            new_par = a_sum[0]/a_total;
        } else {
            new_par = 0.5;
        }
        if (new_par <= 0 || new_par >= 1){
            fprintf(stderr, "alpha[%i]=%f\n", i, new_par);
            fprintf(stderr, "alpha[0]=%f, alpha[1]=%f\n", a_sum[0], a_sum[1]);
            return(-1);
        }
        // only calculate delt for unfixed droplets
        fl = bflg_get(mdl->fix_flag, i);
        if (!fl) mdl->mp->_par_diff += fabs(new_par - mdl->mp->alpha[i]);
        if (fl) mdl->mp->alpha[i] = 1.0;
        else mdl->mp->alpha[i] = new_par;
    }

    // rho
    if (mdl->has_rna){
        /*
        f_t rho_tot[2] = {-INFINITY,-INFINITY};
        rho_tot[0] = logsumexpd2(mdl->mp->_rho_sum, rho_nrow);
        rho_tot[1] = logsumexpd2(mdl->mp->_rho_sum + rho_nrow, rho_nrow);
        */
        f_t rho_tot[2] = {0,0};
        for (i = 0; i < rho_nrow; ++i){
            for (j = 0; j < 2; ++j){
                uint32_t rix = CMI(i, j, rho_nrow);
                rho_tot[j] += mdl->mp->_rho_sum[rix];
            }
        }
                
        // check for errors
        if ( (isnan(rho_tot[0])) || (isnan(rho_tot[1])) )
            return err_msg(-1, 0, "mdl_m_est: rho prob. sum is nan: [%f,%f]", 
                    rho_tot[0], rho_tot[1]);
        if (isinf(rho_tot[0]) || isinf(rho_tot[1]))
            return err_msg(-1, 0, "mdl_m_est: rho prob. sum is 0: [%f,%f]", 
                    rho_tot[0], rho_tot[1]);

        for (i = 0; i < rho_nrow; ++i){
            for (j = 0; j < 2; ++j){
                uint32_t rix = CMI(i, j, rho_nrow);
                // new_par = exp(mdl->mp->_rho_sum[rix] - rho_tot[j]);
                new_par = mdl->mp->_rho_sum[rix] / rho_tot[j];
                mdl->mp->_par_diff += fabs(new_par - mdl->mp->rho[rix]);
                mdl->mp->rho[rix] = new_par;
            }
        }
    }

    // sigma
    if (mdl->has_atac){
        uint32_t sig_nrow = 2;
        f_t sig_tot[2];
        /*
        sig_tot[0] = logsumexpd2(mdl->mp->_sigma_sum, sig_nrow);
        sig_tot[1] = logsumexpd2(mdl->mp->_sigma_sum + sig_nrow, sig_nrow);
        */
        for (i = 0; i < 2; ++i){
            sig_tot[i] = mdl->mp->_sigma_sum[CMI(0,i,sig_nrow)] + 
                mdl->mp->_sigma_sum[CMI(1,i,sig_nrow)];
        }
        if (isnan(sig_tot[0]) || isnan(sig_tot[1]))
            return err_msg(-1, 0, "mdl_max_sigma: sigma prob. sum is nan: [%f,%f]", 
                    sig_tot[0], sig_tot[1]);
        if (isinf(sig_tot[0]) || isinf(sig_tot[1]))
            return err_msg(-1, 0, "mdl_max_sigma: sigma prob. sum is inf: [%f,%f]", 
                    sig_tot[0], sig_tot[1]);
        // new_par = exp(mdl->mp->_sigma_sum[ CMI(1, 0, sig_nrow) ] - sig_tot[0]);
        new_par = mdl->mp->_sigma_sum[ CMI(1, 0, sig_nrow) ] / sig_tot[0];
        mdl->mp->_par_diff += fabs(new_par - mdl->mp->sigma[0]);
        mdl->mp->sigma[0] = new_par;

        // new_par = exp(mdl->mp->_sigma_sum[ CMI(1, 1, sig_nrow) ] - sig_tot[1]);
        new_par = mdl->mp->_sigma_sum[ CMI(1, 1, sig_nrow) ] / sig_tot[1];
        mdl->mp->_par_diff += fabs(new_par - mdl->mp->sigma[1]);
        mdl->mp->sigma[1] = new_par;
    }

    return(0);
}

int mdl_delta(mdl_t *mdl, f_t *delta){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_delta: arguments are NULL");

    f_t n_pars = 0.0;
    n_pars += 3.0; // lambda
    n_pars += (f_t)mdl->test_bcs->n; // alpha (for unfixed droplets)
    if (mdl->has_rna) n_pars += 6.0 * (f_t)mdl->mp->G; // rho
    if (mdl->has_atac) n_pars += 2.0; // sigma

    // *delta = exp(log(mdl->mp->_par_diff) - log(n_pars));
    *delta = mdl->mp->_par_diff / n_pars;

    return(0);
}

/*******************************************************************************
 * EM functions
 ******************************************************************************/

int mdl_init_par_dat(mdl_t *mdl, bam_data_t *bam_data){

    if (mdl_m_init(mdl, 1.0) < 0) return(-1);

    int lsret = 0;
    uint32_t i;
    uint32_t G = mdl->mp->G;
    uint32_t D = mdl->mp->D;
    uint32_t rho_nrow = 3 * G;

    uint32_t sig_nrow = 2;

    // set pi as fixed
    if (mdl_m_pi_fix(mdl) < 0) return(-1);

    mdl->mp->lambda[0] = 0.9;
    mdl->mp->lambda[1] = 0.09;
    mdl->mp->lambda[2] = 0.01;

    mdl->mp->tau = TAU;

    khint_t k;
    for (i = 0; i < D; ++i){
        int fl;

        // if no barcode
        fl = bflg_get(mdl->absent_bc, i);
        if (fl == 1) continue;

        char *bc_key = str_map_str(mdl->all_bcs, i);
        k = kh_get(kh_bc_dat, bam_data->bc_data, bc_key);
        if (k == kh_end(bam_data->bc_data)) continue;

        fl = bflg_get(mdl->fix_flag, i);
        uint32_t c_ix = fl ? 0 : 1;

        // bam data
        bc_data_t *bc_dat = kh_val(bam_data->bc_data, k);
        khash_t(khrmn) *mols = bc_dat->rna_mols;
        khash_t(khaf) *frags = bc_dat->atac_frags;

        khint_t k_rna, k_atac;
        if (bam_data->has_rna){
            for (k_rna = kh_begin(mols); k_rna != kh_end(mols); ++k_rna){
                if (!kh_exist(mols, k_rna)) continue;
                rna_mol_t *mol = kh_val(mols, k_rna);

                size_t n_genes = ml_size(&mol->gl);
                if (n_genes > 1) n_genes = 0; // skip multigene UMIs
                if (n_genes == 1){
                    ml_node_t(seq_gene_l) *gn = ml_begin(&mol->gl);
                    seq_gene_t seq_gene = ml_node_val(gn);
                    int32_t g_ix = seq_gene.gene_id;
                    uint8_t spl = seq_gene.splice;
                    uint32_t f_ix = g_ix + (spl * mdl->mp->G);
                    uint32_t rix = CMI(f_ix, c_ix, rho_nrow);
                    // mdl->mp->_rho_sum[rix] = logsum2expd(mdl->mp->_rho_sum[rix], log(1.0));
                    mdl->mp->_rho_sum[rix] += 1.0;
                }
            }
        }
        if (bam_data->has_atac){
            for (k_atac = kh_begin(frags); k_atac != kh_end(frags); ++k_atac){
                if (!kh_exist(frags, k_atac)) continue;
                atac_frag_t *frag = kh_val(frags, k_atac);

                size_t n_peaks = mv_size(&frag->pks);
                uint8_t pk = n_peaks == 0 ? 0 : 1; // 0: outside peak, 1: inside peak

                uint32_t six = CMI(pk, c_ix, sig_nrow);

                // mdl->mp->_sigma_sum[six] = logsum2expd(mdl->mp->_sigma_sum[six], log(1.0));
                mdl->mp->_sigma_sum[six] += 1.0;
            }
        }
    }

    if (mdl_m_est(mdl) < 0) return(-1);

    return(0);
}
int mdl_init_par_uni(mdl_t *mdl){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_init_par_uni: arguments are NULL");

    uint32_t i, j;
    // lambda
    for (i = 0; i < 3; ++i)
        mdl->mp->lambda[i] = 1./3.;

    // set pi as fixed
    if (mdl_m_pi_fix(mdl) < 0) return(-1);

    // set gamma ambient prob.
    if (mdl_m_set_gamma_amb(mdl) < 0) return(-1);

    // alpha
    for (i = 0; i < mdl->mp->D; ++i){
        int fl;
        fl = bflg_get(mdl->absent_bc, i);
        if (fl){
            mdl->mp->alpha[i] = 1.0;
            continue;
        }
        fl = bflg_get(mdl->fix_flag, i);
        if (fl) mdl->mp->alpha[i] = 1.0;
        else mdl->mp->alpha[i] = 0.5;
    }

    // rho
    uint32_t ng = 3 * mdl->mp->G;
    double rho_init = 1. / (double)ng;
    for (i = 0; i < 2; ++i){
        for (j = 0; j < ng; ++j)
            mdl->mp->rho[CMI(j,i,ng)] = rho_init;
    }

    // sigma
    mdl->mp->sigma[0] = 0.5;
    mdl->mp->sigma[1] = 0.5;

    return(0);
}

int mdl_init_par_rand(mdl_t *mdl, unsigned int seed){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_init_par_rand: arguments are NULL");

    srand(seed);

    uint32_t i, j;
    // lambda
    f_t lambda_tot = 0;
    for (i = 0; i < 3; ++i){
        mdl->mp->lambda[i] = frandin();
        lambda_tot += mdl->mp->lambda[i];
    }
    for (i = 0; i < 3; ++i)
        mdl->mp->lambda[i] /= lambda_tot;

    // set pi as fixed
    if (mdl_m_pi_fix(mdl) < 0) return(-1);

    // set gamma ambient prob.
    if (mdl_m_set_gamma_amb(mdl) < 0) return(-1);

    // alpha
    for (i = 0; i < mdl->mp->D; ++i){
        int fl = bflg_get(mdl->absent_bc, i);
        if (fl){
            mdl->mp->alpha[i] = 1.0;
            continue;
        }
        fl = bflg_get(mdl->fix_flag, i);
        if (fl) mdl->mp->alpha[i] = 1.0;
        else mdl->mp->alpha[i] = frandin();
    }

    // rho
    uint32_t ng = 3 * mdl->mp->G;
    for (i = 0; i < 2; ++i){
        f_t rho_tot = 0;
        for (j = 0; j < ng; ++j){
            mdl->mp->rho[CMI(j,i,ng)] = 1.0 + (f_t)(rand() % 10000);
            rho_tot += mdl->mp->rho[CMI(j,i,ng)];
        }
        for (j = 0; j < ng; ++j)
            mdl->mp->rho[CMI(j,i,ng)] /= rho_tot;
    }

    // sigma
    mdl->mp->sigma[0] = frandin();
    mdl->mp->sigma[1] = frandin();

    return(0);
}

typedef struct thrd_args {
    mdl_t *mdl;
    bam_data_t *bam_data;
    int *ixs;
    uint32_t ix_len;
} thrd_args;

void *mdl_thrd_fx(void *arg){
    thrd_args *t = (thrd_args *)arg;
    if (mdl_e_hs(t->mdl, t->bam_data, t->ixs, t->ix_len) < 0) return((void *)-1);
    if (mdl_m_sum(t->mdl, t->bam_data, t->ixs, t->ix_len) < 0) return((void *)-1);
    return((void *)0);
}

int mdl_thrd_call(mdl_t *mdl, bam_data_t *bam_data, int *ixs, uint32_t ix_len){
    if (mdl == NULL || bam_data == NULL)
        return err_msg(-1, 0, "mdl_thrd_call: arguments are NULL");

    uint32_t i;
    if (mdl->threads < 1)
        return err_msg(-1, 0, "mdl_thrd_call: threads=%u is less than 1", mdl->threads);

    uint32_t mt1 = mdl->threads - 1;
    double step = (double)ix_len / (double)mdl->threads;
    step = ceil(step);
    uint32_t stepi = (uint32_t)ROUND_2_INT(step);
    uint32_t step_last = ix_len - ((mt1) * stepi);

    thrd_args *targs = malloc(mdl->threads * sizeof(thrd_args));
    for (i = 0; i < mdl->threads; ++i){
        targs[i].mdl = mdl;
        targs[i].bam_data = bam_data;
        targs[i].ixs = ixs + (stepi * i);
        targs[i].ix_len = i == (mt1) ? step_last : stepi;
    }

    pthread_t *ids = malloc(mdl->threads * sizeof(pthread_t));

    int err;
    for (i = 0; i < mdl->threads; ++i){
        err = pthread_create(ids + i, NULL, mdl_thrd_fx, &targs[i]);
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

int mdl_em(mdl_t *mdl, bam_data_t *bam_data, obj_pars *objs){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_em_step: arguments are NULL");

    f_t q = 0, q_prev = -INFINITY;
    f_t delta, prior = 10;
    uint32_t i, D = mdl->mp->D, ix_len = D;
    int *ixs = malloc(D * sizeof(int));
    for (i = 0; i < ix_len; ++i) ixs[i] = i;
    permute(ixs, D);

    if (objs->verbose) log_msg("running EM");
    for (i = 0; i < mdl->max_iter; ++i){
        mdl_check_pars(mdl);

        if (mdl_m_init(mdl, prior) < 0) return(-1); 

        if (mdl_thrd_call(mdl, bam_data, ixs, ix_len) < 0) return(-1);
        /*
        if (objs->verbose) log_msg("getting conditional probabilities");
        if (mdl_e_hs(mdl, bam_data, ixs, ix_len) < 0) return(-1);

        if (objs->verbose) log_msg("getting maximum likelihood (sum)");
        if (mdl_m_sum(mdl, bam_data, ixs, ix_len) < 0) return(-1);
        */

        if (mdl_m_est(mdl) < 0) return(-1);

        if (mdl_data_llk(mdl, &q) < 0) return(-1);
        if ( (q - q_prev) < 0 ) err_msg(0, 1, "llk decreased: delta=%f", q - q_prev);
        if (objs->verbose) log_msg("Q=%f; delta=%f", q, q - q_prev);
        q_prev = q;

        if (mdl_delta(mdl, &delta) < 0) return(-1);
        log_msg("iteration %u: delta=%f", i+1, delta);
        mdl_check_pars(mdl);
        if (delta < mdl->eps) break;
    }
    if (objs->verbose) log_msg("finished EM");

    free(ixs);

    return(0);
}

int mdl_data_llk(mdl_t *mdl, f_t *llk){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_llk: mdl is null");
    if (mdl->all_bcs == NULL)
        return err_msg(-1, 0, "mdl_llk: all_bcs is null");
    if (mdl->all_bcs->n < 0)
        return err_msg(-1, 0, "mdl_llk: all_bcs n is less than 0");

    uint32_t i, n_bcs = (uint32_t)mdl->all_bcs->n;
    f_t llkt = 0;
    int fl;
    for (i = 0; i < n_bcs; ++i){
        fl = bflg_get(mdl->absent_bc, i);
        if (fl) continue;
        llkt += mdl->c_probs[i].p_x;
        assert(!isnan(llkt));
    }
    *llk = llkt;
    return(0);
}

int mdl_check_pars(mdl_t *mdl){

    f_t lt1 = 1 - feps;
    f_t ut1 = 1 + feps;

    int cond;
    int lsret = 0;
    uint32_t i, j, M1 = mdl->mp->M + 1;
    mdl_pars_t *mp = mdl->mp;

    // check pi
    f_t pi_sum = 0;
    for (i = 1; i < M1; ++i){
        int s1 = mdl->hs_ix[CMI(i, 1, mdl->_nrow_hs)];
        pi_sum += mdl->mp->pi[s1];
    }
    assert(pi_sum >= lt1 && pi_sum <= ut1);

    pi_sum = 0;
    for (; i < mdl->_nrow_hs; ++i){
        int s1 = mdl->hs_ix[CMI(i, 1, mdl->_nrow_hs)];
        int s2 = mdl->hs_ix[CMI(i, 2, mdl->_nrow_hs)];
        pi_sum += (mdl->mp->pi[s1] * mdl->mp->pi[s2]) / mdl->mp->pi_d_sum;
    }
    assert(pi_sum >= lt1 && pi_sum <= ut1);

    // check c_probs
    for (i = 0; i < mdl->mp->D; ++i){
        int fl = bflg_get(mdl->fix_flag, i);
        uint32_t n_hs = fl == 1 ? 1 : mdl->_nrow_hs;

        // f_t c_sum = logsumexpd2(mdl->c_probs[i].cp_hs, n_hs);
        // c_sum = exp(c_sum);
        f_t c_sum = 0;
        for (j = 0; j < n_hs; ++j){
            c_sum += mdl->c_probs[i].cp_hs[j];
        }
        cond = c_sum >= lt1 && c_sum <= ut1;
        if (!cond){
            fprintf(stderr, "c_prob[%u]=", i);
            for (j = 0; j < n_hs; ++j)
                fprintf(stderr, ",%f", mdl->c_probs[i].cp_hs[j]);
        }
        assert(cond);
    }

    // check lambda
    f_t lambda_tot = 0, pi_tot = 0, rho_tot[2] = {0,0};
    for (i = 0; i < 3; ++i)
        lambda_tot += mp->lambda[i];
    cond = lambda_tot >= lt1 && lambda_tot <= ut1;
    if (!cond){
        for (i = 0; i < 3; ++i)
            fprintf(stderr, "lambda[%u] = %f\n", i, mp->lambda[i]);
    }
    assert(lambda_tot >= lt1 && lambda_tot <= ut1);

    // check pi
    for (i = 0; i < mp->M; ++i)
        pi_tot += mp->pi[i];
    assert(pi_tot >= lt1 && pi_tot <= ut1);

    // alpha
    for (i = 0; i < mp->D; ++i){
        cond = mp->alpha[i] >= 0 && mp->alpha[i] <= 1;
        if (!cond)
            fprintf(stderr, "alpha[%u]=%f\n", i, mp->alpha[i]);
        assert(mp->alpha[i] >= 0 && mp->alpha[i] <= 1);
    }

    if (mdl->has_rna){
        uint32_t ng = mp->G * 3;
        for (i = 0; i < 2; ++i){
            for (j = 0; j < ng; ++j){
                rho_tot[i] += mp->rho[CMI(j, i, ng)];
                assert( mp->rho[CMI(j, i, ng)] >= 0 && mp->rho[CMI(j, i, ng)] <= 1 );
            }
        }
        for (i = 0; i < 2; ++i){
            cond = rho_tot[i] >= lt1 && rho_tot[i] <= ut1;
            if (!cond)
                fprintf(stderr, "rho_tot[%u] = %f\n", i, rho_tot[i]);
            assert(cond);
        }
    }

    if (mdl->has_atac){
        for (i = 0; i < 2; ++i){
            cond = mp->sigma[i] >= 0 && mp->sigma[i] <= 1;
            if (!cond)
                fprintf(stderr, "sigma[%u] = %f\n", i, mp->sigma[i]);
            assert(cond);
        }
    }

    for (i = 0; i < mp->V; ++i)
        assert(mp->gamma[CMI(mp->M, i, M1)] >= 0 && mp->gamma[CMI(mp->M, i, M1)] <= 1);

    return(0);
}

/*******************************************************************************
 ******************************************************************************/

int mdl_get_llk(mdl_t *mdl){
    if (mdl == NULL)
        return 0;

    mdl_fit_t *mf = mdl->mf;

    int lsret = 0;
    uint32_t i, ncomb = mdl->_nrow_hs, M = mdl->samples->n;
    uint32_t n_bc = (uint32_t)mdl->test_bcs->n;
    mf->bc_llks = calloc(ncomb * n_bc, sizeof(f_t));
    mf->best_sng_llk = calloc(n_bc, sizeof(f_t));
    mf->sec_sng_llk = calloc(n_bc, sizeof(f_t));
    mf->best_dbl_llk = calloc(n_bc, sizeof(f_t));
    mf->best_sng_ix = calloc(n_bc, sizeof(uint32_t));
    mf->sec_sng_ix = calloc(n_bc, sizeof(uint32_t));
    mf->best_dbl_ix = calloc(n_bc, sizeof(uint32_t));
    mf->pp = calloc(n_bc * 3, sizeof(f_t));

    f_t lsum[3];

    uint32_t bc_i;
    for (bc_i = 0; bc_i < n_bc; ++bc_i){
        char *bc = str_map_str(mdl->test_bcs, bc_i);
        if (bc == NULL)
            return err_msg(-1, 0, "mdl_get_llk: test bc %u not found", bc_i);
        int all_bc_ix = str_map_ix(mdl->all_bcs, bc);
        if (all_bc_ix < 0)
            return err_msg(-1, 0, "mdl_get_llk: test bc %s not found in all bcs", bc);

        c_probs_t cps = mdl->c_probs[all_bc_ix];

        for (i = 0; i < 3; ++i)
            lsum[i] = -INFINITY;

        // fill in bc_llks, get PP
        uint32_t hs_ix = 0;
        for (hs_ix = 0; hs_ix < ncomb; ++hs_ix){
            int h = mdl->hs_ix[CMI(hs_ix, 0, mdl->_nrow_hs)];
            f_t ll = mdl->c_probs[all_bc_ix].lp_hs[hs_ix];
            mf->bc_llks[CMI(bc_i, hs_ix, n_bc)] = ll;
            lsum[h] = logsum2expd(lsum[h], ll);
        }
        f_t lsumall = logsumexpd2(lsum, 3);
        for (i = 0; i < 3; ++i)
            mf->pp[CMI(bc_i, i, n_bc)] = exp(lsum[i] - lsumall);

        f_t sllk[2] = {-INFINITY, -INFINITY};
        uint32_t si[2] = {-1,-1};
        for (hs_ix = 1; hs_ix < M + 1; ++hs_ix){
            if (mf->bc_llks[CMI(bc_i, hs_ix, n_bc)] > sllk[0]){
                sllk[0] = mf->bc_llks[CMI(bc_i, hs_ix, n_bc)];
                si[0] = hs_ix;
                continue;
            }
            if (mf->bc_llks[CMI(bc_i, hs_ix, n_bc)] > sllk[1]){
                sllk[1] = mf->bc_llks[CMI(bc_i, hs_ix, n_bc)];
                si[1] = hs_ix;
            }
        }
        f_t dllk[2] = {-INFINITY, -INFINITY};
        uint32_t di[2] = {-1,-1};
        for (; hs_ix < ncomb; ++hs_ix){
            if (mf->bc_llks[CMI(bc_i, hs_ix, n_bc)] > dllk[0]){
                dllk[0] = mf->bc_llks[CMI(bc_i, hs_ix, n_bc)];
                di[0] = hs_ix;
                continue;
            }
            if (mf->bc_llks[CMI(bc_i, hs_ix, n_bc)] > dllk[1]){
                dllk[1] = mf->bc_llks[CMI(bc_i, hs_ix, n_bc)];
                di[1] = hs_ix;
            }
        }
        mf->best_sng_llk[bc_i] = sllk[0];
        mf->sec_sng_llk[bc_i] = sllk[1];
        mf->best_sng_ix[bc_i] = si[0];
        mf->sec_sng_ix[bc_i] = si[1];
        mf->best_dbl_llk[bc_i] = dllk[0];
        mf->best_dbl_ix[bc_i] = di[0];
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

    if (objs->verbose) log_msg("initializing model");

    mdl_t *mdl = mdl_alloc();
    if (mdl == NULL)
        return -1;

    // sets all bc, flt_bcs, and test_bcs
    if (mdl_set_bcs(mdl, bam_dat, objs) < 0)
        return -1;

    if (mdl_set_samples(mdl, objs) < 0)
        return -1;

    if (mdl_init_all(mdl, bam_dat, objs) < 0)
        return -1;

    if (mdl_em(mdl, bam_dat, objs) < 0)
        return -1;

    if (objs->verbose) log_msg("getting barcode likelihoods");
    if (mdl_get_llk(mdl) < 0)
        return -1;

    if (objs->verbose) log_msg("writing out likelihoods");
    if (write_llk(mdl, objs->out_fn) < 0)
        return(-1);

    if (objs->verbose) log_msg("writing out parameters");
    if (write_lambda(mdl, objs->out_fn) < 0)
        return -1;
    if (write_alpha(mdl, objs->out_fn) < 0)
        return -1;
    if (bam_dat->has_rna && write_rho(mdl, objs, objs->out_fn) < 0)
        return -1;
    if (bam_dat->has_atac && write_sigma(mdl, objs, objs->out_fn) < 0)
        return -1;

    if (write_res(mdl, bam_dat, objs->out_fn) < 0)
        return -1;

    if (objs->verbose) log_msg("writing out summary results");
    mdl_dstry(mdl);

    if (objs->verbose) log_msg("model fit finished");

    return(0);
}

int write_lambda(mdl_t *mdl, char *fn){
    if (mdl == NULL || fn == NULL)
        return err_msg(-1, 0, "write_lambda: arguments are NULL");

    if (mdl->mp == NULL || mdl->mp->lambda == NULL)
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

int write_alpha(mdl_t *mdl, char *fn){
    if (mdl == NULL)
        return err_msg(-1, 0, "write_alpha: arguments are NULL");

    if (mdl->mp == NULL || mdl->mp->alpha == NULL)
        return err_msg(-1, 0, "write_alpha: model hasn't been initialized");

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    char *alpha_fn = ".alpha.txt.gz";
    char *out_alpha_fn = strcat2((const char*)fn, (const char*)alpha_fn);

    // row names
    char **bc_row_names = str_map_ca(mdl->test_bcs);
    int bc_n = mdl->test_bcs->n;

    // col names
    char **col_names = malloc(sizeof(char *));
    col_names[0] = strdup("Alpha");
    if (col_names[0] == NULL)
        return err_msg(-1, 0, "write_alpha: %s", strerror(errno));

    // alpha array for test barcodes
    f_t *al = malloc(bc_n * sizeof(f_t));
    int i;
    for (i = 0; i < bc_n; ++i){
        char *bc = str_map_str(mdl->test_bcs, i);
        int bc_ix = str_map_ix(mdl->all_bcs, bc);
        al[i] = mdl->mp->alpha[bc_ix];
    }

    // write matrix
    ret = write_matrix_double(out_alpha_fn, al, NULL, NULL, NULL, 
            bc_row_names, bc_n, col_names, 1, 
            delim, nl, decp);
    free(al);
    if (ret < 0)
        return err_msg(-1, 0, "write_alpha: failed to write matrix to file");

    free(col_names[0]);
    free(col_names);
    for (i = 0; i < bc_n; ++i)
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
    uint32_t gs, g, G = mdl->mp->G, GA = 3 * G;
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
    char **col_names = malloc(2 * sizeof(char*));
    col_names[0] = strdup("Ambient");
    if (col_names[0] == NULL)
        return err_msg(-1, 0, "write_rho: %s", strerror(errno));
    col_names[1] = strdup("Cell");
    if (col_names[1] == NULL)
        return err_msg(-1, 0, "write_rho: %s", strerror(errno));

    // write matrix
    ret = write_matrix_double(out_rho_fn, mdl->mp->rho, NULL, NULL, NULL, 
            gene_row_names, GA, col_names, 2, 
            delim, nl, decp);
    if (ret < 0)
        return err_msg(-1, 0, "write_rho: failed to write matrix to file");

    free(col_names[0]);
    free(col_names[1]);
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

    if (mdl->mp == NULL || mdl->mp->sigma == NULL)
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
    int i, ret = 0;
    char *llk_fn = ".llk.txt.gz";
    const char ssep[] = "-";

    char *out_llk_fn = strcat2((const char*)fn, (const char*)llk_fn);
    if (out_llk_fn == NULL)
        return err_msg(-1, 0, "write_llk: %s", strerror(errno));

    // row names
    log_msg("getting row names");
    char **bc_row_names = str_map_ca(mdl->test_bcs);
    if (bc_row_names == NULL)
        return err_msg(-1, 0, "write_llk: failed to write to file");

    int bc_n = mdl->test_bcs->n;

    // col names
    log_msg("getting column names");
    char **col_names = malloc(mdl->_nrow_hs * sizeof(char *));
    if (col_names == NULL)
        return err_msg(-1, 0, "write_llk: %s", strerror(errno));

    uint32_t hs_ix;
    for (hs_ix = 0; hs_ix < mdl->_nrow_hs; ++hs_ix){
        int hd, s1_ix, s2_ix, t_ix[3], t_n;
        if (mdl_get_hst(mdl, hs_ix, &hd, &s1_ix, &s2_ix, t_ix, &t_n) < 0) return(-1);
        char *s1, *s2;
        
        if (hd == 0){
            if ( (col_names[hs_ix] = strdup("Empty")) == NULL )
                return err_msg(-1, 0, "write_llk: %s", strerror(errno));
        } else if (hd == 1) {
            s1 = str_map_str(mdl->samples, s1_ix);
            if (s1 == NULL)
                return err_msg(-1, 0, "write_llk: failed to get sample for index %i", s1_ix);
            if ( (col_names[hs_ix] = strdup(s1)) == NULL )
                return err_msg(-1, 0, "write_llk: %s", strerror(errno));
        } else {
            s1 = str_map_str(mdl->samples, s1_ix);
            if (s1 == NULL)
                return err_msg(-1, 0, "write_llk: failed to get sample for index %i", s1_ix);
            s2 = str_map_str(mdl->samples, s2_ix);
            if (s2 == NULL)
                return err_msg(-1, 0, "write_llk: failed to get sample for index %i", s2_ix);
            const char *sa[2] = {s1, s2};
            if ( (col_names[hs_ix] = cat_strs(sa, 2, ssep)) == NULL )
                return err_msg(-1, 0, "write_llk: failed to cat samples '%s' and '%s'", s1, s2);
        }
    }

    // write matrix
    log_msg("writing matrix");
    ret = write_matrix_double(out_llk_fn, mdl->mf->bc_llks, NULL, NULL, NULL, 
            bc_row_names, bc_n, col_names, mdl->_nrow_hs, 
            delim, nl, decp);
    if (ret < 0)
        return err_msg(-1, 0, "write_llk: failed to write matrix to file");

    log_msg("freeing memory");
    for (hs_ix = 0; hs_ix < mdl->_nrow_hs; ++hs_ix) free(col_names[hs_ix]);
    free(col_names);
    for (i = 0; i < bc_n; ++i) free(bc_row_names[i]);
    free(bc_row_names);
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
    FILE *fp = fopen(out_s_fn, "w");
    if (fp == NULL){
        err_msg(-1, 0, "write_samples: Could not open file %s", out_s_fn);
        free(out_s_fn);
        return -1;
    }

    int hs_ix, fret, hd, s1, s2, t_ix[3], t_n, nhs = (int)mdl->_nrow_hs;
    for (hs_ix = 1; hs_ix < nhs; ++hs_ix){
        if (mdl_get_hst(mdl, hs_ix, &hd, &s1, &s2, t_ix, &t_n) < 0) return(-1);

        // if singlet
        if (hd == 1){ // if singlet
            fret = fputs(mdl->samples->strs[s1], fp);
            fret = fputc(delim, fp);
            fret = fputs(mdl->samples->strs[s1], fp);
            fret = fputc(nl, fp);
        } else { // if doublet
            fret = fputs(mdl->samples->strs[s1], fp);
            fret = fputc(delim, fp);
            fret = fputs(mdl->samples->strs[s2], fp);
            fret = fputc(nl, fp);
        }
    }

    free(out_s_fn);

    // close file
    if ( (fret = fclose(fp)) != 0 )
        return err_msg(-1, 0, "write_samples: could not close file %s: %s", 
                out_s_fn, strerror(errno));

    return 0;
}

int write_res(mdl_t *mdl, bam_data_t *bam_dat, char *fn){
    char *r_fn = ".summary.txt";
    mdl_fit_t *mf = mdl->mf;
    str_map *samples = mdl->samples;
    str_map *bcs = mdl->test_bcs;
    int bcs_n = bcs->n;
    unsigned int decp = 8;
    int fret;
    int delim = '\t';
    int nl = '\n';

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "write_res: failed to create output directory for %s", fn);

    char *out_r_fn = strcat2((const char*)fn, (const char*)r_fn);
    FILE *fp = fopen(out_r_fn, "w");
    if (fp == NULL){
        err_msg(-1, 0, "write_res: Could not open file %s", out_r_fn);
        free(out_r_fn);
        return -1;
    }

    char hdr[] = "Barcode\tn_rna_molecules\tn_atac_molecules\tn_features\tn_peaks"
        "\tn_rna_variants\tn_atac_variants\tatac_pct_mt\trna_pct_mt\tFRIP\tAlpha"
        "\tSinglet_best\tSinglet_best_llk"
        "\tSinglet_scnd\tSinglet_scnd_llk"
        "\tDoublet_sample1\tDoublet_sample2\tDoublet_llk"
        "\tPP0\tPP1\tPP2\n";

    fputs(hdr, fp);

    size_t buf_size = decp + 1000;
    char *pstr = (char *)malloc(buf_size * sizeof(char));
    if (pstr == NULL)
        return err_msg(-1, 0, "write_res: %s", strerror(errno));

    char *s1, *s2;
    int s_ix1, s_ix2, bc_i, pstr_len;
    for (bc_i = 0; bc_i < bcs_n; ++bc_i){
        char *bc_name = str_map_str(bcs, bc_i);
        fret = fputs(bc_name, fp);

        khint_t k_bc = kh_get(kh_bc_dat, bam_dat->bc_data, bc_name);
        if (k_bc == kh_end(bam_dat->bc_data))
            return err_msg(-1, 0, "write_res: barcode %s not found", bc_name);

        bc_data_t *bc_data = kh_val(bam_dat->bc_data, k_bc);
        bc_stats_t *bcc = bc_data->bc_stats;

        int all_ix = str_map_ix(mdl->all_bcs, bc_name);
        assert(all_ix >= 0);
        f_t al = mdl->mp->alpha[all_ix];

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

        fputc(delim, fp);
        double2str_in((double)al, &pstr, &buf_size, 4);
        fputs(pstr, fp);

        // write best singlet
        fret = fputc(delim, fp);

        s_ix1 = mdl->hs_ix[CMI(mf->best_sng_ix[bc_i], 1, mdl->_nrow_hs)];
        s1 = str_map_str(samples, s_ix1);
        fret = fputs(s1, fp);

        // write best singlet llk
        double par = mf->best_sng_llk[bc_i];
        if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
            return err_msg(-1, 0, "write_llk: failed to convert %f to string", par);

        fret = fputc(delim, fp);
        fret = fputs(pstr, fp);

        // write 2nd singlet
        fret = fputc(delim, fp);

        s_ix1 = mdl->hs_ix[CMI(mf->sec_sng_ix[bc_i], 1, mdl->_nrow_hs)];
        s1 = str_map_str(samples, s_ix1);
        fret = fputs(s1, fp);

        // write 2nd llk
        par = mf->sec_sng_llk[bc_i];
        if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
            return err_msg(-1, 0, "write_llk: failed to convert %f to string", par);

        fret = fputc(delim, fp);
        fret = fputs(pstr, fp);

        // write best doublet
        s_ix1 = mdl->hs_ix[CMI(mf->best_dbl_ix[bc_i], 1, mdl->_nrow_hs)];
        s_ix2 = mdl->hs_ix[CMI(mf->best_dbl_ix[bc_i], 2, mdl->_nrow_hs)];
        s1 = str_map_str(samples, s_ix1);
        s2 = str_map_str(samples, s_ix2);
        fret = fputc(delim, fp);
        fret = fputs(s1, fp);
        fret = fputc(delim, fp);
        fret = fputs(s2, fp);

        // write best doublet llk
        par = mf->best_dbl_llk[bc_i];
        if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
            return err_msg(-1, 0, "write_llk: failed to convert %f to string", par);

        fret = fputc(delim, fp);
        fret = fputs(pstr, fp);

        // write log likelihood ratio of sng to doublet
        int i;
        for (i = 0; i < 3; ++i){
            par = mf->pp[CMI(bc_i, i, bcs_n)];
            if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
            return err_msg(-1, 0, "write_llk: failed to convert %f to string", par);
            fret = fputc(delim, fp);
            fret = fputs(pstr, fp);
        }

        fret = fputc(nl, fp);
    }

    free(pstr);
    free(out_r_fn);

    // close file
    if ( (fret = fclose(fp)) != 0 )
        return err_msg(-1, 0, "write_samples: could not close file %s: %s", 
                out_r_fn, strerror(errno));

    return 0;
}

