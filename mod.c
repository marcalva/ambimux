
#include "mod.h"
#include "str_util.h"
#include "rna_data.h"
#include "clopts.h"
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>
#include <math.h>

void ix2indices(int ix, int M, int *s1, int *s2){
    if (ix < M){
        *s1 = ix;
        *s2 = -1;
        return;
    }
    int rix = 0;
    int cmax = M - 1;

    ix = ix - M;

    while (ix >= 0){
        if (ix >= 0 && ix < cmax){
            *s1 = rix;
            *s2 = ix + (M-cmax);
            break;
        }
        ix -= cmax;
        cmax--;
        rix++;
    }
}

mdl_pars_t *init_mdl_pars(){
    mdl_pars_t *gp = calloc(1, sizeof(mdl_pars_t));
    return gp;
}

int set_mdl_pars(mdl_pars_t *gp, uint32_t D, uint32_t G, uint32_t V, uint16_t M){
    if (gp == NULL)
        return err_msg(-1, 0, "set_mdl_pars: arguments are NULL");

    gp->D = D;
    gp->G = G;
    gp->V = V;
    gp->M = M;

    gp->alpha = NULL;
    gp->gamma = NULL;
    gp->beta = NULL;

    return 0;
}

void destroy_mdl_pars(mdl_pars_t *gp){
    if (gp == NULL) return;

    if (gp->alpha) free(gp->alpha);
    if (gp->gamma) free(gp->gamma);
    if (gp->beta) free(gp->beta);
    
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
    free(mf->llr);
    free(mf);
}

mdl_t *mdl_alloc(){
    mdl_t *mdl = (mdl_t *)calloc(1, sizeof(mdl_t));
    if (mdl == NULL){
        err_msg(-1, 0, "mdl_alloc: %s", strerror(errno));
        return(NULL);
    }

    mdl->all_bcs = NULL;
    mdl->flt_bcs = NULL;
    mdl->samples = NULL;

    mdl->bcss = NULL;

    mdl->mp = init_mdl_pars();
    mdl->mf = init_mdl_fit();
    if (mdl->mp == NULL || mdl->mf == NULL)
        return(NULL);

    return(mdl);
}

void mdl_dstry(mdl_t *m){
    if (m == NULL) return;

    if (m->bcss){
        bcs_stats_free(m->bcss);
        free(m->bcss);
    }

    if (m->all_bcs) destroy_str_map(m->all_bcs);
    if (m->flt_bcs) destroy_str_map(m->flt_bcs);
    if (m->test_bcs) destroy_str_map(m->test_bcs);
    if (m->samples) destroy_str_map(m->samples);

    destroy_mdl_pars(m->mp);
    destroy_mdl_fit(m->mf);
    free(m);
}

int mdl_set_bcs(mdl_t *mdl, bam_data_t *bam_dat, obj_pars *objs){
    if (mdl == NULL || bam_dat == NULL || objs == NULL)
        return err_msg(-1, 0, "mdl_set_bcs: arguments are null");

    khint_t k;

    // get barcode stats
    if ( (mdl->bcss = bc_stats_alloc()) == NULL )
        return(-1);
    if (bam_data_fill_stats(bam_dat) < 0)
        return(-1);

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

        // if wl_bc not given, set to all varcodes present in data
        for (k = kh_begin(bam_dat->bc_data); k != kh_end(bam_dat->bc_data); ++k){
            if (!kh_exist(bam_dat->bc_data, k))
                continue;

            char *bc_key = kh_key(bam_dat->bc_data, k);

            if (add2str_map(mdl->all_bcs, bc_key, &found) < 0)
                return(-1);
            if (found == 1)
                return err_msg(-1, 0, "mdl_set_bcs: %s wl barcode duplicate "
                        "found, there is a bug", bc_key);
        }
    }

    // set the filtered barcodes to calculate alpha from
    if (objs->flt_bcs != NULL){
        mdl->flt_bcs = str_map_copy(objs->flt_bcs);
        if (mdl->flt_bcs == NULL)
            return err_msg(-1, 0, "mdl_set_bcs: failed to copy flt_bcs");
    } else {
        if (objs->flt_n_bcs >= kh_size(bam_dat->bc_data)){
            return err_msg(-1, 0, "mdl_set_bcs: the --flt-n barcodes (%u) "
                    "is greater than or equal to "
                    "the number of barcodes present in the data (%u). "
                    "Set to a value less than this.", 
                    objs->flt_n_bcs, kh_size(bam_dat->bc_data));
        }

        if ( (mdl->flt_bcs = init_str_map()) == NULL )
            return(-1);

        int ret;
        uint32_t min_c = bam_data_count_of_n(bam_dat, objs->flt_n_bcs, &ret);
        if (ret < 0)
            return(-1);

        for (k = kh_begin(bam_dat->bc_data); k != kh_end(bam_dat->bc_data); ++k){
            if (!kh_exist(bam_dat->bc_data, k))
                continue;

            char *bc_key = kh_key(bam_dat->bc_data, k);
            bc_data_t *bc_data = kh_val(bam_dat->bc_data, k);

            if (bc_data->bc_stats->counts < min_c)
                continue;

            if (add2str_map(mdl->flt_bcs, bc_key, &found) < 0)
                return(-1);
            if (found == 1)
                return err_msg(-1, 0, "mdl_set_bcs: %s flt barcode duplicate "
                        "found, there is a bug", bc_key);
        }
    }

    // set the output barcodes to calculate llk for
    mdl->test_bcs = init_str_map();
    uint32_t out_mins = (uint32_t)objs->out_min;
    for (k = kh_begin(bam_dat->bc_data); k != kh_end(bam_dat->bc_data); ++k){
        if (!kh_exist(bam_dat->bc_data, k))
            continue;
        char *bc_key = kh_key(bam_dat->bc_data, k);
        bc_data_t *bc_data = kh_val(bam_dat->bc_data, k);
        bc_counts *c = bc_data->bc_stats;
        if (c->rna_counts < out_mins && c->atac_counts < out_mins) 
            continue;
        if (add2str_map(mdl->test_bcs, bc_key, &found) < 0)
            return(-1);
        if (found == 1)
            return err_msg(-1, 0, "mdl_set_bcs: %s already found", c->bc);
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

int mdl_pars_est_alpha(mdl_pars_t *gp, bam_data_t *bam_data, str_map *flt_bcs, 
        double smooth){
    if (gp == NULL || bam_data == NULL || flt_bcs == NULL)
        return err_msg(-1, 0, "mdl_pars_est_alpha: arguments are null");

    // if RNA wasn't provided, the model only contains variants from ATAC.
    if (bam_data->has_rna == 0)
        return(0);

    uint32_t gs = gp->G * N_SPL;
    uint32_t ne = gs * 2;
    gp->alpha = calloc(ne, sizeof(double));
    if (gp->alpha == NULL)
        return err_msg(-1, 0, "mdl_pars_est_alpha: %s", strerror(errno));

    double ttl_umi[2] = {0,0};

    // add smoothing prior to alpha and ttl_umi
    int i, j;
    for (i = 0; i < 2; ++i){
        for (j = 0; j < gs; ++j){
            gp->alpha[j + (i * gs)] = smooth;
            ttl_umi[i] += smooth;
        }
    }

    // add RNA gene counts to alpha and ttl_umi
    khash_t(kh_bc_dat) *bc_hash = bam_data->bc_data;
    khint_t k;
    for (k = kh_begin(bc_hash); k != kh_end(bc_hash); ++k){
        if (!kh_exist(bc_hash, k)) continue;
        bc_data_t *bc_data = kh_val(bc_hash, k);
        char *bc_key = kh_key(bc_hash, k);

        int bc_o = C_CELL; // 0 if contamination, 1 if cellular
        if (str_map_ix(flt_bcs, bc_key) < 0)
            bc_o = C_AMBN;

        khint_t l;
        khash_t(khrmn) *mols = bc_data->rna_mols;
        for (l = kh_begin(mols); l != kh_end(mols); ++l){
            if (!kh_exist(mols, l)) continue;
            rna_mol_t *u = kh_val(mols, l);
            if (u == NULL){
                fprintf(stderr, "UMI is null\n");
                continue;
            }
            // TODO: remove UMIs that map to multiple genes
            seq_gene_t *gene;
            for (gene = u->genes.head; gene != NULL; gene = gene->next){
                int32_t gid = gene->gene_id;
                uint8_t sid = gene->splice;
                uint32_t ix = gid + (sid * gp->G) + (bc_o * gs);
                gp->alpha[ix]++;
                ttl_umi[bc_o]++;
            }
        }
    }

    // get prob by dividing by total
    double ptot[2] = {0,0};
    for (i = 0; i < 2; ++i){
        for (j = 0; j < gs; ++j){
            gp->alpha[j + (i * gs)] /= ttl_umi[i];
            ptot[i] += gp->alpha[j + (i * gs)];
        }
    }

    // bug checking
    if (ptot[0] < 0.99 || ptot[0] > 1.01 || ptot[1] < 0.99 || ptot[1] > 1.01)
        return err_msg(-1, 0, "mdl_pars_est_alpha: totals must be [1,1], not [%f, %f]", 
                ptot[0], ptot[1]);

    return(0);
}

int mdl_pars_add_gamma(mdl_pars_t *gp, float **a, int nv, int ns){
    if (a == NULL)
        return err_msg(-1, 0, "mdl_pars_add_gamma: a is NULL");

    if (gp->M != (uint16_t)ns)
        return err_msg(-1, 0, "mdl_pars_add_gamma: ns and number of samples in a don't match");

    uint16_t ns1 = (uint16_t)ns + 1;

    gp->gamma = calloc(nv * ns1, sizeof(float));
    
    int v;
    for (v = 0; v < nv; ++v){
        int n_miss = 0;
        float d_sum = 0.0;
        int s;
        for (s = 0; s < ns; ++s){
            // if variant is missing, set to -1
            float av = a[v] == NULL ? -1 : a[v][s];
            gp->gamma[CMI(s, v, ns1)] = av;
            if (av >= 0)
                d_sum += av;
            else
                n_miss++;
        }
        float ambp;
        if (ns == n_miss){
            ambp = -1;
        } else {
            ambp = d_sum / (float)(ns - n_miss);
        }
        gp->gamma[CMI(s, v, ns1)] = ambp;
    }

    return 0;
}

int mdl_pars_set_beta(mdl_pars_t *gp){
    if (gp == NULL || gp->M == 0)
        return err_msg(-1, 0, "mdl_pars_set_beta: arguments are NULL");

    uint16_t M = gp->M;
    double p = 2.0 / (M + M*M);
    int tot = M + (M * (M-1) / 2);
    gp->beta = calloc(tot, sizeof(double));
    int i;
    for (i = 0; i < tot; ++i){
        gp->beta[i] = p;
    }
    return 0;
}

int mdl_pars_set_tau(mdl_pars_t *gp, double tau){
    if (gp == NULL)
        return err_msg(-1, 0, "mdl_pars_set_tau: arguments are NULL");

    if (tau < 0 || tau > 1)
        return err_msg(-1, 0, "mdl_pars_set_tau: tau must be between [0,1]");
    gp->tau = tau;
    return 0;
}

int mdl_pars_set_eps(mdl_pars_t *gp, double eps){
    if (gp == NULL)
        return err_msg(-1, 0, "mdl_pars_set_eps: arguments are NULL");

    if (eps < 0 || eps > 1)
        return err_msg(-1, 0, "set_eps: eps must be between [0,1]");
    gp->eps = eps;
    return 1;
}

int p_brnli(int x, double p, double *prob){
    if (x != 0 && x != 1)
        return err_msg(-1, 0, "p_brnli: x must be 0 or 1, not %i", x);

    if (p < 0 || p > 1)
        return err_msg(-1, 0, "p_brnli: p must be between [0,1], not %f", p);

    if (x == 0) *prob = 1 - p;
    else *prob = p;

    return 0;
}

int p_b_gce(int allele, int v, mdl_pars_t *gp, int s, int e, int c, 
        double *prob){
    if (c != C_AMBN && c != C_CELL)
        return err_msg(-1, 0, "p_b_gce: c must be 0 or 1");
    if (e != 0 && e != 1)
        return err_msg(-1, 0, "p_b_gce: e must be 0 or 1");
    if (allele != 0 && allele != 1)
        return err_msg(-1, 0, "p_b_gce: allele must be 0 or 1");
    if (v < 0 || v >= gp->V)
        return err_msg(-1, 0, "p_b_gce: v must be between 0 and %i", gp->V);

    if (e == 1){
        *prob = 0.5;
        return 0;
    }

    double gprb[2] = {0,0};
    if (c == C_AMBN){
        gprb[0] = gp->gamma[CMI((gp->M), (v), (gp->M+1))];
        gprb[1] = gp->gamma[CMI((gp->M), (v), (gp->M+1))];
    } else {
        // if singlet
        if (s < gp->M){ // if singlet
            gprb[0] = gp->gamma[CMI((s), (v), (gp->M+1))];
            gprb[1] = gp->gamma[CMI((s), (v), (gp->M+1))];
        } else { // if doublet
            int s1 = 0, s2 = 0;
            ix2indices(s, gp->M, &s1, &s2);
            gprb[0] = gp->gamma[CMI((s1), (v), (gp->M+1))];
            gprb[1] = gp->gamma[CMI((s2), (v), (gp->M+1))];
        }
    }

    // if any missing, return 1
    int pret;
    if (gprb[0] == -1 || gprb[1] == -1){
        *prob = 1;
        pret = 0;
    } else {
        double p = 0.5 * (gprb[0] + gprb[1]);
        pret = p_brnli(allele, p, prob);
    }

    return pret;
}

int p_e(int e, double eps, double *prob){
    if (e != 0 && e != 1)
        return err_msg(-1, 0, "p_e: e must be 0 or 1");
    if (eps < 0 || eps > 1)
        return err_msg(-1, 0, "p_e: eps must be between [0,1]");

    if (e == 0) *prob = 1 - eps;
    else *prob = eps;

    return 0;
}

int p_f_c(int f, int spl, mdl_pars_t *gp, int c, double *prob){
    if (f < 0 || f >= gp->G)
        return err_msg(-1, 0, "p_f_c: f must be between 0 and %i", gp->G);
    if (spl != SPLICE && spl != UNSPLICE && spl != AMBIG)
        return err_msg(-1, 0, "p_f_c: spl must be 0, 1, or 2");
    if (c != C_AMBN && c != C_CELL)
        return err_msg(-1, 0, "p_f_c: c must be 0 or 1");

    uint32_t gs = gp->G * N_SPL;
    uint32_t ix = f + (spl * gp->G) + (c * gs);
    *prob = gp->alpha[ix];
    return 0;
}

int p_c(int c, double tau, double *prob){
    if (c != C_AMBN && c != C_CELL)
        return err_msg(-1, 0, "p_c: c must be 0 or 1");
    
    if (c == C_AMBN) *prob = tau;
    else *prob = 1 - tau;
    
    return 0;
}

int mdl_llk(mdl_t *mdl, bam_data_t *bam_dat, int verbose){
    if (mdl == NULL)
        return err_msg(-1, 0, "mdl_llk: mdl is NULL");
    if (bam_dat->bc_data == NULL)
        return err_msg(-1, 0, "mdl_llk: bam_dat->bc_data is NULL");

    if (mdl->mp == NULL)
        return err_msg(-1, 0, "mdl_llk: mp is NULL");

    // barcodes to test
    str_map *bcs = mdl->test_bcs;

    int fret = 0;

    uint16_t M = mdl->mp->M; // number of samples
    int max_ix = M + (M * (M-1) / 2); // number of singlets + doublets

    // allocate llk memory
    mdl->mf->bc_llks = calloc(bcs->n * max_ix, sizeof(double)); // n_bcs * ix. col major

    khash_t(kh_bc_dat) *bc_hash = bam_dat->bc_data;

    char dt[20];
    get_time(dt, 20);
    int bc_i;
    for (bc_i = 0; bc_i < bcs->n; ++bc_i){
        char *bc_name = str_map_str(bcs, bc_i);
        int s_ix;

        // progress
        int prog = (int)floor(100. * ((float)(bc_i+1) / (float)bcs->n));
        get_time(dt, 20);
        if (verbose) fprintf(stdout, 
                "\033[A\33[2K\r%s: getting likelihood %i%%\n", dt, prog);

        khint_t k_bc = kh_get(kh_bc_dat, bc_hash, bc_name);
        if (k_bc == kh_end(bc_hash)){
            log_msg("mdl_llk: barcode %s in mdl but not in data", bc_name);
            continue;
        }
        bc_data_t *bc_data = kh_val(bc_hash, k_bc);

        for (s_ix = 0; s_ix < max_ix; ++s_ix){
            mdl->mf->bc_llks[CMI( (bc_i), (s_ix), (bcs->n) )] = 0;

            // Loop through RNA
            if (bam_dat->has_rna){
                khash_t(khrmn) *mols = bc_data->rna_mols;
                khint_t k_rna;
                for (k_rna = kh_begin(mols); k_rna != kh_end(mols); ++k_rna){
                    if (!kh_exist(mols, k_rna)) continue;

                    rna_mol_t *rna_mol = kh_val(mols, k_rna);
                    if (rna_mol == NULL) continue;

                    double tmp_p = 1;
                    double prb_uc[2] = {1,1}; // store prob features
                    int c_stat;
                    for (c_stat = 0; c_stat < 2; ++c_stat){
                        // prior prob of contamination
                        fret = p_c(c_stat, mdl->mp->tau, &tmp_p);
                        if (fret < 0) return -1;
                        prb_uc[c_stat] *= tmp_p;
                        // gene expression probs
                        seq_gene_t *gene;
                        for (gene = rna_mol->genes.head; gene != NULL; gene = gene->next){
                            int32_t g_ix = gene->gene_id;
                            uint8_t sp = gene->splice;
                            fret = p_f_c(g_ix, sp, mdl->mp, c_stat, &tmp_p);
                            if (fret < 0) return -1;
                            prb_uc[c_stat] *= tmp_p;
                            if (tmp_p == 0){
                                log_msg("warning: barcode %s has p_f_c of %f", bc_name, tmp_p);
                            }
                        }
                        // get var probs
                        vac_t *v = rna_mol->vacs.head;
                        for (; v != NULL; v = v->next){
                            uint8_t allele = v->allele;
                            if (allele != 0 && allele != 1)
                                continue;
                            // prior prob of seq error
                            double prb_e0 = 1;
                            double prb_e1 = 1;
                            fret = p_e(0, mdl->mp->eps, &prb_e0);
                            if (fret < 0) return -1;
                            fret = p_e(1, mdl->mp->eps, &prb_e1);
                            if (fret < 0) return -1;

                            fret = p_b_gce(allele, v->vix, mdl->mp, s_ix, 0, c_stat, &tmp_p);
                            if (fret < 0) return -1;
                            prb_e0 *= tmp_p;

                            fret = p_b_gce(allele, v->vix, mdl->mp, s_ix, 1, c_stat, &tmp_p);
                            if (fret < 0) return -1;
                            prb_e1 *= tmp_p;

                            double prb_v = prb_e0 + prb_e1;
                            if (prb_v == 0){
                                log_msg("warning: barcode %s has prb_v of %f", bc_name, prb_v);
                            }
                            prb_uc[c_stat] *= prb_v;
                        }
                    }
                    double prb_u = prb_uc[0] + prb_uc[1];
                    if (prb_u == 0){
                        log_msg("warning: barcode %s has prb_u of %f", bc_name, prb_u);
                    }
                    mdl->mf->bc_llks[CMI( (bc_i), (s_ix), (bcs->n) )] += log(prb_u);
                }
            }
            // Loop through ATAC
            if (bam_dat->has_atac){
                khash_t(khaf) *frags = bc_data->atac_frags;
                khint_t k_atac;
                for (k_atac = kh_begin(frags); k_atac != kh_end(frags); ++k_atac){
                    if (!kh_exist(frags, k_atac)) continue;

                    atac_frag_t *atac_frag = kh_val(frags, k_atac);
                    if (atac_frag == NULL) continue;

                    double tmp_p = 1;
                    double prb_uc[2] = {1,1}; // store prob features
                    int c_stat;
                    for (c_stat = 0; c_stat < 2; ++c_stat){
                        // prior prob of contamination
                        fret = p_c(c_stat, mdl->mp->tau, &tmp_p);
                        if (fret < 0) return -1;
                        prb_uc[c_stat] *= tmp_p;
                        // get var probs
                        vac_t *v = atac_frag->vacs.head;
                        for (; v != NULL; v = v->next){
                            uint8_t allele = v->allele;
                            if (allele != 0 && allele != 1)
                                continue;
                            // prior prob of seq error
                            double prb_e0 = 1;
                            double prb_e1 = 1;
                            fret = p_e(0, mdl->mp->eps, &prb_e0);
                            if (fret < 0) return -1;
                            fret = p_e(1, mdl->mp->eps, &prb_e1);
                            if (fret < 0) return -1;

                            fret = p_b_gce(allele, v->vix, mdl->mp, s_ix, 0, c_stat, &tmp_p);
                            if (fret < 0) return -1;
                            prb_e0 *= tmp_p;

                            fret = p_b_gce(allele, v->vix, mdl->mp, s_ix, 1, c_stat, &tmp_p);
                            if (fret < 0) return -1;
                            prb_e1 *= tmp_p;

                            double prb_v = prb_e0 + prb_e1;
                            if (prb_v == 0){
                                log_msg("warning: barcode %s has prb_v of %f", bc_name, prb_v);
                            }
                            prb_uc[c_stat] *= prb_v;
                        }
                    }
                    double prb_u = prb_uc[0] + prb_uc[1];
                    if (prb_u == 0){
                        log_msg("warning: barcode %s has prb_u of %f", bc_name, prb_u);
                    }
                    mdl->mf->bc_llks[CMI( (bc_i), (s_ix), (bcs->n) )] += log(prb_u);
                }
            }
        }
    }
    return(0);
}

int mdl_get_best_llk(mdl_t *mdl){
    if (mdl == NULL)
        return 0;

    mdl_fit_t *mf = mdl->mf;

    mf->best_sng_llk = calloc(mdl->test_bcs->n, sizeof(double));
    mf->sec_sng_llk = calloc(mdl->test_bcs->n, sizeof(double));
    mf->best_dbl_llk = calloc(mdl->test_bcs->n, sizeof(double));
    mf->best_sng_ix = calloc(mdl->test_bcs->n, sizeof(uint32_t));
    mf->sec_sng_ix = calloc(mdl->test_bcs->n, sizeof(uint32_t));
    mf->best_dbl_ix = calloc(mdl->test_bcs->n, sizeof(uint32_t));
    mf->llr = calloc(mdl->test_bcs->n, sizeof(double));

    int M = mdl->samples->n;
    int max_ix = M + (M * (M-1) / 2);

    int bcs_n = mdl->test_bcs->n;
    
    int bc_i;
    for (bc_i = 0; bc_i < bcs_n; ++bc_i){
        uint32_t s_ix = 0;
        // singlet probs
        mf->best_sng_llk[bc_i] = mf->bc_llks[CMI( (bc_i), (s_ix), (bcs_n) )];
        mf->sec_sng_llk[bc_i] = mf->bc_llks[CMI( (bc_i), (s_ix), (bcs_n) )];
        mf->best_sng_ix[bc_i] = s_ix;
        mf->sec_sng_ix[bc_i] = s_ix;
        for (++s_ix; s_ix < M; ++s_ix){
            double tl = mf->bc_llks[CMI( (bc_i), (s_ix), (bcs_n) )];
            if (tl > mf->best_sng_llk[bc_i]){
                mf->best_sng_llk[bc_i] = tl;
                mf->best_sng_ix[bc_i] = s_ix;
            }
            else if (tl > mf->sec_sng_llk[bc_i]){
                mf->sec_sng_llk[bc_i] = tl;
                mf->sec_sng_ix[bc_i] = s_ix;
            }
        }
        // doublet probs
        mf->best_dbl_llk[bc_i] = mf->bc_llks[CMI( (bc_i), (s_ix), (bcs_n) )];
        mf->best_dbl_ix[bc_i] = s_ix;
        for (++s_ix; s_ix < max_ix; ++s_ix){
            double tl = mf->bc_llks[CMI( (bc_i), (s_ix), (bcs_n) )];
            if (tl > mf->best_dbl_llk[bc_i]){
                mf->best_dbl_llk[bc_i] = tl;
                mf->best_dbl_ix[bc_i] = s_ix;
            }
        }
        mf->llr[bc_i] = mf->best_sng_llk[bc_i] - mf->best_dbl_llk[bc_i];
    }
    return 0;
}


int mdl_est_pars(mdl_t *mdl, bam_data_t *bam_dat, obj_pars *objs){
    if (mdl == NULL || bam_dat == NULL || objs == NULL)
        return err_msg(-1, 0, "mdl_est_pars: arguments are NULL");

    int n_bcs = mdl->test_bcs->n;
    int n_genes = 0;
    if (objs->anno)
        n_genes = objs->anno->gene_ix->n;
    int n_vars = objs->gv->n_v;
    int n_samples = bcf_hdr_nsamples(objs->vcf_hdr);

    if (set_mdl_pars(mdl->mp, n_bcs, n_genes, n_vars, n_samples) < 0)
        return(-1);

    if (mdl_pars_est_alpha(mdl->mp, bam_dat, mdl->flt_bcs, 1.0) < 0)
        return(-1);

    // get valid indices
    int32_t nvar = objs->gv->n_e;
    int32_t *vixs = calloc(nvar, sizeof(int32_t));
    int i;
    for (i = 0; i < nvar; ++i) vixs[i] = i;
    mdl->mp->V = nvar;

    float **gm = ap_array_gt(objs->gv, objs->vcf_hdr, vixs, nvar, "GT");
    free(vixs);
    if (gm == NULL)
        return(-1);

    if (mdl_pars_add_gamma(mdl->mp, gm, nvar, n_samples) < 0)
        return(-1);

    // free gm
    for (i = 0; i < n_vars; ++i){
        free(gm[i]);
    }
    free(gm);
    
    if (mdl_pars_set_beta(mdl->mp) < 0)
        return(-1);

    if (mdl_pars_set_tau(mdl->mp, 0.5) < 0)
        return(-1);

    if (mdl_pars_set_eps(mdl->mp, 0.01) < 0)
        return(-1);

    return(0);
}

int fit_mdl_pars(bam_data_t *bam_dat, obj_pars *objs){
    if (bam_dat == NULL || objs == NULL)
        return err_msg(-1, 0, "fit_mdl_pars: arguments are NULL");

    if (bam_dat->bc_data == NULL)
        return err_msg(-1, 0, "fit_mdl_pars: rna and atac are NULL");

    if (objs->gv == NULL || objs->vcf_hdr == NULL)
        return err_msg(-1, 0, "fit_mdl_pars: variant data  is NULL");

    if (objs->verbose) log_msg("fitting model");

    mdl_t *mdl = mdl_alloc();
    if (mdl == NULL)
        return -1;

    if (objs->verbose) log_msg("setting barcodes");
    // sets all bc, flt_bcs, and test_bcs
    if (mdl_set_bcs(mdl, bam_dat, objs) < 0)
        return -1;

    if (objs->verbose) log_msg("setting samples");
    if (mdl_set_samples(mdl, objs) < 0)
        return -1;

    if (objs->verbose) log_msg("estimating parameters");
    if (mdl_est_pars(mdl, bam_dat, objs) < 0)
        return -1;

    if (objs->verbose) log_msg("getting likelihood");
    if (mdl_llk(mdl, bam_dat, objs->verbose) < 0)
        return -1;

    if (objs->verbose) log_msg("classifying droplets");
    if (mdl_get_best_llk(mdl) < 0)
        return -1;

    /*
    write_llk(mf, "llk_test.txt");
    */

    if (bam_dat->rna && write_alpha(mdl, objs, objs->out_fn) < 0)
        return -1;

    if (write_llk(mdl, objs->out_fn) < 0)
        return -1;

    if (write_samples(mdl, objs->out_fn) < 0)
        return -1;

    if (objs->verbose) log_msg("writing output");
    if (write_res(mdl, bam_dat, objs->out_fn) < 0)
        return -1;

    mdl_dstry(mdl);

    return(0);
}


int write_alpha(mdl_t *mdl, obj_pars *objs, char *fn){
    if (mdl == NULL)
        return err_msg(-1, 0, "write_alpha: arguments are NULL");

    if (mdl->mp == NULL || mdl->mp->alpha == NULL)
        return err_msg(-1, 0, "write_alpha: model hasn't been initialized");

    if (objs->anno == NULL)
        return err_msg(-1, 0, "write_alpha: annotation not provided");

    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;

    BGZF *fp;
    char *alpha_fn = ".alpha.txt.gz";
    char *out_alpha_fn = strcat2((const char*)fn, (const char*)alpha_fn);
    fp = bgzf_open(out_alpha_fn, "wg1");
    if (fp == 0){
        err_msg(-1, 0, "write_alpha: Could not open file %s\n", out_alpha_fn);
        free(out_alpha_fn);
        return -1;
    }

    size_t buf_size = decp + 1000;
    char *pstr = (char *)malloc(buf_size * sizeof(char));
    if (pstr == NULL)
        return err_msg(-1, 0, "write_alpha: %s", strerror(errno));
    int pstr_len;

    char *hdr = "gene\tambient\tcell\n";
    ret = bgzf_write(fp, hdr, strlen(hdr));

    mdl_pars_t *mp = mdl->mp;

    int g, gs, a, G = mp->G, GA = 3 * mp->G;
    char r_apnd[3][20] = {"_SPLICE", "_UNSPLICE", "_AMBIG"};

    for (gs = 0; gs < 3; ++gs){
        for (g = 0; g < G; ++g){
            char *g_str = str_map_str(objs->anno->gene_ix, g);
            ret = bgzf_write(fp, g_str, strlen(g_str));
            ret = bgzf_write(fp, r_apnd[gs], strlen(r_apnd[gs]));
            for (a = 0; a < 2; ++a){
                uint32_t ix = g + (gs * G) + (a * GA);
                double par = mp->alpha[ix];
                if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
                    return err_msg(-1, 0, "write_alpha: failed to convert %f to string", par);

                // write to file
                ret = bgzf_write(fp, &delim, 1);
                ret = bgzf_write(fp, pstr, pstr_len);
            }
            ret = bgzf_write(fp, &nl, 1);
        }
    }

    bgzf_close(fp);
    free(pstr);
    if (ret < 0){
        err_msg(-1, 0, "write_alpha: failed to write to file %s", fn);
        return -1;
    }
    free(out_alpha_fn);

    return 0;
}

int write_llk(mdl_t *mdl, char *fn){
    char nl = '\n';
    char delim = '\t';
    unsigned int decp = 8;
    int ret = 0;
    char *llk_fn = ".llk.txt.gz";
    str_map *bcs = mdl->test_bcs;
    int bcs_n = bcs->n;
    int M = mdl->samples->n;
    uint16_t max_ix = M + (M * (M-1) / 2);
    mdl_fit_t *mf = mdl->mf;

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "write_llk: failed to create output directory for %s", fn);

    BGZF *fp;
    char *out_llk_fn = strcat2((const char*)fn, (const char*)llk_fn);
    fp = bgzf_open(out_llk_fn, "wg1");
    if (fp == 0){
        err_msg(-1, 0, "write_llk: Could not open file %s\n", fn);
        free(out_llk_fn);
        return -1;
    }

    size_t buf_size = decp + 1000;
    char *pstr = (char *)malloc(buf_size * sizeof(char));
    if (pstr == NULL)
        return err_msg(-1, 0, "write_llk: %s", strerror(errno));

    int bc_i, s_ix, pstr_len;
    for (bc_i = 0; bc_i < bcs_n; ++bc_i){
        char *bc_name = str_map_str(bcs, bc_i);
        ret = bgzf_write(fp, bc_name, strlen(bc_name));
        for (s_ix = 0; s_ix < max_ix; ++s_ix){
            // convert to string
            double par = mf->bc_llks[CMI( (bc_i), (s_ix), (bcs_n) )];
            if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
                return err_msg(-1, 0, "write_llk: failed to convert %f to string", par);
            
            // write to file
            ret = bgzf_write(fp, &delim, 1);
            ret = bgzf_write(fp, pstr, pstr_len);
        }
        ret = bgzf_write(fp, &nl, 1);
    }
    bgzf_close(fp);
    free(pstr);
    if (ret < 0){
        err_msg(-1, 0, "write_llk: failed to write to file %s", fn);
        return -1;
    }
    free(out_llk_fn);

    return 0;
}

int write_samples(mdl_t *mdl, char *fn){
    char *s_fn = ".samples.txt";
    str_map *samples = mdl->samples;
    int M = samples->n;
    uint16_t max_ix = M + (M * (M-1) / 2);
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

    int s_ix, fret;
    for (s_ix = 0; s_ix < max_ix; ++s_ix){

        // if singlet
        if (s_ix < M){ // if singlet
            fret = fputs(samples->strs[s_ix], fp);
            fret = fputc(delim, fp);
            fret = fputs(samples->strs[s_ix], fp);
            fret = fputc(nl, fp);
        } else { // if doublet
            int s1 = 0, s2 = 0;
            ix2indices(s_ix, M, &s1, &s2);
            fret = fputs(samples->strs[s1], fp);
            fret = fputc(delim, fp);
            fret = fputs(samples->strs[s2], fp);
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
    int M = samples->n;
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
        "\tn_variants\tFRIP\tSinglet_best\tSinglet_best_llk"
        "\tSinglet_scnd\tSinglet_scnd_llk"
        "\tDoublet_sample1\tDoublet_sample2\tDoublet_llk"
        "\tLog_lk_ratio\n";

    fputs(hdr, fp);

    size_t buf_size = decp + 1000;
    char *pstr = (char *)malloc(buf_size * sizeof(char));
    if (pstr == NULL)
        return err_msg(-1, 0, "write_res: %s", strerror(errno));

    int bc_i, pstr_len;
    for (bc_i = 0; bc_i < bcs_n; ++bc_i){
        char *bc_name = str_map_str(bcs, bc_i);
        fret = fputs(bc_name, fp);

        khint_t k_bc = kh_get(kh_bc_dat, bam_dat->bc_data, bc_name);
        if (k_bc == kh_end(bam_dat->bc_data))
            return err_msg(-1, 0, "write_res: barcode %s not found", bc_name);

        bc_data_t *bc_data = kh_val(bam_dat->bc_data, k_bc);
        bc_counts *bcc = bc_data->bc_stats;

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
        int2strp(bcc->n_var, &pstr, &buf_size);
        fputs(pstr, fp);

        fputc(delim, fp);
        double2str_in((double)bcc->frip, &pstr, &buf_size, 4);
        fputs(pstr, fp);

        // write best singlet
        fret = fputc(delim, fp);
        fret = fputs(str_map_str(samples, mf->best_sng_ix[bc_i]), fp);

        // write best singlet llk
        double par = mf->best_sng_llk[bc_i];
        if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
            return err_msg(-1, 0, "write_llk: failed to convert %f to string", par);

        fret = fputc(delim, fp);
        fret = fputs(pstr, fp);

        // write 2nd singlet
        fret = fputc(delim, fp);
        fret = fputs(str_map_str(samples, mf->sec_sng_ix[bc_i]), fp);

        // write 2nd llk
        par = mf->sec_sng_llk[bc_i];
        if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
            return err_msg(-1, 0, "write_llk: failed to convert %f to string", par);

        fret = fputc(delim, fp);
        fret = fputs(pstr, fp);

        // write best doublet
        int s1 = 0, s2 = 0;
        ix2indices(mf->best_dbl_ix[bc_i], M, &s1, &s2);
        fret = fputc(delim, fp);
        fret = fputs(str_map_str(samples, s1), fp);
        fret = fputc(delim, fp);
        fret = fputs(str_map_str(samples, s2), fp);

        // write best doublet llk
        par = mf->best_dbl_llk[bc_i];
        if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
            return err_msg(-1, 0, "write_llk: failed to convert %f to string", par);

        fret = fputc(delim, fp);
        fret = fputs(pstr, fp);

        // write log likelihood ratio of sng to doublet
        par = mf->llr[bc_i];
        if ( (pstr_len = double2str_in(par, &pstr, &buf_size, decp)) < 0)
            return err_msg(-1, 0, "write_llk: failed to convert %f to string", par);

        fret = fputc(delim, fp);
        fret = fputs(pstr, fp);

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

