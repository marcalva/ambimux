
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "atac_data.h"
#include "str_util.h"
#include "variants.h"
#include "region.h"
#include "counts.h"
#include "overlap.h"
#include "htslib/hts.h"

/*******************************************************************************
 * miscellaneous
 ******************************************************************************/

/* Count number of fragments in bam_atac
 */
/*
void count_frags(bam_atac_t *bam_a, int *n_frag, int *n_reads, int *n_bc);

uint64_t get_n_buckets(bam_atac_t *bam_a);

void print_atac_read(atac_read1_t *ar);

void print_vac_bam(bam_atac_t *b, g_var_t *gv);

void print_frag_dup(bam_atac_t *b);
*/

/*******************************************************************************
 * miscellaneous
 ******************************************************************************/

/*
void count_frags(bam_atac_t *bam_a, int *n_frag, int *n_reads, int *n_bc){
    *n_frag = 0;
    *n_reads = 0;
    *n_bc = 0;
    khint_t k;
    bc_atac_t *b;

    for (k = kh_begin(bam_a->bc_dat); k != kh_end(bam_a->bc_dat); ++k){
        if (!kh_exist(bam_a->bc_dat, k)) continue;
        (*n_bc)++;

        b = kh_val(bam_a->bc_dat, k);

        *n_frag += kh_size(b->frags);

        if (b->pairs != NULL)
            *n_reads += kh_size(b->pairs);
    }
    printf("n_frag = %i; n_reads = %i; n_bc = %i\n", *n_frag, *n_reads, *n_bc);
}

uint64_t get_n_buckets(bam_atac_t *bam_a){
    if (bam_a == NULL || bam_a->bc_dat == NULL)
        return(0);

    uint64_t nb_frag = 0, nb_dup = 0, nb_pair = 0;
    uint64_t nb_bc = kh_n_buckets(bam_a->bc_dat);
    khint_t kb;
    for (kb = kh_begin(bam_a->bc_dat); kb != kh_end(bam_a->bc_dat); ++kb){
        if (!kh_exist(bam_a->bc_dat, kb)) continue;
        bc_atac_t *bc_atac = kh_val(bam_a->bc_dat, kb);
        if (bc_atac == NULL) continue;
        if (bc_atac->frags) nb_frag += (uint64_t)kh_n_buckets(bc_atac->frags);
        if (bc_atac->dups) nb_dup += (uint64_t)kh_n_buckets(bc_atac->dups);
        if (bc_atac->pairs) nb_pair += (uint64_t)kh_n_buckets(bc_atac->pairs);
    }
    printf("nb_frag = %" PRIu64 "; nb_dup = %" PRIu64 "; "
           "nb_pair = %" PRIu64 "; nb_barcodes = %" PRIu64 "\n", 
            nb_frag, nb_dup, nb_pair, nb_bc);
    return(nb_bc);
}

void print_vac_bam(bam_atac_t *b, g_var_t *gv){
    if (b == NULL){
        err_msg(-1, 0, "print_vac_bam: arguments must not be NULL");
        return;
    }

    fprintf(stdout, "printing variant allele calls for atac fragments\n");

    khint_t k;
    for (k = kh_begin(b->bc_dat); k != kh_end(b->bc_dat); ++k){
        if (!kh_exist(b->bc_dat, k)) continue;

        // get frags from bam_atac
        char *barcode = kh_key(b->bc_dat, k);
        bc_atac_t *v = kh_val(b->bc_dat, k);
        if (v == NULL) continue;

        size_t n_f = kh_size(v->frags);
        // loop through frags
        khint_t kf;
        for (kf = kh_begin(v->frags); kf != kh_end(v->frags); ++kf){
            if (!kh_exist(v->frags, kf)) continue;
            g_reg_pair regp = kh_key(v->frags, kf);

            atac_frag_t *f_i = kh_val(v->frags, kf);
            if (f_i == NULL) continue;
            fprintf(stdout, "%s: (%zu frags) ", barcode, n_f);
            fprintf(stdout, "r1: %i:%"PRIhts_pos"-%"PRIhts_pos"", 
                    regp.r1.rid, regp.r1.start, regp.r1.end);
            fprintf(stdout, " r2: %i:%"PRIhts_pos"-%"PRIhts_pos"\n", 
                    regp.r2.rid, regp.r2.start, regp.r2.end);

            // print bases
            seq_base_t *sb = f_i->bases.head;
            fprintf(stdout, "bases:\n");
            while (sb != NULL){
                fprintf(stdout, "\t%i", sb->pos.rid);
                fprintf(stdout, "\t%"PRIhts_pos"", sb->pos.pos);
                fprintf(stdout, "\t%u", sb->base);
                fprintf(stdout, "\t%u\n", sb->qual);
                sb = sb->next;
            }

            fprintf(stdout, "variant alleles:\n");
            // print variant calls
            vac_t *va = f_i->vacs.head;
            size_t va_n = f_i->vacs.n;
            fprintf(stdout, "\t%zu variants\n", va_n);
            while (va != NULL){
                var_t * var = gv_vari(gv, va->vix); 
                if (var == NULL)
                    fprintf(stdout, "\t%i is invalid\n", va->vix);
                bcf1_t *rec = var->b;
                fprintf(stdout, "\t%s", rec->d.id);
                fprintf(stdout, "\t%"PRIhts_pos"", rec->pos);
                fprintf(stdout, "\t%u\n", va->allele);
                va = va->next;
            }
        }
    }
}

void print_frag_dup(bam_atac_t *b){
    uint32_t counts[100] = {0}; // 1-5, 6-10, 11-15, ...
    int i, m = 1;
    khint_t k;
    for (k = kh_begin(b->bc_dat); k != kh_end(b->bc_dat); ++k){
        if (!kh_exist(b->bc_dat, k)) continue;

        // get frags from bam_atac
        bc_atac_t *v = kh_val(b->bc_dat, k);
        if (v == NULL) continue;

        // loop through frags
        khint_t kf;
        for (kf = kh_begin(v->frags); kf != kh_end(v->frags); ++kf){
            if (!kh_exist(v->frags, kf)) continue;

            atac_frag_t *f_i = kh_val(v->frags, kf);
            if (f_i == NULL) continue;

            for (i = 1; i <= 99; ++i){
                if (f_i->s <= (m * i)){
                    break;
                }
                ++counts[i-1];
            }
        }
    }

    for (i = 0; i < 100; ++i){
        fprintf(stdout, "%i-%i: %" PRIu32 "\n", i*m, (i+1)*m, counts[i]);
    }
    
}

void print_dups_n(bam_atac_t *b){
    // get number of duplicates
    uint32_t counts[100] = {0};
    int i;
    khint_t k;
    for (k = kh_begin(b->bc_dat); k != kh_end(b->bc_dat); ++k){
        if (!kh_exist(b->bc_dat, k)) continue;

        // get dups from bam_atac
        bc_atac_t *v = kh_val(b->bc_dat, k);
        if (v == NULL) continue;

        // loop through dups
        khint_t kd;
        for (kd = kh_begin(v->dups); kd != kh_end(v->dups); ++kd){
            if (!kh_exist(v->dups, kd)) continue;

            atac_dups_t *d_i = kh_val(v->dups, kd);
            if (d_i == NULL) continue;

            size_t s = 0;
            if (d_i->size > 100)
                s = 99;
            else
                s = d_i->size - 1;

            ++counts[s];
        }
    }

    for (i = 0; i < 100; ++i){
        fprintf(stdout, "%i: %" PRIu32 "\n", i+1, counts[i]);
    }
    
}

*/

