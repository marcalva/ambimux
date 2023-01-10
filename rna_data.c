
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
    init_seq_blist(&r->bases);
    seq_glist_init(&r->genes);
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
    free_seq_blist(&r->bases);
    seq_glist_free(&r->genes);
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

    int cret;

    seq_blist_t *b_cpy = copy_seq_blist(&r->bases, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    cpy->bases = *b_cpy;
    free(b_cpy);

    seq_glist_t *g_cpy = seq_glist_cpy(&r->genes, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    cpy->genes = *g_cpy;
    free(g_cpy);

    return(cpy);
}

int rna_read1_equal(const rna_read1_t *r1, const rna_read1_t *r2, int cmp_qual){
    if (r1 == NULL || r2 == NULL)
        return err_msg(-1, 0, "rna_read1_equal: r1 or r2 are null");

    if (regioncmp(r1->loc, r2->loc) != 0)
        return(0);

    if (r1->bases.n != r2->bases.n)
        return(0);

    if (r1->genes.n != r2->genes.n)
        return(0);

    seq_base_t *b1 = r1->bases.head, *b2 = r2->bases.head;
    while (b1 != NULL && b2 != NULL){
        g_pos pos1 = b1->pos, pos2 = b2->pos;
        if (poscmp(pos1, pos2) != 0)
            return(0);

        uint8_t base1 = b1->base, base2 = b2->base;
        if (base1 != base2)
            return(0);

        uint8_t qual1 = b1->qual, qual2 = b2->qual;
        if (cmp_qual != 0 && qual1 != qual2)
            return(0);

        b1 = b1->next;
        b2 = b2->next;
    }
    if (b1 != NULL || b2 != NULL)
        return err_msg(-1, 0, "rna_read1_equal: equal n but different num. of "
                " bases, there is a bug");

    seq_gene_t *g1 = r1->genes.head, *g2 = r2->genes.head;
    while (g1 != NULL && g2 != NULL){
        int32_t id1 = g1->gene_id, id2 = g2->gene_id;
        if (id1 != id2)
            return(0);

        uint8_t s1 = g1->splice, s2 = g2->splice;
        if (s1 != s2)
            return(0);

        g1 = g1->next;
        g2 = g2->next;
    }
    if (g1 != NULL || g2 != NULL)
        return err_msg(-1, 0, "rna_read1_equal: equal n but different num. of "
                " genes, there is a bug");

    return(1);
}

int rna_read1_add_gene(rna_read1_t *r, const seq_gene_t *gene){
    if (r == NULL || gene == NULL)
        return err_msg(-1, 0, "rna_read1_add_gene: argument is null");

    if (seq_glist_add_gene(&r->genes, gene) < 0) return(-1);

    return(0);
}

int rna_read1_add_base(rna_read1_t *r, const seq_base_t *base){
    if (r == NULL || base == NULL)
        return err_msg(-1, 0, "rna_read1_add_base: argument is nul");

    if (blist_add_base(&r->bases, base, 1, 1) < 0) return(-1);

    return(0);
}

int rna_read1_match_qual(rna_read1_t *r, const rna_read1_t *cmp){
    if (r == NULL || cmp == NULL)
        return err_msg(-1, 0, "rna_read1_match_qual: argument is null");

    if (seq_blist_match_qual(&r->bases, &cmp->bases) < 0)
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
    rd->rds_per_dup = (size_t *)calloc(rd->m, sizeof(size_t));
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
    int i;
    for (i = 0; i < rd->size; ++i)
        rna_read1_free(&rd->dups[i]);

    free(rd->dups);
    free(rd->rds_per_dup);
}

void rna_dups_dstry(rna_dups_t *rd){
    if (rd == NULL) return;

    rna_dups_free(rd);

    free(rd);
}

int rna_dups_add_read(rna_dups_t *rd, const rna_read1_t *r){
    if (rd == NULL || r == NULL)
        return err_msg(-1, 0, "rna_dups_add_read: arguments are null");

    int i, ret;
    for (i = 0; i < rd->size; ++i){
        int rcmp = rna_read1_equal(&rd->dups[i], r, 0);
        if (rcmp < 0)
            return(-1);
        if (rcmp == 1){
            if (rna_read1_match_qual(&rd->dups[i], r) < 0)
                return(-1);
            ++rd->rds_per_dup[i];
            return(0);
        }
    }
    while (rd->size >= rd->m){
        rd->m = rd->size + 1;
        rd->dups = realloc(rd->dups, rd->m * sizeof(rna_read1_t));
        rd->rds_per_dup = realloc(rd->rds_per_dup, rd->m * sizeof(size_t));
        if (rd->dups == NULL || rd->rds_per_dup == NULL)
            return err_msg(-1, 0, "rna_dups_add_read: %s", strerror(errno));
    }

    // copy read pair
    rna_read1_t *tmp = rna_read1_cpy(r, &ret);
    if (ret < 0)
        return(-1);

    rd->dups[rd->size] = *tmp;
    rd->rds_per_dup[rd->size] = 1;
    free(tmp);

    ++rd->size;
    return(0);
}

/*******************************************************************************
 * rna_mol_t
 ******************************************************************************/

int rna_mol_init(rna_mol_t *rmol){
    if (rmol == NULL) return(-1);
    rmol->dups = rna_dups_alloc();
    if (rmol->dups == NULL)
        return(-1);
    init_g_region(&rmol->loc);
    init_seq_blist(&rmol->bases);
    seq_glist_init(&rmol->genes);
    init_vacs(&rmol->vacs);
    rmol->n_reads = 0;
    rmol->_dedup = 0;
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

    rna_dups_dstry(rmol->dups);
    rmol->dups = NULL;
    free_seq_blist(&rmol->bases);
    seq_glist_free(&rmol->genes);
    free_vacs(&rmol->vacs);
    rmol->n_reads = 0;
    rmol->_dedup = 0;
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

    // copy dups
    cpy->dups = rna_dups_alloc();
    if (cpy->dups == NULL){
        *ret = err_msg(-1, 0, "rna_mol_cpy: %s", strerror(errno));
        return(NULL);
    }
    int i;
    for (i = 0; i < m->dups->size; ++i){
        if (rna_dups_add_read(cpy->dups, m->dups->dups + i) < 0)
            return(NULL);
    }

    cpy->loc = m->loc;

    int cret;

    seq_blist_t *b_cpy = copy_seq_blist(&m->bases, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    cpy->bases = *b_cpy;
    free(b_cpy);

    seq_glist_t *g_cpy = seq_glist_cpy(&m->genes, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    cpy->genes = *g_cpy;
    free(g_cpy);

    cpy->n_reads = m->n_reads;

    cpy->_dedup = m->_dedup;

    return(cpy);
}

int rna_mol_add_read(rna_mol_t *mol, const rna_read1_t *r){
    if (mol == NULL || r == NULL)
        return err_msg(-1, 0, "rna_mol_add_read: 'mol' or 'r' is null");
    return rna_dups_add_read(mol->dups, r);
}

int rna_mol_dedup(rna_mol_t *mol){
    if (mol == NULL)
        return err_msg(-1, 0, "rna_mol_dedup: mol is null");

    if (mol->dups == NULL)
        return err_msg(-1, 0, "rna_mol_dedup: mol->dups is null");

    rna_dups_t *dups = mol->dups;

    // if no dups, return 0 successfully
    if (dups->size == 0)
        return(0);

    int ix_best = 0;
    size_t max_c = 0; // store read count
    size_t rd_w_max = 0; // number of read pairs with max_c read counts

    int i;
    for (i = 0; i < dups->size; ++i){
        // check for bugs
        if (dups->rds_per_dup[i] == 0)
            return err_msg(-1, 0, "rna_mol_dedup: a dup read has 0 "
                    "supporting reads, there is a bug");

        if (dups->rds_per_dup[i] == max_c){
            ++rd_w_max;
        } else if (dups->rds_per_dup[i] > max_c){
            max_c = dups->rds_per_dup[i];
            ix_best = i;
            rd_w_max = 1;
        }
    }

    // check for bugs
    if (max_c == 0 || rd_w_max == 0)
        return err_msg(-1, 0, "atac_rna_dedup: no best duplicate found "
                ", there is a bug");

    // initialize molecule
    rna_read1_t rb = dups->dups[ix_best];

    int cret;

    seq_blist_t *b_cpy = copy_seq_blist(&rb.bases, &cret);
    if (cret < 0) return(-1);
    mol->bases = *b_cpy;
    free(b_cpy);

    seq_glist_t *g_cpy = seq_glist_cpy(&rb.genes, &cret);
    if (cret < 0) return(-1);
    mol->genes = *g_cpy;
    free(g_cpy);

    mol->n_reads = dups->rds_per_dup[ix_best];

    // free dups
    rna_dups_dstry(dups);
    mol->dups = NULL;

    mol->_dedup = 1;

    return(0);
}

int rna_mol_var_call(rna_mol_t *m, g_var_t *gv, contig_map *cmap, 
        uint8_t min_qual){
    if (m == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "rna_mol_var_call: argument is null");
    }

    int ret = seq_blist_call_var(&m->bases, &m->vacs, gv, cmap, 
            min_qual);

    // free the bases
    free_seq_blist(&m->bases);

    return(ret);
}

/*******************************************************************************
 * bc_rna_t
 ******************************************************************************/

/*******************************************************************************
 * bam_rna_t
 ******************************************************************************/

/*******************************************************************************
 * miscellaneous
 ******************************************************************************/

/*
void count_umis(bam_rna_t *bam_r, int *n_umi, int *n_reads, int *n_bc){
    *n_umi = 0;
    *n_reads = 0;
    *n_bc = 0;
    khint_t k;
    bc_rna_t *b;

    for (k = kh_begin(bam_r->bc_rna); k != kh_end(bam_r->bc_rna); ++k){
        if (!kh_exist(bam_r->bc_rna, k)) continue;
        (*n_bc)++;

        b = kh_val(bam_r->bc_rna, k);

        *n_umi += kh_size(b->mols);

        khint_t km;
        for (km = kh_begin(b->mols); km != kh_end(b->mols); ++km){
            if (!kh_exist(b->mols, km)) continue;
            rna_mol_t *m = kh_val(b->mols, km);
            *n_reads += m->n_reads;
        }


    }
}

void print_bam_rna(bam_rna_t *b, gene_anno_t *anno, bcf_hdr_t *vcf_hdr, 
        g_var_t *gv){
    if (b == NULL){
        err_msg(-1, 0, "print_bam_rna: arguments must not be NULL");
        return;
    }

    fprintf(stdout, "printing bam rna\n");

    khint_t k;
    for (k = kh_begin(b->bc_rna); k != kh_end(b->bc_rna); ++k){
        if (!kh_exist(b->bc_rna, k)) continue;

        // get frags from bam_atac
        char *barcode = kh_key(b->bc_rna, k);
        bc_rna_t *v = kh_val(b->bc_rna, k);
        if (v == NULL) continue;

        size_t n_m = kh_size(v->mols);
        // loop through molecules
        khint_t km;
        for (km = kh_begin(v->mols); km != kh_end(v->mols); ++km){
            if (!kh_exist(v->mols, km)) continue;

            char *umi = kh_key(v->mols, km);
            rna_mol_t *f_i = kh_val(v->mols, km);
            if (f_i == NULL) continue;
            fprintf(stdout, "%s: (%zu mols) ", barcode, n_m);
            fprintf(stdout, "region: %i:%"PRIhts_pos"-%"PRIhts_pos"\n", 
                    f_i->loc.rid, f_i->loc.start, f_i->loc.end);
            fprintf(stdout, "\tUMI:%s %zu reads\n", umi, f_i->n_reads);

            // print genes
            seq_gene_t *sg = f_i->genes.head;
            while (sg != NULL){
                char *gene_name = str_map_str(anno->gene_ix, sg->gene_id);
                fprintf(stdout, "\tgene:");
                fprintf(stdout, " gene=%s", gene_name);
                fprintf(stdout, " splice=%u\n", sg->splice);
                sg = sg->next;
            }

            // print bases
            seq_base_t *sb = f_i->bases.head;
            while (sb != NULL){
                const char *v_chr = bcf_hdr_id2name(vcf_hdr, sb->pos.rid);
                fprintf(stdout, "\tbase:");
                fprintf(stdout, " %s (%i)", v_chr, sb->pos.rid);
                fprintf(stdout, " %"PRIhts_pos"", sb->pos.pos);
                fprintf(stdout, " %u", sb->base);
                fprintf(stdout, " %u\n", sb->qual);
                sb = sb->next;
            }

            // print variant calls
            vac_t *va = f_i->vacs.head;
            size_t va_n = f_i->vacs.n;
            fprintf(stdout, "\t%zu variants\n", va_n);
            while (va != NULL){
                var_t *var = gv_vari(gv, va->vix);
                fprintf(stdout, "\tvar:");
                fprintf(stdout, "\t%s", var->b->d.id);
                fprintf(stdout, "\t%"PRIhts_pos"", var->b->pos);
                fprintf(stdout, "\t%u\n", va->allele);
                va = va->next;
            }
        }
    }
}
*/

