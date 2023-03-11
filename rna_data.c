
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
    base_list_init(&r->bl);
    gene_list_init(&r->gl);
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
    base_list_free(&r->bl);
    gene_list_free(&r->gl);
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

    if (base_list_cpy(&cpy->bl, &r->bl) < 0){
        *ret = -1;
        free(cpy);
        return(NULL);
    }

    if (gene_list_cpy(&cpy->gl, &r->gl) < 0){
        *ret = -1;
        free(cpy);
        return(NULL);
    }

    return(cpy);
}

int rna_read1_cmp(const rna_read1_t *r1, const rna_read1_t *r2, int cmp_qual){
    int rc = regioncmp(r1->loc, r2->loc);
    if (rc != 0) return(rc);

    int bsc = ml_size(&r1->bl) - ml_size(&r2->bl);
    if (bsc != 0) return(bsc);

    int gsc = ml_size(&r1->gl) - ml_size(&r2->gl);
    if (gsc != 0) return(gsc);

    int bc = base_list_cmp(r1->bl, r2->bl, cmp_qual);
    if (bc != 0) return(bc);

    int gc = gene_list_cmp(r1->gl, r2->gl);
    if (gc != 0) return(gc);

    return(0);
}

int rna_read1_equal(const rna_read1_t *r1, const rna_read1_t *r2, int cmp_qual){
    if (r1 == NULL || r2 == NULL)
        return err_msg(-1, 0, "rna_read1_equal: r1 or r2 are null");

    if (regioncmp(r1->loc, r2->loc) != 0)
        return(0);

    if (ml_size(&r1->bl) != ml_size(&r2->bl))
        return(0);

    if (ml_size(&r1->gl) != ml_size(&r2->gl))
        return(0);

    if ( base_list_equal(r1->bl, r2->bl, cmp_qual) != 1 )
        return(0);

    if ( gene_list_equal(r1->gl, r2->gl) != 1 )
        return(0);

    return(1);
}

int rna_read1_add_gene(rna_read1_t *r, const seq_gene_t *gene){
    if (r == NULL || gene == NULL)
        return err_msg(-1, 0, "rna_read1_add_gene: argument is null");

    if ( gene_list_insert(&r->gl, *gene, 0, 0) < 0 ) return(-1);

    return(0);
}

int rna_read1_add_base(rna_read1_t *r, seq_base_t base){
    if (r == NULL)
        return err_msg(-1, 0, "rna_read1_add_base: argument is null");

    if (base_list_insert(&r->bl, base, 1, 1) < 0) return(-1);

    return(0);
}

int rna_read1_match_qual(rna_read1_t *r, const rna_read1_t *cmp){
    if (r == NULL || cmp == NULL)
        return err_msg(-1, 0, "rna_read1_match_qual: argument is null");

    if (base_list_match_qual(&r->bl, &cmp->bl) < 0)
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
    size_t i;
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

    if (rd->dups == NULL || rd->rds_per_dup == NULL)
        return err_msg(-1, 0, "rna_dups_add_read: rd->dups or rds_per_dup is null");

    size_t i;
    int ret;
    for (i = 0; i < rd->size; ++i){
        int rcmp = rna_read1_cmp(&rd->dups[i], r, 0);
        if (rcmp == 0){
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

    assert(tmp->loc.rid >= 0);

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
    base_list_init(&rmol->bl);
    vac_list_init(&rmol->vl);
    gene_list_init(&rmol->gl);
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
    base_list_free(&rmol->bl);
    vac_list_free(&rmol->vl);
    gene_list_free(&rmol->gl);
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
        rna_mol_dstry(cpy);
        *ret = err_msg(-1, 0, "rna_mol_cpy: %s", strerror(errno));
        return(NULL);
    }
    size_t ix;
    for (ix = 0; ix < m->dups->size; ++ix){
        if (rna_dups_add_read(cpy->dups, m->dups->dups + ix) < 0){
            *ret = -1;
            rna_mol_dstry(cpy);
            return(NULL);
        }
    }

    cpy->loc = m->loc;

    if (base_list_cpy(&cpy->bl, &m->bl) < 0){
        *ret = -1;
        rna_mol_dstry(cpy);
        return(NULL);
    }

    if (vac_list_cpy(&cpy->vl, &m->vl) < 0){
        *ret = -1;
        rna_mol_dstry(cpy);
        return(NULL);
    }

    if (gene_list_cpy(&cpy->gl, &m->gl) < 0){
        *ret = -1;
        rna_mol_dstry(cpy);
        return(NULL);
    }

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

    size_t i;
    for (i = 0; i < dups->size; ++i){
        // check for bugs
        assert(dups->rds_per_dup[i] > 0);

        if (dups->rds_per_dup[i] == max_c){
            ++rd_w_max;
        } else if (dups->rds_per_dup[i] > max_c){
            max_c = dups->rds_per_dup[i];
            ix_best = i;
            rd_w_max = 1;
        }
    }

    // check for bugs
    assert(max_c > 0 && rd_w_max > 0);

    // initialize molecule
    rna_read1_t rb = dups->dups[ix_best];

    mol->loc = rb.loc;

    if ( base_list_cpy(&mol->bl, &rb.bl) < 0 )
        return err_msg(-1, 0, "rna_mol_dedup: failed to copy base list");

    if ( gene_list_cpy(&mol->gl, &rb.gl) < 0 )
        return err_msg(-1, 0, "rna_mol_dedup: failed to copy gene list");

    mol->n_reads = dups->rds_per_dup[ix_best];

    // free dups
    rna_dups_dstry(dups);
    mol->dups = NULL;

    mol->_dedup = 1;

    return(0);
}

int rna_mol_var_call(rna_mol_t *m, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (m == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "rna_mol_var_call: argument is null");
    }

    int ret = vac_list_call_var(&m->bl, &m->vl, gv, cmap, min_qual);

    /*
    if (1 && ml_size(&m->vl) > 0){
        printf("var list::\n");
        ml_node_t(vac_list) *nv;
        for (nv = ml_begin(&m->vl); nv; nv = ml_node_next(nv)){
            seq_vac_t vi = ml_node_val(nv);
            printf("\tvix=%i,allele=%u,qual=%u\n", vi.vix, vi.allele, vi.qual);
        }
    }
    */

    // free the bases
    base_list_free(&m->bl);

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
            seq_vac_t *va = f_i->vacs.head;
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

