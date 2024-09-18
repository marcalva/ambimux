
#include "rna_data.h"
#include <stdlib.h>
#include <string.h>
#include "kbtree.h"
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
    seq_base_l_init(&r->bl);
    seq_gene_l_init(&r->gl);
    r->umi = 0;
    r->n_reads = 0;
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
    init_g_region(&r->loc);
    seq_base_l_free(&r->bl);
    seq_gene_l_free(&r->gl);
    r->umi = 0;
    r->n_reads = 0;
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

    if (seq_base_l_cpy(&cpy->bl, &r->bl) < 0){
        *ret = -1;
        free(cpy);
        return(NULL);
    }

    if (seq_gene_l_cpy(&cpy->gl, &r->gl) < 0){
        *ret = -1;
        free(cpy);
        return(NULL);
    }

    cpy->umi = r->umi;
    cpy->n_reads = r->n_reads;

    return(cpy);
}

int rna_read1_cmp(const rna_read1_t *r1, const rna_read1_t *r2, int cmp_qual){
    int rc = regioncmp(r1->loc, r2->loc);
    if (rc != 0) return(rc);

    int bsc = ml_size(&r1->bl) - ml_size(&r2->bl);
    if (bsc != 0) return(bsc);

    int gsc = ml_size(&r1->gl) - ml_size(&r2->gl);
    if (gsc != 0) return(gsc);

    int bc = seq_base_l_cmp(r1->bl, r2->bl, cmp_qual);
    if (bc != 0) return(bc);

    int gc = seq_gene_l_cmp(r1->gl, r2->gl);
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

    if ( seq_base_l_equal(r1->bl, r2->bl, cmp_qual) != 1 )
        return(0);

    if ( seq_gene_l_equal(r1->gl, r2->gl) != 1 )
        return(0);

    return(1);
}

int rna_read1_add_gene(rna_read1_t *r, const seq_gene_t *gene){
    if (r == NULL || gene == NULL)
        return err_msg(-1, 0, "rna_read1_add_gene: argument is null");

    if ( seq_gene_l_insert(&r->gl, *gene, 0, 0) < 0 ) return(-1);

    return(0);
}

int rna_read1_add_base(rna_read1_t *r, seq_base_t base){
    if (r == NULL)
        return err_msg(-1, 0, "rna_read1_add_base: argument is null");

    if (seq_base_l_insert(&r->bl, base, 1, 1) < 0) return(-1);

    return(0);
}

int rna_read1_match_qual(rna_read1_t *r, const rna_read1_t *cmp){
    if (r == NULL || cmp == NULL)
        return err_msg(-1, 0, "rna_read1_match_qual: argument is null");

    if (seq_base_l_match_qual(&r->bl, &cmp->bl) < 0)
        return(-1);

    return(0);
}

/*******************************************************************************
 * rna_dups_t
 ******************************************************************************/

int rna_dups_init(rna_dups_t *rd){
    if (rd == NULL) return(0);

    rd->umi = 0;

    mv_init(&rd->rds);

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

    rd->umi = 0;

    size_t i;
    for (i = 0; i < mv_size(&rd->rds); ++i){
        rna_read1_t *r = &mv_i(&rd->rds, i);
        rna_read1_free(r);
    }
    mv_free(&rd->rds);
}

void rna_dups_dstry(rna_dups_t *rd){
    if (rd == NULL) return;

    rna_dups_free(rd);

    free(rd);
}

int rna_dups_add_read(rna_dups_t *rd, const rna_read1_t *r){
    if (rd == NULL || r == NULL)
        return err_msg(-1, 0, "rna_dups_add_read: arguments are null");

    // check if read is empty/unitialized
    if (r->loc.rid == -1)
        return err_msg(-1, 0, "rna_dups_add_read: read is empty (rid == -1)");

    int ret;
    // copy read pair
    rna_read1_t *tmp = rna_read1_cpy(r, &ret);
    if (ret < 0)
        return err_msg(-1, 0, "rna_dups_add_read: failed to copy rna_read1_t");

    // add read pair
    size_t i;
    for (i = 0; i < mv_size(&rd->rds); ++i){
        rna_read1_t *r1 = &mv_i(&rd->rds, i);
        // if read pair is found, set existing pair's base qualities to max
        //  of the two.
        if (rna_read1_equal(r1, tmp, 0) == 1){
            if (rna_read1_match_qual(r1, tmp) < 0) {
                rna_read1_dstry(tmp);
                return err_msg(-1, 0, "rna_dups_add_read: failed to match quality");
            }
            r1->n_reads += 1;
            rna_read1_dstry(tmp);
            tmp = NULL;
            return(0);
        }
    }
    if (mv_push(mv_rr, &rd->rds, *tmp) < 0) {
        rna_read1_dstry(tmp);
        return err_msg(-1, 0, "rna_dups_add_read: failed to add read to dups");
    }
    free(tmp);

    return 0;
}

rna_mol_t *rna_dups_dedup(rna_dups_t *dups, int *ret){
    *ret = 0;
    if (dups == NULL) {
        *ret = err_msg(-1, 0, "rna_dups_dedup: mol is null");
        return NULL;
    }

    uint32_t max_c = 0; // store read count
    uint32_t rd_w_max = 0; // number of read pairs with max_c read counts
    rna_read1_t *rd_best = NULL; // store the best read (highest counts)

    size_t i;
    for (i = 0; i < mv_size(&dups->rds); ++i) {
        rna_read1_t *rd  = &mv_i(&dups->rds, i);
        uint32_t rds_per_dup = rd->n_reads;

        if (rds_per_dup == 0) {
            *ret = err_msg(-1, 0, "rna_dups_dedup: rna_read1_t object in 'dups' has 0 underlying reads");
            return NULL;
        }

        if (rds_per_dup == max_c){
            ++rd_w_max;
        } else if (rds_per_dup > max_c){
            max_c = rds_per_dup;
            rd_best = rd;
            rd_w_max = 1;
        }
    }

    if (rd_w_max == 0) {
        *ret = err_msg(-1, 0, "rna_dups_dedup: no reads with read count");
        return NULL;
    }
    if (rd_best == NULL) {
        *ret = err_msg(-1, 0, "rna_dups_dedup: failed to find best read");
        return NULL;
    }

    // skip if ambiguous
    if (rd_w_max > 1)
        return NULL;

    rna_mol_t *mol = rna_mol_alloc();
    if (mol == NULL) {
        *ret = err_msg(-1, 0, "rna_dups_dedup: failed to allocate rna_mol_t");
        return NULL;
    }

    mol->umi = rd_best->umi;
    mol->loc = rd_best->loc;

    // copy the base list
    if ( seq_base_l_cpy(&mol->bl, &rd_best->bl) < 0 ) {
        *ret = err_msg(-1, 0, "rna_dups_dedup: failed to copy base list");
        rna_mol_dstry(mol);
        return NULL;
    }

    // copy the gene list
    if ( seq_gene_l_cpy(&mol->gl, &rd_best->gl) < 0 ) {
        *ret = err_msg(-1, 0, "rna_dups_dedup: failed to copy gene list");
        rna_mol_dstry(mol);
        return NULL;
    }

    // set the number of supporting reads
    mol->n_reads = max_c;

    return mol;
}

/*******************************************************************************
 * rna_mol_t
 ******************************************************************************/

int rna_mol_init(rna_mol_t *rmol){
    if (rmol == NULL) return(-1);
    rmol->umi = 0;
    init_g_region(&rmol->loc);
    if (seq_base_l_init(&rmol->bl) < 0)
        return err_msg(-1, 0, "rna_mol_init: failed to initialize");
    if (seq_vac_l_init(&rmol->vl) < 0)
        return err_msg(-1, 0, "rna_mol_init: failed to initialize");
    if (seq_gene_l_init(&rmol->gl) < 0)
        return err_msg(-1, 0, "rna_mol_init: failed to initialize");
    rmol->n_reads = 0;
    return(0);
}

rna_mol_t *rna_mol_alloc(){
    rna_mol_t *rmol = (rna_mol_t *)calloc(1, sizeof(rna_mol_t));
    if (rmol == NULL){
        err_msg(-1, 0, "rna_mol_alloc: %s", strerror(errno));
        return(NULL);
    }
    if (rna_mol_init(rmol) < 0) {
        rna_mol_free(rmol);
        return(NULL);
    }
    return(rmol);
}

void rna_mol_free(rna_mol_t *rmol){
    if (rmol == NULL) return;

    init_g_region(&rmol->loc);
    rmol->umi = 0;
    seq_base_l_free(&rmol->bl);
    seq_vac_l_free(&rmol->vl);
    seq_gene_l_free(&rmol->gl);
    rmol->n_reads = 0;
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

    cpy->umi = m->umi;
    cpy->loc = m->loc;

    if (seq_base_l_cpy(&cpy->bl, &m->bl) < 0){
        *ret = -1;
        rna_mol_dstry(cpy);
        return(NULL);
    }

    if (seq_vac_l_cpy(&cpy->vl, &m->vl) < 0){
        *ret = -1;
        rna_mol_dstry(cpy);
        return(NULL);
    }

    if (seq_gene_l_cpy(&cpy->gl, &m->gl) < 0){
        *ret = -1;
        rna_mol_dstry(cpy);
        return(NULL);
    }

    cpy->n_reads = m->n_reads;

    return(cpy);
}

int rna_mol_var_call(rna_mol_t *m, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual){
    if (m == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "rna_mol_var_call: argument is null");
    }

    int ret = seq_vac_l_call_var(&m->bl, &m->vl, gv, cmap, min_qual);
    if (ret < 0)
        err_msg(ret, 0, "rna_mol_var_call: failed to call variants");

    // free the bases
    seq_base_l_free(&m->bl);

    return(ret);
}

