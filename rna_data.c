
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
        return err_msg(-1, 0, "rna_read1_equal: r1 or r2 are NULL");

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
        return err_msg(-1, 0, "rna_read1_add_gene: arguments are NULL");

    if (seq_glist_add_gene(&r->genes, gene) < 0) return(-1);

    return(0);
}

int rna_read1_add_base(rna_read1_t *r, const seq_base_t *base){
    if (r == NULL || base == NULL)
        return err_msg(-1, 0, "rna_read1_add_base: arguments are NULL");

    if (blist_add_base(&r->bases, base, 1, 1) < 0) return(-1);

    return(0);
}

int rna_read1_match_qual(rna_read1_t *r, const rna_read1_t *cmp){
    if (r == NULL || cmp == NULL)
        return err_msg(-1, 0, "rna_read1_match_qual: arguments are NULL");

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
        return err_msg(-1, 0, "rna_dups_add_read: arguments are NULL");

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
        return err_msg(-1, 0, "rna_mol_dedup: mol is NULL");

    if (mol->dups == NULL)
        return err_msg(-1, 0, "rna_mol_dedup: mol->dups is NULL");

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

rna_mol_t *rna_dups_dedup(rna_dups_t *rd, int *ret){
    *ret = 0;
    if (rd == NULL){
        *ret = err_msg(-1, 0, "rna_dups_dedup: rd is NULL");
        return(NULL);
    }

    // if no dups, return NULL successfully
    if (rd->size == 0)
        return(NULL);

    int ix_best = 0;
    size_t max_c = 0; // store read count
    size_t rd_w_max = 0; // number of read pairs with max_c read counts

    int i;
    for (i = 0; i < rd->size; ++i){
        // check for bugs
        if (rd->rds_per_dup[i] == 0){
            *ret = err_msg(-1, 0, "rna_dups_dedup: a dup read has 0 "
                    "supporting reads, there is a bug");
            return(NULL);
        }

        if (rd->rds_per_dup[i] == max_c){
            ++rd_w_max;
        } else if (rd->rds_per_dup[i] > max_c){
            max_c = rd->rds_per_dup[i];
            ix_best = i;
            rd_w_max = 1;
        }
    }

    // check for bugs
    if (max_c == 0 || rd_w_max == 0){
        *ret = err_msg(-1, 0, "atac_rna_dedup: no best duplicate found "
                ", there is a bug");
        return(NULL);
    }

    // initialize molecule
    rna_mol_t *m = rna_mol_alloc();
    rna_read1_t rb = rd->dups[ix_best];

    m->loc = rb.loc;
    
    int cret;

    seq_blist_t *b_cpy = copy_seq_blist(&rb.bases, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    m->bases = *b_cpy;
    free(b_cpy);

    seq_glist_t *g_cpy = seq_glist_cpy(&rb.genes, &cret);
    if (cret < 0){
        *ret = -1;
        return(NULL);
    }
    m->genes = *g_cpy;
    free(g_cpy);

    m->n_reads = rd->rds_per_dup[ix_best];

    return(m);
}

int rna_mol_var_call(rna_mol_t *m, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual){
    if (m == NULL || gv == NULL || cmap == NULL){
        return err_msg(-1, 0, "rna_mol_var_call: arguments must not be NULL");
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

void bc_rna_init(bc_rna_t *br){
    if (br == NULL) return;
    br->dups = kh_init(khrdn);
    br->mols = kh_init(khrmn);
}

bc_rna_t *bc_rna_alloc(){
    bc_rna_t *br = (bc_rna_t *)calloc(1, sizeof(bc_rna_t));
    if (br == NULL){
        err_msg(-1, 0, "bc_rna_alloc: %s", strerror(errno));
        return(NULL);
    }
    bc_rna_init(br);
    return(br);
}

void bc_rna_free_dups(bc_rna_t *br){
    if (br == NULL || br->dups == NULL) return;
    khint_t k;
    for (k = kh_begin(br->dups); k != kh_end(br->dups); ++k){
        if (!kh_exist(br->dups, k)) continue;
        rna_dups_t *d = kh_val(br->dups, k);
        rna_dups_dstry(d);
        char *name = kh_key(br->dups, k);
        free(name);
    }
    kh_destroy(khrdn, br->dups);
    br->dups = NULL;
}

void bc_rna_free_mols(bc_rna_t *br){
    if (br == NULL || br->mols == NULL) return;
    khint_t k;
    for (k = kh_begin(br->mols); k != kh_end(br->mols); ++k){
        if (!kh_exist(br->mols, k)) continue;
        rna_mol_t *m = kh_val(br->mols, k);
        rna_mol_dstry(m);
        char *name = kh_key(br->mols, k);
        free(name);
    }
    kh_destroy(khrmn, br->mols);
    br->mols = NULL;
}

void bc_rna_free(bc_rna_t *br){
    if (br == NULL) return;

    bc_rna_free_dups(br);
    bc_rna_free_mols(br);
}

void bc_rna_dstry(bc_rna_t *br){
    if (br == NULL) return;
    bc_rna_free(br);
    free(br);
}

int bc_rna_add_read(bc_rna_t *br, const rna_read1_t *r, const char *name){
    if (br == NULL || r == NULL)
        return err_msg(-1, 0, "bc_rna_add_read: arguments are NULL");
    if (name == NULL)
        return err_msg(-1, 0, "bc_rna_add_read: name is NULL");

    khash_t(khrdn) *dups = br->dups;
    khint_t k = kh_get(khrdn, dups, (char *)name); // add read by name (UMI barcode)

    // if read isn't present, allocate and add
    int kret = 0;
    if (k == kh_end(dups)){
        char *name_cpy = strdup(name);
        if (name_cpy == NULL)
            return err_msg(-1 , 0, "bc_rna_add_read: %s", strerror(errno));

        k = kh_put(khrdn, dups, name_cpy, &kret);
        if (kret < 0)
            return err_msg(-1, 0, "bc_rna_add_read: failed to add name to dups hash table");

        if ( (kh_val(dups, k) = rna_dups_alloc()) == NULL )
            return(-1);
    }

    rna_dups_t *dup = kh_val(dups, k);
    if (rna_dups_add_read(dup, r) < 0) return(-1);

    return(0);
}

int bc_rna_add_mol(bc_rna_t *br, const rna_mol_t *m, const char *name){
    if (br == NULL || m == NULL || name == NULL)
        return err_msg(-1, 0, "bc_rna_add_mol: arguments are NULL");

    khash_t(khrmn) *mols = br->mols;
    khint_t k;

    // copy m
    int cret;
    rna_mol_t *m_cpy = rna_mol_cpy(m, &cret);
    if (cret < 0)
        return(-1);

    char *name_cpy = strdup(name);
    if (name_cpy == NULL)
            return err_msg(-1, 0, "bc_rna_add_mol: %s", strerror(errno));
    k = kh_put(khrmn, mols, name_cpy, &cret); // use memory from m
    if (cret < 0)
        return err_msg(-1, 0, "bc_rna_add_mol: failed to add UMI to mol "
                "hash table");
    // molecule should not be present
    if (cret == 0)
        return err_msg(-1, 0, "bc_rna_add_mol: another molecule %s found", name_cpy);

    kh_val(mols, k) = m_cpy;

    return(0);
}

int bc_rna_dedup(bc_rna_t *br){
    if (br == NULL)
        return err_msg(-1, 0, "bc_rna_dedup: br is NULL");

    khint_t kd;
    for (kd = kh_begin(br->dups); kd != kh_end(br->dups); ++kd){
        if (!kh_exist(br->dups, kd)) continue;

        rna_dups_t *dup = kh_val(br->dups, kd);
        if (dup == NULL)
            return err_msg(-1, 0, "bc_rna_dedup: a duplicate is NULL");

        char *name = kh_key(br->dups, kd);
        if (name == NULL)
            return err_msg(-1, 0, "bc_rna_dedup: read name is NULL");

        // get best RNA molecule by majority
        int dret;
        rna_mol_t *m = rna_dups_dedup(dup, &dret);
        if (dret < 0) return(-1);

        if (bc_rna_add_mol(br, m, name) < 0) return(-1);

        rna_mol_dstry(m); // destroy since we added a copy
    }
    return(0);
}

int bc_rna_var_call(bc_rna_t *br, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual){
    if (br == NULL)
        return err_msg(-1, 0, "bc_rna_var_call: br is NULL");

    if (br->mols == NULL)
        return err_msg(-1, 0, "bc_rna_var_call: mols is NULL");

    int n_add = 0;
    khint_t km;
    for (km = kh_begin(br->mols); km != kh_end(br->mols); ++km){
        if (!kh_exist(br->mols, km)) continue;

        rna_mol_t *mol = kh_val(br->mols, km);
        if (mol == NULL)
            return err_msg(-1, 0, "bc_rna_var_call: a molecule is NULL");

        int a = rna_mol_var_call(mol, gv, cmap, min_qual);
        if (a < 0)
            return(-1);
        n_add += a;
    }
    return(n_add);
}

/*******************************************************************************
 * bam_rna_t
 ******************************************************************************/

void bam_rna_init(bam_rna_t *bam_r){
    if (bam_r == NULL) return;

    bam_r->bc_rna = kh_init(khrbc);
}

bam_rna_t *bam_rna_alloc(){
    bam_rna_t *bam_r = (bam_rna_t *)calloc(1, sizeof(bam_rna_t));
    if (bam_r == NULL){
        err_msg(-1, 0, "bam_rna_alloc: %s", strerror(errno));
        return(NULL);
    }

    bam_rna_init(bam_r);
    return(bam_r);
}

void bam_rna_free_dups(bam_rna_t *bam_r){
    if (bam_r == NULL) return;

    khint_t k;
    for (k = kh_begin(bam_r->bc_rna); k != kh_end(bam_r->bc_rna); ++k){
        if (!kh_exist(bam_r->bc_rna, k)) continue;
        bc_rna_t *bcr = kh_val(bam_r->bc_rna, k);
        bc_rna_free_dups(bcr);
    }
}

void bam_rna_free(bam_rna_t *bam_r){
    if (bam_r == NULL) return;

    khint_t k;
    for (k = kh_begin(bam_r->bc_rna); k != kh_end(bam_r->bc_rna); ++k){
        if (!kh_exist(bam_r->bc_rna, k)) continue;
        bc_rna_t *bcr = kh_val(bam_r->bc_rna, k);
        bc_rna_dstry(bcr);
        char *bc = kh_key(bam_r->bc_rna, k);
        free(bc);
        kh_val(bam_r->bc_rna, k) = NULL;

    }
    kh_destroy(khrbc, bam_r->bc_rna);
}

void bam_rna_dstry(bam_rna_t *bam_r){
    if (bam_r == NULL) return;

    bam_rna_free(bam_r);
    free(bam_r);
}

int bam_rna_add_read(bam_rna_t *bam_r, const char *bc, const rna_read1_t *r, 
        const char *name){
    if (bam_r == NULL || r == NULL)
        return err_msg(-1, 0, "bam_rna_add_read: arguments are NULL");

    int ret;
    khint_t kbc;
    kbc = kh_get(khrbc, bam_r->bc_rna, (char *)bc);
    // if barcode isn't present, add bc_rna
    if (kbc == kh_end(bam_r->bc_rna)){
        char *bc_cpy = strdup(bc);
        kbc = kh_put(khrbc, bam_r->bc_rna, bc_cpy, &ret);
        if (ret < 0)
            return err_msg(-1, 0, "bam_rna_add_read: failed to add read to bcs");

        if ( (kh_val(bam_r->bc_rna, kbc) = bc_rna_alloc()) == NULL )
            return(-1);
    }
    bc_rna_t *bcr = kh_val(bam_r->bc_rna, kbc);
    ret = bc_rna_add_read(bcr, r, name);
    if (ret < 0)
        return err_msg(-1, 0, "bam_rna_add_read: failed to add read to bam");
    return(0);
}

int bam_rna_dedup(bam_rna_t *bam_r){
    if (bam_r == NULL)
        return err_msg(-1, 0, "bam_rna_dedup: arguments are NULL");

    khint_t kbc;
    for (kbc = kh_begin(bam_r->bc_rna); kbc != kh_end(bam_r->bc_rna); ++kbc){
        if (!kh_exist(bam_r->bc_rna, kbc)) continue;

        bc_rna_t *bc_rna = kh_val(bam_r->bc_rna, kbc);
        if (bc_rna == NULL)
            return err_msg(-1, 0, "bam_rna_dedup: a bc_rna object is NULL, "
                    "there is a bug");

        // bc_rna_dedup
        if (bc_rna_dedup(bc_rna) < 0)
            return(-1);
    }
    return(0);
}

int bam_rna_var_call(bam_rna_t *bam_r, GenomeVar *gv, contig_map *cmap, 
        uint8_t min_qual){
    if (bam_r == NULL)
        return err_msg(-1, 0, "bam_rna_var_call: argument bam_r is NULL");

    int ret, n_add = 0;
    khint_t kbc;
    for (kbc = kh_begin(bam_r->bc_rna); kbc != kh_end(bam_r->bc_rna); ++kbc){
        if (!kh_exist(bam_r->bc_rna, kbc)) continue;
        bc_rna_t *bca = kh_val(bam_r->bc_rna, kbc);
        if (bca == NULL)
            return err_msg(-1, 0, "bam_rna_var_call: barcode is NULL");

        ret = bc_rna_var_call(bca, gv, cmap, min_qual);
        if (ret < 0) return(-1);
        n_add += ret;
    }
    return(n_add);
}

/*******************************************************************************
 * miscellaneous
 ******************************************************************************/

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

void print_bam_rna(bam_rna_t *b, Annotation *anno, bcf_hdr_t *vcf_hdr, 
        GenomeVar *gv){
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
                Var *var = gv_vari(gv, va->vix);
                fprintf(stdout, "\tvar:");
                fprintf(stdout, "\t%s", var->b->d.id);
                fprintf(stdout, "\t%"PRIhts_pos"", var->b->pos);
                fprintf(stdout, "\t%u\n", va->allele);
                va = va->next;
            }
        }
    }
}

