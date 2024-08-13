
#ifndef ATAC_DATA_H
#define ATAC_DATA_H

/* 
 * Data structures and methods for scATAC-seq
 */

#include <stdlib.h>
#include "str_util.h"
#include "htslib/hts.h"
#include "htslib/khash.h"
#include "g_list.h"
#include "kavl.h"
#include "variants.h"
#include "counts.h"

// store query read name as hashed 64 bit integer
typedef uint64_t qshort;

#define kh_qname_hash_func(i) (kh_int64_hash_func((i)))
#define kh_qname_hash_equal(a, b) ((a) == (b))

static inline qshort qname2qshort(const char *s){
    qshort qs = (qshort)(str_int64_hash(s));
    return(qs);
}

// collision check
// KHASH_INIT(khcol, qshort, uint8_t, 1, kh_qname_hash_func, kh_qname_hash_equal);

/*! @typedef
 * @abstract Single read
 *
 * @field loc read region.
 * @field bl list of bases.
 */
typedef struct atac_read1_t {
    g_region loc;
    ml_t(seq_base_l) bl;
} atac_read1_t;

/* @brief Store a pair of ATAC reads belonging to the same fragment.
 *
 * ATAC reads belong to the same fragment if they have the same read name
 * and are paired.
 */
typedef struct atac_rd_pair_t {
    /* First atac read in pair. */
    atac_read1_t r1;

    /* Second atac read in pair. */
    atac_read1_t r2;

    /* @brief Number of reads in pair.
     *
     * Keeps track of how many reads were added to the struct.
     * Can be 0, 1, or 2.
     */
    uint8_t s;
} atac_rd_pair_t;

/*******************************************************************************
 * atac pairs bag
 ******************************************************************************/

#define key1_cmp(key1, key2) (((key2) < (key1) - ((key1) < (key2))))
mt_declare(atac_pair, qshort, atac_rd_pair_t, key1_cmp);

// (internal) node for rna_dups_t object
typedef struct atac_pair_node {
    qshort key; // UMI hash
    atac_rd_pair_t atac_pair;
    KAVL_HEAD(struct atac_pair_node) head;
} atac_pair_node;
#define qshort_cmp(p, q) (((q)->key < (p)->key) - ((p)->key < (q)->key))
KAVL_INIT2(bt_qr_pr, static, struct atac_pair_node, head, qshort_cmp);

int atac_pair_node_init(atac_pair_node *node);
void atac_pair_node_free(atac_pair_node *node);

// container of atac_rd_pair_t objects
typedef struct atac_pair_bag_t {
    atac_pair_node *atac_pairs;
} atac_pair_bag_t;

int atac_pair_bag_init(atac_pair_bag_t *bag);
void atac_pair_bag_free(atac_pair_bag_t *bag);
int atac_pair_bag_add_read(atac_pair_bag_t *bag, const atac_read1_t *ar,
        qshort qname);

typedef struct {
    kavl_itr_t(bt_qr_pr) itr;
    uint8_t next;
} atac_pair_bag_itr;

void atac_pair_bag_itr_first(atac_pair_bag_itr *itr, atac_pair_bag_t *bag);
// return 0 if nothing left and no value, 1 if itr has value after call
int atac_pair_bag_itr_next(atac_pair_bag_itr *itr);
int atac_pair_bag_itr_alive(atac_pair_bag_itr *itr);
qshort *atac_pair_bag_itr_key(atac_pair_bag_itr *itr);
atac_rd_pair_t *atac_pair_bag_itr_val(atac_pair_bag_itr *itr);

/* read pairs keyed by query name */
KHASH_INIT(khap, qshort, atac_rd_pair_t *, 1, kh_qname_hash_func, kh_qname_hash_equal);

/*******************************************************************************
 ******************************************************************************/

/*******************************************************************************
 * atac frags
 ******************************************************************************/

typedef struct atac_read2_t{
    atac_read1_t r1;
    atac_read1_t r2;
    uint32_t s; // number os supporting reads
} atac_read2_t;

// vector of atac_read2_t
mv_declare(ar2_vec, atac_read2_t);

/*! @typedef atac_dups_t
 * Store an array of atac_rd_pair_t objects/
 */
typedef struct atac_dups_t {
    struct _atac_dup1_t {
        atac_rd_pair_t rd; // read pair
        uint32_t n_rd; // number of (duplicate) reads per pair.
    } *dups;
    uint32_t size, m;
} atac_dups_t;

/*******************************************************************************
 * atac dups bag
 ******************************************************************************/

mt_declare(atac_dups, g_reg_pair, atac_dups_t, reg_pair_cmp);

// (internal) node for rna_dups_t object
typedef struct atac_dup_node {
    g_reg_pair key;
    atac_dups_t atac_dup;
    KAVL_HEAD(struct atac_dup_node) head;
} atac_dup_node;
#define adn_cmp(p, q) reg_pair_cmp( (p)->key, (q)->key )
KAVL_INIT2(bt_rg_dp, static, struct atac_dup_node, head, adn_cmp);

int atac_dup_node_init(atac_dup_node *node);
void atac_dup_node_free(atac_dup_node *node);

// container of atac_dups_t objects
typedef struct atac_dup_bag_t {
    atac_dup_node *atac_dups;
} atac_dup_bag_t;

int atac_dup_bag_init(atac_dup_bag_t *bag);
void atac_dup_bag_free(atac_dup_bag_t *bag);
int atac_dup_bag_add_read(atac_dup_bag_t *bag, const atac_rd_pair_t *ap,
        int skip_chim);

typedef struct {
    kavl_itr_t(bt_rg_dp) itr;
    uint8_t next;
} atac_dup_bag_itr;

void atac_dup_bag_itr_first(atac_dup_bag_itr *itr, atac_dup_bag_t *bag);
// return 0 if nothing left and no value, 1 if itr has value after call
int atac_dup_bag_itr_next(atac_dup_bag_itr *itr);
int atac_dup_bag_itr_alive(atac_dup_bag_itr *itr);
g_reg_pair *atac_dup_bag_itr_key(atac_dup_bag_itr *itr);
atac_dups_t *atac_dup_bag_itr_val(atac_dup_bag_itr *itr);

/* duplicates keyed by region */
KHASH_INIT(khad, g_reg_pair, atac_dups_t *, 1, kh_reg_pair_hash, kh_reg_pair_equal);

/*******************************************************************************
 * atac frags
 ******************************************************************************/

/*! @typedef atac_frag_t
 */
typedef struct atac_frag_t {
    ml_t(seq_base_l) bl; // bases
    ml_t(seq_vac_l) vl; // variant calls
    mv_t(int_vec) pks; // peaks
    uint32_t s; // number of supporting reads
} atac_frag_t;

/*******************************************************************************
 * atac frags bag
 ******************************************************************************/

mt_declare(atac_frags, g_reg_pair, atac_frag_t, reg_pair_cmp);

// (internal) node for rna_frags_t object
typedef struct atac_frag_node {
    g_reg_pair key;
    atac_frag_t atac_frag;
    KAVL_HEAD(struct atac_frag_node) head;
} atac_frag_node;
#define adn_cmp(p, q) reg_pair_cmp( (p)->key, (q)->key )
KAVL_INIT2(bt_rg_fr, static, struct atac_frag_node, head, adn_cmp);

int atac_frag_node_init(atac_frag_node *node);
void atac_frag_node_free(atac_frag_node *node);

// container of atac_frags_t objects
typedef struct atac_frag_bag_t {
    atac_frag_node *atac_frags;
} atac_frag_bag_t;

int atac_frag_bag_init(atac_frag_bag_t *bag);
void atac_frag_bag_free(atac_frag_bag_t *bag);
int atac_frag_bag_add(atac_frag_bag_t *bag, const atac_frag_t *af,
        g_reg_pair reg);

typedef struct {
    kavl_itr_t(bt_rg_fr) itr;
    uint8_t next;
} atac_frag_bag_itr;

void atac_frag_bag_itr_first(atac_frag_bag_itr *itr, atac_frag_bag_t *bag);
// return 0 if nothing left and no value, 1 if itr has value after call
int atac_frag_bag_itr_next(atac_frag_bag_itr *itr);
int atac_frag_bag_itr_alive(atac_frag_bag_itr *itr);
g_reg_pair *atac_frag_bag_itr_key(atac_frag_bag_itr *itr);
atac_frag_t *atac_frag_bag_itr_val(atac_frag_bag_itr *itr);

/* atac fragments. Key is region, value is pointer to atac_frag_t */
KHASH_INIT(khaf, g_reg_pair, atac_frag_t *, 1, kh_reg_pair_hash, kh_reg_pair_equal);

#define mlaf_lt(p, q) -1
ml_declare(mlaf, atac_frag_t *, mlaf_lt);

/*******************************************************************************
 * atac_read1_t
 ******************************************************************************/

/* Set members of atac read to 0.
 */
void atac_read_set0(atac_read1_t *ar);

/* Initialize atac read
 *
 * Allocates and zeroes the members. Object must be 
 * freed by caller.
 *
 * @return Pointer to dynamically allocated object, or NULL on failure.
 */
atac_read1_t *atac_read_init();

/* Free an atac_read object's members.
 * This does not destroy the object pointed to by @p ar, but frees 
 * the memory of its members.
 */
void atac_read_free(atac_read1_t *ar);

/* Free an atac_read object and its members.
 * The object pointed to by @p ar is freed.
 */
void atac_read_dstry(atac_read1_t *ar);

/* Create a deep copy of an atac read object.
 *
 * The object returned must be freed by the caller.
 * If @p r is NULL, returns NULL.
 *
 * The int pointed to by ret is modified after the function call. If an 
 * error occured, it is modified to -1. Otherwise, it is set to 0.
 *
 * @param r Pointer to atac read object.
 * @param ret Pointer to int to store return status.
 */
atac_read1_t *atac_read1_dup(const atac_read1_t *r, int *ret);
int atac_read1_cpy(atac_read1_t *dest, const atac_read1_t *src);

/* Compare if two atac reads are equal.
 *
 * Two atac reads are equal if they have the same genomic position and 
 * the same base list (base call quality is not considered).
 *
 * @return 1 if the reads are equal, 0 if they are not equal
 */
int atac_read_equal(atac_read1_t r1, atac_read1_t r2);

/* Add base to atac_read1_t object.
 *
 * The seq_base_t object is copied and added to the seq_blist_t in @p r.
 *
 * @return 0 on success, -1 on error.
 */
int atac_read1_add_base(atac_read1_t *r, seq_base_t base);

/*******************************************************************************
 * atac_rd_pair_t
 ******************************************************************************/

/* Set 0 for read pair
 */
void atac_rd_pair_set0(atac_rd_pair_t *rp);

/* Initialize an atac_rd_pair_t object.
 *
 * @return Pointer to allocated object, NULL on error.
 */
atac_rd_pair_t *atac_rd_pair_init();

/* Free underlying memory, but not object. */
void atac_rd_pair_free(atac_rd_pair_t *rp);
/* Free underlying memory and the object. */
void atac_rd_pair_dstry(atac_rd_pair_t *rp);

/* @brief Duplicate an atac_rd_pair_t object.
 *
 * Create a deep copy of the atac_rd_pair_t object pointed to by rp
* and return a pointer to the new object.
 *
 * @param rp Pointer to atac_rd_pair_t object to copy.
 * @param ret Pointer to int to store return status.
 * 
 * @return - NULL if rp is NULL or if an error occurred. ret is set to -1.
 *         - Otherwise, !NULL pointer to a dynamically allocated atac_rd_pair_t object.
 *           ret is set to 0.
 *         
 * @note The return value must be freed by the caller.
 */
atac_rd_pair_t *atac_rd_pair_dup(const atac_rd_pair_t *rp, int *ret);

/**
 * @brief Copy the contents of an atac_rd_pair_t object to another.
 *
 * This function creates a deep copy of the contents of the `src` atac_rd_pair_t
 * object and stores it in the `dest` atac_rd_pair_t object. The `dest` object
 * must be properly initialized before calling this function.
 *
 * @param dest Pointer to the destination atac_rd_pair_t object where the copy
 *             will be stored.
 * @param src Pointer to the source atac_rd_pair_t object to be copied.
 *
 * @return 0 on success, -1 on error.
 *
 * @retval 0 The copy operation was successful.
 * @retval -1 An error occurred during the copy operation. This can happen if
 *            either `dest` or `src` is NULL, or if an error occurred while
 *            copying the individual atac_read1_t objects.
 *
 * @note This function assumes that `dest` has been properly initialized before
 *       calling it. If `dest` contains any existing data, it will be overwritten
 *       by the copy operation.
 *
 * @see atac_read1_cpy
 */
int atac_rd_pair_cpy(atac_rd_pair_t *dest, const atac_rd_pair_t *src);

/*
 * @brief Add an ATAC read to an atac_rd_pair_t struct.
 *
 * This function adds the given ATAC read to the atac_rd_pair_t struct. If the
 * struct is empty, the read is added to r1. If the struct contains one read,
 * the read is added to r2, and the order of r1 and r2 is adjusted so that r1
 * comes before r2 based on their locations. If the struct already contains two
 * reads, an error is returned.
 *
 * The function makes a copy of the input read using atac_read1_dup before
 * adding it to the struct.
 *
 * @param rp Pointer to the atac_rd_pair_t struct to add the read to.
 * @param ar Pointer to the atac_read1_t struct to add.
 *
 * @return 0 on success, -1 on error.
 */
int atac_rd_pair_add_read(atac_rd_pair_t *rp, const atac_read1_t *ar);

/* Test if two atac read pairs are equal.
 *
 * Expects both arguments are not null.
 * The pairs are equal if they have the same number of reads (0, 1, or 2), 
 * and each of the reads are equal.
 * Neither of the parameters can be null.
 *
 * @return 1 for equality, 0 for inequality, -1  for error.
 */
int atac_rd_pair_equal(const atac_rd_pair_t *rp1, const atac_rd_pair_t *rp2);

/* Set the highest base quality in a read of the two reads 
 *
 * Expects `rp` and `cmp` have the same bases, it is an error otherwise.
 * Compares the quality score at each base between the two reads. If the 
 * base quality is higher in `cmp`, then set the quality in `rp` 
 * to that of `cmp`. Calls the function `seq_base_l_match_qual` to match
 * the qualities of bases for each read pair.
 *
 * @param rp Pointer to atac_rd_pair_t to modify the qualities of.
 * @param cmp Pointer to atac_rd_pair_t to compare the qualities.
 * @return 0 on success, -1 on error.
 *
 * @see seq_base_l_match_qual
 */
int atac_rd_pair_match_qual(atac_rd_pair_t *rp, const atac_rd_pair_t *cmp);

/*******************************************************************************
 * atac_dups_t
 ******************************************************************************/

/* Allocate and initialize and atac dups object.
 * @return Allocated and initialized object, or NULL on error.
 */
int atac_dups_init(atac_dups_t *d);
atac_dups_t *atac_dups_alloc();

/* Free the memory an atac dups object and its memers.
 *
 * All read pairs are freed.
 */
/* Free the underlying memory but not the object. */
void atac_dups_free(atac_dups_t *d);
/* Free the underlying memory and the object. */
void atac_dups_dstry(atac_dups_t *d);

/* Add a read pair to an atac dups object.
 *
 * The read pair is added to atac dups. If the read pair was 
 * found, increment @f rds_per_dup in @p d. If the read pair was not found, 
 * then copy the read pair, add to @f dups, and set the rds_per_dup to 1.
 * Sets the region from @p rp.
 *
 * Expects non-null arguments.
 * w
 * @param d A pointer to the atac_dups_t object to add the read pair to.
 * @param rp A pointer to the read pair object to add.
 * @return 0 on success, -1 on error or if memory allocation fails.
 */
int atac_dups_add_pair(atac_dups_t *d, const atac_rd_pair_t *rp);

/*******************************************************************************
 * atac_frag_t
 ******************************************************************************/

/* Allocate and initialize an empty atac frag object. */
int atac_frag_init(atac_frag_t *f);
atac_frag_t *atac_frag_alloc(); 

/* Free underlying memory and object itself. */
void atac_frag_free(atac_frag_t *f);
void atac_frag_dstry(atac_frag_t *f);

/*
 * @note returned object must be freed by caller.
 */
atac_frag_t *atac_dups_dedup(atac_dups_t *dups, int *ret);

/* call variants from atac fragments
 * Variants are called and placed in @f vacs. The seq_blist_t 
 * in @f bases is destroyed after this call.
 *
 * Expects non-null arguments
 * @param f Pointer to atac_frag_t object.
 * @param gv Pointer to g_var_t object.
 * @param cmap str_map between integer IDs and chromosome names.
 * @param min_qual minimum read phred base quality.
 *
 * @return -1 on error, the number of variants called on success.
 */
int atac_frag_var_call(atac_frag_t *f, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual);

/* Call peaks for a fragment.
 * The pairs must be aligned to the same chromosome.
 *
 * @param f The fragment to get overlapping peaks for.
 * @param reg The g_reg_pair of the fragment.
 * @param pks Pointer to the iregs_t which contains the peaks.
 * @param cmap Pointer to the contig map.
 * @return The number of overlapping peaks added, or -1 on error.
 */
int atac_frag_peak_call(atac_frag_t *f, g_reg_pair reg, iregs_t *pks, 
        str_map *cmap);

#endif // ATAC_DATA_H
