
#ifndef BAM_DAT_H
#define BAM_DAT_H

#include <stdlib.h>
#include <pthread.h>
#include "rna_data.h"
#include "atac_data.h"
#include "variants.h"
#include "region.h"
#include "bc_stats.h"
#include "htslib/sam.h"

/*******************************************************************************
 * bc_data_t
 ******************************************************************************/

/* @brief Store barcode data
 */
typedef struct {
    // RNA
    khash_t(kh_rd) *rna_dups;
    ml_t(ml_rm) *rna_mols;
    // kbtree_t(kb_rm) *rna_mols;

    // ATAC
    khash_t(kh_arp) *atac_prs;
    khash_t(kh_ad) *atac_dups;
    ml_t(ml_af) *atac_frgs;

    bc_stats_t *bc_stats; // barcode stats/counts

    pthread_mutex_t bc_lock; // multitheading lock
} bc_data_t;

// hash table for bc_data_t
// key is barcode string, value is pointer to bc_data_t
KHASH_INIT(kh_bc_dat, char *, bc_data_t *, 1, kh_str_hash_func, kh_str_hash_equal);

/* Initialize empty bc_data_t struct.
 * Hash tables are initialized.
 * Return null on error.
 */
bc_data_t *bc_data_init();

/* @brief Free memory underlying bcdat->rna_dupls.
 */
void bc_data_free_rna_dupls(bc_data_t *bcdat);

/* @brief Free memory underlying bcdat->rna_mlcls.
 */
void bc_data_free_rna_mlcls(bc_data_t *bcdat);

/* @brief Free memory underlying bcdat->atac_pairs.
 */
void bc_data_free_atac_pairs(bc_data_t *bcdat);

/* @brief Free memory underlying bcdat->atac_dupls.
 */
void bc_data_free_atac_dupls(bc_data_t *bcdat);

/* @brief Free memory underlying bcdat->atac_frags.
 */
void bc_data_free_atac_frags(bc_data_t *bcdat);

/* @brief Free memory underlying reads of bcdat.
 */
void bc_data_free_reads(bc_data_t *bcdat);

/* Destroy and free all memory of bcdat
 */
void bc_data_dstry(bc_data_t *bcdat);
/* Free pairs in bcdat struct
 */
void bc_data_free_atac_pairs(bc_data_t *bcdat);

int bc_data_fill_list(bc_data_t *bcdat);

/*******************************************************************************
 * BC RNA
 ******************************************************************************/

/**
 * @brief Add an RNA read to the specified bc_data_t structure.
 *
 * This function adds an RNA read to the rna_dupls tree in the given
 * bc_data_t structure. It initializes the required structures, and performs
 * necessary checks.
 *
 * @param bcdat Pointer to the bc_data_t structure where the RNA read 
 * will be added.
 * @param r Pointer to the rna_read1_t structure representing the RNA read.
 * @param umih The UMI (Unique Molecular Identifier) associated with the RNA read.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bc_data_rna_add_read: argument is null" if either bcdat or r is NULL.
 * - "bc_data_rna_add_read: read region is null" if read is uninitialized.
 * - "bc_data_rna_add_read: failed to initialize rna_dups" if initialization 
 * of rna_dup fails.
 * - "bc_data_rna_add_read: failed to add rna_dups to bcdat" if adding rna_dups 
 * to bcdat fails.
 * - "bc_data_rna_add_read: failed to add read to rna_dups" if adding the read 
 * to rna_dups fails.
 *
 * @note This function is not thread-safe.
 */
int bc_data_rna_add_read(bc_data_t *bcdat, const rna_read1_t *r, umishort umih);

/**
 * @brief De-duplicate RNA reads in the specified bc_data_t structure.
 *
 * This function traverses the `rna_dupls` tree in the given bc_data_t structure,
 * performing deduplication on the RNA reads. Deduplicated molecules are then
 * added to the `rna_mlcls` tree. The bcdat->rna_dupls values are freed after deduplication.
 *
 * @param bcdat Pointer to the bc_data_t structure containing RNA reads to be deduplicated.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bc_data_rna_dedup: bcdat is null" if the input argument is NULL.
 * - "bc_data_rna_dedup: failed to initialize iterator" if the iterator initialization fails.
 * - "bc_data_rna_dedup: failed to deduplicate reads" if deduplication fails.
 * - "bc_data_rna_dedup: failed to add rna_mol to bcdat" if adding deduplicated molecules to the `rna_mlcls` tree fails.
 * - "bc_data_rna_dedup: duplicate RNA BC-UMI found" if a duplicate molecule is encountered while adding to the `rna_mlcls` tree.
 *
 * @note This function is not thread-safe.
 */
int bc_data_rna_dedup(bc_data_t *bcdat);

/**
 * @brief Call RNA variants in the specified bc_data_t structure.
 *
 * This function traverses the `rna_mlcls` tree in the provided bc_data_t structure,
 * and performs variant calling on each RNA molecule. The results are added to the
 * provided g_var_t structure.
 *
 * @param bcdat Pointer to the bc_data_t structure containing RNA molecules.
 * @param gv Pointer to the g_var_t structure where variant calls will be stored.
 * @param cmap Pointer to the str_map structure used for chromosome name to index mapping.
 * @param min_qual Minimum base phred quality in read to consider for calling a variant allele.
 *
 * @return The number of variants added on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bc_data_rna_var_call: argument is null" if bcdat, gv, or cmap is NULL.
 * - "bc_data_rna_var_call: failed to initialize iterator" if the iterator initialization fails.
 * - "bc_data_rna_var_call: failed to call variant" if variant calling on an RNA molecule fails.
 *
 * @note This function is not thread-safe.
 */
int bc_data_rna_var_call(bc_data_t *bcdat, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual);

/*******************************************************************************
 * BC ATAC
 ******************************************************************************/

/**
 * @brief Add an ATAC read to the specified bc_data_t structure.
 *
 * This function adds an ATAC read to the `atac_pairs` tree in the given
 * bc_data_t structure. It initializes the required structures, and performs
 * necessary checks.
 *
 * @param bcdat Pointer to the bc_data_t structure where the ATAC read will be added.
 * @param ar Pointer to the atac_read1_t structure representing the ATAC read.
 * @param qname The query name associated with the ATAC read.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bc_data_atac_add_read: argument is null" if either bcdat or ar is NULL.
 * - "bc_data_atac_add_read: failed to add atac pair to bcdat" if adding atac pair to `atac_pairs` tree fails.
 * - "bc_data_atac_add_read: failed to add read to atac pair" if adding the read to the atac pair fails.
 *
 * @note This function is not thread-safe.
 */
int bc_data_atac_add_read(bc_data_t *bcdat, const atac_read1_t *ar, qshort qname);

int bc_data_atac_add_pair(bc_data_t *bcdat, const atac_rd_pair_t *ap);

/**
 * @brief Identify and add ATAC read duplicates in the specified bc_data_t structure.
 *
 * This function traverses the `atac_pairs` tree in the given bc_data_t structure,
 * identifies duplicate read pairs, and adds them to the `atac_dupls` tree.
 *
 * @param bcdat Pointer to the bc_data_t structure containing ATAC read pairs.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bc_data_atac_get_dups: argument is null" if the bcdat argument is NULL.
 * - "bc_data_atac_get_dups: failed to initialize iterator" if the iterator initialization fails.
 * - "bc_data_atac_get_dups: failed to check if pair is chimeric" if the chimeric check fails.
 * - "bc_data_atac_get_dups: failed to add dups to bcdat" if adding duplicates to `atac_dupls` tree fails.
 * - "bc_data_atac_get_dups: failed to add read to dups" if adding the read pair to the duplicates set fails.
 *
 * @note This function is not thread-safe.
 */
int bc_data_atac_get_dups(bc_data_t *bcdat);

/**
 * @brief De-duplicate ATAC read pairs in the specified bc_data_t structure.
 *
 * This function traverses the `atac_dupls` tree in the given bc_data_t structure,
 * performs deduplication on the ATAC read pairs, and adds the deduplicated fragments
 * to the `atac_frags` tree.
 *
 * @param bcdat Pointer to the bc_data_t structure containing ATAC duplicates.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bc_data_atac_dedup: bcdat is null" if the bcdat argument is NULL.
 * - "bc_data_atac_dedup: failed to initialize iterator" if the iterator initialization fails.
 * - "bc_data_atac_dedup: deduplication failed" if deduplication fails.
 * - "bc_data_atac_dedup: failed to add frag to bcdat" if adding deduplicated fragment to `atac_frags` tree fails.
 *
 * @note This function is not thread-safe.
 */
int bc_data_atac_dedup(bc_data_t *bcdat);


/**
 * @brief Call ATAC variants in the specified bc_data_t structure.
 *
 * This function traverses the `atac_frags` tree in the provided bc_data_t structure,
 * and performs variant calling on each ATAC fragment. The results are added to the
 * provided g_var_t structure.
 *
 * @param bcdat Pointer to the bc_data_t structure containing ATAC fragments.
 * @param gv Pointer to the g_var_t structure where variant calls will be stored.
 * @param cmap Pointer to the str_map structure used for chromosome name to index mapping.
 * @param min_qual Minimum base phred quality in read to consider for calling a variant allele.
 *
 * @return The number of variants added on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bc_data_atac_var_call: argument is null" if bcdat, gv, or cmap is NULL.
 * - "bc_data_atac_var_call: failed to initialize iterator" if the iterator initialization fails.
 * - "bc_data_atac_var_call: failed to call variant" if variant calling on an ATAC fragment fails.
 *
 * @note This function is not thread-safe.
 */
int bc_data_atac_var_call(bc_data_t *bcdat, g_var_t *gv, str_map *cmap, 
        uint8_t min_qual);

/**
 * @brief Call ATAC peaks in the specified bc_data_t structure.
 *
 * This function traverses the `atac_frags` tree in the provided bc_data_t structure,
 * and performs peak calling on each ATAC fragment. The results are added to the
 * provided iregs_t structure.
 *
 * @param bcdat Pointer to the bc_data_t structure containing ATAC fragments.
 * @param pks Pointer to the iregs_t structure where peak calls will be stored.
 * @param cmap Pointer to the str_map structure used for chromosome name to index mapping.
 *
 * @return The number of peaks added on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bc_data_atac_peak_call: argument is null" if bcdat, pks, or cmap is NULL.
 * - "bc_data_atac_peak_call: failed to initialize iterator" if the iterator initialization fails.
 * - "bc_data_atac_peak_call: failed to call peaks" if peak calling on an ATAC fragment fails.
 *
 * @note This function is not thread-safe.
 */
int bc_data_atac_peak_call(bc_data_t *bcdat, iregs_t *pks, str_map *cmap);

/*******************************************************************************
 * bam_data_t
 ******************************************************************************/

/*! @typedef
 * @abstract Structure to hold the atac and rna pileup
 */
typedef struct {
    khash_t(kh_bc_dat) *bc_data; // hash table for barcode data

    // flags
    uint8_t has_rna;
    uint8_t has_atac;
    uint8_t has_stats;

    str_map *bcs; // barcodes

    pthread_mutex_t bam_lock; // multitheading lock
} bam_data_t;

/**
 * @brief Initialize a bam_data_t structure.
 *
 * This function allocates and initializes a bam_data_t structure, including its
 * hash tables, string map, and mutex. It ensures that all necessary components
 * are properly allocated and initialized, and handles error cases by cleaning
 * up any allocated resources before returning.
 *
 * @return Pointer to the initialized bam_data_t structure, or NULL on error.
 *
 * Possible error messages:
 * - "bam_data_init: [system error message]" if memory allocation fails.
 * - "bam_data_init: failed to initialize mutex errno=[error code]" if mutex initialization fails.
 *
 * @note This function ensures that allocated resources are properly freed in case
 * of initialization failure to prevent memory leaks.
 */
bam_data_t *bam_data_init();

/**
 * @brief Free resources associated with the bc_data hash table in bam_data_t structure.
 *
 * This function iterates through the `bc_data` hash table in the provided
 * bam_data_t structure and frees each entry's resources. The function then
 * destroys the hash table and sets the pointer to NULL.
 *
 * @param bam_data Pointer to the bam_data_t structure whose bc_data hash table 
 * resources will be freed.
 *
 * This function performs the following steps:
 * 1. Checks if `bam_data` is NULL; if so, it returns immediately.
 * 2. Checks if the `bc_data` hash table is NULL; if so, it returns immediately.
 * 3. Iterates over the hash table, freeing the keys and values.
 * 4. Destroys the hash table and sets the `bc_data` pointer to NULL.
 */
void bam_data_free_bcs(bam_data_t *bam_data);

/* @brief Free all memory in and including @p bam_dat.
 */
void bam_data_dstry(bam_data_t *bam_dat);

/* @brief Free the ATAC read pairs in bam_data.
 */
void bam_data_atac_free_pairs(bam_data_t *bam_data);

/**
 * @brief Fill barcodes in the bam_data_t structure.
 *
 * This function fills the `bcs` field in the provided bam_data_t structure with
 * barcodes. If a `str_map` is provided in the `bcs` argument, it will be copied. 
 * If the `bcs` argument is NULL, barcodes will be filled from the `bc_data` hash
 * table in the `bam_data_t` structure. The function handles initialization, 
 * copying, and error checks to ensure proper memory management.
 *
 * @param bam_data Pointer to the `bam_data_t` structure whose `bcs` field will be filled.
 * @param bcs Pointer to the `str_map` structure to be copied, or NULL to fill from `bc_data`.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bam_data_fill_bcs: argument is NULL" if the `bam_data` argument is NULL.
 * - "bam_data_fill_bcs: bam_data->bc_data is null" if the `bc_data` hash table in `bam_data` is NULL.
 * - "bam_data_fill_bcs: failed to add barcode" if adding a barcode to the string map fails.
 */
int bam_data_fill_bcs(bam_data_t *b, str_map *bcs);

/*******************************************************************************
 * BAM RNA
 ******************************************************************************/

/**
 * @brief Add an RNA read to the specified bam_data_t structure.
 *
 * This function adds an RNA read to the `bc_data` associated with 
 * the provided barcode (`bc`). If the barcode is not already present 
 * in the hash table, it is added along with an initialized `bc_data_t` structure.
 * The function ensures thread-safety by using mutex locks.
 *
 * @param bam_data Pointer to the `bam_data_t` structure where the RNA read will be added.
 * @param bc The barcode string associated with the RNA read.
 * @param r Pointer to the `rna_read1_t` structure representing the RNA read.
 * @param umih The UMI (Unique Molecular Identifier) associated with the RNA read.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bam_data_rna_add_read: argument is null" if `bam_data`, `bc`, or `r` is NULL.
 * - "bam_data_rna_add_read: bam_data->bc_data is null" if the `bc_data` hash table is NULL.
 * - "bam_data_rna_add_read: failed mutex lock" if locking the mutex fails.
 * - "bam_data_rna_add_read: failed mutex unlock" if unlocking the mutex fails.
 * - "bam_data_rna_add_read: %s" where %s represents strerror(errno) if memory allocation fails.
 * - "bam_data_rna_add_read: failed to add read to bcs" if adding the barcode to the hash table fails.
 * - "bam_data_rna_add_read: failed to add barcode to bcs" if initializing `bc_data` for the new barcode fails.
 * - "bam_data_rna_add_read: failed to add read to bam" if adding the RNA read to the `bc_data` fails.
 *
 * @note This function is thread-safe and uses mutex locks to guard against 
 * race conditions when accessing shared resources.
 */
int bam_data_rna_add_read(bam_data_t *bam_data, const char *bc, 
        const rna_read1_t *r, umishort umih);

/**
 * @brief De-duplicate RNA reads in the specified bam_data_t structure.
 *
 * This function traverses the `bc_data` hash table in the given `bam_data_t` 
 * structure and performs deduplication on the RNA data associated with each 
 * barcode. After deduplication, it frees the memory used by the `rna_dupls` 
 * tree to save space.
 *
 * @param bam_data Pointer to the `bam_data_t` structure containing RNA data.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bam_data_rna_dedup: bam_data is null" if the `bam_data` argument is NULL.
 * - "bam_data_rna_dedup: bam_data->bc_data is null" if the `bc_data` hash table is NULL.
 * - "bam_data_rna_dedup: barcode data bc_dat is null" if a barcode's associated `bc_data` is NULL.
 * - "bam_data_rna_dedup: failed to deduplicate reads" if deduplication of RNA reads fails for a barcode.
 */
int bam_data_rna_dedup(bam_data_t *bam_data);

/**
 * @brief Call RNA variants in the specified bam_data_t structure.
 *
 * This function traverses the `bc_data` hash table in the provided `bam_data_t`
 * structure and performs variant calling on the RNA data associated with each
 * barcode. The results are stored in the provided `g_var_t` structure.
 *
 * @param bam_data Pointer to the `bam_data_t` structure containing RNA data.
 * @param gv Pointer to the `g_var_t` structure where variant information will be stored.
 * @param cmap Pointer to the `str_map` structure used for chromosome name to index mapping.
 * @param min_qual Minimum base quality score to consider for calling a variant.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bam_data_rna_var_call: argument is null" if `bam_data`, `gv`, or `cmap` is NULL.
 * - "bam_data_rna_var_call: bam_data->bc_data is null" if the `bc_data` hash table is NULL.
 * - "bam_data_rna_var_call: barcode data bc_dat is null" if a barcode's associated `bc_data` is NULL.
 * - "bam_data_rna_var_call: failed to call variants" if variant calling on RNA data fails for a barcode.
 */
int bam_data_rna_var_call(bam_data_t *bam_data, g_var_t *gv, 
        str_map *cmap, uint8_t min_qual);

/*******************************************************************************
 * BAM ATAC
 ******************************************************************************/

/**
 * @brief Add an ATAC read to the specified bam_data_t structure.
 *
 * This function adds an ATAC read to the `bc_data` associated with the provided 
 * barcode (`bc`). If the barcode is not already present in the hash table, it is 
 * added along with an initialized `bc_data_t` structure. The function ensures 
 * thread-safety by using mutex locks.
 *
 * @param bam_data Pointer to the `bam_data_t` structure where the ATAC read will be added.
 * @param bc The barcode string associated with the ATAC read.
 * @param r Pointer to the `atac_read1_t` structure representing the ATAC read.
 * @param qname The query name associated with the ATAC read.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bam_data_atac_add_read: argument is null" if `bam_data`, `bc`, or `r` is NULL.
 * - "bam_data_atac_add_read: 'bc_data' is null" if the `bc_data` hash table is NULL.
 * - "bam_data_atac_add_read: failed bam mutex lock" if locking the `bam_lock` mutex fails.
 * - "bam_data_atac_add_read: %s" where %s represents strerror(errno) if memory allocation for the barcode copy fails.
 * - "bam_data_atac_add_read: failed to add barcode to bcs" if adding the barcode to the hash table fails.
 * - "bam_data_atac_add_read: failed to add barcode to bcs" if initializing `bc_data` for the new barcode fails.
 * - "bam_data_atac_add_read: failed bam mutex unlock" if unlocking the `bam_lock` mutex fails.
 * - "bam_data_atac_add_read: failed bc mutex lock" if locking the `bc_lock` mutex of the related `bc_data` fails.
 * - "bam_data_atac_add_read: failed bc mutex unlock" if unlocking the `bc_lock` mutex of the related `bc_data` fails.
 * - "bam_data_atac_add_read: failed to add read to bam" if adding the ATAC read to the `bc_data` fails.
 *
 * @note This function is thread-safe and uses mutex locks to guard against 
 * race conditions when accessing shared resources.
 */
int bam_data_atac_add_read(bam_data_t *bam_data, const char *bc, 
        const atac_read1_t *r, qshort qname);

/**
 * @brief Identify and mark duplicate ATAC reads in the specified bam_data_t structure.
 *
 * This function traverses the `bc_data` hash table in the provided `bam_data_t`
 * structure and identifies duplicates in the ATAC data associated with each
 * barcode. It calls `bc_data_atac_get_dups` for each barcode's associated `bc_data_t`.
 *
 * @param bam_data Pointer to the `bam_data_t` structure containing the ATAC data.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bam_data_atac_get_dups: bam_data is null" if the `bam_data` argument is NULL.
 * - "bam_data_atac_get_dups: bam_data->bc_data is null" if the `bc_data` hash table is NULL.
 * - "bam_data_atac_get_dups: barcode data bc_dat is null" if a barcode's associated `bc_data` is NULL.
 * - "bam_data_atac_get_dups: failed to get dups" if identifying duplicates for the ATAC data fails for a barcode.
 */
int bam_data_atac_get_dups(bam_data_t *bam_data);

/**
 * @brief De-duplicate ATAC reads in the specified bam_data_t structure.
 *
 * This function traverses the `bc_data` hash table in the provided `bam_data_t`
 * structure and performs deduplication on the ATAC data associated with each
 * barcode. It calls `bc_data_atac_dedup` for each barcode's associated `bc_data_t`.
 *
 * @param bam_data Pointer to the `bam_data_t` structure containing ATAC data.
 *
 * @return 0 on success, -1 on error with an appropriate error message.
 *
 * Possible error messages:
 * - "bam_data_atac_dedup: bam_data is null" if the `bam_data` argument is NULL.
 * - "bam_data_atac_dedup: bam_data->bc_data is null" if the `bc_data` hash table is NULL.
 * - "bam_data_atac_dedup: barcode data bc_dat is null" if a barcode's associated `bc_data` is NULL.
 * - "bam_data_atac_dedup: failed to deduplicate reads" if deduplication of ATAC reads fails for a barcode.
 */
int bam_data_atac_dedup(bam_data_t *bam_data);

/**
 * @brief Call variants in ATAC bam_data_t object.
 *
 * If no fragments or necessary parameters are present, this function does nothing.
 * It iterates over the barcode data contained in bam_data and attempts to call
 * variants using the bc_data_atac_var_call function for each barcode.
 *
 * @param bam_data Pointer to bam_data_t object.
 * @param gv Pointer to g_var_t object containing the variant information.
 * @param cmap Pointer to str_map object that maps chromosome names to indices.
 * @param min_qual Minimum base Phred quality in read to consider for calling a variant allele.
 *
 * @return -1 on error, or the number of variants called on success.
 * 
 * Errors:
 * - Returns -1 if any of bam_data, gv, or cmap are NULL.
 * - Returns -1 if bam_data->bc_data is NULL.
 * - Returns -1 if bc_data_atac_var_call fails for any barcode.
 * 
 * Appropriate error messages are provided for each failure case.
 */
int bam_data_atac_var_call(bam_data_t *bam_data, g_var_t *gv, 
        str_map *cmap, uint8_t min_qual);

/**
 * @brief Call peaks in ATAC bam_data_t object.
 *
 * If no fragments or necessary parameters are present, this function does nothing.
 * Iterates over the barcode data contained in bam_data and attempts to call
 * peaks using the bc_data_atac_peak_call function for each barcode.
 *
 * @param bam_data Pointer to bam_data_t object.
 * @param pks Pointer to iregs_t object containing the peaks.
 * @param cmap Pointer to str_map object that maps chromosome names to indices.
 *
 * @return -1 on error, or the number of peaks called on success.
 * 
 * Errors:
 * - Returns -1 if any of bam_data, pks, or cmap are NULL.
 * - Returns -1 if bam_data->bc_data is NULL.
 * - Returns -1 if bc_data_atac_peak_call fails for any barcode.
 */
int bam_data_atac_peak_call(bam_data_t *bam_data, iregs_t *pks, 
        str_map *cmap);

/*******************************************************************************
 * bc_stats
 ******************************************************************************/

/**
 * @brief Fill RNA statistics for a given barcode data.
 *
 * This function fills the statistics related to RNA molecules contained in bc_data_t. 
 * It iterates through RNA molecules, counting genes and variant calls, 
 * as well as calculating mitochondrial count and fraction of genes.
 *
 * @param bcc Pointer to bc_stats_t structure where statistics will be stored.
 * @param rna_hdr Pointer to sam_hdr_t structure containing RNA header information.
 * @param bc_dat Pointer to bc_data_t structure containing the barcode data.
 *
 * @return Returns 0 on success, or -1 on error.
 * 
 * Errors:
 * - Returns -1 if any of bcc, rna_hdr, or bc_dat are NULL.
 * - Returns -1 if iterating through RNA molecules fails.
 * - Returns -1 if adding to gene or variant ID hash tables fails.
 * - Returns -1 if chromosome name for RNA molecule is not found.
 * - Returns -1 if checking for mitochondrial chromosome fails.
 */
int bc_data_rna_fill_stats(bc_stats_t *bcc, const sam_hdr_t *rna_hdr, bc_data_t *bc_dat);

/**
 * @brief Fill ATAC statistics for a given barcode data.
 *
 * This function fills the statistics related to ATAC fragments contained in bc_data_t.
 * It iterates through ATAC fragments, counting peaks and variant calls,
 * as well as calculating mitochondrial count and fraction of reads in peaks.
 *
 * @param bcc Pointer to bc_stats_t structure where statistics will be stored.
 * @param atac_hdr Pointer to sam_hdr_t structure containing ATAC header information.
 * @param bc_dat Pointer to bc_data_t structure containing the barcode data.
 *
 * @return Returns 0 on success, or -1 on error.
 * 
 * Errors:
 * - Returns -1 if any of bcc, atac_hdr, or bc_dat are NULL.
 * - Returns -1 if iterating through ATAC fragments fails.
 * - Returns -1 if adding to peak or variant ID hash tables fails.
 * - Returns -1 if chromosome name for ATAC fragment is not found.
 * - Returns -1 if checking for mitochondrial chromosome fails.
 */
int bc_data_atac_fill_stats(bc_stats_t *bcc, const sam_hdr_t *atac_hdr, bc_data_t *bc_dat);

/**
 * @brief Fill statistics for both RNA and ATAC data in bam_data_t.
 *
 * This function iterates over barcode data in bam_data_t and fills the statistics 
 * for RNA and ATAC data. Allocates memory for bc_stats_t and updates the bc_stats 
 * field in each bc_data_t structure. It handles errors by freeing the allocated 
 * memory and returning an error code.
 *
 * @param bam_data Pointer to bam_data_t structure.
 * @param rna_hdr Pointer to sam_hdr_t structure containing RNA header information.
 * @param atac_hdr Pointer to sam_hdr_t structure containing ATAC header information.
 *
 * @return Returns 0 on success, or -1 on error.
 * 
 * Errors:
 * - Returns -1 if bam_data is NULL.
 * - Returns -1 if bam_data->bc_data is NULL.
 * - Returns -1 if a barcode data entry is NULL.
 * - Returns -1 if memory allocation for bc_stats_t fails.
 * - Returns -1 if filling RNA stats fails, with appropriate error message.
 * - Returns -1 if filling ATAC stats fails, with appropriate error message.
 */
int bam_data_fill_stats(bam_data_t *bam_data, 
        const sam_hdr_t *rna_hdr, const sam_hdr_t *atac_hdr);

typedef struct {
    char *bc;
    uint32_t count;
} bc_count_t;

static inline
int cmp_bc_count_rev(const void *bc1, const void *bc2){
    bc_count_t *bcc1 = (bc_count_t *)bc1;
    bc_count_t *bcc2 = (bc_count_t *)bc2;

    if (bcc1->count > bcc2->count) return(-1);
    if (bcc1->count < bcc2->count) return(1);
    return(0);
}

/**
 * @brief Retrieve the count of the top N-th barcode based on its statistics.
 *
 * This function retrieves the count of the N-th top barcode in terms of its statistics.
 * The function allocates memory for an array of `bc_count_t`, sorts the barcodes
 * based on their counts, and returns the count of the top N-th barcode.
 *
 * @param bam_data Pointer to bam_data_t structure containing barcode data.
 * @param top_n Index indicating which top N-th count to retrieve.
 * @param ret Pointer to an int to store return status (0 for success, -1 for error).
 *
 * @return Returns the count of the top N-th barcode on success, or 0 on error.
 *
 * Errors:
 * - Sets *ret to -1 and returns 0 if bam_data is NULL.
 * - Sets *ret to -1 and returns 0 if bam_data->bc_data is NULL.
 * - Sets *ret to -1 and returns 0 if memory allocation for bc_array fails.
 * - Sets *ret to -1 and returns 0 if top_n is out of range.
 * - Sets *ret to -1 and returns 0 if a barcode key or value is NULL.
 * - Sets *ret to -1 and returns 0 if the number of barcodes counted is inconsistent.
 *
 * Appropriate error messages are logged for each failure case.
 */
uint32_t bam_data_count_of_n(bam_data_t *bam_data, int top_n, int *ret);

#endif // BAM_DAT_H
