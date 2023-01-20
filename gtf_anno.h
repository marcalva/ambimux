
#ifndef ANNO_H
#define ANNO_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/hts.h"
#include "kbtree.h"
#include "str_util.h"
#include "bins.h"

/* gtf file parsing */

enum spl {SPLICE, UNSPLICE, AMBIG, SNA};
#define N_SPL 3

/* columns */
enum {GTF_SEQNAME, 
    GTF_SOURCE, 
    GTF_FEAT, 
    GTF_START, 
    GTF_END, 
    GTF_SCORE, 
    GTF_STRAND, 
    GTF_ATTR};

#define FEAT "transcript"
#define GENE_NAME "gene_name"
#define GENE_ID "gene_id"
#define GENE_TYPE "gene_type"
#define TX_NAME "transcript_name"
#define TX_ID "transcript_id"
#define TX_TYPE "transcript_type"
#define EXON_ID "exon_id"
#define GENE "gene"
#define TX "transcript"
#define EXON "exon"

/* all @p beg and @p end are [beg,end) 
 * all coordinates are 0-based coordinates */

/* exon_t structure */
typedef struct exon_t {
    int beg;
    int end;
    struct exon_t *next; // for linked list. NULL if none.
} exon_t;

/* isoform_t structure */
typedef struct isoform_t {
    char *id;
    int beg;
    int end;
    struct exon_t *exons; // list of exons. NULL if empty.
    uint32_t exons_n; // number of exons in exons field.
} isoform_t;

#define iso_cmp(p, q) (strcmp(((p).id), ((q).id)))
KBTREE_INIT(kb_iso, isoform_t, iso_cmp);

/* gene_t structure */
typedef struct gene_t {
    char *id; // gene_id attribute in GTF
    char *name; // gene_name attribute in GTF
    char *type; // gene_type attribute in GTF
    int beg; // 0-based start (inclusive)
    int end; // 0-based start (exclusive)
    char strand;
    int chrm; // chromosome index
    int bin;
    kbtree_t(kb_iso) *bt_isoforms;
    int isoforms_n; // number of isoforms in @field isoforms
    struct gene_t *next; // NULL if empty.
} gene_t;

/* Chromosome structure */
typedef struct chr_genes_t {
    gene_t *bins[MAX_BIN]; // array of pointers to gene_t objects. NULL if empty.
    uint16_t genes_n[MAX_BIN]; // Number of genes in each element of bins.
} chr_genes_t;

/* main struct to hold gene-transcript-exon info
 * Given a chromosome string, get the index in @field chrms with kh_get and kh_val.
 */
typedef struct {
    str_map *chrm_ix; // string to index in chrms array
    str_map *gene_ix; // gene to index. g->ix has same memory as key here
    khash_t(str_int) *gene_chrm; // gene to chrm ix. key is same memory as g->id
    khash_t(str_int) *gene_bin; // gene to bin. key is same memory as g->id
    chr_genes_t **chrms; // array of pointers to chr_genes_t objects.
    int chrms_m; // allocated max number of elements in chrms
} gene_anno_t;

/* Structure to hold gtf information
 * The fields attr_tag and attr_val hold linked lists of the 
 * tag and value members in attribute field. These fields 
 * are populated from @field attribute after calling parse_gtf_attr.
 */
typedef struct {
    kstring_t chrname;
    kstring_t source;
    kstring_t feature;
    int start;
    int end;
    int score;
    char strand;
    int frame;
    kstring_t attribute; // holds the space delimited attribute field from a GTF line.
    kstr_node *attr_tag; // tag strings of each attribute
    kstr_node *attr_val; // value strings of each attribute. Match @field attr_tag.
    int n_attr;
} gtf1_t;


/****************
 * Functions
 ****************/

/****************************
 * gene_anno_t structure
 ****************************/

/* Initialize exon object */
exon_t *init_exon();

/* Destroy exon object.
 *
 * @return pointer to next member of @p e.
 */
exon_t *destroy_exon(exon_t *e);

/* Initialize isoform object */
int isoform_init(isoform_t *iso);
isoform_t *isoform_alloc();

/* Destroy isoform object */
void destroy_isoform(isoform_t *iso);

/* Initialize gene_t object */
gene_t *init_gene();

/* destroy gene object
 *
 * @param g pointer to gene_t object to destroy.
 * @return pointer to @p g->next member
 */
gene_t *destroy_gene(gene_t *g);

/* Initialize anno object */
gene_anno_t *init_anno();

/* Initialize chr_genes_t object
 * return NULL if memory couldn't be allocated
 */
chr_genes_t *init_chr_genes();

/* destroy annotation object a and free memory. */
void destroy_anno(gene_anno_t *a);

/* Add chromosome to the annotation object.
 *
 * @param a pointer to annotation object.
 * @param c character array of chromosome name.
 * @return index of the chromosome in @p a->chrm_ix, or -1 on error.
 */
int add_chrm(gene_anno_t *a, char *c);

/* Return gene object in a given gene name
 *
 * Returns the gene given the chrom index, the bin, and the gene ID.
 *
 * If not found, returns NULL. If found, returns the pointer to it.
 */
gene_t *gene_from_name_chrm(gene_anno_t *a, int chrm_ix, char *gene_id);
gene_t *gene_from_name(gene_anno_t *a, char *gene_id);

/****************************
 * GTF processing
 ****************************/

/* Initialize gtf line object */
gtf1_t *init_gtf1();

/* Clears the memory in gtf1_t g and reset values */
void clear_gtf1(gtf1_t *g);

/* Clear memory and free gtf1_t g object */
void destroy_gtf1(gtf1_t *g);

/* Add gene, isoform, or exon from GTF line to anno
 *
 * @param a pointer to annotation object
 * @param gl pointer to gtf1_t object
 * @return 0 on success, -1 on failure.
 */
int gtf_gene2anno(gene_anno_t *a, gtf1_t *gl);
int gtf_iso2anno(gene_anno_t *a, gtf1_t *gl);
int gtf_exon2anno(gene_anno_t *a, gtf1_t *gl);
int gtf1_to_anno(gene_anno_t *a, gtf1_t *gl);

/* Parse GTF line in string and place data into gtf1_t
 *
 * The function makes successive calls to strtok_r on the 
 * string @p line. The corresponding GTF fields are populated .
 *
 * @param line string that contains the GTF line.
 * @param g pointer to gtf1_t to populate data.
 * @return 0 on success, -1 on error.
 */
int parse_gtf1(kstring_t *line, gtf1_t *g);

/* Parse GTF line attributes.
 *
 * The attributes in @p g are stored in a single string. 
 * This function tokenizes the string and stores the key-value 
 * pairs as strings.
 *
 * @return 0 on success
 */
int parse_gtf_attr(gtf1_t *g);

/* Get attribute value of key from gtf line.
 * Searches for the first occurence of key, and returns the 
 * corresponding value in the gtf line.
 *
 * @param g pointer to gtf1_t object
 * @param key char array of key of GTF attribute
 * @return char array of the first occurence of the key's value.
 *  NULL if the attribute is not found.
 */
char *get_attr_val(gtf1_t *g, char *key);

/* Test if key-value pair is present in gtf line
 *
 * Returns 1 if found, 0 if not found.
 */
int has_key_val(gtf1_t *g, char *key, char *val);

/* Read GTF annotations into gene_anno_t object.
 *
 * The GTF file must have attributes 'gene_id' and 'transcript_id'. The 
 * exons of a transcript must be listed after the transcript, and the 
 * transcripts must be listed after the genes. The basic parameter allows 
 * to filter the GTF to include only transcripts/exons that are tagged as 
 * basic.
 *
 * @param file char array containing file name of GTF.
 * @param basic specify whether to use only isoforms that are tagged as basic (basic=1).
 * @return pointer to gene_anno_t object, or NULL if failure.
 *
 * Returned object must be freed with destroy_anno()
 */
gene_anno_t *read_from_gtf(char *file, int basic);

/* Write out gene summary data
 *
 * Takes the genes in gene_ix and writes summary data in that order.
 * This writes out chr, start, end, strand, type, name and key.
 *
 * @return 0 on success, -1 on error.
 */
int write_gene_data(BGZF *fp, gene_anno_t *anno, str_map *gene_ix);

/****************************
 * Region overlap
 ****************************/

/* Get features that overlap given region [beg, end).
 *
 * @param a gene_anno_t object to retrieve features from.
 * @param ref Reference sequence name (chromosome) in character array.
 * @param beg 0-based start position of region.
 * @param end 0-based position the base pair after the end. Open end position.
 * @param stranded 0 if given region is unstranded, any other value if the region is stranded.
 * @param strand Character of strand, one of '+' or '-'. Ignored if @p stranded is 0.
 * @param genes array of pointers to gene_t objects, for overlapping genes.
 * @param genes_len Pointer to integer that contains the length of the genes array.
 * @param p Minimum fraction of the region that overlaps.
 *
 * @return number of overlapping features returned, -1 on error.
 *
 * The function will reallocate the @p genes array to fit the number of overlapping genes. 
 * If @p ref chromosome is not found in @p a, then 0 is returned and @p genes is unchanged.
 * If no genes overlap, then 0 is returned and @p genes is unchanged.
 * The object @p genes must be allocated before the function is called, and freed after the call.
 */
int feats_from_region_p(const gene_anno_t *a, const char* ref, 
        int32_t beg, int32_t end, uint8_t stranded, char strand, 
        gene_t ***genes, int *genes_len, double p);

/* Get features that overlap completely with @p set to 1 in the above company */
static inline int feats_from_region(const gene_anno_t *a, const char* ref, int32_t beg, 
        int32_t end, uint8_t stranded, char strand, gene_t ***genes, int *genes_len){
    return feats_from_region_p(a, ref, beg, end, stranded, strand, genes, genes_len, 1.0);
}

/****************************
 * Summary functions
 ****************************/

int n_feat(gene_anno_t *a, int *n_gene, int *n_iso, int *n_exon);

int test_anno(char *file);

#endif //ANNO_H

