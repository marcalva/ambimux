
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

/* Exon structure */
typedef struct Exon {
    int beg;
    int end;
    struct Exon *next; // for linked list. NULL if none.
} Exon;

/* Isoform structure */
typedef struct Isoform {
    int beg;
    int end;
    struct Exon *exons; // list of exons. NULL if empty.
    uint32_t exons_n; // number of exons in exons field.
} Isoform;

/* hash table of char* key Isoform* values */
KHASH_INIT(iso, char*, Isoform*, 1, kh_str_hash_func, kh_str_hash_equal);

/* Gene structure */
typedef struct Gene {
    char *id; // gene_id attribute in GTF
    char *name; // gene_name attribute in GTF
    char *type; // gene_type attribute in GTF
    int beg; // 0-based start (inclusive)
    int end; // 0-based start (exclusive)
    char strand;
    int chrm; // chromosome index
    int bin;
    khash_t(iso) *isoforms; // hash table of isoforms
    int isoforms_n; // number of isoforms in @field isoforms
    struct Gene *next; // NULL if empty.
} Gene;

/* Chromosome structure */
typedef struct Chrm {
    Gene *bins[MAX_BIN]; // array of pointers to Gene objects. NULL if empty.
    uint16_t genes_n[MAX_BIN]; // Number of genes in each element of bins.
} Chrm;

/* main struct to hold gene-transcript-exon info
 * Given a chromosome string, get the index in @field chrms with kh_get and kh_val.
 */
typedef struct {
    str_map *chrm_ix; // string to index in chrms array
    str_map *gene_ix; // gene to index. g->ix has same memory as key here
    khash_t(str_int) *gene_chrm; // gene to chrm ix. key is same memory as g->id
    khash_t(str_int) *gene_bin; // gene to bin. key is same memory as g->id
    Chrm **chrms; // array of pointers to Chrm objects.
    int chrms_m; // allocated max number of elements in chrms
} Annotation;

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
} gtf_line;


/****************
 * Functions
 ****************/

/****************************
 * Annotation structure
 ****************************/

/* Initialize exon object */
Exon *init_exon();

/* Destroy exon object.
 *
 * @return pointer to next member of @p e.
 */
Exon *destroy_exon(Exon *e);

/* Initialize isoform object */
Isoform *init_isoform();

/* Destroy isoform object */
void destroy_isoform(Isoform *iso);

/* Initialize Gene object */
Gene *init_gene();

/* destroy gene object
 *
 * @param g pointer to Gene object to destroy.
 * @return pointer to @p g->next member
 */
Gene *destroy_gene(Gene *g);

/* Initialize anno object */
Annotation *init_anno();

/* Initialize Chrm object
 * return NULL if memory couldn't be allocated
 */
Chrm *init_Chrm();

/* destroy annotation object a and free memory. */
void destroy_anno(Annotation *a);

/* Add chromosome to the annotation object.
 *
 * @param a pointer to annotation object.
 * @param c character array of chromosome name.
 * @return index of the chromosome in @p a->chrm_ix, or -1 on error.
 */
int add_chrm(Annotation *a, char *c);

/* Return gene object in a given gene name
 *
 * Returns the gene given the chrom index, the bin, and the gene ID.
 *
 * If not found, returns NULL. If found, returns the pointer to it.
 */
Gene *gene_from_name_chrm(Annotation *a, int chrm_ix, char *gene_id);
Gene *gene_from_name(Annotation *a, char *gene_id);

/****************************
 * GTF processing
 ****************************/

/* Initialize gtf line object */
gtf_line *init_gtf_line();

/* Clears the memory in gtf_line g and reset values */
void clear_gtf_line(gtf_line *g);

/* Clear memory and free gtf_line g object */
void destroy_gtf_line(gtf_line *g);

/* Add gene, isoform, or exon from GTF line to anno
 *
 * @param a pointer to annotation object
 * @param gl pointer to gtf_line object
 * @return 0 on success, -1 on failure.
 */
int gtf_gene2anno(Annotation *a, gtf_line *gl);
int gtf_iso2anno(Annotation *a, gtf_line *gl);
int gtf_exon2anno(Annotation *a, gtf_line *gl);
int gtf_line2anno(Annotation *a, gtf_line *gl);

/* Parse GTF line in string and place data into gtf_line
 *
 * The function makes successive calls to strtok_r on the 
 * string @p line. The corresponding GTF fields are populated .
 *
 * @param line string that contains the GTF line.
 * @param g pointer to gtf_line to populate data.
 * @return 0 on success, -1 on error.
 */
int parse_gtf_line(kstring_t *line, gtf_line *g);

/* Parse GTF line attributes.
 *
 * The attributes in @p g are stored in a single string. 
 * This function tokenizes the string and stores the key-value 
 * pairs as strings.
 *
 * @return 0 on success
 */
int parse_gtf_attr(gtf_line *g);

/* Get attribute value of key from gtf line.
 * Searches for the first occurence of key, and returns the 
 * corresponding value in the gtf line.
 *
 * @param g pointer to gtf_line object
 * @param key char array of key of GTF attribute
 * @return char array of the first occurence of the key's value.
 *  NULL if the attribute is not found.
 */
char *get_attr_val(gtf_line *g, char *key);

/* Test if key-value pair is present in gtf line
 *
 * Returns 1 if found, 0 if not found.
 */
int has_key_val(gtf_line *g, char *key, char *val);

/* Read GTF annotations into Annotation object.
 *
 * The GTF file must have attributes 'gene_id' and 'transcript_id'. The 
 * exons of a transcript must be listed after the transcript, and the 
 * transcripts must be listed after the genes. The basic parameter allows 
 * to filter the GTF to include only transcripts/exons that are tagged as 
 * basic.
 *
 * @param file char array containing file name of GTF.
 * @param basic specify whether to use only isoforms that are tagged as basic (basic=1).
 * @return pointer to Annotation object, or NULL if failure.
 *
 * Returned object must be freed with destroy_anno()
 */
Annotation *read_from_gtf(char *file, int basic);

/* Write out gene summary data
 *
 * Takes the genes in gene_ix and writes summary data in that order.
 * This writes out chr, start, end, strand, type, name and key.
 *
 * @return 0 on success, -1 on error.
 */
int write_gene_data(BGZF *fp, Annotation *anno, str_map *gene_ix);

/****************************
 * Region overlap
 ****************************/

/* Get features that overlap given region [beg, end).
 *
 * @param a Annotation object to retrieve features from.
 * @param ref Reference sequence name (chromosome) in character array.
 * @param beg 0-based start position of region.
 * @param end 0-based position the base pair after the end. Open end position.
 * @param stranded 0 if given region is unstranded, any other value if the region is stranded.
 * @param strand Character of strand, one of '+' or '-'. Ignored if @p stranded is 0.
 * @param genes array of pointers to Gene objects, for overlapping genes.
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
int feats_from_region_p(const Annotation *a, const char* ref, 
        int64_t beg, int64_t end, uint8_t stranded, char strand, 
        Gene ***genes, int *genes_len, double p);

/* Get features that overlap completely with @p set to 1 in the above company */
static inline int feats_from_region(const Annotation *a, const char* ref, hts_pos_t beg, 
        hts_pos_t end, uint8_t stranded, char strand, Gene ***genes, int *genes_len){
    return feats_from_region_p(a, ref, beg, end, stranded, strand, genes, genes_len, 1.0);
}

/****************************
 * Summary functions
 ****************************/

int n_feat(Annotation *a, int *n_gene, int *n_iso, int *n_exon);

int test_anno(char *file);

#endif //ANNO_H

