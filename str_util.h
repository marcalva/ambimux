
#ifndef STR_UTIL_H
#define STR_UTIL_H

#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include "htslib/kstring.h"
#include "htslib/khash.h"
#include "htslib/bgzf.h"
#include "htslib/hts.h"

/* log message
 * 
 * prints a message like printf but with the date and time 
 * prepended.
 */
void log_msg(const char *fmt, ...);

// error message
// @param ret integer for function to return
// @param type 0 for error, 1 for warning
// @param fmt format string
int err_msg(int ret, int type, const char *fmt, ...);

int write_fail();

static inline void print_status_plp(const char* k, int i, const char *chr, int pos){
    log_msg("processed %" PRIu64 " %s through %s:%i", i, k, chr, pos);
}

static inline uint64_t str_int64_hash(const char *s){
    if (s == NULL)
        return(0);

    uint64_t h = 5381;
	h = (uint64_t)*s;
	if (h) for (; *s; ++s) h = ( (h << 5) + h) + ((uint64_t)*s);
	return h;
}

static inline uint64_t str_int64_hash2(const char *s){
    if (s == NULL)
        return(0);

    uint64_t h = 0;

	h = (uint64_t)*s;
	if (h) for (; *s; ++s)
        h = ((uint64_t)*s) + (h << 6) + (h << 16) - h;

    return h;
}

/*****************************
 * linked list of kstring 
 *****************************/

typedef struct kstr_node {
    kstring_t str;
    struct kstr_node *next;
} kstr_node;

/* Initiailize a kstr_node object.
 * Object must be freed */
static inline kstr_node *init_kstr_node(){
    kstr_node *n = (kstr_node*)calloc(1, sizeof(kstr_node));

    if (n == NULL){
        err_msg(-1, 0, "init_kstr_node: %s", strerror(errno));
        return NULL;
    }

    ks_initialize(&(n->str));
    n->next = NULL;
    return n;
}

/* destroy kstr node 
 *
 * Frees the kstring_t @p node->str and calls free on the @p node.
 *
 * @return @p node->next
 */
kstr_node *destroy_kstr_node(kstr_node *node);


/*****************************
 * str_map
 *
 * Provides an ordered map with hash table. 
 * The value of a key string in the hash table 
 * gives the index of that key string.
 *****************************/


KHASH_INIT(str_int, char*, int, 1, kh_str_hash_func, kh_str_hash_equal);

/* str_map
 * index to string, or string to index
 * The char arrays in @p strs and @p ixs point to the same memory location.
 * Only one free is necessary.
 *
 * @field strs array (length len) of char arrays. Maps index to string.
 * @field ixs str_int hash table. Maps string to index.
 * @field n number of elements with data.
 * @field len number of valid elements in array (total added). Can be greater than n if elements are deleted.
 * @field m allocated length of strs.
 */
typedef struct {
    char **strs;
    khash_t(str_int) *ixs;
    int n;
    int m;
} str_map;

/* init str_map */
str_map *init_str_map();
str_map *init_str_map_sized(int size);

/* initialize from strings in @p strs of length @p n
 * This copies each string.
 * Returns a pointer to a newly allocated str_map.
 * Returns NULL on error.
 */
str_map *init_str_map_array(char **strs, int n);

/* destroy str_map */
void destroy_str_map(str_map *sm);

/* hash function helpers */

/* Return 1 if str key is in hash, 0 if not found */
int key_in_hash(khash_t(str_int) *hash, char *str);

/* Get the int value of str key in hash 
 * return -1 if not found. Return index if found*/
int hash_val(khash_t(str_int) *hash, char *str);

/* Add str to hash, placing the next integer index as value
 * Creates a copy of str, so must be freed.
 * If no_dup > 0 and str found in hash, return -1 
 * If no_dup = 0 and str found in hash, return its index value.
 * If not found, return its newly placed integer index value*/
int add_hash_ix(khash_t(str_int) *hash, const char *str, int no_dup);

/* Free memory allocated from string copies in add_hash_ix
 * Frees khash object */
void free_hashix(khash_t(str_int) *hash);

/* */

/* add string to str_map.
 * The string memory is copied.
 * Return index in sm->strs, or -1 on error. 
 * If found, update @p found to 1, otherwise return it at 0.
 */
int add2str_map(str_map *sm, const char *str, int *found);

/* add strings from another str map */
int add_from_str_map(str_map *sm_dst, str_map *sm_src);

/* get index from string
 * return index of character string if found, otherwise return -1
 */
int str_map_ix(str_map *sm, char *str);

/* get string from index 
 * return pointer to char if ix is valid, otherwise return NULL
 * */
char *str_map_str(str_map *sm, int ix);

/* delete string from str_map
 * this is an expensive operation.
 * removes the string and decrements all elements after.
 * return -1 if not found or error, 0 on success
 */
int str_map_del(str_map *sm, char *str);

/* copy str_map
 */
str_map *str_map_copy(str_map *sm);

/* resize array to sm->n to free memory
 */
static inline void str_map_resize(str_map *sm){
    sm->m = sm->n;
    sm->strs = (char **)realloc(sm->strs, sm->m * sizeof(char *));
}

/* write str_map to file 
 * @return total number of bytes written, -1 on error
 * */
int write_str_map(str_map *sm, char *fn);

/* read str_map from file */
str_map * read_str_map(const char *fn);

/* return a copy of the character arrays in str_map */
char **str_map_ca(str_map *sm);



/*****************************
 *****************************/

/* Return a string after applying sprintf to the arguments.
 *
 * Dynamically allocates a character array to fit the string, 
 * calls sprintf, and returns the string.
 *
 * @param format The format string (see printf).
 * @param ... Variant argumens (see printf).
 *
 * @return NULL on error, the string on success.
 * 
 * @note The string must be freed by the caller.
 */
char *aloc_sprintf(const char *format, ...);

/* Format a number to a string
 *
 * Convert int or double to char array string.
 * The returned array must be freed by the caller.
 *
 * @param x int to format to char array
 * @param len length of returned char array, not including NULL.
 * @return char array with formatted string
 */
char *int2str(int x, size_t *len);
char *double2str(double x, size_t *len);
int double2str_in(double x, char **pstr, size_t *buf_size, int decp);

int float2str_ip(float x, char **pstr, size_t *buf_size, int decp);

/* place string representation of an int into a char array
 *
 * convert the int to a string, and place in an array 
 * pointed to by the pointer pointed to by strp. Places 
 * the length of the string in @p len
 *
 * @param x int to format to char array.
 * @param strp pointer to char array.
 * @param strp_size pointer to size_t. Stores the allocated size of the 
 * char array pointed to by strp.
 * @return the string length, including the null byte.
 *
 * On error, returns -1, and the pointer pointed to by strp may be NULL 
 * if a memory error is encountered.
 */
int int2strp(int x, char **strp, size_t *strp_size);
int htspos2strp(hts_pos_t x, char **strp, size_t *strp_size);

/* Concatenate two strings by allocating memory.
 * Returned char array must be freed.
 */
char *strcat2(const char *str1, const char *str2);

void get_time(char dt[], int len);

char *cat_strs(const char **sa, int n, const char *sep);

static inline void destroy_strint_array(char **array, int len){
    if (array != NULL){
        int i;
        for (i = 0; i < len; ++i) free(array[i]);
        free(array);
    }
}

/* Read lines into array of char arrays.
 *
 * @param fn File name
 * @return Array of char arrays, or NULL if failure.
 */
char **read_lines(const char *fn, int *len);

/* split line by delimiter
 * update @p fields
 * The *p tokens parameter is allocated and updated after the call, 
 * and must be freed after. Each element of *tokens points to the 
 * char array in @p s, so these values aren't copied.
 *
 * @param s char array to split into tokens
 * @param tokens pointer to array of char arrays that contains the tokens.
 * @param len pointer to int that contains the number of tokens
 * @param m pointer to int that gives the allocated size of @p *tokens
 */
int split_line(char *s, char ***tokens, char *delim, int *len, int *m);

/* Create directories if not present from a file path. 
 *
 * @param fpath string of file path.
 * @param mode ownership properties of directories (0755).
 * @return 0 on success, -1 on error.
 * */
int mkpath(char *fpath, mode_t mode);

#endif // STR_UTIL_H
