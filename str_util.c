
#include "str_util.h"
#include "htslib/kstring.h"
#include "htslib/hts.h"
#include <sys/stat.h>
#include <string.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include <stdarg.h>
#include <time.h>

void get_time(char dt[], int len){
    time_t now;
    time(&now);
    struct tm *local = localtime(&now);
    strftime(dt, len, "%Y-%m-%d %H:%M:%S", local);
}

void log_msg(const char *fmt, ...){

    char dt[20];
    get_time(dt, 20);
    fprintf(stdout, "%s: ", dt);

    va_list args;
    va_start(args, fmt);
    vfprintf(stdout, fmt, args);
    va_end(args);

    fprintf(stdout, "\n");
    fflush(stdout);
}

int err_msg(int ret, int type, const char *fmt, ...){

    char *type_str;
    if (type == 0) type_str = "error";
    else type_str = "warning";

    fflush(stdout);

    fprintf(stderr, "%s: ", type_str);

    va_list args;
    va_start(args, fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);

    fprintf(stderr, "\n");

    return ret;
}


kstr_node *destroy_kstr_node(kstr_node *node){
    if (node == NULL) return NULL;
    kstr_node *n = node->next;
    ks_free(&(node->str));
    free(node);
    return n;
}

/*****************************
* str_map
******************************/


/* init str_map */
str_map *init_str_map(){
    str_map *sm = calloc(1, sizeof(str_map));

    if (sm == NULL){ // if no memory
        err_msg(-1, 0, "init_str_map: %s", strerror(errno));
        return NULL;
    }

    sm->strs = NULL;
    sm->ixs = kh_init(str_int);

    if (sm->ixs == NULL){
        err_msg(-1, 0, "init_str_map: %s", strerror(errno));
        return NULL;
    }

    sm->n = 0;
    sm->m = 0;
    return sm;
}

str_map *init_str_map_sized(int size){
    str_map *sm = calloc(1, sizeof(str_map));

    if (sm == NULL){ // if no memory
        err_msg(-1, 0, "init_str_map_sized: %s", strerror(errno));
        return NULL;
    }

    sm->strs = calloc(size, sizeof(char *));
    if (sm->strs == NULL){
        err_msg(-1, 0, "init_str_map_sized: %s", strerror(errno));
        return NULL;
    }

    sm->ixs = kh_init(str_int);
    kh_resize(str_int, sm->ixs, size);
    sm->n = 0;
    sm->m = size;
    return sm;
}

str_map *init_str_map_array(char **strs, int n){

    str_map *sm;

    if ( (sm = init_str_map_sized(n)) == NULL){
        err_msg(-1, 0, "init_str_map_array: %s", strerror(errno));
        return NULL;
    } else {
        int i, found;
        for (i = 0; i < n; ++i){
            if (add2str_map(sm, strs[i], &found) < 0) return NULL;
        }
    }

    return sm;
}

/* destroy str_map */
void destroy_str_map(str_map *sm){
    if (sm == NULL) return;
    int i;
    for (i = 0; i < sm->n; ++i){
        if (sm->strs[i]) free(sm->strs[i]);
    }
    if (sm->strs) free(sm->strs);
    kh_destroy(str_int, sm->ixs);
    free(sm);
}

/* test if key in hash */
int key_in_hash(khash_t(str_int) *hash, char *str){
    khint_t k;
    k = kh_get(str_int, hash, str);
    if ( k == kh_end(hash) ) return(0);
    return(1);
}

/* test if key in hash */
int hash_val(khash_t(str_int) *hash, char *str){
    khint_t k;
    k = kh_get(str_int, hash, str);
    if ( k == kh_end(hash) ) return(-1);
    return(kh_val(hash, k));
}

/* Add index to hash
 * If no_dup > 0 and str found in hash, return -1 */
int add_hash_ix(khash_t(str_int) *hash, const char *str, int no_dup){

    khint_t k;
    k = kh_get(str_int, hash, (char *)str);
    if ( no_dup && (k != kh_end(hash)) )
        return err_msg(-1, 0, "add_hash_ix: duplicate %s found", str);

    if ( k == kh_end(hash) ) {
        char *str_cpy = (char*)calloc(strlen(str)+1, sizeof(char));
        if (str_cpy == NULL)
            return err_msg(-1, 0, "add_hash_ix: %s", strerror(errno));
        strcpy(str_cpy, str);
        int ret;
        k = kh_put(str_int, hash, str_cpy, &ret);
        if (ret == -1)
            return err_msg(-1, 0, "add_hash_ix: failed to add %s to hash table", str);
        kh_val(hash, k) = (int)kh_size(hash) - 1;
    }

    return(kh_val(hash, k));
}

void free_hashix(khash_t(str_int) *hash){
    khint_t k;
    for (k = kh_begin(hash); k != kh_end(hash); k++){
        if ( !kh_exist(hash, k)) continue;
        char *s = kh_key(hash, k);
        free(s);
    }
    kh_destroy(str_int, hash);
}

int add2str_map(str_map *sm, const char *str, int *found){
    if (sm->m <= 8){
        sm->m = 8;
        sm->strs = (char **)realloc(sm->strs, sm->m * sizeof(char *));
        if (sm->strs == NULL)
            return err_msg(-1, 0, "add2str_map: %s", strerror(errno));
    }
    while (sm->n >= sm->m){
        if ( sm->m > (INT_MAX / 2))
            return err_msg(-1, 0, "add2str_map: int overflow. Cannot store %i elements. "
                                  "Current size %i", sm->m, sm->n);
        sm->m = 2 * (sm->m);
        sm->strs = (char **)realloc(sm->strs, sm->m * sizeof(char *));
        if (sm->strs == NULL)
            return err_msg(-1, 0, "add2str_map: %s", strerror(errno));
    }

    khint_t k;
    k = kh_get(str_int, sm->ixs, (char *)str);
    if (k == kh_end(sm->ixs)){
        *found = 0;

        char *str_cpy = strdup(str);

        if (str_cpy == NULL)
            return err_msg(-1, 0, "add2str_map: %s", strerror(errno));

        int ret;
        k = kh_put(str_int, sm->ixs, str_cpy, &ret);
        if (ret == -1)
            return err_msg(-1, 0, "add2str_map: failed to add %s to hash table", str);

        kh_val(sm->ixs, k) = sm->n;
        sm->strs[sm->n] = str_cpy;
        (sm->n)++;
    }
    else *found = 1;
    return(kh_val(sm->ixs, k));
}

int add_from_str_map(str_map *sm_dst, str_map *sm_src){
    int i, found;
    for (i = 0; i < sm_src->n; ++i){
        char *str = str_map_str(sm_src, i);
        if (add2str_map(sm_dst, str, &found) < 0) return -1;
    }
    return 0;
}

int str_map_ix(str_map *sm, char *str){
    khint_t k;
    k = kh_get(str_int, sm->ixs, str);
    if (k == kh_end(sm->ixs)) return -1;
    return(kh_val(sm->ixs, k));
}

char *str_map_str(str_map *sm, int ix){
    if (ix < 0 || ix >= sm->n) return NULL;
    return(sm->strs[ix]);
}

int str_map_del(str_map *sm, char *str){
    khint_t k;
    k = kh_get(str_int, sm->ixs, str);
    if (k == kh_end(sm->ixs)) return -1;
    char *key = kh_key(sm->ixs, k);
    int ix = kh_val(sm->ixs, k);
    if (ix == -1) return -1;
    free(key);
    kh_del(str_int, sm->ixs, k);
    int i;
    sm->n--;
    for (i = ix; i < sm->n; ++i){
        sm->strs[i] = sm->strs[i+1];
        k = kh_get(str_int, sm->ixs, sm->strs[i]);
        kh_val(sm->ixs, k) = i;
    }
    return 0;
}

str_map *str_map_copy(str_map *sm){
    str_map *cpy = init_str_map();
    if (cpy == NULL) return NULL;

    int i, n = sm->n, found;
    for (i = 0; i < n; ++i){
        char *key = str_map_str(sm, i);
        if (key == NULL) continue;
        if (add2str_map(cpy, key, &found) < 0) return NULL;
    }
    return cpy;
}

int write_str_map(str_map *sm, char *fn){
    BGZF *fp = bgzf_open(fn, "wg1");
    if (fp == 0)
        return err_msg(-1, 0, "write_str_map: could not open file %s\n", fn);
    
    int ret;
    int k, total = 0;
    for (k = 0; k < sm->n; ++k){
        ret = bgzf_write(fp, sm->strs[k], strlen(sm->strs[k]));
        ret = bgzf_write(fp, "\n", 1);
        if (ret < 0) err_msg(-1, 0, "write_str_map: failed write");
        else total += ret;
    }
    bgzf_close(fp);
    return total;
}

str_map * read_str_map(const char *fn){
    if (fn == NULL)
        return(NULL);

    BGZF *fp = bgzf_open(fn, "r");
    if (fp == NULL){
        err_msg(-1, 0, "read_str_map: failed to open file %s", fn);
        return NULL;
    }
    int ret;
    kstring_t line = KS_INITIALIZE;

    str_map *sm = init_str_map();
    if (sm == NULL) return NULL;

    int found = 0;
    while ((ret = bgzf_getline(fp, '\n', &line)) >= 0){
        if (ret < -1){
            err_msg(-1, 0, "read_str_map: failed to read from file %s", fn);
            return NULL;
        }
        else if (ret == -1) break;
        else if (ret == 0) continue;
        
        add2str_map(sm, line.s, &found);
        if (found == 1){
            err_msg(-1, 0, "read_str_map: duplicate index key %s found in %s", line.s, fn);
            return NULL;
        }
    }
    ks_free(&line);
    bgzf_close(fp);

    return sm;
}

char **str_map_ca(str_map *sm){
    int n = sm->n;
    char **ca = malloc(n * sizeof(char *));
    if (ca == NULL){
        err_msg(-1, 0, "str_map_ca: %s", strerror(errno));
        return NULL;
    }
    int i;
    for (i = 0; i < n; ++i){
        char *str = str_map_str(sm, i);
        ca[i] = strdup(str);
        if (ca[i] == NULL){
            err_msg(-1, 0, "str_map_ca: %s", strerror(errno));
            return NULL;
        }
    }
    return ca;
}

/*****************************
 *****************************/

char *aloc_sprintf(const char *format, ...){
    int r_size;

    va_list args;
    va_start(args, format);

    // copy
    va_list cpy;
    va_copy(cpy, args);
    
    // get size
    r_size = vsnprintf(NULL, 0, format, cpy);

    va_end(cpy);

    if (r_size < 0){
        err_msg(-1, 0, "aloc_sprintf: failed to copy string");
        return(NULL);
    }

    char *str = (char *)calloc(r_size + 1, sizeof(char));

    if (str == NULL){
        err_msg(-1, 0, "aloc_sprintf: %s", strerror(errno));
        return NULL;
    }

    r_size = vsprintf(str, format, args);

    va_end(args);

    return(str);
}

char *int2str(int x, size_t *len){
    size_t size = 30;
    char *str = malloc(size * sizeof(char));
    if (str == NULL){
        err_msg(-1, 0, "int2str: %s", strerror(errno));
        return NULL;
    }

    while ((*len = snprintf(str, size, "%i", x)) >= size){
        size *= 2;
        str = realloc(str, size * sizeof(char));
        if (str == NULL){
            err_msg(-1, 0, "int2str: %s", strerror(errno));
            return NULL;
        }
        *len = snprintf(str, size, "%i", x);
    }

    return str;
}

int int2strp(int x, char **strp, size_t *strp_size){
    int xlen;

    if ( (xlen = snprintf(*strp, *strp_size, "%i", x)) >= (int)*strp_size ){
        *strp_size = xlen + 1;
        *strp = realloc(*strp, *strp_size * sizeof(char));
        if (*strp == NULL)
            return err_msg(-1, 0, "int2strp: %s", strerror(errno));
        xlen = snprintf(*strp, *strp_size, "%i", x);
    }
    if (xlen < 0)
        return err_msg(-1, 0, "int2strp: %s", strerror(errno));

    return xlen;
}

int htspos2strp(hts_pos_t x, char **strp, size_t *strp_size){
    int xlen;
    if ( (xlen = snprintf(*strp, *strp_size, "%" PRIhts_pos, x)) >= (int)*strp_size ){
        *strp_size = xlen + 1;
        *strp = realloc(*strp, *strp_size * sizeof(char));
        if (*strp == NULL)
            return err_msg(-1, 0, "htspos2strp: %s", strerror(errno));
        xlen = snprintf(*strp, *strp_size, "%" PRIhts_pos, x);
    }
    if (xlen < 0)
        return err_msg(-1, 0, "htspos2strp: %s", strerror(errno));

    return(xlen);
}

char *double2str(double x, size_t *len){
    size_t size = 30;
    char *str = malloc(size * sizeof(char));
    if (str == NULL){
        err_msg(-1, 0, "double2str: %s", strerror(errno));
        return NULL;
    }

    while ( (*len = snprintf(str, size, "%f", x)) >= size){
        size *= 2;
        str = realloc(str, size * sizeof(char));
        if (str == NULL){
            err_msg(-1, 0, "double2str: %s", strerror(errno));
            return NULL;
        }
    }

    return str;
}

int double2str_in(double x, char **pstr, size_t *buf_size, int decp){
    int pstr_len = 0;
    while ( (pstr_len = snprintf(*pstr, *buf_size, "%.*g", decp, x)) >= (int)*buf_size){
        *buf_size *= 2;
        void *tmp = realloc(*pstr, *buf_size * sizeof(char));
        if (tmp != NULL) *pstr = tmp;
        else return err_msg(-1, 0, "write_res: %s", strerror(errno));
    }
    return pstr_len;
}

int float2str_ip(float x, char **pstr, size_t *buf_size, int decp){
    if (decp < 0)
        return err_msg(-1, 0, "float2str_ip: decp=%i must be nonnegative", decp);
    if (pstr == NULL)
        return err_msg(-1, 0, "float2str_ip: pstr is null");

    if (*pstr == NULL){
        *buf_size = 32;
        *pstr = malloc(*buf_size * sizeof(char));
        if (*pstr == NULL)
            return err_msg(-1, 0, "float2str_ip: %s", strerror(errno));
    }

    int pstr_len = 0;
    while ( (pstr_len = snprintf(*pstr, *buf_size, "%.*g", decp, x)) >= (int)*buf_size){
        *buf_size = pstr_len + 1;
        void *tmp = realloc(*pstr, *buf_size * sizeof(char));
        if (tmp != NULL) *pstr = tmp;
        else return err_msg(-1, 0, "float2str_ip: %s", strerror(errno));
    }
    return pstr_len;
}

// TODO finished writing this function
int str2float(char *s, float *x){
    if (s == NULL)
        return err_msg(-1, 0, "str2float: string s is null");

    char *p_end = NULL;

    int save_errno = errno;
    errno = 0;

    double val = strtod(s, &p_end);

    if ((int)val == 0 && errno > 0)
        return err_msg(-1, 0, "str2float: %s", strerror(errno));

    *x = val;
    errno = save_errno;
    return(0);
}

char *strcat2(const char *str1, const char *str2){
    size_t s1 = strlen(str1), s2 = strlen(str2);
    size_t s3 = s1 + s2 + 1;
    char *str3 = (char *)malloc(s3 * sizeof(char));
    if (str3 == NULL){
        err_msg(-1, 0, "strcat2: %s", strerror(errno));
        return NULL;
    }
    memcpy(str3, str1, s1);
    memcpy(str3 + s1, str2, s2 + 1);
    return str3;
}

// concatenate strings
// Returned array must be freed
char *cat_strs(const char **sa, int n, const char *sep){
    int len = 0;
    int i;
    for (i = 0; i < n; i++){
        len = len + strlen(sa[i]);
    }
    len = len + ( (n-1)*strlen(sep) ) + 1;

    char *s = (char*)calloc(len, sizeof(char));
    if (s == NULL){
        err_msg(-1, 0, "cat_strs: %s", strerror(errno));
        return NULL;
    }

    strcpy(s, sa[0]);
    for (i = 1; i < n; i++){
        strcat(s, sep);
        strcat(s, sa[i]);
    }
    
    return(s);
}

int get_line(FILE *fp, char **line, size_t *m){
    int n = 0;
    if (*line == NULL){
        *m = 1;
        *line = malloc((*m) * sizeof(char));
        if (*line == NULL)
            return err_msg(-1, 0, "get_line: %s", strerror(errno));
    }
    char c;
    while ( (c = fgetc(fp)) != EOF && c != '\n' && !ferror(fp) ){
        while (n >= (int)*m){
            *m = (*m)<<1;
            *line = (char *)realloc(*line, (*m) * sizeof(char));
            if (*line == NULL)
                return err_msg(-1, 0, "get_line: %s", strerror(errno));
        }
        (*line)[n++] = c;
    }
    while (n >= (int)*m){
        *m = (*m) + 1;
        *line = (char *)realloc(*line, (*m) * sizeof(char));
        if (*line == NULL)
            return err_msg(-1, 0, "get_line: %s", strerror(errno));
    }
    (*line)[n] = '\0';
    return n;
}

char **read_lines(const char *fn, int *len){
    *len = 0;
    int lines_alloc = 1;
    char **lines = (char **)calloc(lines_alloc, sizeof(char *));
    if (lines == NULL){
        err_msg(-1, 0, "read_lines: %s", strerror(errno));
        return NULL;
    }

    FILE *fp = NULL;
    fp = fopen(fn, "r");
    if (fp == NULL){
        err_msg(-1, 0, "read_lines: failed to open file %s: %s", fn, strerror(errno));
        return NULL;
    }

    char *line = NULL;
    size_t m = 0;
    int len2;
    while ( (len2 = get_line(fp, &line, &m)) > 0){
        while (*len >= lines_alloc){
            lines_alloc = 2 * lines_alloc;
            lines = (char **)realloc(lines, lines_alloc * sizeof(char *));
            if (lines == NULL){
                err_msg(-1, 0, "read_lines: %s", strerror(errno));
                return NULL;
            }
        }
        char *line_cpy = strdup(line);
        if (line_cpy == NULL){
            err_msg(-1, 0, "read_lines: %s", strerror(errno));
            return NULL;
        }
        lines[*len] = line_cpy;
        (*len)++;
    }
    lines_alloc = *len;
    lines = (char **)realloc(lines, lines_alloc * sizeof(char *));
    free(line);
    fclose(fp);
    return lines;
}

int split_line(char *s, char ***tokens, char *delim, int *len, int *m){
    int init_len = 8;
    *len = 0;
    if (*tokens == NULL){
        *m = init_len;
        *tokens = (char **)calloc(*m, sizeof(char *));
        if (*tokens == NULL)
            return err_msg(-1, 0, "split_line: %s", strerror(errno));
    }
    char *token = NULL;
    char *rest = NULL;

    token = strtok_r(s, delim, &rest);
    while (token != NULL){
        while (*len >= *m){
            *m *= 2;
            *tokens = (char **)realloc(*tokens, *m * sizeof(char *));
            if (*tokens == NULL)
                return err_msg(-1, 0, "split_line: %s", strerror(errno));
        }
        (*tokens)[(*len)++] = token;
        token = strtok_r(NULL, delim, &rest);
    }
    return 0;
}

int mkpath(char *fpath, mode_t mode){
    if (fpath == NULL) return -1;
    char *p = strchr(fpath + 1, '/');
    for ( ; p; p = strchr(p + 1, '/')){
        *p = '\0';
        if (mkdir(fpath, mode) == -1){
            if (errno != EEXIST){
                *p = '/';
                return -1;
            }
        }
        *p = '/';
    }
    return 0;
}

