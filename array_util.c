
#include "array_util.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "htslib/bgzf.h"
#include "htslib/kstring.h"
#include "str_util.h"

int write_matrix_double(char *fn, double *array1_cm, double *array1_rm, 
        double **array2_cf, double **array2_rf, 
        char **rownames, int nrow, char **colnames, int ncol, 
        char delim, char nl, int decp){

    int nnull = 0;
    if (array1_cm != NULL) nnull++;
    if (array1_rm != NULL) nnull++;
    if (array2_cf != NULL) nnull++;
    if (array2_rf != NULL) nnull++;

    if (nnull != 1)
        return err_msg(-1, 0, "write_matrix_double: only one of arrays must be given");

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "write_matrix_double: failed to create directory for %s", fn);

    int ret = 0;

    BGZF *fp;
    fp = bgzf_open(fn, "wg1");
    if (fp == 0){
        err_msg(-1, 0, "write_matrix_double: failed to open file %s", fn);
        return -1;
    }

    // write column names if set
    if (colnames != NULL){
        int c;
        for (c = 0; c < ncol; ++c){
            if (c || rownames != NULL) ret = bgzf_write(fp, &delim, 1);
            ret = bgzf_write(fp, colnames[c], strlen(colnames[c]));
            if (ret < 0) return err_msg(-1, 0, "write_matrix_double: failed to write to file %s", fn);
        }
        ret = bgzf_write(fp, &nl, 1);
        if (ret < 0) return err_msg(-1, 0, "write_matrix_double: failed to write to file %s", fn);
    }

    size_t buf_size = decp + 1000;
    char *pstr = (char *)malloc(buf_size * sizeof(char));
    if (pstr == NULL)
        return err_msg(-1, 0, "write_matrix_double: %s", strerror(errno));

    int i, j, pstr_len;
    for (i = 0; i < nrow; ++i){
        if (rownames != NULL) ret = bgzf_write(fp, rownames[i], strlen(rownames[i]));
        for (j = 0; j < ncol; ++j){
            // get par
            double par;
            if (array1_cm != NULL)      {par = array1_cm[CMI(i,j,nrow)];}
            else if (array1_rm != NULL) {par = array1_rm[RMI(i,j,ncol)];}
            else if (array2_cf != NULL) {par = array2_cf[j][i];}
            else if (array2_rf != NULL) {par = array2_rf[i][j];}
            else { return -1; }

            // convert par to string
            while ( (pstr_len = snprintf(pstr, buf_size, "%.*g", decp, par)) >= (int)buf_size){
                buf_size *= 2;
                void *tmp = realloc(pstr, buf_size * sizeof(char));
                if (tmp != NULL) pstr = tmp;
                else{
                    return err_msg(-1, 0, "write_matrix_double: %s", strerror(errno));
                }
            }
            if (pstr_len < 0)
                return err_msg(-1, 0, "write_matrix_double: %s", strerror(errno));

            // write to file
            if (j || rownames != NULL) ret = bgzf_write(fp, &delim, 1);
            ret = bgzf_write(fp, pstr, pstr_len);
        }
        ret = bgzf_write(fp, &nl, 1);
    }
    free(pstr);
    if (ret < 0){
        err_msg(-1, 0, "write_matrix_double: failed to write to file %s", fn);
        return -1;
    }
    bgzf_close(fp);
    return 0;
}

int read_matrix_double(char *fn, double ***array, 
        int rowcol, int header,
        char ***rownames, int *nrow, char ***colnames, int *ncol, 
        char *delims, char newline){

    BGZF *fp = bgzf_open(fn, "r");
    if (fp == 0){
        err_msg(-1, 0, "read_matrix_double: failed to open file %s", fn);
        return -1;
    }

    rowcol--; // 0-based

    *nrow = 0;
    *ncol = 0;

    // set up array
    int am = 1; // max length of array
    *array = realloc(*array, am * sizeof(double *));
    if (*array == NULL)
        return err_msg(-1, 0, "read_matrix_double: %s", strerror(errno));

    // row and column names
    str_map *rownames_tmp = init_str_map();
    str_map *colnames_tmp = init_str_map();
    if (rownames_tmp == NULL || colnames_tmp == NULL) return -1;

    // for tokenizing
    int slen = 0, sm = 0, found;
    char **tokens = NULL; // free

    kstring_t kstr = KS_INITIALIZE;
    int nc_set = 0; // is number of columns set?
    int ai = 0; // array index 
    int nc = 0, nl = 0, ret;
    while ((ret = bgzf_getline(fp, newline, &kstr)) >= 0){
        if (ret < -1){
            err_msg(-1, 0, "read_matrix_double: failed to read file %s", fn);
            return -1;
        }
        if (ret == 0) continue;

        // tokenize
        if (split_line(kstr.s, &tokens, delims, &slen, &sm) < 0) return -1;

        // header can contain rowname header or not
        if (nl == 0 && header){ // header line
            int i;
            for (i = 0; i < slen; ++i){
                if ( add2str_map(colnames_tmp, tokens[i], &found) < 0 ) return -1;
            }
            nc = slen;
        } else { // data line
            // get number of columns
            int line_err = 0;
            if (nc_set == 0){
                if (header){ // if nc was set by 0
                    if (rowcol >= 0){ // if row names, length can be nc or nc - 1
                        if ( (slen != nc) && (slen != (nc + 1))) line_err = 1;
                        else if (slen == nc){ // remove column field of row names
                            --nc;
                            char *todel = str_map_str(colnames_tmp, rowcol);
                            str_map_del(colnames_tmp, todel);
                        }
                    }
                }
                else {
                    if (rowcol >= 0) nc = slen - 1; // if row names
                    else nc = slen; // if no row names
                }
                nc_set = 1;
            }
            if (rowcol >= 0 && slen != (nc + 1)) line_err = 1;
            else if (rowcol < 0 && slen != nc) line_err = 1;
            if (line_err){
                err_msg(-1, 0, "read_matrix_double: number of columns (%i) in line %i "
                        "must be the same (%i) in file %s", slen, nl+1, nc, fn);
                return -1;
            }

            // resize array if needed
            while ((ai+1) > am){
                am *= 2;
                void *tmp = realloc(*array, am * sizeof(double *));
                if (tmp != NULL) *array = tmp;
                else {
                    return err_msg(-1, 0, "read_matrix_double: %s", strerror(errno));
                }
            }

            // check rowcol
            if (rowcol >= slen){
                err_msg(-1, 0, "read_matrix_double: rowcol (%i) must be from 0 to %i", 
                        rowcol+1, slen);
                return -1;
            }

            // store data
            int ci = 0;
            (*array)[ai] = malloc(nc * sizeof(double));
            if ((*array)[ai] == NULL)
                return err_msg(-1, 0, "read_matrix_double: %s", strerror(errno));

            int i;
            for (i = 0; i < slen; ++i){
                if (i == rowcol){
                    add2str_map(rownames_tmp, tokens[i], &found); // if row name
                    if (found)
                        return err_msg(-1, 0, "read_matrix_double: duplicate row name %s found in line %i",
                                tokens[i], nl+1);
                }
                else { // if value
                    // convert to double

                    char *p_end = NULL;

                    int save_errno = errno;
                    errno = 0;

                    double val = strtod(tokens[i], &p_end);

                    if ((int)val == 0 && errno > 0)
                        return err_msg(-1, 0, "read_matrix_double: %s", strerror(errno));

                    (*array)[ai][ci++] = val;

                    errno = save_errno;
                }
            }
            ++ai;
        }
        ++nl;
    }
    *array = realloc(*array, ai * sizeof(double *));
    *nrow = ai;
    *ncol = nc;

    // store row and column names
    if (header) *colnames = str_map_ca(colnames_tmp);
    if (colnames == NULL) return -1;

    if (rowcol >= 0) *rownames = str_map_ca(rownames_tmp);
    if (rownames == NULL) return -1;

    free(tokens);
    ks_free(&kstr);
    destroy_str_map(rownames_tmp);
    destroy_str_map(colnames_tmp);
    bgzf_close(fp);

    return 0;
}

int write_matrix_float(char *fn, float *array1_cm, float *array1_rm, 
        float **array2_cf, float **array2_rf, 
        char **rownames, int nrow, char **colnames, int ncol, 
        char delim, char nl, int decp){

    int nnull = 0;
    if (array1_cm != NULL) nnull++;
    if (array1_rm != NULL) nnull++;
    if (array2_cf != NULL) nnull++;
    if (array2_rf != NULL) nnull++;

    if (nnull != 1)
        return err_msg(-1, 0, "write_matrix_float: only one of arrays must be given");

    if (mkpath(fn, 0755) == -1)
        return err_msg(-1, 0, "write_matrix_float: failed to create directory for %s", fn);

    int ret = 0;

    BGZF *fp;
    fp = bgzf_open(fn, "wg1");
    if (fp == 0){
        err_msg(-1, 0, "write_matrix_float: failed to open file %s", fn);
        return -1;
    }

    // write column names if set
    if (colnames != NULL){
        int c;
        for (c = 0; c < ncol; ++c){
            if (c || rownames != NULL) ret = bgzf_write(fp, &delim, 1);
            ret = bgzf_write(fp, colnames[c], strlen(colnames[c]));
        }
        ret = bgzf_write(fp, &nl, 1);
        if (ret < 0){
            err_msg(-1, 0, "write_matrix_float: failed to write to file %s", fn);
            return -1;
        }
    }

    size_t buf_size = decp + 1000;
    char *pstr = (char *)malloc(buf_size * sizeof(char));
    if (pstr == NULL)
        return err_msg(-1, 0, "write_matrix_float: %s", strerror(errno));

    int i, j, pstr_len;
    for (i = 0; i < nrow; ++i){
        if (rownames != NULL) ret = bgzf_write(fp, rownames[i], strlen(rownames[i]));
        for (j = 0; j < ncol; ++j){
            // get par
            float par;
            if (array1_cm != NULL)      {par = array1_cm[CMI(i,j,nrow)];}
            else if (array1_rm != NULL) {par = array1_rm[RMI(i,j,ncol)];}
            else if (array2_cf != NULL) {par = array2_cf[j][i];}
            else if (array2_rf != NULL) {par = array2_rf[i][j];}
            else { return -1; }

            // convert par to string
            while ( (pstr_len = snprintf(pstr, buf_size, "%.*g", decp, par)) >= (int)buf_size){
                buf_size *= 2;
                void *tmp = realloc(pstr, buf_size * sizeof(char));
                if (tmp != NULL) pstr = tmp;
                else{
                    return err_msg(-1, 0, "write_matrix_float: %s", strerror(errno));
                }
            }

            // write to file
            if (j || rownames != NULL) ret = bgzf_write(fp, &delim, 1);
            ret = bgzf_write(fp, pstr, pstr_len);
        }
        ret = bgzf_write(fp, &nl, 1);
    }
    free(pstr);
    if (ret < 0){
        err_msg(-1, 0, "write_matrix_float: failed to write to file %s", fn);
        return -1;
    }
    bgzf_close(fp);
    return 0;
}

int read_matrix_float(char *fn, float ***array, 
        int rowcol, int header,
        char ***rownames, int *nrow, char ***colnames, int *ncol, 
        char *delims, char newline){

    BGZF *fp = bgzf_open(fn, "r");
    if (fp == 0){
        err_msg(-1, 0, "read_matrix_float: failed to open file %s", fn);
        return -1;
    }

    rowcol--; // 0-based

    *nrow = 0;
    *ncol = 0;

    // set up array
    int am = 1; // max length of array
    *array = realloc(*array, am * sizeof(float *));
    if (*array == NULL)
        return err_msg(-1, 0, "read_matrix_float: %s", strerror(errno));

    // row and column names
    str_map *rownames_tmp = init_str_map();
    str_map *colnames_tmp = init_str_map();
    if (rownames_tmp == NULL || colnames_tmp == NULL) return -1;

    // for tokenizing
    int slen = 0, sm = 0, found;
    char **tokens = NULL; // free

    kstring_t kstr = KS_INITIALIZE;
    int nc_set = 0; // is number of columns set?
    int ai = 0; // array index 
    int nc = 0, nl = 0, ret;
    while ((ret = bgzf_getline(fp, newline, &kstr)) >= 0){
        if (ret < -1){
            err_msg(-1, 0, "read_matrix_float: failed to read file %s", fn);
            return -1;
        }
        if (ret == 0) continue;

        // tokenize
        if (split_line(kstr.s, &tokens, delims, &slen, &sm) < 0) return -1;

        // header can contain rowname header or not
        if (nl == 0 && header){ // header line
            int i;
            for (i = 0; i < slen; ++i){
                if ( add2str_map(colnames_tmp, tokens[i], &found) < 0 ) return -1;
            }
            nc = slen;
        } else { // data line
            // get number of columns
            int line_err = 0;
            if (nc_set == 0){
                if (header){ // if nc was set by 0
                    if (rowcol >= 0){ // if row names, length can be nc or nc - 1
                        if ( (slen != nc) && (slen != (nc + 1))) line_err = 1;
                        else if (slen == nc){ // remove column field of row names
                            --nc;
                            char *todel = str_map_str(colnames_tmp, rowcol);
                            str_map_del(colnames_tmp, todel);
                        }
                    }
                }
                else {
                    if (rowcol >= 0) nc = slen - 1; // if row names
                    else nc = slen; // if no row names
                }
                nc_set = 1;
            }
            if (rowcol >= 0 && slen != (nc + 1)) line_err = 1;
            else if (rowcol < 0 && slen != nc) line_err = 1;
            if (line_err){
                err_msg(-1, 0, "read_matrix_float: number of columns (%i) in line %i "
                        "must be the same (%i) in file %s", slen, nl+1, nc, fn);
                return -1;
            }

            // resize array if needed
            while ((ai+1) > am){
                am *= 2;
                void *tmp = realloc(*array, am * sizeof(float *));
                if (tmp != NULL) *array = tmp;
                else {
                    return err_msg(-1, 0, "read_matrix_float: %s", strerror(errno));
                }
            }

            // check rowcol
            if (rowcol >= slen){
                err_msg(-1, 0, "read_matrix_float: rowcol (%i) must be from 0 to %i", 
                        rowcol+1, slen);
                return -1;
            }

            // store data
            int ci = 0;
            (*array)[ai] = malloc(nc * sizeof(float));
            if ((*array)[ai] == NULL)
                return err_msg(-1, 0, "read_matrix_float: %s", strerror(errno));

            int i;
            for (i = 0; i < slen; ++i){
                if (i == rowcol){
                    add2str_map(rownames_tmp, tokens[i], &found); // if row name
                    if (found)
                        return err_msg(-1, 0, "read_matrix_float: duplicate row name %s found in line %i",
                                tokens[i], nl+1);
                }
                else { // if value
                    // convert to float

                    char *p_end = NULL;

                    int save_errno = errno;
                    errno = 0;

                    float val = strtod(tokens[i], &p_end);

                    if ((int)val == 0 && errno > 0)
                        return err_msg(-1, 0, "read_matrix_float: %s", strerror(errno));

                    (*array)[ai][ci++] = val;

                    errno = save_errno;
                }
            }
            ++ai;
        }
        ++nl;
    }
    *array = realloc(*array, ai * sizeof(float *));
    *nrow = ai;
    *ncol = nc;

    // store row and column names
    if (header) *colnames = str_map_ca(colnames_tmp);
    if (colnames == NULL) return -1;

    if (rowcol >= 0) *rownames = str_map_ca(rownames_tmp);
    if (rownames == NULL) return -1;

    free(tokens);
    ks_free(&kstr);
    destroy_str_map(rownames_tmp);
    destroy_str_map(colnames_tmp);
    bgzf_close(fp);

    return 0;
}
