
#ifndef ARRAY_UTILS_H
#define ARRAY_UTILS_H

// [i,j] is [row,column]. I is number of rows, J is number of columns
// column-major indexing
#define CMI(i,j,I) ((I*j)+i) // column-major indexing
#define RMI(i,j,J) ((J*i)+j) // row-major indexing

/* write a matrix of doubles to file 
 *
 * Write a matrix to file. 
 * Only one of @p array1_cm array1_rm array2_cf array2_rf must be given.
 * array1_cm is a one-dimensional continuous array in column-major format, 
 * whereas array1_rm is in row-major format
 * array2_cf 
 *
 * @param fn file name
 * @param array1_cm one-dim continuous array in column-major format (CMI)
 * @param array1_rm one-dim continuous array in row-major format (RMI)
 * @param array2_cf two-dim array. First index is column, second is row
 * @param array2_rf two-dim array. First index is row, second is column
 * @param n1 number of elements in first index (rows)
 * @param n2 number of elements in second index (columns)
 * @param delim delimiter
 * @param nl newline character
 * @param decp number of decimal places to round to
 * @param rownames row names of length @p n1. Set to NULL to ignore
 * @param colnames column names of length @p n2. Set to NULL to ignore
 * @return 0 on success, -1 on error
 */
int write_matrix_double(char *fn, double *array1_cm, double *array1_rm, 
        double **array2_cf, double **array2_rf, 
        char **rownames, int nrow, char **colnames, int ncol, 
        char delim, char nl, int decp);

int write_matrix_float(char *fn, float *array1_cm, float *array1_rm, 
        float **array2_cf, float **array2_rf, 
        char **rownames, int nrow, char **colnames, int ncol, 
        char delim, char nl, int decp);

/* read a matrix of doubles to file
 *
 * The parameters array, rownames, and colnames should point to NULL.
 *
 * @param fn File name
 * @param array Pointer to 2D matrix array
 * @param rowcol Column containing row names (1-based). Set to 0 for no row names
 * @param header Set to 1 for header, 0 if no header
 * @param rownames Pointer to array of char arrays to to store row names
 * @param nrow Pointer to int. Stores the number of rows
 * @param colnames Pointer to array of char arrays to to store column names
 * @param ncol Pointer to int. Stores the number of columns
 * @param delims An array of char to split lines by
 * @param newline The newline character
 * @return 0 on success, -1 on error
 *
 * The memory in array, rownames, and colnames must be freed by the caller.
 */
int read_matrix_double(char *fn, double ***array, 
        int rowcol, int header,
        char ***rownames, int *nrow, char ***colnames, int *ncol, 
        char *delims, char newline);

int read_matrix_float(char *fn, float ***array, 
        int rowcol, int header,
        char ***rownames, int *nrow, char ***colnames, int *ncol, 
        char *delims, char newline);

#endif // ARRAY_UTILS_H

