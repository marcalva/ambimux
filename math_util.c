
#include "math_util.h"
#include <math.h>
#include <errno.h>
#include <assert.h>
#include "str_util.h"

float maxf(float *x, uint32_t n){
    assert(x != NULL);
    assert(n > 0);
    uint32_t i;
    float max_val = x[0];
    for (i = 1; i < n; ++i){
        if (x[i] > max_val) max_val = x[i];
    }
    return max_val;
}

double maxd(double *x, uint32_t n){
    assert(x != NULL);
    assert(n > 0);
    uint32_t i;
    double max_val = x[0];
    for (i = 1; i < n; ++i){
        if (x[i] > max_val) max_val = x[i];
    }
    return max_val;
}

double maxdi(double *x, uint32_t n, uint32_t *i){
    assert(x != NULL);
    assert(n > 0);
    *i = 0;
    double max_val = x[*i];
    uint32_t j;
    for (j = 1; j < n; ++j){
        if (x[j] > max_val) {
            max_val = x[j];
            *i = j;
        }
    }
    return max_val;
}

float logsumexpf(float *x, uint32_t n, int *ret){
    *ret = 0;
    if (x == NULL){
        *ret = err_msg(-1, 0, "logsumexpf: x is null");
        return(0);
    }
    if (n == 0){
        *ret = err_msg(-1, 0, "logsumexpf: n = 0");
        return(0);
    }
    uint32_t i;
    // check input
    for (i = 0; i < n; ++i){
        if (isnan(x[i])){
            *ret = err_msg(-1, 0, "logsumexpf: x[%i]=%e", x[i], i);
            return(0);
        }
    }
    float sum_val = 0, max_val = maxf(x, n);
    for (i = 0; i < n; ++i)
        sum_val += expf(x[i] - max_val);
    return(logf(sum_val) + max_val);
}

double logsumexpd(double *x, uint32_t n, int *ret){
    *ret = 0;
    if (x == NULL){
        *ret = err_msg(-1, 0, "logsumexpd: x is null");
        return(0);
    }
    if (n == 0){
        *ret = err_msg(-1, 0, "logsumexpd: n = 0");
        return(0);
    }
    uint32_t i;
    // check input
    for (i = 0; i < n; ++i){
        if (isnan(x[i])){
            *ret = err_msg(-1, 0, "logsumexpd: x[%i]=%e", x[i], i);
            return(0);
        }
    }
    double sum_val = 0, max_val = maxd(x, n);

    // if max is nan, there is an error
    if (isnan(max_val)){
        *ret = err_msg(-1, 0, "logsumexpd: max val is nan");
        return(0);
    }

    // if all values are 0, return log(0).
    // if one value is infinity, return log(infinity).
    if (isinf(max_val)) return(max_val);

    for (i = 0; i < n; ++i)
        sum_val += exp(x[i] - max_val);
    assert(sum_val >= 0);
    assert(isfinite(sum_val));
    double r = log(sum_val) + max_val;
    return(r);
}

double logsumexpd2(double *x, uint32_t n){
    uint32_t i;
    double max_val = x[0];
    for (i = 1; i < n; ++i){
        if (x[i] > max_val) max_val = x[i];
    }

    // if all values are 0, return log(0).
    // if one value is infinity, return log(infinity).
    if (isinf(max_val)) return(max_val);

    double sum_val = 0;
    for (i = 0; i < n; ++i)
        sum_val += exp(x[i] - max_val);

    double r = log(sum_val) + max_val;
    return(r);
}

uint32_t psr(uint32_t x, uint32_t mod){
    x += 0x9e3779b9;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = ((x >> 16) ^ x) * 0x45d9f3b;
    x = (x >> 16) ^ x;
    x = x % mod;
    return(x);
}

void permute(int *array, uint32_t n){
    uint32_t i = 0, j;
    for (i = 0; i < n; ++i){
        int tmp = array[i];
        j = psr(i+1, (uint32_t)n);
        array[i] = array[j];
        array[j] = tmp;
    }
}

