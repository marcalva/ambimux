
#ifndef MATH_UTIL_H
#define MATH_UTIL_H

#include <math.h>
#include <assert.h>
#include <stdint.h>
#include "str_util.h"

#define ROUND_2_INT(f) ((int)(f >= 0.0 ? (f + 0.5) : (f - 0.5)))

float maxf(float *x, uint32_t n);

/* calculate log sum exp
 * Calculate the log of the sum of small or large numbers given in log scale.
 * @param x array of floats giving log numbers.
 * @param n length of x
 * @return the log(sum(exp(x))) value
 */
float logsumexpf(float *x, uint32_t n, int *ret);

static inline float logsum2expf(float x1, float x2, int *ret){
    float a[2];
    a[0] = x1;
    a[1] = x2;
    return logsumexpf(a, 2, ret);
}

double maxd(double *x, uint32_t n);
double maxdi(double *x, uint32_t n, uint32_t *i);
double logsumexpd(double *x, uint32_t n, int *ret);

double logsumexpd2(double *x, uint32_t n);

static inline double logsum2expd(double x1, double x2){
    double a[2];
    a[0] = x1;
    a[1] = x2;
    return logsumexpd2(a, 2);
}

// generate a pseudorandom positive number
uint32_t psr(uint32_t x, uint32_t mod);

// permute an integer array of n elements
// deterministic output using psr
void permute(int *array, uint32_t n);

// log( (x)/(1-x) )
#define logit(x) log( (x) / (1 - (x)) )

#define logis(x) (1 / (1 + exp(-(x))))

#endif // MATH_UTIL_H
