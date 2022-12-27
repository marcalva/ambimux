
#include "bins.h"
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <inttypes.h>

// Binning scheme with 6 levels
/* Bin 0:         512 Mb
 * Bin 1-8        64 Mb
 * Bin 9-72       8 Mb
 * Bin 73-584     1 Mb
 * Bin 585-4680   128 Kb
 * Bin 4681-37448 16 Kb
 *
 * 32,449 total bins
 */

// calculate bin given an alignment with [beg, end)
int reg2bin(int beg, int end)
{
    --end;
    if (beg>>14 == end>>14) return ((1<<15))/7 + (beg>>14);
    if (beg>>17 == end>>17) return ((1<<12))/7 + (beg>>17);
    if (beg>>20 == end>>20) return ((1<<9))/7 + (beg>>20);
    if (beg>>23 == end>>23) return ((1<<6))/7 + (beg>>23);
    if (beg>>26 == end>>26) return ((1<<3))/7 + (beg>>26);
    return 0;
}

int reg2bins(int beg, int end, uint16_t list[MAX_BIN]){
    int i = 0,k;
    --end;
    list[i++] = 0;
    for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) list[i++] = k;
    for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) list[i++] = k;
    for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) list[i++] = k;
    for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) list[i++] = k;
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) list[i++] = k;
    return i;
}

