
#ifndef BINS_H
#define BINS_H

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

#define MAX_BIN (((1<<18)-1)/7)

// calculate bin given an alignment with [beg, end)
int reg2bin(int beg, int end);

/* Return the bins that the region overlaps with. 
 * Place the bin numbers in list.
 * Returns the number of elements in list.
 */
int reg2bins(int beg, int end, uint16_t list[MAX_BIN]);

#endif // BINS_H
