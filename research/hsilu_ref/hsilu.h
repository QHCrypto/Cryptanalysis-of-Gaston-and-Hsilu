#ifndef HSILU_H
#define HSILU_H

#include <stdint.h>

#define NROWS 5
#define NXOR 13

#define NROUNDS 12
#define END_ROUND 12
#define BEGIN_ROUND (END_ROUND - NROUNDS)

// Hsilu parameters from Table 9 of "Multiple Rows Mixers and Hsilu".
// The entries are indexed by the generated words x[5] ... x[17].
const int n1[NXOR]={0, 3, 1, 0, 5, 5, 7,11,12,12,13,14,15};
const int n2[NXOR]={4, 4, 2, 1, 2, 3, 6, 8, 9, 6,10,11, 9};
const int m[NXOR] ={47,51,54,30, 4,55, 0,14,59, 4,10,45,53};

// The row rotations are applied to x[13] ... x[17], respectively.
const int r[NROWS]={0,57,13,31,20};

typedef struct {
  // Hsilu's external state is five 64-bit words (320 bits).
  uint64_t x[NROWS];
} hsilu_state_t;

static inline
uint64_t rotate_left64(uint64_t x, unsigned int r) {
    return (x << r) | (x >> (-r & 63));
}

void hsilu(hsilu_state_t *A);

#endif //HSILU_H
