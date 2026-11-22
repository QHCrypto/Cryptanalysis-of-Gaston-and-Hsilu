//
// Created by Rusydi Makarim on 25/05/23.
//

#ifndef GASTON_H
#define GASTON_H

#include <stdint.h>

#define GASTON_NROWS 5

#define GASTON_NROUNDS 12
#define GASTON_END_ROUND 12
#define GASTON_BEGIN_ROUND (GASTON_END_ROUND - GASTON_NROUNDS)

// parameters for theta
#define GASTON_t0 25
#define GASTON_t1 32
#define GASTON_t2 52
#define GASTON_t3 60
#define GASTON_t4 63
#define GASTON_r 1
#define GASTON_s 18
#define GASTON_u 23

// rho-east rotation offsets
#define GASTON_e0 0
#define GASTON_e1 60
#define GASTON_e2 22
#define GASTON_e3 27
#define GASTON_e4 4

// rho-west rotation offsets
#define GASTON_w0 0
#define GASTON_w1 56
#define GASTON_w2 31
#define GASTON_w3 46
#define GASTON_w4 43

typedef struct {
  uint64_t x[GASTON_NROWS];
} gaston_state_t;

static inline
uint64_t rotate_left64(uint64_t x, uint64_t r) {
    return (x << r) | (x >> (-r & 63));
}

void gaston(gaston_state_t *A);

#endif //GASTON_H
