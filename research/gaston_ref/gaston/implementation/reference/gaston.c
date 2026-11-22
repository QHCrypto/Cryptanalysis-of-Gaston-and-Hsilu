//
// Created by Rusydi Makarim on 25/05/23.
//

#include "gaston.h"


void gaston_round(gaston_state_t *A, uint64_t rc) {
    uint64_t P = 0, Q = 0;

    // rho-east
    A->x[0] = rotate_left64(A->x[0], GASTON_e0);
    A->x[1] = rotate_left64(A->x[1], GASTON_e1);
    A->x[2] = rotate_left64(A->x[2], GASTON_e2);
    A->x[3] = rotate_left64(A->x[3], GASTON_e3);
    A->x[4] = rotate_left64(A->x[4], GASTON_e4);

    // theta
    P = A->x[0];
    P ^= A->x[1];
    P ^= A->x[2];
    P ^= A->x[3];
    P ^= A->x[4];
    P ^= rotate_left64(P, GASTON_r);
    Q = rotate_left64(A->x[0], GASTON_t0);
    Q ^= rotate_left64(A->x[1], GASTON_t1);
    Q ^= rotate_left64(A->x[2], GASTON_t2);
    Q ^= rotate_left64(A->x[3], GASTON_t3);
    Q ^= rotate_left64(A->x[4], GASTON_t4);
    Q ^= rotate_left64(Q, GASTON_s);
    P ^= Q;
    P = rotate_left64(P, GASTON_u);
    A->x[0] ^= P;
    A->x[1] ^= P;
    A->x[2] ^= P;
    A->x[3] ^= P;
    A->x[4] ^= P;

    // rho-west
    A->x[0] = rotate_left64(A->x[0], GASTON_w0);
    A->x[1] = rotate_left64(A->x[1], GASTON_w1);
    A->x[2] = rotate_left64(A->x[2], GASTON_w2);
    A->x[3] = rotate_left64(A->x[3], GASTON_w3);
    A->x[4] = rotate_left64(A->x[4], GASTON_w4);

    // iota
    A->x[0] ^= rc;

    // chi
    P = A->x[0];
    Q = A->x[1];
    A->x[0] ^= (A->x[2] & ~A->x[1]);
    A->x[1] ^= (A->x[3] & ~A->x[2]);
    A->x[2] ^= (A->x[4] & ~A->x[3]);
    A->x[3] ^= (P & ~A->x[4]);
    A->x[4] ^= (Q & ~P);
}

void gaston(gaston_state_t *A) {
    const uint64_t gaston_rc[] = {
            0x00000000000000F0, 0x00000000000000E1, 0x00000000000000D2,
            0x00000000000000C3, 0x00000000000000B4, 0x00000000000000A5,
            0x0000000000000096, 0x0000000000000087, 0x0000000000000078,
            0x0000000000000069, 0x000000000000005A, 0x000000000000004B
    };

    for(int i = GASTON_BEGIN_ROUND; i < GASTON_END_ROUND; ++i) {
        gaston_round(A, gaston_rc[i]);
    }
}
