#include "hsilu.h"


void hsilu_round(hsilu_state_t *A, uint64_t rc) {
    // Keep the generated MRM words private; the permutation state is only
    // the five 64-bit words defined by the specification.
    uint64_t x[NROWS + NXOR] = {};
    for (int i = 0; i < NROWS; ++i) {
        x[i] = A->x[i];
    }

    // iota
    x[0] ^= rc;

    // chi
    uint64_t P = x[0];
    uint64_t Q = x[1];
    x[0] ^= (x[2] & ~x[1]);
    x[1] ^= (x[3] & ~x[2]);
    x[2] ^= (x[4] & ~x[3]);
    x[3] ^= (P & ~x[4]);
    x[4] ^= (Q & ~P);

    //theta

    for(int i = 5; i < 18; ++i) {
        x[i] = x[n1[i-5]] ^ rotate_left64(x[n2[i-5]], m[i-5]);
    }

    // rho
    for (int i = 0; i < NROWS; ++i) {
        A->x[i] = rotate_left64(x[NXOR + i], r[i]);
    }
}

void hsilu(hsilu_state_t *A) {
    const uint64_t hsilu_rc[] = {
            0x00000000000000F0, 0x00000000000000E1, 0x00000000000000D2,
            0x00000000000000C3, 0x00000000000000B4, 0x00000000000000A5,
            0x0000000000000096, 0x0000000000000087, 0x0000000000000078,
            0x0000000000000069, 0x000000000000005A, 0x000000000000004B
    };

    for(int i = BEGIN_ROUND; i < END_ROUND; ++i) {
        hsilu_round(A, hsilu_rc[i]);
    }
}
