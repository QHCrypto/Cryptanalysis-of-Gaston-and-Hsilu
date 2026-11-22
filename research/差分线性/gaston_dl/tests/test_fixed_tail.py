#!/usr/bin/env python3
"""Independent value enumeration for the conditional DDT/Walsh prefix evaluator."""
import itertools
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "tools"))
from fixed_tail_search import (ExactPrefix, MASK, S, ZERO, column, dot, embed,
                               linear, state_chi, sxor, transpose)


def enumerate_active_values(delta, gamma, rc):
    positions = [j for j in range(64) if column(delta, j)]
    outside = MASK ^ sum(1 << j for j in positions)
    total = 0
    for xs in itertools.product(range(32), repeat=len(positions)):
        z = beta = ZERO
        for j, x in zip(positions, xs):
            z = sxor(z, embed(j, S[x]))
            beta = sxor(beta, embed(j, S[x] ^ S[x ^ column(delta, j)]))
        d = linear(beta)
        # Boolean derivative of chi: h_i = gamma_(i-1)*d_(i+1)
        #                                   xor gamma_(i-2)*d_(i-1).
        h = tuple((gamma[(i - 1) % 5] & d[(i + 1) % 5])
                  ^ (gamma[(i - 2) % 5] & d[(i - 1) % 5]) for i in range(5))
        H = transpose(h)
        if any(v & outside for v in H):
            continue
        y = list(linear(z))
        y[0] ^= rc
        # Directly evaluate the derivative at y rather than using its constant.
        phase = dot(gamma, sxor(state_chi(y), state_chi(sxor(y, d))))
        total += (-1) ** phase
    return total, 32 ** len(positions)


def main():
    rng = random.Random(20261015)
    gammas = [embed(7, 19), (0x0200200200222020, 0x2002002002220200, 0, 0, 0),
              tuple(rng.getrandbits(64) for _ in range(5))]
    deltas = [embed(11, 8), sxor(embed(11, 8), embed(15, 8)),
              sxor(embed(3, 17), embed(40, 5))]
    checked = 0
    for gamma in gammas:
        evaluator = ExactPrefix(gamma)
        for delta, rc in itertools.product(deltas, (0, 0x5a)):
            result = evaluator.evaluate(delta, rc)
            assert (result["numerator"], result["denominator"]) == enumerate_active_values(delta, gamma, rc)
            checked += 1
    print(f"PASS {checked} single/two-column comparisons with independent active-value enumeration")


if __name__ == "__main__":
    main()
