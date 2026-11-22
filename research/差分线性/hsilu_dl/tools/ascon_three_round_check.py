#!/usr/bin/env python3
"""A structural certificate for a three-round Ascon DL correlation of +1.

Coordinates are (word, LSB-based bit). All three diffusion layers are retained.
The cone is an overapproximation: absent coordinates certify independence;
present coordinates are not claimed to be essential variables.
Reference: https://ascon.isec.tugraz.at/specification.html
S box: https://ascon.isec.tugraz.at/images/sbox_instructions.c
"""

import argparse
import json
from pathlib import Path

import numpy as np

MASK = (1 << 64) - 1
ROT = ((19, 28), (61, 39), (1, 6), (10, 17), (7, 41))
RC = (0x69, 0x5A, 0x4B)
TABLE = (4, 11, 31, 20, 26, 21, 9, 2, 27, 5, 8, 18, 29, 3, 6, 28,
         30, 19, 7, 14, 0, 13, 17, 24, 16, 12, 1, 25, 22, 10, 15, 23)


def sbox(state):
    x = list(state)
    x[0] ^= x[4]
    x[4] ^= x[3]
    x[2] ^= x[1]
    y = [x[i] ^ ((~x[(i + 1) % 5]) & x[(i + 2) % 5]) for i in range(5)]
    y[1] ^= y[0]
    y[0] ^= y[4]
    y[3] ^= y[2]
    y[2] = ~y[2]
    return [v & MASK for v in y]


def ror(x, n):
    return ((x >> n) | (x << (64 - n))) & MASK


def linear(state):
    return [x ^ ror(x, a) ^ ror(x, b) for x, (a, b) in zip(state, ROT)]


def permutation(state):
    x = [v.copy() if hasattr(v, 'copy') else v for v in state]
    for rc in RC:
        x[2] ^= rc
        x = linear(sbox(x))
    return x


def table_permutation(state):
    """Independent column-table implementation for a small implementation check."""
    x = list(state)
    for rc in RC:
        x[2] ^= rc
        y = [0] * 5
        for j in range(64):
            column = sum(((x[i] >> j) & 1) << (4 - i) for i in range(5))
            for i in range(5):
                y[i] |= ((TABLE[column] >> (4 - i)) & 1) << j
        x = linear(y)
    return x


def bit_positions(state):
    return [[i, j] for i, word in enumerate(state) for j in range(64) if word >> j & 1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples-log2', type=int, default=22)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    for v in range(32):
        out = sbox([(v >> (4 - i)) & 1 for i in range(5)])
        assert sum((out[i] & 1) << (4 - i) for i in range(5)) == TABLE[v]
    deps = [[i for i in range(5) if any(
        ((TABLE[v] ^ TABLE[v ^ (1 << (4 - i))]) >> (4 - o)) & 1
        for v in range(32))] for o in range(5)]
    cones = []
    for lane in range(5):
        x = [int(i == lane) for i in range(5)]
        trace = []
        for r in range(3):
            y = []
            for indices in deps:
                word = 0
                for i in indices:
                    word |= x[i]
                y.append(word)
            x = y
            trace.append({'layer': f'S{r + 1}', 'potentially_dependent_bits':
                          sum(v.bit_count() for v in x)})
            if r == 2 and lane == 0:
                core_cone = list(x)
            x = [v | ror(v, a) | ror(v, b) for v, (a, b) in zip(x, ROT)]
            trace.append({'layer': f'L{r + 1}', 'potentially_dependent_bits':
                          sum(v.bit_count() for v in x)})
        cones.append({'input_delta': [lane, 0], 'trace': trace,
                      'certified_independent_output_bits': bit_positions([MASK ^ v for v in x])})
    assert [0, 14] in cones[0]['certified_independent_output_bits']
    assert all(not (core_cone[0] >> j & 1) for j in (14, 33, 42))

    seed = 2026101901
    rng = np.random.default_rng(seed)
    for _ in range(128):
        state = [int(v) for v in rng.integers(0, 1 << 64, size=5, dtype=np.uint64)]
        assert permutation(state) == table_permutation(state)
    n = 1 << args.samples_log2
    ones = 0
    for offset in range(0, n, 1 << 16):
        count = min(1 << 16, n - offset)
        state = list(rng.integers(0, 1 << 64, size=(5, count), dtype=np.uint64))
        other = [v.copy() for v in state]
        other[0] ^= np.uint64(1)
        a, b = permutation(state), permutation(other)
        ones += int(np.count_nonzero(((a[0] ^ b[0]) >> 14) & 1))
    assert ones == 0
    result = {
        'definition': 'Three complete Ascon rounds: RC -> S -> L, repeated three times',
        'coordinates': '(word 0..4, bit 0..63), LSB=0; right rotations',
        'sample_round_constants': [hex(c) for c in RC],
        'delta_words_hex': ['0x1', '0x0', '0x0', '0x0', '0x0'],
        'lambda_words_hex': ['0x4000', '0x0', '0x0', '0x0', '0x0'],
        'equivalent_core_mask_before_last_L': [[0, 14], [0, 33], [0, 42]],
        'correlation': 1,
        'status': 'Strict structural certificate, valid for any round constants',
        'proof': 'S-box dependencies are exact. OR propagation overapproximates dependencies. '
                 'Output (0,14) is absent from the cone of input (0,0). Its derivative is identically zero.',
        'sbox_input_dependencies_by_output_row': deps,
        'single_bit_input_cones_at_column_zero': cones,
        'sampling': {'generator': 'NumPy PCG64', 'seed': seed, 'pairs': n,
                     'parity_zero': n - ones, 'parity_one': ones,
                     'empirical_correlation': (n - 2 * ones) / n},
        'cross_checks': {'sbox_table_inputs': 32, 'table_vs_bitslice_permutations': 128},
        'sources': ['https://ascon.isec.tugraz.at/specification.html',
                    'https://ascon.isec.tugraz.at/images/sbox_instructions.c'],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps({'output': str(args.output), 'correlation': 1,
                      'sample': result['sampling'], 'cones': cones}, ensure_ascii=False))


if __name__ == '__main__':
    main()
