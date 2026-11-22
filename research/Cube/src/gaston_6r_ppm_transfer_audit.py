"""Audit the 2025/1029 transfer to six-round Gaston.

The total degree is at most 64, so a 64D cube is key-independent. For a
63D cube the superpoly is affine in all 128 key bits. Its key-linear part
comes solely from total-degree-64 terms, and therefore from the quadratic
round skeleton H=(Q o L)^6, without IV or round constants.

For the cube N0[0]=0, N1=0, reconstruct all key columns from two representative
key bits and all codimension-one nonce facets, using rotation covariance.
Full-size costs are counted, not executed. Full-size matrix rank is unknown.
"""
import json
import math
import random
from pathlib import Path

import numpy as np

from engine import Engine
from permutations import GE, GT, GW, GASTON_KAT, RC, gaston_linear
from structure import polar, rank

ROOT = Path(__file__).resolve().parents[1]
IV = 0x00001000808C0001


def rol(a, shift, n):
    shift %= n
    return ((a << shift) | (a >> ((-shift) % n))) & ((1 << n) - 1)


def evaluate(nonce0, nonce1, key, n, rounds, *, skeleton=False, iv=IV):
    nonce0 = np.asarray(nonce0, dtype=np.uint64)
    mask = np.uint64((1 << n) - 1)
    state = [np.full_like(nonce0, (0 if skeleton else iv) & int(mask)),
             np.full_like(nonce0, key[0]), np.full_like(nonce0, key[1]),
             nonce0.copy(), np.full_like(nonce0, nonce1)
             if np.isscalar(nonce1) else np.asarray(nonce1, dtype=np.uint64)]

    def rotate(x, s):
        s %= n
        return ((x << np.uint64(s)) | (x >> np.uint64((-s) % n))) & mask

    for rc in RC[-rounds:]:
        a = [rotate(x, s) for x, s in zip(state, GE)]
        p, q = np.zeros_like(nonce0), np.zeros_like(nonce0)
        for x, s in zip(a, GT):
            p ^= x
            q ^= rotate(x, s)
        effect = rotate(p ^ rotate(p, 1) ^ q ^ rotate(q, 18), 23)
        a = [rotate(x ^ effect, s) for x, s in zip(a, GW)]
        if skeleton:
            state = [a[(i + 1) % 5] & a[(i + 2) % 5] for i in range(5)]
        else:
            a[0] ^= np.uint64(rc) & mask
            state = [a[i] ^ ((mask ^ a[(i + 1) % 5]) & a[(i + 2) % 5])
                     for i in range(5)]
    return state


def packed_sum(state, n, words=2):
    return sum(int(np.bitwise_xor.reduce(state[i])) << (n * i)
               for i in range(words))


def rotate_output(value, s, n, words=2):
    return sum(rol((value >> (n * i)) & ((1 << n) - 1), s, n) << (n * i)
               for i in range(words))


class FacetStream:
    """Streaming truncated Boolean transform in binary input order.

    A block on ell variables stores its full sum and ell facets with one
    coordinate fixed to zero. Merging left/right uses ell+1 vector XORs.
    A 2^n-point stream needs 2^(n+1)-n-2 vector XORs and O(n^2) vectors.
    """
    def __init__(self):
        self.stack = []
        self.xors = 0
        self.peak_vectors = 0
        self.points = 0

    def push(self, value):
        block, level = [value], 0
        self.points += 1
        self.peak_vectors = max(self.peak_vectors,
                                1 + sum(len(x) for x in self.stack if x is not None))
        while level < len(self.stack) and self.stack[level] is not None:
            left = self.stack[level]
            self.stack[level] = None
            left_full = left[0]
            for i in range(len(left)):
                left[i] ^= block[i]
                self.xors += 1
            left.append(left_full)
            block = left
            level += 1
        if level == len(self.stack):
            self.stack.append(block)
        else:
            self.stack[level] = block

    def finish(self, n):
        assert self.points == 1 << n
        assert len(self.stack) == n + 1 and all(v is None for v in self.stack[:n])
        assert self.xors == 2 * (1 << n) - n - 2
        return self.stack[n]


def apply_columns(columns, key):
    result = 0
    while key:
        bit = key & -key
        result ^= columns[bit.bit_length() - 1]
        key ^= bit
    return result


def small_instance(n, rounds):
    assert n == 1 << rounds
    nonce = np.arange(1 << n, dtype=np.uint64)
    h0 = evaluate(nonce, 0, (0, 0), n, rounds, skeleton=True)
    h1 = evaluate(nonce, 0, (1, 0), n, rounds, skeleton=True)
    h2 = evaluate(nonce, 0, (0, 1), n, rounds, skeleton=True)
    differences = [[h1[i] ^ h0[i] for i in range(2)],
                   [h2[i] ^ h0[i] for i in range(2)]]
    stream = FacetStream()
    for index in range(1 << n):
        value = sum(int(differences[a][w][index]) << (n * (2 * a + w))
                    for a in range(2) for w in range(2))
        stream.push(value)
    combined = stream.finish(n)
    assert combined[0] == 0  # The full nD cube of a key difference is zero.
    representatives = [[(combined[t + 1] >> (2 * n * a)) & ((1 << (2 * n)) - 1)
                        for t in range(n)] for a in range(2)]
    for a in range(2):
        for t in range(n):
            selected = (nonce & np.uint64(1 << t)) == 0
            assert representatives[a][t] == packed_sum(
                [word[selected] for word in differences[a]], n)
    columns = [rotate_output(representatives[a][(-j) % n], j, n)
               for a in range(2) for j in range(n)]

    points = np.arange(1 << (n - 1), dtype=np.uint64) << np.uint64(1)
    constant = packed_sum(evaluate(points, 0, (0, 0), n, rounds), n)
    for j, expected in enumerate(columns):
        key = ((1 << j), 0) if j < n else (0, 1 << (j - n))
        actual = packed_sum(evaluate(points, 0, key, n, rounds), n) ^ constant
        assert actual == expected
    # Different IV: the constant may change but the highest-degree key matrix cannot.
    other_constant = packed_sum(evaluate(points, 0, (0, 0), n, rounds, iv=0), n)
    for j in [0, n // 2, n, 2 * n - 1]:
        key = ((1 << j), 0) if j < n else (0, 1 << (j - n))
        assert packed_sum(evaluate(points, 0, key, n, rounds, iv=0), n) ^ other_constant == columns[j]
    rng = random.Random(20261020 + n)
    keys = list(range(1 << (2 * n))) if n == 4 else [rng.getrandbits(2 * n) for _ in range(32)]
    for key in keys:
        actual = packed_sum(evaluate(points, 0, (key & ((1 << n) - 1), key >> n), n, rounds), n)
        assert actual == constant ^ apply_columns(columns, key)
    full_cube = packed_sum(evaluate(nonce, 0, (0, 0), n, rounds), n)
    for key in keys[:8]:
        assert packed_sum(evaluate(nonce, 0, (key & ((1 << n) - 1), key >> n), n, rounds), n) == full_cube
    return {
        'word_bits': n, 'rounds': rounds, 'cube_dimension': n - 1,
        'key_bits': 2 * n, 'feature_bits': 2 * n,
        'all_key_columns_cross_checked': len(columns),
        'rank_one_output_word': rank([x & ((1 << n) - 1) for x in columns]),
        'rank_two_output_words': rank(columns),
        'random_or_exhaustive_key_checks': len(keys),
        'keys_exhaustive': n == 4,
        'different_iv_key_columns_checked': 4,
        'full_dimension_cube_key_independence_checks': 8,
        'stream_vector_xors': stream.xors,
        'stream_peak_packed_vectors': stream.peak_vectors,
        'columns_hex': [f'{v:x}' for v in columns],
    }


def full_word_reduced_round_checks():
    engine = Engine()
    try:
        for source, target in GASTON_KAT:
            result = evaluate(np.array([source[3]], dtype=np.uint64), source[4],
                              (source[1], source[2]), 64, 12, iv=source[0])
            assert [int(x[0]) for x in result] == target
        # Actual 64-bit words, four rounds, 15D cube. This verifies the
        # highest-total-degree extraction independently against the C++ engine.
        coordinates = [7, 8, 14, 24, 25, 48, 50, 54, 66, 69, 85, 90, 102, 113, 118]
        index = np.arange(1 << len(coordinates), dtype=np.uint64)
        n0, n1 = np.zeros_like(index), np.zeros_like(index)
        for i, p in enumerate(coordinates):
            contribution = ((index >> np.uint64(i)) & np.uint64(1)) << np.uint64(p % 64)
            if p < 64:
                n0 ^= contribution
            else:
                n1 ^= contribution
        base = engine.cube('gaston', 4, IV, (0, 0), coordinates)
        h0 = packed_sum(evaluate(n0, n1, (0, 0), 64, 4, skeleton=True), 64)
        checks = [0, 1, 17, 63, 64, 65, 97, 127]
        for j in checks:
            key = ((1 << j), 0) if j < 64 else (0, 1 << (j - 64))
            observed = engine.cube('gaston', 4, IV, key, coordinates)
            expected = (observed[0] ^ base[0]) | ((observed[1] ^ base[1]) << 64)
            value = packed_sum(evaluate(n0, n1, key, 64, 4, skeleton=True), 64) ^ h0
            assert value == expected
            s = j % 64
            representative = (1, 0) if j < 64 else (0, 1)
            rotated_n0, rotated_n1 = rol(n0, -s, 64), rol(n1, -s, 64)
            r0 = packed_sum(evaluate(rotated_n0, rotated_n1, (0, 0), 64, 4, skeleton=True), 64)
            r1 = packed_sum(evaluate(rotated_n0, rotated_n1, representative, 64, 4, skeleton=True), 64)
            assert rotate_output(r0 ^ r1, s, 64) == expected
        return {'word_bits': 64, 'rounds': 4, 'cube_dimension': 15,
                'original_paper_known_answer_vectors': len(GASTON_KAT),
                'independent_cpp_linear_coefficient_checks': len(checks),
                'rotation_reconstruction_checks': len(checks),
                'nonce_coordinates': coordinates}
    finally:
        engine.close()


def first_round_structure():
    nonce = [gaston_linear([0, 0, 0, 1 << j, 0]) for j in range(64)]
    key = [gaston_linear([0, (1 << j) if j < 64 else 0,
                         (1 << (j - 64)) if j >= 64 else 0, 0, 0]) for j in range(128)]
    columns = [sum(polar(k, u) << (320 * i) for i, u in enumerate(nonce[1:])) for k in key]
    witness = polar(nonce[1], nonce[2])
    assert witness
    return {'cube_nonce_word': 0, 'cube_coordinates': list(range(1, 64)),
            'nonzero_quadratic_nonce_pair': [1, 2], 'quadratic_coefficient_hex': f'{witness:080x}',
            'first_round_mixed_key_coefficient_rank': rank(columns),
            'meaning': 'First round is not affine on this cube. Rank concerns first-round mixed coefficients only, not the final superpoly.'}


def quartic_tail():
    """Exact formal quartics from boundary Z=L(S4), before round-five chi.

    These counts are after Boolean ANF cancellation, BEFORE any prefix-trail
    screening. A quartic in Z has degree at most 64 in original inputs.
    """
    n = 64
    q = [[{(1 << (n * ((i + 1) % 5) + j)) |
           (1 << (n * ((i + 2) % 5) + j))} for j in range(n)] for i in range(5)]
    a = [[q[i][(j - GE[i]) % n] for j in range(n)] for i in range(5)]
    p, v = [set() for _ in range(n)], [set() for _ in range(n)]
    for i in range(5):
        for j in range(n):
            p[j] ^= a[i][j]
            v[j] ^= a[i][(j - GT[i]) % n]
    effect = [p[(j - 23) % n] ^ p[(j - 24) % n] ^
              v[(j - 23) % n] ^ v[(j - 41) % n] for j in range(n)]
    b = [[a[i][(j - GW[i]) % n] ^ effect[(j - GW[i]) % n]
          for j in range(n)] for i in range(5)]
    records, all_classes = [], set()
    for i in range(2):
        terms = set()
        for u in b[(i + 1) % 5][0]:
            for v in b[(i + 2) % 5][0]:
                if u | v in terms:
                    terms.remove(u | v)
                else:
                    terms.add(u | v)
        quartics = [tuple(j for j in range(5 * n) if (m >> j) & 1)
                    for m in terms if m.bit_count() == 4]
        classes = {min(tuple(sorted(n * (j // n) + (j % n + s) % n for j in t))
                       for s in range(n)) for t in quartics}
        all_classes.update(classes)
        records.append({'output_word': i, 'bit_zero_quartics': len(quartics),
                        'rotation_classes': len(classes)})
    return {
        'prefix_boundary': 'Z=L(S4), before fifth-round chi; every Z bit has total degree at most 16.',
        'formal_highest_tail': 'Q(L(Q(Z)))',
        'word_counts': records,
        'two_word_union_rotation_classes': len(all_classes),
        'prefix_division_property_screen_completed': False,
        'key_linear_coefficient_rule': 'For [N^I k_j] with |I|=63, distribute the 63 nonce variables disjointly as 15+16+16+16 among four prefix factors; the size-15 factor supplies k_j.',
        'why_no_all_key_matching_needed': 'The resulting superpoly is affine in the full key; solve a linear system instead of evaluating a nonlinear feature table on 2^128 keys.',
        'scope': 'Formal counts only. Not counts of nonzero full-cube coefficients, not a completed MILP screen, and not the cost basis for the rotation-interpolation schedule.',
    }


def main():
    small = [small_instance(n, r) for n, r in [(4, 2), (8, 3), (16, 4)]]
    full_word = full_word_reduced_round_checks()
    n, points = 64, 1 << 64
    report = {
        'model': 'Six-round Gaston in the existing Ascon-AEAD128-style research mode; (IV,K0,K1,N0,N1), first two output words observable',
        'iv_hex': f'{IV:016x}',
        'round_constants': RC[-6:],
        'strict_degree_conclusions': {
            'total_key_nonce_degree_upper_bound': 64,
            '64D_cube_key_degree_upper_bound': 0,
            '63D_cube_key_degree_upper_bound': 1,
            '64D_cube_key_recovery_via_feature_matching': 'No information: cube sum is independent of every key bit.',
            '63D_matrix_key_bits': 128,
            '63D_matrix_key_linear_part': 'Total-degree-64 component; IV and round constants do not affect it.',
        },
        'first_round_structure': first_round_structure(),
        'four_factor_method_transfer': quartic_tail(),
        'rotation_identity': 'For I={1,...,63} in N0: A_(a,j)=ROL_j B_(a,-j mod 64), where B_(a,t)=XOR_(N0[t]=0) (H(K_a=1,N0)+H(K=0,N0)), H=(Q o L)^6.',
        'coefficient_extraction': 'B for both representative key bits is computed in one 2^64-point stream; each point evaluates H at key 0 and the two representatives.',
        'streaming_transform': {
            'record': 'Full sum plus one zero-fixed facet sum per processed coordinate.',
            'merge': '(a_full,a_facets),(b_full,b_facets) -> (a_full XOR b_full, a_facets XOR b_facets, a_full)',
            'packed_vector_bits': 256,
            'vector_xors': 2 * points - n - 2,
            'peak_stack_vectors_bound': n * (n + 1) // 2 + 1,
            'core_stack_bits': 256 * (n * (n + 1) // 2 + 1),
            'working_memory_conservative_bound_bits': 1 << 20,
            'memory_excludes': 'Language object overhead; bound is for a packed implementation, not NumPy demo arrays.',
        },
        'full_size_cost_schedule_not_executed': {
            'matrix_skeleton_evaluations': 3 * points,
            'constant_term_full_permutation_evaluations': points // 2,
            'online_full_permutation_queries': points // 2,
            'total_evaluations_upper_bound_in_six_round_units': 4 * points,
            'total_evaluations_log2': 66,
            'additional_64bit_xors': 14 * points - 4 * n - 8,
            'additional_64bit_xors_log2': math.log2(14 * points - 4 * n - 8),
            'matrix_reindex_and_solve': 'Polynomial cost in 128, negligible beside 2^64.',
            'data_log2': 63,
            'candidate_space': 'Exactly 2^(128-rho) keys satisfy the affine cube equations; rho=rank(A).',
            'verification_upper_bound': 'O(2^(128-rho)) full-key checks on known observations.',
            'rho': None,
            'old_direct_interpolation_total_evaluations': 130 * (1 << 63),
            'old_direct_interpolation_total_log2': math.log2(130 * (1 << 63)),
            'evaluation_count_speedup': 130 / 8,
        },
        'small_instance_experiments': small,
        'actual_word_size_reduced_round_checks': full_word,
        'limitations': [
            'The six-round 128x128 coefficient matrix has not been evaluated and its rank is unknown.',
            '2^66 is a proved evaluation schedule for obtaining equations and online cube sums, not an unconditional full-key recovery claim.',
            'Reduced-word matrix ranks do not prove anything about the actual six-round rank.',
            'A 64-bit feature observes at most 64 independent equations in 128 key bits; this differs from the previous Hsilu local-key feature.',
            'The schedule uses total-degree truncation and rotation, not a completed Gaston division-property/MILP quartic screen.',
            'No literature novelty claim is made for the symmetry-based interpolation schedule.',
        ],
    }
    path = ROOT / 'results/gaston_6r_ppm_transfer_audit.json'
    path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({k: v for k, v in report.items() if k != 'small_instance_experiments'}, indent=2))
    print(json.dumps([{k: v for k, v in row.items() if k != 'columns_hex'} for row in small], indent=2))


if __name__ == '__main__':
    main()
