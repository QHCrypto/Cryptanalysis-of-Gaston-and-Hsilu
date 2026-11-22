"""Reproduce checks and summarize the actual six-round coefficient attempt.

No unresolved contribution is replaced by zero. Proper-period contributions
are computed using formal key derivatives, so they are linear in the key even
before they are added to the unresolved full-period contribution.
"""
import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / 'src'))
from analyze_repair import z_to_t
from gaston_6r_rank_followup import rotated_matrix, sum_matrices
from permutations import GE, GT, GW


def read(name):
    return json.loads((HERE / name).read_text())


def fold_to_t(record):
    m = record['m']
    sums = [[int(x, 16) for x in row] for row in record['folded_sums_hex']]
    return [[[z_to_t(sums[a+1][w] ^ sums[0][w], m) >> i & 1
              for a in range(2)] for w in range(2)] for i in range(m)]


def orbit_count(m):
    blocks = 64 // m
    numerator = (1 << 63) + sum(1 << (m * math.gcd(blocks, s))
                                for s in range(1, blocks))
    assert numerator % blocks == 0
    return numerator // blocks


def highest_quartic_tail_counts():
    """Exact formal cancellation at Z=L(S4), before prefix substitution."""
    n = 64
    q = [[{(1 << (n*((i+1) % 5)+j)) | (1 << (n*((i+2) % 5)+j))}
          for j in range(n)] for i in range(5)]
    a = [[q[i][(j-GE[i]) % n] for j in range(n)] for i in range(5)]
    p, v = [set() for _ in range(n)], [set() for _ in range(n)]
    for i in range(5):
        for j in range(n):
            p[j] ^= a[i][j]
            v[j] ^= a[i][(j-GT[i]) % n]
    effect = [p[(j-23) % n] ^ p[(j-24) % n] ^
              v[(j-23) % n] ^ v[(j-41) % n] for j in range(n)]
    b = [[a[i][(j-GW[i]) % n] ^ effect[(j-GW[i]) % n]
          for j in range(n)] for i in range(5)]
    result = []
    for w in range(2):
        outputs = []
        for j in range(n):
            terms = set()
            for u in b[(w+1) % 5][j]:
                for v in b[(w+2) % 5][j]:
                    if u & v:
                        continue  # Degree below 64 after the 4-round prefix.
                    term = u | v
                    if term in terms:
                        terms.remove(term)
                    else:
                        terms.add(term)
            outputs.append(terms)
        counts = []
        for i in range(4):
            terms = set()
            for j in range(n):
                if j & i == i:
                    terms ^= outputs[j]
            counts.append(len(terms))
        result.append(counts)
    return result


def main():
    checks = []
    for rounds in (4, 5):
        c = read(f'certificate_r{rounds}.json')
        reps = [int(v, 16) for v in c['representative_columns_hex']]
        expected = [[[z_to_t((v >> (64*w)) & ((1 << 64)-1)) >> i & 1
                      for v in reps] for w in range(2)] for i in range(4)]
        for suffix in ('', '_dual'):
            o = read(f'orbit_r{rounds}_m4{suffix}.json')
            assert o['coefficient_certificate'] and o['enumeration_complete']
            assert o['points_covered'] == 1 << ((1 << rounds)-1)
            assert o['low_coefficient_matrices'] == expected
            checks.append(f'Actual-word {rounds}R: orbit{suffix} agrees with direct full cubes.')

    prior = json.loads((ROOT / 'results/gaston_6r_ppm_transfer_audit.json').read_text())
    for row in prior['small_instance_experiments']:
        n, r = row['word_bits'], row['rounds']
        cols = [int(x, 16) for x in row['columns_hex']]
        even = sum_matrices([rotated_matrix(cols, s, n) for s in range(n)], (1 << n)-1)
        for a in range(2 if r < 4 else 1):
            expected = [[z_to_t((even[a*n] >> (n*w)) & ((1 << n)-1), n) >> i & 1
                         for w in range(2)] for i in range(4)]
            z = read(f'zdd_r{r}_k{a}.json')
            assert z['coefficient_certificate'] and z['low_coefficients'] == expected
        checks.append(f'ZDD vs exhaustive Boolean cubes: word width {n}, rounds {r}.')
    assert (read('zdd_r3_keyfirst_k0.json')['low_coefficients'] ==
            read('zdd_r3_k0.json')['low_coefficients'])

    short = read('orbit_r6_proper_periods_m4_dual.json')
    short_m1 = read('orbit_r6_proper_periods_m1_dual.json')
    assert short['enumeration_complete'] and short['points_covered'] == 1 << 32
    assert short['rounds'] == 6 and short['word_bits'] == 64 and short['mode'] == 'dual'
    assert not short['coefficient_certificate']
    contribution = fold_to_t(short)
    assert contribution[:1] == fold_to_t(short_m1)
    checks.append('Actual 6R proper-period contribution: independent m=1 and m=4 agree on M0.')
    for name in ('orbit_r6_m1_incomplete.json', 'orbit_r6_m4_dual_incomplete.json',
                 'zdd_r6_k0.json', 'zdd_r6_keyfirst_k0.json'):
        assert not read(name)['coefficient_certificate']

    report = {
        'target': 'Actual 64-bit-word six-round Gaston; N1=0 and parity(N0)=0; 128 output bits.',
        'target_coefficient_certificate_complete': False,
        'target_matrix_rank': None,
        'target_matrix_nontrivial_rank_lower_bound_proved': False,
        'second_cube_selected_from_actual_target_kernel': False,
        'exact_formula': {
            'top_skeleton': 'H=(Q o L)^6 without IV, round constants or affine chi terms.',
            'matrix': 'A_U=sum_(N in U) J_key H(0,N), where J is the ordinary formal derivative before Boolean key reduction.',
            'justification': 'The relevant term contains 63 distinct nonce directions and one key direction, hence has total degree64. Terms with two or more key factors or smaller total degree cannot contribute. Formal key differentiation therefore gives exactly the same full-cube matrix.',
            'periodic_key': 'For m a power of two, use g_m(z)=1+z^m+...+z^(64-m)=(z+1)^(64-m) in either key word.',
            'low_coefficients': 'A_U g_m=g_m p(z); one m-bit output block gives p modulo z^m+1=t^m.',
            'orbit_sum': 'The key direction is fixed by R_m. Sum the tangent output over each nonce orbit exactly, including every smaller orbit. An orbit with l members contributes its R_m norm, not just one output value.',
            'scope': 'These identities reduce an exact calculation. They do not supply the missing full-period sums.',
        },
        'rotation_orbits': [
            {'m': m, 'coefficient_bits': 4*m,
             'all_even_nonce_orbits': orbit_count(m),
             'full_period_64_orbits': ((1 << 63)-(1 << 32)) // (64//m),
             'log2_orbit_count': math.log2(orbit_count(m))}
            for m in (1, 2, 4, 8)
        ],
        'completed_actual_six_round_part': {
            'nonce_set': 'N0=x || x for all 32-bit x; N1=0. This is exactly the set of nonce words with period dividing32.',
            'nonce_points': 1 << 32,
            'operator': 'B=sum_(N in short-period set) J_key H(0,N); this auxiliary matrix is linear and block circulant.',
            'polynomial_mod_t4': [['1', '1+t^2'], ['1+t+t^3', '1+t+t^3']],
            'coefficient_matrices_B0_to_B3': contribution,
            'full_target_relation': 'P(t)=B(t)+C(t); C is the unresolved contribution from even-parity nonce words of least period64.',
            'unknown_C0_to_C3': None,
            'warning': 'The rank of B or the values above do not establish the rank or any coefficient of the target matrix P. C can cancel B.',
            'm4_record': 'orbit_r6_proper_periods_m4_dual.json',
            'm1_record': 'orbit_r6_proper_periods_m1_dual.json',
        },
        'formal_quartic_tail': {
            'boundary': 'Z=L(S4)',
            'counts_by_output_word_and_hasse_index_0_to_3': highest_quartic_tail_counts(),
            'meaning': 'Exact formal degree4 supports after tail/output-mask cancellation; no prefix coefficient or rank follows from these counts.',
        },
        'zdd_attempts': [read('zdd_r6_k0.json'), read('zdd_r6_keyfirst_k0.json')],
        'verification': checks,
        'remaining_obstruction': 'No exact parity calculation has completed for the full-period64 contribution. The finite enumeration and ZDD attempts here do not yield the actual 6R kernel.',
    }
    path = ROOT / 'results/gaston_6r_low_coefficients_attempt.json'
    path.write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps({
        'path': str(path),
        'target_rank': None,
        'checks': checks,
        'proper_period_coefficients': contribution,
        'formal_quartic_counts': report['formal_quartic_tail']['counts_by_output_word_and_hasse_index_0_to_3'],
    }, indent=2))


if __name__ == '__main__':
    main()
