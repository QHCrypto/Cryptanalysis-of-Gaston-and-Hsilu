"""Reproducible algebra checks with explicit path-to-native bit mappings.

X1 is the first round input, X2 the second round input, X3 its output.
No sponge initialization or padding is assumed here.
"""
import ast
import collections
import ctypes
import itertools
import json
from pathlib import Path
import random
import subprocess

HERE = Path(__file__).resolve().parent
DATA = {}
for line in (HERE / 'hsilu_2round_best_attack.txt').read_text().splitlines():
    key, value = line.split(': ', 1)
    try:
        DATA[key] = ast.literal_eval(value)
    except (ValueError, SyntaxError):
        DATA[key] = value
MASK = (1 << 64) - 1
PARAM = DATA['mrm_parameters']
U = DATA['free_variables_by_lane'][0]
G = DATA['guessed_by_lane'][0]
S = DATA['selected_outputs_by_lane'][0]
UMASK = sum(1 << j for j in U)

def rot(x, n):
    return ((x << n) | (x >> ((64 - n) % 64))) & MASK

def linear(a):
    x = list(a)
    for j in range(5, 18):
        x.append(x[PARAM['n1'][j]] ^ rot(x[PARAM['n2'][j]], PARAM['m'][j]))
    return [rot(x[13 + i], PARAM['s'][i]) for i in range(5)]

def rnd(a, rc):
    a = list(a)
    a[0] ^= rc
    return linear([a[i] ^ ((~a[(i + 1) % 5]) & a[(i + 2) % 5]) for i in range(5)])

def forward(a, constants=(0x5a, 0x4b)):
    for rc in constants:
        a = rnd(a, rc)
    return a

def rank(rows):
    pivots = {}
    for row in rows:
        while row:
            p = row.bit_length() - 1
            if p in pivots:
                row ^= pivots[p]
            else:
                pivots[p] = row
                break
    return len(pivots)

def solve(rows, rhs, n):
    pivots = {}
    for row, b in zip(rows, rhs):
        while row:
            p = row.bit_length() - 1
            if p in pivots:
                row ^= pivots[p][0]
                b ^= pivots[p][1]
            else:
                pivots[p] = (row, b)
                break
        if not row and b:
            return None
    result = 0
    for p in sorted(pivots):
        row, b = pivots[p]
        result |= (b ^ ((row & result).bit_count() & 1)) << p
    return result

def affine_data(a):
    base = forward(a)[0]
    cols = []
    for j in U:
        b = a.copy()
        b[0] ^= 1 << j
        cols.append(forward(b)[0] ^ base)
    qmask = 0
    for k, l in itertools.combinations(range(len(U)), 2):
        b = a.copy()
        b[0] ^= (1 << U[k]) | (1 << U[l])
        qmask |= forward(b)[0] ^ base ^ cols[k] ^ cols[l]
    rows = [sum(((c >> j) & 1) << k for k, c in enumerate(cols)) for j in S]
    return base, cols, qmask, rows

def validate_convention(name, mapping, sample_count=256):
    global U, G, S, UMASK
    path_u = DATA['free_variables_by_lane'][0]
    path_s = DATA['selected_outputs_by_lane'][0]
    U = [mapping(j) for j in path_u]
    G = [mapping(j) for j in DATA['guessed_by_lane'][0]]
    S = [mapping(j) for j in path_s]
    UMASK = sum(1 << j for j in U)
    condition_bit = mapping(17)
    active_columns = [linear([1 << j, 0, 0, 0, 0]) for j in U]
    condition_active_support = [path_u[k] for k, col in enumerate(active_columns) if (col[1] >> condition_bit) & 1]
    # Exact second-degree coefficients follow solely from these columns.
    # Constants cancel in chi's second finite difference.
    quadratic_by_pair = {}
    for k, l in itertools.combinations(range(len(U)), 2):
        c, d = active_columns[k], active_columns[l]
        q = linear([(c[(i+1)%5] & d[(i+2)%5]) ^ (d[(i+1)%5] & c[(i+2)%5]) for i in range(5)])[0]
        if q:
            quadratic_by_pair[f'{path_u[k]},{path_u[l]}'] = [j for j in range(64) if q >> mapping(j) & 1]
    nonlinear = set(itertools.chain.from_iterable(quadratic_by_pair.values()))
    exact_affine = [j for j in range(64) if j not in nonlinear]
    before_hist = collections.Counter()
    after_hist = collections.Counter()
    all22_hist = collections.Counter()
    condition_support_hist = collections.Counter()
    samples = []
    witness = None
    for trial in range(sample_count):
        a = [random.getrandbits(64) for _ in range(5)]
        a[0] &= ~UMASK
        a[1] &= ~UMASK
        a[4] |= UMASK
        mid = rnd(a, 0x5a)
        coefficients = []
        for j in G:
            b = a.copy()
            b[0] ^= 1 << j
            if ((rnd(b, 0x5a)[1] ^ mid[1]) >> condition_bit) & 1:
                coefficients.append(j)
        before = affine_data(a)
        before_hist[rank(before[3])] += 1
        condition_support_hist[len(coefficients)] += 1
        if (mid[1] >> condition_bit) & 1:
            if not coefficients:
                samples.append({'trial': trial, 'condition_unsatisfiable': True})
                continue
            a[0] ^= 1 << coefficients[0]
        assert ((rnd(a, 0x5a)[1] >> condition_bit) & 1) == 0
        base, cols, qmask, rows = affine_data(a)
        affine_positions = [j for j in range(64) if not (qmask >> mapping(j)) & 1]
        assert affine_positions == exact_affine
        rows22 = [sum(((c >> mapping(j)) & 1) << k for k, c in enumerate(cols)) for j in DATA['final_linearizable_x3_by_lane'][0]]
        after_hist[rank(rows)] += 1
        all22_hist[rank(rows22)] += 1
        sample = {'trial': trial, 'x2_condition_background_support': [j for j in DATA['guessed_by_lane'][0] if mapping(j) in coefficients],
                  'rank_selected': rank(rows), 'rank_all22': rank(rows22),
                  'nonlinear_selected': [j for j in path_s if (qmask >> mapping(j)) & 1]}
        samples.append(sample)
        if witness is None and not sample['nonlinear_selected'] and rank(rows) == len(U):
            planted_u = random.getrandbits(len(U))
            planted = a.copy()
            planted[0] ^= sum(((planted_u >> k) & 1) << j for k, j in enumerate(U))
            target = forward(planted)[0]
            solution = solve(rows, [((target ^ base) >> j) & 1 for j in S], len(U))
            candidate = a.copy()
            candidate[0] ^= sum(((solution >> k) & 1) << j for k, j in enumerate(U))
            actual = forward(candidate)[0]
            assert actual == target
            assert ((rnd(candidate, 0x5a)[1] >> condition_bit) & 1) == 0
            witness = {'input_base_words_hex': [f'{w:016x}' for w in a],
                       'planted_message_hex': f'{planted[0]:016x}', 'target_hex': f'{target:016x}',
                       'solved_message_hex': f'{candidate[0]:016x}', 'candidate_digest_hex': f'{actual:016x}',
                       'selected_bits_match': True, 'full_64_bits_match': True, 'rank': rank(rows),
                       'input_candidate_words_hex': [f'{w:016x}' for w in candidate],
                       'x2_candidate_words_hex': [f'{w:016x}' for w in rnd(candidate, 0x5a)],
                       'output_candidate_words_hex': [f'{w:016x}' for w in forward(candidate)]}
    return {'name': name, 'condition_active_support': condition_active_support,
            'exact_affine_positions': exact_affine, 'exact_nonlinear_selected': [j for j in path_s if j in nonlinear],
            'quadratic_coefficients_by_active_pair': quadratic_by_pair,
            'sample_count': sample_count, 'rank_before_x2_condition': dict(before_hist),
            'rank_after_x2_condition': dict(after_hist), 'rank_all22_after_x2_condition': dict(all22_hist),
            'x2_background_support_size_histogram': dict(condition_support_hist),
            'witness': witness, 'samples': samples}

def reduced_target_experiment(case):
    """Independent 32-bit target on 21 selected + 11 other output positions.

    Capacity is chosen freely here, not generated by a sponge prefix.
    This is explicitly not an arbitrary 64-bit preimage computation.
    """
    global U, G, S, UMASK
    U = [63-j for j in DATA['free_variables_by_lane'][0]]
    G = [63-j for j in DATA['guessed_by_lane'][0]]
    S = [63-j for j in DATA['selected_outputs_by_lane'][0]]
    UMASK = sum(1 << j for j in U)
    background = [int(x,16) for x in case['witness']['input_base_words_hex']]
    selected_path = DATA['selected_outputs_by_lane'][0]
    extra_path = [j for j in range(64) if j not in selected_path][:11]
    target_positions = selected_path + extra_path
    target_mask = sum(1 << (63-j) for j in target_positions)
    target = random.getrandbits(64) & target_mask
    rank_counts = collections.Counter()
    candidate_count = 0
    for trial in range(1, 100001):
        a = background.copy()
        a[0] = random.getrandbits(64) & ~UMASK
        # Coefficient of path m[24] in X2[1][17] is identically one.
        if rnd(a, 0x5a)[1] >> (63-17) & 1:
            a[0] ^= 1 << (63-24)
        base = forward(a)[0]
        columns = []
        for j in U:
            b = a.copy()
            b[0] ^= 1 << j
            columns.append(forward(b)[0] ^ base)
        rows = [sum(((c >> j) & 1) << k for k,c in enumerate(columns)) for j in S]
        rk = rank(rows)
        rank_counts[rk] += 1
        if rk != 21:
            continue
        solution = solve(rows, [((target ^ base) >> j) & 1 for j in S], 21)
        a[0] ^= sum(((solution >> k) & 1) << j for k,j in enumerate(U))
        digest = forward(a)[0]
        assert all(((digest ^ target) >> j) & 1 == 0 for j in S)
        candidate_count += 1
        if (digest & target_mask) == target:
            assert rnd(a, 0x5a)[1] >> (63-17) & 1 == 0
            return {'status': 'found', 'scope': 'independent 32-bit projected target; chosen capacity',
                    'trial_count': trial, 'candidate_count': candidate_count,
                    'rank_histogram': dict(rank_counts), 'path_target_positions': target_positions,
                    'native_target_mask_hex': f'{target_mask:016x}',
                    'native_target_value_hex': f'{target:016x}',
                    'native_input_words_hex': [f'{w:016x}' for w in a],
                    'native_full_digest_hex': f'{digest:016x}', 'matched_bits': 32}
    return {'status': 'not_found_within_limit', 'trial_count': 100000}

def main():
    random.seed(20261013)
    wrapper = HERE / 'verify_reference_wrapper.cpp'
    wrapper.write_text('#include "hsilu.h"\nvoid hsilu_round(hsilu_state_t*, uint64_t);\nextern "C" void verify_round(uint64_t* p, uint64_t c) { hsilu_state_t a; for(int i=0;i<5;i++)a.x[i]=p[i]; hsilu_round(&a,c); for(int i=0;i<5;i++)p[i]=a.x[i]; }\n')
    libpath = HERE / 'verify_reference_wrapper.so'
    subprocess.run(['g++', '-std=c++17', '-O2', '-shared', '-fPIC', '-I'+str(HERE.parents[1]/'hsilu_ref'), str(wrapper), str(HERE.parents[1]/'hsilu_ref/hsilu.cpp'), '-o', str(libpath)], check=True)
    lib = ctypes.CDLL(str(libpath))
    lib.verify_round.argtypes = [ctypes.POINTER(ctypes.c_uint64), ctypes.c_uint64]
    for _ in range(100):
        a = [random.getrandbits(64) for _ in range(5)]
        rc = random.getrandbits(8)
        c = (ctypes.c_uint64 * 5)(*a)
        lib.verify_round(c, rc)
        assert list(c) == rnd(a, rc)

    cases = [validate_convention('native_lsb_j', lambda j: j),
             validate_convention('path_j_to_native_63_minus_j', lambda j: 63-j),
             validate_convention('path_j_to_native_minus_j_mod64', lambda j: -j%64)]
    reduced = reduced_target_experiment(cases[1])
    # Verify concrete witnesses again using the C++ reference, not just Python.
    for c in cases:
        if c['witness']:
            w = c['witness']
            state = (ctypes.c_uint64 * 5)(*[int(x,16) for x in w['input_candidate_words_hex']])
            for rc in (0x5a, 0x4b):
                lib.verify_round(state, rc)
            assert list(state) == [int(x,16) for x in w['output_candidate_words_hex']]
    if reduced['status'] == 'found':
        state = (ctypes.c_uint64 * 5)(*[int(x,16) for x in reduced['native_input_words_hex']])
        for rc in (0x5a, 0x4b):
            lib.verify_round(state, rc)
        assert state[0] == int(reduced['native_full_digest_hex'],16)
    out = {'assumptions': {'round_constants': ['0x5a', '0x4b'],
                           'state_meanings': 'X1=input, X2=after round 1, X3=after round 2',
                           'scope': 'local two-round permutation, digest=row 0'},
           'reference_random_tests': 100, 'capacity_conditions': 42, 'cases': cases,
           'reduced_target_experiment': reduced,
           'x2_condition_exact_message_expression_path_indexing':
           'm[24] + (1+c4[14])*m[14] + c1[1]*m[1] + c1[61]*m[61] + k(c,rc) = 0'}
    (HERE / 'path_algebra_validation_results.json').write_text(json.dumps(out, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({**{k:v for k,v in out.items() if k!='cases'}, 'cases': [
        {k:v for k,v in c.items() if k not in ['samples', 'quadratic_coefficients_by_active_pair']} for c in cases
    ]}, ensure_ascii=False, indent=2))

if __name__ == '__main__':
    main()
