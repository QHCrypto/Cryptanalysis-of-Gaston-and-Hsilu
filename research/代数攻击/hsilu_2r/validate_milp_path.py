"""Validate an exported MILP path against actual two-round Hsilu polynomials.

The full quadratic ANF in active variables is checked, not just degree tags.
Random capacity samples do not establish global condition reachability/rank.
"""
import argparse
import ast
import collections
import ctypes
import json
from pathlib import Path
import random
import subprocess
import verify_path_algebra as v

def read_path(path):
    result = {}
    for line in path.read_text().splitlines():
        if ': ' not in line:
            continue
        key, value = line.split(': ', 1)
        try:
            result[key] = ast.literal_eval(value)
        except (SyntaxError, ValueError):
            result[key] = value
    return result

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('path', type=Path)
    ap.add_argument('--samples', type=int, default=16)
    args = ap.parse_args()
    if args.samples <= 0:
        ap.error('--samples must be positive')
    data = read_path(args.path)
    assert data['mrm_parameters'] == v.DATA['mrm_parameters']
    assert data['rounds'] == 2 and data['rate_lanes'] == [0]
    v.U = [63-j for j in data['free_variables_by_lane'][0]]
    v.G = [63-j for j in data['guessed_by_lane'][0]]
    v.S = [63-j for j in data['selected_outputs_by_lane'][0]]
    v.UMASK = sum(1 << j for j in v.U)
    assert set(v.U).isdisjoint(v.G) and sorted(v.U + v.G) == list(range(64))
    constraints = []
    all_conditions = []
    for value in (0, 1):
        for stage, rows in data[f'conditions_eq_{value}'].items():
            for row, js in rows.items():
                for j in js:
                    item = (stage, row, 63-j, value)
                    all_conditions.append(item)
                    if stage == 'X2' or row == 0:
                        constraints.append(item)
    assert len(all_conditions) == data['condition_count']
    # assert len(v.U) + len(all_conditions) <= 64
    assert len(v.S) == len(v.U)

    ref = v.HERE.parents[1] / 'hsilu_ref'
    subprocess.run(['g++', '-std=c++17', '-O2', '-shared', '-fPIC', '-I'+str(ref),
                    str(v.HERE/'verify_reference_wrapper.cpp'), str(ref/'hsilu.cpp'),
                    '-o', str(v.HERE/'verify_reference_wrapper.so')], check=True)
    lib = ctypes.CDLL(str(v.HERE/'verify_reference_wrapper.so'))
    lib.verify_round.argtypes = [ctypes.POINTER(ctypes.c_uint64), ctypes.c_uint64]

    def condition_values(state):
        mid = v.rnd(state, 0x5a)
        return [((state if stage == 'X1' else mid)[row] >> j) & 1
                for stage, row, j, _ in constraints]

    random.seed(20261013)
    histogram = collections.Counter()
    constraint_rank_histogram = collections.Counter()
    unsat = 0
    witness = None
    for sample in range(args.samples):
        state = [random.getrandbits(64) for _ in range(5)]
        state[0] &= ~v.UMASK
        for stage, row, j, value in all_conditions:
            if stage == 'X1' and row > 0:
                state[row] = (state[row] & ~(1 << j)) | (value << j)
        baseline = condition_values(state)
        rows = [0] * len(constraints)
        for k, j in enumerate(v.G):
            changed = state.copy()
            changed[0] ^= 1 << j
            for r, (x, y) in enumerate(zip(baseline, condition_values(changed))):
                rows[r] |= (x ^ y) << k
        # First-round state is exactly affine in all message bits.
        for j in v.U:
            changed = state.copy()
            changed[0] ^= 1 << j
            assert condition_values(changed) == baseline, 'Condition depends on an active variable'
        rhs = [x ^ item[3] for x, item in zip(baseline, constraints)]
        solution = v.solve(rows, rhs, len(v.G))
        constraint_rank_histogram[v.rank(rows)] += 1
        if solution is None:
            unsat += 1
            continue
        state[0] ^= sum(((solution >> k) & 1) << j for k, j in enumerate(v.G))
        assert condition_values(state) == [item[3] for item in constraints]
        base, columns, quadratic, output_rows = v.affine_data(state)
        bad = [63-j for j in v.S if quadratic >> j & 1]
        assert not bad, f'Selected outputs are nonlinear: {bad}'
        rank = v.rank(output_rows)
        histogram[rank] += 1
        # Compare a random complete active assignment to the compiled reference.
        assignment = random.getrandbits(len(v.U))
        candidate = state.copy()
        candidate[0] ^= sum(((assignment >> k) & 1) << j for k, j in enumerate(v.U))
        native = (ctypes.c_uint64 * 5)(*candidate)
        for rc in (0x5a, 0x4b):
            lib.verify_round(native, rc)
        assert list(native) == v.forward(candidate)
        predicted = base
        for k, column in enumerate(columns):
            if assignment >> k & 1:
                predicted ^= column
        assert all((predicted ^ native[0]) >> j & 1 == 0 for j in v.S)
        if witness is None and rank == len(v.U):
            recovered = v.solve(output_rows, [((native[0] ^ base) >> j) & 1 for j in v.S], len(v.U))
            assert recovered == assignment
            witness = {'scope': 'planted target at chosen capacity/background, not arbitrary XOF preimage',
                       'native_input_hex': [f'{w:016x}' for w in candidate],
                       'first_block_target_hex': f'{native[0]:016x}', 'rank': rank}

    assert sum(histogram.values()) > 0, 'No consistent sampled capacity; feasibility unresolved'
    output = {'input_path': str(args.path.resolve()), 'sample_count': args.samples,
              'sampled_unsatisfiable_capacities': unsat,
              'background_condition_rank_histogram': dict(constraint_rank_histogram),
              'output_rank_histogram': dict(histogram), 'selected_output_quadratic_terms': 'all zero in every consistent sample',
              'cpp_reference_checks': sum(histogram.values()), 'witness': witness,
              'limitations': 'No full 128-bit XOF preimage, no proof of condition independence or sponge prefix reachability.'}
    filename = args.path.with_suffix('.algebra_validation.json')
    filename.write_text(json.dumps(output, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(output, ensure_ascii=False, indent=2))
    print('Saved:', filename)

if __name__ == '__main__':
    main()
