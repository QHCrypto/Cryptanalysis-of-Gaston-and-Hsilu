"""Verify one concrete, background-only linearization of the eight X3 conditions.

This is a construction for the saved capacity sample, not an optimum proof and
not a complete XOF preimage algorithm. All arithmetic is over GF(2).
"""
import collections
import argparse
import json
import random
import analyze_three_round_path as a


def evaluate(poly, x):
    return sum((x & monomial) == monomial for monomial in poly) & 1


def substitute_background(t, particular, basis):
    for k, v in enumerate(basis):
        if t >> k & 1:
            particular ^= v
    return particular


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--case', choices=['first', 'cheap'], default='first')
    args = parser.parse_args()
    detail = json.loads((a.HERE / f'{args.case}_capacity_algebra_details.json').read_text())
    if args.case == 'first':
        sample = json.loads((a.HERE / 'three_round_condition_analysis.json').read_text())['samples'][0]
    else:
        sample = json.loads((a.HERE / 'cheap_capacity_linearization_sample.json').read_text())
    forms = [int(v, 16) for v in sample['cover_forms_hex']]
    polys = detail['X3_residual_polynomials_monomial_masks']
    n = len(detail['background_basis_hex'])
    _, kernel, rank = a.system(forms, [0] * len(forms), n)
    assert rank == len(forms)
    # The polar form is independent of the affine coset. These checks prove
    # that EVERY guess branch has affine, rather than quadratic, conditions.
    for p in polys:
        for i, v in enumerate(kernel):
            for w in kernel[:i]:
                assert evaluate(p, v ^ w) ^ evaluate(p, v) ^ evaluate(p, w) ^ evaluate(p, 0) == 0
    lifts = [a.system(forms, [int(i == j) for j in range(rank)], n)[0] for i in range(rank)]
    state = [0] + [int(v, 16) for v in detail['native_capacity_words_hex']]
    for active in a.U:
        changed = state.copy()
        changed[0] ^= 1 << (63 - active)
        delta = [x ^ y for x, y in zip(a.ref.rnd(state, 0x69), a.ref.rnd(changed, 0x69))]
        assert all((delta[i] >> (63 - j) & 1) == 0 for i, j in sample['cover_original_X2_positions'])
    original_part = int(detail['background_particular_hex'], 16)
    original_basis = [int(v, 16) for v in detail['background_basis_hex']]
    rng = random.Random(20261013)
    trials = 8192
    consistent = 0
    backgrounds = 0
    rank_hist = collections.Counter()
    output_rank_hist = collections.Counter()
    selected_base_hist = collections.Counter()
    witness = None
    for trial in range(trials):
        guess = rng.getrandbits(rank)
        base = substitute_background(guess, 0, lifts)
        rhs = [evaluate(p, base) for p in polys]
        rows = [sum((evaluate(p, base ^ v) ^ c) << k for k, v in enumerate(kernel)) for p, c in zip(polys, rhs)]
        coefficient_rank = a.ref.rank(rows)
        rank_hist[coefficient_rank] += 1
        solution, free, _ = a.system(rows, rhs, len(kernel))
        if solution is None:
            continue
        consistent += 1
        for assignment in range(1 << len(free)):
            v = substitute_background(assignment, solution, free)
            t = substitute_background(v, base, kernel)
            assert all(evaluate(p, t) == 0 for p in polys)
            backgrounds += 1
            g = substitute_background(t, original_part, original_basis)
            message = sum(((g >> k) & 1) << (63 - j) for k, j in enumerate(a.G))
            # Check every active assignment against the exact round function.
            outputs = []
            for active in range(4):
                native = message ^ sum(((active >> k) & 1) << (63 - j) for k, j in enumerate(a.U))
                x2 = a.ref.rnd([native] + state[1:], 0x69)
                x3 = a.ref.rnd(x2, 0x5a)
                assert all((x2[i] >> (63 - j) & 1) == value for i, j, value in a.COND['X2'])
                assert all((x3[i] >> (63 - j) & 1) == value for i, j, value in a.COND['X3'])
                outputs.append(a.ref.rnd(x3, 0x4b)[0])
            assert all((outputs[0] ^ outputs[1] ^ outputs[2] ^ outputs[3]) >> (63-j) & 1 == 0 for j in a.SELECT)
            out_rows = [((outputs[0] ^ outputs[1]) >> (63-j) & 1) | (((outputs[0] ^ outputs[2]) >> (63-j) & 1) << 1) for j in a.SELECT]
            output_rank_hist[a.ref.rank(out_rows)] += 1
            selected_base_hist[sum(((outputs[0] >> (63-j)) & 1) << k for k, j in enumerate(a.SELECT))] += 1
            if witness is None:
                witness = {'guessed_linear_parts_hex': hex(guess), 'linear_system_rows': rows,
                           'linear_system_rhs': rhs, 'solution_hex': hex(solution),
                           'native_background_message_hex': f'{message:016x}'}
    out = {'scope': 'one saved capacity; exact linearization, sampled guess branches',
           'background_dimension_before_guesses': n, 'independent_guess_count': rank,
           'background_dimension_after_guesses': len(kernel), 'condition_count': len(polys),
           'all_branches_affine': True, 'guesses_independent_of_active_variables': True,
           'sampled_branches': trials, 'consistent_branches': consistent,
           'generated_backgrounds': backgrounds, 'linear_system_rank_histogram': dict(rank_hist),
           'selected_output_rank_histogram': dict(output_rank_hist),
           'selected_output_at_zero_active_histogram': dict(selected_base_hist),
           'backgrounds_per_branch': backgrounds / trials,
           'guess_X2_positions': sample['cover_original_X2_positions'], 'witness': witness}
    name = 'quadratic_linearization_validation.json' if args.case == 'first' else 'cheap_quadratic_linearization_validation.json'
    (a.HERE / name).write_text(json.dumps(out, indent=2) + '\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
