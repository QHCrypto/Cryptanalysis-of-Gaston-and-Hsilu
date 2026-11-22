#!/usr/bin/env python3
"""Evaluate the complete K06 prefix/endpoints with ordinary DATF.

No B1/B2 conditions are imposed: all original 320 input bits are uniform.
Replacement cuts are immediately before chi, after L and iota. The first L
is omitted and its constant absorbed into the uniform starting state.
"""
import argparse
import hashlib
import json
import math
import random
import time
from fractions import Fraction
from pathlib import Path

from datf import (ZERO, ONE, OrdinaryDATF, BudgetExceeded, add, mul, evaluate,
                  indices, support)
from linear_seed_scan import EAST, WEST, TWIN, rol, linear, column
from fixed_tail_search import RC, state_chi, dot


def decode(state):
    return tuple(int(w, 16) for w in state['words'])


def rotate_row(row, shift):
    return [row[(j - shift) % 64] for j in range(64)]


def symbolic_linear(state):
    east = [rotate_row(row, e) for row, e in zip(state, EAST)]
    p = [add(*(row[j] for row in east)) for j in range(64)]
    q = [add(*(east[i][(j - TWIN[i]) % 64] for i in range(5))) for j in range(64)]
    effect = [add(p[(j - 23) % 64], p[(j - 24) % 64],
                  q[(j - 23) % 64], q[(j - 41) % 64]) for j in range(64)]
    return [rotate_row([add(row[j], effect[j]) for j in range(64)], w)
            for row, w in zip(east, WEST)]


def symbolic_chi_pair(a, d):
    base, diff = [], []
    for i in range(5):
        k, l = (i + 1) % 5, (i + 2) % 5
        base.append([add(a[i][j], a[l][j], mul(a[k][j], a[l][j])) for j in range(64)])
        diff.append([add(d[i][j], d[l][j], mul(a[k][j], d[l][j]),
                         mul(d[k][j], a[l][j]), mul(d[k][j], d[l][j])) for j in range(64)])
    return base, diff


def construct(delta, rounds, start, max_bias_nodes=1000000, engine_class=OrdinaryDATF):
    engine = engine_class(max_bias_nodes)
    a = [[engine.fresh(kind='input', row=i, column=j) for j in range(64)] for i in range(5)]
    d = [[ONE if delta[i] >> j & 1 else ZERO for j in range(64)] for i in range(5)]
    stages = []
    for r in range(rounds):
        if r:
            a, d = symbolic_linear(a), symbolic_linear(d)
            for j in range(64):
                if RC[start + r] >> j & 1:
                    a[0][j] = add(a[0][j], ONE)
            for i in range(5):
                for j in range(64):
                    a[i][j] = engine.transition(a[i][j], 'alpha', r, i, j)
                    d[i][j] = engine.transition(d[i][j], 'beta', r, i, j)
        a, d = symbolic_chi_pair(a, d)
        stages.append({'chi_layers': r + 1, 'dag_variables': len(engine.definitions),
                       'alpha_variables': sum(m['kind'] == 'alpha' for m in engine.metadata),
                       'beta_variables': sum(m['kind'] == 'beta' for m in engine.metadata),
                       'difference_terms': sum(len(p) for row in d for p in row),
                       'max_bit_difference_terms': max(len(p) for row in d for p in row)})
    return engine, a, d, stages


def parity_polynomial(difference, mask):
    return add(*(difference[i][j] for i in range(5) for j in range(64) if mask[i] >> j & 1))


def integer_core(state, rounds, start):
    state = state_chi(state)
    for r in range(1, rounds):
        state = linear(state)
        state = (state[0] ^ RC[start + r], *state[1:])
        state = state_chi(state)
    return state


def verify_dag(engine, base, diff, delta, rounds, start, samples=32):
    rng = random.Random(20261016 + rounds * 100 + start)
    for _ in range(samples):
        state = tuple(rng.getrandbits(64) for _ in range(5))
        input_bits = [(state[i] >> j) & 1 for i in range(5) for j in range(64)]
        assignment = engine.evaluate_dag(input_bits)
        actual_base = tuple(sum(evaluate(base[i][j], assignment) << j for j in range(64)) for i in range(5))
        actual_diff = tuple(sum(evaluate(diff[i][j], assignment) << j for j in range(64)) for i in range(5))
        expected_base = integer_core(state, rounds, start)
        expected_pair = integer_core(tuple(x ^ y for x, y in zip(state, delta)), rounds, start)
        assert actual_base == expected_base
        assert actual_diff == tuple(x ^ y for x, y in zip(expected_base, expected_pair))


def value_record(value):
    return {'numerator': str(value.numerator), 'denominator': str(value.denominator),
            'log2_abs': math.log2(abs(value.numerator)) - math.log2(value.denominator) if value else None}


def make_jobs(source):
    jobs = []
    for tail_id in ('D3_chi2_to_3', 'D3_chi2_to_4'):
        tail = next(r for r in source['reports'] if r['id'] == tail_id)
        assert tail['best']['orbit'] == 'K06' and tail['best']['rotation'] == 22
        for prefix in (True, False):
            jobs.append({'id': f'K06_prefix3_for_R{tail["R"]}' if prefix else f'K06_full_R{tail["R"]}',
                         'rounds': 3 if prefix else tail['R'], 'start_round_index': tail['start_round_index'],
                         'parent_total_rounds': tail['R'], 'tail_id': tail_id,
                         'input_difference': tail['best']['d0'],
                         'output_mask': source['groups'][tail['connection_group']]['gamma'] if prefix else tail['output_mask'],
                         'previous_restricted_score': tail['best']['restricted_signed_prefix_sum'] if prefix else tail['best']['restricted_splice_score']})
    return sorted(jobs, key=lambda j: (j['rounds'], j['parent_total_rounds']))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--max-bias-nodes', type=int, default=1000000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source = json.loads(args.source.read_text())
    jobs = make_jobs(source)
    results = []
    for job in jobs:
        start_time = time.perf_counter()
        engine, base, diff, stages = construct(decode(job['input_difference']), job['rounds'],
                                              job['start_round_index'], args.max_bias_nodes)
        verify_dag(engine, base, diff, decode(job['input_difference']), job['rounds'], job['start_round_index'])
        target = parity_polynomial(diff, decode(job['output_mask']))
        trace = []
        result = {**job, 'target_terms': len(target), 'target_variables': support(target).bit_count(),
                  'stages': stages, 'dag_verified_samples': 32,
                  'round_constants_after_first_chi': [hex(RC[job['start_round_index'] + r]) for r in range(1, job['rounds'])]}
        try:
            correlation = engine.estimate(target, label=job['id'], trace=trace)
            result.update({'status': 'completed', 'correlation': value_record(correlation),
                           'bias': value_record(correlation / 2),
                           'interpretation': 'ordinary-DATF independence-model estimate; not a proof about the true correlation'})
        except BudgetExceeded as e:
            result.update({'status': 'resource_limit', 'reason': str(e), 'correlation': None, 'bias': None})
        result.update({'trace': trace, 'engine_statistics': dict(engine.stats),
                       'zero_witnesses': engine.zero_witnesses, 'elapsed_seconds': time.perf_counter() - start_time})
        with (args.output / (job['id'] + '.json')).open('x') as f:
            json.dump({'result': result, 'target_polynomial': [list(indices(t)) for t in sorted(target)],
                       'dag': engine.serialized_dag()}, f, indent=2)
            f.write('\n')
        results.append(result)
        print(job['id'], result['status'], result['correlation'], 'nodes', engine.stats['bias_nodes'], flush=True)
        (args.output / 'results.json').write_text(json.dumps({
            'method': 'ordinary DATF / original isolated-variable Algorithm 1',
            'first_L_omitted': True, 'input_distribution': 'all 320 independent uniform input bits; no path conditions',
            'transitional_cut': 'immediately before each chi after the first: L and iota included before substitution',
            'transitional_rules': 'replace expressions with more than one variable; reuse identical/complementary beta definitions',
            'assumptions': ['Current alpha variables are treated as independent uniform variables when Bias is evaluated.',
                            'Nonisolated variables are evaluated under independent marginal distributions.',
                            'Factors peeled at different isolated-variable substitution steps are multiplied as independent.'],
            'enhanced_DATF_features': False,
            'resource_limit': {'max_bias_nodes_per_job': args.max_bias_nodes, 'uncomputed_is_not_zero': True},
            'source': {'path': str(args.source.resolve()), 'sha256': hashlib.sha256(args.source.read_bytes()).hexdigest()},
            'results': results}, ensure_ascii=False, indent=2) + '\n')
    print(args.output.resolve())


if __name__ == '__main__':
    main()
