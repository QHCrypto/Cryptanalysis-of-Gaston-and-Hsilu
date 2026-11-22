#!/usr/bin/env python3
"""Finish limited audits by exact sums over disjoint variable components.

Recompute each affected endpoint with the original Shannon-based ordinary DATF.
This changes neither the estimator nor the saved correlation. Preserve the
original limited audit alongside the completed independent factor audit.
"""
import argparse
import gzip
import json
import time
from fractions import Fraction
from pathlib import Path

from datf import ZERO, add, components, indices, isolated, support
from datf_audit import weighted_quadratic_correlation
from datf_k06 import construct, decode, parity_polynomial, value_record
from datf_paper_parallel import atomic_json


def audit_components(engine, target, trace):
    polynomial, product, factors = target, Fraction(1), []
    for item in trace:
        lone = isolated(polynomial)
        rest = polynomial - frozenset(1 << v for v in indices(lone))
        marginals = {v: engine.marginal(v) for v in indices(support(rest))}
        checked, seen, pieces = Fraction(1), 0, []
        for part in components(rest):
            variables = support(part)
            assert not seen & variables
            seen |= variables
            rho, branches = weighted_quadratic_correlation(part, marginals)
            pieces.append({'polynomial': [list(indices(t)) for t in sorted(part)],
                           'correlation': value_record(rho), 'branches': branches})
            checked *= rho
        assert seen == support(rest)
        assert checked == Fraction(int(item['factor_numerator']), int(item['factor_denominator']))
        factors.append({'step': item['step'], 'correlation': value_record(checked),
                        'biased_marginals': {str(v): value_record(c) for v, c in marginals.items() if c},
                        'disjoint_components': pieces})
        product *= checked
        if not product:
            break
        original_inputs = [v for v in indices(lone) if engine.definitions[v] is None]
        if original_inputs:
            assert item['zero_reason'].startswith('isolated uniform input variable')
            factors[-1]['isolated_uniform_original_inputs'] = original_inputs
            product = Fraction(0)
            break
        polynomial = add(*(engine.definitions[v] for v in indices(lone))) if lone else ZERO
    return {'status': 'PASS', 'correlation': value_record(product), 'factors': factors,
            'method': 'exact weighted quadratic sum of disjoint variable components',
            'scope': 'ordinary-DATF independent factor arithmetic; not true permutation correlation'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    args = parser.parse_args()
    records = []
    for summary_path in sorted((args.results/'reassessment/jobs').glob('*.summary.json')):
        summary = json.loads(summary_path.read_text())
        if summary['factor_audit_status'] == 'PASS':
            continue
        assert summary['status'] == 'completed'
        started = time.perf_counter()
        path = Path(summary['artifact'])
        with gzip.open(path, 'rt') as f:
            payload = json.load(f)
        saved = payload['result']
        assert saved['factor_audit']['status'] == 'audit_resource_limit'
        engine, _, diff, _ = construct(decode(saved['input_difference']), saved['rounds'],
                                      saved['start_round_index'], max_bias_nodes=4000000)
        target = parity_polynomial(diff, decode(saved['output_mask']))
        trace = []
        rho = engine.estimate(target, trace=trace)
        assert value_record(rho) == saved['correlation']
        audit = audit_components(engine, target, trace)
        assert audit['correlation'] == saved['correlation']
        saved['original_limited_factor_audit'] = saved['factor_audit']
        saved['factor_audit'] = audit
        control = {'status': 'PASS', 'id': saved['id'], 'correlation': value_record(rho),
                   'backend': 'original ordinary DATF with Shannon Bias',
                   'engine_statistics': dict(engine.stats),
                   'elapsed_seconds': time.perf_counter()-started}
        saved['supplemental_original_backend_control'] = control
        temp = path.with_name(path.name+'.partial')
        with gzip.open(temp, 'wt') as f:
            json.dump(payload, f, separators=(',', ':'))
            f.write('\n')
        temp.replace(path)
        summary['factor_audit_status'] = 'PASS'
        atomic_json(summary_path, summary)
        records.append(control)
        print('AUDITED', saved['id'], round(control['elapsed_seconds'], 3), flush=True)
    receipt = args.results/'supplemental_factor_audits.json'
    if records:
        previous = json.loads(receipt.read_text()) if receipt.exists() else []
        atomic_json(receipt, previous + records)


if __name__ == '__main__':
    main()
