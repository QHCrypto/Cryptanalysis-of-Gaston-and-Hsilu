#!/usr/bin/env python3
"""Independent algebraic audit of ordinary DATF's model-zero factors.

This does not change DATF's transitions or estimates. Quadratic elimination
checks the same independence factor by a second exact method. Biased variables
are conditioned explicitly; the remaining quadratic sum is uniform.
"""
import argparse
import itertools
import json
from fractions import Fraction
from pathlib import Path

from datf import ONE, add, mul, indices, support, isolated, evaluate, condition
from datf_k06 import construct, decode, parity_polynomial, value_record


def quadratic_correlation(p):
    assert all(t.bit_count() <= 2 for t in p)
    factor = Fraction(1)
    while True:
        edge = next((t for t in sorted(p) if t.bit_count() == 2), None)
        if edge is None:
            return Fraction(0) if support(p) else factor * (-1 if 0 in p else 1)
        u, v = list(indices(edge))
        ub, vb = 1 << u, 1 << v
        a = frozenset(t ^ ub for t in p if t & ub and not t & vb)
        b = frozenset(t ^ vb for t in p if t & vb and not t & ub)
        c = frozenset(t for t in p if not t & edge)
        # Sum over u,v: E(-1)^[uv+u*a+v*b+c] = (1/2)(-1)^[c+a*b].
        p = add(c, mul(a, b))
        factor /= 2


def rank(vectors):
    basis = {}
    for v in vectors:
        while v:
            pivot = v.bit_length() - 1
            if pivot in basis:
                v ^= basis[pivot]
            else:
                basis[pivot] = v
                break
    return len(basis)


def balance_certificate(p):
    """Find h with f(z+h)=f(z)+1, by the kernel of the polar matrix."""
    variables = list(indices(support(p)))
    rows = {v: 0 for v in variables}
    for term in p:
        if term.bit_count() == 2:
            u, v = list(indices(term))
            rows[u] ^= 1 << v
            rows[v] ^= 1 << u
    basis = {}
    for row in rows.values():
        while row:
            pivot = row.bit_length() - 1
            if pivot in basis:
                row ^= basis[pivot]
            else:
                basis[pivot] = row
                break
    for free in variables:
        if free in basis:
            continue
        shift = 1 << free
        for pivot, row in sorted(basis.items()):
            if (row & shift).bit_count() & 1:
                shift |= 1 << pivot
        assert all(not ((row & shift).bit_count() & 1) for row in rows.values())
        if evaluate(p, shift) ^ (0 in p):
            return {'variables': len(variables), 'polar_rank': len(basis),
                    'shift_variable_ids': list(indices(shift)),
                    'identity': 'f(z XOR h) = f(z) XOR 1 in the independent-variable model'}
    return None


def weighted_quadratic_correlation(p, marginals):
    """Independent audit: condition biased bits, then eliminate quadratic pairs."""
    biased = sorted(v for v in indices(support(p)) if marginals[v])
    # These jobs have at most two biased bits in a global factor. Never silently
    # turn an impractical audit into a zero or a uniform approximation.
    if len(biased) > 16:
        raise RuntimeError('independent audit has more than 16 biased variables')
    total, branches = Fraction(0), []
    for values in itertools.product((0, 1), repeat=len(biased)):
        reduced, weight = p, Fraction(1)
        for v, bit in zip(biased, values):
            weight *= (1 + (-1) ** bit * marginals[v]) / 2
            reduced = condition(reduced, v, bit)
        if not weight:
            continue
        rho = quadratic_correlation(reduced)
        certificate = balance_certificate(reduced) if not rho else None
        assert rho or certificate
        total += weight * rho
        branches.append({'assignment': {str(v): bit for v, bit in zip(biased, values)},
                         'weight': value_record(weight), 'uniform_correlation': value_record(rho),
                         'balance_certificate': certificate})
    return total, branches


def audit_job(job):
    engine, _, d, _ = construct(decode(job['input_difference']), job['rounds'], job['start_round_index'])
    e = parity_polynomial(d, decode(job['output_mask']))
    trace = []
    rho = engine.estimate(e, label=job['id'], trace=trace)
    assert value_record(rho) == job['correlation']
    records = []
    for item in trace:
        lone = isolated(e)
        factor_poly = e - frozenset(1 << v for v in indices(lone))
        probabilities = {v: engine.marginal(v) for v in indices(support(factor_poly))}
        checked, branches = weighted_quadratic_correlation(factor_poly, probabilities)
        assert checked == Fraction(int(item['factor_numerator']), int(item['factor_denominator']))
        record = {'step': item['step'], 'variables': len(probabilities), 'terms': len(factor_poly),
                  'independent_quadratic_correlation': value_record(checked),
                  'all_variable_marginals_uniform_in_model': all(c == 0 for c in probabilities.values()),
                  'biased_variable_marginals': {str(v): value_record(c) for v, c in probabilities.items() if c},
                  'conditioned_branches': branches}
        # First-cut beta variables are affine forms in the same 12 original
        # input bits. Count explicitly how much independence this cut discards.
        first_betas = [v for v in probabilities if engine.metadata[v]['kind'] == 'beta'
                       and engine.metadata[v]['stage'] == 1]
        direct_inputs = [v for v in probabilities if engine.metadata[v]['kind'] == 'input']
        if first_betas:
            vectors = [add(engine.definitions[v], ONE) if 0 in engine.definitions[v]
                       else engine.definitions[v] for v in first_betas]
            assert all(all(t.bit_count() == 1 for t in p) for p in vectors)
            linear_rows = [support(p) for p in vectors] + [1 << v for v in direct_inputs]
            record['first_cut_dependency'] = {
                'beta_variables': len(first_betas), 'direct_original_variables': len(direct_inputs),
                'modeled_as_independent_coordinates': len(linear_rows),
                'actual_affine_rank': rank(linear_rows),
                'affine_relations_lost_by_independence': len(linear_rows) - rank(linear_rows)}
        records.append(record)
        e = add(*(engine.definitions[v] for v in indices(lone)))
    return {'id': job['id'], 'factor_audits': records,
            'conclusion': 'The reported zeros are exact for ordinary DATF independence factors; no claim of zero true Gaston correlation.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, required=True)
    args = parser.parse_args()
    data = json.loads((args.results / 'results.json').read_text())
    audits = [audit_job(j) for j in data['results'] if j['status'] == 'completed']
    (args.results / 'zero_audit.json').write_text(json.dumps(audits, ensure_ascii=False, indent=2) + '\n')
    for a in audits:
        print(a['id'], 'independently checked factors',
              [f['independent_quadratic_correlation'] for f in a['factor_audits']], flush=True)


if __name__ == '__main__':
    main()
