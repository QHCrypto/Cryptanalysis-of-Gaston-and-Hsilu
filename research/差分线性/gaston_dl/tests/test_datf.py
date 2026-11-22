#!/usr/bin/env python3
"""Validate ordinary DATF algebra, original bias loop, and Gaston integration."""
import itertools
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from datf import (ZERO, ONE, OrdinaryDATF, add, mul, variable, support,
                  evaluate, isolated, indices)
from datf_k06 import (construct, decode, symbolic_linear, symbolic_chi_pair,
                      parity_polynomial, verify_dag, make_jobs)
from datf_audit import (quadratic_correlation, weighted_quadratic_correlation,
                        balance_certificate, audit_job)
from linear_seed_scan import linear, embed, chi


def brute(poly, weights):
    answer = Fraction(0)
    variables = sorted(weights)
    for values in itertools.product((0, 1), repeat=len(variables)):
        assignment, probability = 0, Fraction(1)
        for v, x in zip(variables, values):
            assignment |= x << v
            probability *= (1 + (-1) ** x * weights[v]) / 2
        answer += (-1) ** evaluate(poly, assignment) * probability
    return answer


def main():
    rng = random.Random(2026101601)
    # Exact GF(2) multilinear arithmetic on independent Boolean values.
    for _ in range(40):
        a, b = (frozenset(i for i in range(64) if rng.randrange(8) == 0) for _ in range(2))
        for x in range(64):
            assert evaluate(add(a, b), x) == evaluate(a, x) ^ evaluate(b, x)
            assert evaluate(mul(a, b), x) == evaluate(a, x) & evaluate(b, x)
        assert mul(a, a) == a
    # Weighted independent Bias, checked against full product-distribution sums.
    for _ in range(40):
        model = OrdinaryDATF()
        for i in range(6):
            model.fresh()
            model.marginals[i] = rng.choice([Fraction(0), Fraction(1, 2), Fraction(-1, 2), Fraction(1), Fraction(-1)])
        poly = frozenset(i for i in range(64) if rng.randrange(8) == 0)
        assert model.independent_correlation(poly) == brute(poly, model.marginals)
    # Compare the independent quadratic elimination and its balanced-shift
    # certificate with exhaustive enumeration, independent of DATF's Bias().
    for _ in range(80):
        poly = frozenset(i for i in range(256) if i.bit_count() <= 2 and rng.randrange(3) == 0)
        actual = quadratic_correlation(poly)
        assert actual == brute(poly, {i: Fraction(0) for i in range(8)})
        if not actual:
            certificate = balance_certificate(poly)
            assert certificate
            h = sum(1 << v for v in certificate['shift_variable_ids'])
            assert all(evaluate(poly, x) ^ evaluate(poly, x ^ h) == 1 for x in range(256))
    for _ in range(40):
        poly = frozenset(i for i in range(64) if i.bit_count() <= 2 and rng.randrange(3) == 0)
        marginals = {i: rng.choice([Fraction(0), Fraction(1, 2), Fraction(-1, 2)])
                     for i in range(6)}
        actual, _ = weighted_quadratic_correlation(poly, marginals)
        assert actual == brute(poly, marginals)
    # Reproduce original Algorithm 1 in the published comparison Example 1.
    # Ordinary estimate = 0; exact expanded rho = 1/4 (bias = 1/8).
    model = OrdinaryDATF()
    x = [model.fresh() for _ in range(4)]
    y0 = model.fresh(add(mul(x[0], x[1]), x[2]), 'beta')
    y1 = model.fresh(add(mul(x[1], x[2]), y0), 'beta')
    z = add(mul(x[2], x[3]), y1)
    assert model.estimate(z) == 0
    exact = add(mul(x[2], x[3]), mul(x[1], x[2]), mul(x[0], x[1]), x[2])
    assert brute(exact, {i: Fraction(0) for i in range(4)}) == Fraction(1, 4)
    # All 32 x 32 single-chi difference/mask combinations, direct 32-value
    # truth-table check. These are exact, including signs and rho=2*epsilon.
    for delta in range(32):
        model = OrdinaryDATF()
        base = [[model.fresh() for _ in range(64)] for _ in range(5)]
        difference = [[ONE if j == 0 and (delta >> i) & 1 else ZERO for j in range(64)] for i in range(5)]
        _, diff = symbolic_chi_pair(base, difference)
        for gamma in range(32):
            poly = add(*(diff[i][0] for i in range(5) if gamma >> i & 1))
            expected = Fraction(sum((-1) ** ((gamma & (chi(v) ^ chi(v ^ delta))).bit_count() & 1)
                                    for v in range(32)), 32)
            assert model.estimate(poly) == expected
    # Every forward-L basis vector, independent integer and polynomial paths.
    for i, j in itertools.product(range(5), range(64)):
        state = [[ZERO] * 64 for _ in range(5)]
        state[i][j] = ONE
        actual = symbolic_linear(state)
        value = tuple(sum((actual[k][b] == ONE) << b for b in range(64)) for k in range(5))
        assert value == linear(embed(j, 1 << i))
    # A previously proved nonzero two-chi control: rho=1/2.
    source = json.loads((ROOT / 'results/pq_kernel_six_all/search.json').read_text())
    model, base, diff, _ = construct(embed(1, 16), 2, 7)
    assert model.estimate(parity_polynomial(diff, decode(source['groups']['join_2']['gamma']))) == Fraction(1, 2)
    # Exact DAG evaluation versus the actual permutation, all four requested
    # jobs and their actual parent-round constants, before any approximation.
    saved = json.loads((ROOT / 'results/datf_k06/results.json').read_text())
    assert [{k: r[k] for k in j} for j, r in zip(make_jobs(source), saved['results'])] == make_jobs(source)
    for job in saved['results']:
        model, base, diff, _ = construct(decode(job['input_difference']), job['rounds'], job['start_round_index'])
        verify_dag(model, base, diff, decode(job['input_difference']), job['rounds'], job['start_round_index'], samples=64)
        audit_job(job)
    result = {'status': 'PASS', 'single_chi_cases': 1024, 'forward_L_basis_checks': 320,
              'weighted_bias_full_enumerations': 40, 'quadratic_full_enumerations': 80,
              'weighted_quadratic_full_enumerations': 40,
              'original_algorithm_example': {'ordinary_rho': '0', 'exact_expanded_rho': '1/4'},
              'nonzero_two_chi_control_rho': '1/2', 'dag_full_state_checks': 256,
              'all_four_jobs_independently_audited': True}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
