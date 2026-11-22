#!/usr/bin/env python3
"""Check saved 1+1+n jobs, exact prefix controls, and serialized DATF DAGs."""
import argparse
import gzip
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from datf import OrdinaryDATF, evaluate, variable
from datf_fixed_tail import make_jobs, evaluate_job, uniform_balance_certificate
from datf_k06 import decode, integer_core, value_record
from fixed_tail_search import ExactPrefix, RC, dot


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path, default=ROOT / 'results/datf_one_one_tail')
    args = parser.parse_args()
    data = json.loads((args.results / 'results.json').read_text())
    source = json.loads(Path(data['source']['path']).read_text())
    jobs = make_jobs(source, data['top_per_tail'])
    assert len(jobs) == len(data['results'])
    assert sum(j['was_best_by_round'] and j['kind'] == 'full_endpoint' for j in jobs) == 5
    rng = random.Random(2026101607)
    # A translation may only move uniform coordinates. Check every returned
    # certificate against the complete truth table, independently of Bias().
    certificates = 0
    assert uniform_balance_certificate(variable(0),{0}) is None
    for _ in range(80):
        poly = frozenset(i for i in range(64) if i.bit_count()<=2 and rng.randrange(3)==0)
        biased = {i for i in range(6) if rng.randrange(3)==0}
        cert = uniform_balance_certificate(poly,biased)
        if cert:
            shift = sum(1<<v for v in cert['shift_variable_ids'])
            assert not any(shift>>v&1 for v in biased)
            assert all(evaluate(poly,x)^evaluate(poly,x^shift)==1 for x in range(64))
            certificates += 1
    exact_cache = {}
    prefix_checks = dag_checks = 0
    for job, result in zip(jobs, data['results']):
        assert {k: result[k] for k in job} == job
        assert job['start_round_index'] == 12-job['parent_total_rounds']
        assert result['remaining_constants'] == [f'0x{RC[job["start_round_index"]+r]:02x}'
                                                  for r in range(1,job['rounds'])]
        assert job['parent_total_rounds'] == 1+1+job['r_l']
        assert result['status'] in ('completed','resource_limit')
        if result['status'] == 'completed':
            assert result['factor_audit']['status'] in ('PASS','audit_resource_limit')
            if result['factor_audit']['status']=='PASS':
                assert result['factor_audit']['correlation'] == result['correlation']
            rho = Fraction(int(result['correlation']['numerator']), int(result['correlation']['denominator']))
            assert result['bias'] == value_record(rho/2)
        else:
            assert result['correlation'] is result['bias'] is None
        if job['kind'] == 'prefix_control':
            gamma = decode(job['output_mask'])
            if gamma not in exact_cache:
                exact_cache[gamma] = ExactPrefix(gamma)
            exact = exact_cache[gamma].evaluate(decode(job['input_difference']), RC[job['start_round_index']+1])
            exact = value_record(Fraction(exact['numerator'],exact['denominator']))
            assert exact == job['exact_prefix_correlation']
            if result['status'] == 'completed':
                assert result['matches_exact_prefix'] == (exact == result['correlation'])
            prefix_checks += 1
        with gzip.open(args.results/(job['id']+'.json.gz'),'rt') as stream:
            payload = json.load(stream)
        assert payload['result'] == result
        engine = OrdinaryDATF()
        for entry in payload['dag']:
            expression = entry['expression']
            assert entry['id'] == len(engine.definitions)
            if expression is not None:
                assert all(all(v < entry['id'] for v in t) for t in expression)
                expression = frozenset(sum(1<<v for v in term) for term in expression)
            engine.definitions.append(expression)
        target = frozenset(sum(1<<v for v in term) for term in payload['target_polynomial'])
        for _ in range(4):
            state = tuple(rng.getrandbits(64) for _ in range(5))
            assignment = engine.evaluate_dag([(state[i]>>j)&1 for i in range(5) for j in range(64)])
            left = integer_core(state, job['rounds'],job['start_round_index'])
            right = integer_core(tuple(x^y for x,y in zip(state,decode(job['input_difference']))),
                                 job['rounds'],job['start_round_index'])
            assert evaluate(target,assignment) == dot(tuple(x^y for x,y in zip(left,right)),decode(job['output_mask']))
            dag_checks += 1
    # Forced exhaustion must stay unknown, never masquerade as a zero estimate.
    result, _ = evaluate_job(jobs[0], max_bias_nodes=1, samples=1)
    assert result['status'] == 'resource_limit' and result['correlation'] is None
    print(json.dumps({'status':'PASS','jobs_checked':len(jobs),'exact_prefix_controls_recomputed':prefix_checks,
                      'serialized_DAG_parity_checks':dag_checks,'resource_limit_is_unknown':True,
                      'weighted_balance_certificates_truth_checked':certificates,
                      'factor_audits_completed':sum(r.get('factor_audit',{}).get('status')=='PASS' for r in data['results']),
                      'prefix_model_disagreements':sum(r.get('matches_exact_prefix') is False for r in data['results'])},indent=2))


if __name__ == '__main__':
    main()
