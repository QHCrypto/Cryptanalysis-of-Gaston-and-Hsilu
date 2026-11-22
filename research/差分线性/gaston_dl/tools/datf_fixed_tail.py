#!/usr/bin/env python3
"""Ordinary DATF for the saved 1+1+n Gaston endpoints and two-round controls.

Uses the unchanged ordinary DATF engine and Gaston boundaries from datf_k06.
Every tail's saved candidates are ordered by the previous splice score; the
default evaluates the first candidate of each nonempty tail, not a new search.
"""
import argparse
import gzip
import hashlib
import json
import signal
import time
from fractions import Fraction
from pathlib import Path

from datf import BudgetExceeded, ZERO, add, evaluate, indices, isolated, support
from datf_k06 import construct, decode, parity_polynomial, value_record, verify_dag
from datf_audit import weighted_quadratic_correlation
from fixed_tail_search import RC


def make_jobs(source, top_per_tail=1):
    winners = {r['tail_id'] for r in source['best_by_R']}
    reports = sorted(source['reports'], key=lambda r: (r['R'], r['id'] not in winners, r['id']))
    jobs = []
    for tail in reports:
        assert tail['r_d'] == tail['r_m'] == 1
        assert tail['R'] == tail['r_l'] + 2 and tail['first_L_omitted']
        for rank, candidate in enumerate(tail['top_candidates'][:top_per_tail], 1):
            prefix_rho = Fraction(candidate['prefix']['numerator'], candidate['prefix']['denominator'])
            score = prefix_rho / (1 << tail['tail_w_squared'])
            for prefix in (True, False):
                jobs.append({'id': f'{tail["id"]}_c{rank:02d}_' + ('prefix2' if prefix else f'full_R{tail["R"]}'),
                             'tail_id': tail['id'], 'candidate_rank': rank,
                             'kind': 'prefix_control' if prefix else 'full_endpoint',
                             'rounds': 2 if prefix else tail['R'], 'parent_total_rounds': tail['R'],
                             'r_d': 1, 'r_m': 1, 'r_l': tail['r_l'],
                             'start_round_index': tail['start_round_index'],
                             'was_best_by_round': tail['id'] in winners and rank == 1,
                             'input_difference': candidate['input_difference'],
                             'output_mask': tail['join_mask'] if prefix else tail['output_mask'],
                             'join_mask': tail['join_mask'],
                             'tail_squared_weight': tail['tail_w_squared'],
                             'exact_prefix_correlation': value_record(prefix_rho),
                             'previous_splice_score': value_record(score)})
    return jobs


def uniform_balance_certificate(poly, biased_variables):
    """Prove zero by a translation that fixes all biased coordinates.

    This is an independent audit shortcut, not a DATF estimation rule.
    """
    variables = list(indices(support(poly)))
    rows = {v: 0 for v in variables}
    for term in poly:
        assert term.bit_count() <= 2
        if term.bit_count() == 2:
            u, v = list(indices(term))
            rows[u] ^= 1 << v
            rows[v] ^= 1 << u
    constraints = list(rows.values()) + [1 << v for v in biased_variables]
    basis = {}
    for row in constraints:
        while row:
            pivot = row.bit_length()-1
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
        assert all((row & shift).bit_count()%2 == 0 for row in constraints)
        if evaluate(poly,shift) ^ (0 in poly):
            return {'shift_variable_ids': list(indices(shift)),
                    'fixed_biased_variables': sorted(biased_variables),
                    'identity': 'f(z XOR h)=f(z) XOR 1, with h zero on every biased coordinate'}
    return None


def audit_trace(engine, target, trace):
    """Audit exact factor arithmetic only; independence remains an assumption."""
    polynomial, product, factors = target, Fraction(1), []
    for item in trace:
        lone = isolated(polynomial)
        rest = polynomial - frozenset(1 << v for v in indices(lone))
        marginals = {v: engine.marginal(v) for v in indices(support(rest))}
        biased = {v for v,c in marginals.items() if c}
        certificate = uniform_balance_certificate(rest,biased) if int(item['factor_numerator']) == 0 else None
        if certificate:
            checked, branches = Fraction(0), []
        else:
            if len(biased) > 12:
                return {'status':'audit_resource_limit','factors':factors,
                        'reason':'no uniform balancing shift; exhaustive audit exceeds 12 biased coordinates'}
            checked, branches = weighted_quadratic_correlation(rest, marginals)
        assert checked == Fraction(int(item['factor_numerator']), int(item['factor_denominator']))
        factors.append({'step': item['step'], 'correlation': value_record(checked),
                        'variables': len(marginals), 'terms': len(rest),
                        'biased_marginals': {str(v): value_record(c) for v, c in marginals.items() if c},
                        'branches': branches, 'uniform_balance_certificate':certificate})
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
            'scope': 'independent verification of ordinary-DATF factor arithmetic, not true permutation correlation'}


def evaluate_job(job, max_bias_nodes=1000000, samples=32, progress=False):
    began = time.perf_counter()
    delta = decode(job['input_difference'])
    engine, base, diff, stages = construct(delta, job['rounds'], job['start_round_index'], max_bias_nodes)
    verify_dag(engine, base, diff, delta, job['rounds'], job['start_round_index'], samples=samples)
    target = parity_polynomial(diff, decode(job['output_mask']))
    trace = []
    old_handler = None
    if progress:
        def heartbeat(signum, frame):
            print('PROGRESS',job['id'],'Bias nodes',engine.stats.get('bias_nodes',0),
                  'seconds',round(time.perf_counter()-began,1),flush=True)
        old_handler = signal.signal(signal.SIGALRM, heartbeat)
        signal.setitimer(signal.ITIMER_REAL,25,25)
    result = {**job, 'dag_verified_samples': samples, 'stages': stages,
              'target_terms': len(target), 'target_variables': support(target).bit_count(),
              'remaining_constants': [f'0x{RC[job["start_round_index"]+r]:02x}' for r in range(1, job['rounds'])]}
    try:
        rho = engine.estimate(target, label=job['id'], trace=trace)
        result.update(status='completed', correlation=value_record(rho), bias=value_record(rho/2))
    except BudgetExceeded as exc:
        result.update(status='resource_limit', correlation=None, bias=None, reason=str(exc))
    finally:
        if progress:
            signal.setitimer(signal.ITIMER_REAL,0)
            signal.signal(signal.SIGALRM,old_handler)
    result.update(trace=trace, engine_statistics=dict(engine.stats),
                  estimation_seconds=time.perf_counter()-began)
    if result['status'] == 'completed':
        if progress:
            print('ESTIMATED',job['id'],result['correlation'],'begin factor audit',flush=True)
        result['factor_audit'] = audit_trace(engine, target, trace)
        if result['factor_audit']['status'] == 'PASS':
            assert result['factor_audit']['correlation'] == result['correlation']
        if job['kind'] == 'prefix_control':
            # A disagreement would be a model discrepancy to report, not a
            # reason to overwrite either the exact value or the DATF estimate.
            result['matches_exact_prefix'] = result['correlation'] == job['exact_prefix_correlation']
    result['elapsed_seconds_including_audit'] = time.perf_counter()-began
    payload = {'result': result, 'target_polynomial': [list(indices(t)) for t in sorted(target)],
               'dag': engine.serialized_dag()}
    return result, payload


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--top-per-tail', type=int, default=1)
    parser.add_argument('--max-bias-nodes', type=int, default=1000000)
    parser.add_argument('--resume',action='store_true',help='continue an interrupted directory with identical inputs and limits')
    args = parser.parse_args()
    if args.top_per_tail < 1 or args.max_bias_nodes < 1:
        parser.error('limits must be positive')
    source = json.loads(args.source.read_text())
    jobs = make_jobs(source, args.top_per_tail)
    args.output.mkdir(parents=True, exist_ok=args.resume)
    document = {'method': 'ordinary DATF / original isolated-variable Algorithm 1',
                'enhanced_DATF_features': False, 'first_L_omitted': True,
                'input_distribution': 'all 320 independent uniform input bits; no path conditions',
                'transitional_cut': 'after L and iota, immediately before each chi after the first',
                'top_per_tail': args.top_per_tail,
                'candidate_scope': 'saved candidates from fixed_paper_tail; no new endpoint search',
                'assumptions': ['independent variable marginals in Bias',
                                'independence of factors at different isolated-variable substitution steps'],
                'resource_limit': {'max_bias_nodes_per_job': args.max_bias_nodes, 'uncomputed_is_not_zero': True},
                'source': {'path': str(args.source.resolve()), 'sha256': hashlib.sha256(args.source.read_bytes()).hexdigest()},
                'skipped_empty_tails': [r['id'] for r in source['reports'] if not r['top_candidates']],
                'results': []}
    if args.resume:
        previous = json.loads((args.output/'results.json').read_text())
        assert {k:v for k,v in previous.items() if k!='results'} == {k:v for k,v in document.items() if k!='results'}
        for job,result in zip(jobs,previous['results']):
            assert {k:result[k] for k in job} == job
            assert (args.output/(job['id']+'.json.gz')).exists()
        document['results'] = previous['results']
    for job in jobs[len(document['results']):]:
        print('START', job['id'], flush=True)
        result, payload = evaluate_job(job, args.max_bias_nodes, progress=True)
        # Exact variable definitions remain available without writing large
        # uncompressed JSON files for every long-tail endpoint.
        with gzip.open(args.output / (job['id'] + '.json.gz'), 'xt', encoding='utf-8') as stream:
            json.dump(payload, stream, separators=(',', ':'))
            stream.write('\n')
        document['results'].append(result)
        (args.output / 'results.json').write_text(json.dumps(document, ensure_ascii=False, indent=2)+'\n')
        print('DONE', job['id'], result['status'], result['correlation'],
              'nodes',result['engine_statistics'].get('bias_nodes',0), flush=True)
    print(args.output.resolve())


if __name__ == '__main__':
    main()
