#!/usr/bin/env python3
"""Re-evaluate eight representative pool winners at their full endpoints.

The splice score and the full ordinary-DATF estimate are different quantities.
This small diagnostic does not re-rank every candidate in the matching pools.
"""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

from datf_fixed_tail import evaluate_job


SELECTION = [
    ('D1', 1, 'L1_pool'), ('D1', 2, 'L1_pool'), ('D2', 1, 'L1_pool'),
    ('D1', 1, 'D3_chi2_to_3'), ('D1', 2, 'D3_chi3_to_4'),
    ('D2', 1, 'D3_chi2_to_3'), ('D2', 2, 'D3_chi3_to_4'),
    ('D1', 2, 'D2_repeat_2'),
]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--max-bias-nodes', type=int, default=1000000)
    p.add_argument('--resume', action='store_true')
    args = p.parse_args()
    source = json.loads(args.source.read_text())
    lookup = {(s['family'], s['middle_rounds'], s['tail']): s for s in source['summaries']}
    jobs = []
    for family, middle, tail in SELECTION:
        b = lookup[family, middle, tail]['best']
        if b is None:
            continue
        jobs.append({'id': f'{family}_M{middle}_{tail}', 'kind': 'full_endpoint',
                     'rounds': b['R'], 'start_round_index': b['start_round_index'],
                     'input_difference': {'words': b['input_difference']},
                     'output_mask': {'words': b['output_mask']},
                     'matching_winner': b})
    document = {'method': 'ordinary DATF Algorithm 1 on full endpoints',
                'first_L_omitted': True,
                'selection': 'First ranked match in eight prespecified nonempty cells; not every tie or every matched endpoint.',
                'source_sha256': hashlib.sha256(args.source.read_bytes()).hexdigest(),
                'max_bias_nodes': args.max_bias_nodes, 'results': []}
    args.output.mkdir(parents=True, exist_ok=args.resume)
    if args.resume:
        previous = json.loads((args.output/'results.json').read_text())
        assert {k:v for k,v in previous.items() if k != 'results'} == {k:v for k,v in document.items() if k != 'results'}
        for job, result in zip(jobs, previous['results']):
            assert {k:result[k] for k in job} == job
            assert (args.output/(job['id']+'.json.gz')).exists()
        document['results'] = previous['results']
    for job in jobs[len(document['results']):]:
        print('START', job['id'], flush=True)
        result, payload = evaluate_job(job, args.max_bias_nodes, samples=32, progress=True)
        with gzip.open(args.output/(job['id']+'.json.gz'), 'xt', encoding='utf-8') as f:
            json.dump(payload, f, separators=(',', ':'))
            f.write('\n')
        document['results'].append(result)
        (args.output/'results.json').write_text(json.dumps(document, ensure_ascii=False, indent=2)+'\n')
        print('DONE', job['id'], result['status'], result['correlation'], flush=True)


if __name__ == '__main__':
    main()
