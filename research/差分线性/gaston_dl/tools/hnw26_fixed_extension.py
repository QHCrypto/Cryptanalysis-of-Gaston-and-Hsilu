#!/usr/bin/env python3
"""Evaluate the two eight-bit fixed-tail extensions with HNW26, with sharding."""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
import math
import subprocess
import time
from collections import Counter
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def save(path, data):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n')
    tmp.replace(path)


def frac(x):
    return Fraction(int(x['numerator']), int(x['denominator']))


def number(x):
    return {'numerator': str(x.numerator), 'denominator': str(x.denominator),
            'log2_abs': math.log2(abs(x.numerator)) - math.log2(x.denominator) if x else None}


def chi(x):
    return sum((((x >> i) & 1) ^ ((1 ^ ((x >> ((i + 1) % 5)) & 1)) &
                                ((x >> ((i + 2) % 5)) & 1))) << i for i in range(5))


DDT = [[sum(chi(x) ^ chi(x ^ a) == b for x in range(32)) for b in range(32)] for a in range(32)]


def column(words, j):
    return sum(((w >> j) & 1) << i for i, w in enumerate(words))


def validate(data, row, begin, end):
    candidate = row['candidate']
    delta = [int(w, 16) for w in candidate['input_difference']['words']]
    assert data['status'] == 'completed'
    assert data['branch_begin'] == begin and data['branch_end'] == end
    assert data['completed_branches'] == end - begin
    assert data['rounds'] == row['rounds'] and data['start_round_index'] == row['start_round_index']
    assert [int(w, 16) for w in data['input_words'].split(';')] == delta
    assert [int(w, 16) for w in data['output_words'].split(';')] == [int(w, 16) for w in row['output_mask']['words']]
    assert [int(w, 16) for w in data['prefix_mask_words'].split(';')] == [int(w, 16) for w in candidate['join_mask']['words']]
    assert data['expected_first_DDT_branches'] == 65536
    if data['mode'] == 'full':
        assert data['first_linear_killed_value_masks'] == 1984
    assert data['certified_fast_second_round']
    support = {j: [b for b in range(32) if DDT[column(delta, j)][b]] for j in range(64) if column(delta, j)}
    mass = total = prefix_total = Fraction(0)
    nonzero = prefix_nonzero = 0
    extinction = Counter()
    prefix_branches = []
    for index, br in enumerate(data['branches'], begin):
        assert index == br['branch_index']
        words = [int(w, 16) for w in br['first_output_difference'].split(';')]
        # Independently reconstruct both the DDT probability and the mixed
        # radix ordinal, proving there are no omissions or duplicated branches.
        probability = Fraction(1)
        ordinal = 0
        for j in range(64):
            a, b = column(delta, j), column(words, j)
            count = DDT[a][b]
            assert count > 0
            if j in support:
                ordinal = ordinal * len(support[j]) + support[j].index(b)
                probability *= Fraction(count, 32)
            else:
                assert b == 0 and count == 32
        assert ordinal == index
        assert probability == frac(br['probability'])
        rho, pre = frac(br['conditional_model_correlation']), frac(br['prefix_conditional_model_correlation'])
        assert abs(rho) <= 1 and abs(pre) <= 1
        assert probability * rho == frac(br['contribution'])
        assert probability * pre == frac(br['prefix_contribution'])
        mass += probability
        total += probability * rho
        prefix_total += probability * pre
        nonzero += rho != 0
        prefix_nonzero += pre != 0
        if pre:
            prefix_branches.append(br)
        if br['absorbed_by_round']:
            extinction[str(br['absorbed_by_round'])] += 1
    assert len(data['branches']) == end - begin
    assert mass == frac(data['probability_mass'])
    assert total == frac(data['correlation']) == 2 * frac(data['bias'])
    assert prefix_total == frac(data['prefix_correlation'])
    assert nonzero == data['nonzero_branch_correlations']
    assert prefix_nonzero == data['nonzero_prefix_branch_correlations']
    assert dict(extinction) == data['absorbed_histogram']
    result = {k: v for k, v in data.items() if k != 'branches'}
    result['nonzero_prefix_branches'] = prefix_branches
    return result


def run_shard(row, mode, begin, end, out):
    path = out / 'jobs' / f"r{row['rounds']}_{mode}_{begin:05d}_{end:05d}.json.gz"
    if not path.exists():
        tmp = path.with_suffix('.partial.json')
        command = [str(ROOT / 'build/hnw_endpoint'), mode, str(row['rounds']), str(row['start_round_index']),
                   *row['candidate']['input_difference']['words'], *row['output_mask']['words'],
                   str(tmp), str(begin), str(end), *row['candidate']['join_mask']['words']]
        with path.with_suffix('.log').open('w') as log:
            subprocess.run(command, check=True, stdout=log, stderr=log)
        raw = tmp.read_bytes()
        data = json.loads(raw)
        result = validate(data, row, begin, end)
        compressed = path.with_suffix('.tmp')
        with gzip.open(compressed, 'wb') as f:
            f.write(raw)
        compressed.replace(path)
        tmp.unlink()
    else:
        with gzip.open(path, 'rt') as f:
            result = validate(json.load(f), row, begin, end)
    result['artifact'] = str(path)
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source', type=Path, default=ROOT / 'results/fixed_tail_extension')
    p.add_argument('--workers', type=int, default=24)
    p.add_argument('--shard-size', type=int, default=1024)
    args = p.parse_args()
    assert args.workers > 0 and args.shard_size > 0
    source = args.source.resolve()
    out = source / 'hnw26'
    (out / 'jobs').mkdir(parents=True, exist_ok=True)
    original = json.loads((source / 'results.json').read_text())
    ids = ['D1_M1_D3_chi1_to_3', 'D1_M1_D3_chi1_to_4']
    rows = []
    for rid in ids:
        row = next(r for r in original['rows'] if r['id'] == rid)
        selected = {k: v for k, v in row.items() if k != 'candidates'}
        selected['candidate'] = row['candidates'][0]
        assert selected['candidate']['input_bit_weight'] == 8
        rows.append(selected)
    sources = [source / 'results.json', ROOT / 'src/hnw_endpoint.cpp', ROOT / 'build/hnw_endpoint',
               ROOT / 'include/gaston.hpp', ROOT / 'include/estimator.hpp', Path(__file__)]
    config = {'rows': rows, 'rd': 1, 'rl': 0, 'first_L_omitted': True,
              'reference_round_constants_retained': True, 'shard_size': args.shard_size,
              'modes': ['full', 'zero'], 'prefix_readout_rounds': 2,
              'expected_first_DDT_branches_per_endpoint_and_mode': 65536,
              'method': 'HNW26 local-column approximation; full first-DDT sum; multi-column product readout',
              'sources': [{'path': str(f), 'sha256': hashlib.sha256(f.read_bytes()).hexdigest()} for f in sources]}
    config_file = out / 'configuration.json'
    if config_file.exists():
        assert json.loads(config_file.read_text()) == config, 'Configuration changed; use another output directory.'
    else:
        save(config_file, config)
    backend = json.loads(subprocess.check_output([str(ROOT / 'build/hnw_endpoint'), '--selftest'], text=True))
    assert backend['status'] == 'PASS'
    save(out / 'backend_validation.json', backend)
    began = time.perf_counter()
    tasks = [(row, mode, begin, min(begin + args.shard_size, 65536))
             for row in rows for mode in ('full', 'zero') for begin in range(0, 65536, args.shard_size)]
    finished = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        jobs = {pool.submit(run_shard, *task, out): task for task in tasks}
        for future in concurrent.futures.as_completed(jobs):
            result = future.result()
            finished.append(result)
            progress = {'completed_shards': len(finished), 'total_shards': len(tasks),
                        'completed_branch_evaluations': sum(r['completed_branches'] for r in finished),
                        'elapsed_seconds': time.perf_counter() - began}
            save(out / 'progress.json', progress)
            print(json.dumps(progress), flush=True)
    combined = []
    for row in rows:
        rr = {'row_id': row['id'], 'rounds': row['rounds'], 'tail_id': row['tail_id'],
              'start_round_index': row['start_round_index'], 'endpoint_id': row['candidate']['endpoint_id'],
              'DATF': row['candidate']['full_DATF']['correlation'],
              'exact_prefix': row['candidate']['exact_prefix'], 'HNW': {}}
        for mode in ('full', 'zero'):
            shards = sorted((s for s in finished if s['rounds'] == row['rounds'] and s['mode'] == mode),
                            key=lambda s: s['branch_begin'])
            assert shards[0]['branch_begin'] == 0 and shards[-1]['branch_end'] == 65536
            assert all(a['branch_end'] == b['branch_begin'] for a, b in zip(shards, shards[1:]))
            mass = sum((frac(s['probability_mass']) for s in shards), Fraction(0))
            total = sum((frac(s['correlation']) for s in shards), Fraction(0))
            prefix = sum((frac(s['prefix_correlation']) for s in shards), Fraction(0))
            assert mass == 1
            hist = Counter()
            for s in shards:
                hist.update(s['absorbed_histogram'])
            rr['HNW'][mode] = {'correlation': number(total), 'bias': number(total / 2),
                               'prefix_correlation': number(prefix), 'prefix_bias': number(prefix / 2),
                               'probability_mass': number(mass), 'completed_branches': 65536,
                               'nonzero_branch_correlations': sum(s['nonzero_branch_correlations'] for s in shards),
                               'nonzero_prefix_branch_correlations': sum(s['nonzero_prefix_branch_correlations'] for s in shards),
                               'absorbed_histogram': dict(hist),
                               'nonzero_prefix_branches': [b for s in shards for b in s['nonzero_prefix_branches']],
                               'shards': shards}
        combined.append(rr)
    result = {'status': 'completed', 'rows': combined, 'first_L_omitted': True,
              'arithmetic': 'GMP exact dyadic rational; no floating zero threshold',
              'scope': 'Model correlations of the complete endpoints, not exact full-permutation correlations',
              'elapsed_seconds': time.perf_counter() - began}
    save(out / 'results.json', result)
    save(out / 'validation.json', {'status': 'PASS', 'branch_evaluations': 4 * 65536,
         'independent_DDT_probability_and_ordinal_checks': 4 * 65536,
         'signed_sums_and_probability_masses': 'PASS', 'full_branch_coverage': 'PASS',
         'prefix_comparisons_to_exact': [frac(r['HNW'][m]['prefix_correlation']) == frac(r['exact_prefix']['correlation'])
                                        for r in combined for m in ('full', 'zero')],
         'backend_selftest': backend, 'unknown_or_truncated_results': False})
    print('RESULT', out / 'results.json', flush=True)


if __name__ == '__main__':
    main()
