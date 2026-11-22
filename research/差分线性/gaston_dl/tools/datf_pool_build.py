#!/usr/bin/env python3
"""Build explicit, rotation-reduced Gaston differential and linear trail pools.

Pool rotation is a combinatorial reduction of DDT/LAT trails, not a claim
that DATF or the fixed-constant permutation is invariant under rotation.
"""
import argparse
import gzip
import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path

from fixed_tail_search import S, lat, column, embed, rotate, linear, transpose
from pq_kernel_six_search import decode, validate_pq
from six_six_prefix_search import assemble, ddt_table, second_round_rank_certificate
from linear_two_pool import POOL_KEY, build_two_pool, persist_pool, render_report


def words(state):
    return [f'{w:016x}' for w in state]


def state(words_):
    return tuple(int(w, 16) for w in words_)


def weight(a):
    return sum(w.bit_count() for w in a)


def bounded_representatives(cap):
    if cap not in (1, 2):
        raise ValueError('Explicit bounded pools currently support caps 1 and 2 only.')
    seeds = {embed(0, 1 << i) for i in range(5)}
    if cap == 2:
        seeds.update(embed(0, (1 << i) | (1 << j))
                     for i, j in itertools.combinations(range(5), 2))
        for j, r, s in itertools.product(range(1, 64), range(5), range(5)):
            seeds.add(tuple(a | b for a, b in zip(embed(0, 1 << r), embed(j, 1 << s))))
    return sorted({min(rotate(a, s) for s in range(64)) for a in seeds})


def single_pool(kind, cap, limit_boundary):
    """Complete nonzero local transitions for one bounded endpoint.

The canonical pair compares the bounded endpoint first, then the other end.
Stabilizers of weight-two states are handled explicitly.
"""
    ddt = ddt_table()
    table = ddt if kind == 'differential' else [[lat(a, b) for b in range(32)] for a in range(32)]
    bound_input = ((kind == 'differential' and limit_boundary != 'inner') or
                   (kind == 'linear' and limit_boundary == 'inner'))
    result = []
    for bounded in bounded_representatives(cap):
        columns = [j for j in range(64) if column(bounded, j)]
        stabilizers = [s for s in range(64) if rotate(bounded, s) == bounded]
        choices = []
        for j in columns:
            a = column(bounded, j)
            choices.append([(b, table[a][b] if bound_input else table[b][a])
                            for b in range(1, 32)
                            if (table[a][b] if bound_input else table[b][a])])
        for entries in itertools.product(*choices):
            other = assemble(columns, [b for b, _ in entries])
            if limit_boundary == 'both' and weight(other) > cap:
                continue
            if other != min(rotate(other, s) for s in stabilizers):
                continue
            a, b = (bounded, other) if bound_input else (other, bounded)
            coefficient = math.prod(v for _, v in entries)
            den_exp = 5 * len(columns)
            abs_weight = den_exp - (abs(coefficient).bit_length() - 1)
            assert abs(coefficient) & (abs(coefficient) - 1) == 0
            orbit = len({(rotate(a, s), rotate(b, s)) for s in range(64)})
            result.append({'rounds': 1, 'input': words(a), 'output': words(b),
                           'weight': abs_weight, 'sign': 1 if coefficient > 0 else -1,
                           'orbit_size': orbit, 'input_bit_weight': weight(a),
                           'output_bit_weight': weight(b)})
    result.sort(key=lambda r: (r['weight'], r['input'], r['output']))
    prefix = 'D1' if kind == 'differential' else 'L1'
    for i, r in enumerate(result):
        r['id'] = f'{prefix}_{i:05d}'
    return result


def pq_pool(patterns):
    ddt = ddt_table()
    result = []
    for p in patterns['orbits']:
        validate_pq(decode(p['theta_input']))
        d0, b1, a2 = (decode(p[k]) for k in ('first_input', 'first_output', 'second_input'))
        assert d0 == b1 and linear(b1) == a2
        assert second_round_rank_certificate(d0, a2)['ddt_product_is_exact_for_this_two_round_trail']
        positions = [j for j in range(64) if column(a2, j)]
        opts = [[b for b in range(1, 32) if ddt[column(a2, j)][b]] for j in positions]
        assert len(positions) == 6 and all(len(o) == 4 for o in opts)
        for i, bs in enumerate(itertools.product(*opts)):
            b2 = assemble(positions, bs)
            assert all(ddt[column(a2, j)][b] == 8 for j, b in zip(positions, bs))
            result.append({'id': f'D2_{p["id"]}_{i:04d}', 'rounds': 2,
                           'input': words(d0), 'first_output': words(b1),
                           'second_input': words(a2), 'output': words(b2),
                           'weight': 24, 'sign': 1, 'orbit_size': 64,
                           'pq_orbit': p['id'], 'second_output_index': i})
    return result


def write_pool(path, records):
    with gzip.open(path, 'wt', encoding='utf-8') as f:
        for record in records:
            f.write(json.dumps(record, separators=(',', ':')) + '\n')


def read_pool(path):
    with gzip.open(path, 'rt', encoding='utf-8') as f:
        return [json.loads(line) for line in f]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cap', type=int, choices=(1, 2), default=2)
    p.add_argument('--boundary', choices=('outer', 'inner', 'both'), required=True)
    p.add_argument('--patterns', type=Path, default=Path('results/pq_kernel_six_all/patterns.json'))
    p.add_argument('--tails', type=Path, default=Path('results/fixed_paper_tail/search.json'))
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    d1 = single_pool('differential', args.cap, args.boundary)
    l1 = single_pool('linear', args.cap, args.boundary)
    d2 = pq_pool(json.loads(args.patterns.read_text()))
    tails = []
    for t in json.loads(args.tails.read_text())['reports']:
        tails.append({'id': t['id'], 'rounds': t['r_l'], 'weight_squared': t['tail_w_squared'],
                      'input': words(tuple(t['tail_rounds'][0]['input'])),
                      'output': words(decode(t['output_mask'])),
                      'join_mask': words(decode(t['join_mask'])), 'trail': t['tail_rounds'],
                      'rotation': 0})
        assert transpose(state(tails[-1]['input'])) == state(tails[-1]['join_mask'])
    pools = {'differential_one': d1, 'differential_pq_two': d2, 'linear_one': l1}
    metadata = {'cap': args.cap, 'bounded_boundary': args.boundary,
                'requested_cap4_endpoint_count': sum(math.comb(320, k) for k in range(1, 5)),
                'cap_endpoint_count': sum(math.comb(320, k) for k in range(1, args.cap + 1)),
                'first_L_omitted': True, 'pools': {}, 'paper_tails': tails,
                'rotation_scope': 'Canonical DDT/LAT trail representatives. Matching expands differential rotations only; paper linear tails remain unrotated.',
                'constant_caveat': 'Relative-shift coverage is not automatically complete fixed-round-constant endpoint coverage.',
                'pq_scope': 'All 11 paired P=Q=0 six-bit orbits, unique maximum-probability first predecessor, all 4^6 second outputs.',
                'sources': [{'path': str(f.resolve()), 'sha256': hashlib.sha256(f.read_bytes()).hexdigest()}
                            for f in (args.patterns, args.tails)]}
    for name, records in pools.items():
        path = args.output / (name + '.jsonl.gz')
        write_pool(path, records)
        metadata['pools'][name] = {'file': path.name, 'representatives': len(records),
                                   'expanded_records': sum(r['orbit_size'] for r in records),
                                   'weight_histogram': dict(sorted(Counter(r['weight'] for r in records).items())),
                                   'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}
        print(name, metadata['pools'][name], flush=True)
    l2, l2_audit = build_two_pool()
    metadata['pools'][POOL_KEY] = persist_pool(args.output, l2, l2_audit)
    metadata['linear_two_scope'] = l2_audit['scope']
    render_report(args.output, l2_audit)
    print(POOL_KEY, metadata['pools'][POOL_KEY], flush=True)
    (args.output / 'pools.json').write_text(json.dumps(metadata, ensure_ascii=False, indent=2) + '\n')


if __name__ == '__main__':
    main()
