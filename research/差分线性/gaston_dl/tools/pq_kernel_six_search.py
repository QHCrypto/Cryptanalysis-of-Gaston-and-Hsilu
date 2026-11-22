#!/usr/bin/env python3
"""Exhaust six-bit P=Q=0 paired-column states and splice fixed paper tails.

The state A is at theta input, after rho_east. Both A and rho_t(A)
must have exactly three active columns of two bits each. Only a common
64-bit rotation is quotiented out. Search results are restricted signed
contributions, not the complete three-chi prefix or full-R correlation.
"""
import argparse
import collections
import hashlib
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path

from linear_seed_scan import EAST, WEST, TWIN, active, rol
from six_six_prefix_search import (
    RC, S, ZERO, assemble, chi_probability, column, ddt_table, derivative_mask,
    dot, embed, encode, flat, fraction_record, lat, linear, match_outputs,
    rotate, second_round_rank_certificate, state_chi, transpose, write_json, xor)


def decode(state):
    return tuple(int(w, 16) for w in state['words'])


def canonical(state):
    return min(rotate(state, s) for s in range(64))


def pair_columns(state):
    return [{'column': j, 'rows': [i for i in range(5) if state[i] >> j & 1]}
            for j in range(64) if column(state, j)]


def validate_pq(state):
    shifted = tuple(rol(w, t) for w, t in zip(state, TWIN))
    for s in (state, shifted):
        pairs = pair_columns(s)
        assert len(pairs) == 3 and all(len(p['rows']) == 2 for p in pairs)
        assert xor(s) == 0
    return shifted


def enumerate_cycles():
    """A simple degree-2 graph on three P and three Q vertices is a 6-cycle.

    Anchor one P vertex at zero, walk five row-labelled edges, and derive
    the sixth row from the closing Q coordinate. Enumerate all five-row
    sequences. Distinct t_i forbid parallel edges / two-edge cycles.
    """
    inverse = {t: i for i, t in enumerate(TWIN)}
    counts = collections.Counter()
    for rows in itertools.product(range(5), repeat=5):
        if any(rows[i] == rows[i + 1] for i in range(4)):
            continue
        q0 = TWIN[rows[0]]
        p1 = (q0 - TWIN[rows[1]]) % 64
        q1 = (p1 + TWIN[rows[2]]) % 64
        p2 = (q1 - TWIN[rows[3]]) % 64
        q2 = (p2 + TWIN[rows[4]]) % 64
        if len({0, p1, p2}) != 3 or len({q0, q1, q2}) != 3:
            continue
        last = inverse.get(q2)
        if last is None or last in (rows[0], rows[4]):
            continue
        state = [0] * 5
        for i, j in zip((*rows, last), (0, p1, p1, p2, p2, 0)):
            state[i] |= 1 << j
        validate_pq(state)
        counts[canonical(state)] += 1
    return counts


def enumerate_column_pairs():
    """Independent enumeration: choose 3 P-column pairs with XOR Q=0.

    One P column can always be rotated to zero. Choose its row pair and
    the next nonzero P column's row pair. Lookup all possible third
    pairs by their two-bit Q signature; require 0 < p1 < p2.
    """
    pairs = list(itertools.combinations(range(5), 2))
    atoms = collections.defaultdict(list)
    for p in range(64):
        for i, j in pairs:
            q = (1 << ((p + TWIN[i]) % 64)) ^ (1 << ((p + TWIN[j]) % 64))
            atoms[q].append((p, i, j))
    counts = collections.Counter()
    for i, j in pairs:
        q0 = (1 << TWIN[i]) ^ (1 << TWIN[j])
        for p in range(1, 64):
            for k, l in pairs:
                q1 = (1 << ((p + TWIN[k]) % 64)) ^ (1 << ((p + TWIN[l]) % 64))
                for p2, u, v in atoms.get(q0 ^ q1, []):
                    if p2 <= p:
                        continue
                    state = [0] * 5
                    for pos, r, s in ((0, i, j), (p, k, l), (p2, u, v)):
                        state[r] ^= 1 << pos
                        state[s] ^= 1 << pos
                    validate_pq(state)
                    counts[canonical(state)] += 1
    return counts


def vortex_states():
    result = {}
    for i, j, k in itertools.combinations(range(5), 3):
        state = [0] * 5
        state[i] = (1 << TWIN[j]) ^ (1 << TWIN[k])
        state[j] = (1 << TWIN[k]) ^ (1 << TWIN[i])
        state[k] = (1 << TWIN[i]) ^ (1 << TWIN[j])
        result[canonical(state)] = [i, j, k]
    return result


def make_patterns(old_pattern, ddt):
    cycles, pairs = enumerate_cycles(), enumerate_column_pairs()
    assert cycles.keys() == pairs.keys()
    vortices = vortex_states()
    assert vortices.keys() <= cycles.keys()
    old_b1 = tuple(int(w, 16) for w in old_pattern['input_words'])
    old_a = tuple(rol(w, e) for w, e in zip(old_b1, EAST))
    records = []
    for index, a in enumerate(sorted(cycles), 1):
        qa = validate_pq(a)
        b1 = tuple(rol(w, -e) for w, e in zip(a, EAST))
        a2 = tuple(rol(w, wshift) for w, wshift in zip(a, WEST))
        assert linear(b1) == a2
        assert active(b1) == active(a2) == 6
        assert len({rotate(a, s) for s in range(64)}) == 64
        assert cycles[a] == 6 and pairs[a] == 3
        # Best predecessors are unique. Any other predecessor has the same
        # active-column support, hence the same conditional rank certificate.
        for j in range(64):
            b = column(b1, j)
            if b:
                best = max(row[b] for row in ddt)
                assert [i for i in range(32) if ddt[i][b] == best] == [b]
        certificate = second_round_rank_certificate(b1, a2)
        assert certificate['ddt_product_is_exact_for_this_two_round_trail']
        assert certificate['second_derivative_rank'] == 12
        assert chi_probability(b1, b1, ddt) == Fraction(1, 1 << 12)
        records.append({
            'id': f'K{index:02d}', 'theta_input': encode(a), 'after_t_shifts': encode(qa),
            'P_pairs': pair_columns(a), 'Q_pairs': pair_columns(qa),
            'first_input': encode(b1), 'first_output': encode(b1), 'second_input': encode(a2),
            'vortex_rows': vortices.get(a), 'orbit_size': 64,
            'is_previous_orbit': canonical(old_a) == a,
            'previous_seed_rotation': next((s for s in range(64) if rotate(a, s) == old_a), None),
            'rank_certificate': certificate,
            'cycle_representations_with_P_anchor_zero': cycles[a],
            'pair_representations_with_P_anchor_zero': pairs[a]})
    return records


def as_fraction(value):
    return Fraction(value['numerator'], value['denominator'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tails', type=Path, required=True)
    parser.add_argument('--previous-pattern', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    tails = json.loads(args.tails.read_text())
    ddt = ddt_table()
    patterns = make_patterns(json.loads(args.previous_pattern.read_text()), ddt)
    write_json(args.output / 'patterns.json', {
        'scope': 'Exactly six theta-input bits, three two-bit P columns and three two-bit Q columns.',
        'bit_index': '(row,bit), zero-based, LSB=0; left rotation adds to bit index modulo 64',
        't': list(TWIN), 'east': list(EAST), 'west': list(WEST),
        'equivalence': 'Only common cyclic rotation of all five rows; rows are not permuted.',
        'canonicalization': 'Lexicographically smallest tuple of five unsigned 64-bit row words.',
        'completeness': 'A degree-two simple bipartite graph on 3+3 vertices must be a six-cycle; all row-labelled cycles are enumerated. Independent enumeration of three P-column pairs with Q-signature XOR zero returns the same states.',
        'orbits': patterns, 'orbit_count': len(patterns), 'state_count': len(patterns) * 64})
    groups, cache = {}, {}
    for tail in tails['reports']:
        gamma = decode(tail['join_mask'])
        if gamma in cache:
            continue
        gid = f'join_{len(groups)}'
        table = [[flat(derivative_mask(gamma, linear(embed(j, b)))) for b in range(32)]
                 for j in range(64)]
        scans = []
        for p in patterns:
            for shift in range(64):
                b1 = rotate(decode(p['first_output']), shift)
                a2 = rotate(decode(p['second_input']), shift)
                positions = [j for j in range(64) if column(a2, j)]
                options = [[b for b in range(32) if ddt[column(a2, j)][b]] for j in positions]
                assert all(len(o) == 4 for o in options)
                connections = []
                for values in match_outputs(positions, options, table):
                    b2 = assemble(positions, values)
                    a3 = linear(b2)
                    # Direct evaluation of all 32 local values independently
                    # verifies the middle derivative is constant and its sign.
                    local = [sum((-1) ** ((column(gamma, j) &
                               (S[x] ^ S[x ^ column(a3, j)])).bit_count() & 1)
                               for x in range(32)) for j in range(64)]
                    assert all(abs(v) == 32 for v in local)
                    sign = math.prod(v // 32 for v in local)
                    assert sign == (-1) ** dot(gamma, state_chi(a3))
                    assert chi_probability(a2, b2, ddt) == Fraction(1, 1 << 12)
                    connections.append({'b2': encode(b2), 'a3': encode(a3), 'middle_correlation': sign})
                positive = sum(c['middle_correlation'] == 1 for c in connections)
                negative = len(connections) - positive
                scans.append({'orbit': p['id'], 'rotation': shift, 'positive': positive, 'negative': negative,
                              'connections': connections,
                              'restricted_signed_prefix_sum': fraction_record(Fraction(positive - negative, 1 << 24))})
        groups[gid] = {'gamma': encode(gamma), 'states_tested': len(scans),
                       'implicit_second_DDT_combinations': len(scans) * 4096,
                       'feasible_trails': sum(len(s['connections']) for s in scans),
                       'nonempty_input_count': sum(bool(s['connections']) for s in scans),
                       'noncancelling_input_count': sum(s['positive'] != s['negative'] for s in scans),
                       'scans': scans}
        cache[gamma] = gid
        print(gid, 'trails', groups[gid]['feasible_trails'],
              'noncancelling inputs', groups[gid]['noncancelling_input_count'], flush=True)
    reports = []
    for tail in tails['reports']:
        gamma = decode(tail['join_mask'])
        gid = cache[gamma]
        new_R = tail['R'] + 1
        signed_tail = Fraction(1)
        previous = gamma
        assert RC[12 - new_R + 3:] == RC[tail['start_round_index'] + 2:]
        for i, rr in enumerate(tail['tail_rounds']):
            a, b = tuple(rr['input']), tuple(rr['output'])
            assert transpose(a) == previous
            c = math.prod(Fraction(lat(column(a, j), column(b, j)), 32) for j in range(64))
            signed_tail *= c * (-1) ** ((a[0] & RC[12 - new_R + 3 + i]).bit_count() & 1)
            previous = b
        assert signed_tail ** 2 == Fraction(1, 1 << tail['tail_w_squared'])
        nonzero = [s for s in groups[gid]['scans'] if s['positive'] != s['negative']]
        best = max(nonzero, key=lambda s: abs(s['positive'] - s['negative'])) if nonzero else None
        best_record = None
        if best:
            p = next(p for p in patterns if p['id'] == best['orbit'])
            best_record = {k: best[k] for k in ('orbit', 'rotation', 'positive', 'negative', 'restricted_signed_prefix_sum')}
            best_record.update({
                'd0': encode(rotate(decode(p['first_input']), best['rotation'])),
                'b1': encode(rotate(decode(p['first_output']), best['rotation'])),
                'a2': encode(rotate(decode(p['second_input']), best['rotation'])),
                'restricted_splice_score': fraction_record(as_fraction(best['restricted_signed_prefix_sum']) * signed_tail ** 2)})
        report = {'id': tail['id'], 'R': new_R, 'r_d': 2, 'r_m': 1, 'r_l': tail['r_l'],
                  'start_round_index': 12 - new_R, 'connection_group': gid,
                  'tail_w_squared': tail['tail_w_squared'], 'tail_signed_LAT_product': fraction_record(signed_tail),
                  'tail_rounds': tail['tail_rounds'], 'output_mask': tail['output_mask'],
                  'feasible_trails': groups[gid]['feasible_trails'],
                  'noncancelling_input_count': len(nonzero), 'best': best_record}
        reports.append(report)
        print(report['id'], 'R', new_R, 'best', best_record and
              (best_record['orbit'], best_record['rotation'], best_record['restricted_splice_score']), flush=True)
    best_by_round = []
    for R in sorted({r['R'] for r in reports}):
        candidates = [r for r in reports if r['R'] == R and r['best']]
        if candidates:
            best_by_round.append(max(candidates, key=lambda r: abs(as_fraction(r['best']['restricted_splice_score'])))['id'])
    data = {
        'sources': [{'path': str(p.resolve()), 'sha256': hashlib.sha256(p.read_bytes()).hexdigest()}
                    for p in (args.tails, args.previous_pattern, Path(__file__))],
        'first_L_omitted': True, 'differential_active_counts': [6, 6],
        'scope': 'All paired P=Q=0 six-bit theta states, all rotations, unique best first-DDT predecessor, all second-DDT outputs, 11 fixed published tails.',
        'orbit_count': len(patterns), 'state_count': len(patterns) * 64,
        'p_first': fraction_record(Fraction(1, 1 << 12)),
        'p_second_conditioned_on_first': fraction_record(Fraction(1, 1 << 12)),
        'p_two_round_trail': fraction_record(Fraction(1, 1 << 24)),
        'groups': groups, 'reports': reports, 'best_by_round': best_by_round,
        'checks': ['Cycle and independent P-column-pair enumeration agree.',
                   'All 11 representatives have rank 12 on unaffected first-output columns; common rotations preserve rank.',
                   'All connecting local middle derivatives enumerated on all 32 inputs.',
                   'Fixed tail transpose, LAT weights and round constants rechecked.'],
        'limits': ['Only the six-bit paired-column P=Q=0 subset; not the general E=0 kernel or all possible 6->6 column patterns.',
                   'Input predecessor is locally optimal for each fixed B1 event. Its active support fixes the rank certificate, so other predecessors cannot increase this restricted event sum.',
                   'Only fixed B1 and middle DLCT +/-1 branches are summed. Other first outputs and conditional contributions of DLCT-zero branches are not evaluated.',
                   'The tail uses one fixed trail pair. These scores are neither full-prefix/full-R correlations nor bounds on them.']}
    write_json(args.output / 'search.json', data)
    # Choose one best example from each nonempty join mask, plus a second
    # B2 branch of the improved five-round candidate to test both contributions.
    jobs, job_records = [], []
    seen = set()
    for r in reports:
        if not r['best'] or r['connection_group'] in seen:
            continue
        seen.add(r['connection_group'])
        b = r['best']
        gr = groups[r['connection_group']]
        scan = next(s for s in gr['scans'] if (s['orbit'], s['rotation']) == (b['orbit'], b['rotation']))
        take = 2 if r['id'] == 'D3_chi2_to_3' else 1
        for n, c in enumerate(scan['connections'][:take]):
            jid = f'{r["id"]}_{b["orbit"]}_rot{b["rotation"]}_b{n}'
            states = [b['d0'], b['b1'], c['b2'], gr['gamma']]
            jobs.append(f'{jid} {r["start_round_index"]} {c["middle_correlation"]} '
                        + ' '.join(w for s in states for w in s['words']))
            job_records.append({'id': jid, 'report': r['id'], 'orbit': b['orbit'],
                                'rotation': b['rotation'], 'connection_index': n})
    (args.output / 'conditional_jobs.txt').write_text('\n'.join(jobs) + '\n')
    write_json(args.output / 'conditional_job_map.json', job_records)
    print('Saved', args.output.resolve(), flush=True)


if __name__ == '__main__':
    main()
