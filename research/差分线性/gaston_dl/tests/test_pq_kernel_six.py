#!/usr/bin/env python3
"""Independent checks of complete paired-PQ enumeration and saved splices."""
import itertools
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'tools'))
from pq_kernel_six_search import decode, make_patterns
from six_six_prefix_search import (
    MASK, ZERO, active_mask, assemble, column, ddt_table, derivative_mask,
    embed, flat, linear, rank, rotate, state_chi, sxor)


def main():
    directory = ROOT / 'results/pq_kernel_six_all'
    patterns = json.loads((directory / 'patterns.json').read_text())
    data = json.loads((directory / 'search.json').read_text())
    old = json.loads((ROOT / 'results/six_six_fixed_tail/search.json').read_text())
    source = ROOT.parents[1] / 'gaston_active_pattern/result/diff_in6.json'
    ddt = ddt_table()
    fresh = make_patterns(json.loads(source.read_text()), ddt)
    assert fresh == patterns['orbits'] and len(fresh) == 11
    assert sum(p['vortex_rows'] is not None for p in fresh) == 10
    assert len({rotate(decode(p['theta_input']), s) for p in fresh for s in range(64)}) == 704
    # Use only forward L and explicit chi derivatives, independent of L^T
    # and the derivative-coefficient formula used in the search certificate.
    for p in fresh:
        d0, a2 = decode(p['first_input']), decode(p['second_input'])
        inactive = MASK ^ active_mask(d0)
        columns = []
        for j in range(64):
            if inactive >> j & 1:
                for i in range(5):
                    y = linear(embed(j, 1 << i))
                    columns.append(flat(sxor(sxor(state_chi(y), state_chi(sxor(y, a2))), state_chi(a2))))
        assert rank(columns) == 12
    # Exhaustive direct 4^6 checks for every input with at least one match,
    # plus rotation zero of every orbit for each mask (including dense/no-hit masks).
    brute_cases = 0
    for gr in data['groups'].values():
        gamma = decode(gr['gamma'])
        for s in gr['scans']:
            if not s['connections'] and s['rotation'] != 0:
                continue
            p = next(p for p in fresh if p['id'] == s['orbit'])
            a2 = rotate(decode(p['second_input']), s['rotation'])
            positions = [j for j in range(64) if column(a2, j)]
            options = [[b for b in range(32) if ddt[column(a2, j)][b]] for j in positions]
            brute = set()
            for values in itertools.product(*options):
                b2 = assemble(positions, values)
                if derivative_mask(gamma, linear(b2)) == ZERO:
                    brute.add(b2)
            assert brute == {decode(c['b2']) for c in s['connections']}
            brute_cases += 1
    # Recover every old record after converting its rotation convention.
    previous = next(p for p in fresh if p['is_previous_orbit'])
    for oldgr in old['groups'].values():
        newgr = next(g for g in data['groups'].values() if g['gamma'] == oldgr['gamma'])
        for oldscan in oldgr['rotations']:
            shift = (oldscan['rotation'] + previous['previous_seed_rotation']) % 64
            newscan = next(s for s in newgr['scans'] if s['orbit'] == previous['id'] and s['rotation'] == shift)
            expected = {(decode(c['b2']), c['middle_correlation']) for c in oldscan['connections']}
            actual = {(decode(c['b2']), c['middle_correlation']) for c in newscan['connections']}
            assert actual == expected
            assert newscan['restricted_signed_prefix_sum'] == oldscan['restricted_signed_prefix_sum']
    expected = {5: -61, 6: -83, 7: -118, 8: -142, 9: -166}
    for rid in data['best_by_round']:
        r = next(r for r in data['reports'] if r['id'] == rid)
        assert r['best']['restricted_splice_score']['log2_abs'] == expected[r['R']]
    print(f'PASS 11 complete orbits / 704 states; 11 independent forward-L rank checks; '
          f'{brute_cases} exhaustive 4^6 comparisons; all 256 prior rotation/mask records reproduced; R=5..9 scores')


if __name__ == '__main__':
    main()
