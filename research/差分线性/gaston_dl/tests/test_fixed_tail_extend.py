#!/usr/bin/env python3
"""Independent counting identities, complete saved witnesses, actual round boundaries."""
import itertools
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools'))
from fixed_tail_search import S, dot, checked_round
from linear_seed_scan import column, embed, linear, transpose, xor
from datf_pool_build import state
from fixed_tail_extend import connector, local_weight, coefficient, theta_effect

ROOT = Path(__file__).resolve().parents[1]/'results/fixed_tail_extension'


def packed(a):
    return sum(w << (64*i) for i, w in enumerate(a))


def truth_syndrome(delta, gamma):
    values = []
    for j in range(64):
        a, m = column(delta, j), column(gamma, j)
        f0 = (m & (S[0] ^ S[a])).bit_count() & 1
        values.append(sum((((m & (S[1<<k] ^ S[a ^ (1<<k)])).bit_count() & 1) ^ f0) << k
                          for k in range(5)))
    return packed(tuple(sum(((v >> i) & 1) << j for j, v in enumerate(values)) for i in range(5)))


def main():
    groups = json.loads((ROOT/'groups.json').read_text())
    audits, checked_seeds, checked_candidates = [], 0, 0
    for group in groups:
        raw = json.loads((ROOT/(group['id']+'_search.json')).read_text())
        gamma = state(group['gamma'])
        syndromes = [truth_syndrome(linear(embed(j, 1<<i)), gamma)
                     for i in range(5) for j in range(64)]
        basis = {}
        for s in syndromes:
            while s:
                k = s.bit_length()-1
                if k in basis:
                    s ^= basis[k]
                else:
                    basis[k] = s
                    break
        assert len(basis) == raw['connector_constraint_rank']
        pairs = Counter(syndromes[a] ^ syndromes[b]
                        for a in range(320) for b in range(a))
        n1, n2 = syndromes.count(0), pairs[0]
        # Each distinct triple has three singleton/pair partitions.
        triple_partition_count = sum(pairs[s] for s in syndromes) - n1*319
        assert triple_partition_count % 3 == 0
        n3 = triple_partition_count//3
        # Each distinct quadruple has three pair/pair partitions. Intersecting
        # pairs are counted once for each zero-syndrome pair and remaining bit.
        disjoint_pair_count = sum(n*(n-1)//2 for n in pairs.values())-n2*318
        assert disjoint_pair_count % 3 == 0
        n4 = disjoint_pair_count//3
        assert [0, n1, n2, n3, n4] == raw['one_round_kernel_counts_by_bit_weight']
        for rd, key in ((1, 'one_round_seeds'), (2, 'two_round_seeds')):
            for r in raw[key]:
                b1, b2 = state(r['first_output']), state(r['second_output'])
                b = b1 if rd == 1 else b2
                assert truth_syndrome(linear(b), gamma) == 0
                if rd == 2:
                    local_weight(linear(b1), b2)
                checked_seeds += 1
        for tail in group['tails']:
            assert tail['rounds'] == len(tail['trail']) >= 3
            assert transpose(state(tail['input'])) == gamma
            assert state(tail['input']) == tuple(tail['trail'][0]['input'])
            assert state(tail['output']) == tuple(tail['trail'][-1]['output'])
            total = 0
            for i, r in enumerate(tail['trail']):
                total += checked_round(tuple(r['input']), tuple(r['output']))['w_squared']
                if i:
                    assert tuple(tail['trail'][i-1]['output']) == transpose(tuple(r['input']))
            assert total == tail['weight_squared']
        audits.append({'group_id': group['id'], 'independent_rank': len(basis),
                       'independent_weight_1_to_4_counts': [n1, n2, n3, n4]})
    candidates = json.loads((ROOT/'candidates.json').read_text())
    for row in candidates['rows']:
        assert row['rounds'] == row['r_d']+1+row['r_l']
        assert row['start_round_index']+row['rounds'] == 12
        for c in row['candidates']:
            delta = state(c['input_difference']['words'])
            gamma = state(c['join_mask']['words'])
            middle = state(c['middle_input_difference']['words'])
            rd = row['r_d']
            weight = 0
            if rd:
                b1 = state(c['first_output']['words'])
                weight += local_weight(delta, b1)
                if rd == 2:
                    a2 = state(c['second_input']['words'])
                    b2 = state(c['second_output']['words'])
                    assert a2 == linear(b1)
                    weight += local_weight(a2, b2)
                    assert middle == linear(b2)
                else:
                    assert middle == linear(b1)
            else:
                assert middle == delta
            assert weight == c['differential_weight']
            assert weight+row['tail_squared_weight'] == c['single_characteristic_weight']
            # Direct local truth sums, independent of the affine coefficient rule.
            value = 1
            for j in range(64):
                a, m = column(middle, j), column(gamma, j)
                s = sum((-1)**((m & (S[x]^S[x^a])).bit_count() & 1) for x in range(32))
                assert abs(s) == 32
                value *= s//32
            assert value == c['middle_correlation']
            checked_candidates += 1
    result = {'status': 'PASS', 'group_audits': audits, 'saved_search_seeds_checked': checked_seeds,
              'selected_characteristics_checked': checked_candidates,
              'tail_trails_checked': 7, 'first_L_omitted': True,
              'scope': 'Independent truth-table syndromes, combinatorial counts, all retained DDT/DLCT witnesses and tail boundaries; not an exact full correlation proof.'}
    if (ROOT/'results.json').exists():
        final = json.loads((ROOT/'results.json').read_text())
        assert final['unique_full_jobs'] == 893 and final['resource_limited'] == 0
        assert final['unique_model_nonzero'] == 2
        dense = next(r for r in final['rows'] if r['id']=='D1_M1_D3_chi1_to_3')['candidates'][0]
        effect = theta_effect(state(dense['first_output']['words']))
        assert int(effect['E'], 16) == 0 and int(effect['P'], 16) and int(effect['Q'], 16)
        assert dense['input_difference'] == dense['first_output']
        assert len(dense['input_difference']['bits']) == 8
        assert len(dense['middle_input_difference']['bits']) == 8
        exact = dense['exact_prefix']
        assert exact['correlation']['numerator']=='1' and exact['correlation']['denominator']=='65536'
        assert exact['branch_count']==65536 and exact['positive_contributions']==1
        assert exact['negative_contributions']==0 and exact['agrees_with_DATF'] is False
        assert dense['prefix_DATF']['correlation']['numerator']=='0'
        result.update(exact_eight_bit_prefix='1/65536',
                      ordinary_DATF_prefix='0; explicitly reported model discrepancy',
                      theta_effect=effect)
    (ROOT/'validation.json').write_text(json.dumps(result, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(result, ensure_ascii=False))


if __name__ == '__main__':
    main()
