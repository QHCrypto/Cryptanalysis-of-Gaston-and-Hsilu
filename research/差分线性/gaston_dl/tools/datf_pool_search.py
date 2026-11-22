#!/usr/bin/env python3
"""Search explicit trail pools with ordinary DATF middle segments.

One-round middles against all paper tails are exhaustive. Two-round middles
use a deterministic, documented batch of differential representatives, with
all rotations of each selected representative. The one-round linear pool is
searched in exact weight strata for one middle round, and with a separate
representative/pair budget for two middle rounds.
The two-round 2->2 linear pool follows the same strategy with its own budget.
"""
import argparse
import collections
import csv
import hashlib
import json
import subprocess
import time
from fractions import Fraction
from pathlib import Path

from datf_pool_build import read_pool, state, words
from fixed_tail_search import RC, linear, transpose, rotate
from datf_k06 import value_record
from linear_two_pool import POOL_KEY, FAMILY


def spread(record):
    a = linear(state(record['output']))
    return ((a[0] | a[1] | a[2] | a[3] | a[4]).bit_count(), sum(w.bit_count() for w in a), record['id'])


def round_robin(groups):
    result = []
    for index in range(max(map(len, groups), default=0)):
        for group in groups:
            if index < len(group):
                result.append(group[index])
    return result


def ordered_differentials(d1, d2, previous):
    # Keep every maximum-probability one-round trail first, then interleave
    # weights 3 and 4. This ordering only selects the bounded middle-2 batch.
    one = sorted([d for d in d1 if d['weight'] == 2], key=spread)
    one += round_robin([sorted([d for d in d1 if d['weight'] == w], key=spread) for w in (3, 4)])
    favored = set()
    for group in previous['groups'].values():
        for scan in group['scans']:
            for c in scan['connections']:
                b = tuple(int(w, 16) for w in c['b2']['words'])
                favored.add((scan['orbit'], rotate(b, -scan['rotation'])))
    grouped = collections.defaultdict(list)
    for d in d2:
        grouped[d['pq_orbit']].append(d)
    groups = [sorted(grouped[k], key=lambda d: ((d['pq_orbit'], state(d['output'])) not in favored, *spread(d)))
              for k in sorted(grouped)]
    return one + round_robin(groups)


def linear_candidates(pools, meta):
    """Canonical representatives; only the differential side is rotated."""
    result = []
    for key, family in [('linear_one', 'L1'), (POOL_KEY, FAMILY)]:
        group = []
        for l in pools.get(key, []):
            join = words(transpose(state(l['input'])))
            if 'join_mask' in l:
                assert join == l['join_mask']
            group.append({**l, 'family': family, 'weight_squared': 2*l['weight'], 'join_mask': join})
        group.sort(key=lambda l:(l['weight_squared'],
                                (lambda a:(a[0]|a[1]|a[2]|a[3]|a[4]).bit_count())(state(l['join_mask'])),
                                sum(w.bit_count() for w in state(l['join_mask'])), l['id']))
        result.extend(group)
    return result + [{**l, 'family': 'paper'} for l in meta['paper_tails']]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pools', type=Path, required=True)
    p.add_argument('--previous', type=Path, default=Path('results/pq_kernel_six_all/search.json'))
    p.add_argument('--middle2-representatives', type=int, default=64)
    p.add_argument('--linear1-pair-budget', type=int, default=200000)
    p.add_argument('--linear1-middle2-representatives', type=int, default=64)
    p.add_argument('--linear2-pair-budget', type=int, default=200000)
    p.add_argument('--linear2-middle2-representatives', type=int, default=64)
    p.add_argument('--export-only', action='store_true', help='Prepare all pools and native input files without running a search')
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    if min(args.middle2_representatives, args.linear1_pair_budget, args.linear1_middle2_representatives,
           args.linear2_pair_budget, args.linear2_middle2_representatives) < 1:
        p.error('budgets must be positive')
    args.output.mkdir(parents=True, exist_ok=False)
    meta = json.loads((args.pools / 'pools.json').read_text())
    pools = {k: read_pool(args.pools / v['file']) for k, v in meta['pools'].items()}
    ds = ordered_differentials(pools['differential_one'], pools['differential_pq_two'], json.loads(args.previous.read_text()))
    ls = linear_candidates(pools, meta)
    with (args.output / 'differentials.txt').open('w') as f:
        for d in ds:
            f.write(f'{d["id"]} D{d["rounds"]} {d["weight"]} {d["orbit_size"]} ' + ' '.join(d['input']+d['output'])+'\n')
    with (args.output / 'linears.txt').open('w') as f:
        for l in ls:
            f.write(f'{l["id"]} {l["family"]} {l["rounds"]} {l["weight_squared"]} ' + ' '.join(l['join_mask']+l['output'])+'\n')
    command = ['build/datf_pool_match', str(args.output/'differentials.txt'), str(args.output/'linears.txt'),
               str(args.output/'match'), str(args.middle2_representatives), str(args.linear1_pair_budget),
               str(args.linear1_middle2_representatives), str(args.linear2_pair_budget),
               str(args.linear2_middle2_representatives)]
    (args.output/'input_manifest.json').write_text(json.dumps({
        'status': 'exported_only' if args.export_only else 'prepared',
        'differential_representatives': len(ds),
        'linear_families': dict(collections.Counter(l['family'] for l in ls)),
        'pool_metadata': meta, 'command': command}, ensure_ascii=False, indent=2)+'\n')
    if args.export_only:
        print('EXPORTED', args.output, dict(collections.Counter(l['family'] for l in ls)), flush=True)
        return
    began = time.perf_counter()
    subprocess.run(command, check=True)
    lookup_d = {d['id']:d for d in ds}; lookup_l = {l['id']:l for l in ls}
    found = collections.defaultdict(list)
    for h in csv.DictReader((args.output/'match_hits.tsv').open(), delimiter='\t'):
        d, l = lookup_d[h['diff_id']], lookup_l[h['linear_id']]
        shift, middle = int(h['rotation']), int(h['middle'])
        rd, rl = d['rounds'], l['rounds']; R = rd+middle+rl
        delta = rotate(state(d['input']), shift)
        dout = rotate(state(d['output']), shift)
        score = Fraction(int(h['sign']), 1<<int(h['weight']))
        record = {'differential_id': d['id'], 'linear_id': l['id'], 'differential_rotation': shift,
                  'r_d':rd, 'r_m':middle, 'r_l':rl, 'R':R, 'start_round_index':12-R,
                  'input_difference':words(delta), 'differential_output':words(dout),
                  'middle_input_difference':words(linear(dout)), 'middle_output_mask':l['join_mask'],
                  'output_mask':l['output'], 'differential_probability':value_record(Fraction(1,1<<d['weight'])),
                  'linear_squared_correlation':value_record(Fraction(1,1<<l['weight_squared'])),
                  'middle_DATF_correlation':value_record(Fraction(int(h['sign']),1<<int(h['middle_weight']))),
                  'score':value_record(score), 'score_weight':int(h['weight']),
                  'middle_second_constant':f'{RC[11-rl]:02x}' if middle==2 else None,
                  'middle_constant_phase_mask':h['rc_phase_mask'],
                  'interpretation':'individual trail-pair splice score; not full endpoint correlation or a signed hull sum'}
        found[h['family'],middle,h['tail']].append(record)
    summaries=[]
    for row in csv.DictReader((args.output/'match_summary.tsv').open(), delimiter='\t'):
        family,middle,tail=row['family'],int(row['middle']),row['tail']
        dcount=meta['pools']['differential_one' if family=='D1' else 'differential_pq_two']['expanded_records']
        pool_key = {'L1_pool': 'linear_one', FAMILY+'_pool': POOL_KEY}.get(tail)
        lcount = len(pools[pool_key]) if pool_key else 1
        exhaustive=row['exhausted']=='1'; certified=row['optimum_certified']=='1'
        candidates=found[family,middle,tail]
        summaries.append({'family':family,'middle_rounds':middle,'tail':tail,
                          'candidate_pair_scope':dcount*lcount,'evaluated_pairs':int(row['evaluated']),
                          'nonzero_evaluated_pairs':int(row['nonzero']),
                          'exhaustive':exhaustive,'optimal_score_certified':certified,
                          'scope_status':'exhaustive' if exhaustive else 'score_optimum_certified' if certified else 'bounded_search',
                          'score_upper_bound':value_record(Fraction(1,1<<int(row['upper_bound_weight']))),
                          'bound_basis':'exact ascending-weight-stratum search' if pool_key and middle==1 else '|rho_middle|<=1 and minimum pool weights',
                          'best':candidates[0] if candidates else None,'top_candidates':candidates})
    result={'method':'differential/linear trail pools + ordinary DATF Algorithm 1 for middle 1 or 2',
            'first_L_omitted':True, 'single_round_cap':meta['cap'],'bounded_boundary':meta['bounded_boundary'],
            'enhanced_DATF':False, 'middle_Bias_backend':'exact uniform quadratic elimination; same factors and substitutions as ordinary Algorithm 1',
            'middle_boundary':'delta_M=L(delta_D_output), gamma_M=L^T(alpha_L_first_chi); these two L layers are included, not omitted',
            'rotation':'differential representatives expanded through their complete cyclic orbit; linear representatives fixed',
            'constant_dependence':'For this one/two-chi ordinary-DATF implementation, middle-2 constants affect only the recorded sign phase. Absolute individual scores are invariant under simultaneous rotation. This does not establish invariance of signed hull sums or true full correlations.',
            'middle2_batch':{'representatives_per_differential_family':args.middle2_representatives,
                             'priority':'D1: all weight-2 trails then interleaved weight-3/4, ordered by middle-input activity. D2: round-robin 11 PQ classes, prior one-middle connectors first, then middle-input activity.',
                             'all_rotations_of_selected_representatives':True},
            'linear1_pair_budget_per_family_and_middle':args.linear1_pair_budget,
            'linear1_middle1_search':'exact ascending weight strata; all lower-weight strata are exhausted before returning a nonzero best score',
            'linear1_middle2_representatives':min(args.linear1_middle2_representatives,len(pools['linear_one'])),
            'linear1_middle2_priority':'squared LAT weight, active columns of L^T(input mask), bit weight of L^T(input mask), identifier',
            'linear2_pair_budget_per_family_and_middle':args.linear2_pair_budget,
            'linear2_middle1_search':'exact ascending weight strata; same stopping rule as the one-round linear pool',
            'linear2_middle2_representatives':min(args.linear2_middle2_representatives,len(pools.get(POOL_KEY, []))),
            'linear2_middle2_priority':'squared LAT weight, active columns of L^T(input mask), bit weight of L^T(input mask), identifier',
            'limits':['Every score is for an individual differential trail and one identical linear trail pair.',
                      'No signed aggregation over differential or linear trail pools is used for ranking.',
                      'Two-middle bounded searches do not establish absence of better unevaluated candidates.',
                      'An upper-bound certificate refers only to this individual-score objective and the stated pools.',
                      'Ordinary DATF uses independence assumptions; zero is not a proof of zero true correlation.'],
            'pool_metadata':meta,'summaries':summaries,'elapsed_seconds':time.perf_counter()-began,
            'command':command,'sources':[{'path':str((args.pools/'pools.json').resolve()),'sha256':hashlib.sha256((args.pools/'pools.json').read_bytes()).hexdigest()}]}
    (args.output/'results.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n')
    print('SAVED',args.output,'seconds',result['elapsed_seconds'],flush=True)


if __name__=='__main__':
    main()
