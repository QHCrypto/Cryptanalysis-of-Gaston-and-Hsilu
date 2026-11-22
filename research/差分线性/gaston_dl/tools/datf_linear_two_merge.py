#!/usr/bin/env python3
"""Audit grouped coverage, expand all nonzero trail pairs and rank endpoints."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from pathlib import Path

from datf_pool_build import read_pool
from datf_paper_parallel import atomic_json
from datf_paper_merge import candidate_record,endpoint_id


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    args=p.parse_args();root=args.results
    config=json.loads((root/'configuration.json').read_text())
    done=json.loads((root/'search_completed.json').read_text())
    for src in config['sources']:
        assert hashlib.sha256(Path(src['path']).read_bytes()).hexdigest()==src['sha256']
    ds={};family_order=collections.Counter()
    for i,line in enumerate((root/'differentials.txt').read_text().splitlines()):
        t=line.split();family=t[1]
        ds[t[0]]={'id':t[0],'family':family,'rounds':int(family[-1]),'weight':int(t[2]),
                  'orbit':int(t[3]),'input':t[4:9],'output':t[9:14],'order':i,
                  'family_order':family_order[family]}
        family_order[family]+=1
    ls={r['id']:{**r,'weight':r['weight_squared'],'gamma':r['join_mask']}
        for r in read_pool(Path(config['sources'][1]['path']))}
    groups={g['id']:g for g in json.loads((root/'linear_groups.json').read_text())}
    all_ids=[lid for g in groups.values() for lid in g['members']]
    assert len(all_ids)==len(set(all_ids))==len(ls)==16000
    for g in groups.values():
        assert all(ls[lid]['gamma']==g['join_mask'] for lid in g['members'])
        assert all(ls[lid]['rounds']==2 for lid in g['members'])
    limit=config['middle2_representatives_per_family'];workers=config['workers']
    def selected(d,m):return m==1 or not limit or d['family_order']<limit
    totals=collections.defaultdict(collections.Counter);found=collections.defaultdict(list)
    assert {s['shard'] for s in done['shards']}==set(range(workers))
    assert len(done['shards'])==workers
    for s in done['shards']:
        shard=s['shard'];assert s['shards']==workers
        counts=collections.Counter();seen=set()
        with (root/'shards'/f'shard_{shard:02d}_hits.tsv').open() as f:
            for h in csv.DictReader(f,delimiter='\t'):
                h={k:int(v) if k in ('middle','rotation','sign','middle_weight') else v for k,v in h.items()}
                d=ds[h['diff_id']];key=(h['family'],h['middle'])
                assert d['order']%workers==shard and d['family']==h['family']
                assert 0<=h['rotation']<d['orbit'] and h['middle'] in (1,2)
                assert selected(d,h['middle']) and h['sign'] in (-1,1)
                signature=(h['diff_id'],h['rotation'],h['middle'],h['group'])
                assert signature not in seen;seen.add(signature);counts[key]+=1
                for lid in groups[h['group']]['members']:
                    found[key].append({**h,'linear_id':lid,'weight':d['weight']+ls[lid]['weight']+h['middle_weight']})
        assert len(s['rows'])==4
        for r in s['rows']:
            key=(r['family'],r['middle'])
            subset=[d for d in ds.values() if d['order']%workers==shard and d['family']==r['family'] and selected(d,r['middle'])]
            states=sum(d['orbit'] for d in subset)
            assert r['representatives']==len(subset) and r['states']==states
            assert r['group_evaluations']==states*len(groups) and r['logical_pairs']==states*len(ls)
            assert r['nonzero_groups']==counts[key] and r['nonzero_pairs']==100*counts[key]
            totals[key].update({k:v for k,v in r.items() if k not in ('family','middle')})
    (root/'matches').mkdir(exist_ok=True);(root/'ranked_endpoints').mkdir(exist_ok=True)
    rows=[];jobs={}
    for family in ('D1','D2'):
        for m in (1,2):
            key=family,m;name=f'{family}_M{m}_L2_2to2_pool'
            matches=sorted(found[key],key=lambda h:(h['weight'],ds[h['diff_id']]['order'],h['rotation'],h['linear_id']))
            assert len(matches)==totals[key]['nonzero_pairs']
            endpoints={};best=None
            with gzip.open(root/'matches'/f'{name}.jsonl.gz','wt') as f:
                for rank,h in enumerate(matches,1):
                    record=candidate_record(h,ds[h['diff_id']],ls[h['linear_id']])
                    job={'kind':'full_endpoint','rounds':record['R'],'start_round_index':record['start_round_index'],
                         'input_difference':{'words':record['input_difference']},'output_mask':{'words':record['output_mask']}}
                    eid=endpoint_id(job);job['id']=eid;jobs[eid]=job
                    if rank==1:best={**record,'endpoint_id':eid}
                    if eid not in endpoints:endpoints[eid]={'endpoint_id':eid,'first_match_rank':rank,'candidate':record}
                    f.write(json.dumps({**h,'rank':rank,'endpoint_id':eid},separators=(',',':'))+'\n')
            with gzip.open(root/'ranked_endpoints'/f'{name}.jsonl.gz','wt') as f:
                for e in endpoints.values():f.write(json.dumps(e,separators=(',',':'))+'\n')
            exhaustive=totals[key]['representatives']==family_order[family]
            rows.append({'id':name,'family':family,'middle_rounds':m,'tail':'L2_2to2_pool',
                         'r_d':int(family[-1]),'r_m':m,'r_l':2,'R':int(family[-1])+m+2,
                         'differential_representatives_total':family_order[family],
                         'differential_representatives_evaluated':totals[key]['representatives'],
                         'differential_states_evaluated':totals[key]['states'],
                         'linear_representatives':len(ls),'distinct_join_masks':len(groups),
                         'evaluated_pairs':totals[key]['logical_pairs'],
                         'group_evaluations':totals[key]['group_evaluations'],
                         'nonzero_groups':totals[key]['nonzero_groups'],'nonzero_matches':len(matches),
                         'unique_endpoints':len(endpoints),'exhaustive':exhaustive,
                         'optimal_score_certified':bool(exhaustive and best),
                         'best':best,'score_weight_histogram':dict(sorted(collections.Counter(h['weight'] for h in matches).items())),
                         'matches_file':f'matches/{name}.jsonl.gz','ranked_endpoints_file':f'ranked_endpoints/{name}.jsonl.gz'})
    with gzip.open(root/'endpoint_jobs.jsonl.gz','wt') as f:
        for job in jobs.values():f.write(json.dumps(job,separators=(',',':'))+'\n')
    result={'method':'ordinary DATF middle 1/2; exact reuse by identical join mask; all 16000 L2_2to2 trails',
            'first_L_omitted':True,'enhanced_DATF':False,'workers':workers,
            'middle2_representatives_per_family':limit,'rows':rows,'search_seconds':done['elapsed_seconds'],
            'total_logical_pairs':sum(r['evaluated_pairs'] for r in rows),
            'total_group_evaluations':sum(r['group_evaluations'] for r in rows),
            'total_nonzero_matches':sum(r['nonzero_matches'] for r in rows),'unique_full_endpoint_jobs':len(jobs),
            'pool_metadata':json.loads((root/'pool_metadata.json').read_text()),
            'ranking':'absolute splice score, original differential priority, differential rotation, linear id',
            'validation':{'shard_coverage':'PASS','no_duplicate_group_hits':'PASS','complete_linear_group_partition':'PASS',
                          'per_row_coverage_and_nonzero_counts':'PASS','source_hashes':'PASS'},
            'limits':['Individual splice score is not a full endpoint correlation or a hull sum.',
                      'Two-middle rows cover only the declared differential batch; all linear representatives are included.',
                      'Ordinary DATF zero is not a proof of zero true permutation correlation.']}
    atomic_json(root/'search_results.json',result)
    for r in rows:print(r['id'],r['evaluated_pairs'],r['nonzero_matches'],(r['best'] or {}).get('score'),flush=True)


if __name__=='__main__':main()
