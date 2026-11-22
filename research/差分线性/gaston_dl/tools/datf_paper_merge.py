#!/usr/bin/env python3
"""Audit exhaustive shard coverage and retain every nonzero paper-tail match."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
from fractions import Fraction
from pathlib import Path

from datf_k06 import value_record
from datf_pool_build import state,words
from datf_paper_parallel import atomic_json
from fixed_tail_search import rotate,linear,RC


def endpoint_id(job):
    key=[job['rounds'],job['start_round_index'],job['input_difference']['words'],job['output_mask']['words']]
    return hashlib.sha256(json.dumps(key,separators=(',',':')).encode()).hexdigest()


def candidate_record(h,d,l):
    rotation=h['rotation'];rd=d['rounds'];rm=h['middle'];rl=l['rounds'];R=rd+rm+rl
    delta=words(rotate(state(d['input']),rotation));out=rotate(state(d['output']),rotation)
    return {'differential_id':d['id'],'linear_id':l['id'],'differential_rotation':rotation,
            'r_d':rd,'r_m':rm,'r_l':rl,'R':R,'start_round_index':12-R,
            'input_difference':delta,'output_mask':l['output'],
            'differential_output':words(out),'middle_input_difference':words(linear(out)),
            'middle_output_mask':l['gamma'],
            'differential_probability':value_record(Fraction(1,1<<d['weight'])),
            'linear_squared_correlation':value_record(Fraction(1,1<<l['weight'])),
            'middle_DATF_correlation':value_record(Fraction(h['sign'],1<<h['middle_weight'])),
            'score':value_record(Fraction(h['sign'],1<<h['weight'])),'score_weight':h['weight'],
            'middle_second_constant':f'{RC[11-rl]:02x}' if rm==2 else None}


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True)
    p.add_argument('--pools',type=Path,default=Path('results/datf_pools/pools'))
    args=p.parse_args();root=args.results
    config=json.loads((root/'configuration.json').read_text())
    completed=json.loads((root/'search_completed.json').read_text())
    assert len(completed['shards'])==config['workers']==24
    for f in config['sources']:assert hashlib.sha256(Path(f['path']).read_bytes()).hexdigest()==f['sha256']
    ds={};ls={}
    for i,line in enumerate(Path(config['sources'][0]['path']).read_text().splitlines()):
        r=line.split();ds[r[0]]={'id':r[0],'family':r[1],'rounds':int(r[1][-1]),'weight':int(r[2]),
                               'orbit':int(r[3]),'input':r[4:9],'output':r[9:14],'order':i}
    for line in Path(config['sources'][1]['path']).read_text().splitlines():
        r=line.split()
        if r[1]=='paper':ls[r[0]]={'id':r[0],'rounds':int(r[2]),'weight':int(r[3]),'gamma':r[4:9],'output':r[9:14]}
    assert len(ds)==57766 and len(ls)==11
    totals=collections.Counter();nonzeros=collections.Counter();found=collections.defaultdict(list)
    shard_ids=set()
    for summary in completed['shards']:
        shard=summary['shard'];assert shard not in shard_ids;shard_ids.add(shard)
        expected=[d for d in ds.values() if d['order']%config['workers']==shard]
        for family in ('D1','D2'):
            subset=[d for d in expected if d['family']==family]
            assert summary[family+'_representatives']==len(subset)
            assert summary[family+'_states']==sum(d['orbit'] for d in subset)
        counts=collections.Counter()
        with (root/'shards'/f'shard_{shard:02d}_hits.tsv').open() as f:
            for row in csv.DictReader(f,delimiter='\t'):
                h={k:(int(v) if k in ('middle','rotation','weight','sign','middle_weight') else v) for k,v in row.items()}
                d=ds[h['diff_id']];assert d['order']%config['workers']==shard
                assert h['family']==d['family'] and 0<=h['rotation']<d['orbit'] and h['tail'] in ls
                assert h['sign'] in (-1,1) and h['middle'] in (1,2)
                assert h['weight']==d['weight']+ls[h['tail']]['weight']+h['middle_weight']
                key=(h['family'],h['middle'],h['tail']);counts[key]+=1;found[key].append(h)
        assert len(summary['rows'])==44
        for r in summary['rows']:
            key=(r['family'],r['middle'],r['tail'])
            assert r['evaluated']==summary[r['family']+'_states']
            assert counts[key]==r['nonzero'];totals[key]+=r['evaluated'];nonzeros[key]+=r['nonzero']
    assert shard_ids==set(range(24))
    (root/'matches').mkdir(exist_ok=True);(root/'ranked_endpoints').mkdir(exist_ok=True)
    all_endpoints={};rows=[]
    for family in ('D1','D2'):
        expected_states=sum(d['orbit'] for d in ds.values() if d['family']==family)
        assert expected_states==({'D1':812800,'D2':2883584}[family])
        for middle in (1,2):
            for tail,l in ls.items():
                key=(family,middle,tail);name=f'{family}_M{middle}_{tail}'
                assert totals[key]==expected_states
                matches=sorted(found[key],key=lambda h:(h['weight'],ds[h['diff_id']]['order'],h['rotation']))
                assert len(matches)==nonzeros[key]
                assert len({(h['diff_id'],h['rotation']) for h in matches})==len(matches)
                endpoints={};best=None
                with gzip.open(root/'matches'/f'{name}.jsonl.gz','wt') as f:
                    for rank,h in enumerate(matches,1):
                        d=ds[h['diff_id']];record=candidate_record(h,d,l)
                        job={'kind':'full_endpoint','rounds':record['R'],'start_round_index':record['start_round_index'],
                             'input_difference':{'words':record['input_difference']},'output_mask':{'words':record['output_mask']}}
                        eid=endpoint_id(job);job['id']=eid
                        if eid in all_endpoints:assert all_endpoints[eid]==job
                        else:all_endpoints[eid]=job
                        if rank==1:best={**record,'endpoint_id':eid}
                        if eid not in endpoints:endpoints[eid]={'endpoint_id':eid,'first_match_rank':rank,'candidate':record}
                        f.write(json.dumps({**h,'rank':rank,'endpoint_id':eid},separators=(',',':'))+'\n')
                with gzip.open(root/'ranked_endpoints'/f'{name}.jsonl.gz','wt') as f:
                    for e in endpoints.values():f.write(json.dumps(e,separators=(',',':'))+'\n')
                rows.append({'id':name,'family':family,'middle_rounds':middle,'tail':tail,
                             'r_d':int(family[-1]),'r_m':middle,'r_l':l['rounds'],'R':int(family[-1])+middle+l['rounds'],
                             'evaluated_pairs':totals[key],'nonzero_matches':len(matches),'unique_endpoints':len(endpoints),
                             'exhaustive':True,'best':best,'matches_file':f'matches/{name}.jsonl.gz',
                             'ranked_endpoints_file':f'ranked_endpoints/{name}.jsonl.gz'})
    with gzip.open(root/'endpoint_jobs.jsonl.gz','wt') as f:
        for job in all_endpoints.values():f.write(json.dumps(job,separators=(',',':'))+'\n')
    old=json.loads(Path('results/datf_pools/search_final/results.json').read_text())
    old_lookup={(s['family'],s['middle_rounds'],s['tail']):s for s in old['summaries'] if s['tail']!='L1_pool'}
    for r in rows:
        s=old_lookup[r['family'],r['middle_rounds'],r['tail']]
        if r['middle_rounds']==1:
            assert r['nonzero_matches']==s['nonzero_evaluated_pairs']
            assert (r['best'] or {}).get('score')==(s['best'] or {}).get('score')
        else:
            assert r['nonzero_matches']>=s['nonzero_evaluated_pairs']
            if s['best']:assert r['best']['score_weight']<=s['best']['score_weight']
    result={'method':'ordinary DATF middle 1/2, exhaustive paper-tail matching','workers':24,
            'first_L_omitted':True,'single_round_cap':2,'all_nonzero_matches_retained':True,
            'ranking':'absolute individual splice score, then original differential order, then rotation',
            'pool_metadata':json.loads((args.pools/'pools.json').read_text()),'rows':rows,
            'unique_full_endpoint_jobs':len(all_endpoints),'search_seconds':completed['elapsed_seconds'],
            'total_logical_pairs':sum(totals.values()),'total_nonzero_matches':sum(nonzeros.values()),
            'validation':{'shard_coverage':'PASS','per_row_full_pool_counts':'PASS','no_duplicate_matches':'PASS',
                          'all_22_old_middle1_rows_reproduced':'PASS','all_22_middle2_rows_superset_of_old_batch':'PASS'}}
    atomic_json(root/'search_results.json',result)
    print('MERGED',len(rows),'rows',sum(nonzeros.values()),'nonzero matches',len(all_endpoints),'unique full endpoints',flush=True)


if __name__=='__main__':main()
