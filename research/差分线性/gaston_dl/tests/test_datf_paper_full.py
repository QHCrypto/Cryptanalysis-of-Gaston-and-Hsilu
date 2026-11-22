#!/usr/bin/env python3
"""Audit full matching coverage and every required ranked endpoint result."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from datf_paper_merge import endpoint_id
from datf_paper_reassess import read_jsonl,row_state
from datf_paper_parallel import atomic_json
from datf_pool_build import state,read_pool
from datf_k06 import construct,parity_polynomial
from fixed_tail_search import linear,rotate,column,S,RC


def payload_result(path):
    # Read the leading result object, without materializing the much larger DAG.
    decoder=json.JSONDecoder();text=''
    with gzip.open(path,'rt') as f:
        while True:
            chunk=f.read(65536)
            if not chunk:raise AssertionError('incomplete endpoint artifact: '+str(path))
            text+=chunk
            assert text.startswith('{"result":')
            try:return decoder.raw_decode(text[len('{"result":'):])[0]
            except json.JSONDecodeError:pass


def main():
    p=argparse.ArgumentParser();p.add_argument('--results',type=Path,default=ROOT/'results/datf_pools/paper_full24')
    p.add_argument('--search-only',action='store_true');args=p.parse_args();root=args.results
    search=json.loads((root/'search_results.json').read_text());config=json.loads((root/'configuration.json').read_text())
    completed=json.loads((root/'search_completed.json').read_text())
    assert config['workers']==24 and len(completed['shards'])==24
    assert {s['shard'] for s in completed['shards']}==set(range(24))
    ds={}
    for i,line in enumerate(Path(config['sources'][0]['path']).read_text().splitlines()):
        a=line.split();ds[a[0]]={'family':a[1],'w':int(a[2]),'orbit':int(a[3]),'input':state(a[4:9]),'output':state(a[9:14]),'order':i}
    meta=search['pool_metadata'];tails={t['id']:t for t in meta['paper_tails']}
    pooldir=ROOT/'results/datf_pools/pools'
    for pool_key in ('differential_one','differential_pq_two'):
        desc=meta['pools'][pool_key];path=pooldir/desc['file']
        assert hashlib.sha256(path.read_bytes()).hexdigest()==desc['sha256']
        for original in read_pool(path):
            d=ds[original['id']]
            assert d['input']==state(original['input']) and d['output']==state(original['output'])
            assert d['w']==original['weight'] and d['orbit']==original['orbit_size']
    expected={'D1':812800,'D2':2883584}
    for family in expected:
        assert sum(d['orbit'] for d in ds.values() if d['family']==family)==expected[family]
        assert sum(s[family+'_states'] for s in completed['shards'])==expected[family]
    assert search['total_logical_pairs']==2*11*sum(expected.values())==81320448
    jobs={j['id']:j for j in read_jsonl(root/'endpoint_jobs.jsonl.gz')}
    assert len(jobs)==search['unique_full_endpoint_jobs']
    for eid,j in jobs.items():assert eid==endpoint_id(j) and j['start_round_index']==12-j['rounds']
    dlct=[[Fraction(sum((-1)**((mask&(S[x]^S[x^delta])).bit_count()%2) for x in range(32)),32)
           for mask in range(32)] for delta in range(32)]
    groups={};middle2=[];one_checks=0;all_count=0
    assert len(search['rows'])==44
    for r in search['rows']:
        assert r['exhaustive'] and r['evaluated_pairs']==expected[r['family']]
        matches=read_jsonl(root/r['matches_file']);entries=read_jsonl(root/r['ranked_endpoints_file']);groups[r['id']]=entries
        assert len(matches)==r['nonzero_matches'] and len(entries)==r['unique_endpoints']
        assert matches==sorted(matches,key=lambda h:(h['weight'],ds[h['diff_id']]['order'],h['rotation']))
        assert len({(h['diff_id'],h['rotation']) for h in matches})==len(matches)
        firsts={}
        tail=tails[r['tail']];gamma=state(tail['join_mask'])
        for rank,h in enumerate(matches,1):
            assert h['rank']==rank and h['family']==r['family'] and h['middle']==r['r_m'] and h['tail']==r['tail']
            d=ds[h['diff_id']];assert 0<=h['rotation']<d['orbit']
            delta=rotate(d['input'],h['rotation']);job=jobs[h['endpoint_id']]
            assert state(job['input_difference']['words'])==delta and job['output_mask']['words']==tail['output']
            assert job['rounds']==r['R']
            assert h['weight']==d['w']+tail['weight_squared']+h['middle_weight']
            firsts.setdefault(h['endpoint_id'],rank)
            dm=linear(rotate(d['output'],h['rotation']))
            if r['r_m']==1:
                c=Fraction(1)
                for j in range(64):c*=dlct[column(dm,j)][column(gamma,j)]
                assert c==Fraction(h['sign'],1<<h['middle_weight']);one_checks+=1
            else:middle2.append((dm,gamma,r['r_l'],h['sign'],h['middle_weight']))
        assert [(e['endpoint_id'],e['first_match_rank']) for e in entries]==list(firsts.items())
        all_count+=len(matches)
    assert all_count==search['total_nonzero_matches']
    rng=random.Random(20261016101)
    sampled=rng.sample(middle2,min(64,len(middle2)))
    for dm,gamma,rl,sign,w in sampled:
        engine,_,diff,_=construct(dm,2,10-rl)
        assert engine.estimate(parity_polynomial(diff,gamma))==Fraction(sign,1<<w)
    validation={'status':'PASS','workers':24,'rows':44,'logical_pairs':search['total_logical_pairs'],
                'all_nonzero_matches_checked':all_count,'independent_single_chi_DLCT_checks':one_checks,
                'original_ordinary_DATF_middle2_checks':len(sampled),'unique_full_endpoint_jobs':len(jobs)}
    if not args.search_only:
        review=json.loads((root/'reassessment/results.json').read_text())
        receipt=json.loads((root/'reassessment/completed.json').read_text())
        assert review['complete'] and receipt['all_rows_terminal'] and receipt['workers']==24
        assert len(review['worker_pids_observed'])>=24
        known={};audits=collections.Counter()
        for f in (root/'reassessment/jobs').glob('*.summary.json'):
            c=json.loads(f.read_text());eid=c['id'];assert eid in jobs
            full=payload_result(Path(c['artifact']))
            assert {k:full[k] for k in jobs[eid]}==jobs[eid]
            for k in ('status','correlation','bias','worker_pid','max_bias_nodes'):assert c[k]==full[k]
            if c['status']=='completed':
                assert full['dag_verified_samples']>=4
                assert full['factor_audit']['status']=='PASS'
                assert full['factor_audit']['correlation']==full['correlation']
                rho=Fraction(int(full['correlation']['numerator']),int(full['correlation']['denominator']))
                assert abs(rho)<=1
                assert Fraction(int(full['bias']['numerator']),int(full['bias']['denominator']))==rho/2
                audits['PASS']+=1
            known[eid]=c
        for r in review['rows']:
            checked=row_state(r,groups[r['id']],known)
            assert checked['status']==r['status']!='running'
            assert checked['first_nonzero']==r['first_nonzero']
            assert checked['verified_prefix_endpoints']==r['verified_prefix_endpoints']
            if checked['status']=='found_nonzero':assert checked['full_DATF']==r['full_DATF']
            # Independent termination audit, without relying on the scheduler's
            # row_state routine: walk the actual ranked unique endpoint list.
            if r['status']=='no_nonzero_matches':
                assert not groups[r['id']] and r['nonzero_matches']==0
                continue
            winner=r['first_nonzero'];required=0;reached=False
            for e in groups[r['id']]:
                saved=known[e['endpoint_id']]
                assert saved['status']=='completed'
                numerator=int(saved['correlation']['numerator']);required+=1
                if winner and e['endpoint_id']==winner['endpoint_id']:
                    assert numerator!=0 and e['first_match_rank']==winner['first_match_rank']
                    reached=True;break
                assert numerator==0
            if winner:assert reached
            else:assert required==r['unique_endpoints']
        validation.update(full_endpoint_artifacts_checked=len(known),independent_factor_audits=dict(audits),
                          every_required_rank_prefix_verified=True,
                          row_status_counts=dict(collections.Counter(r['status'] for r in review['rows'])))
    atomic_json(root/('validation_search.json' if args.search_only else 'validation.json'),validation)
    print(json.dumps(validation),flush=True)


if __name__=='__main__':main()
