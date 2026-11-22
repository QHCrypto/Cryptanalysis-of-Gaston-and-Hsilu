#!/usr/bin/env python3
"""Independent middle checks and full-DATF reassessment of all best-score ties."""
import argparse
import collections
import concurrent.futures
import csv
import gzip
import json
import multiprocessing
import random
import time
from fractions import Fraction
from pathlib import Path

from datf_k06 import construct,parity_polynomial,value_record
from datf_paper_parallel import atomic_json
from datf_paper_reassess import evaluate_endpoint
from datf_pool_build import state
from fixed_tail_search import S,linear,rotate


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results',type=Path,required=True)
    p.add_argument('--workers',type=int,default=16)
    args=p.parse_args();root=args.results;began=time.perf_counter()
    result=json.loads((root/'search_results.json').read_text())
    config=json.loads((root/'configuration.json').read_text())
    groups={g['id']:state(g['join_mask']) for g in json.loads((root/'linear_groups.json').read_text())}
    ds={};families=collections.defaultdict(list)
    for line in (root/'differentials.txt').read_text().splitlines():
        t=line.split();ds[t[0]]={'family':t[1],'orbit':int(t[3]),'output':state(t[9:14])};families[t[1]].append(t[0])
    dlct=[[Fraction(sum((-1)**((b&(S[x]^S[x^a])).bit_count()&1) for x in range(32)),32)
           for b in range(32)] for a in range(32)]
    checked=0
    for f in sorted((root/'shards').glob('*_hits.tsv')):
        for h in csv.DictReader(f.open(),delimiter='\t'):
            if int(h['middle'])!=1:continue
            delta=linear(rotate(ds[h['diff_id']]['output'],int(h['rotation'])));gamma=groups[h['group']]
            rho=Fraction(1)
            for j in range(64):
                a=sum(((v>>j)&1)<<i for i,v in enumerate(delta))
                b=sum(((v>>j)&1)<<i for i,v in enumerate(gamma))
                rho*=dlct[a][b]
            assert rho==Fraction(int(h['sign']),1<<int(h['middle_weight']))
            checked+=1
    rng=random.Random(202610162200);zero_checks=0
    # Check sampled zeros with the original Python ordinary-DATF implementation.
    # Full DAG begins at the middle's first chi; r_l=2 makes its start index 8.
    for family in ('D1','D2'):
        chosen=families[family][:config['middle2_representatives_per_family'] or None]
        for _ in range(16):
            did=rng.choice(chosen);d=ds[did];shift=rng.randrange(d['orbit']);gamma=rng.choice(list(groups.values()))
            delta=linear(rotate(d['output'],shift))
            eng,_,diff,_=construct(delta,2,8)
            rho=eng.estimate(parity_polynomial(diff,gamma))
            assert not rho
            zero_checks+=1
    audit={'status':'PASS','all_nonzero_middle1_groups_direct_DLCT_checked':checked,
           'middle2_zero_samples_original_python_DATF':zero_checks,
           'source_and_group_coverage_audit':result['validation'],
           'scope':'Exact one-chi DLCT and ordinary-DATF model implementation checks, not true full permutation correlations.'}
    atomic_json(root/'matching_validation.json',audit)
    folder=root/'endpoint_check';folder.mkdir(exist_ok=True)
    selected=[]
    for row in result['rows']:
        if not row['best']:continue
        w=row['best']['score_weight']
        with gzip.open(root/row['ranked_endpoints_file'],'rt') as f:
            for line in f:
                e=json.loads(line)
                if e['candidate']['score_weight']!=w:break
                selected.append({'row_id':row['id'],**e})
    results={};jobs=[]
    for e in selected:
        eid=e['endpoint_id'];b=e['candidate'];cached=folder/(eid+'.summary.json')
        if cached.exists():results[eid]=json.loads(cached.read_text());continue
        jobs.append({'id':eid,'kind':'full_endpoint','rounds':b['R'],'start_round_index':b['start_round_index'],
                     'input_difference':{'words':b['input_difference']},'output_mask':{'words':b['output_mask']}})
    ctx=multiprocessing.get_context('spawn')
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers,mp_context=ctx) as ex:
        futures=[ex.submit(evaluate_endpoint,j,str(folder),1000000,8) for j in jobs]
        for f in concurrent.futures.as_completed(futures):
            r=f.result();results[r['id']]=r
    entries=[{**e,'result':results[e['endpoint_id']]} for e in selected]
    completed=[e for e in entries if e['result']['status']=='completed']
    nonzero=[e for e in completed if int(e['result']['correlation']['numerator'])]
    def magnitude(e):
        r=e['result']['correlation'];return abs(Fraction(int(r['numerator']),int(r['denominator'])))
    best=max(nonzero,key=magnitude) if nonzero else None
    # Independently recompute the best full endpoint with the original backend.
    original=None
    if best:
        b=best['candidate'];eng,_,diff,_=construct(state(b['input_difference']),b['R'],b['start_round_index'])
        rho=eng.estimate(parity_polynomial(diff,state(b['output_mask'])))
        assert value_record(rho)==best['result']['correlation']
        original={'status':'PASS','endpoint_id':best['endpoint_id'],'correlation':value_record(rho),
                  'backend':'original ordinary DATF with Shannon Bias'}
    document={'scope':'Every unique full endpoint tied for the best individual splice score in each nonempty matched row; not all matched endpoints.',
              'first_L_omitted':True,'selected_endpoints':len(entries),'completed_endpoints':len(completed),
              'nonzero_endpoints':len(nonzero),'model_zero_endpoints':len(completed)-len(nonzero),
              'resource_limited_endpoints':len(entries)-len(completed),'best_full_DATF_within_selection':best,
              'original_backend_control':original,'entries':entries,
              'elapsed_seconds':time.perf_counter()-began}
    atomic_json(folder/'results.json',document)
    print(json.dumps({k:v for k,v in document.items() if k not in ('entries','best_full_DATF_within_selection')},ensure_ascii=False),flush=True)
    if best:print('BEST_FULL',best['candidate'],best['result'],flush=True)


if __name__=='__main__':main()
