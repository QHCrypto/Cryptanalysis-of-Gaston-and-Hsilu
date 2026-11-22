#!/usr/bin/env python3
import argparse
import collections
import json
import random
import subprocess
import sys
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from datf_pool_build import read_pool,state,weight
from datf_k06 import construct,parity_polynomial
from fixed_tail_search import S,RC,lat,column,rotate,linear
from six_six_prefix_search import ddt_table


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--results',type=Path,default=ROOT/'results/datf_pools/search_final')
    parser.add_argument('--pools',type=Path,default=ROOT/'results/datf_pools/pools')
    args=parser.parse_args()
    pooldir=args.pools
    meta=json.loads((pooldir/'pools.json').read_text())
    data={k:read_pool(pooldir/v['file']) for k,v in meta['pools'].items()}
    assert sum(d['orbit_size'] for d in data['differential_one'])==812800
    assert sum(d['orbit_size'] for d in data['linear_one'])==814080
    assert len(data['differential_pq_two'])==11*4096
    ddt=ddt_table();checked=0
    for kind,key in [('d','differential_one'),('l','linear_one')]:
        pairs=set()
        for d in data[key]:
            a,b=state(d['input']),state(d['output'])
            assert weight(a if kind=='d' else b)<=2
            assert (a,b) not in pairs;pairs.add((a,b))
            score=Fraction(1)
            for j in range(64):score*=Fraction(ddt[column(a,j)][column(b,j)] if kind=='d' else lat(column(a,j),column(b,j)),32)
            assert score==Fraction(d['sign'],1<<d['weight'])
            checked+=1
    rng=random.Random(20261017);jobs=[];expected=[]
    def job(delta,mask,r,start):
        eng,_,diff,_=construct(delta,r,start)
        expected.append(eng.estimate(parity_polynomial(diff,mask)))
        jobs.append(str(r)+' '+' '.join(f'{w:x}' for w in (*delta,*mask,RC[start+1]))+'\n')
    for i in range(96):
        delta=[0]*5;mask=[0]*5
        for _ in range(1+i%6):delta[rng.randrange(5)]^=1<<rng.randrange(64)
        for _ in range(1+i%8):mask[rng.randrange(5)]^=1<<rng.randrange(64)
        job(tuple(delta),tuple(mask),1 if i<24 else 2,rng.randrange(10))
    result=json.loads((args.results/'results.json').read_text()) if (args.results/'results.json').exists() else None
    winners=[]
    if result:
        for s in result['summaries']:
            if s['best']:
                b=s['best'];winners.append(b)
                delta,mask=state(b['middle_input_difference']),state(b['middle_output_mask'])
                job(delta,mask,b['r_m'],12-b['r_l']-b['r_m'])
                assert expected[-1]==Fraction(int(b['middle_DATF_correlation']['numerator']),int(b['middle_DATF_correlation']['denominator']))
                assert Fraction(int(b['score']['numerator']),int(b['score']['denominator']))==expected[-1]*Fraction(int(b['differential_probability']['numerator']),int(b['differential_probability']['denominator']))*Fraction(int(b['linear_squared_correlation']['numerator']),int(b['linear_squared_correlation']['denominator']))
        old=json.loads((ROOT/'results/pq_kernel_six_all/search.json').read_text())
        counts={r['id']:r['feasible_trails'] for r in old['reports']}
        for s in result['summaries']:
            if s['family']=='D2' and s['middle_rounds']==1 and s['tail']!='L1_pool':
                assert s['exhaustive'] and s['nonzero_evaluated_pairs']==counts[s['tail']]
    native=subprocess.check_output([str(ROOT/'build/datf_middle_cli')],input=''.join(jobs),text=True)
    assert len(native.splitlines())==len(expected)
    for row,want in zip(native.splitlines(),expected):
        sign,w,_=row.split();assert Fraction(int(sign),1<<int(w))==want
    # Check the two-chi ordinary model's phase formula at every remaining RC,
    # including simultaneous rotations; this is not a true-correlation test.
    rotation_jobs=[];groups=[]
    for b in [x for x in winners if x['r_m']==2][:8]:
        delta,mask=state(b['middle_input_difference']),state(b['middle_output_mask'])
        group=[]
        for shift in range(64):
            for rc in (0,RC[10],RC[8]):
                group.append(len(rotation_jobs))
                rotation_jobs.append('2 '+' '.join(f'{w:x}' for w in (*rotate(delta,shift),*rotate(mask,shift),rc))+'\n')
        groups.append(group)
    rotation_lines=subprocess.check_output([str(ROOT/'build/datf_middle_cli')],input=''.join(rotation_jobs),text=True).splitlines()
    rotation_values=[Fraction(int(row.split()[0]),1<<int(row.split()[1])) for row in rotation_lines]
    for group in groups:assert len({abs(rotation_values[i]) for i in group})==1
    validation={'status':'PASS','single_round_pool_trails_checked':checked,
                'compiled_vs_original_python_evaluations':len(jobs),
                'nonzero_oracle_cases':sum(bool(x) for x in expected),
                'simultaneous_rotation_constant_checks':len(rotation_jobs),
                'old_PQ_one_middle_match_counts_reproduced':11 if result else 0,
                'scope':'factor arithmetic, pool coverage, boundaries and model invariance; no true full-correlation claim'}
    if result:(args.results/'validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation))


if __name__=='__main__':main()
