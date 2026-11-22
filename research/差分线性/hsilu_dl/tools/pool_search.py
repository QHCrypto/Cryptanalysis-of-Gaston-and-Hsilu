#!/usr/bin/env python3
"""Hsilu linear two-round tails and exact-two-chi prefix matching."""
import argparse
import csv
import itertools
import json
import math
from fractions import Fraction
from pathlib import Path
from hsilu import *

def chi5(x):return sum((((x>>i)&1)^((1^((x>>((i+1)%5))&1))&((x>>((i+2)%5))&1)))<<i for i in range(5))
LAT=[[sum(1 if ((a&x).bit_count()+(b&chi5(x)).bit_count())%2==0 else -1 for x in range(32)) for b in range(32)] for a in range(32)]
def local_score(a,b):
    result=Fraction(1)
    for j in range(64):
        u=sum(((a[i]>>j)&1)<<i for i in range(5));v=sum(((b[i]>>j)&1)<<i for i in range(5))
        result*=Fraction(LAT[u][v],32)
    assert result
    return {'sign':1 if result>0 else -1,'weight':int(-math.log2(abs(result))),
            'correlation':str(result),'input':encode(a),'output':encode(b)}
def decode(a):return tuple(int(w,16) for w in a['words'])
def completions(fixed,input_free):
    cols=[j for j in range(64) if any(w>>j&1 for w in fixed)]
    patterns=[sum(((fixed[i]>>j)&1)<<i for i in range(5)) for j in cols]
    choices=[[v for v in range(1,32) if (LAT[v][u] if input_free else LAT[u][v])] for u in patterns]
    for values in itertools.product(*choices):
        yield tuple(sum(((v>>i)&1)<<j for j,v in zip(cols,values)) for i in range(5))
def build(directory):
    bridges=json.loads((directory/'linear_bridges.json').read_text())['bridges']
    heads=[];tails=[]
    for index,bridge in enumerate(bridges):
        b1,a2=decode(bridge['alpha']),decode(bridge['beta'])
        assert transpose(a2)==b1
        last=[local_score(a2,b2) for b2 in completions(a2,False)]
        tails.append({'id':index,'first_output':encode(b1),'second_input':encode(a2),'last_chi_choices':last})
        for a1 in completions(b1,True):
            local=local_score(a1,b1)
            heads.append({'id':f'H{len(heads):04d}','bridge_id':index,'join_mask':encode(transpose(a1)),
                          'first_chi':local})
    with (directory/'heads.txt').open('x') as f:
        for h in heads:f.write(f'{h["id"]} {h["first_chi"]["weight"]} '+ ';'.join(h['join_mask']['words'])+'\n')
    (directory/'linear_pool.json').write_text(json.dumps({'heads':heads,'bridges':tails},indent=2)+'\n')
    print('heads',len(heads),'complete tails',sum(sum(h['bridge_id']==i for h in heads)*len(t['last_chi_choices']) for i,t in enumerate(tails)))
def rank(directory):
    pool=json.loads((directory/'linear_pool.json').read_text());heads={h['id']:h for h in pool['heads']}
    best={};nonzero=0
    with (directory/'prefix_matches.csv').open() as f:
        for row in csv.DictReader(f):
            nonzero+=1;head=heads[row['head_id']];delta=embed(int(row['input_col']),int(row['input_pattern']))
            prefix=Fraction(int(row['prefix_numerator']),32)
            for last in pool['bridges'][head['bridge_id']]['last_chi_choices']:
                weight=head['first_chi']['weight']+last['weight'];score=prefix/Fraction(1<< (2*weight))
                mask=decode(last['output']);key=(delta,mask)
                if key not in best or abs(score)>abs(Fraction(best[key]['splice_score'])):
                    best[key]={'rounds':4,'input_difference':encode(delta),'core_output_mask':last['output'],
                               'exact_prefix_correlation':str(prefix),'linear_weight':weight,
                               'splice_score':str(score),'head':head,'last_chi':last}
    ranked=sorted(best.values(),key=lambda r:(-abs(Fraction(r['splice_score'])),r['input_difference']['words'],r['core_output_mask']['words']))
    invt=inverse_function(transpose)
    for i,r in enumerate(ranked):
        r['id']=f'R4_pool_{i:06d}'
        if i<128:r['full_output_mask']=encode(invt(decode(r['core_output_mask'])))
    summary={'prefix_nonzero_matches':nonzero,'unique_endpoints':len(ranked),
             'scope':'single-active-column first difference, exact two-chi prefix, enumerated 3->2 linear bridges',
             'warning':'prefix correlation is exact; multiplying by a single tail squared correlation is a splice heuristic, not a full correlation',
             'top':ranked[:128]}
    (directory/'pool_candidates.json').write_text(json.dumps(summary,indent=2)+'\n')
    # Diverse initial selection: retain at most two per input difference and
    # at most one per output mask, to avoid spending all samples on one family.
    selected=[];inputs={};outputs=set()
    for r in ranked:
        a=tuple(r['input_difference']['words']);b=tuple(r['core_output_mask']['words'])
        if inputs.get(a,0)>=2 or b in outputs:continue
        selected.append(r);inputs[a]=inputs.get(a,0)+1;outputs.add(b)
        if len(selected)==16:break
    for r in selected:
        r['full_output_mask']=encode(invt(decode(r['core_output_mask'])))
    (directory/'pool_selection.json').write_text(json.dumps(selected,indent=2)+'\n')
    with (directory/'pool_jobs.txt').open('x') as f:
        for r in selected:f.write(f'{r["id"]} 4 '+ ';'.join(r['input_difference']['words'])+' '+';'.join(r['core_output_mask']['words'])+'\n')
    print('nonzero prefixes',nonzero,'unique endpoints',len(ranked),'selected',len(selected))
    print([(r['id'],r['splice_score']) for r in selected[:4]])
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['build','rank']);p.add_argument('--output',type=Path,required=True);a=p.parse_args()
    (build if a.action=='build' else rank)(a.output)
