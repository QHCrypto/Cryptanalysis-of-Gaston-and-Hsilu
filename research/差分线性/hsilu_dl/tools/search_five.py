#!/usr/bin/env python3
"""Search 3-chi prefixes attached to fixed 2-chi linear tails, total 5 chi.

The prefix uses RC indices 7,8,9 (not the last three-round constant sequence).
All 31 one-column differences and 64 positions are considered; no fixed-RC
rotation equivalence is assumed. Ordinary DATF zeros remain model results.
"""
import argparse
import concurrent.futures
import json
import time
from fractions import Fraction
from pathlib import Path
from hsilu import *
from datf import BudgetExceeded, add

def initialize(path):
    global POOL, HEADS, LAST
    POOL=json.loads(Path(path).read_text());HEADS=POOL['heads'];LAST={}
    for b in POOL['bridges']:
        weight=min(x['weight'] for x in b['last_chi_choices'])
        LAST[b['id']]=[x for x in b['last_chi_choices'] if x['weight']==weight]

def work(job):
    col,pattern=job;delta=embed(col,pattern);e,a,d=construct(delta,3,start=7)
    candidates=[];limited=0;nonzero=0
    for head in HEADS:
        mask=tuple(int(w,16) for w in head['join_mask']['words'])
        target=add(*(d[i][j] for i in range(5) for j in range(64) if mask[i]>>j&1))
        e.max_bias_nodes=e.stats['bias_nodes']+20000
        try:rho=e.estimate(target)
        except BudgetExceeded:limited+=1;continue
        if rho:
            nonzero+=1
            weight=head['first_chi']['weight']+LAST[head['bridge_id']][0]['weight']
            candidates.append((abs(rho)/2**(2*weight),head['id'],str(rho),weight))
    candidates.sort(reverse=True)
    return {'input_col':col,'input_pattern':pattern,'nonzero':nonzero,'resource_limited':limited,
            'top_heads':[{'head_id':hid,'prefix_datf':rho,'linear_weight':weight,'abs_splice_score':str(score)}
                         for score,hid,rho,weight in candidates[:4]]}

def main():
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,default=4)
    p.add_argument('--patterns',default='all');a=p.parse_args();out=a.output
    target=out/'five_prefixes.jsonl'
    if target.exists():raise FileExistsError(target)
    patterns=list(range(1,32)) if a.patterns=='all' else [int(v) for v in a.patterns.split(',')]
    jobs=[(col,pattern) for pattern in patterns for col in range(64)];records=[];began=time.monotonic()
    with target.open('x') as stream,concurrent.futures.ProcessPoolExecutor(max_workers=a.workers,initializer=initialize,initargs=(str(out/'linear_pool.json'),)) as executor:
        for row in executor.map(work,jobs,chunksize=1):
            stream.write(json.dumps(row)+'\n');stream.flush();records.append(row)
            if len(records)%32==0:print('FIVE',len(records),'/',len(jobs),'nonzero',sum(r['nonzero'] for r in records),'seconds',round(time.monotonic()-began,1),flush=True)
    initialize(out/'linear_pool.json');heads={h['id']:h for h in HEADS};endpoints={}
    for row in records:
        for h in row['top_heads']:
            head=heads[h['head_id']]
            for last in LAST[head['bridge_id']]:
                delta=embed(row['input_col'],row['input_pattern']);key=(delta,tuple(last['output']['words']))
                score=Fraction(h['prefix_datf'])/2**(2*h['linear_weight'])
                if key not in endpoints or abs(score)>abs(Fraction(endpoints[key]['splice_score'])):
                    endpoints[key]={'rounds':5,'input_difference':encode(delta),'core_output_mask':last['output'],
                                    'prefix_datf':h['prefix_datf'],'splice_score':str(score),
                                    'linear_weight':h['linear_weight'],'head':head,'last_chi':last}
    ranked=sorted(endpoints.values(),key=lambda r:(-abs(Fraction(r['splice_score'])),r['input_difference']['words'],r['core_output_mask']['words']))
    invt=inverse_function(transpose)
    for i,r in enumerate(ranked):
        r['id']=f'R5_pool_{i:06d}'
        if i<128:r['full_output_mask']=encode(invt(tuple(int(w,16) for w in r['core_output_mask']['words'])))
    document={'scope':'single active input column, all 31 patterns × 64 positions, 384 fixed linear-tail heads',
              'prefix_chi_round_indices':[7,8,9],'prefix_rc_after_first':['78','69'],
              'input_candidates':len(records),'head_evaluations':len(records)*len(HEADS),
              'nonzero_prefix_DATF':sum(r['nonzero'] for r in records),
              'resource_limited':sum(r['resource_limited'] for r in records),
              'candidate_retention':'top four heads per input; minimum-weight last-chi completions',
              'zero_is_not_proof':True,'elapsed_seconds':time.monotonic()-began,'top':ranked[:128]}
    (out/'five_candidates.json').write_text(json.dumps(document,indent=2)+'\n')
    selected=[];seen={}
    for r in ranked:
        key=tuple(r['input_difference']['words'])
        if seen.get(key,0)>=2:continue
        selected.append(r);seen[key]=seen.get(key,0)+1
        if len(selected)==16:break
    for r in selected:r['full_output_mask']=encode(invt(tuple(int(w,16) for w in r['core_output_mask']['words'])))
    (out/'five_selection.json').write_text(json.dumps(selected,indent=2)+'\n')
    with (out/'five_jobs.txt').open('x') as f:
        for r in selected:f.write(f'{r["id"]} 5 '+ ';'.join(r['input_difference']['words'])+' '+';'.join(r['core_output_mask']['words'])+'\n')
    print('DONE FIVE',len(ranked),'unique candidates, selected',len(selected),flush=True)

if __name__=='__main__':main()
