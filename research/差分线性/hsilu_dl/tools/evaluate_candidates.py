#!/usr/bin/env python3
import argparse
import json
import time
from pathlib import Path
from hsilu import *
from datf import BudgetExceeded, add
from datf_k06 import value_record

p=argparse.ArgumentParser()
p.add_argument('--source',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args()
if a.output.exists():raise FileExistsError(a.output)
jobs=json.loads(a.source.read_text());results=[];cache={}
for job in jobs:
    R=job['rounds'];delta=tuple(int(w,16) for w in job['input_difference']['words'])
    mask=tuple(int(w,16) for w in job['core_output_mask']['words'])
    start=time.monotonic();key=(R,delta)
    if key not in cache:cache[key]=construct(delta,R)
    e,base,diff=cache[key]
    target=add(*(diff[i][j] for i in range(5) for j in range(64) if mask[i]>>j&1))
    e.max_bias_nodes=e.stats['bias_nodes']+500000;trace=[]
    row={'id':job['id'],'rounds':R,'input_difference':job['input_difference'],
         'core_output_mask':job['core_output_mask'],'target_terms':len(target)}
    try:
        value=e.estimate(target,label=job['id'],trace=trace)
        row.update(status='completed',correlation=value_record(value))
    except BudgetExceeded as exc:row.update(status='resource_limit',correlation=None,reason=str(exc))
    row.update(trace=trace,elapsed_seconds=time.monotonic()-start)
    results.append(row)
    a.output.write_text(json.dumps({'method':'ordinary DATF, unchanged Gaston-project engine; Hsilu symbolic linear layer',
                                    'zero_is_not_proof':True,'results':results},indent=2)+'\n')
    print(row['id'],row['status'],'log2_abs',None if row['correlation'] is None else row['correlation']['log2_abs'],flush=True)
