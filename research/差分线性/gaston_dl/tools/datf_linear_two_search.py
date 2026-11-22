#!/usr/bin/env python3
"""Match all L2 2->2 trails, sharing ordinary DATF by identical join mask.

Middle 1 exhausts both differential pools; middle 2 uses a stated prefix of
each priority-ordered differential pool (0 means all). All linears are covered.
"""
import argparse
import collections
import concurrent.futures
import hashlib
import json
import shutil
import subprocess
import time
from pathlib import Path

from datf_pool_build import read_pool
from datf_paper_parallel import atomic_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,default=Path('results/datf_pools/linear_two_ready'))
    p.add_argument('--pools',type=Path,default=Path('results/datf_pools/pools'))
    p.add_argument('--middle2-representatives',type=int,default=64)
    p.add_argument('--workers',type=int,default=24)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    if args.workers<1 or args.middle2_representatives<0:p.error('invalid workers or batch size')
    args.output.mkdir(parents=True,exist_ok=args.resume)
    (args.output/'shards').mkdir(exist_ok=True)
    pool_path=args.pools/'linear_two_2to2.jsonl.gz'
    records=read_pool(pool_path)
    grouped=collections.defaultdict(list)
    for r in records:grouped[tuple(r['join_mask'])].append(r['id'])
    groups=[{'id':f'G{i:03d}','join_mask':list(mask),'members':sorted(ids)}
            for i,(mask,ids) in enumerate(sorted(grouped.items()))]
    assert len(groups)==160 and len(records)==16000
    assert all(len(g['members'])==100 for g in groups)
    source_paths=[args.inputs/'differentials.txt', pool_path, Path('build/datf_linear_two_shard')]
    config={'workers':args.workers,'middle2_representatives_per_family':args.middle2_representatives,
            'linear_representatives':len(records),'distinct_join_masks':len(groups),
            'first_L_omitted':True,'middle1_scope':'all differential representatives and rotations',
            'middle2_scope':'first N priority-ordered representatives per differential family, all their rotations; N=0 means all',
            'linear_scope':'all 16000 L2_2to2 representatives; exact reuse of the 160 distinct middle output masks',
            'sources':[{'path':str(f.resolve()),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in source_paths]}
    if args.resume:
        assert json.loads((args.output/'configuration.json').read_text())==config
    else:
        atomic_json(args.output/'configuration.json',config)
        atomic_json(args.output/'linear_groups.json',groups)
        atomic_json(args.output/'pool_metadata.json',json.loads((args.pools/'pools.json').read_text()))
        shutil.copyfile(args.inputs/'differentials.txt',args.output/'differentials.txt')
        (args.output/'groups.txt').write_text(''.join(g['id']+' '+str(len(g['members']))+' '+' '.join(g['join_mask'])+'\n' for g in groups))
    began=time.time()
    def run(i):
        prefix=args.output/'shards'/f'shard_{i:02d}'
        if Path(str(prefix)+'_summary.json').exists():return i
        command=['build/datf_linear_two_shard',str(args.output/'differentials.txt'),str(args.output/'groups.txt'),
                 str(prefix),str(i),str(args.workers),str(args.middle2_representatives)]
        with Path(str(prefix)+'.log').open('a') as f:subprocess.run(command,stdout=f,stderr=subprocess.STDOUT,check=True)
        return i
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        pending={ex.submit(run,i) for i in range(args.workers)}
        while pending:
            done,pending=concurrent.futures.wait(pending,timeout=20,return_when=concurrent.futures.FIRST_COMPLETED)
            for f in done:print('DONE',f.result(),flush=True)
            shards=[json.loads(f.read_text()) for f in (args.output/'shards').glob('*_progress.json')]
            progress={'live_shards':len(pending),'elapsed_seconds':time.time()-began,
                      **{k:sum(s[k] for s in shards) for k in ('representatives_done','representatives_total','middle1_states','middle2_states')}}
            atomic_json(args.output/'progress.json',progress)
            print('PROGRESS',json.dumps(progress),flush=True)
    atomic_json(args.output/'search_completed.json',{'elapsed_seconds':time.time()-began,
                'shards':[json.loads(f.read_text()) for f in sorted((args.output/'shards').glob('*_summary.json'))]})


if __name__=='__main__':main()
