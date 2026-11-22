#!/usr/bin/env python3
"""Exhaust all paper-tail matches in independent native process shards."""
import argparse
import concurrent.futures
import hashlib
import json
import os
import subprocess
import time
from pathlib import Path


def atomic_json(path, doc):
    temp=path.with_suffix(path.suffix+'.partial')
    temp.write_text(json.dumps(doc,ensure_ascii=False,indent=2)+'\n')
    temp.replace(path)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,default=Path('results/datf_pools/search_final'))
    p.add_argument('--workers',type=int,default=24)
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--resume',action='store_true')
    args=p.parse_args()
    if args.workers<1:p.error('workers must be positive')
    args.output.mkdir(parents=True,exist_ok=args.resume)
    (args.output/'shards').mkdir(exist_ok=True)
    sources=[args.inputs/'differentials.txt',args.inputs/'linears.txt']
    config={'workers':args.workers,'scope':'D1,D2 x 11 paper tails x middle 1,2; all nonzero matches retained',
            'sources':[{'path':str(f.resolve()),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in sources],
            'binary_sha256':hashlib.sha256(Path('build/datf_paper_shard').read_bytes()).hexdigest()}
    if args.resume:assert json.loads((args.output/'configuration.json').read_text())==config
    else:atomic_json(args.output/'configuration.json',config)
    began=time.time()
    def run(shard):
        prefix=args.output/'shards'/f'shard_{shard:02d}'
        if Path(str(prefix)+'_summary.json').exists():return {'shard':shard,'cached':True}
        command=['build/datf_paper_shard',*map(str,sources),str(prefix),str(shard),str(args.workers)]
        with Path(str(prefix)+'.log').open('a') as log:
            process=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT)
            atomic_json(Path(str(prefix)+'_process.json'),{'pid':process.pid,'command':command,'started_unix':time.time()})
            code=process.wait()
        if code:raise RuntimeError(f'shard {shard} failed: exit {code}')
        return {'shard':shard,'cached':False}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as ex:
        pending={ex.submit(run,i):i for i in range(args.workers)}
        while pending:
            done,_=concurrent.futures.wait(pending,timeout=25,return_when=concurrent.futures.FIRST_COMPLETED)
            for f in done:print('DONE',f.result(),flush=True);del pending[f]
            states=[]
            for f in sorted((args.output/'shards').glob('*_progress.json')):
                states.append(json.loads(f.read_text()))
            progress={'coordinator_pid':os.getpid(),'workers':args.workers,'live_shards':len(pending),
                      'elapsed_seconds':time.time()-began,'representatives_done':sum(s['representatives_done'] for s in states),
                      'representatives_total':sum(s['representatives_total'] for s in states),
                      'states_done':sum(s['D1_states']+s['D2_states'] for s in states)}
            atomic_json(args.output/'progress.json',progress)
            print('PROGRESS',json.dumps(progress),flush=True)
    atomic_json(args.output/'search_completed.json',{'workers':args.workers,'elapsed_seconds':time.time()-began,
                'shards':[json.loads(f.read_text()) for f in sorted((args.output/'shards').glob('*_summary.json'))]})


if __name__=='__main__':main()
