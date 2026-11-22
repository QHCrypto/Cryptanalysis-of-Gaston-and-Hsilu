#!/usr/bin/env python3
"""24-process ranked full-endpoint reassessment until each row finds nonzero.

Rows retain every match, deduplicate identical full endpoints, and require every
earlier unique endpoint to have completed as zero before accepting a winner.
Speculative parallel work is allowed, but completion order never changes rank.
"""
import argparse
import concurrent.futures
import gzip
import hashlib
import json
import multiprocessing
import os
import time
from pathlib import Path

from datf import BudgetExceeded,indices,support
from datf_k06 import construct,decode,parity_polynomial,value_record,verify_dag
from datf_fixed_tail import audit_trace
from datf_quadratic import QuadraticOrdinaryDATF
from datf_paper_parallel import atomic_json
from fixed_tail_search import RC


def read_jsonl(path):
    with gzip.open(path,'rt') as f:return [json.loads(line) for line in f]


def evaluate_endpoint(job,folder,node_limit,samples):
    began=time.time();pid=os.getpid();folder=Path(folder)
    atomic_json(folder/(job['id']+'.running.json'),{'pid':pid,'started_unix':began,'node_limit':node_limit})
    engine,base,diff,stages=construct(decode(job['input_difference']),job['rounds'],job['start_round_index'],
                                    node_limit,engine_class=QuadraticOrdinaryDATF)
    verify_dag(engine,base,diff,decode(job['input_difference']),job['rounds'],job['start_round_index'],samples=samples)
    target=parity_polynomial(diff,decode(job['output_mask']));trace=[]
    result={**job,'worker_pid':pid,'max_bias_nodes':node_limit,'dag_verified_samples':samples,'stages':stages,
            'target_terms':len(target),'target_variables':support(target).bit_count(),
            'remaining_constants':[f'0x{RC[job["start_round_index"]+r]:02x}' for r in range(1,job['rounds'])],
            'backend':'ordinary Algorithm 1 with exact independent quadratic Bias; unchanged assumptions'}
    try:
        rho=engine.estimate(target,label=job['id'],trace=trace)
        result.update(status='completed',correlation=value_record(rho),bias=value_record(rho/2))
        result['factor_audit']=audit_trace(engine,target,trace)
    except BudgetExceeded as exc:
        result.update(status='resource_limit',correlation=None,bias=None,reason=str(exc))
    result.update(trace=trace,engine_statistics=dict(engine.stats),elapsed_seconds=time.time()-began)
    payload={'result':result,'target_polynomial':[list(indices(t)) for t in sorted(target)],'dag':engine.serialized_dag()}
    path=folder/(job['id']+'.json.gz');temp=folder/(job['id']+'.json.gz.partial')
    with gzip.open(temp,'wt') as f:json.dump(payload,f,separators=(',',':'));f.write('\n')
    temp.replace(path)
    compact={k:result[k] for k in ('id','status','correlation','bias','worker_pid','max_bias_nodes','elapsed_seconds')}
    compact['factor_audit_status']=result.get('factor_audit',{}).get('status')
    compact['artifact']=str(path.resolve())
    atomic_json(folder/(job['id']+'.summary.json'),compact)
    return compact


def row_state(row,entries,known):
    prefix=0
    for e in entries:
        r=known.get(e['endpoint_id'])
        if not r or r['status']!='completed':break
        prefix+=1
        if int(r['correlation']['numerator']):
            return {'status':'found_nonzero','verified_prefix_endpoints':prefix,'first_nonzero':e,
                    'full_DATF':r['correlation'],'bias':r['bias'],'result_artifact':r['artifact']}
    if prefix==len(entries):return {'status':'no_nonzero_matches' if not entries else 'all_endpoints_model_zero',
                                   'verified_prefix_endpoints':prefix,'first_nonzero':None}
    return {'status':'running','verified_prefix_endpoints':prefix,'first_nonzero':None}


def next_endpoint(entries,known,inflight):
    for e in entries:
        eid=e['endpoint_id'];r=known.get(eid)
        if r and r['status']=='completed':
            if int(r['correlation']['numerator']):return None
            continue
        if eid not in inflight:return eid
    return None


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--results',type=Path,required=True)
    p.add_argument('--workers',type=int,default=24);p.add_argument('--max-bias-nodes',type=int,default=1000000)
    p.add_argument('--dag-samples',type=int,default=4);p.add_argument('--resume',action='store_true')
    args=p.parse_args();root=args.results;folder=root/'reassessment';folder.mkdir(exist_ok=args.resume)
    (folder/'jobs').mkdir(exist_ok=True)
    search=json.loads((root/'search_results.json').read_text());rows=search['rows']
    rankings={r['id']:read_jsonl(root/r['ranked_endpoints_file']) for r in rows}
    jobs={j['id']:j for j in read_jsonl(root/'endpoint_jobs.jsonl.gz')}
    code_paths=[Path(__file__),Path('tools/datf.py'),Path('tools/datf_k06.py'),Path('tools/datf_quadratic.py')]
    config={'workers':args.workers,'initial_max_bias_nodes':args.max_bias_nodes,'dag_samples':args.dag_samples,
            'search_sha256':hashlib.sha256((root/'search_results.json').read_bytes()).hexdigest(),
            'code_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in code_paths},
            'objective':'first nonzero full ordinary-DATF endpoint in each row, or exhaustive model-zero result',
            'unknown_is_not_zero':True,'first_L_omitted':True}
    if args.resume:assert json.loads((folder/'configuration.json').read_text())==config
    else:atomic_json(folder/'configuration.json',config)
    known={}
    for f in (folder/'jobs').glob('*.summary.json'):
        r=json.loads(f.read_text());assert r['id'] in jobs and Path(r['artifact']).exists();known[r['id']]=r
    began=time.time();pending={};inflight=set();worker_pids=set();last_update=0;next_row=0
    def checkpoint(complete=False):
        states=[]
        for row in rows:
            entries=rankings[row['id']];state=row_state(row,entries,known)
            tested=[known[e['endpoint_id']] for e in entries if e['endpoint_id'] in known]
            first=known.get(entries[0]['endpoint_id']) if entries else None
            states.append({**row,**state,'evaluated_unique_endpoints':sum(r['status']=='completed' for r in tested),
                           'first_ranked_endpoint_result':first})
        doc={'workers':args.workers,'worker_pids_observed':sorted(worker_pids),
             'elapsed_seconds_this_run':time.time()-began,'complete':complete,
             'completed_endpoint_jobs':sum(r['status']=='completed' for r in known.values()),
             'resource_limited_jobs':sum(r['status']=='resource_limit' for r in known.values()),
             'inflight':len(inflight),'rows':states}
        atomic_json(folder/'results.json',doc)
        counts={s:sum(r['status']==s for r in states) for s in ('running','found_nonzero','all_endpoints_model_zero','no_nonzero_matches')}
        print('REEVALUATION',json.dumps({**counts,'completed_endpoints':doc['completed_endpoint_jobs'],
              'inflight':len(inflight),'seconds':round(time.time()-began,1)}),flush=True)
        return states
    context=multiprocessing.get_context('spawn')
    with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers,mp_context=context,max_tasks_per_child=32) as ex:
        while True:
            no_candidate_round=0
            while len(pending)<args.workers and no_candidate_round<len(rows):
                row=rows[next_row%len(rows)];next_row+=1
                if row_state(row,rankings[row['id']],known)['status']!='running':
                    no_candidate_round+=1;continue
                eid=next_endpoint(rankings[row['id']],known,inflight)
                if eid is None:no_candidate_round+=1;continue
                no_candidate_round=0
                limit=known[eid]['max_bias_nodes']*4 if eid in known else args.max_bias_nodes
                future=ex.submit(evaluate_endpoint,jobs[eid],str(folder/'jobs'),limit,args.dag_samples)
                pending[future]=eid;inflight.add(eid)
            if not pending:break
            done,_=concurrent.futures.wait(pending,timeout=10,return_when=concurrent.futures.FIRST_COMPLETED)
            for f in done:
                eid=pending.pop(f);inflight.remove(eid);r=f.result();assert r['id']==eid
                known[eid]=r;worker_pids.add(r['worker_pid'])
            if time.time()-last_update>15:checkpoint();last_update=time.time()
    states=checkpoint(complete=True)
    assert all(r['status']!='running' for r in states)
    atomic_json(folder/'completed.json',{'rows':len(states),'workers':args.workers,
                'all_rows_terminal':True,'completed_endpoint_jobs':sum(r['status']=='completed' for r in known.values()),
                'worker_pids_observed':sorted(worker_pids),'elapsed_seconds':time.time()-began})


if __name__=='__main__':main()
