import argparse
import hashlib
import json
import os
import random
import resource
import signal
import sys
import time
from pathlib import Path
from algorithm3 import Algorithm3, ResourceLimit, ROOT
from datf import evaluate, indices
from datf_k06 import construct, parity_polynomial, verify_dag, value_record

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--manifest',type=Path,required=True)
    ap.add_argument('--index',type=int,required=True)
    ap.add_argument('--output',type=Path,required=True)
    ap.add_argument('--seconds',type=int,default=120)
    ap.add_argument('--max-terms',type=int,default=1000000)
    ap.add_argument('--max-products',type=int,default=100000000)
    ap.add_argument('--max-nodes',type=int,default=200000)
    ap.add_argument('--input-shannon',action='store_true')
    args=ap.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,(4<<30,4<<30))
    def alarm(signum,frame):raise ResourceLimit('hard wall time limit')
    signal.signal(signal.SIGALRM,alarm)
    signal.alarm(args.seconds+15)
    job=json.loads(args.manifest.read_text())['jobs'][args.index]
    started=time.perf_counter()
    result={**job,'method':'original DATF Algorithm 3','status':'running','correlation':None,
            'bias':None,'limits':{'seconds':args.seconds,'max_terms':args.max_terms,
            'max_products':args.max_products,'max_nodes':args.max_nodes,'address_space_bytes':4<<30},
            'first_L_omitted':True,'remaining_alpha_assumed_jointly_uniform':True}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    def save(a=None):
        if a: result.update(phase=a.phase,statistics=dict(a.stats))
        result['elapsed_seconds']=time.perf_counter()-started
        result['peak_rss_kib']=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        tmp=args.output.with_suffix('.tmp');tmp.write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n');tmp.replace(args.output)
    alg=None
    try:
        delta=tuple(int(w,16) for w in job['delta_core']);mask=tuple(int(w,16) for w in job['lambda_core'])
        rounds,start=job['rounds'],job['start_round_index']
        if job['algorithm']=='Gaston':
            e,base,diff,stages=construct(delta,rounds,start)
            verify_dag(e,base,diff,delta,rounds,start,samples=4)
        else:
            sys.path.insert(0,str(ROOT/'hsilu_dl/tools'))
            import hsilu
            e,base,diff=hsilu.construct(delta,rounds,start=start)
            rng=random.Random(1011)
            for _ in range(4):
                state=tuple(rng.getrandbits(64) for _ in range(5))
                values=e.evaluate_dag([(w>>j)&1 for w in state for j in range(64)])
                base0=hsilu.core(state,rounds,start)
                base1=hsilu.core(tuple(x^y for x,y in zip(state,delta)),rounds,start)
                assert tuple(sum(evaluate(p,values)<<j for j,p in enumerate(row)) for row in base)==base0
                assert tuple(sum(evaluate(p,values)<<j for j,p in enumerate(row)) for row in diff)==tuple(x^y for x,y in zip(base0,base1))
        result['dag_verified_full_states']=4
        p=parity_polynomial(diff,mask)
        alg=Algorithm3(e,args.max_terms,args.max_products,args.max_nodes,args.seconds,save)
        ans=alg.estimate_input_shannon(p) if args.input_shannon else alg.estimate(p)
        # Expansion is checked using independently drawn formal free variables;
        # beta values are computed from the original, unexpanded DAG.
        if not args.input_shannon:
            rng=random.Random(20261021)
            for _ in range(16):
                values=0
                for v,m in enumerate(e.metadata):
                    bit=evaluate(e.definitions[v],values) if m['kind']=='beta' else rng.getrandbits(1)
                    values|=bit<<v
                assert evaluate(p,values)==evaluate(alg.final_polynomial,values)
        result.update(status='completed',correlation=value_record(ans),bias=value_record(ans/2),
                      expansion_verified_formal_assignments=0 if args.input_shannon else 16,zero_certificate=alg.zero_certificate)
        # Retain the complete expanded ANF so completed results are reviewable.
        expanded=args.output.with_suffix('.anf.json')
        expanded.write_text(json.dumps({'terms':None if args.input_shannon else [hex(t) for t in sorted(alg.final_polynomial)],
                                       'exact_Shannon_branches':alg.branches,
                                       'variable_metadata':e.metadata},ensure_ascii=False)+'\n')
        result['expanded_anf_file']=str(expanded)
        result['expanded_anf_sha256']=hashlib.sha256(expanded.read_bytes()).hexdigest()
    except (ResourceLimit,MemoryError,RecursionError) as exc:
        result.update(status='resource_limit',reason=type(exc).__name__+': '+str(exc),correlation=None,bias=None)
    except Exception as exc:
        result.update(status='error',reason=repr(exc))
        raise
    finally:
        signal.alarm(0)
        result['evaluation_order']='exact input Shannon then full beta expansion per branch' if args.input_shannon else 'full beta expansion then Bias'
        if alg and args.input_shannon:result['branches']=alg.branches
        save(alg)
    print(job['id'],result['status'],result['correlation'],round(result['elapsed_seconds'],3),flush=True)

if __name__=='__main__':main()
