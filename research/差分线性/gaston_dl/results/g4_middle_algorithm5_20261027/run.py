"""Re-evaluate G4a/G4b ORIGINAL middle boundaries using JOC Algorithms 3/5."""
import gzip
import hashlib
import html
import json
import math
import random
import resource
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
G = HERE.parents[1]
ROOT = G.parent
sys.path.insert(0, str(G/'tools'))
sys.path.insert(0, str(ROOT/'DATF算法3复核_20261021'))
sys.path.insert(0, str(G/'results/splice_long_20261027'))
from algorithm3 import Algorithm3, ResourceLimit
from algorithm5 import Algorithm5
from datf import evaluate, support, indices, variable, add, mul, ONE, ZERO, OrdinaryDATF
from datf_k06 import construct, verify_dag, parity_polynomial, value_record
from fixed_tail_search import S, column, embed, linear, transpose, lat, RC
from hnw_middle import one, two

def read(p): return json.loads(p.read_text())
def save(p, obj): p.write_text(json.dumps(obj, ensure_ascii=False, indent=2)+'\n')
def state(a): return tuple(int(w,16) if isinstance(w,str) else w for w in a)
def rec(q): return value_record(q)
def frac(d): return Fraction(int(d['numerator']),int(d['denominator']))
def pool(name, identifier):
    with gzip.open(G/'results/datf_pools/pools'/name, 'rt') as f:
        for line in f:
            row=json.loads(line)
            if row['id']==identifier:return row
    raise KeyError(identifier)

def validate_weighted_backend():
    # Direct nonuniform enumeration independently checks BiasP, including
    # quadratic/nonlinear, constant, deterministic and disconnected cases.
    rng=random.Random(20261027);cases=0
    for n in range(1,8):
        e=OrdinaryDATF()
        for i in range(n):e.fresh()
        for trial in range(20):
            a=Algorithm5(e)
            a.distributions={i:rng.choice([Fraction(0),Fraction(1,2),Fraction(-1,4),Fraction(1),Fraction(-1)]) for i in range(n)}
            f=frozenset(t for t in range(1<<n) if rng.random()<.15)
            expected=Fraction(0)
            for assignment in range(1<<n):
                probability=math.prod((1+(-1 if assignment>>v&1 else 1)*a.distributions[v])/2 for v in range(n))
                expected+=(-1 if evaluate(f,assignment) else 1)*probability
            assert a.bias_p_correlation(f)==expected
            cases+=1
    # The layerwise marginal of beta2 = beta1*x is 3/4, whereas full
    # substitution beta1=x*y would produce 1/2. This catches the old hybrid.
    e=OrdinaryDATF();x=e.fresh();y=e.fresh()
    b=e.fresh(mul(x,y),kind='beta',stage=1)
    c=e.fresh(mul(b,x),kind='beta',stage=2)
    a=e.fresh(mul(x,y),kind='alpha',stage=1)
    alg=Algorithm5(e);alg.initialize_distributions()
    assert alg.distributions[2]==Fraction(1,2)
    assert alg.distributions[3]==Fraction(3,4)
    assert alg.distributions[4]==Fraction(1,2)
    # Once beta2 is substituted, beta1*x is nonisolated and BiasP uses
    # the stored independent marginals (3/4), not full substitution (1/2).
    assert Algorithm5(e).estimate(b)==Fraction(1,2)
    assert Algorithm5(e).estimate(c)==Fraction(3,4)
    return dict(weighted_exhaustive_cases=cases,layerwise_not_recursive_marginal='PASS',
                nonuniform_alpha='PASS',isolated_substitution='PASS')

def main():
    resource.setrlimit(resource.RLIMIT_AS,(4<<30,4<<30))
    (HERE/'raw').mkdir(exist_ok=True)
    sources=[G/'results/hnw26_comparison/configuration.json',
             G/'results/datf_pools/linear_two_match/hnw26/results.json',
             ROOT/'路径分段汇总_20261027/results.json',
             ROOT/'文献/Differential-Linear Cryptanalysis from an Algebraic Perspective.pdf']
    c_a=next(x['candidate'] for x in read(sources[0])['jobs'] if x['id']=='D1_M2_L1_pool')
    c_b=read(sources[1])['candidate']
    prior={r['code']:r for r in read(sources[2])['rows']}
    validation=validate_weighted_backend()
    rows=[];jobs=[]
    for code,c in [('G4a',c_a),('G4b',c_b)]:
        delta=state(c['middle_input_difference']);mask=state(c['middle_output_mask'])
        assert linear(state(c['differential_output']))==delta
        probability=math.prod(Fraction(sum(S[x]^S[x^column(state(c['input_difference']),j)]==column(state(c['differential_output']),j) for x in range(32)),32) for j in range(64))
        assert probability==frac(c['differential_probability'])
        tail=pool('linear_one.jsonl.gz' if code=='G4a' else 'linear_two_2to2.jsonl.gz',c['linear_id'])
        assert transpose(state(tail['input']))==mask
        trail=tail.get('trail',[dict(input=tail['input'],output=tail['output'])])
        q=Fraction(1)
        for i,t in enumerate(trail):
            q*=math.prod(Fraction(lat(column(state(t['input']),j),column(state(t['output']),j)),32) for j in range(64))
            if i:assert state(trail[i-1]['output'])==transpose(state(t['input']))
        assert q*q==frac(c['linear_squared_correlation'])
        assert state(tail['output'])==state(c['output_mask'])
        assert c['input_difference']==prior[code]['delta_core']
        assert c['output_mask']==prior[code]['lambda_core']
        start=c['start_round_index']+c['r_d'];rounds=c['r_m']
        job=dict(id=code+'_middle',algorithm='Gaston',rounds=rounds,start_round_index=start,
                 delta_core=c['middle_input_difference'],lambda_core=c['middle_output_mask'])
        jobs.append(job)
        started=time.perf_counter();e,base,diff,stages=construct(delta,rounds,start)
        verify_dag(e,base,diff,delta,rounds,start,samples=16)
        target=parity_polynomial(diff,mask)
        construction_seconds=time.perf_counter()-started
        dagpath=HERE/'raw'/f'{code}_dag.json'
        save(dagpath,dict(job=job,target_terms=[hex(t) for t in sorted(target)],dag=e.serialized_dag()))
        row=dict(code=code,split=[c['r_d'],rounds,c['r_l']],candidate=c,middle_job=job,
                 construction_and_16_dag_checks_seconds=construction_seconds,
                 middle_bit_weights=dict(input=sum(w.bit_count() for w in delta),output=sum(w.bit_count() for w in mask)),
                 p=rec(probability),q_squared=rec(q*q),sampling=prior[code]['sampling'],
                 previous_full_algorithm3=prior[code]['datf_algorithm3'],dag_file=str(dagpath))
        for name,cls in [('algorithm3',Algorithm3),('algorithm5',Algorithm5)]:
            started=time.perf_counter();alg=cls(e,seconds=120,max_nodes=500000,max_terms=1000000)
            result=dict(job=job,method='JOC DATF '+name,status='running',correlation=None)
            try:
                value=alg.estimate(target)
                result.update(status='completed',correlation=rec(value),bias=rec(value/2),
                              splice_single_trail=rec(probability*value*q*q))
                if name=='algorithm3':
                    result['expanded_polynomial']=[hex(t) for t in sorted(alg.final_polynomial)]
                    result['zero_certificate']=alg.zero_certificate
                    rng=random.Random(20261027)
                    for _ in range(32):
                        values=0
                        for v,m in enumerate(e.metadata):
                            bit=evaluate(e.definitions[v],values) if m['kind']=='beta' else rng.getrandbits(1)
                            values|=bit<<v
                        assert evaluate(target,values)==evaluate(alg.final_polynomial,values)
                    result['expansion_checks']=32
                else:
                    result['trace']=alg.trace
                    result['marginal_correlations']={str(v):rec(x) for v,x in alg.distributions.items()}
            except (ResourceLimit,MemoryError,RecursionError) as exc:
                result.update(status='resource_limit',reason=str(exc),correlation=None)
            result['elapsed_seconds']=time.perf_counter()-started
            result['statistics']=dict(alg.stats) if name=='algorithm3' else alg.statistics()
            result['phase']=alg.phase
            save(HERE/'raw'/f'{code}_{name}.json',result)
            row[name]={k:v for k,v in result.items() if k not in ['expanded_polynomial','marginal_correlations','trace']}
            print(code,name,result['status'],result['correlation'],result['elapsed_seconds'],flush=True)
        # Independent local operator: one chi exact DLCT; two chi HNW zero model.
        value=one(delta,mask) if rounds==1 else two(delta,mask)[0]
        row['independent_middle_crosscheck']=dict(method='exact one-chi DLCT' if rounds==1 else 'HNW zero model, affine/quadratic evaluation',correlation=rec(value))
        for name in ['algorithm3','algorithm5']:
            assert row[name]['status']=='completed'
            assert frac(row[name]['correlation'])==value
        # All four first-DDT branches, for THIS one fixed linear trail only.
        j=next(j for j in range(64) if column(state(c['input_difference']),j))
        a=column(state(c['input_difference']),j);branches=[];summed=Fraction(0)
        for b in range(32):
            count=sum(S[x]^S[x^a]==b for x in range(32))
            if not count:continue
            d=linear(embed(j,b));cm=one(d,mask) if rounds==1 else two(d,mask)[0]
            p=Fraction(count,32);summed+=p*cm*q*q
            branches.append(dict(first_output_pattern=b,probability=rec(p),middle_crosscheck=rec(cm)))
        row['fixed_linear_trail_all_first_DDT_branches']=dict(branches=branches,splice_sum=rec(summed),
            note='Independent local/HNW middle evaluation; this is not a sum over other linear trails.')
        rows.append(row)
    save(HERE/'manifest.json',dict(jobs=jobs))
    result=dict(title='G4a / G4b 原始中间段：DATF 算法 3 与折中算法 5',
                date='2026-10-27',scope='Only original 2/1-chi middle segments; no new full-four-round evaluation',
                first_L_omitted=True,first_constant_absorbed_into_uniform_middle_input=True,
                transitional_cut='before each chi after the first, after L and iota; reference-code DAG simplifications',
                algorithm5='Paper page 19, Algorithm 5; all immediate-predecessor alpha/beta marginals use BiasP; no enhanced DATF',
                quantity='signed correlation = twice paper bias',rows=rows)
    save(HERE/'results.json',result)
    validation.update(status='PASS',dag_checks=32,algorithm3_expansion_checks=64,
                      outer_DDT_and_tail_LAT_and_boundary_checks='PASS',
                      algorithms_3_5_and_independent_middle_crosscheck='PASS',
                      new_full_endpoint_jobs=0,whole_endpoint_resource_limits_preserved=True)
    save(HERE/'validation.json',validation)
    sources += [HERE/'algorithm5.py',HERE/'run.py',ROOT/'DATF算法3复核_20261021/algorithm3.py',G/'tools/datf.py',G/'tools/datf_k06.py',G/'tools/datf_quadratic.py',G/'tools/fixed_tail_search.py',G/'tools/linear_seed_scan.py',G/'results/splice_long_20261027/hnw_middle.py']
    sources += [G/'results/datf_pools/pools'/name for name in ['linear_one.jsonl.gz','linear_two_2to2.jsonl.gz']]
    save(HERE/'source_manifest.json',dict(files=[dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in sources]))
    print('PASS',validation,flush=True)

if __name__=='__main__':main()
