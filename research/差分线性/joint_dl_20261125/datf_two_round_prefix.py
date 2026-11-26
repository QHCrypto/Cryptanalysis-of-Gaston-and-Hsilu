"""Two-round prefix weighting with original DATF Algorithm 3 on the tail.

For each fixed tail input difference d, Algorithm 3 fully expands the beta
variables of chi-L-chi. The last-layer alpha variables are independent uniform
formal variables in that model. Averaging them first gives the squared-LAT
expansion, with full masks pulled through L^T. This identity holds for every d.
Affine elimination then computes the COMPLETE weighted prefix sum; it does not
sample or retain only selected prefix differences. Prefix weights use the same
differential Markov convention as the manuscript's HNW calculation.
"""
from pathlib import Path
import sys
import json
import random
import hashlib
import time
from collections import Counter
from fractions import Fraction
from itertools import product

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT/'DATF算法3复核_20261021'))
sys.path.insert(0, str(ROOT/'gaston_dl/tools'))
sys.path.insert(0, str(ROOT/'hsilu_dl/tools'))
from algorithm3 import Algorithm3
import datf as D
import datf_k06 as G
import hsilu as H
import core as C
from joint_estimator import affine_character, span_basis, fraction_record


def local_tables():
    """Compute both local tables through the original DATF Algorithm-3 backend."""
    engine = D.OrdinaryDATF()
    x = [engine.fresh(row=i,column=0) for i in range(5)]
    sx = [D.add(x[i],x[(i+2)%5],D.mul(x[(i+1)%5],x[(i+2)%5])) for i in range(5)]
    walsh = [[Fraction(0) for _ in range(32)] for _ in range(32)]
    act = [[Fraction(0) for _ in range(32)] for _ in range(32)]
    for b in range(32):
        component = D.add(*(sx[i] for i in range(5) if (b>>i)&1))
        for a in range(32):
            p = D.add(component,*(x[i] for i in range(5) if (a>>i)&1))
            walsh[b][a] = Algorithm3(engine).estimate(p)
            direct = Fraction(sum(1-2*(((a&u).bit_count()^(b&C.sbox(u)).bit_count())&1)
                                  for u in range(32)),32)
            assert walsh[b][a] == direct
    for d in range(32):
        bits = [(d>>i)&1 for i in range(5)]
        diff = []
        for i in range(5):
            k,l = (i+1)%5,(i+2)%5
            diff.append(D.add(D.ONE if bits[i]^bits[l]^(bits[k]&bits[l]) else D.ZERO,
                              x[k] if bits[l] else D.ZERO,x[l] if bits[k] else D.ZERO))
        for b in range(32):
            p = D.add(*(diff[i] for i in range(5) if (b>>i)&1))
            act[d][b] = Algorithm3(engine).estimate(p)
            direct = Fraction(sum(1-2*((b&(C.sbox(u)^C.sbox(u^d))).bit_count()&1)
                                  for u in range(32)),32)
            assert act[d][b] == direct
    supports = [[a for a in range(32) if walsh[b][a]] for b in range(32)]
    bases = [span_basis(a^aa[0] for a in aa) for aa in supports]
    for b in range(32):
        assert sum(c*c for c in walsh[b]) == 1
        for d in range(32):
            value = 0 if any((v&d).bit_count()&1 for v in bases[b]) else 1-2*((supports[b][0]&d).bit_count()&1)
            assert act[d][b] == value
    return walsh,act,supports,bases


def last_layer_terms(name,walsh,supports):
    LT = C.lgt if name=='Gaston' else C.lht
    _,mask = C.ENDPOINTS[name]
    cols = [j for j in range(64) if C.column(mask,j)]
    terms = []
    for choices in product(*(supports[C.column(mask,j)] for j in cols)):
        alpha = C.ZERO
        weight = Fraction(1)
        for j,a in zip(cols,choices):
            alpha = C.add(alpha,C.embed(j,a))
            weight *= walsh[C.column(mask,j)][a]**2
        terms.append((alpha,LT(alpha),weight))
    assert sum(w for _,_,w in terms) == 1
    return terms


def tail_formula(d,terms,act):
    result = Fraction(0)
    for _,b,w in terms:
        value = w
        for j in range(64):
            value *= act[C.column(d,j)][C.column(b,j)]
            if not value:break
        result += value
    return result


def prefix_sum(name,walsh,supports,bases):
    start = time.perf_counter()
    L,LT = (C.lg,C.lgt) if name=='Gaston' else (C.lh,C.lht)
    delta,mask = C.ENDPOINTS[name]
    j0 = next(j for j in range(64) if C.column(delta,j))
    d0 = C.column(delta,j0)
    first = Counter(C.sbox(x)^C.sbox(x^d0) for x in range(32))
    terms = last_layer_terms(name,walsh,supports)
    total = Fraction(0)
    records = []
    prefixes = []
    for d1,count in sorted(first.items()):
        d2 = L(C.embed(j0,d1))
        offset = L(C.chi(d2))
        probability = Fraction(count,32)
        # Every local quadratic derivative has an affine image. The dimension
        # below is the log2 of this branch's full second-round difference set.
        rank = sum(len(span_basis(C.sbox(x)^C.sbox(x^C.column(d2,j))^C.sbox(C.column(d2,j))
                                  for x in range(32))) for j in range(64))
        prefixes.append({'first_column_difference':d1,'probability':str(probability),
                         'second_layer_difference_image_dimension':rank})
        def adjoint(mask):
            return C.pack(C.polar_transpose(LT(mask),d2))
        for alpha,b3,weight in terms:
            ks = [C.embed(j,v) for j in range(64) for v in bases[C.column(b3,j)]]
            origin = C.add(*(C.embed(j,supports[C.column(b3,j)][0]) for j in range(64)))
            value,rk,reason = affine_character([adjoint(k) for k in ks],
                                               [C.dot(k,offset) for k in ks],
                                               adjoint(origin),C.dot(origin,offset))
            contribution = probability*weight*value
            total += contribution
            records.append({'first_column_difference':d1,'first_probability':str(probability),
                            'last_chi_input_mask':[f'{v:016x}' for v in alpha],
                            'last_chi_squared_weight':str(weight),'prefix_average':str(value),
                            'rank':rk,'reason':reason,'contribution':str(contribution)})
    return {'target':name,'delta':[f'{v:016x}' for v in delta],
            'lambda':[f'{v:016x}' for v in mask],'correlation':fraction_record(total),
            'bias':fraction_record(total/2),'first_prefix_branches':prefixes,
            'last_mask_branches':len(terms),'all_summed_cases':len(records),
            'nonzero_cases':sum(bool(Fraction(r['contribution'])) for r in records),
            'branches':records,'seconds':time.perf_counter()-start}


def validate_tail(name,walsh,act,supports):
    """Run full original two-round DATF Algorithm 3, not a renamed HNW call."""
    rng = random.Random(2026112507+(name=='Hsilu'))
    L = C.lg if name=='Gaston' else C.lh
    delta,mask = C.ENDPOINTS[name]
    j0 = next(j for j in range(64) if C.column(delta,j))
    d0 = C.column(delta,j0)
    first = sorted({C.sbox(x)^C.sbox(x^d0) for x in range(32)})
    terms = last_layer_terms(name,walsh,supports)
    cases = [('zero',C.ZERO)]
    for d1 in first:
        d2 = L(C.embed(j0,d1))
        for _ in range(24):
            y = tuple(rng.getrandbits(64) for _ in range(5))
            d = L(C.add(C.chi(y),C.chi(C.add(y,d2))))
            cases.append((f'prefix_{d1}',d))
    for _ in range(16):
        cases.append(('unrestricted',tuple(rng.getrandbits(64) for _ in range(5))))
    for _ in range(8):
        cases.append(('single_column',C.embed(rng.randrange(64),rng.randrange(1,32))))
    records=[]
    dag_checks=0
    for index,(kind,d) in enumerate(cases):
        if name=='Gaston':engine,base,diff,_=G.construct(d,2,10)
        else:engine,base,diff=H.construct(d,2,start=10)
        # First chi input -> chi3 -> L -> C_0x4b -> chi4.
        assert sum(m['kind']=='alpha' for m in engine.metadata)==320
        p=G.parity_polynomial(diff,mask)
        alg=Algorithm3(engine,seconds=30)
        value=alg.estimate(p)
        expected=tail_formula(d,terms,act)
        assert value==expected,(name,d,value,expected)
        assert alg.stats['expanded_degree']<=2
        if index%24==0:
            x=tuple(rng.getrandbits(64) for _ in range(5))
            assignment=engine.evaluate_dag([(v>>j)&1 for v in x for j in range(64)])
            direct=C.dot(mask,C.add(C.permutation(x,name,2),C.permutation(C.add(x,d),name,2)))
            assert D.evaluate(p,assignment)==direct
            dag_checks+=1
        records.append({'kind':kind,'difference':[f'{v:016x}' for v in d],
                        'algorithm3_correlation':str(value),'fullmask_formula_correlation':str(expected),
                        'expanded_degree':alg.stats['expanded_degree'],
                        'expanded_terms':alg.stats['expanded_terms']})
    return {'target':name,'cases':len(cases),'nonzero_cases':sum(bool(Fraction(r['algorithm3_correlation'])) for r in records),
            'full_state_dag_checks':dag_checks,'status':'PASS','records':records}


if __name__=='__main__':
    started=time.perf_counter()
    walsh,act,supports,bases=local_tables()
    print('DATF Algorithm 3 local tables: 2048 exact truth-table checks PASS',flush=True)
    results=[];checks=[]
    for name in ['Gaston','Hsilu']:
        r=prefix_sum(name,walsh,supports,bases)
        results.append(r)
        print(name,'weighted two-round-prefix DATF',r['correlation'],flush=True)
        check=validate_tail(name,walsh,act,supports)
        checks.append(check)
        print(name,'full two-round DATF checks',check['cases'],'nonzero',check['nonzero_cases'],'PASS',flush=True)
    doc={'method':'Two-round differential-prefix weighting plus original DATF Algorithm 3 on chi-L-chi tail',
         'quantity':'signed correlation; bias is half of correlation',
         'connection_model':'independent uniform tail base state, as in the manuscript HNW comparison',
         'prefix_model':'full two-round differential Markov distribution, integrated exactly as a mixture of affine images',
         'tail_identity':'Fully expand DATF beta variables; average last-layer independent uniform alpha variables; squared-LAT expansion and complete L^T masks are an equivalent evaluation order for this two-round DATF model.',
         'scope':'Exact rational model evaluation, not a proof of the actual fixed-permutation correlation',
         'local_algorithm3_truth_table_checks':2048,'results':results,'tail_validation':checks,
         'seconds':time.perf_counter()-started}
    sources=[Path(__file__),HERE/'core.py',HERE/'joint_estimator.py',
             ROOT/'DATF算法3复核_20261021/algorithm3.py',ROOT/'gaston_dl/tools/datf.py',
             ROOT/'gaston_dl/tools/datf_k06.py',ROOT/'gaston_dl/tools/datf_quadratic.py',
             ROOT/'hsilu_dl/tools/hsilu.py']
    doc['sources']=[{'path':str(p.relative_to(ROOT.parent)),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]
    (HERE/'datf_prefix_fullmask_results.json').write_text(json.dumps(doc,indent=2)+'\n')
