#!/usr/bin/env python3
import json
import random
import sys
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'tools'))
from datf import OrdinaryDATF,evaluate
from datf_audit import quadratic_correlation
from datf_quadratic import QuadraticOrdinaryDATF,uniform_quadratic
from datf_k06 import construct,parity_polynomial

rng=random.Random(20261016099)
for trial in range(200):
    n=1+trial%10
    terms=[0]+[1<<i for i in range(n)]+[(1<<i)|(1<<j) for i in range(n) for j in range(i+1,n)]
    p=frozenset(t for t in terms if rng.randrange(2))
    actual=Fraction(sum((-1)**evaluate(p,a) for a in range(1<<n)),1<<n)
    assert uniform_quadratic(p)==actual==quadratic_correlation(p)
    original=OrdinaryDATF();fast=QuadraticOrdinaryDATF()
    probs={i:rng.choice([Fraction(0),Fraction(0),Fraction(1,2),Fraction(-1,4),Fraction(1)]) for i in range(n)}
    original.marginal=probs.__getitem__;fast.marginal=probs.__getitem__
    assert original.independent_correlation(p)==fast.independent_correlation(p)

sources=[ROOT/'results/datf_pools/search_final/endpoints/results.json',ROOT/'results/datf_k06/results.json',ROOT/'results/datf_one_one_tail/results.json']
checked=0;nonzero=0
for source in sources:
    doc=json.loads(source.read_text())
    for job in doc['results']:
        if job.get('status')!='completed':continue
        delta=tuple(int(w,16) for w in job['input_difference']['words'])
        mask=tuple(int(w,16) for w in job['output_mask']['words'])
        engine,_,diff,_=construct(delta,job['rounds'],job['start_round_index'],engine_class=QuadraticOrdinaryDATF)
        rho=engine.estimate(parity_polynomial(diff,mask))
        assert rho==Fraction(int(job['correlation']['numerator']),int(job['correlation']['denominator'])),job['id']
        checked+=1;nonzero+=bool(rho)
print(json.dumps({'status':'PASS','exhaustive_small_quadratics':200,'weighted_shannon_comparisons':200,'saved_ordinary_DATF_endpoints':checked,'nonzero_saved_endpoints':nonzero}))
