#!/usr/bin/env python3
"""Independent arithmetic, selection, and local-model equivalence checks."""
import argparse
import collections
import gzip
import itertools
import json
from fractions import Fraction
from functools import lru_cache
from pathlib import Path

from fixed_tail_search import S, column, embed, linear, transpose
from hnw26_compare import ROOT, fraction, save


def reference_zero(candidate):
    """Independent demand-driven Fraction oracle, not the C++ vector updater."""
    delta=tuple(int(x,16) for x in candidate['input_difference'])
    target=tuple(int(x,16) for x in candidate['output_mask'])
    R=candidate['R']
    ddt=[[sum(S[x]^S[x^a]==b for x in range(32)) for b in range(32)] for a in range(32)]
    dlct=[[Fraction(sum((-1)**((v&(S[x]^S[x^a])).bit_count()%2) for x in range(32)),32) for v in range(32)] for a in range(32)]
    kernel=[]
    for v in range(32):
        terms=[]
        for u in range(32):
            n=sum((-1)**(((u&x).bit_count()+(v&S[x]).bit_count())%2) for x in range(32))
            if n:terms.append((u,Fraction(n*n,1024)))
        kernel.append(terms)
    active=[j for j in range(64) if column(delta,j)]
    choices=[[(b,Fraction(n,32)) for b,n in enumerate(ddt[column(delta,j)]) if n] for j in active]
    total=Fraction(0)
    @lru_cache(None)
    def linear_terms(j,u):
        mask=transpose(embed(j,u))
        return [(k,column(mask,k)) for k in range(64) if column(mask,k)]
    for branch in itertools.product(*choices):
        difference=[0]*5;prob=Fraction(1)
        for j,(b,p) in zip(active,branch):
            prob*=p
            for i in range(5):difference[i]^=((b>>i)&1)<<j
        after=linear(difference)
        @lru_cache(None)
        def q(r,j,v):
            if not v:return Fraction(1)
            if r==1:return Fraction((-1)**((v&column(difference,j)).bit_count()%2))
            if r==2:return dlct[column(after,j)][v]
            value=Fraction(0)
            for u,weight in kernel[v]:
                product=weight
                for k,mask in linear_terms(j,u):
                    product*=q(r-1,k,mask)
                    if not product:break
                value+=product
            return value
        value=Fraction(1)
        for j in range(64):
            value*=q(R,j,column(target,j))
            if not value:break
        total+=prob*value
    return total


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--results',type=Path,default=ROOT/'results/hnw26_comparison')
    args=parser.parse_args();root=args.results
    data=json.loads((root/'results.json').read_text());config=json.loads((root/'configuration.json').read_text())
    assert data['status']=='completed'
    source=json.loads((ROOT/'results/datf_pools/merged/merged_results.json').read_text())
    original={r['id']:r for r in source['rows']}
    ddt=[[sum(S[x]^S[x^a]==b for x in range(32)) for b in range(32)] for a in range(32)]
    checks=0;correlations=[];probability_branches=0
    for row in data['results']:
        r=original[row['id']];assert row['candidate']==(r['selected'] or r['best'])
        assert row['datf_correlation']==(r['selected_correlation'] or r['first_correlation'])
        b=row['candidate'];delta=[int(w,16) for w in b['input_difference']]
        expected=1
        for j in range(64):expected*=sum(n>0 for n in ddt[column(delta,j)])
        for mode,saved in row['hnw'].items():
            job=json.loads(Path(saved['artifact']).read_text());seen=set()
            assert job['completed_branches']==expected
            for branch in job['branches']:
                words=tuple(int(w,16) for w in branch['first_output_difference'].split(';'))
                assert words not in seen;seen.add(words);p=Fraction(1)
                for j in range(64):p*=Fraction(ddt[column(delta,j)][column(words,j)],32)
                assert p>0 and p==fraction(branch['probability']);probability_branches+=1
            assert fraction(job['correlation'])==sum((fraction(x['contribution']) for x in job['branches']),Fraction(0))
            assert fraction(job['bias'])*2==fraction(job['correlation']);checks+=1
        if all(m in row['hnw'] for m in ('zero','full')):
            assert fraction(row['hnw']['zero']['correlation'])==fraction(row['hnw']['full']['correlation'])
        if row['category']=='DATF_nonzero':
            oracle=reference_zero(b);assert oracle==fraction(row['hnw'].get('zero',row['hnw'].get('full'))['correlation'])
            correlations.append({'id':row['id'],'independent_fraction_oracle':str(oracle)})
    # Exact structural reason that full and zero-value-mask variants coincide
    # for these endpoints. This is a statement about the factorized HNW model.
    masks=[[transpose(embed(j,a)) for a in range(32)] for j in range(64)]
    witnesses=[]
    for j in range(64):
        for a in range(1,32):
            for b in range(32):
                k=next((k for k in range(64) if column(masks[j][a],k) and column(masks[j][b],k) in (0,column(masks[j][a],k))),None)
                assert k is not None
                witnesses.append([j,a,b,k])
    inputs={tuple(r['candidate']['input_difference']) for r in data['results']};first_checks=[]
    for words in sorted(inputs):
        delta=[int(w,16) for w in words];inactive={j for j in range(64) if column(delta,j)==0}
        for j in range(64):
            for a in range(1,32):assert any(column(masks[j][a],k) for k in inactive)
        first_checks.append({'input_words':list(words),'nonzero_value_masks_checked':64*31})
    with gzip.open(root/'linear_mask_witnesses.json.gz','wt') as f:json.dump(witnesses,f,separators=(',',':'))
    save(root/'model_equivalence.json',{'status':'PASS','scope':'These endpoints under the HNW local-column factorization; not true Gaston correlation.',
          'first_conditional_linear_layer':first_checks,'uniform_marginal_annihilation_checks':len(witnesses),
          'witness_format':['output_column','base_value_mask_a','difference_mask_b','input_column_k'],
          'witness_file':'linear_mask_witnesses.json.gz',
          'argument':[
              'At the first L, every nonzero local base-value mask touches an inactive first-round input column, whose conditional base value remains uniform and whose difference is zero. The corresponding product factor is zero.',
              'With only a=0 input coordinates, both unprimed and primed state marginals are uniform. Applying a bijective local chi preserves both uniform marginals.',
              'For every a!=0,b, L-transposed masks have a column with u0!=0 and u1=0 or u1=u0. Its factor is a nontrivial mask of the unprimed or primed uniform marginal and vanishes.',
              'Thus all a!=0 coordinates vanish before every later chi. Round constants multiply only these coordinates by signs and have no effect in this model. The a=0 update reduces to LAT squared, proving the two computed variants coincide for these endpoints.']})
    save(root/'validation.json',{'status':'PASS','endpoints':len(data['results']),'native_evaluations':checks,
          'first_DDT_probabilities_independently_checked':probability_branches,
          'independent_exact_Fraction_oracles':correlations,
          'full_and_zero_agree_on_all_selected_endpoints':True,
          'structural_equivalence_certificate':'model_equivalence.json',
          'backend_validation':json.loads((root/'backend_validation.json').read_text())})
    print('PASS',len(data['results']),'endpoints',checks,'evaluations',probability_branches,'DDT branch probabilities',flush=True)


if __name__=='__main__':main()
