#!/usr/bin/env python3
"""Revalidate existing Hsilu SAT trails and join two-chi prefixes to L2 tails.

When the pulled second-DDT constraint matrix has the same rank after removing
the first active columns, the two-round characteristic probability is EXACT.
The middle derivative may then be identically constant on that characteristic.
Neither statement makes a claim about cancellation with other characteristics.
"""
import json
from fractions import Fraction
from pathlib import Path
from hsilu import *
from pool_search import chi5

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT.parents[1]/'Hsilu_trail/pysat_kissat'
OUT=ROOT/'results/sparse_diff'
DDT=[[sum((chi5(x)^chi5(x^a))==b for x in range(32)) for b in range(32)] for a in range(32)]
def column(a,j):return sum(((a[i]>>j)&1)<<i for i in range(5))
def rank(rows):
    basis={}
    for value in rows:
        while value:
            pivot=value.bit_length()-1
            if pivot in basis:value^=basis[pivot]
            else:basis[pivot]=value;break
    return len(basis)
def rank_certificate(x0,x1):
    equations=[]
    for j in range(64):
        d=column(x1,j)
        for i in range(5):
            a=(((d>>((i+2)%5))&1)<<((i+1)%5))^(((d>>((i+1)%5))&1)<<((i+2)%5))
            if a:equations.append(pack(transpose(embed(j,a))))
    fixed=pack(tuple(x0[0]|x0[1]|x0[2]|x0[3]|x0[4] for _ in range(5)))
    full=rank(equations);outside=rank([r&~fixed for r in equations])
    return {'full_constraint_rank':full,'outside_first_active_columns_rank':outside,
            'exact_product_probability_certified':full==outside,
            'argument':'Outside first active columns, first-chi outputs are independent uniform. Equal ranks make every reachable RHS equiprobable regardless of the conditioned active columns.'}
def middle(delta,gamma):
    coefficients=tuple((gamma[(i+4)%5]&delta[(i+1)%5])^(gamma[(i+3)%5]&delta[(i+4)%5]) for i in range(5))
    if any(coefficients):return 0
    constant=tuple(delta[i]^delta[(i+2)%5]^(delta[(i+1)%5]&delta[(i+2)%5]) for i in range(5))
    return -1 if sum((a&b).bit_count() for a,b in zip(gamma,constant))%2 else 1
def main():
    OUT.mkdir(exist_ok=False);seen={};valid=0;rejected=[]
    for path in sorted(SOURCE.rglob('*.trail.txt')):
        lines=path.read_text().splitlines();states={}
        try:
            for j,line in enumerate(lines):
                if line in ('x0','y0','x1','y1','x2'):
                    states[line]=tuple(int(w.replace('.','0'),2) for w in lines[j+1:j+6])
            if not all(k in states for k in ('x0','y0','x1','y1')):continue
            assert linear(states['y0'])==states['x1']
            if 'x2' in states:assert linear(states['y1'])==states['x2']
            weights=[]
            for r in (0,1):
                counts=[DDT[column(states['x'+str(r)],j)][column(states['y'+str(r)],j)] for j in range(64)]
                assert all(counts)
                weights.append(sum(5-(n.bit_length()-1) for n in counts))
            key=min(tuple(tuple(rol(x,s) for x in states[k]) for k in ('x0','y0','x1','y1')) for s in range(64))
            valid+=1
            if key not in seen:seen[key]={'sources':[],'weights':weights}
            seen[key]['sources'].append(str(path))
        except (AssertionError,ValueError,IndexError):rejected.append(str(path))
    prefixes=[]
    for i,(states,record) in enumerate(sorted(seen.items())):
        x0,y0,x1,y1=states;cert=rank_certificate(x0,x1)
        assert cert['full_constraint_rank']==record['weights'][1]
        prefixes.append({'id':f'D2_{i:02d}',**record,'states':{k:encode(v) for k,v in zip(('x0','y0','x1','y1'),states)},'certificate':cert})
    pool=json.loads((ROOT/'results/baseline/linear_pool.json').read_text())
    heads=[(h,tuple(int(w,16) for w in h['join_mask']['words'])) for h in pool['heads']]
    best={};matches=0
    for prefix in prefixes:
        states={k:tuple(int(w,16) for w in v['words']) for k,v in prefix['states'].items()}
        if not prefix['certificate']['exact_product_probability_certified']:continue
        weight=sum(prefix['weights']);before=linear(states['y1'])
        for rotation in range(64):
            delta=tuple(rol(x,rotation) for x in before)
            input_delta=tuple(rol(x,rotation) for x in states['x0'])
            for head,gamma in heads:
                sign=middle(delta,gamma)
                if not sign:continue
                matches+=1
                for last in pool['bridges'][head['bridge_id']]['last_chi_choices']:
                    w=weight+2*(head['first_chi']['weight']+last['weight'])
                    score=Fraction(sign,1<<w);key=(input_delta,tuple(last['output']['words']))
                    if key not in best or abs(score)>abs(Fraction(best[key]['splice_score'])):
                        best[key]={'rounds':5,'input_difference':encode(input_delta),'core_output_mask':last['output'],
                                   'prefix_id':prefix['id'],'rotation':rotation,'differential_weight':weight,
                                   'restricted_prefix_contribution':str(Fraction(sign,1<<weight)),
                                   'middle_constant_correlation':sign,'splice_score':str(score),
                                   'head':head,'last_chi':last}
    ranked=sorted(best.values(),key=lambda r:(-abs(Fraction(r['splice_score'])),r['input_difference']['words'],r['core_output_mask']['words']))
    invt=inverse_function(transpose)
    for i,r in enumerate(ranked):
        r['id']=f'R5_sparse_{i:06d}'
        if i<128:r['full_output_mask']=encode(invt(tuple(int(w,16) for w in r['core_output_mask']['words'])))
    summary={'source_trails':valid,'rejected_sources':rejected,'distinct_prefix_orbits':len(prefixes),'prefixes':prefixes,
             'middle_nonzero_matches':matches,'unique_endpoints':len(ranked),'top':ranked[:128],
             'limitation':'Exact restricted-prefix contributions can cancel in the unrestricted prefix. Tail multiplication is a heuristic, not a proved full correlation.'}
    (OUT/'candidates.json').write_text(json.dumps(summary,indent=2)+'\n')
    selected=[];seen_inputs={};seen_outputs=set()
    for r in ranked:
        a=tuple(r['input_difference']['words']);b=tuple(r['core_output_mask']['words'])
        if seen_inputs.get(a,0)>=2 or b in seen_outputs:continue
        selected.append(r);seen_inputs[a]=seen_inputs.get(a,0)+1;seen_outputs.add(b)
        if len(selected)==16:break
    for r in selected:r['full_output_mask']=encode(invt(tuple(int(w,16) for w in r['core_output_mask']['words'])))
    (OUT/'selection.json').write_text(json.dumps(selected,indent=2)+'\n')
    with (OUT/'jobs.txt').open('x') as f:
        for r in selected:f.write(f'{r["id"]} 5 '+ ';'.join(r['input_difference']['words'])+' '+';'.join(r['core_output_mask']['words'])+'\n')
    print('source trails',valid,'prefix orbits',len(prefixes),'rank certificates',sum(p['certificate']['exact_product_probability_certified'] for p in prefixes))
    print('matches',matches,'endpoints',len(ranked),'selected',len(selected))
    print([(r['id'],r['splice_score'],r['prefix_id']) for r in selected[:4]])
if __name__=='__main__':main()
