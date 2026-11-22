#!/usr/bin/env python3
"""Exhaust the exact nullspaces for every support passing both rank checks."""
from functools import reduce
import hashlib
import json
from operator import or_
from pathlib import Path
from audit_kernel_structure import EAST,WEST,rol,maps
from verify_lossy_model import COST,column

ROOT=Path(__file__).resolve().parents[1]
matrix_path=ROOT/'results/support_reference_matrices_20261024.txt'
raw=[int(value,16) for value in matrix_path.read_text().split()]
assert len(raw)==640
tables=[raw[:320],raw[320:]]


def nullspace(columns):
    basis={};kernel=[]
    for j,value in enumerate(columns):
        coefficients=1<<j
        while value:
            pivot=value.bit_length()-1
            if pivot in basis:
                value^=basis[pivot][0];coefficients^=basis[pivot][1]
            else:
                basis[pivot]=(value,coefficients)
                break
        if not value:
            kernel.append(coefficients)
    return kernel


def chi_weight(words):
    return sum(COST[column(words,z).bit_count()] for z in range(64))


def support(words):
    return reduce(or_,words)


def candidates(mask,side):
    positions=[z for z in range(64) if mask>>z&1]
    columns=[tables[side][5*z+i] for z in positions for i in range(5)]
    basis=nullspace(columns)
    # Abort rather than silently truncate if a large nullspace appears.
    if len(basis)>24:
        raise RuntimeError('Unexpected large nullspace needs separate enumeration')
    valid=[]
    for selector in range(1,1<<len(basis)):
        v=0
        for j,b in enumerate(basis):
            if selector>>j&1:v^=b
        syndrome=0
        for j,c in enumerate(columns):
            if v>>j&1:syndrome^=c
        assert syndrome==0
        words=[sum(((v>>(5*j+i))&1)<<z for j,z in enumerate(positions)) for i in range(5)]
        theta=[rol(w,(-WEST[i] if side==0 else EAST[i])) for i,w in enumerate(words)]
        assert maps(theta)[3]==0
        if support(words)!=mask:
            continue
        other=[rol(w,(-1 if side==0 else 1)*(EAST[i]+WEST[i])) for i,w in enumerate(words)]
        cost=2*support(other).bit_count()+chi_weight(words) if side==0 else chi_weight(other)
        valid.append({'middle_words':[f'{w:016x}' for w in words],
                      'other_words':[f'{w:016x}' for w in other],
                      'bit_weight':sum(w.bit_count() for w in words),
                      'other_active_columns':support(other).bit_count(),'partial_lossy_cost':cost})
    return {'rank':len(columns)-len(basis),'nullity':len(basis),
            'kernel_basis':[f'{v:016x}' for v in basis],
            'nonzero_vectors_exhausted':(1<<len(basis))-1,
            'full_support_candidates':valid}


def main():
    branches={}
    for k,name in ((6,'six'),(7,'seven'),(8,'eight')):
        probe_path=ROOT/f'results/support_rank_{name}_20261024.json'
        probe=json.loads(probe_path.read_text())
        records=[];minimum=None;best=None;full=0
        for mask_hex in probe['both_orbit_masks']:
            mask=int(mask_hex,16)
            assert mask.bit_count()==k
            left,right=candidates(mask,0),candidates(mask,1)
            record={'mask':mask_hex,'left':left,'right':right,'minimum_lossy_weight':None}
            if left['full_support_candidates'] and right['full_support_candidates']:
                full+=1
                a=min(left['full_support_candidates'],key=lambda v:v['partial_lossy_cost'])
                b=min(right['full_support_candidates'],key=lambda v:v['partial_lossy_cost'])
                weight=a['partial_lossy_cost']+b['partial_lossy_cost']
                record['minimum_lossy_weight']=weight
                if minimum is None or weight<minimum:
                    minimum=weight;best={'mask':mask_hex,'left':a,'right':b,'lossy_weight':weight}
            records.append(record)
        branches[str(k)]={'rank_candidate_orbits':len(records),'both_full_support_orbits':full,
                          'minimum_lossy_weight':minimum,'minimum_witness':best,
                          'below_44_excluded':minimum is None or minimum>=44,'records':records,
                          'rank_enumeration_sha256':hashlib.sha256(probe_path.read_bytes()).hexdigest()}
    report={'status':'PASS','scope':'Exact lossy double-kernel support-space analysis; actual DDT validity not assumed',
            'matrix_sha256':hashlib.sha256(matrix_path.read_bytes()).hexdigest(),
            'objective_identity':'min_L [2*cols(T^-1*a1)+chi_cost(a1)] + min_R chi_cost(T*b1), with both middle supports exactly S',
            'branches':branches}
    output=ROOT/'results/support_space_analysis_20261024.json'
    output.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:{name:v[name] for name in ('rank_candidate_orbits','both_full_support_orbits','minimum_lossy_weight','below_44_excluded')} for k,v in branches.items()},indent=2))


if __name__=='__main__':main()
