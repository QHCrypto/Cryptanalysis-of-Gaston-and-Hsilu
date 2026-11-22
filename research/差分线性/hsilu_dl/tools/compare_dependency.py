#!/usr/bin/env python3
"""Compare normalized three-chi cores and certify constant-derivative masks.

For d_a(X)=F(X)^F(X^a), a constant parity lambda.d_a(X) implies
lambda is orthogonal to every d_a(X)^d_a(0). Finite rank 320 therefore
rigorously excludes EVERY nonzero constant-parity mask for that a. The
samples here are witnesses, not a statistical confidence argument.
"""
import argparse
import hashlib
import json
from pathlib import Path
from hsilu import core, chi, pack, embed, rol, MASK, encode
from audit_structure import dependency, cone_chi
from linear_seed_scan import linear as gaston_linear, EAST, WEST, TWIN

ROOT=Path(__file__).resolve().parents[1]
SEED=2026101840


def mix(z):
    z=((z^(z>>30))*0xbf58476d1ce4e5b9)&MASK
    z=((z^(z>>27))*0x94d049bb133111eb)&MASK
    return z^(z>>31)


def sample(k):
    return tuple(mix((SEED+(5*k+i+1)*0x9e3779b97f4a7c15)&MASK) for i in range(5))


def gaston_core(s):
    s=chi(s)
    for rc in (0x5a,0x4b):
        s=gaston_linear(s)
        s=chi((s[0]^rc,*s[1:]))
    return s


def union(values):
    result=0
    for v in values:result|=v
    return result


def gaston_cone(s):
    stages=[]
    for r in range(3):
        if r:
            a=tuple(rol(v,n) for v,n in zip(s,EAST))
            p=union(a);q=union(rol(v,n) for v,n in zip(a,TWIN))
            effect=rol(p|rol(p,1)|q|rol(q,18),23)
            s=tuple(rol(v|effect,n) for v,n in zip(a,WEST))
        s=cone_chi(s);stages.append(encode(s))
    return stages


def insert(value,pivots):
    while value:
        p=value.bit_length()-1
        if p not in pivots:
            pivots[p]=value
            return True
        value^=pivots[p]
    return False


def independent_rank(rows):
    # Recheck each certificate by elimination from the opposite end.
    pivots={}
    for v in rows:
        while v:
            low=v&-v
            if low not in pivots:
                pivots[low]=v
                break
            v^=pivots[low]
    return len(pivots)


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    states=[sample(k) for k in range(1024)]
    result={'status':'PASS','core_boundary':'chi -> L -> iota_5a -> chi -> L -> iota_4b -> chi',
            'bit_numbering':'lane 0..4, LSB-index 0..63','seed':SEED,
            'witness_generator':'SplitMix64 finalizer of seed+(5*k+i+1)*golden_ratio mod 2^64',
            'rank_definition':'span over GF(2) of d_delta(X) xor d_delta(0)',
            'inference':'A nonzero constant output parity must lie in the orthogonal complement of this span.',
            'scope':'All 320 one-bit differences at the first-chi input; all nonzero 320-bit output masks; fixed last-three-round constants',
            'not_claimed':'No conclusion about multi-bit input differences, other RC schedules, or a bound away from |rho|=1.',
            'cone_counts_by_lane':[], 'algorithms':{}}
    for lane in range(5):
        delta=embed(0,1<<lane)
        result['cone_counts_by_lane'].append({'lane':lane,
             'Hsilu':[r['bit_weight'] for r in dependency(delta,3)],
             'Gaston':[r['bit_weight'] for r in gaston_cone(delta)]})
    for name,F in [('Hsilu',lambda s:core(s,3)),('Gaston',gaston_core)]:
        base=[pack(F(x)) for x in states];f0=pack(F((0,)*5));records=[]
        for lane in range(5):
            for col in range(64):
                delta=embed(col,1<<lane);d0=pack(F(delta))^f0
                target_rank=318 if name=='Hsilu' and lane==3 else 320
                pivots={};witnesses=[];vectors=[]
                for k,(x,y) in enumerate(zip(states,base)):
                    v=pack(F(tuple(a^b for a,b in zip(x,delta))))^y^d0
                    if insert(v,pivots):witnesses.append(k);vectors.append(v)
                    if len(pivots)==target_rank:break
                assert len(pivots)==target_rank,(name,lane,col,len(pivots))
                assert independent_rank(vectors)==target_rank
                known=[]
                if target_rank==318:
                    # Upper bound is a separate all-input dependency proof,
                    # not a conclusion from finite rank deficiency.
                    cone=dependency(delta,3)[-1]['words']
                    for outlane in (0,1):
                        j=(col+29)%64
                        assert not (int(cone[outlane],16)>>j&1)
                        assert not (d0>>(64*outlane+j)&1)
                        known.append([outlane,j])
                records.append({'input_bit':[lane,col],'affine_rank':target_rank,
                    'witness_indices':witnesses,'known_zero_derivative_output_bits':known,
                    'constant_parity_mask_dimension':320-target_rank})
            print(name,'lane',lane,'completed',flush=True)
        histogram={str(r):sum(x['affine_rank']==r for x in records) for r in sorted({x['affine_rank'] for x in records})}
        result['algorithms'][name]={'rank_histogram':histogram,'records':records}
    sources=[Path(__file__),ROOT/'tools/hsilu.py',ROOT/'tools/audit_structure.py',
             ROOT.parent/'gaston_dl/tools/linear_seed_scan.py']
    result['source_sha256']={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with args.output.open('x') as stream:json.dump(result,stream,ensure_ascii=False,indent=2);stream.write('\n')
    print('PASS', {name:data['rank_histogram'] for name,data in result['algorithms'].items()})


if __name__=='__main__':main()
