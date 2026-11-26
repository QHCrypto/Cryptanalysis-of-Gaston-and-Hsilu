"""HNW two-round prefix weighting and a separate joint-mask extension.

We integrate the prefix difference law by affine elimination, rather than
enumerating every difference. The tail is chi3-L-chi4 with fresh uniform
base values. For each fixed prefix difference, HNW's column readout multiplies
the expectations for the final active columns. This is different from taking
one expectation for their joint parity. Both are integrated against exactly
the same prefix distribution below.
"""
from joint_estimator import *
from pathlib import Path

def compare(name,prefix='exact',split_output=True):
    started = time.perf_counter()
    L,LT = (lg,lgt) if name=='Gaston' else (lh,lht)
    delta,target = ENDPOINTS[name]
    j0 = next(j for j in range(64) if column(delta,j))
    free = ((1 << 320)-1)^pack(embed(j0,31))
    outcols = [j for j in range(64) if column(target,j)]
    total = Fraction(0)
    records=[]
    for choices in product(*(SUPPORT[column(target,j)] for j in outcols)):
        weight = Fraction(1)
        pulled = []
        for j,a in zip(outcols,choices):
            weight *= Fraction(LAT[column(target,j)][a]**2,1024)
            pulled.append(LT(embed(j,a)))
        if not split_output:
            pulled = [add(*pulled)]
        # Split mode represents PRODUCT of chi3 autocorrelations, one for each
        # final output column. Joint mode represents the autocorrelation of
        # the XOR of the pulled masks. Each product is then averaged jointly
        # over the entire prefix difference distribution.
        constraints = [embed(j,v) for b in pulled for j in range(64)
                       for v in BASES[column(b,j)]]
        origin = ZERO
        for b in pulled:
            for j in range(64):
                origin = add(origin,embed(j,SUPPORT[column(b,j)][0]))
        mean = Fraction(0)
        nonzero=[]
        for x0 in range(32):
            d1 = embed(j0,sbox(x0)^sbox(x0^column(delta,j0)))
            d2 = L(d1)
            def adj(k):
                a = polar_transpose(LT(k),d2)
                return pack(a) if prefix=='fresh' else pack(LT(a)) & free
            y = ZERO if prefix=='fresh' else add(L(embed(j0,sbox(x0))),(RC[9],0,0,0,0))
            off = L(add(chi(y),chi(add(y,d2))))
            val,rank,reason = affine_character([adj(k) for k in constraints],
                                              [dot(k,off) for k in constraints],
                                              adj(origin),dot(origin,off))
            mean += val/32
            if val:
                nonzero.append({'first_column_value':x0,'value':str(val),'rank':rank})
        total += weight*mean
        records.append({'last_input_column_patterns':choices,'weight':str(weight),
                        'conditional_mean':str(mean),'nonzero':nonzero})
    return {'target':name,'prefix':prefix,'split_final_column_readout':split_output,
            'correlation':fraction_record(total),'branches':records,
            'seconds':time.perf_counter()-started,
            'scope':'HNW Eq. (5) with two-layer prefix and Algorithm 1 column-product tail readout, zero value masks; exact arithmetic for stated model.'
                    if split_output else 'Joint-mask extension, not original HNW26: two-round differential-prefix weighting with zero value masks and complete L^T masks in the tail; exact arithmetic for the stated approximation, not a fixed-permutation correlation certificate.'}

if __name__=='__main__':
    records=[]
    for name in ENDPOINTS:
        for prefix in ['fresh','exact']:
            for split in [True,False]:
                r=compare(name,prefix,split)
                if not split:
                    assert r['correlation']==estimate(name,prefix)['correlation']
                records.append(r)
                print(name,prefix,'split' if split else 'joint',r['correlation'],flush=True)
    Path('hnw_prefix_comparison.json').write_text(json.dumps(records,indent=2)+'\n')
