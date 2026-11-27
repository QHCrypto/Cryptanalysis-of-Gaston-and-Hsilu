"""HNW two-round prefix weighting with basic and joint final-layer readouts.

We integrate the prefix difference law by affine elimination, rather than
enumerating every difference. The tail is chi3-L-chi4 with fresh uniform
base values. The basic Algorithm 1 readout multiplies final-column expectations.
The joint readout follows the authors' PRESENT getBias_Opt2/getBias_Opt4
implementation, combining the active-column masks before the transpose and
the final squared-LAT evaluation. Both use the same prefix distribution below.
"""
from joint_estimator import *
from pathlib import Path

HNW_IMPLEMENTATION = {
    'repository': 'https://github.com/hukaisdu/Round-based-approximation',
    'commit': '783e8edc560b848d0f2ba63f11838aecd23da40f',
    'file': 'PRESENT/present.py',
    'functions': ['getBias_Opt2', 'getBias_Opt4'],
    'source_url': 'https://github.com/hukaisdu/Round-based-approximation/blob/783e8edc560b848d0f2ba63f11838aecd23da40f/PRESENT/present.py',
    'adaptation': 'Joint final-layer treatment adapted to chi5 and the Gaston/Hsilu mixers; affine elimination evaluates the two-round differential-prefix sum.',
}

def compare(name,prefix='exact',split_output=False):
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
            'scope':'HNW Eq. (5) with two-round prefix and basic Algorithm 1 column-product tail readout, zero value masks; exact arithmetic for the stated model.'
                    if split_output else 'HNW two-round differential-prefix weighting with the joint final-layer treatment from the authors\' public PRESENT implementation, adapted to Gaston/Hsilu; zero value masks and complete L^T masks; exact arithmetic for the stated model, not a fixed-permutation correlation certificate.'}

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
