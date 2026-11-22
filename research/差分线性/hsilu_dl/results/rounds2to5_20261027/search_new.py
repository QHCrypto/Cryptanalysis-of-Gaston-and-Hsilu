"""Search new five-chi splits (2,2,1) and (1,2,2), without full endpoints."""
import json,math,sys,time
from pathlib import Path
from fractions import Fraction
from middle2 import two,affine,ROOT
from hsilu import linear,transpose,embed,rol,words,encode
from pool_search import LAT,decode,chi5
H=Path(__file__).resolve().parent;BASE=ROOT/'hsilu_dl/results/baseline'
DDT=[[sum(chi5(x)^chi5(x^a)==b for x in range(32)) for b in range(32)] for a in range(32)]
def save(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
def rec(q):return dict(numerator=str(q.numerator),denominator=str(q.denominator),log2_abs=math.log2(abs(q)) if q else None)
def main():
    started=time.perf_counter();pool=json.loads((BASE/'linear_pool.json').read_text())
    prefixes=json.loads((ROOT/'hsilu_dl/results/sparse_diff/candidates.json').read_text())['prefixes']
    heads=[]
    for lane in range(5):
        v=1<<lane
        for u in range(1,32):
            if LAT[u][v]:
                heads.append(dict(id=f'L1_{u}_{v}',input=words(embed(0,u)),output=words(embed(0,v)),
                    gamma=transpose(embed(0,u)),q2=Fraction(LAT[u][v]**2,1024)))
    kept=[];tested=nonzero=0
    for pref in prefixes:
        assert pref['certificate']['exact_product_probability_certified']
        before=linear(decode(pref['states']['y1']));din=decode(pref['states']['x0']);p=Fraction(1,2**sum(pref['weights']))
        for rot in range(64):
            delta=tuple(rol(w,rot) for w in before)
            for h in heads:
                cm=two(delta,h['gamma']);tested+=1
                if not cm:continue
                nonzero+=1;score=p*cm*h['q2']
                kept.append(dict(family='D2_M2_L1',split=[2,2,1],input_difference=words(tuple(rol(w,rot) for w in din)),
                    output_mask=h['output'],middle_input=words(delta),middle_mask=words(h['gamma']),
                    p=rec(p),cm_hnw_zero=rec(cm),q2=rec(h['q2']),score=rec(score),prefix_id=pref['id'],rotation=rot,tail={k:v for k,v in h.items() if k not in ['gamma','q2']}))
        print('D2_M2_L1',pref['id'],'tested',tested,'nonzero',nonzero,'seconds',round(time.perf_counter()-started,1),flush=True)
    kept.sort(key=lambda r:-abs(Fraction(int(r['score']['numerator']),int(r['score']['denominator']))))
    doc=dict(status='completed',family='D2_M2_L1',tested=tested,nonzero=nonzero,linear_heads=len(heads),
             scope='10 certified two-round differential prefixes × 64 rotations × all one-bit output one-chi linear trails at column 0',
             top=kept[:128],elapsed_seconds=time.perf_counter()-started)
    save(H/'search_221.json',doc)
    # Enumerate middle deltas induced by EVERY first-column DDT output,
    # independently of the first input pattern. This factors the search.
    retained=[];tests=hits=0;start=time.perf_counter()
    final={b['id']:min(b['last_chi_choices'],key=lambda z:(z['weight'],z['output']['words'])) for b in pool['bridges']}
    hs=[(h,decode(h['join_mask']),final[h['bridge_id']],Fraction(1,2**(2*(h['first_chi']['weight']+final[h['bridge_id']]['weight'])))) for h in pool['heads']]
    for col in range(64):
        for b in range(1,32):
            delta=linear(embed(col,b))
            for h,gamma,last,q2 in hs:
                cm=two(delta,gamma);tests+=1
                if not cm:continue
                hits+=1
                # For a fixed DDT output choose a maximum-probability nonzero
                # input, retaining one representative (tie: fewer bits).
                a=min(range(1,32),key=lambda a:(-DDT[a][b],a.bit_count(),a))
                p=Fraction(DDT[a][b],32);score=p*cm*q2
                retained.append(dict(family='D1_M2_L2',split=[1,2,2],input_difference=words(embed(col,a)),output_mask=last['output']['words'],
                     middle_input=words(delta),middle_mask=words(gamma),p=rec(p),cm_hnw_zero=rec(cm),q2=rec(q2),score=rec(score),
                     differential_output=words(embed(col,b)),head=h,last_chi=last))
        if col%8==7:print('D1_M2_L2','columns',col+1,'tested',tests,'nonzero',hits,'seconds',round(time.perf_counter()-start,1),flush=True)
    retained.sort(key=lambda r:-abs(Fraction(int(r['score']['numerator']),int(r['score']['denominator']))))
    save(H/'search_122.json',dict(status='completed',tested=tests,nonzero=hits,
         scope='64 columns × 31 nonzero first-chi output patterns × 384 existing linear-tail heads; best first DDT input and minimum-weight last completion retained',
         top=retained[:128],elapsed_seconds=time.perf_counter()-start))
    print('DONE new five-round matches',nonzero,hits,flush=True)
if __name__=='__main__':main()
