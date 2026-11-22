"""Prepend a chi to each existing two-chi linear tail. Exact affine LAT matching.
One-column DDT front, exact one-chi DLCT middle, three-chi single linear trail.
"""
import json,sys,time,math
from pathlib import Path
from fractions import Fraction
from search_new import H,BASE,DDT,save,rec
from hsilu import linear,transpose,embed,chi,words,encode
from pool_search import LAT,decode

def xor(a,b):return tuple(x^y for x,y in zip(a,b))
def support(v):
    us=[u for u in range(32) if LAT[u][v]];off=us[0];piv={}
    for u in us:
        z=u^off
        while z:
            k=z.bit_length()-1
            if k in piv:z^=piv[k]
            else:piv[k]=z;break
    basis=list(piv.values());span={0}
    for z in basis:span|={x^z for x in list(span)}
    assert {off^x for x in span}==set(us)
    assert len({abs(LAT[u][v]) for u in us})==1
    return off,basis,int(round(math.log2(32/abs(LAT[us[0]][v]))))
SUP=[support(v) for v in range(32)]

def setup(h):
    B=decode(h['join_mask']);u0=(0,)*5;ubs=[];weight=0
    for j in range(64):
        v=sum((B[i]>>j&1)<<i for i in range(5));off,bs,w=SUP[v]
        u0=xor(u0,embed(j,off));ubs.extend(embed(j,b) for b in bs);weight+=w
    gs=[transpose(u0)]+[transpose(b) for b in ubs]
    gforms=[sum((g[i]>>j&1)<<k for k,g in enumerate(gs)) for i in range(5) for j in range(64)]
    return B,u0,ubs,gforms,weight

def solve(delta,forms):
    piv={}
    for j in range(64):
        d=[delta[i]>>j&1 for i in range(5)]
        if not any(d):continue
        for i in range(5):
            eq=(forms[((i-1)%5)*64+j] if d[(i+1)%5] else 0)^(forms[((i-2)%5)*64+j] if d[(i-1)%5] else 0)
            while eq&~1:
                k=eq.bit_length()-1
                if k in piv:eq^=piv[k]
                else:piv[k]=eq;break
            else:
                if eq:return None
    # Bit 0 represents the fixed constant 1, not a free variable.
    assignment=1
    for k in sorted(piv):
        if ((piv[k]^(1<<k))&assignment).bit_count()&1:assignment|=1<<k
    return assignment

def main():
    t=time.perf_counter();pool=json.loads((BASE/'linear_pool.json').read_text())
    final={b['id']:min(b['last_chi_choices'],key=lambda z:(z['weight'],z['output']['words'])) for b in pool['bridges']}
    heads=sorted(pool['heads'],key=lambda h:2*(setup(h)[-1]+h['first_chi']['weight']+final[h['bridge_id']]['weight']))
    ds=[(col,b,linear(embed(col,b))) for col in range(64) for b in range(1,32)]
    total=hits=0;kept=[];perhead=[]
    for hi,h in enumerate(heads):
        B,u0,ubs,forms,w=setup(h);last=final[h['bridge_id']];qw=2*(w+h['first_chi']['weight']+last['weight']);found=0
        for col,b,delta in ds:
            total+=1;s=solve(delta,forms)
            if s is None:continue
            hits+=1;found+=1
            u=u0
            for k,ub in enumerate(ubs,1):
                if s>>k&1:u=xor(u,ub)
            gamma=transpose(u);assert solve(delta,[((gamma[i]>>j)&1) for i in range(5) for j in range(64)])==1
            sign=-1 if sum((g&c).bit_count() for g,c in zip(gamma,chi(delta)))%2 else 1
            a=min(range(1,32),key=lambda a:(-DDT[a][b],a.bit_count(),a));p=Fraction(DDT[a][b],32);q2=Fraction(1,2**qw)
            r=dict(split=[1,1,3],head_id=h['id'],input_difference=words(embed(col,a)),output_mask=last['output']['words'],differential_output=words(embed(col,b)),
                   middle_input=words(delta),middle_mask=words(gamma),prepended_chi=dict(input=words(u),output=words(B),weight=w),head=h,last_chi=last,
                   p=rec(p),cm=rec(Fraction(sign)),q2=rec(q2),score=rec(p*sign*q2))
            kept.append(r)
        perhead.append(dict(id=h['id'],q_squared_weight=qw,hits=found))
        if found or hi%16==15:print('L3 affine',hi+1,h['id'],'q2weight',qw,'hits',found,'totalhits',hits,'sec',round(time.perf_counter()-t,2),flush=True)
    kept.sort(key=lambda r:-abs(Fraction(int(r['score']['numerator']),int(r['score']['denominator']))))
    unique=[];seen=set()
    for r in kept:
        k=(tuple(r['input_difference']),tuple(r['output_mask']))
        if k not in seen:seen.add(k);unique.append(r)
    save(H/'search_113.json',dict(status='completed',scope='all 384 old two-chi heads extended backwards by all nonzero LAT supports, via affine feasibility; all 64 × 31 one-column DDT outputs; one deterministic free-variable assignment per feasible match',tested=total,nonzero=hits,distinct_endpoints=len(seen),perhead=perhead,top=unique[:128],elapsed_seconds=time.perf_counter()-t))
    print('DONE',total,hits,len(seen),flush=True)
if __name__=='__main__':main()
