"""Exact evaluation of the HNW26 zero-value-mask model for 1/2 chi layers.

For two chi layers this computes precisely
 sum_b Pr[D_delta chi(X)=b] prod_j DLCT((L b)_j, gamma_j).
Because chi is quadratic, its derivative is affine in X. Each last-layer
DLCT is either zero (linear equations on Lb) or a quadratic sign. Solve the
equations and evaluate that quadratic sign exactly. This changes only the
summation order of the first-DDT-enumerated HNW zero model, not the model.
No full endpoint and no enhanced DATF computation is performed here.
"""
import json,sys,math,time,subprocess,random
from fractions import Fraction
from pathlib import Path
H=Path(__file__).resolve().parent;G=H.parents[1]
sys.path.insert(0,str(G/'tools'))
from fixed_tail_search import S,column,embed,linear,transpose,state_chi,lat
from datf_quadratic import uniform_quadratic

def bits(n):
    while n:
        b=n&-n;yield b;n^=b
def record(q):return dict(numerator=str(q.numerator),denominator=str(q.denominator),log2_abs=math.log2(abs(q.numerator))-math.log2(q.denominator) if q else None)
def fraction(q):return Fraction(int(q['numerator']),int(q['denominator']))
def dlct(a,v):return Fraction(sum((-1)**((v&(S[x]^S[x^a])).bit_count()&1) for x in range(32)),32)
def one(delta,gamma):return math.prod(dlct(column(delta,j),column(gamma,j)) for j in range(64))
def two(delta,gamma):
    # Affine encoding: bit 0 is the constant, bit k+1 is input variable k.
    base=linear(state_chi(delta))
    d=[[((base[i]>>j)&1) for j in range(64)] for i in range(5)]
    for k in range(5):
        for j in range(64):
            pattern=(((delta[(k+1)%5]>>j)&1)<<((k-1)%5))^(((delta[(k-1)%5]>>j)&1)<<((k-2)%5))
            if not pattern:continue
            image=linear(embed(j,pattern));var=1<<(k*64+j+1)
            for i in range(5):
                for b in bits(image[i]):d[i][b.bit_length()-1]^=var
    pivots={}
    for j in range(64):
        v=column(gamma,j)
        for k in range(5):
            eq=(d[(k+1)%5][j] if (v>>((k-1)%5))&1 else 0)^(d[(k-1)%5][j] if (v>>((k-2)%5))&1 else 0)
            while eq&~1:
                pivot=1<<(eq.bit_length()-1)
                if pivot in pivots:eq^=pivots[pivot]
                else:pivots[pivot]=eq;break
            else:
                if eq==1:return Fraction(0),dict(inconsistent_constraints=True,rank=len(pivots))
    substitutions={}
    def subst(f):
        out=f&1
        for b in bits(f&~1):out^=substitutions[b]
        return out
    for k in range(1,321):
        bit=1<<k;substitutions[bit]=subst(pivots[bit]^bit) if bit in pivots else bit
    d=[[subst(f) for f in row] for row in d]
    def poly(f):return set([0] if f&1 else [])|set(bits(f&~1))
    def multiply(a,b):
        out=set()
        for x in poly(a):
            for y in poly(b):
                term=x|y
                if term in out:out.remove(term)
                else:out.add(term)
        return out
    p=set()
    for j in range(64):
        for i in range(5):
            if gamma[i]>>j&1:
                p.symmetric_difference_update(poly(d[i][j]^d[(i+2)%5][j]))
                p.symmetric_difference_update(multiply(d[(i+1)%5][j],d[(i+2)%5][j]))
    value=uniform_quadratic(frozenset(p))/2**len(pivots)
    return value,dict(rank=len(pivots),quadratic_terms=len(p),inconsistent_constraints=False)

def full_zero_certificate(delta):
    inactive={j for j in range(64) if not column(delta,j)}
    killed=0
    for j in range(64):
        for a in range(1,32):
            mask=transpose(embed(j,a))
            if any(column(mask,k) for k in inactive):killed+=1
    return dict(nonzero_base_value_masks_killed=killed,total=64*31,
        full_equals_zero_certified_for_two_chi=killed==64*31)

def main():
    jobs=json.loads((H/'manifest.json').read_text())['jobs'];testdir=H/'raw/checks';testdir.mkdir(exist_ok=True)
    # All one-chi LAT-square identities, independent of the DATF worker.
    for a in range(32):
        for v in range(32):
            assert dlct(a,v)==sum(Fraction(lat(u,v)**2,1024)*((-1)**((u&a).bit_count()&1)) for u in range(32))
    rng=random.Random(20261027);checks=[]
    for i in range(8):
        delta=[0]*5
        for j in rng.sample(range(64),1+i%2):
            a=rng.choice([1,4,8,13,31])
            for k in range(5):delta[k]^=((a>>k)&1)<<j
        gamma=[0]*5
        for j in rng.sample(range(64),1+i%3):
            a=rng.randrange(1,32)
            for k in range(5):gamma[k]^=((a>>k)&1)<<j
        if i==7:gamma=[int(w,16) for w in next(j for j in jobs if j['rounds']==2)['lambda_core']]
        out=testdir/f'case{i}.json'
        subprocess.run([str(G/'build/hnw_endpoint'),'zero','2','8',*[f'{w:016x}' for w in delta+gamma],str(out)],check=True,stdout=subprocess.DEVNULL)
        native=json.loads(out.read_text());fast,stats=two(delta,gamma)
        assert fast==fraction(native['correlation']),(i,fast,native['correlation'])
        checks.append(dict(case=i,native_first_DDT_branches=native['completed_branches'],correlation=record(fast),statistics=stats))
    for j in jobs:
        started=time.perf_counter();delta=tuple(int(w,16) for w in j['delta_core']);gamma=tuple(int(w,16) for w in j['lambda_core'])
        if j['rounds']==1:q=one(delta,gamma);stats={};cert=dict(full_equals_zero_certified_for_one_chi=True)
        else:q,stats=two(delta,gamma);cert=full_zero_certificate(delta)
        result=dict(**j,status='completed',method='HNW26 zero-value-mask model with exact first-DDT summation',
            backend='local DLCT product' if j['rounds']==1 else 'affine constraints and exact quadratic sign sum',
            correlation=record(q),statistics=stats,certificate=cert,elapsed_seconds=time.perf_counter()-started)
        (H/'raw'/f"{j['id']}_hnw.json").write_text(json.dumps(result,indent=2)+'\n')
    (H/'hnw_validation.json').write_text(json.dumps(dict(status='PASS',one_chi_LAT_identities=1024,native_two_chi_crosschecks=checks),indent=2)+'\n')
    print('HNW completed',len(jobs),'middle jobs; 8 native first-DDT-enumeration cross-checks passed.')
if __name__=='__main__':main()
