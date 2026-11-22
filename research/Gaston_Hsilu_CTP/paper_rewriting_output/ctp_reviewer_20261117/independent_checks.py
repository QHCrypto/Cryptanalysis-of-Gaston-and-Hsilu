from pathlib import Path
import re, math, json, random, collections

TEX=(Path(__file__).resolve().parents[3] / 'Gaston_Hsilu_CTP/verification_inputs/paper_parameters.tex').read_text()
MASK=(1<<64)-1
RC=[0xf0,0xe1,0xd2,0xc3,0xb4,0xa5,0x96,0x87,0x78,0x69,0x5a,0x4b]
def rot(x,n):
    n%=64
    return ((x<<n)|(x>>(64-n)))&MASK
def chi(x): return [x[i]^((~x[(i+1)%5])&x[(i+2)%5]) for i in range(5)]
def lh(x):
    z=list(x)
    for a,b,r in [(0,4,47),(3,4,51),(1,2,54),(0,1,30),(5,2,4),(5,3,55),(7,6,0),(11,8,14),(12,9,59),(12,6,4),(13,10,10),(14,11,45),(15,9,53)]: z.append(z[a]^rot(z[b],r))
    return [rot(z[13+i],r) for i,r in enumerate([0,57,13,31,20])]
def effect(a):
    p=0;q=0
    for x,t in zip(a,[25,32,52,60,63]): p^=x;q^=rot(x,t)
    return p^rot(p,1)^q^rot(q,18)
def lg(x):
    a=[rot(v,r) for v,r in zip(x,[0,60,22,27,4])]; e=rot(effect(a),23)
    return [rot(v^e,r) for v,r in zip(a,[0,56,31,46,43])]
def perm(x,r,kind='h'):
    x=list(x)
    for c in RC[-r:]:
        if kind=='g': x=lg(x)
        x[0]^=c;x=chi(x)
        if kind=='h': x=lh(x)
    return x
def stages(x,r):
    out=[]
    for c in RC[-r:]:
        x=list(x);x[0]^=c;x=lh(chi(x));out.append(x)
    return out
def prefix(m0,m1,r,kind='h'):
    x=perm([0x0000080000cc0003,0,0,0,0],r,kind)
    for m in [m0,m1]:
        x[0]^=m;x=perm(x,r,kind)
    return x
def solve(rows,rhs,n):
    a=[x|(b<<n) for x,b in zip(rows,rhs)];p=[];m=(1<<n)-1
    for j in range(n):
        k=next((k for k in range(len(p),len(a)) if (a[k]>>j)&1),None)
        if k is None: continue
        i=len(p);a[i],a[k]=a[k],a[i]
        for k in range(len(a)):
            if k!=i and (a[k]>>j)&1:a[k]^=a[i]
        p.append(j)
    if any((x&m)==0 and x>>n for x in a):return None,[],len(p)
    base=sum(((a[i]>>n)&1)<<j for i,j in enumerate(p))
    basis=[(1<<j)|sum(((a[i]>>j)&1)<<k for i,k in enumerate(p)) for j in range(n) if j not in p]
    return base,basis,len(p)
def comb(x,basis):
    z=0
    for j,v in enumerate(basis):
        if (x>>j)&1:z^=v
    return z
def embed(x,positions): return sum(((x>>i)&1)<<j for i,j in enumerate(positions))
def cols(x): return [(sum(((v>>j)&1)<<i for i,v in enumerate(x))) for j in range(64)]
def state(label):
    line=next(l for l in TEX.splitlines() if l.startswith('$'+label+'$ &'))
    return [int(v,16) for v in re.findall(r'\\hex\{([0-9a-f]+)\}',line)]

def main():
    report={};rng=random.Random(20261117)
    box=[sum(v<<i for i,v in enumerate(chi([(a>>i)&1 for i in range(5)]))) for a in range(32)]
    ddt=[[sum(box[x]^box[x^a]==b for x in range(32)) for b in range(32)] for a in range(32)]
    weights=[];active=[]
    for i in range(3):
        x,y=state(f'x_{i}'),state(f'y_{i}')
        weights.append(sum(-math.log2(ddt[a][b]/32) for a,b in zip(cols(x),cols(y))))
        active.append(sum(a!=0 for a in cols(x)))
        if i<2:assert lh(y)==state(f'x_{i+1}')
    assert weights==[16,14,50] and active==[8,7,22]
    report['hsilu_trail']={'weights':weights,'active_columns':active,'linear_connections':'PASS'}
    assert all(-math.log2(ddt[a][b]/32)>=2+math.ceil((a.bit_count()-1)/2) for a in range(1,32) for b in range(32) if ddt[a][b])
    report['chi_input_weight_inequality']='checked every valid nonzero DDT transition'
    gx=[0,0x0000008000000100,0,0,1]
    gy=[0x0400000000002000,0x1004000008000020,0x0000100002000000,0x0800010000000000,0x0100802000000000]
    assert lg(gx)==gy
    assert effect([rot(v,r) for v,r in zip(gx,[0,60,22,27,4])])==0x0040000800000000
    assert [sum(a!=0 for a in cols(s)) for s in [gx,gy]]==[3,13]
    assert all(a==0 or a.bit_count()==1 for s in [gx,gy] for a in cols(s))
    report['gaston_weight32_witness']='PASS'
    def f(k):
        if k<=8:return [0,18,17,13,13,12,6,7,7][k]
        for hi,v in [(24,6),(33,5),(38,4),(47,3),(56,2),(64,1)]:
            if k<=hi:return v
    def g(k):
        if k<=5:return [0,57,48,39,34,25][k]
        for hi,v in [(11,6),(12,5),(16,3),(17,2),(64,1)]:
            if k<=hi:return v
    triples=[(u,v,w) for u in range(1,65) for v in range(1,65) for w in range(1,65) if v>=f(u) and u>=g(v) and w>=f(v) and v>=g(w) and u+v>=12 and v+w>=12]
    def outside(u,v):return u+v>=16 and (v!=6 or u>=13)
    report['profile_triple_minima']={str(pattern):min(sum(t) for t in triples if (not pattern[0] or outside(t[0],t[1])) and (not pattern[1] or outside(t[1],t[2]))) for pattern in [(1,0),(0,1),(1,1)]}
    report['both_kernel_below44_middle_counts']=sorted(set(t[1] for t in triples if sum(t)<22))
    report['middle9_below44_triples']=[t for t in triples if sum(t)<22 and t[1]==9]
    U={0,4,5,10,11,22,23,27,33,34,40,41,45,51,56,57,63}
    J=[0,1,6,17,18,19,24,29,30,41,42,47,52,53,59,60,63]
    gr=[sum(1<<((-j-d)%64) for d in [23,24,48,2]) for j in range(64) if j not in U]
    base,basis,rank=solve(gr,[0]*47,64);assert rank==47
    jm=sum(1<<((-j)%64) for j in J)
    def beforechi2(x):
        x=lg(x);x[0]^=RC[-2];x=lg(chi(x));x[0]^=RC[-1];return x[0]&jm
    for _ in range(8):
        s=prefix(rng.getrandbits(64),rng.getrandbits(64),2,'g')
        y=beforechi2(s);d=[]
        for v in basis:d.append(beforechi2([s[0]^v]+s[1:])^y)
        for i,v in enumerate(basis):
            for j in range(i):assert beforechi2([s[0]^v^basis[j]]+s[1:])==y^d[i]^d[j]
    report['gaston_shared_forms']={'rank':rank,'second_derivative_checks':8*136}
    A2=[0,1,2,4,6,9,10,14,19,23,26,28,32,38,40,45,47,49,51,55,57,59]
    J2=[0,1,2,5,6,7,10,14,17,19,25,33,34,38,42,43,50,52,56,57,61,62]
    m2=sum(1<<j for j in J2)
    for _ in range(12):
        x=[rng.getrandbits(64) for i in range(5)]
        for j in set(A2)-{10,55}:x[1]&=~(1<<j)
        for j in (set(A2)-{14,59})|{44}:x[4]|=1<<j
        if (stages(x,2)[0][2]>>54)&1:x[0]^=1<<27
        assert (stages(x,2)[0][2]>>54)&1==0
        y=perm(x,2)[0]&m2
        ds=[(perm([x[0]^(1<<j)]+x[1:],2)[0]&m2)^y for j in A2]
        for i,j in enumerate(A2):
            for k in range(i):assert perm([x[0]^(1<<j)^(1<<A2[k])]+x[1:],2)[0]&m2==y^ds[i]^ds[k]
    report['hsilu2_local_structure']={'conditioned_states':12,'second_derivative_checks':12*231,'legal_prefix_density':'not verified'}
    s=prefix(0x63c92946c6b0ac2c,0x101be7b317c13396,3)
    assert all((s[1]>>j)&1==0 and (s[4]>>j)&1==1 for j in [12,18])
    cond=[(0,[21,42,46,48,52],[19,25]),(1,[7,13,26],[20,35,39,45]),(2,[],[57,63]),(3,[26,35,41,45],[21,27,33,42,46,52]),(4,[57,63],[7,13,26,32])]
    c1=[(i,j,b) for i,zero,one in cond for b,js in [(0,zero),(1,one)] for j in js]
    c2=[(0,53,0),(0,48,1),(1,15,1),(1,52,1)]
    guesses=[(i,j) for i,js in enumerate([[4,11,18,22,31,53],[2,6,8,58,61],[15,34,48,53],[4,9,22,49,59],[5,54,55,58,59]]) for j in js]
    bg=[j for j in range(64) if j not in [12,18]]
    def hs(x):return stages([x]+s[1:],3)
    z=hs(0)[0];cd=[[(hs(1<<j)[0][i]^z[i]) for i in range(5)] for j in bg]
    rows=[sum(((cd[k][i]>>j)&1)<<k for k in range(62)) for i,j,b in c1]
    rhs=[b^((z[i]>>j)&1) for i,j,b in c1]
    b0,B,rk=solve(rows,rhs,62);assert rk==32 and b0 is not None
    z=hs(embed(b0,bg))[0]
    gr=[sum((((hs(embed(b0^v,bg))[0][i]^z[i])>>j)&1)<<k for k,v in enumerate(B)) for i,j in guesses]
    _,BB,grank=solve(gr,[0]*25,30);assert grank==23
    valid=0;branches=0
    for _ in range(32):
        tv=rng.getrandbits(30)
        bv=embed(b0^comb(tv,B),bg)
        BBx=[embed(comb(v,B),bg) for v in BB]
        zs=hs(bv)[1]
        rows2=[sum((((hs(bv^v)[1][i]^zs[i])>>j)&1)<<k for k,v in enumerate(BBx)) for i,j,b in c2]
        rhs2=[b^((zs[i]>>j)&1) for i,j,b in c2]
        a,bs,r=solve(rows2,rhs2,7)
        branches+=1
        if a is None:continue
        for v in range(1<<len(bs)):
            word=bv^comb(a^comb(v,bs),BBx)
            ys=[]
            for av in range(4):
                ss=hs(word^embed(av,[12,18]));ys.append(ss[-1][0]&((1<<5)|(1<<15)))
                assert all(((ss[0][i]>>j)&1)==b for i,j,b in c1)
                assert all(((ss[1][i]>>j)&1)==b for i,j,b in c2)
            assert ys[0]^ys[1]^ys[2]^ys[3]==0
            valid+=1
    report['hsilu3_concrete_structure']={'condition_rank':rk,'guess_rank':grank,'branches_sampled':branches,'backgrounds_checked':valid,'active_assignments_checked':4*valid}
    report['arithmetic']={'gaston_candidate_exponent':111-17*math.log2(.75),'hsilu2_candidate_exponent':106,'hsilu3_candidate_exponent':126,'hsilu_batch_exponent':math.log2((12*3**64-10*2**64)/7)}
    (Path(__file__).resolve().parent/'check_results.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))
if __name__=='__main__':main()
