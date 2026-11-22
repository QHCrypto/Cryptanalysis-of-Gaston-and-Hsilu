from pathlib import Path
import sys,random,json,statistics,collections,math
sys.path.insert(0,str(Path(__file__).resolve().parent))
from independent_checks import prefix,stages,solve,comb,embed

rng=random.Random(20261117)
bg=[j for j in range(64) if j not in [12,18]]
cond=[(0,[21,42,46,48,52],[19,25]),(1,[7,13,26],[20,35,39,45]),(2,[],[57,63]),(3,[26,35,41,45],[21,27,33,42,46,52]),(4,[57,63],[7,13,26,32])]
c1=[(i,j,b) for i,zero,one in cond for b,js in [(0,zero),(1,one)] for j in js]
c2=[(0,53,0),(0,48,1),(1,15,1),(1,52,1)]
guesses=[(i,j) for i,js in enumerate([[4,11,18,22,31,53],[2,6,8,58,61],[15,34,48,53],[4,9,22,49,59],[5,54,55,58,59]]) for j in js]
counts=collections.Counter();etas=[]
for _ in range(4096):
    counts['prefixes']+=1
    s=prefix(rng.getrandbits(64),rng.getrandbits(64),3)
    if not all((s[1]>>j)&1==0 and (s[4]>>j)&1==1 for j in [12,18]):continue
    counts['capacity_pass']+=1
    def hs(x):return stages([x]+s[1:],3)
    z=hs(0)[0];cd=[[(hs(1<<j)[0][i]^z[i]) for i in range(5)] for j in bg]
    rows=[sum(((cd[k][i]>>j)&1)<<k for k in range(62)) for i,j,b in c1]
    rhs=[b^((z[i]>>j)&1) for i,j,b in c1]
    b0,B,rk=solve(rows,rhs,62)
    if b0 is None or rk!=32:continue
    counts['consistent_rank32']+=1
    z=hs(embed(b0,bg))[0]
    gr=[sum((((hs(embed(b0^v,bg))[0][i]^z[i])>>j)&1)<<k for k,v in enumerate(B)) for i,j in guesses]
    _,BB,grank=solve(gr,[0]*25,30)
    if grank!=23:continue
    counts['retained']+=1
    BBx=[embed(comb(v,B),bg) for v in BB]
    backgrounds=0
    for trial in range(256):
        tv=rng.getrandbits(30);bv=embed(b0^comb(tv,B),bg);zs=hs(bv)[1]
        rows2=[sum((((hs(bv^v)[1][i]^zs[i])>>j)&1)<<k for k,v in enumerate(BBx)) for i,j,b in c2]
        rhs2=[b^((zs[i]>>j)&1) for i,j,b in c2]
        a,bs,r=solve(rows2,rhs2,7)
        if a is not None:backgrounds+=1<<len(bs)
    etas.append(backgrounds/256)
out={'counts':dict(counts),'sampled_branches_per_retained_prefix':256,'eta_mean':statistics.mean(etas),'eta_min':min(etas),'eta_max':max(etas),'eta_histogram':dict(sorted(collections.Counter(etas).items())),'scope':'Independent sampling of legal prefixes and branches; not an exhaustive yield guarantee or preimage experiment.'}
print(json.dumps(out,indent=2))
open(Path(__file__).resolve().parent/'h3_prefix_yield.json','w').write(json.dumps(out,indent=2)+'\n')
