from pathlib import Path
import sys, json, math
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent))
from independent_checks import RC, lg, lh, chi

def rot(x,n):
    n%=64
    if not n:return x
    return (x<<np.uint64(n))|(x>>np.uint64(64-n))
def L(x,kind):
    if kind=='g':
        a=[rot(v,r) for v,r in zip(x,[0,60,22,27,4])]
        p=a[0]^a[1]^a[2]^a[3]^a[4]
        q=rot(a[0],25)^rot(a[1],32)^rot(a[2],52)^rot(a[3],60)^rot(a[4],63)
        e=rot(p^rot(p,1)^q^rot(q,18),23)
        return [rot(v^e,r) for v,r in zip(a,[0,56,31,46,43])]
    z=list(x)
    for a,b,r in [(0,4,47),(3,4,51),(1,2,54),(0,1,30),(5,2,4),(5,3,55),(7,6,0),(11,8,14),(12,9,59),(12,6,4),(13,10,10),(14,11,45),(15,9,53)]:z.append(z[a]^rot(z[b],r))
    return [rot(z[13+i],r) for i,r in enumerate([0,57,13,31,20])]
def sandwich(x,kind):
    for i,c in enumerate(RC[-4:]):
        x=list(x);x[0]=x[0]^np.uint64(c)
        x=[x[j]^((~x[(j+1)%5])&x[(j+2)%5]) for j in range(5)]
        if i<3:x=L(x,kind)
    return x
def parity(x):
    for n in [32,16,8,4,2,1]:x=x^(x>>np.uint64(n))
    return x&np.uint64(1)
rg=np.random.default_rng(20261117)
check=rg.integers(0,2**64,size=(5,3),dtype=np.uint64)
for kind,scalar in [('g',lg),('h',lh)]:
    got=L(check,kind)
    for j in range(3):assert [int(v[j]) for v in got]==scalar([int(v[j]) for v in check])
out={}
for kind,delta,mask in [('g',[2,0,0,0,0],[0,1,0,0x0040000000000000,0]),('h',[0,0,0,1<<18,0],[1<<9,0,4,0,1<<9])]:
    n=2**22;ones=0
    for _ in range(n//65536):
        x=rg.integers(0,2**64,size=(5,65536),dtype=np.uint64)
        y=[a^np.uint64(d) for a,d in zip(x,delta)]
        a=sandwich(x,kind);b=sandwich(y,kind)
        p=np.zeros(65536,dtype=np.uint64)
        for u,v,m in zip(a,b,mask):p^=(u^v)&np.uint64(m)
        ones+=int(np.sum(parity(p)))
    cor=1-2*ones/n;se=math.sqrt((1-cor*cor)/n)
    out[kind]={'N':n,'odd_count':ones,'correlation':cor,'normal95':[cor-1.96*se,cor+1.96*se]}
print(json.dumps(out,indent=2))
open(Path(__file__).resolve().parent/'dl_results.json','w').write(json.dumps(out,indent=2)+'\n')
