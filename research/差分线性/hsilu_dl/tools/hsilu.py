"""Hsilu Table 9; standalone integer and symbolic linear maps."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'gaston_dl/tools'))
from datf import OrdinaryDATF, ONE, ZERO, add
from datf_k06 import symbolic_chi_pair

MASK=(1<<64)-1
N1=(0,3,1,0,5,5,7,11,12,12,13,14,15)
N2=(4,4,2,1,2,3,6,8,9,6,10,11,9)
SH=(47,51,54,30,4,55,0,14,59,4,10,45,53)
ROT=(0,57,13,31,20)
RC=(0xf0,0xe1,0xd2,0xc3,0xb4,0xa5,0x96,0x87,0x78,0x69,0x5a,0x4b)
def rol(x,n):
    n%=64
    return ((x<<n)|(x>>((-n)%64)))&MASK
def linear(a):
    x=list(a)
    for u,v,n in zip(N1,N2,SH):x.append(x[u]^rol(x[v],n))
    return tuple(rol(x[13+i],ROT[i]) for i in range(5))
def transpose(a):
    x=[0]*18
    for i in range(5):x[13+i]=rol(a[i],-ROT[i])
    for j in range(12,-1,-1):
        x[N1[j]]^=x[5+j];x[N2[j]]^=rol(x[5+j],-SH[j])
    return tuple(x[:5])
def chi(a):return tuple((a[i]^((~a[(i+1)%5])&a[(i+2)%5]))&MASK for i in range(5))
def core(a,R,start=None):
    if start is None:start=12-R
    a=chi(a)
    for rc in RC[start+1:start+R]:
        a=linear(a);a=chi((a[0]^rc,*a[1:]))
    return a
def pack(a):return sum(w<<(64*i) for i,w in enumerate(a))
def unpack(a):return tuple((a>>(64*i))&MASK for i in range(5))
def embed(j,v):return tuple(((v>>i)&1)<<j for i in range(5))
def words(a):return [f'{w:016x}' for w in a]
def active(a):return (a[0]|a[1]|a[2]|a[3]|a[4]).bit_count()
def encode(a):return {'words':words(a),'bits':[[i,j] for i in range(5) for j in range(64) if a[i]>>j&1],'active_columns':active(a),'bit_weight':sum(w.bit_count() for w in a)}
def inverse_function(f):
    cols=[pack(f(unpack(1<<j))) for j in range(320)]
    rows=[sum(((cols[j]>>i)&1)<<j for j in range(320))|(1<<(320+i)) for i in range(320)]
    for j in range(320):
        k=next(k for k in range(j,320) if rows[k]>>j&1)
        rows[j],rows[k]=rows[k],rows[j]
        for k in range(320):
            if k!=j and rows[k]>>j&1:rows[k]^=rows[j]
    mask=(1<<320)-1
    assert all((r&mask)==1<<j for j,r in enumerate(rows))
    invrows=[r>>320 for r in rows]
    def inverse(a):
        a=pack(a)
        return unpack(sum(((a&r).bit_count()%2)<<j for j,r in enumerate(invrows)))
    return inverse
def symbolic_linear(a):
    x=list(a)
    for u,v,n in zip(N1,N2,SH):x.append([add(x[u][j],x[v][(j-n)%64]) for j in range(64)])
    return [[x[13+i][(j-ROT[i])%64] for j in range(64)] for i in range(5)]
def construct(delta,R,max_nodes=100000,start=None):
    if start is None:start=12-R
    assert 0<=start and start+R<=12
    e=OrdinaryDATF(max_nodes)
    a=[[e.fresh(kind='input',row=i,column=j) for j in range(64)] for i in range(5)]
    d=[[ONE if delta[i]>>j&1 else ZERO for j in range(64)] for i in range(5)]
    for r in range(R):
        if r:
            a,d=symbolic_linear(a),symbolic_linear(d)
            for j in range(8):
                if RC[start+r]>>j&1:a[0][j]=add(a[0][j],ONE)
            for i in range(5):
                for j in range(64):
                    a[i][j]=e.transition(a[i][j],'alpha',r,i,j)
                    d[i][j]=e.transition(d[i][j],'beta',r,i,j)
        a,d=symbolic_chi_pair(a,d)
    return e,a,d
