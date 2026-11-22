"""Sample GF(2) rank of the 17 affine target equations, without dependencies.

The five c words are fixed first-round state contributions. This is a
structural diagnostic, not a full Gaston-XOF attack or a complexity proof.
"""
import random, collections
MASK=(1<<64)-1
U={0,4,5,10,11,22,23,27,33,34,40,41,45,51,56,57,63}
S=sorted(set(range(64))-U)
J=[0,1,6,17,18,19,24,29,30,41,42,47,52,53,59,60,63]
D=[[0,23,24,48,2],[19,20,51,5],[45,46,33,51],[50,51,46,0],[27,28,26,44]]
def rot(x,d): return ((x>>d)|(x<<(64-d)))&MASK if d else x
def u_of(m): return rot(m,23)^rot(m,24)^rot(m,48)^rot(m,2)
def f(m,c):
 u=u_of(m)
 y=[m^u^c[0],rot(u,56)^c[1],rot(u,31)^c[2],rot(u,46)^c[3],rot(u,43)^c[4]]
 x=[y[i]^((~y[(i+1)%5]&MASK)&y[(i+2)%5]) for i in range(5)]
 o=0
 for i in range(5):
  for d in D[i]:o^=rot(x[i],d)
 return sum(((o>>j)&1)<<r for r,j in enumerate(J))
rows=[sum(1<<t for t in ((s+23)%64,(s+24)%64,(s+48)%64,(s+2)%64)) for s in S]
piv=[]
for col in range(64):
 k=next((k for k in range(len(piv),len(rows)) if rows[k]>>col&1),None)
 if k is None:continue
 rows[len(piv)],rows[k]=rows[k],rows[len(piv)]
 for k2 in range(len(rows)):
  if k2!=len(piv) and rows[k2]>>col&1:rows[k2]^=rows[len(piv)]
 piv.append(col)
 if len(piv)==47:break
assert len(piv)==47
null=[]
for free in set(range(64))-set(piv):
 v=1<<free
 for k,p in enumerate(piv):
  if rows[k]>>free&1:v|=1<<p
 assert all((u_of(v)>>s&1)==0 for s in S)
 null.append(v)
def rank(cols):
 basis={}
 for v in cols:
  while v:
   p=v.bit_length()-1
   if p not in basis:
    basis[p]=v;break
   v^=basis[p]
 return len(basis)
rng=random.Random(1234567)
for ci in range(4):
 c=[0]*5 if ci==0 else [rng.getrandbits(64) for _ in range(5)]
 counts=collections.Counter()
 for _ in range(3000):
  m=rng.getrandbits(64)
  base=f(m,c)
  counts[rank([f(m^n,c)^base for n in null])]+=1
 print(ci,sorted(counts.items()))
