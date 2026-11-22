import random,json,subprocess
from fractions import Fraction
from search_new import H,ROOT,save,rec
from middle2 import two
from hsilu import words,embed,linear,transpose,core,chi
from pool_search import chi5
# Independent DDT enumeration across initial active columns; local DLCT lookup.
DDT=[[sum(chi5(x)^chi5(x^a)==b for x in range(32)) for b in range(32)] for a in range(32)]
DLCT=[[sum((-1)**((v&(chi5(x)^chi5(x^a))).bit_count()%2) for x in range(32)) for v in range(32)] for a in range(32)]
def brute(delta,gamma):
 active=[(j,sum((delta[i]>>j&1)<<i for i in range(5))) for j in range(64) if any(w>>j&1 for w in delta)]
 states=[((0,)*5,Fraction(1))]
 for j,a in active:
  states=[(tuple(x^y for x,y in zip(st,embed(j,b))),p*Fraction(n,32)) for st,p in states for b,n in enumerate(DDT[a]) if n]
 out=Fraction(0)
 for st,p in states:
  d=linear(st)
  for j in range(64):
   a=sum((d[i]>>j&1)<<i for i in range(5));v=sum((gamma[i]>>j&1)<<i for i in range(5));p*=Fraction(DLCT[a][v],32)
   if not p:break
  out+=p
 return out
rng=random.Random(20261027);cases=[]
for k in range(120):
 d=embed(rng.randrange(64),rng.randrange(1,32));g=embed(rng.randrange(64),rng.randrange(1,32))
 if k%3==0:d=tuple(x^y for x,y in zip(d,embed(rng.randrange(64),rng.randrange(1,32))))
 if k%4==0:g=transpose(g)
 a=two(d,g);b=brute(d,g);assert a==b,(k,a,b)
 cases.append(dict(delta=words(d),mask=words(g),rho=rec(a)))
# Check exact mapping from input x to affine derivative independently.
from middle2 import affine
for k in range(20):
 d=tuple(rng.getrandbits(64) for _ in range(5));x=tuple(rng.getrandbits(64) for _ in range(5));assignment=1|sum(w<<(64*i+1) for i,w in enumerate(x))
 actual=linear(tuple(a^b for a,b in zip(chi(x),chi(tuple(a^b for a,b in zip(x,d))))))
 forms=affine(d);got=tuple(sum(((forms[64*i+j]&assignment).bit_count()%2)<<j for j in range(64)) for i in range(5));assert got==actual
# Native HNW zero recurrence, two nonzero and two zero controls.
controls=[next(c for c in cases if int(c['rho']['numerator'])==0),next(c for c in cases if int(c['rho']['numerator'])!=0)]
for i,c in enumerate(controls):
 p=H/'raw'/f'backend_control_{i}.json';r=subprocess.run([str(ROOT/'hsilu_dl/build/hnw_endpoint'),'zero','2','8',*c['delta'],*c['mask'],str(p)],capture_output=True,text=True,check=True)
 v=json.loads(p.read_text());q=v['correlation'];assert Fraction(int(q['numerator']),int(q['denominator']))==Fraction(int(c['rho']['numerator']),int(c['rho']['denominator']))
save(H/'middle2_validation.json',dict(status='passed',affine_derivative_checks=20,independent_DDT_DLCT_sums=120,native_hnw_controls=2,nonzero_controls=sum(int(c['rho']['numerator'])!=0 for c in cases),cases=cases))
print('middle2 validated',len(cases),'cases',flush=True)
