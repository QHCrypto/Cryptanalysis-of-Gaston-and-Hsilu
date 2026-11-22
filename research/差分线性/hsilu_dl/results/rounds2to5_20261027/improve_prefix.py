"""Seek simultaneous same-sign first-DDT contributions in best-weight L3 support."""
import json,itertools,time,subprocess
from fractions import Fraction
from search_113 import setup,xor
from search_new import H,BASE,ROOT,DDT,save,rec
from hsilu import embed,linear,transpose,chi,words,encode
from pool_search import decode
from middle2 import two

def equations(d,forms,sign):
 eq=[]
 for j in range(64):
  for k in range(5):
   z=(forms[((k-1)%5)*64+j] if d[(k+1)%5]>>j&1 else 0)^(forms[((k-2)%5)*64+j] if d[(k-1)%5]>>j&1 else 0)
   if z:eq.append(z)
 z=1 if sign==-1 else 0
 for i,w in enumerate(chi(d)):
  while w:
   b=w&-w;j=b.bit_length()-1;z^=forms[i*64+j];w^=b
 eq.append(z)
 return eq

def solve(eq):
 piv={}
 for z in eq:
  while z&~1:
   k=z.bit_length()-1
   if k in piv:z^=piv[k]
   else:piv[k]=z;break
  else:
   if z:return None
 s=1
 for k in sorted(piv):
  if ((piv[k]^(1<<k))&s).bit_count()%2:s|=1<<k
 return s

def main():
 t=time.perf_counter();pool=json.loads((BASE/'linear_pool.json').read_text());h=next(h for h in pool['heads'] if h['id']=='H0064');B,u0,ubs,forms,w=setup(h)
 last=min(pool['bridges'][h['bridge_id']]['last_chi_choices'],key=lambda z:(z['weight'],z['output']['words']))
 q2=Fraction(1,2**(2*(w+h['first_chi']['weight']+last['weight'])));results=[];tests=0
 for col in range(64):
  for lane in range(5):
   a=1<<lane;bs=[b for b in range(32) if DDT[a][b]];best=None
   eqs={(b,sign):equations(linear(embed(col,b)),forms,sign) for b in bs for sign in [1,-1]}
   for k in range(len(bs),0,-1):
    for subset in itertools.combinations(bs,k):
     for sign in [1,-1]:
      tests+=1;s=solve([e for b in subset for e in eqs[b,sign]])
      if s is None:continue
      u=u0
      for i,ub in enumerate(ubs,1):
       if s>>i&1:u=xor(u,ub)
      g=transpose(u);model=two(embed(col,a),g)
      if best is None or abs(model)>abs(best['rho']):best=dict(rho=model,u=u,g=g,subset=subset,sign=sign)
   if best is not None:
    native=subprocess.run([str(ROOT/'hsilu_dl/build/hsilu_dl'),'prefix-eval',str(col),str(a),';'.join(words(best['g'])),'78'],capture_output=True,text=True,check=True)
    exact=Fraction(int(native.stdout),32)
    results.append(dict(id=f'optimized_c{col}_a{a}',split=[0,2,3],input_difference=words(embed(col,a)),output_mask=last['output']['words'],middle_mask=words(best['g']),
     exact_prefix2=rec(exact),hnw_prefix2=rec(best['rho']),q2=rec(q2),score=rec(exact*q2),prepended_chi=dict(input=words(best['u']),output=words(B),weight=w),head=h,last_chi=last,matched_branches=list(best['subset']),matched_sign=best['sign']))
  if col%16==15:print('same-sign',col+1,'columns','found',len(results),flush=True)
 results.sort(key=lambda c:-abs(Fraction(int(c['score']['numerator']),int(c['score']['denominator']))))
 save(H/'search_023.json',dict(status='completed',scope='H0064 backwards linear extension; all 320 one-bit input differences; simultaneous same-sign feasibility for every nonempty subset of four first-DDT outputs; one free assignment per feasible system',tests=tests,found=len(results),top=results,elapsed_seconds=time.perf_counter()-t))
 print('best',[(r['id'],r['score']) for r in results[:5]],flush=True)
if __name__=='__main__':main()
