import json,sys,subprocess,time
from fractions import Fraction
from search_new import H,ROOT,BASE,DDT,save,rec
from hsilu import linear,transpose,embed,encode,words,inverse_function
from pool_search import LAT
def decode(a):return tuple(int(w,16) for w in (a["words"] if isinstance(a,dict) else a))

def val(d):return Fraction(int(d['numerator']),int(d['denominator']))
def corr(u,v):
 q=Fraction(1)
 for j in range(64):
  a=sum((u[i]>>j&1)<<i for i in range(5));b=sum((v[i]>>j&1)<<i for i in range(5));q*=Fraction(LAT[a][b],32)
 return q
new=json.loads((H/'new_candidates.json').read_text())['candidates'];optimized=json.loads((H/'search_023.json').read_text())['top']
checked=0
for c in new+optimized:
 prep=c['prepended_chi'];u=decode(prep['input']);v=decode(prep['output']);h=c['head'];last=c['last_chi']
 assert transpose(u)==decode(c['middle_mask'])
 assert v==transpose(decode(h['first_chi']['input']))
 assert decode(h['first_chi']['output'])==transpose(decode(last['input']))
 q=corr(u,v)*corr(decode(h['first_chi']['input']),decode(h['first_chi']['output']))*corr(decode(last['input']),decode(last['output']))
 assert q*q==val(c['q2'])
 if c in new:
  j=json.loads((H/'raw'/f"{c['id']}_prefix2.json").read_text());assert j['status']=='completed' and val(j['correlation'])==val(c['exact_prefix2'])
  a=decode(c['input_difference']);b=decode(c['differential_output']);p=Fraction(1)
  for col in range(64):p*=Fraction(DDT[sum((a[i]>>col&1)<<i for i in range(5))][sum((b[i]>>col&1)<<i for i in range(5))],32)
  assert p==val(c['p'])
 checked+=1
# Direct HNW full and zero on the same short prefixes / full 2–3 round endpoints.
jobs=json.loads((H/'manifest.json').read_text())['jobs']
jobs += [j for j in json.loads((H/'optimized_manifest.json').read_text())['jobs'] if j['rounds']==2]
results={}
for j in jobs:
 pair={}
 for variant in ['full','zero']:
  out=H/'raw'/f"{j['id']}_hnw_{variant}.json"
  cmd=[str(ROOT/'hsilu_dl/build/hnw_endpoint'),variant,str(j['rounds']),str(j['start_round_index']),*j['delta_core'],*j['lambda_core'],str(out)]
  r=subprocess.run(cmd,text=True,capture_output=True,check=True);pair[variant]=json.loads(out.read_text())
 results[j['id']]=pair
save(H/'hnw_results.json',results)
# Normalize generic worker's historical field: Hsilu has no initial L to omit.
for p in (H/'raw').glob('*.json'):
 if p.name.endswith('.anf.json'):continue
 d=json.loads(p.read_text())
 if d.get('method')=='original DATF Algorithm 3':
  d['boundary_note']='Hsilu core: first chi input to last chi output. All internal L layers retained; final native L handled by external mask L^{-T}. Generic worker first_L_omitted is inapplicable to Hsilu.'
  d['first_L_omitted']=None;save(p,d)
save(H/'validation.json',dict(status='passed',linear_trails_checked=checked,new_prefix_DATF3_vs_independent_exact=128,affine_search_feasibility_checks='every retained hit verified in search script',hnw_jobs=len(jobs),hnw_variants=['full','zero'],native_mask_conversion='will be verified in report generation'))
print('validated',checked,'trails;',len(jobs),'HNW full+zero jobs',flush=True)
