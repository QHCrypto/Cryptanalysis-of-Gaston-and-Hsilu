import json,sys,subprocess,concurrent.futures,math
from fractions import Fraction
from search_new import H,ROOT,save,rec
from hsilu import words,encode
from pool_search import decode,LAT
WORKER=ROOT/'DATF算法3复核_20261021/worker.py'
def main():
 candidates=json.loads((H/'search_113.json').read_text())['top'];rows=[];jobs=[]
 for i,c in enumerate(candidates):
  d=tuple(int(x,16) for x in c['input_difference']);col=next(j for j in range(64) if any(w>>j&1 for w in d));a=sum((w>>col&1)<<k for k,w in enumerate(d))
  res=subprocess.run([str(ROOT/'hsilu_dl/build/hsilu_dl'),'prefix-eval',str(col),str(a),';'.join(c['middle_mask']),'78'],capture_output=True,text=True,check=True)
  exact=Fraction(int(res.stdout),32);q2=Fraction(int(c['q2']['numerator']),int(c['q2']['denominator']))
  id=f'new5_{i:03d}';c=dict(c,id=id,exact_prefix2=rec(exact),prefix_splice=rec(exact*q2));rows.append(c)
  jobs.append(dict(id=id+'_prefix2',algorithm='Hsilu',rounds=2,start_round_index=7,delta_core=c['input_difference'],lambda_core=c['middle_mask'],scope='middle of (0,2,3)'))
 save(H/'new_candidates.json',dict(candidates=rows,scope='128 best distinct endpoints from affine matching; one selected three-round tail per endpoint'))
 # Full Algorithm3 feasibility attempts for three strongest NONZERO full-prefix scores,
 # as well as exact single-chi middle verification for their original 1+1+3 cuts.
 best=sorted(rows,key=lambda c:-abs(Fraction(int(c['prefix_splice']['numerator']),int(c['prefix_splice']['denominator']))))[:3]
 for c in best:
  jobs.append(dict(id=c['id']+'_whole5',algorithm='Hsilu',rounds=5,start_round_index=7,delta_core=c['input_difference'],lambda_core=c['output_mask'],scope='whole feasibility'))
  jobs.append(dict(id=c['id']+'_middle1',algorithm='Hsilu',rounds=1,start_round_index=8,delta_core=c['middle_input'],lambda_core=c['middle_mask'],scope='middle of (1,1,3)'))
 save(H/'new_manifest.json',dict(jobs=jobs,top_ids=[c['id'] for c in best]))
 def run(pair):
  i,j=pair;out=H/'raw'/f'{j["id"]}.json'
  if out.exists():return j['id']+' cached'
  r=subprocess.run([sys.executable,str(WORKER),'--manifest',str(H/'new_manifest.json'),'--index',str(i),'--output',str(out),'--seconds','30'],text=True,capture_output=True)
  (H/'raw'/f'{j["id"]}.log').write_text(r.stdout+r.stderr)
  if r.returncode:raise RuntimeError(r.stderr)
  return r.stdout.strip()
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as e:
  for s in e.map(run,enumerate(jobs)):print(s,flush=True)
 print('new candidate exact prefix2 nonzero',sum(int(c['exact_prefix2']['numerator'])!=0 for c in rows),flush=True)
 print('top',[(c['id'],encode(tuple(int(x,16) for x in c['input_difference']))['bits'],c['prefix_splice']) for c in best],flush=True)
if __name__=='__main__':main()
