import json,sys,subprocess,time,concurrent.futures
from pathlib import Path
from search_new import H,BASE,ROOT,save
from hsilu import words,embed
WORKER=ROOT/'DATF算法3复核_20261021/worker.py'
def main():
 jobs=[];mapping={}
 def add(id,R,start,d,m,scope):
  key=(R,start,tuple(d),tuple(m))
  if key in mapping:return mapping[key]
  mapping[key]=id;jobs.append(dict(id=id,algorithm='Hsilu',rounds=R,start_round_index=start,delta_core=d,lambda_core=m,scope=scope));return id
 for ln in (BASE/'initial_jobs.txt').read_text().splitlines():
  id,r,d,m=ln.split();add(id,int(r),12-int(r),d.split(';'),m.split(';'),'whole')
 reps=json.loads((ROOT/'hsilu_dl/results/research_summary.json').read_text())
 links={}
 for c in reps['four_candidates']:
  links[c['id']]=add(c['id']+'_prefix2',2,8,c['input_difference']['words'],c['head']['join_mask']['words'],'middle of (0,2,2)')
 for c in reps['five_candidates']:
  links[c['id']]=add(c['id']+'_prefix3',3,7,c['input_difference']['words'],c['head']['join_mask']['words'],'middle of (0,3,2)')
 manifest={'jobs':jobs,'candidate_to_job':links};save(H/'manifest.json',manifest)
 def run(pair):
  i,j=pair;out=H/'raw'/f'{j["id"]}.json'
  if out.exists():return j['id']+' cached'
  cmd=[sys.executable,str(WORKER),'--manifest',str(H/'manifest.json'),'--index',str(i),'--output',str(out),'--seconds','30']
  r=subprocess.run(cmd,text=True,capture_output=True)
  (H/'raw'/f'{j["id"]}.log').write_text(r.stdout+r.stderr)
  if r.returncode:raise RuntimeError(j['id']+r.stderr)
  return r.stdout.strip()
 (H/'raw').mkdir(exist_ok=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=3) as e:
  for s in e.map(run,enumerate(jobs)):print(s,flush=True)
if __name__=='__main__':main()
