"""Collect all explicitly reported endpoints, without re-running the search space."""
import csv
import gzip
import hashlib
import json
from pathlib import Path

OUT=Path(__file__).resolve().parent
ROOT=OUT.parent
G=ROOT/'gaston_dl/results'
H=ROOT/'hsilu_dl/results'
jobs=[];by_key={};sources={}
def read(p):
    sources[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    return json.loads(p.read_text())
def words(s):
    if isinstance(s,dict):s=s['words']
    if isinstance(s,str):s=s.split(';')
    return [f'{int(w,16):016x}' for w in s]
def add(identifier,algo,R,start,delta,mask,previous,source,**extra):
    d,m=words(delta),words(mask)
    key=(algo,R,start,tuple(d),tuple(m))
    if key in by_key:
        row=by_key[key];row['aliases'].append(identifier)
        row['sources'].append(source)
        return
    row=dict(id=identifier,algorithm=algo,rounds=R,start_round_index=start,
             delta_core=d,lambda_core=m,previous_datf=previous,
             aliases=[identifier],sources=[source],**extra)
    by_key[key]=row;jobs.append(row)

p=ROOT/'研究汇总_20261019/results_catalog.json'
for x in read(p)['endpoints']:
    add(x['code'],x['algorithm'],x['rounds'],12-x['rounds'],x['delta_core'],x['lambda_core'],x['datf'],str(p),sampling=x['sampling'],previous_hnw=x['hnw'])
p=G/'datf_k06/results.json'
for x in read(p)['results']:
    add(x['id'],'Gaston',x['rounds'],x['start_round_index'],x['input_difference'],x['output_mask'],x['correlation'],str(p))
p=G/'hnw26_comparison/results.json'
for x in read(p)['results']:
    c=x['candidate'];add(x['id'],'Gaston',c['R'],c['start_round_index'],c['input_difference'],c['output_mask'],x['datf_correlation'],str(p))
p=H/'research_summary.json'
d=read(p)
for x in d['four_DATF']+d['five_DATF']:
    add(x['id'],'Hsilu',x['rounds'],12-x['rounds'],x['input_difference'],x['core_output_mask'],x['correlation'],str(p))
p=G/'extend_121_tail/results.json';d=read(p)
for x in d['results']:
    add('G4a_'+x['id'],'Gaston',d['rounds'],d['start_round_index'],d['source_candidate']['input_difference'],x['extension']['final_output'],x['DATF']['correlation'],str(p))
p=G/'datf_one_one_tail/summary.csv'
sources[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
for x in csv.DictReader(p.open()):
    identifier='early_'+x['tail_id']+'_rank'+x['candidate_rank']
    add(identifier,'Gaston',int(x['R']),int(x['start_round_index']),x['input_words'],x['output_words'],x['full_DATF_rho'],str(p))
    add(identifier+'_prefix2','Gaston',2,int(x['start_round_index']),x['input_words'],x['join_words'],x['prefix_DATF_rho'],str(p))
p=G/'fixed_tail_extension/results.json';d=read(p)
for x in d['rows']:
    c=x['best_nonzero_candidate']
    if c:add('forward_'+x['id'],'Gaston',x['rounds'],x['start_round_index'],c['input_difference'],x['output_mask'],c['full_DATF']['correlation'],str(p))
p=G/'fixed_tail_extension/hnw26/results.json'
for x in read(p)['rows']:
    for folder,identifier in [('full_jobs',x['endpoint_id']),('prefix_jobs',x['exact_prefix']['id'])]:
        q=G/'fixed_tail_extension'/folder/(identifier+'.json.gz')
        sources[str(q)]=hashlib.sha256(q.read_bytes()).hexdigest()
        d=json.loads(gzip.open(q,'rt').read())['result']
        add('forward_'+x['row_id']+'_'+folder,'Gaston',d['rounds'],d['start_round_index'],d['input_difference'],d['output_mask'],d['correlation'],str(q))

# The two-to-two pool's 176 retained full endpoints are a separate, reproducible
# batch; the 24,515 historical screening evaluations are not silently conflated
# with the set of representative results above.
p=G/'datf_pools/linear_two_match/endpoint_check/results.json'
for x in read(p)['entries']:
    c=x['candidate'];add('L2pool_'+x['endpoint_id'][:12],'Gaston',c['R'],c['start_round_index'],c['input_difference'],c['output_mask'],x['result']['correlation'],str(p))

p=H/'three_round_datf/results.json'
for x in read(p)['results']:
    add(x['id'],'Hsilu',x['rounds'],12-x['rounds'],x['input_difference'],x['core_output_mask'],x['correlation'],str(p))
p=G/'extend_121_tail/results.json';d=read(p)
add('G4a_shifted_constant_window','Gaston',4,d['start_round_index'],d['source_candidate']['input_difference'],d['source_candidate']['output_mask'],d['shifted_old_endpoint_DATF']['correlation'],str(p))

(OUT/'manifest.json').write_text(json.dumps({'scope':'reported representative endpoints, all 32 Hsilu retained endpoints, 26 HNW comparison endpoints, K06, early and forward/backward extensions, 176 retained L2 endpoints; deduplicated by full endpoint AND constant window','excludes':'complete historical search space and 24515 screening evaluations, not re-run','jobs':jobs,'source_sha256':sources},ensure_ascii=False,indent=2)+'\n')
print('unique jobs:',len(jobs),'source records:',sum(len(j['aliases']) for j in jobs))

if __name__=='__main__':pass
