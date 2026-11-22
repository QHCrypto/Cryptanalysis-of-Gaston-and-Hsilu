"""Extract the 13 agreed representative endpoints and their archived cuts."""
import json
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
def read(path): return json.loads((ROOT / path).read_text())
def words(ws): return [f'{int(x, 16):016x}' for x in ws]

catalog = read('研究汇总_20261019/results_catalog.json')
endpoints = {e['code']: e for e in catalog['endpoints']}
jobs = []
splits = {}
def middle(code, suffix, algo, rounds, start, delta, mask):
    job = dict(id=code+suffix, algorithm=algo, rounds=rounds,
               start_round_index=start, delta_core=words(delta), lambda_core=words(mask))
    jobs.append(job)
    return job['id']

paper = read('gaston_dl/results/fixed_paper_tail/search.json')
tail = next(r for r in paper['reports'] if r['id']=='D3_chi2_to_3')
c = next(c for c in tail['top_candidates'] if c['input_column']==1 and c['input_pattern']==16)
branches = []
for i,b in enumerate(c['prefix']['branches']):
    job = middle('G4c',f'_b{i}', 'Gaston',1,9,
                 b['middle_chi_input_difference']['words'], tail['join_mask']['words'])
    branches.append(dict(job_id=job, probability=f"{b['ddt_count']}/{b['denominator']}",
                         differential_output=words(b['first_chi_output_difference']['words'])))
splits['G4c'] = dict(split=[1,1,2], linear_squared_weight=38, branches=branches,
    source='gaston_dl/results/fixed_paper_tail/search.json', tail_id=tail['id'],
    aggregation='sum of two distinct retained DDT events; each probability 1/4; not one trail of probability 1/2')

pool = read('gaston_dl/results/datf_pools/paper_full24/search_results.json')
c = next(r for r in pool['rows'] if r['id']=='D1_M1_D3_chi2_to_3')['best']
assert words(c['input_difference']) == endpoints['G4d']['delta_core']
assert words(c['output_mask']) == endpoints['G4d']['lambda_core']
job = middle('G4d','_b0','Gaston',1,9,c['middle_input_difference'],c['middle_output_mask'])
splits['G4d'] = dict(split=[1,1,2],linear_squared_weight=38,
    branches=[dict(job_id=job,probability='1/4',differential_output=c['differential_output'])],
    source='gaston_dl/results/datf_pools/paper_full24/search_results.json',
    tail_id=c['linear_id'], aggregation='single retained DDT event')

h = read('hsilu_dl/results/research_summary.json')['five_candidates'][0]
assert h['id']=='R5_pool_000000'
assert h['input_difference']['words']==endpoints['H5']['delta_core']
assert h['core_output_mask']['words']==endpoints['H5']['lambda_core']
job = middle('H5','_prefix3','Hsilu',3,7,h['input_difference']['words'],h['head']['join_mask']['words'])
splits['H5'] = dict(split=[0,3,2],linear_squared_weight=2*h['linear_weight'],
    branches=[dict(job_id=job,probability='1')],
    source='hsilu_dl/results/research_summary.json',tail_id=h['head']['id'],
    aggregation='no differential outer prefix; original search used a three-chi prefix and a fixed two-chi linear tail')

(OUT/'middle_manifest.json').write_text(json.dumps({'jobs':jobs},ensure_ascii=False,indent=2)+'\n')
(OUT/'cuts.json').write_text(json.dumps(splits,ensure_ascii=False,indent=2)+'\n')
print('Prepared',len(endpoints),'representative endpoints;',len(jobs),'middle jobs.')
