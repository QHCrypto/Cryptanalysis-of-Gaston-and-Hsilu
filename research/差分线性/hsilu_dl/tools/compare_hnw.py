#!/usr/bin/env python3
import json
import subprocess
from pathlib import Path
from hsilu import embed,encode

root=Path(__file__).resolve().parents[1]
base=root/'results/baseline';out=root/'results/hnw';out.mkdir(exist_ok=True)
four=json.loads((base/'pool_selection.json').read_text())
five=json.loads((base/'five_selection.json').read_text())
jobs=[{'id':'R3_structural','rounds':3,'input_difference':encode(embed(0,8)),
       'core_output_mask':encode(embed(29,1))},
      *[r for r in four if r['id'] in ('R4_pool_000000','R4_pool_000027')],five[0]]
summary=[]
for job in jobs:
    for mode in ('zero','full'):
        path=out/(job['id']+'_'+mode+'.json')
        if path.exists():raise FileExistsError(path)
        args=[root/'build/hnw_endpoint',mode,str(job['rounds']),str(12-job['rounds']),
              *job['input_difference']['words'],*job['core_output_mask']['words'],path]
        subprocess.run(args,check=True)
        data=json.loads(path.read_text())
        summary.append({'id':job['id'],'rounds':job['rounds'],'mode':mode,
                        'correlation':data['correlation'],'completed_branches':data['completed_branches'],
                        'path':str(path.relative_to(root))})
        (out/'summary.json').write_text(json.dumps({'method':'HNW26 round-wise independence approximation, exact dyadic arithmetic',
          'not_exact_permutation_correlation':True,'results':summary},indent=2)+'\n')
