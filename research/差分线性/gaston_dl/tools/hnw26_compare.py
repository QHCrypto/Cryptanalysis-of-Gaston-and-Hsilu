#!/usr/bin/env python3
"""Re-evaluate fixed saved Gaston endpoints with HNW26 rd=1, rl=0."""
import argparse
import concurrent.futures
import hashlib
import json
import subprocess
import time
from fractions import Fraction
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


def save(path,data):
    temp=path.with_name(path.name+'.partial')
    temp.write_text(json.dumps(data,ensure_ascii=False,indent=2)+'\n')
    temp.replace(path)


def fraction(record):
    return Fraction(int(record['numerator']),int(record['denominator']))


def run_job(job,mode,output):
    folder=output/'jobs';folder.mkdir(exist_ok=True)
    dest=folder/(job['id']+'.'+mode+'.json')
    if dest.exists():
        result=json.loads(dest.read_text())
        assert result['status']=='completed'
    else:
        temp=dest.with_name(dest.name+'.partial')
        b=job['candidate']
        args=[str(ROOT/'build/hnw_endpoint'),mode,str(b['R']),str(b['start_round_index']),
              *b['input_difference'],*b['output_mask'],str(temp)]
        with (folder/(job['id']+'.'+mode+'.log')).open('w') as log:
            subprocess.run(args,check=True,stdout=log,stderr=log)
        result=json.loads(temp.read_text());temp.replace(dest)
    b=job['candidate']
    assert result['rounds']==b['R'] and result['start_round_index']==b['start_round_index']
    assert result['input_words'].split(';')==b['input_difference']
    assert result['output_words'].split(';')==b['output_mask']
    assert result['mode']==mode and result['rd']==1 and result['rl']==0
    assert result['completed_branches']==result['expected_first_DDT_branches']==len(result['branches'])
    assert sum((fraction(x['probability']) for x in result['branches']),Fraction(0))==1
    assert sum((fraction(x['contribution']) for x in result['branches']),Fraction(0))==fraction(result['correlation'])
    for x in result['branches']:
        assert fraction(x['probability'])*fraction(x['conditional_model_correlation'])==fraction(x['contribution'])
    assert abs(fraction(result['correlation']))<=1
    return {k:v for k,v in result.items() if k!='branches'}|{'artifact':str(dest.resolve())}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=ROOT/'results/hnw26_comparison')
    parser.add_argument('--workers',type=int,default=24)
    parser.add_argument('--modes',nargs='+',choices=['zero','full'],default=['zero','full'])
    parser.add_argument('--selection',choices=['all_rows','four_groups'],default='all_rows')
    args=parser.parse_args();output=args.output.resolve();output.mkdir(parents=True,exist_ok=True)
    source=ROOT/'results/datf_pools/merged/merged_results.json'
    data=json.loads(source.read_text())
    positives=[];zeros=[]
    for r in data['rows']:
        if r['selected_correlation'] and fraction(r['selected_correlation']):positives.append(r)
        elif r['best'] and r['first_correlation'] is not None and not fraction(r['first_correlation']):zeros.append(r)
    if args.selection=='four_groups':
        grouped={}
        for r in zeros:
            key=(r['family'],r['middle'])
            if key not in grouped or r['best']['score_weight']<grouped[key]['best']['score_weight']:grouped[key]=r
        zeros=list(grouped.values())
    order=['D1_M1_L1_pool','D1_M2_L1_pool','D1_M1_D3_chi2_to_3']
    positives.sort(key=lambda r:order.index(r['id']))
    jobs=[]
    for r in positives+zeros:
        b=r['selected'] or r['best']
        jobs.append({'id':r['id'],'category':'DATF_nonzero' if r in positives else 'DATF_zero',
                     'original_split':r['split'],'hnw_split':f'1+{r["R"]-1}+0',
                     'candidate':b,'datf_correlation':r['selected_correlation'] or r['first_correlation'],
                     'source_coverage':r['coverage'],'source_status':r['status']})
    assert len(positives)==3 and len(zeros)==(23 if args.selection=='all_rows' else 4)
    files=[source,ROOT/'src/hnw_endpoint.cpp',ROOT/'include/gaston.hpp',ROOT/'include/estimator.hpp',ROOT/'build/hnw_endpoint']
    config={'selection':args.selection,'modes':args.modes,'workers':args.workers,'jobs':jobs,
            'sources':[{'path':str(f),'sha256':hashlib.sha256(f.read_bytes()).hexdigest()} for f in files],
            'method':'HNW26 first-order, exact first DDT, rd=1 and rl=0; local block factorization',
            'full':'1024 coordinates per column for (base state, difference); first-round conditional base distribution retained; reference constants retained',
            'zero':'32 difference-mask coordinates per column; base-state-mask components discarded; constants have no effect in this model',
            'output_extension':'For masks spanning multiple columns, multiply local output components within each DDT branch, then sum signed branch contributions.',
            'first_L_omitted':True}
    cfg=output/'configuration.json'
    if cfg.exists():assert json.loads(cfg.read_text())==config,'Configuration changed; use a new output directory.'
    else:save(cfg,config)
    began=time.time();results={j['id']:dict(j,hnw={}) for j in jobs}
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        pending={pool.submit(run_job,j,m,output):(j['id'],m) for j in jobs for m in args.modes}
        for future in concurrent.futures.as_completed(pending):
            rid,mode=pending[future];result=future.result();results[rid]['hnw'][mode]=result
            done=sum(len(r['hnw']) for r in results.values())
            save(output/'progress.json',{'done':done,'total':len(pending),'elapsed_seconds':time.time()-began})
            print('DONE',done,'/',len(pending),rid,mode,str(fraction(result['correlation'])),flush=True)
    final={'status':'completed','method':config['method'],'selection':args.selection,'modes':args.modes,
           'elapsed_seconds':time.time()-began,'workers':args.workers,'results':list(results.values()),
           'existing_experiments':data['experiments']}
    save(output/'results.json',final)
    print('COMPLETED',len(jobs),'endpoints',round(final['elapsed_seconds'],3),'seconds',flush=True)


if __name__=='__main__':main()
