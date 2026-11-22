"""Reassess archived Gaston splice candidates; never evaluate full endpoints."""
import json, hashlib, sys, math
from pathlib import Path
from fractions import Fraction
HERE=Path(__file__).resolve().parent
G=HERE.parents[1]
ROOT=G.parent
sys.path.insert(0,str(G/'tools'))
from fixed_tail_search import column, linear, transpose, lat, S
def read(p):return json.loads(p.read_text())
def ws(a):
    if isinstance(a,dict):a=a['words']
    return [f'{int(w,16) if isinstance(w,str) else w:016x}' for w in a]
def state(a):return tuple(int(w,16) for w in ws(a))
def rec(q):return dict(numerator=str(q.numerator),denominator=str(q.denominator),log2_abs=math.log2(abs(q.numerator))-math.log2(q.denominator) if q else None)
def fraction(q):return Fraction(int(q['numerator']),int(q['denominator']))
def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2)+'\n')

def main():
    meta=read(G/'results/datf_pools/pools/pools.json');tails={r['id']:r for r in meta['paper_tails']}
    jobs={};rows=[];inventory=[]
    for t in tails.values():
        q=Fraction(1)
        for r in t['trail']:
            q*=math.prod(Fraction(lat(column(r['input'],j),column(r['output'],j)),32) for j in range(64))
        assert q*q==Fraction(1,2**t['weight_squared'])
        assert state(t['join_mask'])==transpose(t['trail'][0]['input'])
        assert state(t['output'])==tuple(t['trail'][-1]['output'])
        for a,b in zip(t['trail'],t['trail'][1:]):assert tuple(a['output'])==transpose(b['input'])
        inventory.append(dict(id=t['id'],rounds=t['rounds'],q_squared=rec(q*q),input=t['input'],output=t['output'],join_mask=t['join_mask']))
    def add(id,source,family,rd,rm,tail,delta,branches,extra=None):
        t=tails[tail];R=rd+rm+t['rounds']
        if R<5:return
        bs=[]
        for b in branches:
            middle=ws(b['middle']);mask=t['join_mask'];start=12-R+rd
            key=json.dumps([rm,start,middle,mask]);jid=hashlib.sha256(key.encode()).hexdigest()[:16]
            jobs[jid]=dict(id=jid,algorithm='Gaston',rounds=rm,start_round_index=start,delta_core=middle,lambda_core=mask)
            bs.append(dict(job_id=jid,probability=rec(b['p']),expected_old_middle=b.get('expected')))
        rows.append(dict(id=id,source=source,family=family,rounds=R,split=[rd,rm,t['rounds']],start_round_index=12-R,
            tail_id=tail,input_difference=ws(delta),output_mask=t['output'],q_squared=rec(Fraction(1,2**t['weight_squared'])),
            branches=bs,aggregation='single retained trail' if len(bs)==1 else 'sum of retained differential events sharing endpoints and tail',extra=extra))
    src='datf_pools/paper_full24/search_results.json'
    for row in read(G/'results'/src)['rows']:
        c=row['best']
        if c:
            add(row['id'],src,'pool_single',c['r_d'],c['r_m'],c['linear_id'],c['input_difference'],
                [dict(middle=c['middle_input_difference'],p=fraction(c['differential_probability']),expected=c['middle_DATF_correlation'])],
                extra=dict(original_candidate=c,archived_search_exhaustive=row['exhaustive']))
    src='fixed_paper_tail/search.json'
    for r in read(G/'results'/src)['reports']:
        if not r['top_candidates']:continue
        c=r['top_candidates'][0]
        add('cluster1_'+r['id'],src,'single_column_cluster',1,1,r['id'],c['input_difference'],
            [dict(middle=b['middle_chi_input_difference'],p=Fraction(b['ddt_count'],b['denominator']),expected=b['middle_uniform_dlct']) for b in c['prefix']['branches']],extra=dict(candidate=c))
    src='pq_kernel_six_all/search.json';data=read(G/'results'/src)
    for r in data['reports']:
        b=r['best']
        if not b:continue
        gr=data['groups'][r['connection_group']]
        scan=next(s for s in gr['scans'] if (s['orbit'],s['rotation'])==(b['orbit'],b['rotation']))
        add('cluster2_'+r['id'],src,'pq_cluster',2,1,r['id'],b['d0'],
            [dict(middle=c['a3'],p=Fraction(1,2**24),expected=c['middle_correlation']) for c in scan['connections']],
            extra=dict(orbit=b['orbit'],rotation=b['rotation'],first_output=b['b1'],second_input=b['a2'],connections=scan['connections']))
    src='fixed_tail_extension/results.json'
    for r in read(G/'results'/src)['rows']:
        if not r['candidates']:continue
        c=r['candidates'][0]
        if c['differential_weight']>64:continue
        tail=r['id'].split('_',2)[2]
        add('extended_'+r['id'],src,'fixed_tail_backward',r['r_d'],r['r_m'],tail,c['input_difference'],
            [dict(middle=c['middle_input_difference'],p=Fraction(1,2**c['differential_weight']),expected=c['middle_correlation'])],
            extra=dict(candidate=c))
    (HERE/'raw').mkdir(exist_ok=True)
    save(HERE/'manifest.json',dict(jobs=list(jobs.values())))
    save(HERE/'candidates.json',dict(rows=rows,linear_inventory=inventory,
         scope='Best archived candidate per search family plus explicitly archived same-endpoint clusters; no new exhaustive search.',
         first_L_omitted=True))
    print('Prepared',len(rows),'candidates;',len(jobs),'distinct middle jobs; max rounds',max(r['rounds'] for r in rows))
if __name__=='__main__':main()
