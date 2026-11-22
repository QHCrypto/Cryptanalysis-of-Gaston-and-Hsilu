"""Evaluate only the four middle jobs. Whole endpoint evidence is reused.

Usage: python3 evaluate_middle.py [--reuse]
One-chi HNW is exactly the local DLCT product: its two-wise LAT operator
reduces to sum_u LAT(u,v)^2 (-1)^(u.delta). The factorization across the
64 independent uniform input columns is exact for a single chi layer.
"""
import argparse
import hashlib
import json
import math
import subprocess
import sys
import time
from fractions import Fraction
from pathlib import Path

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
sys.path.insert(0,str(ROOT/'gaston_dl/tools'))
from fixed_tail_search import S, column, embed, lat, linear, transpose

def record(q):
    return dict(numerator=str(q.numerator),denominator=str(q.denominator),
                log2_abs=math.log2(abs(q.numerator))-math.log2(q.denominator) if q else None)
def fraction(q): return Fraction(int(q['numerator']),int(q['denominator']))
def save(path,obj): path.write_text(json.dumps(obj,ensure_ascii=False,indent=2)+'\n')
def dlct(a,v):
    return Fraction(sum(1 if ((v&(S[x]^S[x^a])).bit_count()&1)==0 else -1
                        for x in range(32)),32)
def one_chi(delta,mask):
    return math.prod(dlct(column(delta,j),column(mask,j)) for j in range(64))

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--reuse',action='store_true');args=ap.parse_args()
    jobs=json.loads((OUT/'middle_manifest.json').read_text())['jobs']
    worker=ROOT/'DATF算法3复核_20261021/worker.py'
    # Exhaustive check of the one-chi HNW operator against direct enumeration.
    checked=0
    for a in range(32):
        for v in range(32):
            via_lat=sum(Fraction(lat(u,v)**2,1024)*((-1)**((u&a).bit_count()&1)) for u in range(32))
            assert via_lat==dlct(a,v)
            checked+=1
    for i,j in enumerate(jobs):
        p=OUT/'raw'/f"{j['id']}_DATF.json"
        if not args.reuse or not p.exists():
            subprocess.run([sys.executable,str(worker),'--manifest',str(OUT/'middle_manifest.json'),
                '--index',str(i),'--output',str(p),'--seconds','120'],check=True)
        if j['rounds']==1:
            started=time.perf_counter()
            delta=tuple(int(w,16) for w in j['delta_core'])
            mask=tuple(int(w,16) for w in j['lambda_core'])
            q=one_chi(delta,mask)
            result=dict(**j,status='completed',method='HNW26 one-chi two-wise LAT operator',
                backend='exact local DLCT product; analytically sums all first-DDT branches',
                mode='full and zero coincide for one chi',correlation=record(q),
                elapsed_seconds=time.perf_counter()-started,local_operator_checks=checked)
            assert fraction(json.loads(p.read_text())['correlation'])==q
            save(OUT/'raw'/f"{j['id']}_HNW.json",result)
        else:
            for mode in ['zero','full']:
                p=OUT/'raw'/f"{j['id']}_HNW_{mode}.json"
                if not args.reuse or not p.exists():
                    subprocess.run([str(ROOT/'hsilu_dl/build/hnw_endpoint'),mode,str(j['rounds']),
                        str(j['start_round_index']),*j['delta_core'],*j['lambda_core'],str(p)],
                        check=True,timeout=180)
    sp=OUT/'raw/H5_prefix3_DATF_shannon.json'
    if not args.reuse or not sp.exists():
        subprocess.run([sys.executable,str(worker),'--manifest',str(OUT/'middle_manifest.json'),
            '--index','3','--output',str(sp),'--seconds','30','--input-shannon'],check=True)
    # Every first-round DDT branch is checked, even where the archived search
    # kept only the nonzero one-chi connections. No full four-round evaluation.
    cuts=json.loads((OUT/'cuts.json').read_text())
    audits={}
    for code,a,col in [('G4c',16,1),('G4d',4,0)]:
        j=next(j for j in jobs if j['id']==cuts[code]['branches'][0]['job_id'])
        mask=tuple(int(w,16) for w in j['lambda_core']);branches=[];total=Fraction(0)
        for b in range(32):
            count=sum(S[x]^S[x^a]==b for x in range(32))
            if not count: continue
            beta=embed(col,b);mid=linear(beta);p=Fraction(count,32);q=one_chi(mid,mask)
            branches.append(dict(first_output_pattern=b,probability=record(p),
                middle_input_difference=[f'{w:016x}' for w in mid],middle_correlation=record(q)))
            total+=p*q
        assert sum(fraction(b['probability']) for b in branches)==1
        assert len([b for b in branches if fraction(b['middle_correlation'])])==len(cuts[code]['branches'])
        audits[code]=dict(branches=branches,signed_sum=record(total))
    assert fraction(audits['G4c']['signed_sum'])==Fraction(1,2)
    assert fraction(audits['G4d']['signed_sum'])==-Fraction(1,4)

    # Recheck the archived two-chi tails: all LAT entries and boundary maps.
    paper=json.loads((ROOT/'gaston_dl/results/fixed_paper_tail/search.json').read_text())
    tail=next(r for r in paper['reports'] if r['id']=='D3_chi2_to_3')
    ts=tail['tail_rounds'];signed=Fraction(1)
    for r in ts:
        signed*=math.prod(Fraction(lat(column(r['input'],j),column(r['output'],j)),32) for j in range(64))
    assert signed*signed==Fraction(1,2**38)
    assert tuple(ts[0]['output'])==transpose(ts[1]['input'])
    assert tuple(int(w,16) for w in tail['join_mask']['words'])==transpose(ts[0]['input'])

    sys.path.insert(0,str(ROOT/'hsilu_dl/tools'));import hsilu
    h=json.loads((ROOT/'hsilu_dl/results/research_summary.json').read_text())['five_candidates'][0]
    rs=[h['head']['first_chi'],h['last_chi']];hm=[];signed=Fraction(1)
    for r in rs:
        x=tuple(int(w,16) for w in r['input']['words']);y=tuple(int(w,16) for w in r['output']['words'])
        hm.append((x,y));signed*=math.prod(Fraction(lat(column(x,j),column(y,j)),32) for j in range(64))
    assert signed*signed==Fraction(1,2**10)
    assert hm[0][1]==hsilu.transpose(hm[1][0])
    assert tuple(int(w,16) for w in h['head']['join_mask']['words'])==hsilu.transpose(hm[0][0])
    main=json.loads((OUT/'raw/H5_prefix3_DATF.json').read_text())
    independent=json.loads((OUT/'raw/H5_prefix3_DATF_shannon.json').read_text())
    assert main['status']==independent['status']=='completed'
    assert fraction(main['correlation'])==fraction(independent['correlation'])==0
    for mode in ['full','zero']:
        x=json.loads((OUT/'raw'/f'H5_prefix3_HNW_{mode}.json').read_text())
        assert x['status']=='completed' and fraction(x['probability_mass'])==1
        assert x['completed_branches']==8 and fraction(x['correlation'])==0
    save(OUT/'validation.json',dict(status='PASS',one_chi_operator_coordinates_checked=checked,
        first_round_branch_audits=audits,tail_checks=['Gaston LAT product and linear masks','Hsilu LAT product and linear masks'],
        H5_algorithm3_independent_order_check='full expansion and exact input-Shannon order agree',
        H5_hnw_modes='full and zero agree, each covers all 8 first-DDT branches',
        new_whole_low_correlation_jobs=0))
    print('PASS:',checked,'one-chi operator identities; cuts, tails, H5 cross-checks.')

if __name__=='__main__': main()
