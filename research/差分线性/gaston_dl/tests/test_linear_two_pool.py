#!/usr/bin/env python3
"""Independent forward-matrix coverage and LAT/constant/matcher checks."""
import argparse
import collections
import csv
import hashlib
import itertools
import json
import subprocess
import sys
import tempfile
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from datf_pool_build import read_pool, state, words
from datf_pool_search import linear_candidates
from fixed_tail_search import RC, S, column, embed, linear, rotate, transpose
from linear_seed_scan import MASK, active
from linear_two_pool import POOL_KEY, FAMILY


def pack(s):
    return sum(w << (64*i) for i,w in enumerate(s))


def unpack(z):
    return tuple((z >> (64*i)) & MASK for i in range(5))


def exhaustive_forward_matrix_scan():
    # Derive every entry of L^T from the independent forward linear function.
    images = [pack(linear(embed(j, 1 << i))) for i in range(5) for j in range(64)]
    basis = [sum(((v >> k) & 1) << p for p,v in enumerate(images)) for k in range(320)]
    table = []
    for j in range(64):
        values = [0]*32
        for a in range(1, 32):
            low = a & -a
            values[a] = values[a ^ low] ^ basis[64*(low.bit_length()-1)+j]
        table.append(values)
    found, checked = set(), 0
    for j,k in itertools.combinations(range(64), 2):
        for a,b in itertools.product(range(1,32), repeat=2):
            z = table[j][a] ^ table[k][b]
            columns = (z | (z>>64) | (z>>128) | (z>>192) | (z>>256)) & MASK
            if columns.bit_count() == 2:
                found.add((z, pack(tuple(x|y for x,y in zip(embed(j,a),embed(k,b))))))
            checked += 1
    return found, checked


def check_native_matcher(pools, meta):
    ds = []
    for key, wanted, shift in [('differential_one','D1_00016',1),
                                ('differential_pq_two','D2_K09_1258',29)]:
        d = next(x for x in pools[key] if x['id']==wanted)
        ds.append({**d, 'input': words(rotate(state(d['input']),shift)),
                   'output': words(rotate(state(d['output']),shift)), 'orbit_size': 2})
    candidates = linear_candidates(pools, meta)
    ls = [next(l for l in candidates if l['id']=='L1_01977')]
    ls += [next(l for l in candidates if l['family']==FAMILY and l['weight']==w) for w in (4,5,6)]
    ls += [next(l for l in candidates if l['id']=='D3_chi3_to_4')]
    # Adapter-only control: re-use a known nonzero middle mask in the L2
    # family to exercise aggregation. This is not installed as a real trail.
    ls.append({**ls[0], 'id':'test_L2_adapter_control', 'family':FAMILY,
               'rounds':2, 'weight_squared':8})
    jobs, keys = [], []
    for d,l,m,shift in itertools.product(ds,ls,(1,2),range(2)):
        delta = linear(rotate(state(d['output']),shift))
        rc = RC[11-l['rounds']]
        jobs.append(str(m)+' '+' '.join(f'{v:x}' for v in (*delta,*state(l['join_mask']),rc))+'\n')
        keys.append((d['id'],l['id'],m,shift))
    output = subprocess.check_output([str(ROOT/'build/datf_middle_cli')], input=''.join(jobs), text=True)
    oracle = {}
    for key,line in zip(keys,output.splitlines(),strict=True):
        sign,w,_ = line.split()
        oracle[key] = Fraction(int(sign), 1<<int(w))
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        (tmp/'d.txt').write_text(''.join(f'{d["id"]} D{d["rounds"]} {d["weight"]} {d["orbit_size"]} '
                                +' '.join(d['input']+d['output'])+'\n' for d in ds))
        (tmp/'l.txt').write_text(''.join(f'{l["id"]} {l["family"]} {l["rounds"]} {l["weight_squared"]} '
                                +' '.join(l['join_mask']+l['output'])+'\n' for l in ls))
        subprocess.run([str(ROOT/'build/datf_pool_match'), str(tmp/'d.txt'), str(tmp/'l.txt'),
                        str(tmp/'match'), '2','10000','100','10000','100'],
                       check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        summaries = list(csv.DictReader((tmp/'match_summary.tsv').open(),delimiter='\t'))
        hits = list(csv.DictReader((tmp/'match_hits.tsv').open(),delimiter='\t'))
        assert len(summaries)==12  # 2 differential families x 2 middles x 3 tail groups
        for s in summaries:
            d = next(d for d in ds if 'D'+str(d['rounds'])==s['family'])
            selected = [l for l in ls if (l['family']+'_pool' if l['family']!='paper' else l['id'])==s['tail']]
            values = [abs(oracle[d['id'],l['id'],int(s['middle']),shift])
                      / (1 << (d['weight']+l['weight_squared']))
                      for l,shift in itertools.product(selected, range(2))]
            best = max(values)
            assert (Fraction(1,1<<int(s['best_weight'])) if s['best_weight'] else 0)==best
            if not best:
                assert s['exhausted']=='1' and int(s['evaluated'])==2*len(selected)
        for h in hits:
            d = next(d for d in ds if d['id']==h['diff_id'])
            l = next(l for l in ls if l['id']==h['linear_id'])
            want = oracle[d['id'],l['id'],int(h['middle']),int(h['rotation'])]
            assert Fraction(int(h['sign']),1<<int(h['middle_weight']))==want
            assert Fraction(int(h['sign']),1<<int(h['weight']))==want/(1<<(d['weight']+l['weight_squared']))
        assert any(h['tail']==FAMILY+'_pool' for h in hits)
    return len(oracle), len(summaries)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--pools',type=Path,default=ROOT/'results/datf_pools/pools')
    args=p.parse_args()
    meta=json.loads((args.pools/'pools.json').read_text())
    pools={key:read_pool(args.pools/desc['file']) for key,desc in meta['pools'].items()}
    records=pools[POOL_KEY]
    audit=json.loads((args.pools/(POOL_KEY+'_audit.json')).read_text())
    scanned,count=exhaustive_forward_matrix_scan()
    expected=set()
    for r in audit['bridges']:
        a=state(r['first_output']['words']); b=state(r['second_input']['words'])
        expected.update((pack(rotate(a,s)),pack(rotate(b,s))) for s in range(64))
    assert scanned==expected and len(scanned)==640
    assert count==1937376
    old=json.loads((ROOT/'results/linear_seeds/two_to_two.json').read_text())
    assert {(state(x['alpha']['words']),state(x['beta']['words'])) for x in old['seeds']}=={
        (state(x['first_output']['words']),state(x['second_input']['words'])) for x in audit['bridges']}
    # Independently evaluate the complete 32-point Walsh sum for each LAT entry.
    lat=[[sum(1 if ((a&x).bit_count()+(b&S[x]).bit_count())%2==0 else -1 for x in range(32))
          for b in range(32)] for a in range(32)]
    seen, groups = set(), collections.defaultdict(set)
    for r in records:
        a1,b1,a2,b2=map(lambda k:state(r[k]),('input','first_output','second_input','output'))
        signature=(a1,b1,a2,b2)
        assert signature not in seen
        seen.add(signature); groups[r['bridge_id']].add((a1,b2))
        assert transpose(a2)==b1 and transpose(a1)==state(r['join_mask'])
        assert all(active(x)==2 for x in signature)
        corr=Fraction(1)
        for pair, round_record in zip(((a1,b1),(a2,b2)),r['trail']):
            a,b=pair
            local=Fraction(1)
            for j in range(64):
                local*=Fraction(lat[column(a,j)][column(b,j)],32)
            assert local==Fraction(round_record['sign'],1<<(round_record['w_squared']//2))
            corr*=local
        assert corr==Fraction(r['sign'],1<<r['weight'])
        assert r['weight_squared']==2*r['weight']
        for s, sign in enumerate(r['sign_by_start_round']):
            phase=sum(((v & RC[s+i]).bit_count()) for i,v in enumerate((a1[0],a2[0])))%2
            assert sign==r['sign']*((-1)**phase)
        assert r['orbit_size']==64
    for bridge in audit['bridges']:
        b1,a2=state(bridge['first_output']['words']),state(bridge['second_input']['words'])
        possibilities=[]
        for fixed,before in [(b1,True),(a2,False)]:
            cols=[j for j in range(64) if column(fixed,j)]
            opts=[]
            for j in cols:
                v=column(fixed,j)
                opts.append([m for m in range(32) if (lat[m][v] if before else lat[v][m])])
            possibilities.append({tuple(x|y for x,y in zip(embed(cols[0],u),embed(cols[1],v)))
                                  for u,v in itertools.product(*opts)})
        assert len(possibilities[0])==16 and len(possibilities[1])==100
        assert groups[bridge['id']]==set(itertools.product(*possibilities))
    assert len(records)==16000
    assert collections.Counter(r['weight'] for r in records)=={4:640,5:5120,6:10240}
    oracle_count, summary_count=check_native_matcher(pools,meta)
    # Existing pool data stay byte-for-byte identical to the saved search snapshot.
    old_meta=json.loads((ROOT/'results/datf_pools/search_final/results.json').read_text())['pool_metadata']
    for key in ('differential_one','differential_pq_two','linear_one'):
        assert hashlib.sha256((args.pools/meta['pools'][key]['file']).read_bytes()).hexdigest()==old_meta['pools'][key]['sha256']
    validation={'status':'PASS','independent_forward_matrix_two_column_patterns':count,
                'expanded_bridges':len(scanned),'complete_trails_checked':len(records),
                'nonzero_LAT_rounds_checked':2*len(records),'constant_start_signs_checked':11*len(records),
                'native_middle_oracle_cases':oracle_count,'native_matcher_summary_groups':summary_count,
                'old_seed_classes_reproduced':10,'old_pool_files_hash_unchanged':True,
                'scope':'Coverage and individual trail correctness; not a full-correlation validation.'}
    (args.pools/(POOL_KEY+'_validation.json')).write_text(json.dumps(validation,indent=2)+'\n')
    print(json.dumps(validation))


if __name__=='__main__':
    main()
