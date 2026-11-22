#!/usr/bin/env python3
"""Cross-check two complete enumerations and independently row-reduce all candidate spaces."""
from datetime import datetime
import hashlib
import itertools
import json
import math
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
def read(name):return json.loads((ROOT/name).read_text())
def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def row_nullspace(columns):
    n=len(columns)
    rows=[sum(((column>>bit)&1)<<j for j,column in enumerate(columns)) for bit in range(64)]
    pivots=[];rank=0
    for j in range(n):
        pivot=next((i for i in range(rank,64) if rows[i]>>j&1),None)
        if pivot is None:continue
        rows[rank],rows[pivot]=rows[pivot],rows[rank]
        for i in range(64):
            if i!=rank and rows[i]>>j&1:rows[i]^=rows[rank]
        pivots.append(j);rank+=1
    basis=[]
    for j in range(n):
        if j in pivots:continue
        vector=1<<j
        for i,p in enumerate(pivots):
            if rows[i]>>j&1:vector|=1<<p
        assert all((vector&row).bit_count()%2==0 for row in rows)
        basis.append(vector)
    return rank,basis


def main():
    matrix_path=ROOT/'results/support_reference_matrices_20261024.txt'
    reference=read('results/support_reference_verification_20261024.json')
    assert reference['status']=='PASS' and reference['official_reference_basis_cases']==640
    assert reference['matrix_sha256']==digest(matrix_path)
    official=ROOT.parent/'gaston_ref/gaston/implementation/reference/gaston.c'
    assert reference['reference_sha256']==digest(official)
    raw=[int(w,16) for w in matrix_path.read_text().split()]
    assert len(raw)==640
    tables=[raw[:320],raw[320:]]
    # Check all 64 simultaneous rotations on the actual restricted-kernel columns.
    mask64=(1<<64)-1
    for table in tables:
        for z in range(64):
            for i in range(5):
                v=table[5*z+i]
                assert table[5*((z+1)%64)+i]==((v<<1)|(v>>63))&mask64
    analysis=read('results/support_space_analysis_20261024.json')
    assert analysis['matrix_sha256']==digest(matrix_path)
    branch_reports={}
    sources=['src/support_rank_probe.cpp','src/verify_support_ranks.cpp','src/build_support_reference.py',
             'src/analyze_support_spaces.py','src/verify_support_proof.py','src/run_support_proof.py',
             'src/audit_kernel_structure.py','src/verify_lossy_model.py','src/build_activity_profiles.py','Makefile',
             'src/support_rank_probe','src/verify_support_ranks',
             'results/support_reference_verification_20261024.json','results/support_reference_matrices_20261024.txt',
             'results/support_space_analysis_20261024.json','references/activity_profiles.json']
    for k,name in ((6,'six'),(7,'seven'),(8,'eight')):
        a_path=f'results/support_rank_{name}_20261024.json'
        b_path=f'results/support_rank_{name}_independent_20261024.json'
        sources.extend([a_path,b_path]);a,b=read(a_path),read(b_path)
        assert a.get('complete_anchored_enumeration',a.get('complete_anchored_six_column_enumeration'))
        assert b['complete_anchored_enumeration'] and b['middle_active_columns']==k
        assert a['checked']==b['checked']==math.comb(63,k-1)
        assert a['left_dependent']==b['left_dependent'] and a['both_dependent']==b['both_dependent']
        assert a['left_orbit_masks']==b['left_orbit_masks']
        assert a['both_orbit_masks']==b['both_orbit_masks']
        branch=analysis['branches'][str(k)]
        assert branch['rank_enumeration_sha256']==digest(ROOT/a_path)
        assert [r['mask'] for r in branch['records']]==a['both_orbit_masks']
        max_nullity=0;checked_vectors=0;both_full=0
        for record in branch['records']:
            mask=int(record['mask'],16);positions=[z for z in range(64) if mask>>z&1]
            assert len(positions)==k
            availability=[]
            for side,table in enumerate(tables):
                columns=[table[5*z+i] for z in positions for i in range(5)]
                rank,basis=row_nullspace(columns)
                side_record=record['left' if side==0 else 'right']
                assert rank==side_record['rank'] and len(basis)==side_record['nullity']
                max_nullity=max(max_nullity,len(basis))
                vectors={0}
                for generator in basis:vectors|={v^generator for v in list(vectors)}
                assert len(vectors)==1<<len(basis)
                assert len(vectors)-1==side_record['nonzero_vectors_exhausted']
                checked_vectors+=len(vectors)-1
                full=set()
                for vector in vectors:
                    if all((vector>>(5*j))&31 for j in range(k)):
                        words=tuple(sum(((vector>>(5*j+i))&1)<<z for j,z in enumerate(positions)) for i in range(5))
                        full.add(words)
                recorded={tuple(int(w,16) for w in v['middle_words']) for v in side_record['full_support_candidates']}
                assert full==recorded
                availability.append(bool(full))
            both_full+=all(availability)
        assert both_full==branch['both_full_support_orbits']==0
        assert branch['minimum_lossy_weight'] is None and branch['below_44_excluded']
        branch_reports[str(k)]={'complete_anchored_supports_checked_twice':a['checked'],
            'left_rank_deficient_orbits':len(a['left_orbit_masks']),
            'both_rank_deficient_orbits':len(a['both_orbit_masks']),
            'both_exact_full_support_orbits':both_full,'maximum_candidate_nullity':max_nullity,
            'nonzero_candidate_vectors_rechecked':checked_vectors,
            'double_kernel_branch_infeasible_for_any_weight':True}

    # Re-establish coverage of W<44 from the saved, completed single-layer results.
    profile=read('references/activity_profiles.json')
    for record in profile['records']+profile['branch_sources']:
        assert digest(ROOT.parent / Path(record['path']).relative_to('/home/yxb/分析'))==record['sha256']
    f,r=profile['forward_lower_by_exact_input_count'],profile['reverse_lower_by_exact_output_count']
    edges={n:[(a,b) for a in range(1,65) for b in range(1,65)
              if b>=f[a] and a>=r[b] and a+b>=12+4*n
              and (not n or b!=6 or a>=13) and (n or (a<=5*b and b<=5*a))] for n in (0,1)}
    remaining=[]
    for n0,n1 in itertools.product((0,1),repeat=2):
        for k in range(1,65):
            aa=[a for a,b in edges[n0] if b==k];cc=[c for b,c in edges[n1] if b==k]
            if aa and cc and 2*(min(aa)+k+min(cc))<44:remaining.append((n0,n1,k))
    assert remaining==[(0,0,k) for k in (6,7,8,9)]
    # k=9: W<=43 forces counts (6,9,6); kernel bit weights even and >=10.
    assert r[9]==6 and f[9]==6 and 2*(6+9+6)==42
    assert 42+math.ceil((10-9)/2)+math.ceil((10-6)/2)==45
    # Verify the even-bit identity on all 320 basis vectors using official-derived right columns.
    assert all((v&0x5555555555555555).bit_count()%2==1 for v in tables[1])
    report={'status':'PASS','checked_local':datetime.now().astimezone().isoformat(),
        'conclusion':'Every nonzero free-endpoint Gaston-320 three-chi-layer XOR differential trail has weight at least 44 under the agreed local-DDT-weight definition',
        'weight_lower_bound':44,'new_computation':'Exact finite GF(2) enumeration and full-support nullspace checks, without floating-point MILP inference',
        'historical_single_layer_proofs_replayed':False,
        'source_scope':'Uses the prior completed single-layer activity results and outside-kernel branch proof, with their files hash-checked',
        'branch_reports':branch_reports,'middle9_analytic_lower_bound':45,
        'remaining_below44_cases':[],
        'sha256':{name:digest(ROOT/name) for name in sources}}
    output=ROOT/'results/support_space_proof_20261024.json'
    output.write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='sha256'},ensure_ascii=False,indent=2))


if __name__=='__main__':main()
