"""Paper-style address/factor-expression cost, explicitly excluding evaluation.

Uses the completed 430-class prefix screen. Rotation wrappers rename key
variables when sharing a quartic coefficient expression. The 3->2 numerical
batch evaluator is a separate schedule; its cost is not removed by compact
factor-expression construction.
"""
import json
import math
from pathlib import Path

from hsilu_7r_batch_audit import leading_value, rol

ROOT=Path(__file__).resolve().parents[1]

def rotate_term(term,s,n):
    return tuple(sorted(n*(p//n)+(p%n+s)%n for p in term))

def reference(term,n):
    canonical=min(rotate_term(term,s,n) for s in range(n))
    shift=next(s for s in range(n) if rotate_term(canonical,s,n)==tuple(term))
    return canonical,shift

def validate_rotation_wrappers():
    n=8
    data=json.loads((ROOT/'hsilu_quartic_screen/terms_8.json').read_text())
    coefficients=json.loads((ROOT/'hsilu_quartic_screen/validation_8.json').read_text())
    tables={tuple(data['classes'][row['id']]):int(row['coefficient_truth_table_hex'],16)
            for row in coefficients['details']}
    terms=[data['mapping'][f'observable_{i}']['terms'] for i in range(2)]
    refs=[[reference(tuple(t),n) for t in word] for word in terms]
    for lam in range(1<<n):
        combined=0
        for word in range(2):
            for bit in range(n):
                value=0
                for representative,shift in refs[word]:
                    key=rol(lam,-shift-bit,n)
                    value^=(tables[representative]>>key)&1
                combined|=value<<(n*word+bit)
        direct=0
        for x in range(1<<n):
            direct^=leading_value(x,x&lam,4,n)
        assert combined==direct&((1<<(2*n))-1)
    return {'word_bits':n,'rounds':4,'key_settings':1<<n,
            'nonce_points_per_key':1<<n,'output_bits':2*n,
            'all_rotated_quartic_expression_combinations_match_direct_cube':True}

def count(value,unit):
    return {'value':value,'log2':math.log2(value),'unit':unit}

def construction(q):
    a=math.comb(64,32)*math.comb(32,16)
    return {
        'quartic_representatives':q,
        'pair_list_construction':count(2*q*a,'paper-style logical table-access units'),
        'outer_product_construction':count(2*q*a,'paper-style logical table-access units'),
        'combined_construction':count(4*q*a,'paper-style logical table-access units'),
        'seven_round_equivalent':count(4*q*a/7,'seven-round equivalents, one logical access = one round'),
        'pair_list_memory':count(2*q*a*128,'bits; two lists per quartic, two 64-bit subset addresses per entry'),
        'outer_product_node_memory_upper_bound':count(q*math.comb(64,32)*128,'bits; two subset-address references per outer product, lower order'),
    }

def main():
    screen=json.loads((ROOT/'results/hsilu_7r_quartic_screen.json').read_text())
    assert screen['complete'] and screen['rotation_classes_retained']==430
    assert screen['rotation_classes_excluded']==0
    data=json.loads((ROOT/'hsilu_quartic_screen/terms_64.json').read_text())
    refs=[]
    for w in range(2):
        for term in data['mapping'][f'observable_{w}']['terms']:
            canonical,shift=reference(tuple(term),64)
            assert list(canonical) in data['classes']
            refs.append((w,canonical,shift))
    r16=sum(sum(math.gcd(i,d)==1 for i in range(1,d+1))*math.comb(64//d,16//d)
            for d in [1,2,4,8,16])//64
    report={
        'scope':'Construction of factored symbolic expressions; not all-key numerical evaluation or a full key-recovery complexity',
        'model':screen['model'],
        'validation':validate_rotation_wrappers(),
        'rotation_identity':'H_rot_s(t)(k) = H_t(rot_-s(k)); each shared expression is used with the required key-variable renaming',
        'word_zero_mapping_entries':len(refs),
        'all_output_mapping_entries':64*len(refs),
        'base_coefficient_recovery':count(r16*(1<<32),'five-round Hsilu calls'),
        'base_coefficient_recovery_seven_round_equivalent':count(r16*(1<<32)*5/7,'seven-round equivalents'),
        'base_coefficient_memory_compressed':count(320*r16*(1<<16),'bits'),
        'base_coefficient_memory_fully_materialized':count(320*math.comb(64,16)*(1<<16),'bits; optional rotation expansion to simplify table addressing'),
        'hsilu_430_classes':construction(430),
        'hsilu_no_quartic_orbit_sharing_2042_types':construction(2042),
        'ascon_121_paper_formula_recomputed':construction(121),
        'ascon_paper_reported':{'time_log2':95.96,'memory_bits_log2':104.85},
        'full_time_formula':'T = T_base + T_expression + T_all_key_evaluation + T_matching + L*2^64; L remains the actual signature bucket size',
        'data':count(1<<64,'distinct-nonce queries, each observing 128 bits'),
        'limitations':[
            'No additional pair deduplication or alternative factor pairing is claimed; this follows the paper\'s two lists per quartic convention.',
            'Logical table-access units are the paper\'s expression-construction model, not a count of physical 64-bit loads and stores.',
            'Memory uses two 64-bit subset addresses; table identity and variable-renaming descriptors are shared metadata. Allocator/hash overhead is not included.',
            'Base rotation recovery can optionally materialize all original subset tables. This lower-order stage does not change the dominant expression-memory bound.',
            '430 counts shared quartic types for 128 observable Hsilu bits; 121 is the paper\'s chosen Ascon output expression. Neither count establishes optimality.',
            'The previous shared ternary numerical evaluator remains a separate complete evaluation schedule; compact expressions alone do not prove a cheaper evaluator.'
        ]
    }
    path=ROOT/'results/hsilu_7r_factored_cost.json'
    path.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))

if __name__=='__main__':main()
