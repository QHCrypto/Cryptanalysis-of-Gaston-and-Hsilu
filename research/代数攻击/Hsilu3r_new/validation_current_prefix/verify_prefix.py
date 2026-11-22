#!/usr/bin/env python3
"""Audit the new SVG at the manuscript prefix and a nearby legal prefix.

All path coordinates are MSB-first. Constants 69/5a/4b are included in
the native computations. Original search and analysis programs are unchanged.
"""
from pathlib import Path
import collections, hashlib, json, random, re, sys

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent))
import analyze_hsilu_x3_conditions as alg
sys.path.insert(0,str(ROOT/'代数攻击/Hsilu_3r'))
import analyze_three_round_path as ref
from verify_quadratic_linearization import evaluate,substitute_background
sys.path.insert(0,str(ROOT/'nist_compatibility_20261109'))
import nist_modes as nist
sys.path.insert(0,str(ROOT/'Gaston_Hsilu_CTP/paper_rewriting_output/hsilu_figure_correction'))
import make_conditional_figures as figs

D=json.loads((HERE/'source_path.json').read_text())
GRID=json.loads((HERE/'source_grid.json').read_text())
U=D['free_variables_by_lane']['0']
BG=[j for j in range(64) if j not in U]
SELECT=D['selected_outputs_by_lane']['0']
COND={st:[(int(i),j,v) for key,v in [('conditions_eq_0',0),('conditions_eq_1',1)]
          for i,js in D[key].get(st,{}).items() for j in js] for st in ['X1','X2','X3']}
TEXT=(HERE/'x3_coordinate_linearization.txt').read_text().split('One minimum set of additional X2 guesses:')[1].split('Linearized equations')[0]
GUESS=[tuple(map(int,v)) for v in re.findall(r'X2\[(\d+)\]\[(\d+)\]',TEXT)]

def first_system(state):
    zero=ref.ref.rnd([0]+state[1:],0x69)
    def delta(j):
        return [a^b for a,b in zip(zero,ref.ref.rnd([1<<(63-j)]+state[1:],0x69))]
    cols=[delta(j) for j in BG]
    mid=[[((zero[i]>>(63-j))&1)|sum(((v[i]>>(63-j))&1)<<(k+1) for k,v in enumerate(cols)) for j in range(64)] for i in range(5)]
    rows=[mid[i][j]>>1 for i,j,v in COND['X2']]
    rhs=[(mid[i][j]&1)^v for i,j,v in COND['X2']]
    p,basis,rk=ref.system(rows,rhs,len(BG))
    return mid,rows,rhs,p,basis,rk,{j:delta(j) for j in U}

def exact_x3_polys():
    fixed=alg.collect_fixed_x2(D)
    y={}
    for i in range(5):
        for j in range(64):
            def inp(row):
                a=alg.variable_poly((row,j))
                if row==0 and 0x5a>>(63-j)&1:a^=alg.constant_poly(1)
                return alg.substitute(a,fixed)
            a,b,c=inp(i),inp((i+1)%5),inp((i+2)%5)
            y[i,j]=a^c^alg.multiply_poly(b,c)
    result=[]
    for i,j,v in COND['X3']:
        p=alg.constant_poly(v)
        for pos in alg.hsilu.MRM_SUPPORTS[i,j]:p^=y[pos]
        result.append(p)
    return result

def substitute_x2(poly,forms):
    out=set()
    for term in poly:
        t={0}
        for i,j in term:t=ref.multiply(t,ref.affine_poly(forms[i][j]))
        out^=t
    return out

def independent_indices(rows):
    pivots={};selected=[]
    for index,row in enumerate(rows):
        v=row
        while v:
            k=v.bit_length()-1
            if k in pivots:v^=pivots[k]
            else:pivots[k]=v;selected.append(index);break
    return selected

def audit_original(words):
    s=nist.prefix_state(words,'hsilu',3)
    mid,rows,rhs,p,basis,rk,deltas=first_system(s)
    violations=[{'lane':i,'msb_bit':j,'native_bit':63-j,'required':v,'actual':s[i]>>(63-j)&1}
                for i,j,v in COND['X1'] if (s[i]>>(63-j)&1)!=v]
    dependencies={str(j):[[i,k] for i,k,v in COND['X2'] if deltas[j][i]>>(63-k)&1] for j in U}
    return {'prefix_hex':[f'{v:016x}' for v in words],'state_hex':[f'{v:016x}' for v in s],
            'input_violations':violations,'valid_for_new_path':not violations,
            'X2_rank_at_zero_active':rk,'X2_consistent_at_zero_active':p is not None,
            'X2_conditions_depending_on_active_bits':dependencies,
            'combined_X2_and_25_guesses_rank_at_zero_active':ref.ref.rank(rows+[mid[i][j]>>1 for i,j in GUESS])}

def main():
    grids={st:[[figs.PALETTE[c] for c in row] for row in grid] for st,grid in GRID.items()}
    structural=figs.check(3,grids)
    original_words=[0x63c92946c6b0ac2c,0x101be7b317c13375]
    original=audit_original(original_words)
    chosen=json.loads((HERE/'new_prefix_candidate.json').read_text())
    words=[int(x,16) for x in chosen['prefix_hex']]
    s=nist.prefix_state(words,'hsilu',3)
    assert all((s[i]>>(63-j)&1)==v for i,j,v in COND['X1'])
    mid,rows,rhs,part,basis,rk,deltas=first_system(s)
    assert part is not None and rk==32
    for delta in deltas.values():
        assert all((delta[i]>>(63-j)&1)==0 for i,j,v in COND['X2'])
        assert all((delta[i]>>(63-j)&1)==0 for i,j in GUESS)
    reduced=[[ref.restrict_affine(x,part,basis) for x in row] for row in mid]
    guessforms=[reduced[i][j] for i,j in GUESS]
    guessrows=[x>>1 for x in guessforms]
    independent=independent_indices(guessrows)
    forms=[guessrows[i] for i in independent]
    _,kernel,grank=ref.system(forms,[0]*len(forms),len(basis))
    assert grank==len(forms)
    lifts=[ref.system(forms,[int(k==i) for k in range(grank)],len(basis))[0] for i in range(grank)]
    # Left kernel describes the two dependencies among actual X2 guess values.
    transpose=[sum(((row>>j)&1)<<i for i,row in enumerate(guessrows)) for j in range(len(basis))]
    _,relations,_=ref.system(transpose,[0]*len(transpose),len(GUESS))
    relations_report=[]
    for rel in relations:
        indices=[i for i in range(len(GUESS)) if rel>>i&1]
        value=sum(guessforms[i]&1 for i in indices)&1
        relations_report.append({'guess_indices':indices,'coordinates_msb':[GUESS[i] for i in indices],'rhs':value})
    exact=exact_x3_polys()
    initial=[alg.x3_condition_polynomial(i,j,v,alg.collect_fixed_x2(D)) for i,j,v in COND['X3']]
    corrections=[alg.format_polynomial(a^b) for a,b in zip(exact,initial)]
    # Constants alter at most linear terms: the coordinate cover is unchanged.
    assert all(all(len(m)<=1 for m in a^b) for a,b in zip(exact,initial))
    assert all(any(v in GUESS for v in term) for p in exact for term in p if len(term)==2)
    polys=[substitute_x2(p,reduced) for p in exact]
    for poly in polys:
        for i,u in enumerate(kernel):
            for v in kernel[:i]:
                assert evaluate(poly,u^v)^evaluate(poly,u)^evaluate(poly,v)^evaluate(poly,0)==0
    transformed=[[ref.restrict_affine(x,0,lifts+kernel) for x in row] for row in reduced]
    branched=[substitute_x2(p,transformed) for p in exact]
    free=len(kernel)
    # q(b,w)=rhs(b)+A(b)w. Prepare exact Gray-code updates for full enumeration.
    constants=[];changes=[];neighbors=[]
    for poly in branched:
        constant=0;change=[0]*grank;neighbor=[0]*grank
        for monomial in poly:
            bits=[i for i in range(grank+free) if monomial>>i&1]
            if not bits:constant^=1<<free
            elif len(bits)==1:
                j=bits[0]
                if j<grank:change[j]^=1<<free
                else:constant^=1<<(j-grank)
            else:
                assert len(bits)==2
                i,j=bits
                assert i<grank,'Quadratic term remains in branch variables'
                if j<grank:neighbor[i]^=1<<j;neighbor[j]^=1<<i
                else:change[i]^=1<<(j-grank)
        constants.append(constant);changes.append(change);neighbors.append(neighbor)
    inputs=[f'{grank} {free} {len(polys)}']
    inputs.append(' '.join(map(str,constants)))
    for change,neighbor in zip(changes,neighbors):
        inputs.append(' '.join(map(str,change)));inputs.append(' '.join(map(str,neighbor)))
    (HERE/'branch_system.txt').write_text('\n'.join(inputs)+'\n')
    rng=random.Random(20261111);samples=4096;consistent=backgrounds=0
    branchhist=collections.Counter();outputhist=collections.Counter();witnesses=[];planted=[]
    for trial in range(samples):
        b=rng.getrandbits(grank)
        augmented=[c for c in constants]
        for k,poly in enumerate(branched):
            aug=evaluate(poly,b)<<free
            for j in range(free):aug|=(evaluate(poly,b^(1<<(grank+j)))^evaluate(poly,b))<<j
            augmented[k]=aug
        br=[x&((1<<free)-1) for x in augmented];bv=[x>>free for x in augmented]
        sol,kbasis,r=ref.system(br,bv,free)
        # The solver returns early on a contradiction; compute the full
        # coefficient rank separately for an accurate histogram.
        r=ref.ref.rank(br)
        branchhist[f'{r}:{int(sol is not None)}']+=1
        if sol is None:continue
        consistent+=1
        for assignment in range(1<<len(kbasis)):
            w=substitute_background(assignment,sol,kbasis)
            joined=b|(w<<grank)
            assert all(evaluate(poly,joined)==0 for poly in branched)
            t=substitute_background(b,0,lifts)^substitute_background(w,0,kernel)
            assert all(evaluate(poly,t)==0 for poly in polys)
            g=substitute_background(t,part,basis)
            x0=sum(((g>>k)&1)<<(63-j) for k,j in enumerate(BG))
            outputs=[];states=[]
            for active in range(4):
                native=x0^sum(((active>>k)&1)<<(63-j) for k,j in enumerate(U))
                x2=ref.ref.rnd([native]+s[1:],0x69);x3=ref.ref.rnd(x2,0x5a)
                assert all((x2[i]>>(63-j)&1)==v for i,j,v in COND['X2'])
                assert all((x3[i]>>(63-j)&1)==v for i,j,v in COND['X3'])
                z=ref.ref.rnd(x3,0x4b);outputs.append(z[0]);states.append(z)
            assert all((outputs[0]^outputs[1]^outputs[2]^outputs[3])>>(63-j)&1==0 for j in SELECT)
            orows=[((outputs[0]^outputs[1])>>(63-j)&1)|(((outputs[0]^outputs[2])>>(63-j)&1)<<1) for j in SELECT]
            outputhist[ref.ref.rank(orows)]+=1;backgrounds+=1
            if len(witnesses)<8:witnesses.append({'branch':b,'free_assignment':w,'x0_hex':f'{x0:016x}','branch_rank':r})
            if len(planted)<32:
                active=rng.randrange(4);native=x0^sum(((active>>k)&1)<<(63-j) for k,j in enumerate(U));m=native^s[0]
                if not m:continue
                length=m.bit_length()-1;message=words[0]|words[1]<<64|(m^(1<<length))<<128
                h=nist.xof_bits(message,128+length,'hsilu',3)
                assert h==states[active][0]|nist.permute(states[active],'hsilu',3)[0]<<64
                target=[(outputs[active]^outputs[0])>>(63-j)&1 for j in SELECT]
                op,ok,_=ref.system(orows,target,2);assert op is not None
                aset=[substitute_background(v,op,ok) for v in range(1<<len(ok))]
                assert active in aset
                found=[v for v in aset if (states[v][0]|nist.permute(states[v],'hsilu',3)[0]<<64)==h]
                assert active in found
                planted.append({'active':active,'tail_length':length,'digest_hex':f'{h:032x}'})
    result={'status':'PASS','coordinate_convention':'all path coordinates MSB-first; native bit = 63-j',
            'source_svg':D['source_svg'],'source_sha256':hashlib.sha256((ROOT / Path(D['source_svg']).relative_to('/home/yxb/分析')).read_bytes()).hexdigest(),
            'structural_cells_verified':structural['verified_cells'],'round_constants':['69','5a','4b'],
            'counts':{'active':len(U),'X1':len(COND['X1']),'X2':len(COND['X2']),'X3':len(COND['X3']),
                      'coordinate_guesses':len(GUESS),'budget_sum':len(U)+len(COND['X2'])+len(COND['X3'])+len(GUESS)},
            'original_prefix':original,
            'new_prefix':{'prefix_hex':chosen['prefix_hex'],'state_hex':[f'{v:016x}' for v in s],
                          'increment_from_original_M1':chosen['increment'],'all_input_conditions_satisfied':True,
                          'X2_rank':rk,'X2_consistent':True,'background_dimension_after_X2':len(basis),
                          'guessed_coordinates_msb':GUESS,'independent_guess_indices':independent,
                          'independent_guess_rank':grank,'combined_rank_before_X3':rk+grank,
                          'remaining_variables_before_X3':free,'guess_dependencies':relations_report,
                          'all_branches_affine':True,'guesses_and_conditions_independent_of_active_variables':True},
            'round_constant_corrections_to_analyzer':corrections,
            'sampled_branches':samples,'consistent_sampled_branches':consistent,'sampled_backgrounds':backgrounds,
            'sampled_branch_rank_consistency_histogram':dict(branchhist),'sampled_output_rank_histogram':dict(outputhist),
            'sampled_backgrounds_per_branch':backgrounds/samples,'witnesses':witnesses,'planted_targets':planted,
            'scope':'original prefix rejected; exact all-branch affinity and sampled concrete verification for the displayed replacement prefix; no global prefix-density or full 128-bit search-time claim'}
    (HERE/'prefix_validation.json').write_text(json.dumps(result,indent=2)+'\n')
    algebra={'source_path':D,'background_particular':part,'background_basis':basis,
             'independent_guess_lifts':lifts,'independent_guess_kernel':kernel,
             'X3_polynomials_after_X2':[sorted(p) for p in polys],
             'X3_polynomials_in_branch_and_free_variables':[sorted(p) for p in branched],
             'gray_constants':constants,'gray_changes':changes,'gray_neighbors':neighbors}
    (HERE/'algebra_certificate.json').write_text(json.dumps(algebra,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k not in ['witnesses','planted_targets','scope']},indent=2))

if __name__=='__main__':main()
