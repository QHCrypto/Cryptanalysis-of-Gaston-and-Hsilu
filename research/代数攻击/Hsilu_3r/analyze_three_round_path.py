"""Read the SVG path and validate its actual Hsilu conditions.

Path bit j is native integer bit 63-j. Xr is round input, Yr chi output.
Numerical scope: the last three rounds (69,5a,4b), rate row 0.
"""
import collections
import ctypes
import itertools
import json
from pathlib import Path
import random
import sys
import xml.etree.ElementTree as ET

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'hsilu_2r'))
import verify_path_algebra as ref

def read_svg():
    grid = {name:[[None]*64 for _ in range(5)] for name in ['X1','Y1','X2','Y2','X3','Y3','X4']}
    root = ET.parse(HERE/'hsilu_3round_best_attack.svg').getroot()
    for e in root.findall('{http://www.w3.org/2000/svg}rect'):
        if e.get('width') != '7' or e.get('height') != '7':
            continue
        x,y=int(e.get('x')),int(e.get('y'))
        stage=(y-78)//78; row=((y-78)%78)//8; col=(x-86)//8
        grid[list(grid)[stage]][row][col]=e.get('fill')
    assert all(v for a in grid.values() for row in a for v in row)
    return grid

GRID=read_svg()
U=[j for j,c in enumerate(GRID['X1'][0]) if c=='#ffd84d']
G=[j for j in range(64) if j not in U]
COND={name:[(i,j,{'#d73a49':0,'#2da44e':1}[c]) for i,row in enumerate(GRID[name]) for j,c in enumerate(row) if c in ['#d73a49','#2da44e']] for name in ['X1','X2','X3']}
SELECT=[j for j,c in enumerate(GRID['X4'][0]) if c=='#8ec5ff']

def system(rows, rhs, n):
    pivots={}
    for row,b in zip(rows,rhs):
        while row:
            p=row.bit_length()-1
            if p in pivots:
                row^=pivots[p][0];b^=pivots[p][1]
            else:
                pivots[p]=(row,b);break
        if not row and b:
            return None,[],len(pivots)
    def complete(free, homogeneous):
        answer=free
        for p in sorted(pivots):
            row,b=pivots[p]
            bit=((row&answer).bit_count()&1)^(0 if homogeneous else b)
            answer|=bit<<p
        return answer
    return complete(0,False),[complete(1<<i,True) for i in range(n) if i not in pivots],len(pivots)

def first_affine(state):
    base=ref.rnd(state,0x69)
    cols=[]
    for j in G:
        b=state.copy();b[0]^=1<<(63-j)
        cols.append([x^y for x,y in zip(base,ref.rnd(b,0x69))])
    aff=[[((base[i]>>(63-j))&1) | sum(((c[i]>>(63-j))&1)<<(k+1) for k,c in enumerate(cols)) for j in range(64)] for i in range(5)]
    for j in U:
        b=state.copy();b[0]^=1<<(63-j)
        mid=ref.rnd(b,0x69)
        assert all((base[i]^mid[i])>>(63-z)&1==0 for i,z,value in COND['X2'])
    return aff

def restrict_affine(a, particular, basis):
    return (a&1)^((a>>1&particular).bit_count()&1) | sum(((a>>1&v).bit_count()&1)<<(k+1) for k,v in enumerate(basis))

def linear_symbolic(a):
    # MSB-indexed rows, so left rotation reads j+s.
    x=list(a)
    for j in range(5,18):
        x.append([x[ref.PARAM['n1'][j]][z] ^ x[ref.PARAM['n2'][j]][(z+ref.PARAM['m'][j])%64] for z in range(64)])
    return [[x[13+i][(j+ref.PARAM['s'][i])%64] for j in range(64)] for i in range(5)]

def affine_poly(a):
    s={0} if a&1 else set()
    a>>=1
    while a:
        b=a&-a;s.add(b);a^=b
    return s

def multiply(p,q):
    out=set()
    for a in p:
        for b in q:
            c=a|b
            if c in out:out.remove(c)
            else:out.add(c)
    return out

def second_condition_polys(mid):
    a=[[affine_poly(v) for v in row] for row in mid]
    for j in range(64):
        if 0x5a>>(63-j)&1:a[0][j]^={0}
    y=[[a[i][j] ^ multiply(a[(i+1)%5][j]^{0},a[(i+2)%5][j]) for j in range(64)] for i in range(5)]
    x=linear_symbolic(y)
    return [x[i][j]^({0} if value else set()) for i,j,value in COND['X3']]

SOURCES=linear_symbolic([[{(i,j)} for j in range(64)] for i in range(5)])

def product_edges(mid):
    edges=set();origins=collections.defaultdict(list)
    bycondition=[]
    for row,j,value in COND['X3']:
        condition_edges=set()
        for i,z in SOURCES[row][j]:
            b=mid[(i+1)%5][z]>>1;c=mid[(i+2)%5][z]>>1
            if not b or not c or b==c:continue
            pair=tuple(sorted((b,c)))
            if pair in condition_edges:condition_edges.remove(pair)
            else:condition_edges.add(pair)
            origins[b].append([(i+1)%5,z]);origins[c].append([(i+2)%5,z])
        edges|=condition_edges
        bycondition.append(condition_edges)
    return edges,origins,bycondition

def basis_reduce(v,basis):
    while v:
        p=v.bit_length()-1
        if p not in basis:break
        v^=basis[p]
    return v

def greedy_span_cover(edges,trials=100):
    best=None
    nodes=set(itertools.chain.from_iterable(edges))
    for trial in range(trials):
        basis={};selected=[];remaining=set(edges)
        while remaining:
            candidates=set(itertools.chain.from_iterable(remaining))
            equivalence=collections.defaultdict(set)
            for node in nodes:
                residue=basis_reduce(node,basis)
                if residue:equivalence[residue].add(node)
            choices=[]
            for residue,group in equivalence.items():
                score=sum(a in group or b in group for a,b in remaining)
                if score:choices.append((score+random.random()*(0 if trial==0 else 3),residue,next(iter(group))))
            _,residue,node=max(choices)
            basis[residue.bit_length()-1]=residue;selected.append(node)
            remaining={(a,b) for a,b in remaining if basis_reduce(a,basis) and basis_reduce(b,basis)}
        if best is None or len(selected)<len(best):best=selected
    return best

def main():
    random.seed(20261013)
    rank_hist=collections.Counter();inconsistent=0;samples=[]
    for trial in range(32):
        state=[0]+[random.getrandbits(64) for _ in range(4)]
        for i,j,value in COND['X1']:
            state[i]=(state[i]&~(1<<(63-j)))|(value<<(63-j))
        mid=first_affine(state)
        rows=[mid[i][j]>>1 for i,j,value in COND['X2']]
        rhs=[(mid[i][j]&1)^value for i,j,value in COND['X2']]
        g0,basis,rk=system(rows,rhs,len(G))
        rank_hist[rk]+=1
        if g0 is None:
            inconsistent+=1;continue
        restricted=[[restrict_affine(a,g0,basis) for a in row] for row in mid]
        polys=second_condition_polys(restricted)
        edges,origins,bycondition=product_edges(restricted)
        cover=greedy_span_cover(edges,100 if trial==0 else 10)
        samples.append({'trial':trial,'X2_rank':rk,'remaining_background_dimension':len(basis),
                        'X3_quadratic_monomial_counts':[sum(m.bit_count()==2 for m in p) for p in polys],
                        'X3_degrees':[max((m.bit_count() for m in p),default=0) for p in polys],
                        'distinct_product_edges':len(edges),'span_cover_rank':len(cover),
                        'cover_forms_hex':[hex(x) for x in cover],
                        'cover_original_X2_positions':[origins[x][0] for x in cover]})
        if trial==0:
            detailed={'native_capacity_words_hex':[f'{w:016x}' for w in state[1:]],
                      'background_particular_hex':hex(g0),'background_basis_hex':[hex(x) for x in basis],
                      'X3_residual_polynomials_monomial_masks': [sorted(p) for p in polys],
                      'reduced_X2_affine_forms_hex':[[hex(a) for a in row] for row in restricted]}
            (HERE/'first_capacity_algebra_details.json').write_text(json.dumps(detailed,indent=2)+'\n')
    out={'active_positions':U,'background_count':len(G),'selected_output_positions':SELECT,
         'conditions':COND,'X2_rank_histogram':dict(rank_hist),'inconsistent_linear_conditions':inconsistent,
         'samples':samples}
    (HERE/'three_round_condition_analysis.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k!='samples'},indent=2))
    print('SAMPLES',json.dumps(samples,indent=2))

if __name__=='__main__':main()
