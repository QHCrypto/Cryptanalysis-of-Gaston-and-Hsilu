"""Independent bit-vector/DDT verification, including the original C++ round."""
import ctypes
import json
from pathlib import Path
import re
import subprocess
from itertools import product
ROOT=Path(__file__).resolve().parent
REF=ROOT.parents[1]/'hsilu_ref'
MASK=(1<<64)-1
# Independently tabulated standard chi5, bits in row order, row 0 least significant.
SBOX=(0,9,18,11,5,12,22,15,10,3,24,1,13,4,30,7,20,21,6,23,17,16,2,19,26,27,8,25,29,28,14,31)
DDT=[[0]*32 for _ in range(32)]
for x in range(32):
    for y in range(32):DDT[x^y][SBOX[x]^SBOX[y]]+=1

def cols(words):return [sum(((v>>z)&1)<<i for i,v in enumerate(words)) for z in range(64)]
def words(columns):return [sum(((v>>i)&1)<<z for z,v in enumerate(columns)) for i in range(5)]
def hexwords(values):return [f'{x:016x}' for x in values]
def cost(a,b):
    n=DDT[a][b]
    assert n>0 and n&(n-1)==0
    return 5-(n.bit_length()-1)

source=(REF/'hsilu.h').read_text()
PARAM={name:list(map(int,re.search(rf'\b{name}\[[^\]]+\]\s*=\s*\{{([^}}]+)\}}',source)[1].split(',')))
       for name in ('n1','n2','m','r')}

def linear_bits(state):
    # Explicit bit destinations. This does not call a rotate/word-XOR implementation.
    nodes=[[(v>>z)&1 for z in range(64)] for v in state]
    for a,b,m in zip(PARAM['n1'],PARAM['n2'],PARAM['m']):
        nodes.append([nodes[a][z]^nodes[b][(z-m)%64] for z in range(64)])
    return [sum(nodes[13+i][(z-PARAM['r'][i])%64]<<z for z in range(64)) for i in range(5)]

def main():
    wrapper=ROOT/'reference_wrapper.cpp'
    wrapper.write_text('#include "'+str(REF/'hsilu.cpp')+'"\n'+r'''
extern "C" void reference_linear(const uint64_t* in,uint64_t* out){
    int inv[32];
    for(int x=0;x<32;x++){
        int y=0;for(int i=0;i<5;i++)y|=(((x>>i&1)^((1^(x>>((i+1)%5)&1))&(x>>((i+2)%5)&1)))<<i);
        inv[y]=x;
    }
    hsilu_state_t state={};
    for(int z=0;z<64;z++){
        int y=0;for(int i=0;i<5;i++)y|=(in[i]>>z&1)<<i;
        int x=inv[y];for(int i=0;i<5;i++)state.x[i]|=uint64_t(x>>i&1)<<z;
    }
    hsilu_round(&state,0);
    for(int i=0;i<5;i++)out[i]=state.x[i];
}
''')
    library=ROOT/'reference_linear.so'
    subprocess.run(['g++','-std=c++17','-O2','-shared','-fPIC',str(wrapper),'-o',str(library)],check=True)
    lib=ctypes.CDLL(str(library));array=ctypes.c_uint64*5
    lib.reference_linear.argtypes=[ctypes.POINTER(ctypes.c_uint64)]*2
    def reference(a):
        b=array();lib.reference_linear(array(*a),b);return list(b)
    for j in range(320):
        a=[0]*5;a[j//64]=1<<(j%64)
        assert linear_bits(a)==reference(a)

    core=json.loads((ROOT/'results/best_core.json').read_text())
    b0,a1,b1,a2=([int(x,16) for x in core[key]] for key in ('b0','a1','b1','a2'))
    assert linear_bits(b0)==a1==reference(b0)
    assert linear_bits(b1)==a2==reference(b1)
    a0=words([max(range(32),key=lambda a:DDT[a][b]) for b in cols(b0)])
    b2=words([max(range(32),key=lambda b:DDT[a][b]) for a in cols(a2)])
    states={'a0':a0,'b0':b0,'a1':a1,'b1':b1,'a2':a2,'b2':b2}
    local_weights=[[cost(a,b) for a,b in zip(cols(states[f'a{r}']),cols(states[f'b{r}']))] for r in range(3)]
    weights=list(map(sum,local_weights))
    assert weights==core['weights'] and sum(weights)==core['weight']==80
    # Exhaustively recheck the fixed middle-input extension, independently of C++.
    middle=cols(a1)
    positions=[z for z,a in enumerate(middle) if a]
    alternatives=[]
    for z in positions:
        values=[]
        for b in range(32):
            if DDT[middle[z]][b]:
                column_words=[((b>>i)&1)<<z for i in range(5)]
                values.append(linear_bits(column_words))
        alternatives.append(values)
    forward_weights=[5-(max(row).bit_length()-1) for row in DDT]
    count=0;minimum=1000;minimizers=0
    for selected in product(*alternatives):
        candidate=[0]*5
        for image in selected:
            for i in range(5):candidate[i]^=image[i]
        value=sum(forward_weights[a] for a in cols(candidate))
        count+=1
        if value<minimum:minimum=value;minimizers=1
        elif value==minimum:minimizers+=1
    assert count==16384 and minimum==50
    # Closed-form family: R is a 64-bit left rotation, sum is XOR.
    def rotate(t,k):return ((t<<k)|(t>>((-k)%64)))&MASK
    input_shifts=[[4,48],[13,18,54],[0,44],[13,57],[]]
    output_shifts=[[],[10,54],[],[5,28,30,46,48],[]]
    def evaluate(shifts,t):
        values=[]
        for row in shifts:
            v=0
            for k in row:v^=rotate(t,k)
            values.append(v)
        return values
    for z in range(64):
        assert linear_bits(evaluate(input_shifts,1<<z))==evaluate(output_shifts,1<<z)
    assert evaluate(input_shifts,1<<10)==b0
    result={
        'verified':True,'weight_convention':'sum of exact local -log2(DDT/32), free endpoints',
        'scope':'No global optimality or fixed-permutation probability/realizability claim',
        'weights':weights,'total_weight':sum(weights),
        'active_columns':[sum(bool(v) for v in cols(states[f'a{r}'])) for r in range(3)],
        'input_bit_weights':[sum(v.bit_count() for v in states[f'a{r}']) for r in range(3)],
        'states':{key:hexwords(a) for key,a in states.items()},
        'local_weights':local_weights,
        'transition_ddt_counts':[[DDT[a][b] for a,b in zip(cols(states[f'a{r}']),cols(states[f'b{r}']))] for r in range(3)],
        'checks':{'full_ddt_supported_pairs':sum(bool(n) for row in DDT for n in row),
                  'reference_round_linear_basis_vectors':320,'linear_boundaries':2,'chi_columns':192},
        'final_linear_output':hexwords(linear_bits(b2)),
        'fixed_core_extension':{'exhausted':count,'minimum_third_weight':minimum,'minimizers':minimizers},
        'cancellation_family':{'input_shifts':input_shifts,'output_shifts':output_shifts,
                               'selected_parameter_bit':10,'basis_vectors_verified':64},
    }
    seeds=json.loads((ROOT/'cut_seeds.json').read_text())
    result['generating_zero_nodes']=next(x['zero_nodes'] for x in seeds if x['a1']==core['a1'])
    (ROOT/'results/trail80.json').write_text(json.dumps(result,indent=2)+'\n')
    labels=['a0','b0','a1','b1','a2','b2']
    lines=['Hsilu 3-round differential trail: 80 = 16 + 14 + 50',
           'Row order 0..4; bit z is bit z of each hexadecimal word; left rotations.',
           'a[r] --chi--> b[r] --L--> a[r+1]; first and last differences free.',
           'Exact local DDT weight sum; no full-permutation probability claim.','']
    lines += [key+' = '+', '.join(hexwords(states[key])) for key in labels]
    lines += ['',json.dumps({k:result[k] for k in ('weights','active_columns','generating_zero_nodes','checks')})]
    (ROOT/'results/trail80.txt').write_text('\n'.join(lines)+'\n')
    print(json.dumps({k:result[k] for k in ('verified','weights','total_weight','active_columns','input_bit_weights','generating_zero_nodes','checks','fixed_core_extension')},indent=2))

if __name__=='__main__':main()
