"""Generate polynomial syzygies by setting four MRM word signals to zero.

Both raw minors and gcd-reduced minors are retained: removing a common factor
can increase cyclic Hamming weight. This is a seed generator, not exhaustive
over all low-weight states or over all parameters of each family.
"""
from itertools import combinations
from pathlib import Path
import json
import sys
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'structure_probe_20261117'))
import cancellation_probe as p

def mul(a,b):
    if a.bit_count()<b.bit_count():a,b=b,a
    ans=0
    while b:
        bit=b&-b;ans^=a<<(bit.bit_length()-1);b^=bit
    return ans

def divmod_poly(a,b):
    q=0
    while a and a.bit_length()>=b.bit_length():
        k=a.bit_length()-b.bit_length();q^=1<<k;a^=b<<k
    return q,a

def gcd(a,b):
    while b:a,b=b,divmod_poly(a,b)[1]
    return a

def cyc(a):
    out=0
    while a:out^=a&p.MASK;a>>=64
    return out

def minors(rows):
    dp={0:1}
    for row in rows:
        new={}
        for mask,value in dp.items():
            for j,entry in enumerate(row):
                if not(mask>>j&1):new[mask|(1<<j)]=new.get(mask|(1<<j),0)^mul(value,entry)
        dp=new
    return [dp.get(31^(1<<j),0) for j in range(5)]

def main():
    nodes=[[int(i==j) for i in range(5)] for j in range(5)]
    for a,b,s in zip(p.N1,p.N2,p.M):nodes.append([x^p.rol(y,s) for x,y in zip(nodes[a],nodes[b])])
    seeds={};cuts=0;nonzero=0
    for selected in combinations(range(18),4):
        cuts+=1;raw=minors([nodes[i] for i in selected])
        if not any(raw):continue
        nonzero+=1;g=0
        for a in raw:g=gcd(g,a)
        forms=[raw,[divmod_poly(a,g)[0] for a in raw]]
        current=raw
        for k in range(6):
            quot=[divmod_poly(a,3) for a in current]
            if any(rem for _,rem in quot):break
            current=[q for q,_ in quot];forms.append(current)
        for form in forms:
            state=tuple(cyc(a) for a in form)
            if not any(state):continue
            words=p.circuit(state)
            assert all(words[i]==0 for i in selected)
            a1=p.linear(state)
            if p.reverse_weight(state)+p.weight(a1)>79:continue
            c=min(tuple(p.rol(x,r) for x in a1) for r in range(64))
            seeds.setdefault(c,selected)
    (ROOT/'cut_seeds.txt').write_text(''.join(' '.join(f'{x:016x}' for x in a)+'\n' for a in seeds))
    (ROOT/'cut_seeds.json').write_text(json.dumps([{'a1':[f'{x:016x}' for x in a],
                                                 'zero_nodes':selected}
                                                for a,selected in seeds.items()],indent=2)+'\n')
    stats={'cuts':cuts,'nonzero_minor_sets':nonzero,'distinct_canonical_seeds_core_cost_le79':len(seeds)}
    (ROOT/'results/cut_generation.json').write_text(json.dumps(stats,indent=2)+'\n')
    print(stats,flush=True)

if __name__=='__main__':main()
