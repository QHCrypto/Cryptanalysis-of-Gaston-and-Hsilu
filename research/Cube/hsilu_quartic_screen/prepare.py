"""Exact quartic tails, rotation classes, and reduced-size screening inputs."""
import json
import sys
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from permutations import HN1,HN2,HM,HR

def multiply(a,b):
    out=set()
    for u in a:
        for v in b:
            w=u|v
            if w in out: out.remove(w)
            else: out.add(w)
    return out

def linear(state,n):
    x=list(state)
    for u,v,s in zip(HN1,HN2,HM):
        x.append([x[u][j]^x[v][(j-s)%n] for j in range(n)])
    return [[x[13+i][(j-HR[i])%n] for j in range(n)] for i in range(5)]

def quadratic(state,n):
    return [[multiply(state[(i+1)%5][j],state[(i+2)%5][j]) for j in range(n)] for i in range(5)]

def rotate(term,shift,n):
    return tuple(sorted(n*(p//n)+(p%n+shift)%n for p in term))

def canonical(term,n):
    return min(rotate(term,s,n) for s in range(n))

def prepare(n,rounds):
    source=[[{1<<(n*i+j)} for j in range(n)] for i in range(5)]
    before=quadratic(linear(quadratic(source,n),n),n)
    after=linear(before,n)
    classes=set();mapping={}
    for label,state,words in [('before_final_linear',before,range(5)),('observable',after,range(2))]:
        for i in words:
            terms=[tuple(j for j in range(5*n) if mask>>j&1) for mask in state[i][0] if mask.bit_count()==4]
            terms.sort();cs=[canonical(t,n) for t in terms]
            mapping[f'{label}_{i}']={'terms':terms,'classes':cs,'count':len(terms),'class_count':len(set(cs))}
            if label=='observable':classes.update(cs)
    ordered=sorted(classes)
    base=Path(__file__).resolve().parent
    with (base/f'terms_{n}.txt').open('w') as f:
        f.write(f'{n} {rounds} {len(ordered)}\n')
        for t in ordered:f.write(' '.join(map(str,t))+'\n')
    report={'word_bits':n,'prefix_rounds':rounds,'classes':ordered,'mapping':mapping,
            'meaning':'Exact tail quartics before prefix feasibility screening; short rotation periods included.'}
    (base/f'terms_{n}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(n,len(ordered),{k:(v['count'],v['class_count']) for k,v in mapping.items()})

if __name__=='__main__':
    prepare(4,1)
    prepare(64,5)
